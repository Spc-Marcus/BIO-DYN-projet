# Architecture — pipeline Slurm

Composants du pipeline, fichiers qu'ils lisent et écrivent, et découpage du travail en tâches.

## 1. Vue d'ensemble

```mermaid
graph TD
    subgraph Etape1["Étape 1 — Préparation (frontale)"]
        P1(["01_prepare_chunks.sh<br/>plages h, lX, bande, graines"])
        C1[("params_chunks/chunk_1.txt")]
        C2[("params_chunks/chunk_2.txt")]
        CN[("params_chunks/chunk_N.txt")]
        P1 --> C1 & C2 & CN
    end

    subgraph Etape2["Étape 2 — Calcul (nœuds du cluster)"]
        S(["03_submit_all.sh"])
        JA{{"Job Array Slurm<br/>--array=1-N%40"}}
        W1["Tâche 1<br/>02_run_worker.slurm"]
        W2["Tâche 2<br/>02_run_worker.slurm"]
        WN["Tâche N<br/>02_run_worker.slurm"]
        SIM["simulation.py<br/>un run par ligne du chunk"]
        S -->|"sbatch"| JA
        JA --> W1 & W2 & WN
        W1 & W2 & WN --> SIM
    end

    subgraph Sorties["Sorties"]
        OUT[("out/L{L}_h{h}_lX{lX}_b{bande}_g{graine}/<br/>grille.npy · params.json · grille.png")]
        LOGS[("logs/slurm_JOBID_TASKID.out")]
    end

    subgraph Etape3["Étape 3 — Récupération (frontale)"]
        M(["04_merge_results.py"])
        CSV[("results/final_runs.csv")]
        M --> CSV
    end

    C1 -.->|"lu par"| W1
    C2 -.->|"lu par"| W2
    CN -.->|"lu par"| WN
    SIM --> OUT
    W1 & W2 & WN --> LOGS
    OUT -->|"params.json"| M
    CSV -.->|"runs non terminés : relancer l'étape 2"| S

    style Etape1 fill:#e6f5ff,stroke:#4a90c2
    style Etape2 fill:#e6ffe6,stroke:#4caf50
    style Etape3 fill:#fff5e6,stroke:#f0a500
    style Sorties fill:#f5f5f5,stroke:#999
```

| Fichier | Où il tourne | Rôle |
|---|---|---|
| `01_prepare_chunks.sh` | frontale | Génère une ligne `L h lX bande graine n_pas` par run et les répartit en chunks de `CHUNK_SIZE` lignes |
| `03_submit_all.sh` | frontale | Compte les chunks et soumet un seul Job Array, une tâche par chunk |
| `02_run_worker.slurm` | nœud de calcul | Lance l'un après l'autre les runs de son chunk |
| `simulation.py` | nœud de calcul | Un run : simulation, sauvegarde régulière, reprise |
| `04_merge_results.py` | frontale | Rassemble les `params.json` de tous les runs dans un CSV |

## 2. Découpage en tâches : TASK_ID → chunk

Toutes les combinaisons des plages sont écrites, une par ligne, puis découpées en paquets de `CHUNK_SIZE` lignes (10 par défaut). La tâche `N` du Job Array traite `params_chunks/chunk_N.txt` :

```
CHUNK_ID = TASK_ID
nombre de tâches = ⌈ nombre de runs / CHUNK_SIZE ⌉
```

```mermaid
flowchart LR
    R["Plages<br/>h : start → end pas step<br/>lX : start → end pas step<br/>bande : start → end pas step<br/>graines : 0 à G-1"]
    X["Produit cartésien<br/>une ligne par run"]
    C1["chunk_1.txt<br/>lignes 1 à 10"]
    C2["chunk_2.txt<br/>lignes 11 à 20"]
    CN["chunk_N.txt<br/>lignes restantes"]
    T1["tâche 1"]
    T2["tâche 2"]
    TN["tâche N"]

    R --> X
    X --> C1 --> T1
    X --> C2 --> T2
    X --> CN --> TN
```

Ordre des lignes : h varie le plus lentement, puis lX, puis bande, puis la graine. Avec les valeurs par défaut (19 h × 30 lX × 1 bande × 3 graines = 1710 runs), on obtient 171 chunks, soit `--array=1-171%40` : bien en dessous de la limite de taille des Job Arrays.

## 3. Un run : `simulation.py`

`simulation.py` ne connaît pas Slurm : il reçoit les paramètres d'un run en ligne de commande, ce qui permet aussi de le lancer à la main en local.

```mermaid
flowchart TD
    A["lire_arguments<br/>L, h, lX, bande, graine, n_pas"]

    subgraph Prep["preparer_run"]
        B{"out/run/params.json<br/>existe ?"}
        C["nouveau run<br/>initialiser_grille, générateur créé depuis la graine"]
        D["charger_run<br/>grille, pas_faits, état du générateur"]
        B -->|"non"| C
        B -->|"oui"| D
    end

    E{"pas_faits ≥ n_pas ?"}
    Z(["Run déjà terminé, rien à faire"])

    subgraph Boucle["simulation_peau"]
        MK["masque_iridophores(L, bande)"]
        T["tirer une case (i, j)<br/>tirer_regle"]
        AP["appliquer_regle<br/>voisin_au_hasard, case_a_distance, masque"]
        SV{"PAS_ENTRE_SAUVEGARDES<br/>pas depuis la dernière sauvegarde ?"}
        SR["sauvegarder_run"]
        FIN{"pas_faits < n_pas ?"}
        MK --> T --> AP --> SV
        SV -->|"oui"| SR --> FIN
        SV -->|"non"| FIN
        FIN -->|"oui"| T
    end

    F["sauvegarder_run final"]

    A --> B
    C & D --> E
    E -->|"oui"| Z
    E -->|"non"| MK
    FIN -->|"non"| F
```

## 4. Iridophores

Les iridophores ne sont pas simulés comme des cellules : ils sont représentés par un masque de cases où le noir ne peut pas naître (`masque_iridophores` dans `simulation.py`). La règle (8) y est bloquée, ainsi que la règle (3) si bM devenait non nul.

```mermaid
flowchart LR
    BD["--bande N"] --> MK["masque_iridophores(L, bande)"]
    MK --> NI[("noir_interdit<br/>tableau L × L de booléens")]
    NI --> AP["appliquer_regle<br/>règle (8) bloquée sur les cases True"]
```

Aujourd'hui le masque est une bande horizontale de `bande` lignes au milieu de la grille (Fig. 1b–c), 0 = pas d'iridophores. Pour une autre forme, seule `masque_iridophores` est à modifier : le reste du code ne lit que le masque.

## Principes clés

- **Une tâche = un chunk**, comme dans le pipeline de BenchQBR : le nombre de tâches reste faible même pour de gros balayages, et `MAX_CONCURRENT_CPU` limite le nombre de tâches simultanées.
- **Un run = un dossier** `out/L{L}_h{h}_lX{lX}_b{bande}_g{graine}/`. Son nom ne dépend pas de `n_pas` : relancer avec un `n_pas` plus grand prolonge le run au lieu d'en créer un nouveau.
- **Reprise automatique** : le run est sauvegardé tous les `PAS_ENTRE_SAUVEGARDES` pas, avec l'état du générateur aléatoire. Relancer `03_submit_all.sh` reprend les runs interrompus au dernier point sauvegardé et saute les runs terminés.
- **Chemins relatifs à la racine du projet** : `03_submit_all.sh` passe `--chdir` à Slurm, les tâches écrivent donc toujours dans `out/` et `logs/` à la racine.
- **Ne pas relancer `01_prepare_chunks.sh` tant que des tâches sont en attente** : il efface et réécrit `params_chunks/`, et chaque tâche ne lit son chunk qu'au moment où elle démarre.
