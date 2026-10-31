# Pipeline Slurm — BIO-DYN-projet

Architecture Scatter-Gather pour lancer les simulations Monte Carlo sur le cluster Genouest, en balayant les paramètres **h**, **lX** et la largeur de la **bande d'iridophores** (et plusieurs graines).

## Prérequis

- Accès au cluster avec Slurm
- Environnement conda `biodyn` accessible via `/local/env/envconda.sh` (création une seule fois : `conda create -n biodyn python numpy matplotlib`)
- Se placer à la **racine du projet** pour toutes les commandes

## Utilisation

### Étape 1 — Générer les runs et les découper en chunks

```bash
# Défauts (pseudocode.md) : L=200, h de 2 à 20 pas 1, lX de 1.1 à 4.0 pas 0.1, pas d'iridophores, 3 graines
bash slurm_pipeline/01_prepare_chunks.sh

# Plages personnalisées (début / fin / pas)
bash slurm_pipeline/01_prepare_chunks.sh --h-start 4 --h-end 16 --h-step 4 \
    --lx-start 1.5 --lx-end 3.5 --lx-step 1 --graines 1

# Un seul run
bash slurm_pipeline/01_prepare_chunks.sh --L 100 --h-start 16 --h-end 16 \
    --lx-start 2.5 --lx-end 2.5 --graines 1 --n-pas 1000000000

# Avec iridophores (Fig. 1b–c) : bandes de 1 et 10 lignes
bash slurm_pipeline/01_prepare_chunks.sh --L 100 --h-start 16 --h-end 16 \
    --lx-start 2.5 --lx-end 2.5 --bande-start 1 --bande-end 10 --bande-step 9
```

Crée `params_chunks/chunk_1.txt`, `chunk_2.txt`, … : une ligne `L h lX bande graine n_pas` par run.
La taille des chunks est configurable via `CHUNK_SIZE=10` en tête du script.

**Options :**

| Option | Défaut | Description |
|--------|--------|-------------|
| `--L` | `200` | Taille de la grille. |
| `--h-start` | `2` | Première valeur de h. |
| `--h-end` | `20` | Dernière valeur de h, incluse. |
| `--h-step` | `1` | Pas du balayage de h. |
| `--lx-start` | `1.1` | Première valeur de lX. |
| `--lx-end` | `4.0` | Dernière valeur de lX, incluse. |
| `--lx-step` | `0.1` | Pas du balayage de lX. |
| `--bande-start` | `0` | Première largeur de bande d'iridophores, en lignes (0 = pas d'iridophores). |
| `--bande-end` | `0` | Dernière largeur, incluse. |
| `--bande-step` | `1` | Pas du balayage de la bande. |
| `--graines` | `3` | Nombre de répétitions (graines 0 à N-1). |
| `--n-pas` | `100000 × L²` | Nombre de pas par run. |

### Étape 2 — Soumettre les jobs

```bash
bash slurm_pipeline/03_submit_all.sh
```

Soumet **un seul Job Array CPU**, une tâche par chunk. Surveiller l'avancement :

```bash
squeue -u $USER
```

**Relancer** : refaire la même commande. Les runs terminés sont sautés, les runs interrompus (limite de 24 h, nœud perdu) reprennent depuis leur dernière sauvegarde.

### Étape 3 — Récupérer les paramètres des runs

```bash
python slurm_pipeline/04_merge_results.py
# options : --out-dir DOSSIER  --output CSV
```

Produit `results/final_runs.csv` : une ligne par run (`L, h, lX, bande, graine, n_pas, pas_faits, t, termine, dossier`). Les runs avec `termine = False` sont à relancer (étape 2).

---

## Architecture des scripts

| Fichier | Rôle |
|---------|------|
| `01_prepare_chunks.sh` | Génère les runs des plages et les chunks (`params_chunks/`) |
| `02_run_worker.slurm` | Worker exécuté par chaque tâche Slurm : lance les runs de son chunk |
| `03_submit_all.sh` | Soumission du Job Array (et relance) |
| `04_merge_results.py` | Récupère les paramètres et l'état de tous les runs dans un CSV |

## Paramètres configurables

| Fichier | Variable | Défaut | Description |
|---------|----------|--------|-------------|
| `01_prepare_chunks.sh` | `CHUNK_SIZE` | `10` | Runs par chunk |
| `01_prepare_chunks.sh` | `DEFAULT_*` | voir options | Plages par défaut |
| `03_submit_all.sh` | `MAX_CONCURRENT_CPU` | `40` | Tâches simultanées |

## Structure des sorties

```
params_chunks/
├── chunk_1.txt                            # une ligne "L h lX bande graine n_pas" par run
└── ...

out/
└── L200_h2_lX1.1_b0_g0/                   # un dossier par run (simulation.py)
    ├── grille.npy
    ├── params.json
    └── grille.png

results/
└── final_runs.csv                         # paramètres et état de tous les runs (étape 3)

logs/
└── slurm_<JOB_ID>_<TASK_ID>.out          # Logs Slurm par tâche
```

## Iridophores

Ils ne sont pas simulés comme des cellules : ils sont représentés par un masque de cases où le noir ne peut pas naître (`masque_iridophores` dans `simulation.py`). Aujourd'hui c'est une bande horizontale de `bande` lignes au milieu de la grille (Fig. 1b–c). Pour une autre forme, seule cette fonction est à modifier.

## Mapping TASK_ID → chunk

Le Job Array a une tâche par chunk :

```
CHUNK_ID = TASK_ID
```

Exemple avec les défauts : 19 valeurs de h × 30 valeurs de lX × 1 bande × 3 graines = 1710 runs → 171 chunks de 10 runs → `--array=1-171%40`.
