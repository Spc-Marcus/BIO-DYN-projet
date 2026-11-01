# Diagrammes de séquence — pipeline Slurm

## 1. Flux complet (cas nominal)

```bash
bash slurm_pipeline/01_prepare_chunks.sh --h-start 5 --h-end 15 --h-step 5 --lx-start 0.5 --lx-end 3.5 --lx-step 1
bash slurm_pipeline/03_submit_all.sh
python slurm_pipeline/04_merge_results.py
```

```mermaid
sequenceDiagram
    actor User as Utilisateur
    participant Prep as 01_prepare_chunks.sh
    participant FS as Disque<br/>(racine du projet)
    participant Submit as 03_submit_all.sh
    participant Slurm
    participant Worker as 02_run_worker.slurm<br/>tâche N
    participant Sim as simulation.py
    participant Merge as 04_merge_results.py

    User->>Prep: plages h, lX, bande, graines
    rect rgb(230, 245, 255)
        Note over Prep,FS: Étape 1 — Préparation
        Prep->>FS: efface params_chunks/
        Prep->>FS: écrit chunk_1.txt … chunk_N.txt
        Prep-->>User: N (nombre de chunks)
    end

    User->>Submit: bash slurm_pipeline/03_submit_all.sh
    rect rgb(230, 255, 230)
        Note over Submit,Sim: Étape 2 — Calcul
        Submit->>FS: compte les chunks
        Submit->>Slurm: sbatch --array=1-N%40 02_run_worker.slurm
        Slurm-->>Submit: Submitted batch job JOBID
        loop pour chaque tâche, 40 en même temps au plus
            Slurm->>Worker: démarre la tâche N sur un nœud
            Worker->>FS: lit params_chunks/chunk_N.txt
            loop pour chaque ligne du chunk
                Worker->>Sim: python simulation.py --L --h --lX --bande --graine --n_pas
                Sim->>FS: sauvegardes régulières dans out/run/
                Sim-->>Worker: code de sortie
            end
            Worker->>FS: journal logs/slurm_JOBID_N.out
        end
    end

    User->>Merge: python slurm_pipeline/04_merge_results.py
    rect rgb(255, 245, 230)
        Note over Merge,FS: Étape 3 — Récupération
        Merge->>FS: lit out/*/params.json
        Merge->>FS: écrit results/final_runs.csv
        Merge-->>User: nombre de runs terminés et à relancer
    end
```

## 2. Reprise après la limite de temps

Une tâche est tuée par Slurm au bout de 24 h (`#SBATCH --time`). On relance avec la **même** commande `03_submit_all.sh`, sans repasser par l'étape 1.

```mermaid
sequenceDiagram
    actor User as Utilisateur
    participant Slurm
    participant Worker as 02_run_worker.slurm<br/>tâche 1
    participant Sim as simulation.py
    participant FS as out/L100_h5_lX3.5_b0_g0/

    rect rgb(255, 235, 235)
        Note over Slurm,FS: Premier passage
        Slurm->>Worker: démarre la tâche 1
        Worker->>Sim: runs 1 à 6 du chunk, terminés
        Worker->>Sim: run 7 : h=5, lX=3.5, graine=0
        loop tous les PAS_ENTRE_SAUVEGARDES pas
            Sim->>FS: grille.npy, params.json (pas_faits, état du générateur)
        end
        Slurm-xWorker: limite de 24 h atteinte, tâche tuée
    end

    User->>Slurm: bash slurm_pipeline/03_submit_all.sh (même commande)

    rect rgb(230, 255, 230)
        Note over Slurm,FS: Relance
        Slurm->>Worker: démarre la tâche 1
        Worker->>Sim: runs 1 à 6
        Sim-->>Worker: Run déjà terminé, rien à faire
        Worker->>Sim: run 7
        Sim->>FS: lit params.json et grille.npy
        Note over Sim: reprend à pas_faits<br/>avec le même état du générateur
        Sim->>FS: sauvegarde finale
        Worker->>Sim: runs 8 à 10
    end
```

## 3. États d'un run

```mermaid
stateDiagram-v2
    state "Nouveau" as Nouveau
    state "En cours" as EnCours
    state "Interrompu" as Interrompu
    state "Terminé" as Termine

    [*] --> Nouveau: pas de params.json
    Nouveau --> EnCours: simulation_peau
    EnCours --> EnCours: sauvegarde tous les PAS_ENTRE_SAUVEGARDES pas
    EnCours --> Interrompu: limite de temps, nœud perdu
    Interrompu --> EnCours: relance de 03_submit_all.sh, reprise à pas_faits
    EnCours --> Termine: pas_faits = n_pas
    Termine --> EnCours: relance avec un n_pas plus grand
    Termine --> [*]
```

Dans `results/final_runs.csv` (étape 3), un run **Terminé** a `termine = True` ; un run **Interrompu** a `termine = False` et doit être relancé.
