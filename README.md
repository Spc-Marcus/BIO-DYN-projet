# BIO-DYN-projet

Reproduction de Bullara & De Decker (2015), *Pigment cell movement is not required for generation of Turing patterns in zebrafish skin*, Nat. Commun. 6:6971.

Notes sur l'article : [projet.md](projet.md). Pseudo-code de la simulation : [pseudocode.md](pseudocode.md).

## Fichiers

Base du projet : les fonctions sont déclarées et documentées, pas encore implémentées (elles lèvent `NotImplementedError`).

| Fichier | Rôle |
|---|---|
| `simulation.py` | Un run Monte Carlo : paramètres, simulation, sauvegarde et reprise |
| `slurm_pipeline/` | Lancement sur GenOuest en Job Array : scripts, [README](slurm_pipeline/README.md) et documentation avec diagrammes dans [slurm_pipeline/doc/](slurm_pipeline/doc/) |

## Lancer sur GenOuest

Depuis la racine du projet :

```bash
bash slurm_pipeline/01_prepare_chunks.sh --h-start 2 --h-end 20 --h-step 2   # runs et chunks
bash slurm_pipeline/03_submit_all.sh                                         # Job Array (et relance)
python slurm_pipeline/04_merge_results.py                                    # paramètres et état des runs
```

Exemple complet pas à pas : [slurm_pipeline/doc/exemple.md](slurm_pipeline/doc/exemple.md).
