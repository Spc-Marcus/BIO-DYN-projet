# BIO-DYN-projet

Reproduction de Bullara & De Decker (2015), *Pigment cell movement is not required for generation of Turing patterns in zebrafish skin*, Nat. Commun. 6:6971.

Notes sur l'article : [projet.md](projet.md). Pseudo-code de la simulation : [pseudocode.md](pseudocode.md).

## Fichiers

Base du projet : les fonctions sont déclarées et documentées, pas encore implémentées (elles lèvent `NotImplementedError`).

| Fichier | Rôle |
|---|---|
| `simulation.py` | Un run Monte Carlo : paramètres, simulation, sauvegarde et reprise |
| `run.slurm` | Lance ou relance un run sur GenOuest |

## Lancer sur GenOuest

Depuis la racine du projet :

```bash
mkdir -p logs
sbatch run.slurm --h 16 --lX 2.5                     # nouveau run
sbatch run.slurm --reprise out/L100_h16_lX2.5_g0     # relance, paramètres récupérés
```

Chaque run écrit `grille.npy`, `params.json` et `grille.png` dans son dossier sous `out/`.
