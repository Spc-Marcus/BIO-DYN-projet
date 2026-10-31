# BIO-DYN-projet

Reproduction de Bullara & De Decker (2015), *Pigment cell movement is not required for generation of Turing patterns in zebrafish skin*, Nat. Commun. 6:6971.

Notes sur l'article : [projet.md](projet.md). Pseudo-code de la simulation : [pseudocode.md](pseudocode.md).

## Fichiers

Base du projet : les fonctions sont déclarées et documentées, pas encore implémentées (elles lèvent `NotImplementedError`).

| Fichier | Rôle |
|---|---|
| `simulation.py` | Un run Monte Carlo : paramètres, simulation, sauvegarde et reprise |
| `lancer.sh` | Lance ou relance les runs sur GenOuest (job array, plages de paramètres) |
| `run.slurm` | Une tâche du job array : une combinaison de paramètres |

## Lancer sur GenOuest

Chaque paramètre (`--L`, `--h`, `--lX`, `--graine`, `--n_pas`) prend une valeur seule ou, sauf `--n_pas`, une plage `MIN:MAX:PAS` (bornes incluses). Une tâche est lancée par combinaison.

```bash
bash lancer.sh --h 16 --lX 2.5                                # un seul run
bash lancer.sh --h 2:20:2 --lX 1.1:4.0:0.1 --graine 0:2:1     # 10 x 30 x 3 = 900 runs
```

Pour relancer, refaire la même commande : les runs terminés sont sautés, les runs interrompus reprennent depuis leur dernière sauvegarde.

Chaque run écrit `grille.npy`, `params.json` et `grille.png` dans `out/L{L}_h{h}_lX{lX}_g{graine}/`, et son journal dans `logs/`.
