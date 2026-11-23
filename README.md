# BIO-DYN-projet

Reproduction de Bullara & De Decker (2015), *Pigment cell movement is not required for generation of Turing patterns in zebrafish skin*, Nat. Commun. 6:6971.

Notes sur l'article : [projet.md](projet.md). Pseudo-code de la simulation : [pseudocode.md](pseudocode.md).

## Fichiers

Code du projet : un run Monte Carlo (`simulation.py`), la timeline d'un run (`timeline.py`) et le pipeline de lancement sur GenOuest (`slurm_pipeline/`).

| Fichier | Rôle |
|---|---|
| `simulation.py` | Un run Monte Carlo : paramètres, simulation, sauvegarde et reprise |
| `timeline.py` | GIF et planche « timeline » d'un run lancé avec `--historique` (indépendant de Slurm) |
| `slurm_pipeline/` | Lancement sur GenOuest en Job Array : scripts, [README](slurm_pipeline/README.md) et documentation avec diagrammes dans [slurm_pipeline/doc/](slurm_pipeline/doc/) |

## Lancer un run en local

Avec numpy et matplotlib installés, depuis la racine du projet :

```bash
python simulation.py --L 50 --h 8 --lX 2.5 --n_pas 5000000 --historique   # environ 5 s
python timeline.py out/L50_h8_lX2.5_b0_g0                                  # GIF et planche
```

Le run écrit `grille.npy`, `params.json` et `grille.png` dans `out/L50_h8_lX2.5_b0_g0/`, et avec `--historique` une grille tous les `PAS_ENTRE_SAUVEGARDES` pas dans `historique/`, plus une toutes les 5 u.t. pendant 500 u.t. au départ et après chaque blessure. Relancer la même commande reprend un run interrompu et ne refait pas un run terminé.

## Options du modèle

Les règles sont dans un dictionnaire règle → taux construit à partir des options (`--bX`, `--bM`, `--dX`, `--dM`, `--s`, `--lX`, `--lB`) ; seules les règles de taux non nul sont tirées. Par défaut, c'est le modèle de l'article. Trois extensions, désactivées par défaut :

| Options | Défaut | Rôle |
|---|---|---|
| `--reseau hex` | `carre` | réseau hexagonal (6 voisins), image en briques |
| `--contact`, `--rho0`, `--fX0` | non, 0, 1 | naissances seulement au contact d'une case occupée, en partant d'une fraction `rho0` de cases occupées (dont une part `fX0` de jaunes) |
| `--lB`, `--alpha`, `--tB` | 0, 1, aucun | blessure (9) de taux `lB` et de rayon `alpha·h` ; avec `--tB`, une seule blessure au centre au temps `tB` |

Chaque option différente de son défaut ajoute un suffixe au dossier du run :

```bash
python simulation.py --L 50 --h 8 --lX 2.5 --n_pas 5000000 --reseau hex --contact --rho0 0.01   # out/L50_h8_lX2.5_b0_g0_hex_c_r0.01
```

## Lancer sur GenOuest

Depuis la racine du projet :

```bash
bash slurm_pipeline/01_prepare_chunks.sh                                     # runs de la Fig. 2 et chunks
# options du modèle pour tous les runs : --extra "--reseau hex --contact --rho0 0.01"
bash slurm_pipeline/03_submit_all.sh                                         # Job Array (et relance)
python slurm_pipeline/04_merge_results.py                                    # paramètres et état des runs
bash slurm_pipeline/05_tests_extensions.sh                                   # tests des extensions, avec GIF
```

Valeurs des tests et raisons : [slurm_pipeline/doc/tests_extensions.md](slurm_pipeline/doc/tests_extensions.md).

Exemple complet pas à pas : [slurm_pipeline/doc/exemple.md](slurm_pipeline/doc/exemple.md).
