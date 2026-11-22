# BIO-DYN-projet

Reproduction de Bullara & De Decker (2015), *Pigment cell movement is not required for generation of Turing patterns in zebrafish skin*, Nat. Commun. 6:6971.

Notes sur l'article : [projet.md](projet.md). Pseudo-code de la simulation : [pseudocode.md](pseudocode.md).

## Fichiers

Code du projet : un run Monte Carlo (`simulation.py`), la timeline d'un run (`timeline.py`) et le pipeline de lancement sur GenOuest (`slurm_pipeline/`).

| Fichier | Rôle |
|---|---|
| `simulation.py` | Un run Monte Carlo : paramètres et options du modèle, simulation, sauvegarde et reprise |
| `timeline.py` | GIF et planche « timeline » d'un run lancé avec `--historique` (indépendant de Slurm) |
| `slurm_pipeline/` | Lancement sur GenOuest en Job Array : scripts, [README](slurm_pipeline/README.md) et documentation avec diagrammes dans [slurm_pipeline/doc/](slurm_pipeline/doc/) |

## Lancer un run en local

Avec numpy et matplotlib installés, depuis la racine du projet :

```bash
python simulation.py --L 50 --h 8 --lX 2.5 --n_pas 5000000 --historique   # environ 5 s
python timeline.py out/L50_h8_lX2.5_b0_g0                                  # GIF et planche
```

Le run écrit `grille.npy`, `params.json` et `grille.png` dans `out/L50_h8_lX2.5_b0_g0/`, et une grille tous les `PAS_ENTRE_SAUVEGARDES` pas dans `historique/` avec `--historique`. Relancer la même commande reprend un run interrompu et ne refait pas un run terminé.

## Options du modèle

Les règles (2) à (8) de l'article sont rangées dans un dictionnaire règle → taux construit à partir des options ; seules les règles de taux non nul sont tirées. Trois extensions s'y ajoutent : la blessure (règle (9)), la naissance au contact et le réseau hexagonal. Avec les valeurs par défaut, on retrouve le modèle de l'article et exactement les mêmes grilles qu'avant ces options.

| Option | Défaut | Suffixe du dossier (si ≠ défaut) | Rôle |
|---|---|---|---|
| `--bX`, `--bM` | 1, 0 | `_bX…`, `_bM…` | taux de naissance spontanée du jaune (2), du noir (3) |
| `--dX`, `--dM` | 0, 0 | `_dX…`, `_dM…` | taux de mort du jaune (4), du noir (5) |
| `--s` | 1 | `_s…` | taux d'élimination entre voisins (6) et (7) |
| `--reseau {carre, hex}` | `carre` | `_hex` | réseau carré (4 voisins) ou hexagonal (6 voisins, image en briques) |
| `--contact` | désactivé | `_c` | naissances (2), (3) et (8) seulement si un voisin direct est occupé ; demande `--rho0` > 0 |
| `--rho0` | 0 | `_r…` | fraction de cases occupées au départ |
| `--fX0` | 1 | `_fx…` (si `rho0` > 0) | part de jaunes parmi ces cases |
| `--lB` | 0 | `_B…` | taux de la blessure (9) : la case tirée et toutes celles à distance ≤ R deviennent vides |
| `--alpha` | 1 | `_a…` (s'il y a une blessure) | rayon de la blessure R = alpha·h |
| `--tB` | aucun | `_tB…` | une seule blessure au centre au temps tB, grille d'avant gardée dans `grille_avant_blessure.npy` |

Par exemple :

```bash
python simulation.py --L 50 --h 8 --lX 2.5 --n_pas 5000000 --reseau hex            # out/L50_h8_lX2.5_b0_g0_hex
python simulation.py --L 50 --h 8 --lX 2.5 --n_pas 5000000 --contact --rho0 0.01   # out/L50_h8_lX2.5_b0_g0_c_r0.01
python simulation.py --L 50 --h 8 --lX 2.5 --n_pas 5000000 --tB 1000 --alpha 0.5   # out/L50_h8_lX2.5_b0_g0_a0.5_tB1000.0
```

Un run repris doit avoir les mêmes paramètres que ceux sauvegardés dans son `params.json` (sauf `n_pas` et `--historique`) : sinon `simulation.py` s'arrête avec un message d'erreur. `lambda_fft` (dans `simulation.py`) mesure la longueur d'onde d'une grille, sur les deux réseaux.

## Lancer sur GenOuest

Depuis la racine du projet :

```bash
bash slurm_pipeline/01_prepare_chunks.sh --h-start 2 --h-end 20 --h-step 2   # runs et chunks
# options du modèle pour tous les runs : --extra "--reseau hex --contact --rho0 0.01"
bash slurm_pipeline/03_submit_all.sh                                         # Job Array (et relance)
python slurm_pipeline/04_merge_results.py                                    # paramètres et état des runs
```

Exemple complet pas à pas : [slurm_pipeline/doc/exemple.md](slurm_pipeline/doc/exemple.md).
