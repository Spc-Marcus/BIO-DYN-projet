# BIO-DYN-projet

Reproduction de Bullara & De Decker (2015), *Pigment cell movement is not required for generation of Turing patterns in zebrafish skin*, Nat. Commun. 6:6971.

Notes sur l'article : [projet.md](projet.md). Pseudo-code de la simulation : [pseudocode.md](pseudocode.md).

## Fichiers

Code du projet : un run Monte Carlo (`simulation.py`) et son lancement sur GenOuest (`slurm_pipeline/`), puis la théorie et l'analyse des motifs, qui se lancent en local.

| Fichier | Rôle |
|---|---|
| `simulation.py` | Un run Monte Carlo : les 7 règles du pseudo-code, sauvegarde et reprise. Les taux bX, bM, dX, dM et s sont des constantes en tête du fichier |
| `timeline.py` | GIF et planche « timeline » d'un run lancé avec `--historique` |
| `dispersion.py` | Stabilité linéaire : relation de dispersion, seuil h_T (19), longueur d'onde (20), Fig. 2 et 3 de l'article |
| `edp.py` | Équations de champ moyen (EDP (14)-(15) ou version sans développement de Taylor), Euler ou RK2 |
| `analyse.py` | Longueur d'onde et type de motif des grilles Monte Carlo, comparaison avec la théorie |
| `doc/` | [Théorie et analyse](doc/theorie_et_analyse.md) : méthodes, vérifications, résultats et exemple complet |
| `slurm_pipeline/` | Lancement sur GenOuest en Job Array : scripts, [README](slurm_pipeline/README.md) et documentation avec diagrammes dans [slurm_pipeline/doc/](slurm_pipeline/doc/) |

## Lancer un run en local

Avec numpy et matplotlib installés, depuis la racine du projet :

```bash
python simulation.py --L 50 --h 8 --lX 2.5 --n_pas 5000000 --historique   # environ 5 s
python timeline.py out/L50_h8_lX2.5_b0_g0                                  # GIF et planche
```

Le run écrit `grille.npy`, `params.json` et `grille.png` dans `out/L50_h8_lX2.5_b0_g0/`, et une grille tous les `PAS_ENTRE_SAUVEGARDES` pas dans `historique/` avec `--historique`. Relancer la même commande reprend un run interrompu et ne refait pas un run terminé.

## Théorie et analyse

Toujours en local, depuis la racine du projet (quelques secondes chacun) :

```bash
python dispersion.py              # Fig. 2 et 3, tableau des seuils
python edp.py --h 16 --lX 2.5     # champ moyen : un motif se forme
python analyse.py                 # mesure les runs de out/ et les compare à la théorie
```

Les figures et tableaux sont écrits dans `results/`. Méthodes et résultats : [doc/theorie_et_analyse.md](doc/theorie_et_analyse.md).

## Lancer sur GenOuest

Depuis la racine du projet :

```bash
bash slurm_pipeline/01_prepare_chunks.sh --h-start 2 --h-end 20 --h-step 2   # runs et chunks
bash slurm_pipeline/03_submit_all.sh                                         # Job Array (et relance)
python slurm_pipeline/04_merge_results.py                                    # paramètres et état des runs
python analyse.py                                                            # motifs et longueurs d'onde
```

Exemple complet pas à pas : [slurm_pipeline/doc/exemple.md](slurm_pipeline/doc/exemple.md).
