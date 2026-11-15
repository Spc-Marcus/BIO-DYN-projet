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

Le run écrit `grille.npy`, `params.json` et `grille.png` dans `out/L50_h8_lX2.5_b0_g0/`, et une grille tous les `PAS_ENTRE_SAUVEGARDES` pas dans `historique/` avec `--historique`. Relancer la même commande reprend un run interrompu et ne refait pas un run terminé.

## Lancer sur GenOuest

Depuis la racine du projet :

```bash
bash slurm_pipeline/01_prepare_chunks.sh --h-start 2 --h-end 20 --h-step 2   # runs et chunks
bash slurm_pipeline/03_submit_all.sh                                         # Job Array (et relance)
python slurm_pipeline/04_merge_results.py                                    # paramètres et état des runs
```

Exemple complet pas à pas : [slurm_pipeline/doc/exemple.md](slurm_pipeline/doc/exemple.md).
