# Tests des extensions aux points de la Fig. 2

`05_tests_extensions.sh` lance chaque extension du modèle aux 12 points de la Fig. 2 : seuls h (5, 10, 15) et lX (0.5, 1.5, 2.5, 3.5) varient, comme dans la figure ; tout le reste est celui de l'article (L = 100, n_pas = 10⁹, b_X = s = 1, graine 0). Chaque série se compare donc directement aux runs de la Fig. 2 (`01_prepare_chunks.sh` sans option). Toutes les durées sont en pas (L = 100 : un passage sur toute la grille = 10⁴ pas).

```bash
bash slurm_pipeline/05_tests_extensions.sh
```

5 variantes × 3 h × 4 lX = 60 runs, en 6 tâches de 10 runs (2 à 3 h chacune), tous avec `--historique`. Une fois toutes les tâches finies, `06_timelines.slurm` écrit `timeline.gif` et la frise `timeline.png` dans le dossier de chaque run. Relancer la même commande reprend les runs interrompus.

## Variantes

| Variante | Options | Ce qu'on voit |
|---|---|---|
| Réseau hexagonal | `--reseau hex` | Les motifs de la Fig. 2 sur un réseau à 6 voisins (image en briques) : même séquence, λ ≈ 2h, rayures moins collées aux axes. |
| Contact | `--contact --rho0 0.01` | Au départ, les cellules poussent en fronts autour de 100 graines jaunes, au lieu d'apparaître partout ; à la fin, les mêmes motifs. |
| Blessures | `--lB 2.5e-08` | 4 à 7 trous au hasard pendant le run : chacun se remplit partout à la fois, puis le motif s'y reforme. |
| Blessures et contact | `--lB 2.5e-08 --contact --rho0 0.01` | Les mêmes trous, remplis depuis le bord vers le centre. |
| Blessure unique | `--pas_blessure 500000000` | Un seul trou, au centre, au milieu du run, sur un motif formé depuis longtemps ; la grille d'avant est gardée dans `grille_avant_blessure.npy` pour comparer. |

## Valeurs

**Taille de la blessure.** R = α·h avec α = 1 (`--alpha`, valeur par défaut) : le diamètre 2h vaut environ une période T du motif (T ≈ 2h dans l'article), soit une rayure noire et une jaune. La blessure garde ainsi la même taille par rapport au motif quand h change.

| h | Période T ≈ 2h | Rayon R | Diamètre | Cases vidées | Part de la peau |
|---|---|---|---|---|---|
| 5 | ≈ 10 cases | 5 | 10 cases | 81 | 0.8 % |
| 10 | ≈ 20 cases | 10 | 20 cases | 317 | 3.2 % |
| 15 | ≈ 30 cases | 15 | 30 cases | 709 | 7.1 % |

Même la plus grande reste sous 10 % de la peau : on voit le trou au milieu d'un motif intact. Pour deux périodes, `--alpha 2` (R = 2h, mais 28 % de la peau pour h = 15).

**Taux de blessure lB = 2.5·10⁻⁸.** La règle (9) est tirée avec la probabilité p_B = lB / (Σ + lB) à chaque pas, avec Σ = b_X + 2s + l_X = 3 + l_X :

| lX | p_B par pas | Blessures par run (10⁹ pas) | Run sans blessure | Pas entre deux blessures (moyenne) |
|---|---|---|---|---|
| 0.5 | 7.1·10⁻⁹ | 7.1 | 0.1 % | 1.4·10⁸ |
| 1.5 | 5.6·10⁻⁹ | 5.6 | 0.4 % | 1.8·10⁸ |
| 2.5 | 4.5·10⁻⁹ | 4.5 | 1.1 % | 2.2·10⁸ |
| 3.5 | 3.8·10⁻⁹ | 3.8 | 2.1 % | 2.6·10⁸ |

- Extrêmement rare : une case donnée est dans une blessure de l'ordre d'une fois en 10¹⁰ pas, alors qu'elle est tirée une fois tous les 10⁴ pas.
- Quelques fois dans une vie : 4 à 7 blessures par run, au plus 2 % de runs sans blessure.
- Séparées : 1.4 à 2.6·10⁸ pas entre deux blessures, bien plus que la reconstruction (au plus 10⁷ pas, ci-dessous). Chaque trou se referme avant le suivant.

**Blessure unique au pas 5·10⁸**, le milieu du run : le motif est formé en moins de 5·10⁷ pas ; on blesse donc une peau « adulte », et il reste 5·10⁸ pas pour suivre la suite.

**Contact, rho0 = 0.01** : avec la règle du contact, une peau vide le resterait ; on sème 1 % de jaunes (100 cases), la valeur du prototype de la spec.

## Historique

Avec `--historique`, `simulation.py` garde une grille à chaque sauvegarde (tous les `PAS_ENTRE_SAUVEGARDES` = 10⁷ pas). Pour voir la reconstruction, il en garde aussi une tous les `PAS_ENTRE_GRILLES` = 5·10⁴ pas pendant `PAS_RECONSTRUCTION` = 5·10⁶ pas, au départ et après chaque blessure (constantes en tête de `simulation.py`).

Mesuré à L = 100 sur un motif formé, avec une blessure R = h au centre (un run par cas, h = 10 et 15) :

| Après la blessure | Sans contact | Avec contact |
|---|---|---|
| Trou rempli à ~90 % | 1 à 2·10⁵ pas | 2 à 4·10⁵ pas, par un front venu du bord |
| Autant de noir qu'ailleurs | 1 à 2·10⁶ pas | 1 à 2·10⁶ pas |
| Zone qu'on ne distingue plus du reste | 2 à 5·10⁶ pas (jusqu'à 10⁷ pour lX = 3.5) | idem |

Une grille tous les 5·10⁴ pas montre donc le remplissage en 2 à 8 images et la reconstruction en 40 à 100 images. Un run garde ainsi environ 200 grilles sans blessure et 600 à 900 avec (quelques Mo). Le GIF demande environ 5 Mo de mémoire par grille, d'où les 8 Go de `06_timelines.slurm`.
