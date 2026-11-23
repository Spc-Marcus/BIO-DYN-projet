# Tests des extensions aux points de la Fig. 2

`05_tests_extensions.sh` lance chaque extension du modèle aux 12 points de la Fig. 2 : seuls h (5, 10, 15) et lX (0.5, 1.5, 2.5, 3.5) varient, comme dans la figure ; tout le reste est celui de l'article (L = 100, 10⁹ pas soit t = 10⁵, b_X = s = 1, graine 0). Chaque série se compare donc directement aux runs de la Fig. 2 (`01_prepare_chunks.sh` sans option).

```bash
bash slurm_pipeline/05_tests_extensions.sh
```

7 variantes × 3 h × 4 lX = 84 runs, en 9 tâches de 10 runs (2 à 3 h chacune), tous avec `--historique`. Une fois toutes les tâches finies, `06_timelines.slurm` écrit `timeline.gif` et la frise `timeline.png` dans le dossier de chaque run. Relancer la même commande reprend les runs interrompus.

## Variantes

| Variante | Options | Pourquoi |
|---|---|---|
| Réseau hexagonal | `--reseau hex` | Retrouve-t-on la même séquence de motifs et λ ≈ 2h ? Les rayures devraient moins suivre les axes du réseau. |
| Départ semé | `--rho0 0.01` | Témoin du contact : même départ (1 % de jaunes, soit 100 cases), sans la règle du contact. |
| Contact | `--contact --rho0 0.01` | Les naissances partent des cases occupées : fronts au départ ; mêmes motifs à la fin attendus. |
| Blessures | `--lB 2.5e-08` | Quelques blessures au hasard par run : la zone se remplit, puis le motif se reforme. |
| Blessures et contact | `--lB 2.5e-08 --contact --rho0 0.01` | C'est après une blessure que le contact change le plus : remplissage depuis le bord au lieu de partout à la fois. |
| Blessure unique | `--tB 50000` | Une blessure au centre, à un instant connu, sur un motif formé depuis longtemps ; la grille d'avant est gardée dans `grille_avant_blessure.npy`. |
| Morts dispersées | `--dX d --dM d` | Témoin des blessures : autant de cases vidées en moyenne, mais une par une et partout (règles (4) et (5)). Seule la forme de la blessure change. |

## Valeurs

**Taux de blessure lB = 2.5·10⁻⁸.** La règle (9) est tirée avec la probabilité p_B = lB / (Σ + lB) à chaque pas, avec Σ = b_X + 2s + l_X = 3 + l_X :

| lX | p_B par pas | Blessures par run (10⁹ pas) | Run sans blessure | Intervalle moyen |
|---|---|---|---|---|
| 0.5 | 7.1·10⁻⁹ | 7.1 | 0.1 % | 14 000 u.t. |
| 1.5 | 5.6·10⁻⁹ | 5.6 | 0.4 % | 18 000 u.t. |
| 2.5 | 4.5·10⁻⁹ | 4.5 | 1.1 % | 22 000 u.t. |
| 3.5 | 3.8·10⁻⁹ | 3.8 | 2.1 % | 26 000 u.t. |

- Extrêmement rare : une case donnée est dans une blessure de l'ordre d'une fois en 10⁶ u.t. (p_B × |D_R| par u.t.), alors qu'elle est tirée une fois par u.t.
- Quelques fois dans une vie : 4 à 7 blessures par run, au plus 2 % de runs sans blessure.
- Séparées : 14 000 à 26 000 u.t. entre deux blessures, bien plus que la reconstruction (ci-dessous). Chaque blessure cicatrise avant la suivante.

**Rayon R = h** (`--alpha` 1, la valeur par défaut). Le diamètre 2h ≈ λ efface une rayure noire et une jaune. Pour h = 5, 10, 15, le disque fait 81, 317 et 709 cases, soit 0.8 %, 3.2 % et 7.1 % de la peau : une large zone reste intacte pour comparer.

**Morts dispersées** : d = lB × |D_R| = 2.025·10⁻⁶, 7.925·10⁻⁶ et 1.7725·10⁻⁵ pour h = 5, 10 et 15, soit autant de cases vidées par u.t. que les blessures.

**Blessure unique à t = 50 000** : le motif est formé en moins de 5 000 u.t. ; il est donc stationnaire depuis longtemps, et il reste 50 000 u.t. pour suivre la suite.

**Contact** : rho0 = 0.01, comme dans le prototype de la spec ; sans cases occupées au départ, une peau vide le resterait.

## Historique

Avec `--historique`, `simulation.py` garde une grille à chaque sauvegarde (tous les 10⁷ pas, soit 1 000 u.t. à L = 100). Pour voir la reconstruction, il en garde aussi une toutes les 5 u.t. pendant 500 u.t., au départ et après chaque blessure (`DT_RECONSTRUCTION` et `DUREE_RECONSTRUCTION` en tête de `simulation.py`).

Mesuré à L = 100 sur un motif formé, avec une blessure R = h au centre (un run par cas, h = 10 et 15) :

| Après la blessure | Sans contact | Avec contact |
|---|---|---|
| Zone remplie à ~90 % | 10 à 20 u.t. | 20 à 40 u.t., par un front venu du bord |
| Autant de noir qu'ailleurs | 100 à 200 u.t. | 100 à 200 u.t. |
| Zone qu'on ne distingue plus du reste | 200 à 500 u.t. (jusqu'à 1 000 pour lX = 3.5) | idem |

Une grille toutes les 5 u.t. montre donc le remplissage en 2 à 8 images et la reconstruction en 40 à 100 images. Un run garde ainsi environ 200 grilles sans blessure et 600 à 900 avec (quelques Mo). Le GIF demande environ 5 Mo de mémoire par grille, d'où les 8 Go de `06_timelines.slurm`.
