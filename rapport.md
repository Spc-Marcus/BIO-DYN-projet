<!--
SQUELETTE DU RAPPORT — 5 pages maximum, page de titre non comptée.

Mode d'emploi
- Les passages « > À rédiger » et « > Observations » sont des trous à remplacer par ton texte.
- Chaque image pointe vers figures/<nom> : copie dans figures/ le fichier indiqué en commentaire
  sous l'image (dossiers de run de out/, voir slurm_pipeline/doc/tests_extensions.md).
- Les GIF bougent dans un lecteur Markdown (GitHub, VS Code), pas dans un PDF : en PDF, garde
  la frise (timeline.png) et mets le GIF en lien. Le dossier figures/ doit être envoyé avec le rapport.
- Budget : environ 1 500 mots, 6 figures et 2 tableaux. Longueur visée en tête de chaque partie.
- Les figures 1, 2, 3 et 6 sont à assembler à la main (images côte à côte) : aucun script ne le fait.
- Pour vérifier le nombre de pages : pandoc rapport.md -o rapport.pdf (les GIF doivent alors
  être remplacés par des PNG).
-->

# Motifs de Turing sans mouvement cellulaire : reproduction de Bullara & De Decker (2015) et extensions

**Prénom NOM** — BS1 / DYN — rendu le 29 novembre 2026

Code : `simulation.py` (un run Monte Carlo) et `timeline.py` (GIF et frises).

---

## 1. Contexte et question

<!-- ≈ ½ page -->

> À rédiger (≈ 250 mots) :
> - les rayures du poisson zèbre : cellules pigmentaires jaunes (xantophores, X), noires (mélanophores, M), iridophores ;
> - le mécanisme de Turing : activation à courte portée, inhibition à longue portée, longueur d'onde propre ;
> - les observations de Kondo : X et M s'éliminent au contact, X favorise M à ~10 cellules ; les cellules ne bougent presque pas ;
> - la question de l'article : ces seules règles, sans mouvement des cellules, donnent-elles un vrai motif de Turing ?
> - notre objectif : refaire les Fig. 1 et 2, puis tester trois extensions (blessure, naissance au contact, réseau hexagonal).

## 2. Modèle et méthodes

<!-- ≈ 1 page -->

### 2.1 Modèle de l'article

> À rédiger (≈ 100 mots) : la peau est une grille L × L périodique ; chaque case est vide (S), jaune (X) ou noire (M) ; des règles locales et une règle à distance h.

| Éq. | Règle | Réaction | Taux | Option |
|---|---|---|---|---|
| (2) | naissance de X | S → X | b_X = 1 | `--bX` |
| (3) | naissance de M | S → M | b_M = 0 | `--bM` |
| (4), (5) | mort de X, de M | X → S, M → S | d_X = d_M = 0 | `--dX`, `--dM` |
| (6) | M voisin tue X | M + X → M + S | s = 1 | `--s` |
| (7) | X voisin tue M | X + M → X + S | s = 1 | `--s` |
| (8) | X à distance h fait naître M | X + S → X + M | l_X | `--lX` |

### 2.2 Simulation Monte Carlo

> À rédiger (≈ 200 mots) :
> - un pas = tirer une case, tirer une règle avec P_k = taux_k / somme des taux, l'appliquer si la configuration le permet ;
> - les règles sont dans un dictionnaire construit à partir des options ; seules celles de taux non nul sont tirées ;
> - paramètres de l'article : L = 100, n_pas = 10⁹, peau vide au départ ; les iridophores sont une bande où (3) et (8) sont interdites ;
> - runs lancés sur GenOuest (Job Array Slurm), reprise exacte après une coupure.

### 2.3 Extensions testées

> À rédiger (≈ 150 mots) : pour chaque extension, l'idée biologique et la règle ajoutée.

| Extension | Règle | Valeurs testées |
|---|---|---|
| Blessure | (9) : la case tirée et toutes celles à distance ≤ R = α·h deviennent vides | `--lB 2.5e-08` (4 à 7 blessures par run), R = h ; `--pas_blessure 500000000` (une blessure au centre) |
| Naissance au contact | (2′), (8′) : une naissance demande un voisin occupé | `--contact --rho0 0.01` (départ avec 1 % de jaunes) |
| Réseau hexagonal | 6 voisins à distance 1 au lieu de 4 | `--reseau hex` |

Chaque extension est testée aux 12 points de la Fig. 2 (h = 5, 10, 15 et l_X = 0.5, 1.5, 2.5, 3.5), le reste comme dans l'article.

## 3. Reproduction de l'article

<!-- ≈ 1 page ¼ -->

### 3.1 Fig. 1 : rayures et bande d'iridophores

![Fig. 1 — h = 16, l_X = 2.5, 10⁹ pas : (a) sans iridophores, (b) bande de 1 ligne, (c) bande de 10 lignes](figures/fig1.png)
<!-- les grille.png de out/L100_h16_lX2.5_b0_g0, out/L100_h16_lX2.5_b1_g0, out/L100_h16_lX2.5_b10_g0, côte à côte -->

> Observations :
> - les rayures sont-elles stationnaires ? obliques sans bande, parallèles à la bande avec ?
> - longueur d'onde mesurée (`lambda_fft`) : … cases, à comparer à 2h = 32 ;
> - différences avec la Fig. 1 de l'article ?

### 3.2 Fig. 2 : forme du motif selon l_X et h

![Fig. 2 — grille finale des 12 runs dans le plan (l_X, h)](figures/fig2.png)
<!-- planche des grille.png de out/L100_h{5,10,15}_lX{0.5,1.5,2.5,3.5}_b0_g0 (l_X en abscisse, h en ordonnée) -->

| h \ l_X | 0.5 | 1.5 | 2.5 | 3.5 |
|---|---|---|---|---|
| 15 | motif : … ; λ = … | … | … | … |
| 10 | … | … | … | … |
| 5 | … | … | … | … |

> Observations :
> - la séquence tout jaune → taches noires → rayures → taches jaunes quand l_X augmente se retrouve-t-elle ?
> - λ ≈ 2h ? l_X change-t-il λ ?
> - accord avec les courbes h_T (éq. 19) et λ/2 (éq. 20) de l'article ?

## 4. Extensions

<!-- ≈ 1 page ½ -->

### 4.1 Réseau hexagonal

![Fig. 3 — même run sur réseau carré (gauche) et hexagonal (droite), h = 10](figures/hex.png)
<!-- grille.png de out/L100_h10_lX2.5_b0_g0 et de out/L100_h10_lX2.5_b0_g0_hex (et l_X = 1.5 pour les taches) -->

> Observations :
> - même séquence de motifs qu'en carré ? λ = … (hex) contre … (carré) ;
> - les rayures et les taches suivent-elles moins les axes du réseau ?

### 4.2 Naissance au contact

![Fig. 4 — début du run, sans contact (haut) et avec contact (bas)](figures/contact.png)
<!-- premières images des timeline.gif de out/L100_h10_lX2.5_b0_g0_B2.5e-08 (sans contact) et out/L100_h10_lX2.5_b0_g0_c_r0.01 (avec) :
     les 100 premières grilles couvrent les 5·10⁶ premiers pas -->
<!-- GIF : out/L100_h10_lX2.5_b0_g0_c_r0.01/timeline.gif -->

> Observations :
> - fronts autour des graines au lieu d'apparitions partout ?
> - nombre de pas pour remplir la peau : … avec contact, … sans ;
> - motif final identique à la Fig. 2 ?

### 4.3 Blessures et reconstruction

![Fig. 5 — une blessure qui se referme, sans contact (haut) et avec contact (bas)](figures/blessures.png)
<!-- images des timeline.gif juste avant et après une blessure (une grille tous les 5·10⁴ pas pendant 5·10⁶ pas) :
     out/L100_h10_lX2.5_b0_g0_B2.5e-08 et out/L100_h10_lX2.5_b0_g0_c_r0.01_B2.5e-08 -->

![GIF — blessures au hasard pendant un run](figures/blessures.gif)
<!-- out/L100_h10_lX2.5_b0_g0_B2.5e-08/timeline.gif (en PDF : lien vers le GIF) -->

![Fig. 6 — blessure unique au centre : juste avant, puis pendant la reconstruction](figures/blessure_unique.png)
<!-- frise timeline.png de out/L100_h10_lX2.5_b0_g0_pB500000000 ; grille d'avant : grille_avant_blessure.npy -->

> Observations :
> - nombre de blessures vues par run : … (attendu 4 à 7) ;
> - pas pour remplir le trou : … sans contact, … avec (front venu du bord ?) ;
> - le motif se reforme-t-il comme avant (comparer à grille_avant_blessure.npy) ou autrement (orientation des rayures) ?
> - comparaison avec Yamaguchi 2007 (cellules réapparues dispersées dans la zone) : quelle règle colle le mieux, (2) ou (2′) ?

## 5. Discussion et conclusion

<!-- ≈ ½ page -->

> À rédiger (≈ 250 mots) :
> - ce que montrent la reproduction et les extensions ;
> - limites : champ moyen qualitatif, une case = une cellule sans forme ni taille ni mouvement, une seule couche, iridophores fixes, nombre entier de périodes imposé par L, une seule graine par point ;
> - suites possibles.

## Références

- Bullara, D. & De Decker, Y. (2015). Pigment cell movement is not required for generation of Turing patterns in zebrafish skin. *Nat. Commun.* 6:6971. doi:10.1038/ncomms7971
- Nakamasu, A. et al. (2009). Interactions between zebrafish pigment cells responsible for the generation of Turing patterns. *PNAS* 106:8429. PMC2689028
- Yamaguchi, M., Yoshimoto, E. & Kondo, S. (2007). Pattern regulation in the stripe of zebrafish suggests an underlying dynamic and autonomous mechanism. *PNAS* 104:4790. PMC1820822
- Mahalwar, P. et al. (2016). Heterotypic interactions regulate cell shape and density during color pattern formation in zebrafish. *Biol. Open* 5:1680. doi:10.1242/bio.022251
