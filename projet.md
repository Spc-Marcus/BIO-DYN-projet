# Notes : Bullara & De Decker (2015)

> *Pigment cell movement is not required for generation of Turing patterns in zebrafish skin*, Nat. Commun. 6:6971. DOI [10.1038/ncomms7971](https://doi.org/10.1038/ncomms7971)
> [n] renvoie à la biblio de l'article ; (éq. n) et (Fig. n) à la numérotation de l'article ; **[perso]** = vérifié ou ajouté par moi (sympy/numpy), absent de l'article.

---

## 1. Contexte et question

- Zébrafish : rayures jaunes/noires formées par des **cellules pigmentaires** sous la peau.

| Cellule | Couleur | Notation | Rôle dans le modèle |
|---|---|---|---|
| Xantophore | jaune | X | oui |
| Mélanophore | noir | M | oui |
| Iridophore | argenté | — | seulement comme pré-motif (Fig. 1b–c) |

- **Turing [2]** : un activateur s'auto-amplifie localement, un inhibiteur diffuse plus loin → *activation courte portée + inhibition longue portée* → motif périodique de taille intrinsèque.
  - Cadre RD (éq. 1) : $\partial_t \mathbf{x} = \mathbf{f}(\mathbf{x};\lambda) + \nabla\cdot[\mathbf{D}\nabla\mathbf{x}]$
  - Bifurcation de Turing : une valeur propre réelle $\omega(k)$ du système linéarisé passe par 0 en un $k_T \neq 0$ unique → longueur d'onde $2\pi/k_T$.
- **Pourquoi le poisson pose problème** :
  - motifs faits de **cellules discrètes**, pas de molécules qui diffusent ;
  - cellules **quasi immobiles** [15] ;
  - échelle minuscule : demi-longueur d'onde ≈ 10 cellules (~1 mm), donc une approximation continue discutable.
- **Données de Kondo [13]** :
  - X et M s'éliminent mutuellement **au contact** ;
  - X favorise l'apparition de M **à ~10 cellules**.
  - Aucun transport nécessaire. Nüsslein-Volhard [14, 15] : même schéma avec les iridophores, pas de mouvement à grande échelle.
- **Question** : ces règles seules, avec des cellules immobiles, donnent-elles une vraie bifurcation de Turing ?
- **Pourquoi 2 espèces** : X + M suffisent dans les nageoires [16], dont le motif est continu avec celui du tronc → X–M est le mécanisme cœur.

---

## 2. Modèle discret (section *Model derivation*)

- Peau = **grille** de pas $a$ (taille d'une cellule), bords périodiques.
  - Case $i$ ∈ {X, M, S (vide)}, avec $X_i + M_i + S_i = 1$ (pas de superposition [23]).
  - Pour un informaticien : un **automate cellulaire stochastique**.
- Processus (taux = probabilité par unité de temps) :

| Éq. | Processus | Réaction | Taux |
|---|---|---|---|
| (2) | naissance de X | $`\mathrm{S}_i \to \mathrm{X}_i`$ | $b_X$ |
| (3) | naissance de M | $`\mathrm{S}_i \to \mathrm{M}_i`$ | $b_M$ |
| (4) | mort de X | $`\mathrm{X}_i \to \mathrm{S}_i`$ | $d_X$ |
| (5) | mort de M | $`\mathrm{M}_i \to \mathrm{S}_i`$ | $d_M$ |
| (6) | M voisin tue X | $`\mathrm{M}_{i\pm1} + \mathrm{X}_i \to \mathrm{M}_{i\pm1} + \mathrm{S}_i`$ | $s_M$ |
| (7) | X voisin tue M | $`\mathrm{X}_{i\pm1} + \mathrm{M}_i \to \mathrm{X}_{i\pm1} + \mathrm{S}_i`$ | $s_X$ |
| (8) | X à distance $h$ fait naître M | $`\mathrm{X}_{i\pm h} + \mathrm{S}_i \to \mathrm{X}_{i\pm h} + \mathrm{M}_i`$ | $l_X$ |

- Hypothèses :
  - $s_X, s_M \gg d_X, d_M$ (la mort par compétition est rapide).
  - $h$ ≈ 10 cellules, **non relié a priori** à λ : le fait que $\lambda \approx 2h$ est un résultat, pas une hypothèse.
- **Longue portée** : 3 effets observés.
  - (a) X favorise la *naissance* de M ;
  - (b) X favorise la *survie* de M (Delta-Notch [27]) ;
  - (c) M inhibe la naissance de M.
  - Seul (a) est indispensable pour que $h$ fixe λ → **seul (8) est gardé** (modèle minimal).
- **Cellules immobiles** : c'est l'hypothèse testée. Les M bougent un peu [28, 29], mais rien ne prouve que ce mouvement soit nécessaire.

---

## 3. Intuition du mécanisme

- **(6)–(7) ségrègent** les couleurs : chaque couleur est protégée en groupe → activation courte portée.
  - Sans (8) : domaines qui grossissent sans fin (coarsening), aucune taille préférée.
- **(8)** : un amas X fait naître du M à distance $h$, ce M tue le X à son contact → **X s'inhibe lui-même à distance $h$** → inhibition longue portée sans transport.

```
distance au centre d'un amas X :
0 ────────── ~h ────────── ~2h
[  X stable  ][ M né ici  ][  X  ]
←────────── λ ≈ 2h ──────────→
```

- Nom proposé par les auteurs : **differential growth**. Les cellules ne bougent pas, mais leurs taux de naissance et de mort dépendent du voisinage → redistribution effective des populations.

---

## 4. Simulations Monte Carlo (section *Methods*)

**Algorithme** :
1. grille $N_0 = L^2$, 4 voisins, périodique, initialement tout S ;
2. $P_k = \text{taux}_k / \sum \text{taux}$ ;
3. tirer une case $i$ et un processus $k$ ; appliquer si la configuration le permet :
   - naissance : $i$ doit être S ;
   - (6)/(7) : tirer 1 voisin parmi 4, vérifier les couleurs ;
   - (8) : $i$ doit être S ; tirer un angle $\alpha$, regarder $j = (i_x + h\cos\alpha,\ i_y + h\sin\alpha)$ ; si $j$ = X, alors $i$ devient M ;
4. $t \mathrel{+}= 1/N_0$ à chaque tentative, réussie ou non.

**Paramètres** : $b_M = d_X = d_M = 0$ (processus lents devant les interactions), $s_X = s_M = s$.

| Figure | Paramètres |
|---|---|
| Fig. 1 | $100\times100$, $b_X = s = 1$, $l_X = 2.5$, $h = 16$, $10^9$ pas MC |
| Fig. 2 | $b_X = s = 1$, $h \in \{5,10,15\}$, $l_X \in \{0.5, 1.5, 2.5, 3.5\}$ |

**Résultats** :
- Motifs **stationnaires**.
- **$l_X$ → forme** (Fig. 2) :

| $l_X$ | Motif |
|---|---|
| faible | tout jaune |
| ↑ | taches noires sur fond jaune |
| ↑↑ | rayures |
| ↑↑↑ | taches jaunes sur fond noir |

- **$h$ → taille** : $\lambda \approx 2h$ (cf. expériences). $l_X$ affecte peu λ.
- Robustesse : motifs nets à quelques cellules de période, alors que le bruit détruit souvent les petits motifs en RD stochastique.
- **Pré-motif d'iridophores** (Fig. 1b–c ; bande où (3) et (8) sont interdits) :
  - **oriente** les rayures parallèlement à la bande, même pour une bande de 1 cellule ;
  - λ inchangée ;
  - sans bande (Fig. 1a) : rayures obliques, pour caser un nombre entier de périodes dans la boîte.
  - Cohérent avec [15, 25].

---

## 5. Équations continues

**Équations des moyennes** (1D, issues de l'équation maîtresse) :

$$
\frac{d\langle X_i\rangle}{dt} = b_X\langle S_i\rangle - d_X\langle X_i\rangle - \tfrac12 s_M\big[\langle X_i, M_{i-1}\rangle + \langle X_i, M_{i+1}\rangle\big] \quad (9)
$$

$$
\frac{d\langle M_i\rangle}{dt} = b_M\langle S_i\rangle - d_M\langle M_i\rangle - \tfrac12 s_X\big[\langle M_i, X_{i-1}\rangle + \langle M_i, X_{i+1}\rangle\big] + \tfrac12 l_X\big[\langle S_i, X_{i-h}\rangle + \langle S_i, X_{i+h}\rangle\big] \quad (10)
$$

- $\langle A_i\rangle$ = probabilité que la case $i$ soit A ; $\langle A_i, B_j\rangle$ = probabilité jointe.
- $\tfrac12$ = choix d'un voisin parmi 2 en 1D ($\tfrac14$ en 2D).

**Approximations** :
1. **champ moyen** : $\langle A_i, B_j\rangle \approx \langle A_i\rangle\langle B_j\rangle$, i.e. cases supposées indépendantes ;
2. **continu** : $r = ia$, $x(r) = \langle X_i\rangle$, $m(r) = \langle M_i\rangle$ ;
3. **Taylor** à l'ordre 2 en $a$ (éq. 11–13). La moyenne des voisins = valeur au centre + correction de courbure :

$$
\tfrac12[m(r+a) + m(r-a)] \approx m + \tfrac{a^2}{2}\nabla^2 m, \qquad \tfrac12[x(r+ha) + x(r-ha)] \approx x + \tfrac{h^2a^2}{2}\nabla^2 x
$$

**EDP obtenues** :

$$
\partial_t x = b_X(1-x-m) - (d_X + s_M m)\,x - a^2\frac{s_M x}{2}\nabla^2 m \quad (14)
$$

$$
\partial_t m = (b_M + l_X x)(1-x-m) - (d_M + s_X x)\,m + a^2\Big[\frac{h^2 l_X}{2}(1-x-m) - \frac{s_X m}{2}\Big]\nabla^2 x \quad (15)
$$

- Termes sans $\nabla^2$ = cinétique locale (peau bien mélangée).
- Termes en $\nabla^2$ : **uniquement croisés** ($\nabla^2 m$ dans l'équation de $x$ et inversement), aucune auto-diffusion. Ils ressemblent à de la *cross-diffusion*, mais **ce n'est pas du mouvement** : c'est la trace des interactions non locales.
- Le facteur $h^2$ amplifie l'effet longue portée.

---

## 6. Analyse de stabilité linéaire

**Limite étudiée** : $d_X = d_M = b_M = 0$, $s_X = s_M = s$. Cinétique locale :

$$
f = b_X S - s\,xm, \qquad g = l_X x S - s\,xm, \qquad S = 1-x-m
$$

**États stationnaires** :

| Éq. | $(x, m)$ | Sens | Turing ? |
|---|---|---|---|
| (16) | $(1, 0)$ | tout jaune | non |
| (17) | $\left(\dfrac{b_X}{l_X},\ \dfrac{l_X - b_X}{l_X + s}\right)$ | mélange ; existe si $l_X \ge b_X$ | **oui** |
| (18) | $(0, 1)$ | tout noir | non |

- (16) et (18) sont des artefacts de $d_X = d_M = 0$.
- (17) s'obtient en combinant $b_X S = sxm$ et $l_X S = sm$ → $x = b_X/l_X$.

**Linéarisation autour de (17)** [perso] : $\delta x, \delta m \propto e^{ikr + \omega t}$, $\nabla^2 \to -k^2$, $q = a^2k^2$.

$$
\mathbf{J}(k) = \begin{pmatrix} -\dfrac{l_X(b_X+s)}{l_X+s} & -\dfrac{b_X(l_X+s)}{l_X} \\[2mm] -b_X & -\dfrac{b_X(l_X+s)}{l_X}\end{pmatrix} + q\begin{pmatrix} 0 & \dfrac{s\,x_2}{2} \\[2mm] -\dfrac{s\,m_2(h^2-1)}{2} & 0\end{pmatrix}
$$

(on a utilisé $l_X S_2 = s m_2$)

- **Trace** $= -\dfrac{b_X(l_X+s)^2 + l_X^2(b_X+s)}{l_X(l_X+s)} < 0$, **indépendante de $k$**.
  - → pas d'instabilité oscillante possible (pas de Hopf) ;
  - → l'instabilité passe forcément par $\det < 0$, donc motif **stationnaire** = signature de Turing.
- **Déterminant** : $\det \mathbf{J}(q) = c_0 + c_1 q + c_2 q^2$

$$
c_0 = \frac{b_X s(l_X-b_X)}{l_X}, \qquad c_1 = \frac{b_X s[l_X - h^2(l_X-b_X)]}{2l_X}, \qquad c_2 = \frac{b_X s^2(l_X-b_X)(h^2-1)}{4l_X(l_X+s)}
$$

| Coefficient | Signe | Conséquence |
|---|---|---|
| $c_0$ | > 0 | stable à $k = 0$ (perturbation uniforme) |
| $c_2$ | > 0 si $h > 1$ | stable à grand $k$ (détails fins) |
| $c_1$ | < 0 si $h^2 > l_X/(l_X - b_X)$ | condition nécessaire d'instabilité |

- **Seuil** : $c_1^2 = 4c_0c_2$ ; mode critique $q_T = \sqrt{c_0/c_2}$.

**Résultats de l'article** :

$$
h_T = \sqrt{\frac{2s(l_X-b_X) + l_X(s+l_X) + 2\sqrt{l_X(s+b_X)\,s\,(l_X-b_X)}}{(l_X-b_X)(s+l_X)}} \quad (19)
$$

$$
\lambda = \frac{2\pi}{a k_T} = \sqrt2\,\pi\left(\frac{\sqrt{l_X(s+b_X)\,s\,(l_X-b_X)} - s(l_X-b_X)}{s\,b_X}\right)^{-1/2} \quad (20)
$$

- Motif si $h \ge h_T$.
- **Valeurs** [perso], $b_X = s = 1$ :

| $l_X$ | $h_T$ (19) | $\lambda$ (20) | $\lambda/h_T$ |
|---|---|---|---|
| 1.0193 | 8.00 | 10.5 | 1.31 |
| 1.5 | 2.40 | 5.22 | 2.17 |
| 2.5 | 1.81 | 3.99 | 2.20 |
| 3.5 | 1.61 | 3.42 | 2.13 |

- → $\lambda \approx 2h_T$, sauf quand $l_X \to b_X$ (où $m_2 \to 0$ et $h_T \to \infty$).

**Figures** :
- **Fig. 2** : plan $(l_X, h)$. Bleu = $h_T$, rouge = $\lambda/2$, quasi confondues sauf près de $l_X = 1$. Au-dessus : motifs MC ; en dessous : tout jaune.
- **Fig. 3** : $\omega(k)$ pour $l_X = 1.0193$, $s = a = 1$ (il faut $b_X = 1$, absent de la légende [perso]) :

| $h$ | Courbe | État |
|---|---|---|
| 7 | rouge, $\omega < 0$ partout | stable |
| 8 | bleue, max = 0 en $k_T \approx 0.6$ | seuil |
| 9 | jaune, $\omega > 0$ pour $k \in (0.27, 1.16)$ | motif |

  - $\omega(0) \approx -0.006$ pour les 3 courbes : stable aux perturbations homogènes, ce qui caractérise Turing.
  - Reproduite à l'identique par `dispersion.py`.

**Conclusion** : bifurcation de Turing **sans aucun transport**.

---

## 7. Limites

- **Champ moyen** : il néglige les anticorrélations X–M créées par (6)–(7). L'accord MC ↔ EDP n'est que qualitatif (grille grossière de la Fig. 2).
- **Taylor hors validité** [perso] : il faut $kha \ll 1$, or au seuil $k_T h \approx 0.6 \times 8 \approx 4.8$.
  - Noyau exact en 1D : $\cos(kha)$ au lieu de $1 - (kha)^2/2$ ; en 2D : Bessel $J_0(kha)$.
  - Résultats avec noyau exact (1D, $h$ entier) :

| $l_X$ | $h_T$ Taylor | $h_T$ exact | $\lambda/h$ exact |
|---|---|---|---|
| 1.0193 | 8 | 15 | ≈ 2.5 |
| 1.5 | 2.4 | 4 | ≈ 2.5 |
| 2.5 | 1.8 | 3 | ≈ 2.5 |

  - → l'instabilité existe toujours (conclusion **renforcée**), mais (19) et (20) sont approximatives.
- **1D vs 2D** [perso] : les équations sont dérivées en 1D, les simulations sont en 2D. En 2D, les coefficients passent de $a^2/2$ à $a^2/4$ → $h_T$ inchangé, λ divisée par $\sqrt2$.
- $d_X = d_M = 0$ : irréaliste sur des temps longs.
- Mouvement réel des M ignoré : le modèle montre seulement qu'il n'est **pas nécessaire**.
- Iridophores réduits à un masque statique.

---

## 8. Pistes projet DYN

- **(A) Analyse linéaire** : `dispersion.py` reproduit les Fig. 2–3. Extension : noyau exact en remplaçant les termes en $q$ par $1 - \cos(ka)$ et $1 - \cos(kha)$.
- **(B) Monte Carlo** : grille `0 = S, 1 = X, 2 = M`, algorithme de la section 4.
  - $10^9$ pas est impossible en Python pur → `numba` (`@njit`) ou petite grille ($50\times50$, $h \approx 8$).
- **(C) EDP (14)–(15)** : Laplacien à 5 points, périodique, Euler vs RK2.
  - Les modes à grand $k$ ont des valeurs propres complexes (trace fixe, det ∝ $q^2$) : ce sont des **oscillateurs amortis**.
  - Même piège que `EquaRessort.py` : Euler explicite instable si $\Delta t$ est trop grand ; RK2 plus tolérant.
- **(D) Comparaison** : λ mesurée par FFT 2D sur les motifs MC vs λ Taylor (20) vs λ noyau exact. Balayer $l_X$ pour retrouver la séquence taches → rayures → taches inversées.

---

## Références clés

- [2] Turing, A. *The chemical basis of morphogenesis*. Phil. Trans. R. Soc. B 237, 37–72 (1952).
- [13] Nakamasu, A., Takahashi, G., Kanbe, A. & Kondo, S. *Interactions between zebrafish pigment cells responsible for the generation of Turing patterns*. PNAS 106, 8429–8434 (2009).
- [14] Frohnhöfer, H. G., Krauss, J., Maischein, H. & Nüsslein-Volhard, C. *Iridophores and their interactions with other chromatophores are required for stripe formation in zebrafish*. Development 140, 2997–3007 (2013).
- [15] Mahalwar, P., Walderich, B., Singh, A. P. & Nüsslein-Volhard, C. *Local reorganization of xanthophores fine-tunes and colors the striped pattern of zebrafish*. Science 345, 1362–1364 (2014).
- [16] Patterson, L. B. & Parichy, D. M. *Interactions with iridophores and the tissue environment required for patterning melanophores and xanthophores during zebrafish adult pigment stripe formation*. PLoS Genet. 9, e1003561 (2013).
- [23] Hirata, M. et al. *Pigment cell organization in the hypodermis of zebrafish*. Dev. Dyn. 227, 497–503 (2003).
- [25] Singh, A. P., Schach, U. & Nüsslein-Volhard, C. *Proliferation, dispersal and patterned aggregation of iridophores in the skin prefigure striped colouration of zebrafish*. Nat. Cell Biol. 16, 604–611 (2014).
- [27] Hamada, H. et al. *Involvement of Delta/Notch signaling in zebrafish adult pigment stripe patterning*. Development 141, 318–324 (2014).
- [28] Takahashi, G. & Kondo, S. *Melanophores in the stripes of adult zebrafish do not have the nature to gather, but disperse when they have the space to move*. Pigment Cell Melanoma Res. 21, 677–686 (2008).
- [29] Yamanaka, H. & Kondo, S. *In vitro analysis suggests that difference in cell movement during direct interaction can generate various pigment patterns in vivo*. PNAS 111, 1867–1872 (2014).
