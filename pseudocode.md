# Pseudo-code : simulation Monte Carlo (Bullara & De Decker 2015)

## Données

| Nom | Sens | Valeur Fig. 1 |
|---|---|---|
| `L` | taille de la grille (L × L cases) | 100 |
| `h` | portée de l'interaction à distance (en cases) | 16 |
| `bX`, `bM` | taux de naissance spontanée du jaune, du noir | 1, 0 |
| `dX`, `dM` | taux de mort naturelle du jaune, du noir | 0, 0 |
| `s` | taux d'élimination entre voisins de couleurs différentes | 1 |
| `lX` | taux de naissance du noir grâce à un jaune à distance `h` | 2.5 |
| `n_pas` | nombre de tentatives | 10⁹ |
| `bande` | largeur de la bande d'iridophores (0 = aucune) | 0, 1 ou 10 |

Les paramètres des extensions (blessure, naissance au contact, réseau hexagonal) sont dans la section [Extensions](#extensions).

## Algorithme principal

```text
ALGORITHME SimulationPeau(L, h, bX, bM, dX, dM, s, lX, n_pas, bande)

    // --- Initialisation ---
    grille ← tableau L × L rempli de VIDE              // peau vide au départ
    noir_interdit[ligne] ← VRAI si la ligne est dans la bande d'iridophores

    taux ← { (2): bX, (3): bM, (4): dX, (5): dM, (6): s, (7): s, (8): lX }
                                                        // dictionnaire règle → taux, dans l'ordre (2) à (8)
    retirer de taux les règles de taux nul              // seules les règles actives sont tirées
    proba[k] ← taux[k] / somme(taux)                    // probabilité de choisir la règle k

    // --- Boucle Monte Carlo ---
    RÉPÉTER n_pas FOIS
        (i, j) ← une case tirée au hasard
        règle  ← une règle tirée au hasard selon proba

        SELON règle
            NAISSANCE_JAUNE (2) :
                SI case(i,j) = VIDE ALORS case(i,j) ← JAUNE

            NAISSANCE_NOIR (3) :
                SI case(i,j) = VIDE ET NON noir_interdit[i] ALORS case(i,j) ← NOIR

            MORT_JAUNE (4) :
                SI case(i,j) = JAUNE ALORS case(i,j) ← VIDE

            MORT_NOIR (5) :
                SI case(i,j) = NOIR ALORS case(i,j) ← VIDE

            NOIR_TUE_JAUNE (6) :
                SI case(i,j) = JAUNE ET VoisinAuHasard(i,j) = NOIR
                    ALORS case(i,j) ← VIDE

            JAUNE_TUE_NOIR (7) :
                SI case(i,j) = NOIR ET VoisinAuHasard(i,j) = JAUNE
                    ALORS case(i,j) ← VIDE

            NOIR_NAIT_A_DISTANCE (8) :
                SI case(i,j) = VIDE ET NON noir_interdit[i]
                   ET CaseADistance(i,j,h) = JAUNE
                    ALORS case(i,j) ← NOIR
        FIN SELON
        // si la condition n'est pas remplie, on ne fait rien : la tentative est « ratée »
    FIN RÉPÉTER

    RENVOYER grille
```

## Fonctions auxiliaires

```text
FONCTION VoisinAuHasard(i, j)
    // un des 4 voisins directs : haut, bas, gauche, droite
    (di, dj) ← tiré au hasard parmi (1,0), (-1,0), (0,1), (0,-1)
    RENVOYER case((i + di) mod L, (j + dj) mod L)       // mod L : bords périodiques

FONCTION CaseADistance(i, j, h)
    // un point sur le cercle de rayon h autour de (i, j)
    α ← angle tiré au hasard dans [0, 2π[
    k ← arrondi(i + h·cos α) mod L
    l ← arrondi(j + h·sin α) mod L
    RENVOYER case(k, l)

FONCTION TirerRègle(proba)
    // tirage pondéré : on découpe [0, 1] en segments de longueur proba[k]
    u ← nombre aléatoire dans [0, 1[
    cumul ← 0
    POUR CHAQUE règle k DE proba, dans l'ordre
        cumul ← cumul + proba[k]
        SI u < cumul ALORS RENVOYER k
    RENVOYER la dernière règle                         // u très proche de 1 : erreurs d'arrondi
```

## Remarques

- **Bords périodiques.** Le `mod L` fait que la grille se referme sur elle-même, comme un tore : sortir à droite fait revenir à gauche, comme dans Pac-Man. On évite ainsi les effets de bord.
- **Temps.** Chaque tentative, réussie ou non, avance le temps de 1/L². Une unité de temps correspond donc en moyenne à une tentative par case.
- **Pourquoi « tirer une règle puis vérifier ».** C'est plus simple que de calculer quelles règles sont possibles à chaque instant. Le prix à payer, ce sont beaucoup de tentatives ratées, d'où le grand nombre de pas nécessaires.
- **Iridophores.** Ils ne sont pas simulés comme des cellules. On les modélise seulement par un masque : sur les lignes de la bande, le noir ne peut pas naître, ni par la règle (3), ni par la règle (8).
- **Règles actives.** Une règle de taux nul est retirée du dictionnaire : placée en dernier, elle serait sinon renvoyée par le repli de `TirerRègle`. Avec les taux de l'article il reste (2), (6), (7) et (8).

## Extensions

Trois extensions du modèle, chacune activée par une option de `simulation.py`. Avec les valeurs par défaut, on retrouve exactement l'algorithme principal.

| Nom | Sens | Défaut |
|---|---|---|
| `lB` | taux de la blessure, règle (9) | 0 (pas de blessure) |
| `alpha` | rayon de la blessure en périodes : `R = alpha · h` | 1 |
| `tB` | temps d'une blessure unique au centre de la grille | aucun |
| `contact` | naissances seulement au contact d'une case occupée, règles (2′), (3′) et (8′) | non |
| `rho0`, `fX0` | fraction de cases occupées au départ, part de jaunes parmi elles | 0, 1 |
| `reseau` | réseau carré (z = 4 voisins) ou hexagonal (z = 6) | carré |

### Blessure : règle (9)

La case tirée `i` et toutes les cases `j` à distance `d(i, j) ≤ R` deviennent vides, quel que soit leur état (A ∈ {S, X, M}). `d` est la distance euclidienne en cases, au plus court sur le tore, comme pour la règle (8).

| Éq. | Processus | Réaction | Taux |
|---|---|---|---|
| (9) | blessure | A_j → S_j pour toute case j à distance ≤ R de i | `lB` |

```text
// --- Initialisation, en plus ---
R ← alpha · h                                       // il faut 1 ≤ R < L/2
(DI, DJ) ← Disque(R)                                // décalages de toutes les cases à distance ≤ R
taux[(9)] ← lB                                      // retirée si lB = 0, comme les autres

// --- Dans le SELON ---
BLESSURE (9) :
    POUR CHAQUE (di, dj) DE (DI, DJ)
        case((i + di) mod L, (j + dj) mod L) ← VIDE

// --- Blessure unique (option tB), dans la boucle, avant le tirage du pas ---
SI pas = arrondi(tB · L²) ALORS
    sauvegarder la grille                            // grille_avant_blessure.npy
    vider le disque centré en (L div 2, L div 2), comme la règle (9)
```

La blessure unique n'est pas une règle tirée. Comme la boucle reprend à `pas_faits`, elle n'est faite qu'une fois, même si le run est repris ou prolongé.

### Naissance au contact : règles (2′), (3′) et (8′)

Une case vide ne reçoit une naissance que si au moins un de ses z voisins directs est occupé :

```text
NAISSANCE_JAUNE (2') :
    SI case(i,j) = VIDE ET (NON contact OU ContactOK(i,j)) ALORS case(i,j) ← JAUNE

NAISSANCE_NOIR (3') :
    SI case(i,j) = VIDE ET NON noir_interdit[i] ET (NON contact OU ContactOK(i,j))
        ALORS case(i,j) ← NOIR

NOIR_NAIT_A_DISTANCE (8') :
    SI case(i,j) = VIDE ET NON noir_interdit[i] ET (NON contact OU ContactOK(i,j))
       ET CaseADistance(i,j,h) = JAUNE
        ALORS case(i,j) ← NOIR
```

`ContactOK` ne tire rien ; le tester avant `CaseADistance` évite un tirage d'angle inutile. Avec le contact, une peau vide le resterait : on part d'une fraction `rho0` de cases occupées (il faut `rho0 > 0`).

```text
FONCTION InitialiserGrille(L, rho0, fX0)
    grille ← VIDE partout
    SI rho0 > 0 ALORS
        n0 ← max(1, arrondi(rho0 · L²))
        POUR CHAQUE case c parmi n0 cases tirées au hasard sans remise    // exactement n0
            case(c) ← JAUNE avec la probabilité fX0, NOIR sinon
    RENVOYER grille

FONCTION ContactOK(i, j)
    POUR CHAQUE (di, dj) DES voisins directs du réseau
        SI case((i + di) mod L, (j + dj) mod L) ≠ VIDE ALORS RENVOYER VRAI
    RENVOYER FAUX
```

Avec `rho0 = 0`, rien n'est tiré au départ : la suite des tirages est celle de l'algorithme principal.

### Réseau hexagonal

On garde la matrice L × L, en coordonnées axiales : la case `(i, j)` est au point `x = j + i/2`, `y = (√3/2)·i` du plan. On ajoute une seule diagonale aux 4 voisins du réseau carré : les 6 voisins sont alors tous à distance 1.

```text
  dans la matrice (di, dj)                  dans le plan (ligne i décalée de i/2)

  (-1,-1)  (-1, 0)  (-1,+1)                     (-1, 0)   (-1,+1)
  ( 0,-1)  ( i, j)  ( 0,+1)                ( 0,-1)    ( i, j)    ( 0,+1)
  (+1,-1)  (+1, 0)  (+1,+1)                     (+1,-1)   (+1, 0)

  voisins : les 4 du réseau carré + (-1,+1) et (+1,-1) ; (-1,-1) et (+1,+1) n'en sont pas
```

Seules les fonctions qui dépendent de la géométrie changent :

```text
FONCTION VoisinAuHasard(i, j)
    (di, dj) ← tiré au hasard parmi (1,0), (-1,0), (0,1), (0,-1), (1,-1), (-1,1)
    RENVOYER case((i + di) mod L, (j + dj) mod L)

FONCTION CaseADistance(i, j, h)
    α ← angle tiré au hasard dans [0, 2π[
    (di, dj) ← ArrondiHex(h·cos α, h·sin α)          // case la plus proche du point
    RENVOYER case((i + di) mod L, (j + dj) mod L)

FONCTION ArrondiHex(dx, dy)
    // coordonnées « cube » fractionnaires, de somme nulle
    r ← 2·dy / √3 ;  q ← dx − r/2 ;  s ← −q − r
    arrondir q, r et s à l'entier le plus proche
    recalculer celle des trois qui a le plus grand écart d'arrondi, pour garder q + r + s = 0
    RENVOYER (di, dj) = (r, q)

FONCTION Distance2(di, dj)
    SI réseau carré ALORS RENVOYER di² + dj²
    SINON RENVOYER di² + di·dj + dj²                 // réseau hexagonal

FONCTION Disque(R)
    RENVOYER tous les (di, dj) avec Distance2(di, dj) ≤ R²
```

- **Bords.** Ils restent `mod L` sur `i` et sur `j` : le domaine est un losange périodique à 60°, donc un tore, d'aire (√3/2)·L².
- **Iridophores.** Les lignes `i` restent horizontales dans le plan : la bande reste une bande.
- **Image.** Chaque case occupe 2 pixels d'une ligne, et la ligne `i` est décalée de `i/2` case (image « briques » L × 2L, affichée avec un rapport d'aspect √3) : chaque brique touche exactement ses 6 voisines.
- **Longueur d'onde.** Avec les indices signés `(m_i, m_j)` du pic de la FFT de la matrice, λ = L / √(m_i² + m_j²) sur le réseau carré, et λ = (√3/2)·L / √(m_i² − m_i·m_j + m_j²) sur le réseau hexagonal.

## Balayage de paramètres (cluster)

```text
// gen_params.py : une ligne de paramètres par simulation
POUR lX DANS 30 valeurs entre 1.1 et 4.0
    POUR h DE 2 À 20
        POUR graine DE 0 À 2
            écrire (L=200, h, lX, n_pas = 10⁵ × L², graine) dans params.txt

// sweep.sbatch : chaque tâche du job array traite un paquet de 10 lignes
POUR chaque ligne du paquet attribué à la tâche
    SI le fichier résultat existe déjà ALORS passer à la suivante   // reprise après plantage
    SimulationPeau(paramètres de la ligne) → sauvegarder la grille (.npy)
```
