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

Trois extensions, chacune activée par une option de `simulation.py`. Désactivées (par défaut), on retrouve l'algorithme principal.

### Blessure : règle (9), taux `lB`

La case tirée et toutes les cases à distance ≤ R = alpha·h deviennent vides, quel que soit leur état.

```text
(DI, DJ) ← Disque(R)                     // une fois, au début ; il faut 1 ≤ R < L/2
taux[(9)] ← lB

BLESSURE (9) :
    POUR CHAQUE (di, dj) DE (DI, DJ) : case((i + di) mod L, (j + dj) mod L) ← VIDE

// option tB : une seule blessure au centre, dans la boucle, avant le tirage du pas
SI pas = arrondi(tB · L²) ALORS
    sauvegarder la grille, puis vider le disque centré en (L div 2, L div 2)
```

### Naissance au contact : règles (2′), (3′) et (8′)

Avec `contact`, une naissance demande en plus qu'un voisin direct de la case soit occupé : `ContactOK(i, j)` s'ajoute aux conditions de (2), (3) et (8), avant `CaseADistance` (on ne tire pas d'angle pour rien). Une peau vide le resterait : on part d'une fraction `rho0` de cases occupées.

```text
FONCTION ContactOK(i, j)
    RENVOYER VRAI si un des voisins directs de (i, j) n'est pas VIDE    // ne tire rien

FONCTION InitialiserGrille(L, rho0, fX0)
    grille ← VIDE partout
    SI rho0 > 0 ALORS
        POUR max(1, arrondi(rho0 · L²)) cases tirées sans remise
            case ← JAUNE avec la probabilité fX0, NOIR sinon
    RENVOYER grille
```

### Réseau hexagonal

Même matrice L × L, en coordonnées axiales : la case `(i, j)` est au point `x = j + i/2`, `y = (√3/2)·i`. Les 6 voisins, tous à distance 1, sont les 4 du carré plus `(1, −1)` et `(−1, 1)`, et d²(di, dj) = di² + di·dj + dj².

```text
  dans la matrice (di, dj)                  dans le plan (ligne i décalée de i/2)

  (-1,-1)  (-1, 0)  (-1,+1)                     (-1, 0)   (-1,+1)
  ( 0,-1)  ( i, j)  ( 0,+1)                ( 0,-1)    ( i, j)    ( 0,+1)
  (+1,-1)  (+1, 0)  (+1,+1)                     (+1,-1)   (+1, 0)
```

- `VoisinAuHasard` tire parmi les 6 voisins, `ContactOK` les lit tous, `Disque` utilise la distance hexagonale.
- `CaseADistance` tire α comme en carré, puis prend la case la plus proche du point (h cos α, h sin α) par l'arrondi hexagonal (coordonnées « cube »).
- Image « briques » L × 2L (2 pixels par case, ligne i décalée de i/2 case), affichée avec un aspect √3.
- Longueur d'onde par FFT : λ = (√3/2)·L / √(m_i² − m_i·m_j + m_j²) au lieu de L / √(m_i² + m_j²).

## Balayage de paramètres (cluster)

```text
// 01_prepare_chunks.sh : une ligne de paramètres par simulation (défauts = Fig. 2)
POUR h DANS {5, 10, 15}
    POUR lX DANS {0.5, 1.5, 2.5, 3.5}
        écrire (L = 100, h, lX, bande = 0, graine = 0, n_pas = 10⁹) dans un paquet

// 02_run_worker.slurm : chaque tâche du job array traite un paquet de 10 lignes
POUR chaque ligne du paquet attribué à la tâche
    SI le fichier résultat existe déjà ALORS passer à la suivante   // reprise après plantage
    SimulationPeau(paramètres de la ligne) → sauvegarder la grille (.npy)
```
