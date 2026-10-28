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

## Algorithme principal

```text
ALGORITHME SimulationPeau(L, h, bX, bM, dX, dM, s, lX, n_pas, bande)

    // --- Initialisation ---
    grille ← tableau L × L rempli de VIDE              // peau vide au départ
    noir_interdit[ligne] ← VRAI si la ligne est dans la bande d'iridophores

    taux ← [bX, bM, dX, dM, s, s, lX]                   // une entrée par règle (2) à (8)
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
    POUR k DE 1 À 7
        cumul ← cumul + proba[k]
        SI u < cumul ALORS RENVOYER k
```

## Remarques

- **Bords périodiques.** Le `mod L` fait que la grille se referme sur elle-même, comme un tore : sortir à droite fait revenir à gauche, comme dans Pac-Man. On évite ainsi les effets de bord.
- **Temps.** Chaque tentative, réussie ou non, avance le temps de 1/L². Une unité de temps correspond donc en moyenne à une tentative par case.
- **Pourquoi « tirer une règle puis vérifier ».** C'est plus simple que de calculer quelles règles sont possibles à chaque instant. Le prix à payer, ce sont beaucoup de tentatives ratées, d'où le grand nombre de pas nécessaires.
- **Iridophores.** Ils ne sont pas simulés comme des cellules. On les modélise seulement par un masque : sur les lignes de la bande, le noir ne peut pas naître, ni par la règle (3), ni par la règle (8).

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
