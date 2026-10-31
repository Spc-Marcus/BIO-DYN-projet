"""
Un run Monte Carlo du modèle de Bullara & De Decker (2015), lancé sur
GenOuest par slurm_pipeline/ (algorithme dans pseudocode.md).

La peau est une grille L x L à bords périodiques. Chaque case vaut VIDE,
JAUNE (xantophore) ou NOIR (mélanophore). À chaque pas on tire une case et
une règle selon les taux, puis on applique la règle si la configuration le
permet. Comme dans l'article : bM = dX = dM = 0 et sX = sM = s, il ne reste
donc que les règles (2), (6), (7) et (8).

Iridophores : ils ne sont pas simulés comme des cellules, seulement par un
masque de cases où le noir ne peut pas naître (Fig. 1b-c). Aujourd'hui une
bande horizontale de `bande` lignes (0 = pas d'iridophores) ; une autre
forme se branche en ne modifiant que masque_iridophores.

Le run est sauvegardé régulièrement dans out/L{L}_h{h}_lX{lX}_b{bande}_g{graine}.
Relancer la même commande récupère ce qui y est sauvegardé : un run coupé
par la limite de temps Slurm reprend là où il s'est arrêté, un run terminé
n'est pas refait.

Utilisation :
    python simulation.py --h 16 --lX 2.5
"""

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")   # pas d'écran sur les nœuds du cluster
import matplotlib.pyplot as plt

# 1. États d'une case
VIDE = 0
JAUNE = 1
NOIR = 2

# 2. Règles restantes, dans l'ordre des taux [BX, S, S, lX]
NAISSANCE_JAUNE = 0        # (2) S_i -> X_i
NOIR_TUE_JAUNE = 1         # (6) M_{i±1} + X_i -> M_{i±1} + S_i
JAUNE_TUE_NOIR = 2         # (7) X_{i±1} + M_i -> X_{i±1} + S_i
NOIR_NAIT_A_DISTANCE = 3   # (8) X_{i±h} + S_i -> X_{i±h} + M_i

# 3. Taux fixés comme dans l'article
BX = 1.0   # naissance spontanée du jaune
S = 1.0    # élimination entre voisins de couleurs différentes

# 4. Nombre de pas entre deux sauvegardes (pour pouvoir relancer)
PAS_ENTRE_SAUVEGARDES = 10**7


def lire_arguments():
    """
    Lit les paramètres du run sur la ligne de commande.

    Options :
      --L       taille de la grille, entier (défaut 100) ;
      --h       portée de l'interaction à distance, en cases, entier
                (défaut 16) ;
      --lX      taux de naissance du noir grâce à un jaune à distance h
                (défaut 2.5) ;
      --bande   largeur de la bande d'iridophores en lignes, entier
                (défaut 0 = pas d'iridophores) ;
      --n_pas   nombre total de tentatives, entier (défaut 1000000000) ;
      --graine  graine du générateur aléatoire, entier (défaut 0).

    Renvoie
    -------
    args : argparse.Namespace
    """
    raise NotImplementedError


def preparer_run(args):
    """
    Prépare le run : nouveau, ou repris depuis son dossier s'il existe déjà.

    Le dossier du run est out/L{L}_h{h}_lX{lX}_b{bande}_g{graine} ; il ne
    dépend pas de n_pas.
      - Pas de params.json dans ce dossier : nouveau run, grille vide,
        générateur créé à partir de la graine, 0 pas effectué.
      - Sinon : reprise, la grille, pas_faits et l'état du générateur sont
        récupérés avec charger_run. Le n_pas de la ligne de commande
        remplace celui sauvegardé, ce qui permet de prolonger un run.

    Paramètres
    ----------
    args : argparse.Namespace
        Sortie de lire_arguments.

    Renvoie
    -------
    params : dict
        Clés "L", "h", "lX", "bande", "n_pas", "graine", "sortie".
    grille : np.ndarray d'entiers, forme (L, L)
        Grille de départ.
    pas_faits : int
        Nombre de pas déjà effectués.
    rng : np.random.Generator
        Générateur aléatoire.
    """
    raise NotImplementedError


def initialiser_grille(L):
    """
    Crée la grille de départ : peau entièrement vide.

    Paramètres
    ----------
    L : int
        Taille de la grille (L x L cases).

    Renvoie
    -------
    grille : np.ndarray d'entiers, forme (L, L)
        Tableau rempli de VIDE.
    """
    raise NotImplementedError


def masque_iridophores(L, bande):
    """
    Construit le masque des cases couvertes par les iridophores.

    Sur ces cases le noir ne peut pas naître : la règle (8) y est bloquée,
    comme la règle (3) si bM était non nul (Fig. 1b-c). Pour l'instant le
    masque est une bande horizontale de `bande` lignes au milieu de la
    grille. Pour une autre forme d'iridophores, il suffit de modifier cette
    fonction : le reste du code ne lit que le masque.

    Paramètres
    ----------
    L : int
        Taille de la grille.
    bande : int
        Largeur de la bande en nombre de lignes (0 = pas d'iridophores).

    Renvoie
    -------
    noir_interdit : np.ndarray de booléens, forme (L, L)
        True sur les cases où le noir ne peut pas naître.
    """
    raise NotImplementedError


def tirer_regle(proba, rng):
    """
    Tire une règle au hasard selon les probabilités données.

    Tirage pondéré : on découpe [0, 1[ en segments de longueur proba[k] et
    on renvoie le segment dans lequel tombe un nombre uniforme u.

    Paramètres
    ----------
    proba : np.ndarray, forme (4,)
        Probabilités des règles, [BX, S, S, lX] / somme.
    rng : np.random.Generator
        Générateur aléatoire.

    Renvoie
    -------
    regle : int
        Indice de la règle tirée, entre 0 et 3.
    """
    raise NotImplementedError


def voisin_au_hasard(grille, i, j, rng):
    """
    Renvoie l'état d'un des 4 voisins directs de la case (i, j).

    Le voisin est tiré parmi haut, bas, gauche, droite ; les bords sont
    périodiques (mod L). Utilisé par les règles (6) et (7).

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État courant de la peau.
    i, j : int
        Coordonnées de la case.
    rng : np.random.Generator
        Générateur aléatoire.

    Renvoie
    -------
    etat : int
        VIDE, JAUNE ou NOIR.
    """
    raise NotImplementedError


def case_a_distance(grille, i, j, h, rng):
    """
    Renvoie l'état d'une case située à distance h de la case (i, j).

    On tire un angle alpha uniforme dans [0, 2π[ et on regarde la case
    (arrondi(i + h cos alpha), arrondi(j + h sin alpha)), modulo L.
    Utilisé par la règle (8).

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État courant de la peau.
    i, j : int
        Coordonnées de la case.
    h : int
        Portée de l'interaction, en cases.
    rng : np.random.Generator
        Générateur aléatoire.

    Renvoie
    -------
    etat : int
        VIDE, JAUNE ou NOIR.
    """
    raise NotImplementedError


def appliquer_regle(grille, i, j, regle, h, noir_interdit, rng):
    """
    Applique une règle à la case (i, j) si la configuration le permet.

      - (2) naissance jaune : la case doit être VIDE ;
      - (6) noir tue jaune : case JAUNE et voisin tiré NOIR ;
      - (7) jaune tue noir : case NOIR et voisin tiré JAUNE ;
      - (8) noir naît à distance : case VIDE, hors iridophores
        (noir_interdit[i, j] faux) et case à distance h JAUNE.

    Sinon la tentative est ratée et la grille n'est pas modifiée.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État courant de la peau, modifié sur place.
    i, j : int
        Coordonnées de la case tirée.
    regle : int
        Indice de la règle, entre 0 et 3.
    h : int
        Portée de l'interaction à distance.
    noir_interdit : np.ndarray de booléens, forme (L, L)
        Masque des iridophores (voir masque_iridophores).
    rng : np.random.Generator
        Générateur aléatoire.

    Renvoie
    -------
    reussie : bool
        True si la grille a été modifiée.
    """
    raise NotImplementedError


def simulation_peau(grille, params, pas_faits, rng):
    """
    Boucle Monte Carlo : fait les pas restants, de pas_faits à n_pas.

    Le masque des iridophores est construit une fois avec
    masque_iridophores(L, bande). À chaque pas : tirer une case, tirer une
    règle, l'appliquer si possible. Tous les PAS_ENTRE_SAUVEGARDES pas, le run est sauvegardé avec
    sauvegarder_run pour pouvoir être relancé.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        Grille de départ, modifiée sur place.
    params : dict
        Paramètres du run (voir preparer_run).
    pas_faits : int
        Nombre de pas déjà effectués (0 pour un nouveau run).
    rng : np.random.Generator
        Générateur aléatoire.

    Renvoie
    -------
    grille : np.ndarray, forme (L, L)
        Grille après n_pas pas.

    Notes
    -----
    Si la boucle est trop lente en Python pur, la compiler avec numba
    (@njit).
    """
    raise NotImplementedError


def sauvegarder_run(grille, params, pas_faits, rng):
    """
    Sauvegarde le run dans le dossier params["sortie"] (créé si besoin).

    Fichiers écrits :
      - grille.npy  : la grille ;
      - params.json : les paramètres, pas_faits, le temps t = pas_faits / L²
                      et l'état du générateur aléatoire ;
      - grille.png  : l'image de la grille (jaune, noir, vide en blanc),
                      avec la zone des iridophores marquée si bande > 0.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau.
    params : dict
        Paramètres du run.
    pas_faits : int
        Nombre de pas effectués.
    rng : np.random.Generator
        Générateur aléatoire, dont l'état est sauvegardé pour que la reprise
        continue la même suite de tirages.
    """
    raise NotImplementedError


def charger_run(dossier):
    """
    Récupère un run sauvegardé par sauvegarder_run, pour le relancer.

    Paramètres
    ----------
    dossier : str
        Dossier du run.

    Renvoie
    -------
    params : dict
        Paramètres relus dans params.json.
    grille : np.ndarray, forme (L, L)
        Grille relue dans grille.npy.
    pas_faits : int
        Nombre de pas déjà effectués.
    rng : np.random.Generator
        Générateur remis dans l'état sauvegardé.
    """
    raise NotImplementedError


if __name__ == "__main__":
    # 1. Paramètres : ligne de commande, et run récupéré s'il existe déjà
    args = lire_arguments()
    params, grille, pas_faits, rng = preparer_run(args)
    print(f"Run {params['sortie']} : {pas_faits} / {params['n_pas']} pas déjà faits")

    if pas_faits >= params["n_pas"]:
        print("Run déjà terminé, rien à faire")
    else:
        # 2. Simulation (sauvegardes régulières)
        grille = simulation_peau(grille, params, pas_faits, rng)

        # 3. Sauvegarde finale
        sauvegarder_run(grille, params, params["n_pas"], rng)
        print(f"Run terminé, résultats dans {params['sortie']}")
