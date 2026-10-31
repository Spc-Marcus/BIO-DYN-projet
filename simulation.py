"""
Un run Monte Carlo du modèle de Bullara & De Decker (2015), lançable sur
GenOuest avec run.slurm (algorithme dans pseudocode.md).

La peau est une grille L x L à bords périodiques. Chaque case vaut VIDE,
JAUNE (xantophore) ou NOIR (mélanophore). À chaque pas on tire une case et
une règle selon les taux, puis on applique la règle si la configuration le
permet. Comme dans l'article : bM = dX = dM = 0 et sX = sM = s, il ne reste
donc que les règles (2), (6), (7) et (8).

La grille et les paramètres sont sauvegardés régulièrement : un run coupé
par la limite de temps Slurm peut être relancé là où il s'est arrêté.

Utilisation :
    python simulation.py --h 16 --lX 2.5
    python simulation.py --reprise out/L100_h16_lX2.5_g0
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
      --L       taille de la grille (défaut 100) ;
      --h       portée de l'interaction à distance, en cases (défaut 16) ;
      --lX      taux de naissance du noir grâce à un jaune à distance h
                (défaut 2.5) ;
      --n_pas   nombre total de tentatives (défaut 10⁹) ;
      --graine  graine du générateur aléatoire (défaut 0) ;
      --sortie  dossier de sauvegarde (défaut out/L{L}_h{h}_lX{lX}_g{graine}) ;
      --reprise dossier d'un run à relancer : les autres options sont alors
                ignorées, les paramètres sont relus dans ce dossier.

    Renvoie
    -------
    args : argparse.Namespace
    """
    raise NotImplementedError


def preparer_run(args):
    """
    Prépare un nouveau run, ou récupère un run existant si --reprise est donné.

    Nouveau run : grille vide, générateur créé à partir de la graine,
    0 pas effectué. Reprise : tout est relu avec charger_run.

    Paramètres
    ----------
    args : argparse.Namespace
        Sortie de lire_arguments.

    Renvoie
    -------
    params : dict
        Clés "L", "h", "lX", "n_pas", "graine", "sortie".
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


def appliquer_regle(grille, i, j, regle, h, rng):
    """
    Applique une règle à la case (i, j) si la configuration le permet.

      - (2) naissance jaune : la case doit être VIDE ;
      - (6) noir tue jaune : case JAUNE et voisin tiré NOIR ;
      - (7) jaune tue noir : case NOIR et voisin tiré JAUNE ;
      - (8) noir naît à distance : case VIDE et case à distance h JAUNE.

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

    À chaque pas : tirer une case, tirer une règle, l'appliquer si possible.
    Tous les PAS_ENTRE_SAUVEGARDES pas, le run est sauvegardé avec
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
      - grille.png  : l'image de la grille (jaune, noir, vide en blanc).

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
    # 1. Paramètres : ligne de commande, ou récupérés en cas de reprise
    args = lire_arguments()
    params, grille, pas_faits, rng = preparer_run(args)
    print(f"Run {params['sortie']} : {pas_faits} / {params['n_pas']} pas déjà faits")

    # 2. Simulation (sauvegardes régulières)
    grille = simulation_peau(grille, params, pas_faits, rng)

    # 3. Sauvegarde finale
    sauvegarder_run(grille, params, params["n_pas"], rng)
    print(f"Run terminé, résultats dans {params['sortie']}")
