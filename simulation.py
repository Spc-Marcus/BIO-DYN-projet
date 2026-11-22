"""
Un run Monte Carlo du modèle de Bullara & De Decker (2015), lancé sur
GenOuest par slurm_pipeline/ (algorithme dans pseudocode.md).

La peau est une grille L x L à bords périodiques. Chaque case vaut VIDE,
JAUNE (xantophore) ou NOIR (mélanophore). À chaque pas on tire une case et
une règle selon les taux, puis on applique la règle si la configuration le
permet.

Règles : les règles (2) à (8) de l'article sont rangées dans un
dictionnaire règle -> taux construit à partir des options de la ligne de
commande (--bX, --bM, --dX, --dM, --s, --lX). Seules les règles de taux non
nul y sont gardées, et chacune est appliquée par sa fonction (table
ACTIONS). Par défaut les taux sont ceux de l'article : bM = dX = dM = 0 et
sX = sM = s, il ne reste donc que les règles (2), (6), (7) et (8).

Naissance au contact (--contact) : une case vide ne reçoit une naissance,
règles (2'), (3') et (8'), que si au moins un de ses voisins directs est
occupé. La peau part alors d'une fraction rho0 de cases occupées
(--rho0, --fX0) au lieu d'être vide, sans quoi elle le resterait.

Blessure : règle (9), de taux --lB, qui vide la case tirée et toutes celles
à distance <= R = alpha h (--alpha), comme une ablation laser. Avec --tB,
une seule blessure au centre de la grille au temps tB, après avoir gardé la
grille d'avant dans grille_avant_blessure.npy.

Réseau : carré par défaut, comme dans l'article (4 voisins). Avec
--reseau hex, réseau hexagonal (6 voisins à distance 1) dans la même
matrice L x L, en coordonnées axiales : la case (i, j) est au point
x = j + i/2, y = (√3/2) i du plan, et ses voisins sont les 4 du carré plus
(1, -1) et (-1, 1). Tout ce qui dépend du réseau est dans RESEAUX.

Iridophores : ils ne sont pas simulés comme des cellules, seulement par un
masque de cases où le noir ne peut pas naître (Fig. 1b-c). Aujourd'hui une
bande horizontale de `bande` lignes (0 = pas d'iridophores) ; une autre
forme se branche en ne modifiant que masque_iridophores.

Le run est sauvegardé régulièrement dans out/L{L}_h{h}_lX{lX}_b{bande}_g{graine},
suivi d'un suffixe par option du modèle différente de sa valeur par défaut
(par exemple out/L100_h16_lX2.5_b0_g0_dX0.1_dM0.1). Relancer la même
commande récupère ce qui y est sauvegardé : un run coupé par la limite de
temps Slurm reprend là où il s'est arrêté, un run terminé n'est pas refait,
et un run sauvegardé avec d'autres paramètres n'est pas repris. Le
générateur aléatoire est un random.Random dont l'état est sauvegardé avec
la grille : un run repris donne exactement la même grille finale qu'un run
fait d'une traite.

Durée : environ 1 µs par pas en Python pur, soit 15 à 20 min pour 10⁹ pas
avec L = 100 (mesuré sur une machine de test).

Utilisation :
    python simulation.py --h 16 --lX 2.5
    python simulation.py --h 16 --lX 2.5 --historique        # garde toutes les grilles
    python simulation.py --h 16 --lX 2.5 --dX 0.1 --dM 0.1   # morts naturelles (4) et (5)
    python simulation.py --h 16 --lX 2.5 --contact --rho0 0.01   # naissances au contact
    python simulation.py --h 16 --lX 2.5 --tB 7000 --n_pas 800000000   # blessure au centre à t = 7000
    python simulation.py --h 16 --lX 2.5 --reseau hex   # réseau hexagonal
"""

import argparse
import json
import math
import os
import random

import numpy as np
import matplotlib
matplotlib.use("Agg")   # pas d'écran sur les nœuds du cluster
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# 1. États d'une case
VIDE = 0
JAUNE = 1
NOIR = 2
COULEURS = ["white", "gold", "black"]   # couleur de VIDE, JAUNE, NOIR sur les images

# 2. Règles du modèle, numérotées comme les équations de l'article
NAISSANCE_JAUNE = 2        # (2) S_i -> X_i
NAISSANCE_NOIR = 3         # (3) S_i -> M_i
MORT_JAUNE = 4             # (4) X_i -> S_i
MORT_NOIR = 5              # (5) M_i -> S_i
NOIR_TUE_JAUNE = 6         # (6) M_{i±1} + X_i -> M_{i±1} + S_i
JAUNE_TUE_NOIR = 7         # (7) X_{i±1} + M_i -> X_{i±1} + S_i
NOIR_NAIT_A_DISTANCE = 8   # (8) X_{i±h} + S_i -> X_{i±h} + M_i
BLESSURE = 9               # (9) A_j -> S_j pour toute case j à distance <= R de i

# 3. Taux fixés comme dans l'article (valeurs par défaut de --bX et --s)
BX = 1.0   # naissance spontanée du jaune
S = 1.0    # élimination entre voisins de couleurs différentes

# 4. Valeurs par défaut des options du modèle. Une option différente de sa
#    valeur par défaut ajoute un suffixe au nom du dossier du run (voir
#    nom_dossier) ; un run sauvegardé avant l'ajout d'une option est relu
#    avec la valeur par défaut de cette option.
DEFAUTS = {
    "bX": BX,           # (2) naissance spontanée du jaune
    "bM": 0.0,          # (3) naissance spontanée du noir
    "dX": 0.0,          # (4) mort du jaune
    "dM": 0.0,          # (5) mort du noir
    "s": S,             # (6) et (7) élimination entre voisins, sX = sM = s
    "reseau": "carre",  # réseau carré (4 voisins) ou hexagonal (6 voisins), voir RESEAUX
    "contact": False,   # naissances seulement au contact d'une case occupée
    "rho0": 0.0,        # fraction de cases occupées au départ
    "fX0": 1.0,         # part de jaunes parmi ces cases
    "lB": 0.0,          # (9) blessure
    "alpha": 1.0,       # rayon de la blessure R = alpha h
    "tB": None,         # temps de la blessure unique au centre (None : pas de blessure unique)
}

# 5. Nombre de pas entre deux sauvegardes (pour pouvoir relancer).
#    Avec --historique, chaque sauvegarde est aussi gardée dans historique/
#    (pour une timeline ou un GIF avec timeline.py).
PAS_ENTRE_SAUVEGARDES = 10**7

# 6. Voisins directs sur le réseau carré : haut, bas, gauche, droite. Ne pas
#    changer l'ordre : il fixe le voisin associé à chaque tirage. Sur le
#    réseau hexagonal (coordonnées axiales), les mêmes plus une diagonale.
VOISINS_CARRE = ((1, 0), (-1, 0), (0, 1), (0, -1))
VOISINS_HEX = VOISINS_CARRE + ((1, -1), (-1, 1))
RACINE_3 = math.sqrt(3)


def lire_arguments():
    """
    Lit les paramètres du run sur la ligne de commande.

    Options :
      --L       taille de la grille, entier (défaut 100) ;
      --h       portée de l'interaction à distance, en cases, entier
                (défaut 16) ;
      --lX      taux de naissance du noir grâce à un jaune à distance h,
                règle (8) (défaut 2.5) ;
      --bande   largeur de la bande d'iridophores en lignes, entier
                (défaut 0 = pas d'iridophores) ;
      --n_pas   nombre total de tentatives, entier (défaut 1000000000) ;
      --graine  graine du générateur aléatoire, entier (défaut 0) ;
      --historique  garde une copie de la grille à chaque sauvegarde
                (désactivé par défaut) ;
      --bX, --bM, --dX, --dM, --s  taux des règles (2), (3), (4), (5), et
                (6)-(7) (défauts 1, 0, 0, 0, 1 : ceux de l'article, DEFAUTS) ;
      --reseau  carre (défaut, 4 voisins) ou hex (6 voisins) ;
      --contact naissances (2), (3) et (8) seulement si un voisin direct de
                la case est occupé (désactivé par défaut) ;
      --rho0    fraction de cases occupées au départ, entre 0 et 1
                (défaut 0 = peau vide) ;
      --fX0     part de jaunes parmi ces cases, le reste en noir, entre 0
                et 1 (défaut 1 = que des jaunes) ;
      --lB      taux de la blessure, règle (9) (défaut 0 = pas de blessure) ;
      --alpha   rayon de la blessure en périodes : R = alpha h (défaut 1) ;
      --tB      temps (en unités t = pas / L²) d'une blessure unique au
                centre de la grille (défaut : pas de blessure unique).

    Le programme s'arrête avec un message d'erreur si un taux est négatif,
    si tous les taux sont nuls, si rho0 ou fX0 sort de [0, 1], si
    --contact est donné sans rho0 > 0 (une peau vide le resterait), si une
    blessure n'a pas 1 <= alpha h < L/2 (au moins une case, sans que le
    disque se recouvre sur le tore), ou si la blessure unique n'a pas lieu
    pendant le run (round(tB L²) doit être entre 0 et n_pas - 1).

    Renvoie
    -------
    args : argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="Un run Monte Carlo (Bullara & De Decker 2015)")
    parser.add_argument("--L", type=int, default=100, help="taille de la grille")
    parser.add_argument("--h", type=int, default=16, help="portée de l'interaction à distance")
    parser.add_argument("--lX", type=float, default=2.5, help="(8) naissance du noir à distance h")
    parser.add_argument("--bande", type=int, default=0, help="largeur de la bande d'iridophores")
    parser.add_argument("--n_pas", type=int, default=1000000000, help="nombre total de tentatives")
    parser.add_argument("--graine", type=int, default=0, help="graine du générateur aléatoire")
    parser.add_argument("--historique", action="store_true", help="garde la grille à chaque sauvegarde")
    parser.add_argument("--bX", type=float, default=DEFAUTS["bX"], help="(2) naissance spontanée du jaune")
    parser.add_argument("--bM", type=float, default=DEFAUTS["bM"], help="(3) naissance spontanée du noir")
    parser.add_argument("--dX", type=float, default=DEFAUTS["dX"], help="(4) mort du jaune")
    parser.add_argument("--dM", type=float, default=DEFAUTS["dM"], help="(5) mort du noir")
    parser.add_argument("--s", type=float, default=DEFAUTS["s"], help="(6) et (7) élimination entre voisins")
    parser.add_argument("--reseau", choices=list(RESEAUX), default=DEFAUTS["reseau"], help="réseau de la peau")
    parser.add_argument("--contact", action="store_true", help="naissances seulement au contact d'une case occupée")
    parser.add_argument("--rho0", type=float, default=DEFAUTS["rho0"], help="fraction de cases occupées au départ")
    parser.add_argument("--fX0", type=float, default=DEFAUTS["fX0"], help="part de jaunes parmi ces cases")
    parser.add_argument("--lB", type=float, default=DEFAUTS["lB"], help="(9) blessure")
    parser.add_argument("--alpha", type=float, default=DEFAUTS["alpha"], help="rayon de la blessure R = alpha h")
    parser.add_argument("--tB", type=float, default=DEFAUTS["tB"], help="temps de la blessure unique au centre")
    args = parser.parse_args()

    # 1. Taux positifs ou nuls, dont au moins un non nul
    taux = [args.bX, args.bM, args.dX, args.dM, args.s, args.lX, args.lB]
    if min(taux) < 0 or max(taux) == 0:
        parser.error("les taux doivent être positifs ou nuls, et au moins un non nul")

    # 2. Départ : des fractions, et des cases occupées s'il faut un contact pour naître
    if not (0 <= args.rho0 <= 1 and 0 <= args.fX0 <= 1):
        parser.error("--rho0 et --fX0 doivent être entre 0 et 1")
    if args.contact and args.rho0 == 0:
        parser.error("--contact demande --rho0 > 0 : sinon la peau vide le reste")

    # 3. Blessure : au moins une case, sans recouvrement sur le tore, et pendant le run
    if (args.lB > 0 or args.tB is not None) and not 1 <= args.alpha * args.h < args.L / 2:
        parser.error("avec une blessure, il faut 1 <= alpha h < L/2")
    if args.tB is not None and not 0 <= round(args.tB * args.L**2) < args.n_pas:
        parser.error("la blessure unique doit avoir lieu pendant le run : 0 <= round(tB L²) < n_pas")
    return args


def nom_dossier(params):
    """
    Construit le nom du dossier du run.

    out/L{L}_h{h}_lX{lX}_b{bande}_g{graine}, suivi d'un suffixe par option
    du modèle différente de sa valeur par défaut (DEFAUTS), dans cet ordre :
      - _bX{bX}, _bM{bM}, _dX{dX}, _dM{dM}, _s{s} pour les taux ;
      - _hex pour le réseau hexagonal ;
      - _c pour --contact, _r{rho0} et, si rho0 > 0, _fx{fX0} pour le départ ;
      - _B{lB}, _a{alpha} (s'il y a une blessure) et _tB{tB} pour la blessure.
    Avec les options par défaut, c'est le nom des runs lancés avant ces
    options : ils sont retrouvés et repris.

    Paramètres
    ----------
    params : dict
        Paramètres du run (voir preparer_run).

    Renvoie
    -------
    sortie : str
        Chemin du dossier, relatif à la racine du projet.
    """
    sortie = f"out/L{params['L']}_h{params['h']}_lX{params['lX']}_b{params['bande']}_g{params['graine']}"
    for cle in ("bX", "bM", "dX", "dM", "s"):
        if params[cle] != DEFAUTS[cle]:
            sortie += f"_{cle}{params[cle]}"
    if params["reseau"] != DEFAUTS["reseau"]:
        sortie += f"_{params['reseau']}"
    if params["contact"]:
        sortie += "_c"
    if params["rho0"] != DEFAUTS["rho0"]:
        sortie += f"_r{params['rho0']}"
        if params["fX0"] != DEFAUTS["fX0"]:
            sortie += f"_fx{params['fX0']}"
    if params["lB"] != DEFAUTS["lB"]:
        sortie += f"_B{params['lB']}"
    if (params["lB"] > 0 or params["tB"] is not None) and params["alpha"] != DEFAUTS["alpha"]:
        sortie += f"_a{params['alpha']}"
    if params["tB"] is not None:
        sortie += f"_tB{params['tB']}"
    return sortie


def preparer_run(args):
    """
    Prépare le run : nouveau, ou repris depuis son dossier s'il existe déjà.

    Le dossier du run est donné par nom_dossier ; il ne dépend ni de n_pas
    ni de historique.
      - Pas de params.json dans ce dossier : nouveau run, générateur
        random.Random(graine), grille de départ initialiser_grille (vide si
        rho0 = 0), 0 pas effectué.
      - Sinon : reprise, la grille, pas_faits et l'état du générateur sont
        récupérés avec charger_run. Les paramètres sauvegardés doivent être
        ceux de la ligne de commande, sauf n_pas et historique qui peuvent
        changer, ce qui permet de prolonger un run ou de commencer à garder
        son historique. Sinon le programme s'arrête plutôt que de reprendre
        un run avec d'autres paramètres. Une option absente du params.json
        (run lancé avant l'ajout de l'option) y vaut sa valeur par défaut.

    Paramètres
    ----------
    args : argparse.Namespace
        Sortie de lire_arguments.

    Renvoie
    -------
    params : dict
        Clés "L", "h", "lX", "bande", "n_pas", "graine", "historique",
        celles de DEFAUTS ("bX", "bM", "dX", "dM", "s", "reseau",
        "contact", "rho0", "fX0", "lB", "alpha", "tB") et "sortie".
    grille : np.ndarray uint8, forme (L, L)
        Grille de départ (initialiser_grille, ou relue pour une reprise).
    pas_faits : int
        Nombre de pas déjà effectués.
    rng : random.Random
        Générateur aléatoire.

    Lève
    ----
    SystemExit
        Si le run sauvegardé dans le dossier a d'autres paramètres que la
        ligne de commande.
    """
    # 1. Paramètres de la ligne de commande
    params = {
        "L": args.L,
        "h": args.h,
        "lX": args.lX,
        "bande": args.bande,
        "n_pas": args.n_pas,
        "graine": args.graine,
        "historique": args.historique,
        "bX": args.bX,
        "bM": args.bM,
        "dX": args.dX,
        "dM": args.dM,
        "s": args.s,
        "reseau": args.reseau,
        "contact": args.contact,
        "rho0": args.rho0,
        "fX0": args.fX0,
        "lB": args.lB,
        "alpha": args.alpha,
        "tB": args.tB,
    }
    params["sortie"] = nom_dossier(params)

    # 2. Reprise si le run a déjà été sauvegardé, avec les mêmes paramètres
    if os.path.exists(os.path.join(params["sortie"], "params.json")):
        params_sauvegardes, grille, pas_faits, rng = charger_run(params["sortie"])
        differents = [cle for cle in params
                      if cle not in ("n_pas", "historique", "sortie")
                      and params_sauvegardes.get(cle, DEFAUTS.get(cle)) != params[cle]]
        if differents:
            detail = ", ".join(f"{cle} = {params_sauvegardes.get(cle, DEFAUTS.get(cle))} au lieu de {params[cle]}"
                               for cle in differents)
            raise SystemExit(f"Erreur : {params['sortie']} a été sauvegardé avec d'autres paramètres ({detail})")
        return params, grille, pas_faits, rng

    # 3. Sinon nouveau run (le générateur d'abord : il tire les cases occupées au départ)
    rng = random.Random(args.graine)
    grille = initialiser_grille(args.L, args.rho0, args.fX0, rng)
    return params, grille, 0, rng


def initialiser_grille(L, rho0, fX0, rng):
    """
    Crée la grille de départ : peau vide, ou avec une fraction rho0 de cases
    occupées.

    Avec rho0 > 0 : n0 = max(1, round(rho0 L²)) cases tirées sans remise
    (rng.sample, donc exactement n0), puis chacune JAUNE avec la
    probabilité fX0, NOIR sinon. Avec rho0 = 0 rien n'est tiré : la suite
    des tirages de la simulation est celle d'une peau vide.

    Paramètres
    ----------
    L : int
        Taille de la grille (L x L cases).
    rho0 : float
        Fraction de cases occupées au départ, entre 0 et 1.
    fX0 : float
        Part de jaunes parmi les cases occupées, entre 0 et 1.
    rng : random.Random
        Générateur aléatoire.

    Renvoie
    -------
    grille : np.ndarray uint8, forme (L, L)
        Tableau de VIDE, avec n0 cases JAUNE ou NOIR si rho0 > 0.
    """
    grille = np.full((L, L), VIDE, dtype=np.uint8)
    if rho0 > 0:
        n0 = max(1, round(rho0 * L**2))
        for c in rng.sample(range(L**2), n0):
            grille[c // L, c % L] = JAUNE if rng.random() < fX0 else NOIR
    return grille


def masque_iridophores(L, bande):
    """
    Construit le masque des cases couvertes par les iridophores.

    Sur ces cases le noir ne peut pas naître : les règles (3) et (8) y sont
    bloquées (Fig. 1b-c). Pour l'instant le masque est une bande
    horizontale de `bande` lignes au milieu de la grille : les lignes
    (L - bande) // 2 à (L - bande) // 2 + bande - 1. Pour une autre forme
    d'iridophores, il suffit de modifier cette fonction : le reste du code
    (simulation et image) ne lit que le masque.

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
    noir_interdit = np.zeros((L, L), dtype=bool)
    if bande > 0:
        debut = (L - bande) // 2
        noir_interdit[debut:debut + bande, :] = True
    return noir_interdit


def regles_actives(params):
    """
    Construit le dictionnaire des règles actives et de leur probabilité.

    Chaque règle reçoit le taux donné par son option (DEFAUTS : taux de
    l'article, et pas de blessure), dans l'ordre des équations (2) à (9).
    La blessure unique de --tB n'est pas une règle tirée : elle n'est pas
    dans ce dictionnaire (voir simulation_peau). Les règles de taux nul
    sont retirées : placée en dernier, une règle de taux nul serait sinon
    renvoyée par le repli d'arrondi de tirer_regle. Puis
    P_k = taux_k / somme(taux). Avec les taux par défaut il reste (2), (6),
    (7) et (8), dans cet ordre : les tirages sont les mêmes que lorsque ces
    quatre règles étaient écrites en dur.

    Pour ajouter une règle : sa fonction dans ACTIONS et son taux ici.

    Paramètres
    ----------
    params : dict
        Paramètres du run (voir preparer_run).

    Renvoie
    -------
    proba : dict
        Règle -> probabilité de la choisir à chaque pas, dans l'ordre des
        équations, de somme 1.
    """
    # 1. Taux de chaque règle, donnés par les options
    taux = {
        NAISSANCE_JAUNE: params["bX"],        # (2) --bX
        NAISSANCE_NOIR: params["bM"],         # (3) --bM
        MORT_JAUNE: params["dX"],             # (4) --dX
        MORT_NOIR: params["dM"],              # (5) --dM
        NOIR_TUE_JAUNE: params["s"],          # (6) --s
        JAUNE_TUE_NOIR: params["s"],          # (7) --s
        NOIR_NAIT_A_DISTANCE: params["lX"],   # (8) --lX
        BLESSURE: params["lB"],               # (9) --lB
    }

    # 2. Règles actives seulement
    taux = {regle: t for regle, t in taux.items() if t > 0}

    # 3. Probabilité de choisir chaque règle : P_k = taux_k / somme(taux)
    somme = sum(taux.values())
    return {regle: t / somme for regle, t in taux.items()}


def tirer_regle(proba, rng):
    """
    Tire une règle au hasard selon les probabilités données.

    Tirage pondéré : on découpe [0, 1[ en segments de longueur proba[regle],
    dans l'ordre du dictionnaire, et on renvoie la règle du segment dans
    lequel tombe un nombre uniforme u. Si les arrondis font que u dépasse la
    somme cumulée, la dernière règle est renvoyée.

    Paramètres
    ----------
    proba : dict
        Règle -> probabilité, sortie de regles_actives.
    rng : random.Random
        Générateur aléatoire.

    Renvoie
    -------
    regle : int
        Numéro de la règle tirée (NAISSANCE_JAUNE, ...).
    """
    u = rng.random()
    cumul = 0.0
    for regle, p in proba.items():
        cumul = cumul + p
        if u < cumul:
            return regle
    return regle   # u très proche de 1 : erreurs d'arrondi sur le cumul


def voisin_au_hasard(grille, i, j, voisins, rng):
    """
    Renvoie l'état d'un des voisins directs de la case (i, j).

    Le voisin est tiré uniformément parmi les z voisins du réseau (haut,
    bas, gauche, droite sur le réseau carré, et les deux diagonales
    (1, -1), (-1, 1) en plus sur le réseau hexagonal) avec un seul nombre
    aléatoire ; les bords sont périodiques (mod L). Utilisé par les règles
    (6) et (7).

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État courant de la peau.
    i, j : int
        Coordonnées de la case.
    voisins : tuple de (di, dj)
        Décalages des voisins directs, RESEAUX[...]["voisins"].
    rng : random.Random
        Générateur aléatoire.

    Renvoie
    -------
    etat : int
        VIDE, JAUNE ou NOIR.
    """
    L = grille.shape[0]
    di, dj = voisins[int(rng.random() * len(voisins))]
    return grille[(i + di) % L, (j + dj) % L]


def contact_ok(grille, i, j, voisins):
    """
    Dit si la case (i, j) touche une case occupée : au moins un de ses
    voisins directs n'est pas VIDE (bords périodiques). Ne tire rien.
    Utilisé par les naissances (2'), (3') et (8') avec --contact.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État courant de la peau.
    i, j : int
        Coordonnées de la case.
    voisins : tuple de (di, dj)
        Décalages des voisins directs, RESEAUX[...]["voisins"].

    Renvoie
    -------
    contact : bool
    """
    L = grille.shape[0]
    for di, dj in voisins:
        if grille[(i + di) % L, (j + dj) % L] != VIDE:
            return True
    return False


def case_a_distance(grille, i, j, h, rng):
    """
    Renvoie l'état d'une case située à distance h de la case (i, j), sur le
    réseau carré.

    On tire un angle alpha uniforme dans [0, 2π[ et on regarde la case
    (round(i + h cos alpha), round(j + h sin alpha)), modulo L, avec le
    round de Python (entier le plus proche). Utilisé par la règle (8).

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État courant de la peau.
    i, j : int
        Coordonnées de la case.
    h : int
        Portée de l'interaction, en cases.
    rng : random.Random
        Générateur aléatoire.

    Renvoie
    -------
    etat : int
        VIDE, JAUNE ou NOIR.
    """
    L = grille.shape[0]
    alpha = rng.random() * 2 * math.pi
    k = round(i + h * math.cos(alpha)) % L
    l = round(j + h * math.sin(alpha)) % L
    return grille[k, l]


def vers_image_carre(grille):
    """
    Image d'une grille du réseau carré : la grille elle-même, une case par
    pixel (à afficher avec aspect = 1).

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)

    Renvoie
    -------
    image : np.ndarray, forme (L, L)
    """
    return grille


def d2_carre(di, dj):
    """
    Distance² (en cases) d'un décalage (di, dj) sur le réseau carré.

    Paramètres
    ----------
    di, dj : int
        Décalage en lignes et en colonnes.

    Renvoie
    -------
    d2 : int
        di² + dj².
    """
    return di * di + dj * dj


def d2_hex(di, dj):
    """
    Distance² (en cases) d'un décalage axial (di, dj) sur le réseau
    hexagonal : avec x = dj + di/2 et y = (√3/2) di, x² + y² vaut
    di² + di dj + dj². Les 6 voisins sont à distance 1.

    Paramètres
    ----------
    di, dj : int
        Décalage en coordonnées axiales.

    Renvoie
    -------
    d2 : int
        di² + di dj + dj².
    """
    return di * di + di * dj + dj * dj


def arrondi_hex(dx, dy):
    """
    Décalage axial (di, dj) de la case du réseau hexagonal la plus proche
    du déplacement (dx, dy) dans le plan.

    On repasse en coordonnées axiales fractionnaires (r = di, q = dj), on
    arrondit les trois coordonnées « cube » q, r, s = -q - r, puis on
    recalcule celle qui a le plus grand écart d'arrondi pour garder
    q + r + s = 0 (arrondi hexagonal de Red Blob Games, Hexagonal Grids).

    Paramètres
    ----------
    dx, dy : float
        Déplacement dans le plan, en cases.

    Renvoie
    -------
    di, dj : int
        Décalage axial de la case la plus proche.
    """
    r = 2 * dy / RACINE_3   # di fractionnaire
    q = dx - r / 2          # dj fractionnaire
    s = -q - r              # 3e coordonnée « cube » : q + r + s = 0
    rq, rr, rs = round(q), round(r), round(s)
    eq, er, es = abs(rq - q), abs(rr - r), abs(rs - s)
    if eq > er and eq > es:
        rq = -rr - rs
    elif er > es:
        rr = -rq - rs
    return rr, rq


def case_a_distance_hex(grille, i, j, h, rng):
    """
    Renvoie l'état d'une case située à distance h de la case (i, j), sur le
    réseau hexagonal.

    Comme sur le réseau carré : angle alpha uniforme dans [0, 2π[, point
    (h cos alpha, h sin alpha) dans le plan, puis la case la plus proche
    (arrondi_hex), modulo L. Utilisé par la règle (8).

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État courant de la peau.
    i, j : int
        Coordonnées axiales de la case.
    h : int
        Portée de l'interaction, en cases.
    rng : random.Random
        Générateur aléatoire.

    Renvoie
    -------
    etat : int
        VIDE, JAUNE ou NOIR.
    """
    L = grille.shape[0]
    alpha = rng.random() * 2 * math.pi
    di, dj = arrondi_hex(h * math.cos(alpha), h * math.sin(alpha))
    return grille[(i + di) % L, (j + dj) % L]


def vers_image_hex(grille):
    """
    Image « briques » d'une grille du réseau hexagonal, de forme L x 2L.

    Chaque case occupe 2 pixels d'une ligne et la ligne i est décalée de i
    pixels (i/2 case), modulo 2L : chaque brique touche alors exactement
    ses 6 voisines. À afficher avec aspect = √3 (un pixel fait 1/2 case de
    large et √3/2 de haut). Le haut et le bas de l'image se recollent avec
    un décalage de L/2 cases : une rayure qui sort en bas revient décalée
    en haut.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)

    Renvoie
    -------
    image : np.ndarray, forme (L, 2L), même type que grille
    """
    L = grille.shape[0]
    I = np.arange(L)[:, None]
    colonnes = (2 * np.arange(L)[None, :] + I) % (2 * L)
    image = np.empty((L, 2 * L), dtype=grille.dtype)
    image[I, colonnes] = grille
    image[I, (colonnes + 1) % (2 * L)] = grille
    return image


# 7. Réseaux : tout ce qui dépend de la géométrie de la peau, lu par les
#    règles et par l'image. Le réseau carré est celui de l'article ; le
#    réseau hexagonal garde la matrice L x L en coordonnées axiales.
#      voisins         : décalages (di, dj) des voisins directs ;
#      d2              : distance² d'un décalage (di, dj), pour la blessure ;
#      case_a_distance : case à distance h, pour la règle (8) ;
#      vers_image      : image de la grille pour plt.imshow ;
#      aspect          : rapport hauteur / largeur d'un pixel de l'image.
RESEAUX = {
    "carre": {
        "voisins": VOISINS_CARRE,
        "d2": d2_carre,
        "case_a_distance": case_a_distance,
        "vers_image": vers_image_carre,
        "aspect": 1.0,
    },
    "hex": {
        "voisins": VOISINS_HEX,
        "d2": d2_hex,
        "case_a_distance": case_a_distance_hex,
        "vers_image": vers_image_hex,
        "aspect": RACINE_3,
    },
}


def disque(d2, R):
    """
    Décalages de toutes les cases à distance <= R de la case (0, 0), centre
    compris : le disque D_R de la blessure, calculé une fois par run.

    On parcourt le carré |di|, |dj| <= ceil(2R/√3) + 1, qui contient le
    disque sur les deux réseaux (sur le réseau hexagonal, d2 >= (3/4) di²
    et (3/4) dj²), et on garde les décalages avec d2(di, dj) <= R².
    Le disque autour de (i, j) est alors grille[(i + DI) % L, (j + DJ) % L].

    Paramètres
    ----------
    d2 : fonction
        Distance² d'un décalage, celle du réseau (RESEAUX[...]["d2"]).
    R : float
        Rayon du disque, en cases.

    Renvoie
    -------
    DI, DJ : np.ndarray d'entiers, forme (n,)
        Décalages en lignes et en colonnes des n cases du disque.
    """
    n = int(math.ceil(2 * R / RACINE_3)) + 1
    decalages = [(di, dj) for di in range(-n, n + 1) for dj in range(-n, n + 1) if d2(di, dj) <= R * R]
    DI = np.array([di for di, dj in decalages])
    DJ = np.array([dj for di, dj in decalages])
    return DI, DJ


def naissance_jaune(grille, i, j, modele, rng):
    """
    (2) S_i -> X_i : la case doit être VIDE. Avec --contact, (2') : elle
    doit aussi toucher une case occupée.
    """
    if grille[i, j] == VIDE and (not modele["contact"] or contact_ok(grille, i, j, modele["voisins"])):
        grille[i, j] = JAUNE
        return True
    return False


def naissance_noir(grille, i, j, modele, rng):
    """
    (3) S_i -> M_i : case VIDE et hors iridophores. Avec --contact, (3') :
    elle doit aussi toucher une case occupée.
    """
    if (grille[i, j] == VIDE and not modele["noir_interdit"][i, j]
            and (not modele["contact"] or contact_ok(grille, i, j, modele["voisins"]))):
        grille[i, j] = NOIR
        return True
    return False


def mort_jaune(grille, i, j, modele, rng):
    """(4) X_i -> S_i : la case doit être JAUNE."""
    if grille[i, j] == JAUNE:
        grille[i, j] = VIDE
        return True
    return False


def mort_noir(grille, i, j, modele, rng):
    """(5) M_i -> S_i : la case doit être NOIR."""
    if grille[i, j] == NOIR:
        grille[i, j] = VIDE
        return True
    return False


def noir_tue_jaune(grille, i, j, modele, rng):
    """(6) un noir voisin élimine le jaune : case JAUNE et voisin tiré NOIR."""
    if grille[i, j] == JAUNE and voisin_au_hasard(grille, i, j, modele["voisins"], rng) == NOIR:
        grille[i, j] = VIDE
        return True
    return False


def jaune_tue_noir(grille, i, j, modele, rng):
    """(7) un jaune voisin élimine le noir : case NOIR et voisin tiré JAUNE."""
    if grille[i, j] == NOIR and voisin_au_hasard(grille, i, j, modele["voisins"], rng) == JAUNE:
        grille[i, j] = VIDE
        return True
    return False


def noir_nait_a_distance(grille, i, j, modele, rng):
    """
    (8) un jaune à distance h fait naître un noir : case VIDE, hors
    iridophores, et case à distance h JAUNE. Avec --contact, (8') : la case
    doit aussi toucher une case occupée, testé avant de tirer la case à
    distance pour ne pas tirer d'angle inutile.
    """
    if (grille[i, j] == VIDE and not modele["noir_interdit"][i, j]
            and (not modele["contact"] or contact_ok(grille, i, j, modele["voisins"]))
            and modele["case_a_distance"](grille, i, j, modele["h"], rng) == JAUNE):
        grille[i, j] = NOIR
        return True
    return False


def blessure(grille, i, j, modele, rng):
    """
    (9) blessure : la case (i, j) et toutes celles à distance <= R
    deviennent VIDE, quel que soit leur état (bords périodiques). Le masque
    des iridophores n'est pas touché. Ne tire rien.
    """
    DI, DJ = modele["disque"]
    L = grille.shape[0]
    grille[(i + DI) % L, (j + DJ) % L] = VIDE
    return True


# 8. Table règle -> action : la fonction qui applique chaque règle
ACTIONS = {
    NAISSANCE_JAUNE: naissance_jaune,
    NAISSANCE_NOIR: naissance_noir,
    MORT_JAUNE: mort_jaune,
    MORT_NOIR: mort_noir,
    NOIR_TUE_JAUNE: noir_tue_jaune,
    JAUNE_TUE_NOIR: jaune_tue_noir,
    NOIR_NAIT_A_DISTANCE: noir_nait_a_distance,
    BLESSURE: blessure,
}


def appliquer_regle(grille, i, j, regle, modele, rng):
    """
    Applique une règle à la case (i, j) si la configuration le permet.

    La règle est appliquée par sa fonction, trouvée dans ACTIONS :
      - (2) naissance jaune : la case doit être VIDE ;
      - (3) naissance noire : case VIDE, hors iridophores ;
      - (4) et (5) mort du jaune, du noir : case JAUNE, NOIR ;
      - (6) noir tue jaune : case JAUNE et voisin tiré NOIR ;
      - (7) jaune tue noir : case NOIR et voisin tiré JAUNE ;
      - (8) noir naît à distance : case VIDE, hors iridophores
        (noir_interdit[i, j] faux) et case à distance h JAUNE ;
      - avec --contact, (2'), (3') et (8') : les naissances demandent en
        plus qu'un voisin direct de la case soit occupé (contact_ok) ;
      - (9) blessure : le disque de rayon R autour de la case devient VIDE,
        toujours possible.

    Sinon la tentative est ratée et la grille n'est pas modifiée. Le
    voisin (règles 6 et 7) ou la case à distance (règle 8) n'est tiré que
    si la case a le bon état, comme dans pseudocode.md : le nombre de
    tirages aléatoires dépend donc de la configuration.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État courant de la peau, modifié sur place.
    i, j : int
        Coordonnées de la case tirée.
    regle : int
        Numéro de la règle (NAISSANCE_JAUNE, ...).
    modele : dict
        Données du run lues par les règles, préparées une fois par
        simulation_peau : "h" (portée de l'interaction à distance),
        "noir_interdit" (masque des iridophores, voir masque_iridophores),
        "contact" (option --contact), "voisins" et "case_a_distance"
        (ceux du réseau, voir RESEAUX), "disque" (décalages (DI, DJ) du
        disque de la blessure, voir disque ; None sans blessure).
    rng : random.Random
        Générateur aléatoire.

    Renvoie
    -------
    reussie : bool
        True si la grille a été modifiée.
    """
    return ACTIONS[regle](grille, i, j, modele, rng)


def blessure_unique(grille, params, modele, pas):
    """
    Blessure unique de --tB, faite par simulation_peau au pas round(tB L²),
    avant le tirage de ce pas.

    La grille est d'abord sauvegardée dans grille_avant_blessure.npy (écrit
    sous un nom temporaire puis renommé, comme grille.npy), puis le disque
    de rayon R centré en (L // 2, L // 2) est vidé comme par la règle (9).

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau, modifié sur place.
    params : dict
        Paramètres du run.
    modele : dict
        Données lues par les règles (voir appliquer_regle), dont "disque".
    pas : int
        Pas de la blessure, pour le message du journal.
    """
    dossier = params["sortie"]
    os.makedirs(dossier, exist_ok=True)
    np.save(os.path.join(dossier, "grille_avant_blessure.tmp.npy"), grille)
    os.replace(os.path.join(dossier, "grille_avant_blessure.tmp.npy"),
               os.path.join(dossier, "grille_avant_blessure.npy"))

    L = params["L"]
    blessure(grille, L // 2, L // 2, modele, None)
    print(f"Blessure au centre : pas {pas}, t = {pas / L**2:.0f}", flush=True)


def simulation_peau(grille, params, pas_faits, rng):
    """
    Boucle Monte Carlo : fait les pas restants, de pas_faits à n_pas.

    Une fois au début : les données lues par les règles (modele : le
    masque des iridophores, construit avec masque_iridophores(L, bande),
    l'option --contact, les voisins et la case à distance du réseau, le
    disque de la blessure de rayon R = alpha h) et le dictionnaire des
    règles actives avec leur probabilité (regles_actives). À chaque pas :
    tirer une case (i, j) uniforme, tirer une règle, l'appliquer si
    possible. Avec --tB, au pas round(tB L²) et avant ce tirage, la
    blessure unique au centre (blessure_unique) : la boucle reprenant à
    pas_faits, elle n'est faite qu'une fois, même si le run est repris ou
    prolongé. Tous les PAS_ENTRE_SAUVEGARDES pas, le run est sauvegardé
    avec sauvegarder_run pour pouvoir être relancé, sauf au dernier pas :
    la sauvegarde finale est faite par le programme principal.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        Grille de départ, modifiée sur place.
    params : dict
        Paramètres du run (voir preparer_run).
    pas_faits : int
        Nombre de pas déjà effectués (0 pour un nouveau run).
    rng : random.Random
        Générateur aléatoire.

    Renvoie
    -------
    grille : np.ndarray, forme (L, L)
        Grille après n_pas pas.

    Notes
    -----
    Environ 1 µs par pas en Python pur. Pour aller plus vite, numba ne
    compile pas random.Random : il faudrait passer la boucle et les tirages
    à numba avec np.random, en gardant une façon de sauvegarder l'état du
    générateur pour la reprise.
    """
    # 1. Paramètres du run et données lues par les règles
    L = params["L"]
    reseau = RESEAUX[params["reseau"]]
    modele = {
        "h": params["h"],
        "noir_interdit": masque_iridophores(L, params["bande"]),
        "contact": params["contact"],
        "voisins": reseau["voisins"],
        "case_a_distance": reseau["case_a_distance"],
        "disque": None,
    }
    if params["lB"] > 0 or params["tB"] is not None:
        modele["disque"] = disque(reseau["d2"], params["alpha"] * params["h"])

    # 2. Règles actives et probabilité de choisir chacune
    proba = regles_actives(params)

    # 3. Pas de la blessure unique (-1 sans --tB : jamais atteint)
    pas_blessure = round(params["tB"] * L**2) if params["tB"] is not None else -1

    # 4. Boucle Monte Carlo
    for pas in range(pas_faits, params["n_pas"]):
        # a. Blessure unique au centre, avant le tirage de ce pas
        if pas == pas_blessure:
            blessure_unique(grille, params, modele, pas)

        # b. Tirage d'une case et d'une règle
        i = int(rng.random() * L)
        j = int(rng.random() * L)
        regle = tirer_regle(proba, rng)

        # c. Application de la règle si la configuration le permet
        appliquer_regle(grille, i, j, regle, modele, rng)

        # d. Sauvegarde régulière, pour pouvoir relancer (la finale est faite après la boucle)
        if (pas + 1) % PAS_ENTRE_SAUVEGARDES == 0 and pas + 1 < params["n_pas"]:
            sauvegarder_run(grille, params, pas + 1, rng)

    return grille


def sauvegarder_run(grille, params, pas_faits, rng):
    """
    Sauvegarde le run dans le dossier params["sortie"] (créé si besoin),
    puis affiche une ligne « Sauvegarde : pas_faits / n_pas pas, t = … »
    pour suivre l'avancement dans le journal Slurm.

    grille.npy et params.json sont d'abord écrits sous des noms temporaires
    puis renommés : un job tué pendant la sauvegarde laisse la sauvegarde
    précédente intacte.

    Fichiers écrits :
      - grille.npy  : la grille ;
      - params.json : les paramètres, pas_faits, le temps t = pas_faits / L²
                      et l'état du générateur aléatoire (rng.getstate(),
                      sous la clé "etat_generateur") ;
      - grille.png  : l'image de la grille (vers_image du réseau, couleurs
                      COULEURS : vide en blanc, jaune, noir), avec la zone
                      des iridophores en gris argenté semi-transparent si le
                      masque n'est pas vide ;
      - historique/grille_{pas_faits:012d}.npz : seulement si
                      params["historique"] est vrai. Une copie de la grille
                      à chaque sauvegarde, jamais écrasée, lue par
                      timeline.py. Elle contient deux matrices binaires,
                      "jaune" (grille == JAUNE) et "noir" (grille == NOIR),
                      compactées ligne par ligne avec
                      np.packbits(..., axis=-1) puis écrites avec
                      np.savez_compressed. Le vide n'est pas stocké : c'est
                      ni jaune ni noir.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau.
    params : dict
        Paramètres du run.
    pas_faits : int
        Nombre de pas effectués.
    rng : random.Random
        Générateur aléatoire, dont l'état est sauvegardé pour que la reprise
        continue la même suite de tirages.
    """
    dossier = params["sortie"]
    os.makedirs(dossier, exist_ok=True)
    L = params["L"]
    t = pas_faits / L**2

    # 1. Grille et paramètres, écrits dans des fichiers temporaires puis
    #    renommés : un job tué pendant l'écriture ne laisse pas de fichier à moitié écrit
    np.save(os.path.join(dossier, "grille.tmp.npy"), grille)

    infos = dict(params)
    infos["pas_faits"] = pas_faits
    infos["t"] = t
    infos["etat_generateur"] = rng.getstate()
    with open(os.path.join(dossier, "params.tmp.json"), "w") as f:
        json.dump(infos, f, indent=2)

    os.replace(os.path.join(dossier, "grille.tmp.npy"), os.path.join(dossier, "grille.npy"))
    os.replace(os.path.join(dossier, "params.tmp.json"), os.path.join(dossier, "params.json"))

    # 2. Image de la grille (vers_image du réseau), iridophores en gris argenté par-dessus
    reseau = RESEAUX[params["reseau"]]
    plt.figure(figsize=(6, 6))
    plt.imshow(reseau["vers_image"](grille), cmap=ListedColormap(COULEURS), vmin=0, vmax=2,
               interpolation="nearest", aspect=reseau["aspect"])
    noir_interdit = reseau["vers_image"](masque_iridophores(L, params["bande"]))
    if noir_interdit.any():
        plt.imshow(np.ma.masked_where(~noir_interdit, noir_interdit), cmap=ListedColormap(["silver"]),
                   alpha=0.5, interpolation="nearest", aspect=reseau["aspect"])
    plt.title(f"L={L}, h={params['h']}, lX={params['lX']}, bande={params['bande']}, t = {t:.0f}")
    plt.axis("off")
    plt.savefig(os.path.join(dossier, "grille.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # 3. Historique : deux matrices binaires compactées et compressées
    if params["historique"]:
        os.makedirs(os.path.join(dossier, "historique"), exist_ok=True)
        np.savez_compressed(os.path.join(dossier, "historique", f"grille_{pas_faits:012d}.npz"),
                            jaune=np.packbits(grille == JAUNE, axis=-1),
                            noir=np.packbits(grille == NOIR, axis=-1))

    print(f"Sauvegarde : {pas_faits} / {params['n_pas']} pas, t = {t:.0f}", flush=True)


def charger_run(dossier):
    """
    Récupère un run sauvegardé par sauvegarder_run, pour le relancer.

    Le générateur est recréé puis remis dans l'état sauvegardé (JSON ayant
    changé les tuples de rng.getstate() en listes, ils sont reconvertis).

    Paramètres
    ----------
    dossier : str
        Dossier du run.

    Renvoie
    -------
    params : dict
        Paramètres relus dans params.json, sans pas_faits, t ni
        etat_generateur (renvoyés à part ou utilisés pour rng).
    grille : np.ndarray uint8, forme (L, L)
        Grille relue dans grille.npy.
    pas_faits : int
        Nombre de pas déjà effectués.
    rng : random.Random
        Générateur remis dans l'état sauvegardé.
    """
    with open(os.path.join(dossier, "params.json")) as f:
        params = json.load(f)
    grille = np.load(os.path.join(dossier, "grille.npy"))

    # JSON a transformé les tuples de l'état du générateur en listes
    version, etat, gauss = params.pop("etat_generateur")
    rng = random.Random()
    rng.setstate((version, tuple(etat), gauss))

    pas_faits = params.pop("pas_faits")
    params.pop("t")
    return params, grille, pas_faits, rng


def lambda_fft(grille, reseau):
    """
    Longueur d'onde dominante du motif, en cases, mesurée par FFT 2D du
    champ noir (grille == NOIR, moyenne retirée).

    Le pic du spectre a les indices signés (m_i, m_j) de np.fft.fftfreq.
    Sur le réseau carré, lambda = L / sqrt(m_i² + m_j²). Sur le réseau
    hexagonal la matrice est en coordonnées axiales et il faut repasser
    dans le plan : lambda = (√3/2) L / sqrt(m_i² - m_i m_j + m_j²), exact
    sur une onde plane.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        Grille de VIDE, JAUNE et NOIR.
    reseau : str
        "carre" ou "hex" (params["reseau"]).

    Renvoie
    -------
    lam : float
        Longueur d'onde, en cases.
    """
    L = grille.shape[0]
    f = (grille == NOIR).astype(float)
    F = np.abs(np.fft.fft2(f - f.mean())) ** 2
    F[0, 0] = 0
    ki, kj = np.unravel_index(np.argmax(F), F.shape)
    mi, mj = np.fft.fftfreq(L, 1 / L)[[ki, kj]]
    if reseau == "hex":
        return RACINE_3 / 2 * L / math.sqrt(mi * mi - mi * mj + mj * mj)
    return L / math.hypot(mi, mj)


if __name__ == "__main__":
    # 1. Paramètres : ligne de commande, et run récupéré s'il existe déjà
    args = lire_arguments()
    params, grille, pas_faits, rng = preparer_run(args)
    print(f"Run {params['sortie']} : {pas_faits} / {params['n_pas']} pas déjà faits", flush=True)

    if pas_faits >= params["n_pas"]:
        print("Run déjà terminé, rien à faire", flush=True)
    else:
        # 2. Simulation (sauvegardes régulières)
        grille = simulation_peau(grille, params, pas_faits, rng)

        # 3. Sauvegarde finale
        sauvegarder_run(grille, params, params["n_pas"], rng)
        print(f"Run terminé, résultats dans {params['sortie']}", flush=True)
