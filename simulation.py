"""
Un run Monte Carlo du modèle de Bullara & De Decker (2015), lancé sur
GenOuest par slurm_pipeline/ (algorithme dans pseudocode.md).

Grille L x L à bords périodiques, chaque case VIDE, JAUNE (xantophore) ou
NOIR (mélanophore). À chaque pas, une case et une règle sont tirées au
hasard et la règle est appliquée si la configuration le permet. Les règles
sont dans un dictionnaire règle -> taux construit à partir des options ;
seules celles de taux non nul sont tirées. Par défaut, c'est le modèle de
l'article : règles (2), (6), (7) et (8). Les iridophores sont un masque de
cases où le noir ne peut pas naître (Fig. 1b-c).

Extensions, désactivées par défaut : blessure (--lB, --alpha, --tB),
naissance au contact (--contact, --rho0, --fX0), réseau hexagonal
(--reseau hex).

Le run est sauvegardé dans out/L{L}_h{h}_lX{lX}_b{bande}_g{graine}, plus un
suffixe par option différente de sa valeur par défaut. Relancer la même
commande reprend le run là où il s'est arrêté, avec exactement le même
résultat qu'un run fait d'une traite. Environ 1 µs par pas en Python pur.

Utilisation :
    python simulation.py --h 16 --lX 2.5
    python simulation.py --h 16 --lX 2.5 --historique                  # garde toutes les grilles
    python simulation.py --h 16 --lX 2.5 --reseau hex --contact --rho0 0.01
    python simulation.py --h 16 --lX 2.5 --tB 7000 --n_pas 800000000   # blessure au centre à t = 7000
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

# 2. Règles, numérotées comme les équations de l'article
NAISSANCE_JAUNE = 2        # (2) S_i -> X_i
NAISSANCE_NOIR = 3         # (3) S_i -> M_i
MORT_JAUNE = 4             # (4) X_i -> S_i
MORT_NOIR = 5              # (5) M_i -> S_i
NOIR_TUE_JAUNE = 6         # (6) M_{i±1} + X_i -> M_{i±1} + S_i
JAUNE_TUE_NOIR = 7         # (7) X_{i±1} + M_i -> X_{i±1} + S_i
NOIR_NAIT_A_DISTANCE = 8   # (8) X_{i±h} + S_i -> X_{i±h} + M_i
BLESSURE = 9               # (9) A_j -> S_j pour toute case j à distance <= R de i

# 3. Taux fixés comme dans l'article
BX = 1.0   # naissance spontanée du jaune
S = 1.0    # élimination entre voisins de couleurs différentes

# 4. Valeurs par défaut des options du modèle (celles de l'article, extensions désactivées)
DEFAUTS = {
    "bX": BX, "bM": 0.0, "dX": 0.0, "dM": 0.0, "s": S,   # taux des règles (2) à (7)
    "reseau": "carre",                                   # ou "hex"
    "contact": False, "rho0": 0.0, "fX0": 1.0,           # naissance au contact et départ
    "lB": 0.0, "alpha": 1.0, "tB": None,                 # blessure
}

# 5. Nombre de pas entre deux sauvegardes (pour pouvoir relancer)
PAS_ENTRE_SAUVEGARDES = 10**7

# 6. Voisins directs : réseau carré (ne pas changer l'ordre, il fixe les tirages),
#    réseau hexagonal = les mêmes plus une diagonale
VOISINS_CARRE = ((1, 0), (-1, 0), (0, 1), (0, -1))
VOISINS_HEX = VOISINS_CARRE + ((1, -1), (-1, 1))
RACINE_3 = math.sqrt(3)


def lire_arguments():
    """
    Lit les paramètres du run (voir python simulation.py --help) et vérifie
    les options : --contact demande rho0 > 0, une blessure demande
    1 <= alpha h < L/2, et la blessure unique doit avoir lieu avant n_pas.
    """
    parser = argparse.ArgumentParser(description="Un run Monte Carlo (Bullara & De Decker 2015)")

    # 1. Paramètres du run
    parser.add_argument("--L", type=int, default=100, help="taille de la grille")
    parser.add_argument("--h", type=int, default=16, help="portée de l'interaction à distance")
    parser.add_argument("--lX", type=float, default=2.5, help="(8) naissance du noir à distance h")
    parser.add_argument("--bande", type=int, default=0, help="largeur de la bande d'iridophores")
    parser.add_argument("--n_pas", type=int, default=1000000000, help="nombre total de tentatives")
    parser.add_argument("--graine", type=int, default=0, help="graine du générateur aléatoire")
    parser.add_argument("--historique", action="store_true", help="garde la grille à chaque sauvegarde")

    # 2. Options du modèle
    parser.add_argument("--bX", type=float, default=DEFAUTS["bX"], help="(2) naissance du jaune")
    parser.add_argument("--bM", type=float, default=DEFAUTS["bM"], help="(3) naissance du noir")
    parser.add_argument("--dX", type=float, default=DEFAUTS["dX"], help="(4) mort du jaune")
    parser.add_argument("--dM", type=float, default=DEFAUTS["dM"], help="(5) mort du noir")
    parser.add_argument("--s", type=float, default=DEFAUTS["s"], help="(6) et (7) élimination entre voisins")
    parser.add_argument("--reseau", choices=list(RESEAUX), default=DEFAUTS["reseau"], help="carré ou hexagonal")
    parser.add_argument("--contact", action="store_true", help="naissances seulement au contact d'une case occupée")
    parser.add_argument("--rho0", type=float, default=DEFAUTS["rho0"], help="fraction de cases occupées au départ")
    parser.add_argument("--fX0", type=float, default=DEFAUTS["fX0"], help="part de jaunes parmi ces cases")
    parser.add_argument("--lB", type=float, default=DEFAUTS["lB"], help="(9) blessure")
    parser.add_argument("--alpha", type=float, default=DEFAUTS["alpha"], help="rayon de la blessure R = alpha h")
    parser.add_argument("--tB", type=float, default=DEFAUTS["tB"], help="temps d'une blessure unique au centre")
    args = parser.parse_args()

    # 3. Vérifications
    if args.contact and args.rho0 == 0:
        parser.error("--contact demande --rho0 > 0 (sinon la peau vide le reste)")
    if (args.lB > 0 or args.tB is not None) and not 1 <= args.alpha * args.h < args.L / 2:
        parser.error("avec une blessure, il faut 1 <= alpha h < L/2")
    if args.tB is not None and round(args.tB * args.L**2) >= args.n_pas:
        parser.error("la blessure unique doit avoir lieu pendant le run : tB L² < n_pas")
    return args


def nom_dossier(params):
    """
    Dossier du run : out/L{L}_h{h}_lX{lX}_b{bande}_g{graine}, plus un suffixe
    par option différente de sa valeur par défaut (par exemple _dX0.1, _hex,
    _c, _r0.01). Sans option, c'est le nom des runs d'avant les options.
    """
    sortie = f"out/L{params['L']}_h{params['h']}_lX{params['lX']}_b{params['bande']}_g{params['graine']}"

    # 1. Taux
    for cle in ("bX", "bM", "dX", "dM", "s"):
        if params[cle] != DEFAUTS[cle]:
            sortie += f"_{cle}{params[cle]}"

    # 2. Réseau, contact et départ
    if params["reseau"] != DEFAUTS["reseau"]:
        sortie += f"_{params['reseau']}"
    if params["contact"]:
        sortie += "_c"
    if params["rho0"] != DEFAUTS["rho0"]:
        sortie += f"_r{params['rho0']}"
        if params["fX0"] != DEFAUTS["fX0"]:
            sortie += f"_fx{params['fX0']}"

    # 3. Blessure
    if params["lB"] != DEFAUTS["lB"]:
        sortie += f"_B{params['lB']}"
    if (params["lB"] > 0 or params["tB"] is not None) and params["alpha"] != DEFAUTS["alpha"]:
        sortie += f"_a{params['alpha']}"
    if params["tB"] is not None:
        sortie += f"_tB{params['tB']}"
    return sortie


def preparer_run(args):
    """
    Nouveau run, ou reprise depuis son dossier s'il contient un params.json.
    La reprise s'arrête si les paramètres sauvegardés diffèrent de la ligne
    de commande (sauf n_pas et historique, qui peuvent changer).
    Renvoie params, grille, pas_faits et le générateur.
    """
    # 1. Paramètres de la ligne de commande
    params = dict(vars(args))
    params["sortie"] = nom_dossier(params)

    # 2. Reprise si le run a déjà été sauvegardé (une option absente d'un ancien run vaut son défaut)
    if os.path.exists(os.path.join(params["sortie"], "params.json")):
        sauves, grille, pas_faits, rng = charger_run(params["sortie"])
        for cle in params:
            valeur = sauves.get(cle, DEFAUTS.get(cle))
            if cle not in ("n_pas", "historique", "sortie") and valeur != params[cle]:
                raise SystemExit(f"Erreur : {params['sortie']} a été lancé avec {cle} = {valeur}, pas {params[cle]}")
        return params, grille, pas_faits, rng

    # 3. Sinon nouveau run (le générateur d'abord : il tire les cases du départ)
    rng = random.Random(args.graine)
    grille = initialiser_grille(args.L, args.rho0, args.fX0, rng)
    return params, grille, 0, rng


def initialiser_grille(L, rho0, fX0, rng):
    """
    Grille de départ : vide, ou avec max(1, round(rho0 L²)) cases tirées sans
    remise, chacune JAUNE avec la probabilité fX0 et NOIR sinon.
    """
    grille = np.full((L, L), VIDE, dtype=np.uint8)
    if rho0 > 0:   # avec rho0 = 0, rien n'est tiré
        n0 = max(1, round(rho0 * L**2))
        for c in rng.sample(range(L**2), n0):
            grille[c // L, c % L] = JAUNE if rng.random() < fX0 else NOIR
    return grille


def masque_iridophores(L, bande):
    """
    Masque des iridophores : True sur les cases où le noir ne peut pas naître
    (règles (3) et (8)). Aujourd'hui une bande de `bande` lignes au milieu de
    la grille ; pour une autre forme, il suffit de changer cette fonction.
    """
    noir_interdit = np.zeros((L, L), dtype=bool)
    if bande > 0:
        debut = (L - bande) // 2
        noir_interdit[debut:debut + bande, :] = True
    return noir_interdit


def regles_actives(params):
    """
    Dictionnaire règle -> probabilité des règles de taux non nul, dans
    l'ordre (2) à (9). Pour ajouter une règle : son taux ici et sa fonction
    dans ACTIONS.
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

    # 2. Règles actives seulement (placée en dernier, une règle de taux nul
    #    serait sinon renvoyée par le repli d'arrondi de tirer_regle)
    taux = {regle: t for regle, t in taux.items() if t > 0}

    # 3. Probabilité de choisir chaque règle : P_k = taux_k / somme(taux)
    somme = sum(taux.values())
    return {regle: t / somme for regle, t in taux.items()}


def tirer_regle(proba, rng):
    """
    Tire une règle selon les probabilités : on découpe [0, 1[ en segments de
    longueur proba[regle] et on renvoie celui où tombe un nombre uniforme u.
    """
    u = rng.random()
    cumul = 0.0
    for regle, p in proba.items():
        cumul = cumul + p
        if u < cumul:
            return regle
    return regle   # u très proche de 1 : erreurs d'arrondi sur le cumul


def voisin_au_hasard(grille, i, j, voisins, rng):
    """État d'un voisin direct de (i, j) tiré au hasard, bords périodiques (règles 6 et 7)."""
    L = grille.shape[0]
    di, dj = voisins[int(rng.random() * len(voisins))]
    return grille[(i + di) % L, (j + dj) % L]


def contact_ok(grille, i, j, voisins):
    """Vrai si un voisin direct de (i, j) est occupé (naissances avec --contact). Ne tire rien."""
    L = grille.shape[0]
    for di, dj in voisins:
        if grille[(i + di) % L, (j + dj) % L] != VIDE:
            return True
    return False


def case_a_distance(grille, i, j, h, rng):
    """
    Réseau carré : état de la case (round(i + h cos alpha), round(j + h sin alpha)),
    avec un angle alpha tiré au hasard, bords périodiques (règle 8).
    """
    L = grille.shape[0]
    alpha = rng.random() * 2 * math.pi
    k = round(i + h * math.cos(alpha)) % L
    l = round(j + h * math.sin(alpha)) % L
    return grille[k, l]


def d2_carre(di, dj):
    """Réseau carré : distance² d'un décalage (di, dj)."""
    return di * di + dj * dj


def vers_image_carre(grille):
    """Réseau carré : l'image est la grille elle-même."""
    return grille


def d2_hex(di, dj):
    """Réseau hexagonal : distance² d'un décalage axial (di, dj) ; les 6 voisins sont à distance 1."""
    return di * di + di * dj + dj * dj


def arrondi_hex(dx, dy):
    """
    Réseau hexagonal : décalage axial (di, dj) de la case la plus proche du
    déplacement (dx, dy) dans le plan (arrondi « cube » de Red Blob Games).
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
    """Réseau hexagonal : comme case_a_distance, avec la case la plus proche du point visé."""
    L = grille.shape[0]
    alpha = rng.random() * 2 * math.pi
    di, dj = arrondi_hex(h * math.cos(alpha), h * math.sin(alpha))
    return grille[(i + di) % L, (j + dj) % L]


def vers_image_hex(grille):
    """
    Réseau hexagonal : image « briques » L x 2L, 2 pixels par case et la
    ligne i décalée de i/2 case, à afficher avec aspect = √3.
    """
    L = grille.shape[0]
    I = np.arange(L)[:, None]
    colonnes = (2 * np.arange(L)[None, :] + I) % (2 * L)
    image = np.empty((L, 2 * L), dtype=grille.dtype)
    image[I, colonnes] = grille
    image[I, (colonnes + 1) % (2 * L)] = grille
    return image


# 7. Réseaux : tout ce qui dépend de la géométrie de la peau. En hexagonal, on
#    garde la matrice L x L en coordonnées axiales : la case (i, j) est au point
#    x = j + i/2, y = (√3/2) i du plan.
RESEAUX = {
    "carre": {"voisins": VOISINS_CARRE, "d2": d2_carre, "case_a_distance": case_a_distance,
              "vers_image": vers_image_carre, "aspect": 1.0},
    "hex": {"voisins": VOISINS_HEX, "d2": d2_hex, "case_a_distance": case_a_distance_hex,
            "vers_image": vers_image_hex, "aspect": RACINE_3},
}


def disque(d2, R):
    """Décalages (DI, DJ) de toutes les cases à distance <= R de (0, 0), pour la blessure."""
    n = int(math.ceil(2 * R / RACINE_3)) + 1   # le disque tient dans |di|, |dj| <= n sur les deux réseaux
    decalages = [(di, dj) for di in range(-n, n + 1) for dj in range(-n, n + 1) if d2(di, dj) <= R * R]
    return np.array([d[0] for d in decalages]), np.array([d[1] for d in decalages])


def naissance_jaune(grille, i, j, modele, rng):
    """(2) S -> X si la case est vide ; (2') avec --contact, il faut aussi un voisin occupé."""
    if grille[i, j] == VIDE and (not modele["contact"] or contact_ok(grille, i, j, modele["voisins"])):
        grille[i, j] = JAUNE
        return True
    return False


def naissance_noir(grille, i, j, modele, rng):
    """(3) S -> M si la case est vide et hors iridophores ; (3') avec --contact."""
    if (grille[i, j] == VIDE and not modele["noir_interdit"][i, j]
            and (not modele["contact"] or contact_ok(grille, i, j, modele["voisins"]))):
        grille[i, j] = NOIR
        return True
    return False


def mort_jaune(grille, i, j, modele, rng):
    """(4) X -> S."""
    if grille[i, j] == JAUNE:
        grille[i, j] = VIDE
        return True
    return False


def mort_noir(grille, i, j, modele, rng):
    """(5) M -> S."""
    if grille[i, j] == NOIR:
        grille[i, j] = VIDE
        return True
    return False


def noir_tue_jaune(grille, i, j, modele, rng):
    """(6) un noir voisin élimine le jaune."""
    if grille[i, j] == JAUNE and voisin_au_hasard(grille, i, j, modele["voisins"], rng) == NOIR:
        grille[i, j] = VIDE
        return True
    return False


def jaune_tue_noir(grille, i, j, modele, rng):
    """(7) un jaune voisin élimine le noir."""
    if grille[i, j] == NOIR and voisin_au_hasard(grille, i, j, modele["voisins"], rng) == JAUNE:
        grille[i, j] = VIDE
        return True
    return False


def noir_nait_a_distance(grille, i, j, modele, rng):
    """
    (8) un jaune à distance h fait naître un noir, hors iridophores ; (8') avec
    --contact, le contact est testé avant de tirer l'angle.
    """
    if (grille[i, j] == VIDE and not modele["noir_interdit"][i, j]
            and (not modele["contact"] or contact_ok(grille, i, j, modele["voisins"]))
            and modele["case_a_distance"](grille, i, j, modele["h"], rng) == JAUNE):
        grille[i, j] = NOIR
        return True
    return False


def blessure(grille, i, j, modele, rng):
    """(9) la case (i, j) et toutes celles à distance <= R deviennent vides."""
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
    Applique la règle à la case (i, j) si la configuration le permet, avec
    sa fonction dans ACTIONS. Renvoie True si la grille a changé.
    """
    return ACTIONS[regle](grille, i, j, modele, rng)


def simulation_peau(grille, params, pas_faits, rng):
    """
    Boucle Monte Carlo de pas_faits à n_pas : à chaque pas, une case et une
    règle tirées au hasard. Sauvegarde tous les PAS_ENTRE_SAUVEGARDES pas.
    """
    # 1. Données lues par les règles
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

    # 3. Pas de la blessure unique (--tB), -1 sinon (jamais atteint)
    pas_blessure = round(params["tB"] * L**2) if params["tB"] is not None else -1

    # 4. Boucle Monte Carlo
    for pas in range(pas_faits, params["n_pas"]):
        # a. Blessure unique au centre, avant le tirage de ce pas. La boucle
        #    reprenant à pas_faits, elle n'est faite qu'une fois.
        if pas == pas_blessure:
            os.makedirs(params["sortie"], exist_ok=True)
            np.save(os.path.join(params["sortie"], "grille_avant_blessure.npy"), grille)
            blessure(grille, L // 2, L // 2, modele, rng)
            print(f"Blessure au centre : pas {pas}, t = {pas / L**2:.0f}", flush=True)

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
    Sauvegarde le run dans params["sortie"] : grille.npy, params.json (avec
    pas_faits, t = pas_faits / L² et l'état du générateur), grille.png et,
    avec --historique, historique/grille_<pas>.npz.
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

    # 2. Image de la grille sur son réseau, iridophores en gris argenté par-dessus
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

    # 3. Historique : deux matrices binaires (jaune, noir) compactées et compressées
    if params["historique"]:
        os.makedirs(os.path.join(dossier, "historique"), exist_ok=True)
        np.savez_compressed(os.path.join(dossier, "historique", f"grille_{pas_faits:012d}.npz"),
                            jaune=np.packbits(grille == JAUNE, axis=-1),
                            noir=np.packbits(grille == NOIR, axis=-1))

    print(f"Sauvegarde : {pas_faits} / {params['n_pas']} pas, t = {t:.0f}", flush=True)


def charger_run(dossier):
    """
    Relit un run sauvegardé par sauvegarder_run. Renvoie params (sans
    pas_faits, t ni etat_generateur), grille, pas_faits et le générateur
    remis dans son état.
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
    Longueur d'onde du motif en cases, donnée par le pic de la FFT 2D du champ
    noir. En hexagonal, la formule repasse des coordonnées axiales au plan.
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
