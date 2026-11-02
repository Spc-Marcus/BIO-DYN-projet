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
    python simulation.py --h 16 --lX 2.5 --historique   # garde toutes les grilles
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

# 2. Règles restantes, dans l'ordre des taux [BX, S, S, lX]
NAISSANCE_JAUNE = 0        # (2) S_i -> X_i
NOIR_TUE_JAUNE = 1         # (6) M_{i±1} + X_i -> M_{i±1} + S_i
JAUNE_TUE_NOIR = 2         # (7) X_{i±1} + M_i -> X_{i±1} + S_i
NOIR_NAIT_A_DISTANCE = 3   # (8) X_{i±h} + S_i -> X_{i±h} + M_i

# 3. Taux fixés comme dans l'article
BX = 1.0   # naissance spontanée du jaune
S = 1.0    # élimination entre voisins de couleurs différentes

# 4. Nombre de pas entre deux sauvegardes (pour pouvoir relancer).
#    Avec --historique, chaque sauvegarde est aussi gardée dans historique/
#    (pour une timeline ou un GIF avec timeline.py).
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
      --graine  graine du générateur aléatoire, entier (défaut 0) ;
      --historique  garde une copie de la grille à chaque sauvegarde
                (désactivé par défaut).

    Renvoie
    -------
    args : argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="Un run Monte Carlo (Bullara & De Decker 2015)")
    parser.add_argument("--L", type=int, default=100, help="taille de la grille")
    parser.add_argument("--h", type=int, default=16, help="portée de l'interaction à distance")
    parser.add_argument("--lX", type=float, default=2.5, help="naissance du noir à distance h")
    parser.add_argument("--bande", type=int, default=0, help="largeur de la bande d'iridophores")
    parser.add_argument("--n_pas", type=int, default=1000000000, help="nombre total de tentatives")
    parser.add_argument("--graine", type=int, default=0, help="graine du générateur aléatoire")
    parser.add_argument("--historique", action="store_true", help="garde la grille à chaque sauvegarde")
    return parser.parse_args()


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
        Clés "L", "h", "lX", "bande", "n_pas", "graine", "historique",
        "sortie".
    grille : np.ndarray d'entiers, forme (L, L)
        Grille de départ.
    pas_faits : int
        Nombre de pas déjà effectués.
    rng : random.Random
        Générateur aléatoire.
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
    }
    params["sortie"] = f"out/L{args.L}_h{args.h}_lX{args.lX}_b{args.bande}_g{args.graine}"

    # 2. Reprise si le run a déjà été sauvegardé
    if os.path.exists(os.path.join(params["sortie"], "params.json")):
        params_sauvegardes, grille, pas_faits, rng = charger_run(params["sortie"])
        return params, grille, pas_faits, rng

    # 3. Sinon nouveau run
    grille = initialiser_grille(args.L)
    rng = random.Random(args.graine)
    return params, grille, 0, rng


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
    return np.full((L, L), VIDE, dtype=np.uint8)


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
    noir_interdit = np.zeros((L, L), dtype=bool)
    if bande > 0:
        debut = (L - bande) // 2
        noir_interdit[debut:debut + bande, :] = True
    return noir_interdit


def tirer_regle(proba, rng):
    """
    Tire une règle au hasard selon les probabilités données.

    Tirage pondéré : on découpe [0, 1[ en segments de longueur proba[k] et
    on renvoie le segment dans lequel tombe un nombre uniforme u.

    Paramètres
    ----------
    proba : liste de 4 float
        Probabilités des règles, [BX, S, S, lX] / somme.
    rng : random.Random
        Générateur aléatoire.

    Renvoie
    -------
    regle : int
        Indice de la règle tirée, entre 0 et 3.
    """
    u = rng.random()
    cumul = 0.0
    for k in range(len(proba)):
        cumul = cumul + proba[k]
        if u < cumul:
            return k
    return len(proba) - 1   # u très proche de 1 : erreurs d'arrondi sur le cumul


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
    rng : random.Random
        Générateur aléatoire.

    Renvoie
    -------
    etat : int
        VIDE, JAUNE ou NOIR.
    """
    L = grille.shape[0]
    di, dj = ((1, 0), (-1, 0), (0, 1), (0, -1))[int(rng.random() * 4)]
    return grille[(i + di) % L, (j + dj) % L]


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
    rng : random.Random
        Générateur aléatoire.

    Renvoie
    -------
    reussie : bool
        True si la grille a été modifiée.
    """
    case = grille[i, j]

    if regle == NAISSANCE_JAUNE:
        # (2) S_i -> X_i
        if case == VIDE:
            grille[i, j] = JAUNE
            return True

    elif regle == NOIR_TUE_JAUNE:
        # (6) un noir voisin élimine le jaune
        if case == JAUNE and voisin_au_hasard(grille, i, j, rng) == NOIR:
            grille[i, j] = VIDE
            return True

    elif regle == JAUNE_TUE_NOIR:
        # (7) un jaune voisin élimine le noir
        if case == NOIR and voisin_au_hasard(grille, i, j, rng) == JAUNE:
            grille[i, j] = VIDE
            return True

    elif regle == NOIR_NAIT_A_DISTANCE:
        # (8) un jaune à distance h fait naître un noir, sauf sur les iridophores
        if case == VIDE and not noir_interdit[i, j] and case_a_distance(grille, i, j, h, rng) == JAUNE:
            grille[i, j] = NOIR
            return True

    # La configuration ne permet pas la règle : tentative ratée
    return False


def simulation_peau(grille, params, pas_faits, rng):
    """
    Boucle Monte Carlo : fait les pas restants, de pas_faits à n_pas.

    Le masque des iridophores est construit une fois avec
    masque_iridophores(L, bande). À chaque pas : tirer une case, tirer une
    règle, l'appliquer si possible. Tous les PAS_ENTRE_SAUVEGARDES pas, le
    run est sauvegardé avec sauvegarder_run pour pouvoir être relancé.

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
    Si la boucle est trop lente en Python pur, la compiler avec numba
    (@njit).
    """
    # 1. Paramètres du run
    L = params["L"]
    h = params["h"]
    noir_interdit = masque_iridophores(L, params["bande"])

    # 2. Probabilité de choisir chaque règle : P_k = taux_k / somme(taux)
    taux = [BX, S, S, params["lX"]]
    somme = sum(taux)
    proba = [t / somme for t in taux]

    # 3. Boucle Monte Carlo
    for pas in range(pas_faits, params["n_pas"]):
        # a. Tirage d'une case et d'une règle
        i = int(rng.random() * L)
        j = int(rng.random() * L)
        regle = tirer_regle(proba, rng)

        # b. Application de la règle si la configuration le permet
        appliquer_regle(grille, i, j, regle, h, noir_interdit, rng)

        # c. Sauvegarde régulière, pour pouvoir relancer (la finale est faite après la boucle)
        if (pas + 1) % PAS_ENTRE_SAUVEGARDES == 0 and pas + 1 < params["n_pas"]:
            sauvegarder_run(grille, params, pas + 1, rng)

    return grille


def sauvegarder_run(grille, params, pas_faits, rng):
    """
    Sauvegarde le run dans le dossier params["sortie"] (créé si besoin).

    Fichiers écrits :
      - grille.npy  : la grille ;
      - params.json : les paramètres, pas_faits, le temps t = pas_faits / L²
                      et l'état du générateur aléatoire ;
      - grille.png  : l'image de la grille (jaune, noir, vide en blanc),
                      avec la zone des iridophores marquée si bande > 0 ;
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

    # 2. Image de la grille, iridophores en gris argenté par-dessus
    plt.figure(figsize=(6, 6))
    plt.imshow(grille, cmap=ListedColormap(COULEURS), vmin=0, vmax=2, interpolation="nearest")
    noir_interdit = masque_iridophores(L, params["bande"])
    if noir_interdit.any():
        plt.imshow(np.ma.masked_where(~noir_interdit, noir_interdit),
                   cmap=ListedColormap(["silver"]), alpha=0.5, interpolation="nearest")
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
