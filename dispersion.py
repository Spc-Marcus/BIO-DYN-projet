"""
Analyse de stabilité linéaire du modèle (sections 5 à 7 de projet.md).

Indépendant de Slurm : se lance en local, en quelques secondes.

On linéarise les équations de champ moyen autour de l'état mélangé (17),
dans la limite de l'article (bM = dX = dM = 0, sX = sM = s), avec des
perturbations proportionnelles à exp(ikr + ωt). Le mode k est amplifié si
ω(k) > 0 : c'est l'instabilité de Turing.

Deux noyaux pour les interactions non locales :
  - "taylor" : développement à l'ordre 2 de l'article, qui donne les EDP
               (14)-(15) et les formules (19) et (20) ;
  - "exact"  : moyenne exacte sur les voisins et sur la distance h, sans
               développement (section 7 de projet.md).
Et deux dimensions : 1 (calcul de l'article) ou 2 (comme les simulations).

Sorties, écrites dans results/ :
  - dispersion_fig3.png : ω(k) pour h = 7, 8, 9 (Fig. 3 de l'article) ;
  - diagramme_fig2.png  : seuil h_T et λ/2 en fonction de lX (Fig. 2) ;
  - seuils.csv          : h_T et λ pour quelques lX, noyaux Taylor et exact.

Utilisation :
    python dispersion.py
"""

import csv
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulation import BX, S

# 1. Dossier des figures et tableaux
DOSSIER_RESULTATS = "results"

# 2. Paramètres de la Fig. 3 (bX = 1 est absent de la légende de l'article)
LX_FIG3 = 1.0193
H_FIG3 = [7, 8, 9]   # stable, seuil, motif

# 3. Valeurs de lX du tableau des seuils (projet.md, sections 6 et 7)
LX_TABLEAU = [1.0193, 1.5, 2.5, 3.5]

# 4. Nombres d'onde explorés pour chercher le mode le plus instable :
#    de 0 à π/a (plus petite longueur d'onde représentable : 2 cases)
N_K = 2000


def etat_stationnaire(bX, lX, s):
    """
    Calcule l'état stationnaire mélangé (x2, m2) de l'éq. (17).

    x2 = bX / lX,    m2 = (lX - bX) / (lX + s).

    Paramètres
    ----------
    bX : float
        Taux de naissance spontanée du jaune.
    lX : float
        Taux de naissance du noir à distance h.
    s : float
        Taux d'élimination entre voisins.

    Renvoie
    -------
    (x2, m2) : tuple de float
        Densités de jaune et de noir à l'équilibre.

    Lève
    ----
    ValueError
        Si lX < bX (l'état n'existe pas : m2 serait négatif).
    """
    if lX < bX:
        raise ValueError(f"L'état mélangé (17) n'existe que si lX >= bX (ici lX = {lX}, bX = {bX})")
    return bX / lX, (lX - bX) / (lX + s)


def facteurs_noyau(k, h, a=1.0, noyau="taylor", dim=1):
    """
    Facteurs de Fourier des deux moyennes non locales, pour le mode k.

    Une perturbation exp(ikr) moyennée sur les voisins directs est
    multipliée par A_n(k), et moyennée à distance h par A_h(k) :

      noyau     dim   A_n                    A_h
      taylor    1     1 - (ka)²/2            1 - (kha)²/2
      taylor    2     1 - (ka)²/4            1 - (kha)²/4
      exact     1     cos(ka)                cos(kha)
      exact     2     (1 + cos(ka)) / 2      J0(kha)

    En dimension 2, k est pris le long d'un axe de la grille, et J0
    (fonction de Bessel) est la moyenne de cos(kha cos θ) sur les angles θ,
    calculée numériquement. L'arrondi de la case à distance h sur la grille
    est négligé.

    Paramètres
    ----------
    k : float ou np.ndarray
        Nombre(s) d'onde.
    h : float
        Portée de l'interaction à distance, en cases.
    a : float, optionnel
        Taille d'une case.
    noyau : {"taylor", "exact"}, optionnel
        Développement de l'article ou moyenne exacte.
    dim : {1, 2}, optionnel
        Dimension de l'espace.

    Renvoie
    -------
    (A_n, A_h) : tuple de np.ndarray
        Même forme que k.
    """
    k = np.asarray(k, dtype=float)

    if noyau == "taylor":
        A_n = 1 - (k * a) ** 2 / (2 * dim)
        A_h = 1 - (k * h * a) ** 2 / (2 * dim)
    elif noyau == "exact" and dim == 1:
        A_n = np.cos(k * a)
        A_h = np.cos(k * h * a)
    elif noyau == "exact" and dim == 2:
        A_n = (1 + np.cos(k * a)) / 2
        # J0(kha) : moyenne de cos(kha cos θ) sur un tour
        theta = np.linspace(0, 2 * np.pi, 720, endpoint=False)
        A_h = np.cos(np.multiply.outer(k * h * a, np.cos(theta))).mean(axis=-1)
    else:
        raise ValueError(f"noyau = {noyau!r} et dim = {dim!r} : attendu 'taylor' ou 'exact', et 1 ou 2")

    return A_n, A_h


def matrice_jacobienne(k, bX, lX, s, h, a=1.0, noyau="taylor", dim=1):
    """
    Matrice J(k) du système linéarisé autour de l'état (17).

    Avec A_n et A_h les facteurs de facteurs_noyau et S2 = 1 - x2 - m2 :

        J11 = -bX - s·m2              J12 = -bX - s·x2·A_n
        J21 = -lX·x2 - s·m2·A_n + lX·S2·A_h
        J22 = -(s + lX)·x2

    Avec le noyau de Taylor en dimension 1, on retrouve J(k) = J0 + q·J1 de
    projet.md (section 6), q = a²k².

    Paramètres
    ----------
    k : float ou np.ndarray
        Nombre(s) d'onde de la perturbation.
    bX, lX, s : float
        Paramètres du modèle.
    h : float
        Portée de l'interaction à distance (en cases).
    a : float, optionnel
        Taille d'une case.
    noyau, dim : optionnels
        Voir facteurs_noyau.

    Renvoie
    -------
    J : np.ndarray, forme (..., 2, 2)
        Une matrice 2 x 2 par valeur de k.
    """
    x2, m2 = etat_stationnaire(bX, lX, s)
    S2 = 1 - x2 - m2
    A_n, A_h = facteurs_noyau(k, h, a, noyau, dim)

    J = np.empty(A_n.shape + (2, 2))
    J[..., 0, 0] = -bX - s * m2
    J[..., 0, 1] = -bX - s * x2 * A_n
    J[..., 1, 0] = -lX * x2 - s * m2 * A_n + lX * S2 * A_h
    J[..., 1, 1] = -(s + lX) * x2
    return J


def taux_croissance(k, bX, lX, s, h, a=1.0, noyau="taylor", dim=1):
    """
    Relation de dispersion ω(k) : plus grande partie réelle des valeurs
    propres de J(k).

    ω(k) > 0 signifie que le mode k est amplifié (motif de Turing).

    Paramètres
    ----------
    k : float ou np.ndarray
        Nombre(s) d'onde.
    bX, lX, s, h, a, noyau, dim :
        Voir matrice_jacobienne.

    Renvoie
    -------
    omega : np.ndarray
        Même forme que k.
    """
    J = matrice_jacobienne(k, bX, lX, s, h, a, noyau, dim)
    return np.linalg.eigvals(J).real.max(axis=-1)


def coefficients_determinant(bX, lX, s, h):
    """
    Coefficients du déterminant det J(q) = c0 + c1·q + c2·q², q = a²k²
    (noyau de Taylor, dimension 1).

    c0 = bX·s·(lX - bX) / lX
    c1 = bX·s·[lX - h²(lX - bX)] / (2 lX)
    c2 = bX·s²·(lX - bX)(h² - 1) / [4 lX (lX + s)]

    c0 > 0 : stable à k = 0. c2 > 0 si h > 1 : stable à grand k. Il y a
    instabilité si c1 est assez négatif : c1² > 4·c0·c2.

    Paramètres
    ----------
    bX, lX, s : float
        Paramètres du modèle.
    h : float
        Portée de l'interaction à distance.

    Renvoie
    -------
    (c0, c1, c2) : tuple de float
    """
    c0 = bX * s * (lX - bX) / lX
    c1 = bX * s * (lX - h**2 * (lX - bX)) / (2 * lX)
    c2 = bX * s**2 * (lX - bX) * (h**2 - 1) / (4 * lX * (lX + s))
    return c0, c1, c2


def h_seuil_turing(bX, lX, s):
    """
    Portée critique h_T au-delà de laquelle un motif apparaît, éq. (19).

    Obtenue en annulant le discriminant du déterminant : c1² = 4·c0·c2.
    Elle est la même en dimension 1 et 2.

    Paramètres
    ----------
    bX, lX, s : float
        Paramètres du modèle (lX > bX).

    Renvoie
    -------
    h_T : float
    """
    racine = math.sqrt(lX * (s + bX) * s * (lX - bX))
    numerateur = 2 * s * (lX - bX) + lX * (s + lX) + 2 * racine
    return math.sqrt(numerateur / ((lX - bX) * (s + lX)))


def longueur_onde_turing(bX, lX, s):
    """
    Longueur d'onde du motif au seuil h = h_T, éq. (20), en nombre de cases.

    Formule de l'article, en dimension 1. En dimension 2 la longueur d'onde
    de Taylor est divisée par racine de 2.

    Paramètres
    ----------
    bX, lX, s : float
        Paramètres du modèle (lX > bX).

    Renvoie
    -------
    lam : float
        Longueur d'onde au seuil.
    """
    racine = math.sqrt(lX * (s + bX) * s * (lX - bX))
    return math.sqrt(2) * math.pi / math.sqrt((racine - s * (lX - bX)) / (s * bX))


def mode_le_plus_instable(bX, lX, s, h, a=1.0, noyau="taylor", dim=1):
    """
    Cherche le nombre d'onde qui maximise ω(k), pour k > 0.

    ω est évalué sur N_K valeurs de k entre 0 (exclu) et π/a, puis le
    maximum est affiné par une parabole passant par ses deux voisins.

    Paramètres
    ----------
    bX, lX, s, h, a, noyau, dim :
        Voir matrice_jacobienne.

    Renvoie
    -------
    k_max : float
        Nombre d'onde le plus amplifié.
    omega_max : float
        Taux de croissance correspondant (positif s'il y a un motif).
    """
    k = np.linspace(0, np.pi / a, N_K + 1)[1:]
    omega = taux_croissance(k, bX, lX, s, h, a, noyau, dim)
    i = int(omega.argmax())

    # Affinage : sommet de la parabole passant par les points i-1, i, i+1
    if 0 < i < len(k) - 1:
        y0, y1, y2 = omega[i - 1], omega[i], omega[i + 1]
        courbure = y0 - 2 * y1 + y2
        if courbure < 0:
            decalage = 0.5 * (y0 - y2) / courbure
            return k[i] + decalage * (k[1] - k[0]), y1 - 0.25 * (y0 - y2) * decalage

    return k[i], omega[i]


def longueur_onde_instable(bX, lX, s, h, a=1.0, noyau="taylor", dim=1):
    """
    Longueur d'onde du mode le plus instable, λ = 2π / k_max, en unités de a.

    Contrairement à (20), valable au seuil seulement, elle est calculée pour
    la portée h donnée : c'est elle qu'on compare aux motifs Monte Carlo.

    Paramètres
    ----------
    bX, lX, s, h, a, noyau, dim :
        Voir matrice_jacobienne.

    Renvoie
    -------
    lam : float
        Longueur d'onde, ou nan si aucun mode n'est amplifié.
    """
    k_max, omega_max = mode_le_plus_instable(bX, lX, s, h, a, noyau, dim)
    if omega_max <= 0:
        return float("nan")
    return 2 * math.pi / (k_max * a)


def h_seuil_numerique(bX, lX, s, a=1.0, noyau="exact", dim=1, h_max=60):
    """
    Plus petit h entier pour lequel un mode est amplifié (max_k ω(k) > 0).

    Sert à vérifier (19) et à obtenir le seuil avec le noyau exact, pour
    lequel il n'y a pas de formule.

    Paramètres
    ----------
    bX, lX, s, a, noyau, dim :
        Voir matrice_jacobienne.
    h_max : int, optionnel
        Plus grande valeur de h testée.

    Renvoie
    -------
    h_T : int ou None
        Seuil trouvé, ou None s'il n'y a pas d'instabilité jusqu'à h_max.
    """
    for h in range(2, h_max + 1):
        k_max, omega_max = mode_le_plus_instable(bX, lX, s, h, a, noyau, dim)
        if omega_max > 0:
            return h
    return None


def tracer_dispersion(h_valeurs, bX, lX, s, chemin, a=1.0, noyau="taylor", dim=1):
    """
    Trace ω(k) pour plusieurs portées h (reproduction de la Fig. 3).

    Paramètres
    ----------
    h_valeurs : liste de float
        Portées à comparer (Fig. 3 : 7, 8, 9).
    bX, lX, s : float
        Paramètres du modèle (Fig. 3 : bX = s = 1, lX = 1.0193).
    chemin : str
        Fichier PNG à écrire.
    a, noyau, dim : optionnels
        Voir matrice_jacobienne.
    """
    k = np.linspace(0, 1.4, 500)

    plt.figure(figsize=(8, 5))
    for h, couleur in zip(h_valeurs, ["red", "blue", "orange", "green", "purple"]):
        omega = taux_croissance(k, bX, lX, s, h, a, noyau, dim)
        plt.plot(k, omega, label=f"h = {h}", color=couleur, linewidth=2)
    plt.axhline(0, color="black", linewidth=0.8)

    plt.xlabel("Nombre d'onde k")
    plt.ylabel("Taux de croissance ω(k)")
    plt.title(f"Relation de dispersion (lX = {lX}, bX = {bX}, s = {s}, noyau {noyau}, dimension {dim})")
    plt.ylim(-0.03, 0.01)
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.savefig(chemin, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Graphique sauvegardé avec succès dans : {chemin}")


def tracer_diagramme_phase(bX, s, chemin):
    """
    Trace h_T(lX) (19) et λ(lX)/2 (20) dans le plan (lX, h), comme la Fig. 2,
    avec en plus le seuil du noyau exact (dimensions 1 et 2).

    Au-dessus de la courbe h_T : motifs ; en dessous : tout jaune.

    Paramètres
    ----------
    bX, s : float
        Paramètres du modèle.
    chemin : str
        Fichier PNG à écrire.
    """
    # 1. Formules de l'article (noyau de Taylor)
    lX_valeurs = np.linspace(bX * 1.01, 4.0, 300)
    h_T = [h_seuil_turing(bX, lX, s) for lX in lX_valeurs]
    demi_lambda = [longueur_onde_turing(bX, lX, s) / 2 for lX in lX_valeurs]

    # 2. Seuil numérique du noyau exact, sur moins de points (plus lent)
    lX_exact = np.linspace(bX * 1.02, 4.0, 40)
    h_exact_1d = [h_seuil_numerique(bX, lX, s, noyau="exact", dim=1) for lX in lX_exact]
    h_exact_2d = [h_seuil_numerique(bX, lX, s, noyau="exact", dim=2) for lX in lX_exact]

    plt.figure(figsize=(8, 5))
    plt.plot(lX_valeurs, h_T, label="h_T, éq. (19)", color="blue", linewidth=2)
    plt.plot(lX_valeurs, demi_lambda, label="λ/2, éq. (20)", color="red", linewidth=2)
    plt.step(lX_exact, h_exact_1d, label="h_T noyau exact, dimension 1", color="green", where="mid")
    plt.step(lX_exact, h_exact_2d, label="h_T noyau exact, dimension 2", color="purple", where="mid")

    plt.xlabel("lX")
    plt.ylabel("h (en cases)")
    plt.title(f"Seuil de Turing (bX = {bX}, s = {s}) : motifs au-dessus de h_T")
    plt.ylim(0, 20)
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.savefig(chemin, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Graphique sauvegardé avec succès dans : {chemin}")


def ecrire_tableau_seuils(lX_valeurs, bX, s, chemin):
    """
    Écrit et affiche le tableau des seuils pour quelques valeurs de lX.

    Colonnes : lX, h_T (19), λ (20), λ/h_T, puis pour le noyau exact en
    dimension 1 et 2 : le seuil entier h_T et λ/h à ce seuil.

    Paramètres
    ----------
    lX_valeurs : liste de float
        Valeurs de lX (toutes > bX).
    bX, s : float
        Paramètres du modèle.
    chemin : str
        Fichier CSV à écrire.
    """
    colonnes = ["lX", "hT_19", "lambda_20", "lambda_sur_hT",
                "hT_exact_1d", "lambda_sur_h_exact_1d", "hT_exact_2d", "lambda_sur_h_exact_2d"]
    lignes = []
    for lX in lX_valeurs:
        h_T = h_seuil_turing(bX, lX, s)
        lam = longueur_onde_turing(bX, lX, s)
        ligne = [lX, round(h_T, 2), round(lam, 2), round(lam / h_T, 2)]
        for dim in (1, 2):
            h_exact = h_seuil_numerique(bX, lX, s, noyau="exact", dim=dim)
            lam_exact = longueur_onde_instable(bX, lX, s, h_exact, noyau="exact", dim=dim)
            ligne += [h_exact, round(lam_exact / h_exact, 2)]
        lignes.append(ligne)

    with open(chemin, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(colonnes)
        writer.writerows(lignes)

    print("   ".join(colonnes))
    for ligne in lignes:
        print("   ".join(str(v) for v in ligne))
    print(f"Tableau sauvegardé dans : {chemin}")


if __name__ == "__main__":
    # 1. Création du dossier de sortie s'il n'existe pas
    os.makedirs(DOSSIER_RESULTATS, exist_ok=True)

    # 2. Relation de dispersion (Fig. 3), noyau de Taylor puis noyau exact
    tracer_dispersion(H_FIG3, BX, LX_FIG3, S, os.path.join(DOSSIER_RESULTATS, "dispersion_fig3.png"))
    tracer_dispersion([14, 15, 16], BX, LX_FIG3, S,
                      os.path.join(DOSSIER_RESULTATS, "dispersion_exact.png"), noyau="exact")

    # 3. Diagramme de phase (Fig. 2)
    tracer_diagramme_phase(BX, S, os.path.join(DOSSIER_RESULTATS, "diagramme_fig2.png"))

    # 4. Tableau des seuils
    ecrire_tableau_seuils(LX_TABLEAU, BX, S, os.path.join(DOSSIER_RESULTATS, "seuils.csv"))
