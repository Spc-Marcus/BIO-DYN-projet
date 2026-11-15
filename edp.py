"""
Équations de champ moyen du modèle, intégrées avec Euler ou RK2 (piste (C)
de projet.md).

Indépendant de Slurm : se lance en local, en une minute environ.

Au lieu de tirer des cases au hasard, on suit les densités x (jaune) et m
(noir) de chaque case, avec S = 1 - x - m, en supposant les cases
indépendantes (champ moyen, éq. (9)-(10) de l'article) :

    dx/dt = bX·S - dX·x - sM·x·<m>_voisins
    dm/dt = bM·S - dM·m - sX·m·<x>_voisins + lX·S·<x>_distance h

Deux façons de calculer les moyennes <.> :
  - noyau "exact"  : vraie moyenne sur les 4 voisins et sur le cercle de
                     rayon h (mêmes cases que case_a_distance du Monte
                     Carlo). Les densités restent entre 0 et 1 et le motif
                     se stabilise.
  - noyau "taylor" : développement à l'ordre 2 de l'article,
                     <u> ≈ u + c·∇²u, qui donne exactement les EDP
                     (14)-(15). Elles ne sont valables que tant que la
                     perturbation est petite : au-delà du seuil de Turing
                     les densités sortent de [0, 1] et la solution diverge,
                     l'intégration s'arrête alors d'elle-même.

La grille a un point par case (dx = a = 1), en dimension 1 ou 2, avec des
bords périodiques. Les taux bX, bM, dX, dM et s sont ceux de simulation.py.

Comme pour EquaRessort.py, Euler explicite devient instable si dt est trop
grand (voir pas_temps_critique_euler). Sur cette grille le mode le plus
contraignant a une valeur propre réelle (ω ≈ -3, soit dt < 0,7 environ) :
RK2 a alors la même limite qu'Euler, il n'est plus tolérant que pour des
modes oscillants. Son intérêt ici est sa précision (ordre 2).

L'état de départ est l'état mélangé (17) plus un petit bruit : il faut
donc lX >= bX.

Sorties, écrites dans results/ :
  - edp_<noyau>_L{L}_h{h}_lX{lX}.png : champs x et m à la fin ;
  - edp_<noyau>_L{L}_h{h}_lX{lX}.npz : champs x et m sauvegardés au cours du
                                       temps, et les instants t.

Utilisation :
    python edp.py --h 16 --lX 2.5
    python edp.py --h 16 --lX 2.5 --noyau taylor --T 20
    python edp.py --dim 1 --L 200 --h 9 --lX 1.0193 --methode euler --dt 0.1
"""

import argparse
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulation import BX, BM, DX, DM, S
from dispersion import etat_stationnaire, matrice_jacobienne

# 1. Dossier des figures et des champs sauvegardés
DOSSIER_RESULTATS = "results"

# 2. Nombre d'angles pour construire le cercle de rayon h (noyau exact, 2D)
N_ANGLES = 3600

# 3. Au-delà de cette valeur absolue, on considère que la solution diverge
SEUIL_DIVERGENCE = 10.0


def lire_arguments():
    """
    Lit les paramètres sur la ligne de commande.

    Options :
      --L        nombre de cases par côté, entier (défaut 100) ;
      --dim      dimension de la grille, 1 ou 2 (défaut 2) ;
      --h        portée de l'interaction à distance, entier (défaut 16) ;
      --lX       taux de naissance du noir à distance h (défaut 2.5) ;
      --noyau    "exact" ou "taylor" (défaut "exact") ;
      --methode  "euler" ou "rk2" (défaut "rk2") ;
      --dt       pas de temps (défaut 0.05) ;
      --T        temps total de la simulation (défaut 500) ;
      --graine   graine du bruit initial, entier (défaut 0).

    Renvoie
    -------
    args : argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="Équations de champ moyen (Bullara & De Decker 2015)")
    parser.add_argument("--L", type=int, default=100, help="nombre de cases par côté")
    parser.add_argument("--dim", type=int, default=2, choices=[1, 2], help="dimension de la grille")
    parser.add_argument("--h", type=int, default=16, help="portée de l'interaction à distance")
    parser.add_argument("--lX", type=float, default=2.5, help="naissance du noir à distance h")
    parser.add_argument("--noyau", default="exact", choices=["exact", "taylor"], help="calcul des moyennes")
    parser.add_argument("--methode", default="rk2", choices=["euler", "rk2"], help="schéma d'intégration")
    parser.add_argument("--dt", type=float, default=0.05, help="pas de temps")
    parser.add_argument("--T", type=float, default=500.0, help="temps total")
    parser.add_argument("--graine", type=int, default=0, help="graine du bruit initial")
    return parser.parse_args()


def laplacien(u):
    """
    Laplacien discret avec bords périodiques, pour dx = 1.

    En dimension 1 : u[i+1] + u[i-1] - 2 u[i]. En dimension 2 (5 points) :
    u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1] - 4 u[i, j].

    Paramètres
    ----------
    u : np.ndarray, forme (L,) ou (L, L)
        Champ à dériver.

    Renvoie
    -------
    lap : np.ndarray, même forme que u
    """
    lap = -2 * u.ndim * u
    for axe in range(u.ndim):
        lap = lap + np.roll(u, 1, axis=axe) + np.roll(u, -1, axis=axe)
    return lap


def construire_noyaux(forme, h):
    """
    Construit les deux noyaux de moyenne du noyau "exact", dans l'espace de
    Fourier (les moyennes se calculent ensuite par FFT, voir moyennes).

    Noyau des voisins : poids égaux sur les 2 voisins (dimension 1) ou les
    4 voisins (dimension 2). Noyau à distance h : poids 1/2 en +h et -h
    (dimension 1) ; en dimension 2, N_ANGLES angles régulièrement espacés
    sur le cercle, chacun arrondi à la case la plus proche comme dans
    case_a_distance de simulation.py.

    Paramètres
    ----------
    forme : tuple d'int
        Forme de la grille, (L,) ou (L, L).
    h : int
        Portée de l'interaction à distance, en cases.

    Renvoie
    -------
    (noyau_voisins, noyau_distance) : tuple de np.ndarray complexes
        Transformées de Fourier (np.fft.rfftn) des deux noyaux.
    """
    L = forme[0]
    voisins = np.zeros(forme)
    distance = np.zeros(forme)

    if len(forme) == 1:
        voisins[1 % L] += 0.5
        voisins[-1 % L] += 0.5
        distance[h % L] += 0.5
        distance[-h % L] += 0.5
    else:
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            voisins[di % L, dj % L] += 0.25
        for n in range(N_ANGLES):
            alpha = (n + 0.5) * 2 * math.pi / N_ANGLES
            distance[round(h * math.cos(alpha)) % L, round(h * math.sin(alpha)) % L] += 1 / N_ANGLES

    return np.fft.rfftn(voisins), np.fft.rfftn(distance)


def moyennes(u, h, noyau, noyaux):
    """
    Moyenne d'un champ sur les voisins directs et à distance h.

    Noyau "taylor" : <u> ≈ u + c·∇²u pour les voisins et u + c·h²·∇²u à
    distance h, avec c = 1/(2·dimension), soit 1/2 en dimension 1 (comme
    dans les éq. (11)-(13)) et 1/4 en dimension 2.
    Noyau "exact" : produit de convolution périodique avec les noyaux de
    construire_noyaux, calculé par FFT.

    Paramètres
    ----------
    u : np.ndarray, forme (L,) ou (L, L)
        Champ à moyenner.
    h : int
        Portée de l'interaction à distance.
    noyau : {"exact", "taylor"}
        Façon de calculer les moyennes.
    noyaux : tuple ou None
        Sortie de construire_noyaux (utilisée seulement par "exact").

    Renvoie
    -------
    (moyenne_voisins, moyenne_distance) : tuple de np.ndarray
        Même forme que u.
    """
    if noyau == "taylor":
        c = 1 / (2 * u.ndim)
        lap = laplacien(u)
        return u + c * lap, u + c * h**2 * lap

    # Les noyaux sont symétriques : la convolution est un produit en Fourier
    axes = tuple(range(u.ndim))
    u_fourier = np.fft.rfftn(u)
    moyenne_voisins = np.fft.irfftn(u_fourier * noyaux[0], s=u.shape, axes=axes)
    moyenne_distance = np.fft.irfftn(u_fourier * noyaux[1], s=u.shape, axes=axes)
    return moyenne_voisins, moyenne_distance


def second_membre(x, m, lX, h, noyau, noyaux):
    """
    Calcule les dérivées temporelles (dx/dt, dm/dt) du champ moyen.

        dx/dt = bX·S - dX·x - s·x·<m>_voisins
        dm/dt = bM·S - dM·m - s·m·<x>_voisins + lX·S·<x>_distance

    Avec le noyau "taylor", ce sont exactement les EDP (14)-(15).

    Paramètres
    ----------
    x, m : np.ndarray, forme (L,) ou (L, L)
        Densités de jaune et de noir.
    lX : float
        Taux de naissance du noir à distance h.
    h : int
        Portée de l'interaction à distance.
    noyau, noyaux :
        Voir moyennes.

    Renvoie
    -------
    (dxdt, dmdt) : tuple de np.ndarray, même forme que x
    """
    vide = 1 - x - m
    x_voisins, x_distance = moyennes(x, h, noyau, noyaux)
    m_voisins, m_distance = moyennes(m, h, noyau, noyaux)

    dxdt = BX * vide - DX * x - S * x * m_voisins
    dmdt = BM * vide - DM * m - S * m * x_voisins + lX * vide * x_distance
    return dxdt, dmdt


def pas_euler(x, m, dt, lX, h, noyau, noyaux):
    """
    Avance l'état d'un pas de temps avec la méthode d'Euler explicite.

    Valeur_suivante = Valeur_actuelle + Dérivée * dt.

    Paramètres
    ----------
    x, m : np.ndarray
        État à l'instant t.
    dt : float
        Pas de temps.
    lX, h, noyau, noyaux :
        Voir second_membre.

    Renvoie
    -------
    (x_suivant, m_suivant) : tuple de np.ndarray
    """
    dxdt, dmdt = second_membre(x, m, lX, h, noyau, noyaux)
    return x + dxdt * dt, m + dmdt * dt


def pas_rk2(x, m, dt, lX, h, noyau, noyaux):
    """
    Avance l'état d'un pas de temps avec Runge-Kutta 2 (point milieu).

    Pentes k1 au début de l'intervalle, état estimé au point milieu
    (t + dt/2), pentes k2 au point milieu, puis mise à jour avec k2.

    Paramètres
    ----------
    x, m : np.ndarray
        État à l'instant t.
    dt : float
        Pas de temps.
    lX, h, noyau, noyaux :
        Voir second_membre.

    Renvoie
    -------
    (x_suivant, m_suivant) : tuple de np.ndarray
    """
    # a. Pentes au début de l'intervalle (k1)
    k1_x, k1_m = second_membre(x, m, lX, h, noyau, noyaux)

    # b. Estimation de l'état au point milieu (t + dt/2)
    x_mid = x + 0.5 * dt * k1_x
    m_mid = m + 0.5 * dt * k1_m

    # c. Pentes au point milieu (k2)
    k2_x, k2_m = second_membre(x_mid, m_mid, lX, h, noyau, noyaux)

    # d. Mise à jour
    return x + dt * k2_x, m + dt * k2_m


def condition_initiale(forme, lX, amplitude=0.001, graine=0):
    """
    État de départ : équilibre mélangé (17) plus un petit bruit gaussien.

    Le bruit sert à déclencher l'instabilité de Turing.

    Paramètres
    ----------
    forme : tuple d'int
        Forme de la grille, (L,) ou (L, L).
    lX : float
        Taux de naissance du noir à distance h (lX >= bX).
    amplitude : float, optionnel
        Écart-type du bruit ajouté à x et à m.
    graine : int, optionnel
        Graine du générateur aléatoire.

    Renvoie
    -------
    (x0, m0) : tuple de np.ndarray, de forme `forme`
    """
    rng = np.random.default_rng(graine)
    x2, m2 = etat_stationnaire(BX, lX, S)
    x0 = x2 + amplitude * rng.standard_normal(forme)
    m0 = m2 + amplitude * rng.standard_normal(forme)
    return x0, m0


def pas_temps_critique_euler(lX, h, dim, noyau):
    """
    Estime le plus grand dt pour lequel Euler explicite reste stable.

    Un mode de valeur propre ω est amorti par Euler si |1 + ω·dt| < 1, soit
    dt < -2·Re(ω) / |ω|². On prend le minimum sur les modes amortis
    (Re(ω) < 0) de la grille, avec les valeurs propres de la matrice J(k) de
    dispersion.py. Les modes amplifiés (Re(ω) > 0) ne comptent pas : leur
    croissance est celle du motif.

    Paramètres
    ----------
    lX : float
        Taux de naissance du noir à distance h.
    h : int
        Portée de l'interaction à distance.
    dim : {1, 2}
        Dimension de la grille.
    noyau : {"exact", "taylor"}
        Façon de calculer les moyennes.

    Renvoie
    -------
    dt_max : float

    Notes
    -----
    Pour le noyau "taylor", le laplacien discret a des valeurs propres
    entre 0 et -4·dim : on balaie donc k² de 0 à 4·dim. Pour le noyau
    "exact", k va de 0 à π le long d'un axe : c'est une estimation.
    """
    if noyau == "taylor":
        k = np.sqrt(np.linspace(0, 4 * dim, 2000))
    else:
        k = np.linspace(0, np.pi, 2000)

    omega = np.linalg.eigvals(matrice_jacobienne(k, BX, lX, S, h, noyau=noyau, dim=dim)).ravel()
    amortis = omega[omega.real < 0]
    return float((-2 * amortis.real / np.abs(amortis) ** 2).min())


def integrer(x0, m0, lX, h, dt, T_max, noyau="exact", methode="rk2", n_sauvegardes=50):
    """
    Intègre les équations de t = 0 à T_max.

    L'intégration s'arrête avant T_max si la solution diverge (valeur non
    finie ou supérieure à SEUIL_DIVERGENCE en valeur absolue), ce qui arrive
    avec le noyau "taylor" au-delà du seuil de Turing, ou avec un dt trop
    grand.

    Paramètres
    ----------
    x0, m0 : np.ndarray, forme (L,) ou (L, L)
        Conditions initiales.
    lX : float
        Taux de naissance du noir à distance h.
    h : int
        Portée de l'interaction à distance.
    dt : float
        Pas de temps.
    T_max : float
        Temps total de la simulation.
    noyau : {"exact", "taylor"}, optionnel
        Façon de calculer les moyennes.
    methode : {"euler", "rk2"}, optionnel
        Schéma d'intégration.
    n_sauvegardes : int, optionnel
        Nombre d'instantanés conservés, régulièrement espacés (plus l'état
        initial).

    Renvoie
    -------
    t : np.ndarray, forme (n,)
        Instants sauvegardés, de 0 au dernier instant atteint.
    x_hist, m_hist : np.ndarray, forme (n, L) ou (n, L, L)
        Champs aux instants sauvegardés.
    diverge : bool
        True si l'intégration s'est arrêtée sur une divergence.
    """
    # 1. Préparation
    noyaux = construire_noyaux(x0.shape, h) if noyau == "exact" else None
    pas = pas_euler if methode == "euler" else pas_rk2
    n = int(round(T_max / dt))
    intervalle = max(1, n // n_sauvegardes)

    # 2. Initialisation des tableaux de stockage
    t = [0.0]
    x_hist = [x0]
    m_hist = [m0]
    x, m = x0, m0
    diverge = False

    # 3. Boucle d'intégration
    for i in range(n):
        x, m = pas(x, m, dt, lX, h, noyau, noyaux)

        # a. Arrêt si la solution diverge
        if not (np.abs(x).max() < SEUIL_DIVERGENCE and np.abs(m).max() < SEUIL_DIVERGENCE):
            diverge = True
            print(f"Divergence à t = {(i + 1) * dt:.2f} : intégration arrêtée")
            break

        # b. Enregistrement
        if (i + 1) % intervalle == 0 or i == n - 1:
            t.append((i + 1) * dt)
            x_hist.append(x)
            m_hist.append(m)

    return np.array(t), np.array(x_hist), np.array(m_hist), diverge


def afficher_champs(x, m, titre, chemin):
    """
    Affiche les champs x (jaune) et m (noir) et sauvegarde la figure.

    En dimension 2 : deux images côte à côte, de 0 (blanc) à 1. En
    dimension 1 : les deux courbes en fonction de la position.

    Paramètres
    ----------
    x, m : np.ndarray, forme (L,) ou (L, L)
        Densités de jaune et de noir.
    titre : str
        Titre de la figure.
    chemin : str
        Fichier PNG à écrire.
    """
    if x.ndim == 1:
        plt.figure(figsize=(10, 6))
        plt.plot(x, label="Jaune (x)", color="orange")
        plt.plot(m, label="Noir (m)", color="black")
        plt.xlabel("Position (cases)")
        plt.ylabel("Densité")
        plt.title(titre)
        plt.legend()
        plt.grid(True)
    else:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        for ax, champ, nom, couleurs in zip(axes, (x, m), ("Jaune (x)", "Noir (m)"), ("YlOrBr", "Greys")):
            image = ax.imshow(champ, cmap=couleurs, vmin=0, vmax=1, interpolation="nearest")
            ax.set_title(nom)
            ax.axis("off")
            fig.colorbar(image, ax=ax, fraction=0.046)
        fig.suptitle(titre)

    plt.savefig(chemin, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Graphique sauvegardé avec succès dans : {chemin}")


if __name__ == "__main__":
    # 1. Paramètres
    args = lire_arguments()
    if args.lX < BX:
        raise SystemExit(f"lX = {args.lX} < bX = {BX} : l'état mélangé (17) n'existe pas, la peau reste jaune")
    forme = (args.L,) * args.dim
    os.makedirs(DOSSIER_RESULTATS, exist_ok=True)

    # 2. Conditions initiales : état (17) + bruit
    x0, m0 = condition_initiale(forme, args.lX, graine=args.graine)

    # 3. Pas de temps : rappel de la limite de stabilité d'Euler
    dt_max = pas_temps_critique_euler(args.lX, args.h, args.dim, args.noyau)
    print(f"dt = {args.dt} ({args.methode}) ; Euler explicite est stable si dt < {dt_max:.3f}")

    # 4. Intégration
    t, x_hist, m_hist, diverge = integrer(x0, m0, args.lX, args.h, args.dt, args.T,
                                          noyau=args.noyau, methode=args.methode)
    print(f"Fin à t = {t[-1]:.1f} : x entre {x_hist[-1].min():.3f} et {x_hist[-1].max():.3f}, "
          f"m entre {m_hist[-1].min():.3f} et {m_hist[-1].max():.3f}")

    # 5. Sauvegarde des champs et de la figure sous results/
    nom = f"edp_{args.noyau}_L{args.L}_h{args.h}_lX{args.lX}"
    if args.dim == 1:
        nom = nom + "_1d"
    np.savez_compressed(os.path.join(DOSSIER_RESULTATS, nom + ".npz"), t=t, x=x_hist, m=m_hist)
    afficher_champs(x_hist[-1], m_hist[-1],
                    f"Champ moyen, noyau {args.noyau}, {args.methode} : h = {args.h}, lX = {args.lX}, t = {t[-1]:.0f}",
                    os.path.join(DOSSIER_RESULTATS, nom + ".png"))
