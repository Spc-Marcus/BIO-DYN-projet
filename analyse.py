"""
Analyse des motifs Monte Carlo et comparaison avec la théorie (piste (D)
de projet.md).

Indépendant de Slurm : se lance en local après avoir rapatrié out/, ou sur
la frontale.

Pour chaque run de out/ (dossier contenant params.json et grille.npy,
écrits par simulation.py), on mesure sur la grille finale :
  - les fractions de cases vides, jaunes et noires ;
  - la longueur d'onde du motif, par transformée de Fourier 2D ;
  - le type de motif : tout jaune, taches noires, rayures, taches jaunes
    ou tout noir.
Puis on compare la longueur d'onde mesurée à celles de dispersion.py :
éq. (20) de l'article (noyau de Taylor, au seuil) et mode le plus instable
du noyau exact, en dimension 1 et 2.

Sorties, écrites dans results/ :
  - analyse.csv        : une ligne par run ;
  - analyse_motifs.png : les grilles finales rangées comme la Fig. 2
                         (lX en abscisse, h en ordonnée), première graine ;
  - analyse_lambda.png : longueur d'onde mesurée et théorique en fonction
                         de h.

Utilisation :
    python analyse.py
    python analyse.py --out-dir out --dossier-resultats results
"""

import argparse
import csv
import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from simulation import VIDE, JAUNE, NOIR, COULEURS, BX, S
from dispersion import longueur_onde_turing, longueur_onde_instable

# 1. En dessous de cette fraction, une couleur est considérée comme absente
FRACTION_MIN = 0.02

# 2. Seuils de classer_motif. Caractéristique d'Euler du noir par maille du
#    motif : au-dessus de SEUIL_EULER, taches noires ; en dessous de
#    l'opposé, taches jaunes (un réseau hexagonal parfait donne ±0,87).
#    Anisotropie du spectre : au-dessus de SEUIL_ANISOTROPIE, rayures.
SEUIL_EULER = 0.4
SEUIL_ANISOTROPIE = 0.4

# 3. Colonnes de analyse.csv
COLONNES = ["L", "h", "lX", "bande", "graine", "t", "vide", "jaune", "noir", "motif",
            "lambda_mc", "lambda_sur_h", "lambda_taylor_20", "lambda_exact_1d", "lambda_exact_2d"]


def lire_arguments():
    """
    Lit les options --out-dir (défaut out) et --dossier-resultats
    (défaut results).

    Renvoie
    -------
    args : argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="Analyse des motifs Monte Carlo")
    parser.add_argument("--out-dir", default="out", help="dossier contenant les dossiers de run")
    parser.add_argument("--dossier-resultats", default="results", help="dossier des sorties")
    return parser.parse_args()


def charger_runs(out_dir):
    """
    Charge la grille finale et les paramètres de chaque run de out_dir.

    Un dossier sans params.json ou sans grille.npy est ignoré.

    Paramètres
    ----------
    out_dir : str
        Dossier contenant les dossiers de run.

    Renvoie
    -------
    runs : liste de dict
        Pour chaque run, les clés de params.json (L, h, lX, bande, graine,
        pas_faits, t, ...) et "grille" (np.ndarray, forme (L, L)), triés
        par h, lX, bande puis graine.
    """
    runs = []
    for nom in sorted(os.listdir(out_dir)):
        chemin_params = os.path.join(out_dir, nom, "params.json")
        chemin_grille = os.path.join(out_dir, nom, "grille.npy")
        if not (os.path.exists(chemin_params) and os.path.exists(chemin_grille)):
            continue
        with open(chemin_params) as f:
            run = json.load(f)
        run.pop("etat_generateur", None)
        run["grille"] = np.load(chemin_grille)
        runs.append(run)

    runs.sort(key=lambda r: (r["h"], r["lX"], r["bande"], r["graine"]))
    return runs


def fractions(grille):
    """
    Calcule la proportion de cases de chaque type.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau.

    Renvoie
    -------
    (vide, jaune, noir) : tuple de float
        Proportions, de somme 1.
    """
    return (float((grille == VIDE).mean()),
            float((grille == JAUNE).mean()),
            float((grille == NOIR).mean()))


def champ_signe(grille):
    """
    Convertit la grille en champ +1 (jaune), -1 (noir), 0 (vide).

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau.

    Renvoie
    -------
    champ : np.ndarray de float, forme (L, L)
    """
    return (grille == JAUNE).astype(float) - (grille == NOIR).astype(float)


def spectre_puissance(grille):
    """
    Spectre de puissance 2D du motif.

    Le champ +1 / -1 / 0 de champ_signe, moyenne retirée, est transformé
    par FFT 2D ; la puissance est |FFT|².

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau.

    Renvoie
    -------
    kx, ky : np.ndarray, forme (L, L)
        Composantes du nombre d'onde de chaque mode.
    puissance : np.ndarray, forme (L, L)
        Puissance de chaque mode.
    """
    L = grille.shape[0]
    champ = champ_signe(grille)
    puissance = np.abs(np.fft.fft2(champ - champ.mean())) ** 2
    k_1d = np.fft.fftfreq(L) * 2 * math.pi
    kx, ky = np.meshgrid(k_1d, k_1d, indexing="ij")
    return kx, ky, puissance


def anneau_dominant(grille):
    """
    Repère les modes de Fourier qui portent le motif.

    Les modes sont rangés en anneaux de même |k|, de largeur 2π/L (l'écart
    entre deux nombres d'onde de la grille). On cherche l'anneau de plus
    grande puissance moyenne, puis on garde ses modes et ceux de ses deux
    anneaux voisins.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau.

    Renvoie
    -------
    kx, ky, puissance : np.ndarray à une dimension
        Nombre d'onde et puissance des modes gardés.
    """
    L = grille.shape[0]
    kx, ky, puissance = spectre_puissance(grille)

    # Numéro d'anneau de chaque mode, puis puissance moyenne par anneau
    anneau = np.rint(np.hypot(kx, ky) / (2 * math.pi / L)).astype(int)
    n = L // 2 + 1
    compte = np.maximum(np.bincount(anneau.ravel())[:n], 1)
    radial = np.bincount(anneau.ravel(), puissance.ravel())[:n] / compte
    pic = int(radial[1:].argmax()) + 1

    garde = (np.abs(anneau - pic) <= 1) & (anneau > 0)
    return kx[garde], ky[garde], puissance[garde]


def longueur_onde(grille):
    """
    Longueur d'onde dominante du motif, en cases : λ = 2π / k_pic.

    k_pic est la moyenne de |k|, pondérée par la puissance, sur les modes
    de anneau_dominant : la mesure est ainsi plus fine que l'écart 2π/L
    entre deux anneaux.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau.

    Renvoie
    -------
    lam : float
        Longueur d'onde, ou nan si la grille est uniforme (pas de motif).
    """
    vide, jaune, noir = fractions(grille)
    if jaune < FRACTION_MIN or noir < FRACTION_MIN:
        return float("nan")

    kx, ky, puissance = anneau_dominant(grille)
    k_pic = (np.hypot(kx, ky) * puissance).sum() / puissance.sum()
    return 2 * math.pi / k_pic


def anisotropie(grille):
    """
    Mesure si le motif a une direction privilégiée, entre 0 et 1.

    Sur les modes de anneau_dominant, d'angle θ dans le plan (kx, ky) :
    |somme(P·exp(2iθ))| / somme(P). Des rayures parallèles mettent toute la
    puissance dans une direction : proche de 1. Des taches, ou un
    labyrinthe sans direction, la répartissent sur l'anneau : proche de 0.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau.

    Renvoie
    -------
    anisotropie : float
    """
    kx, ky, puissance = anneau_dominant(grille)
    theta = np.arctan2(ky, kx)
    return float(np.abs((puissance * np.exp(2j * theta)).sum()) / puissance.sum())


def lisser(grille):
    """
    Enlève le bruit d'une case : chaque case prend la couleur majoritaire
    de son voisinage 3 x 3 (bords périodiques).

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau.

    Renvoie
    -------
    noir : np.ndarray de booléens, forme (L, L)
        True où le noir est majoritaire dans le voisinage.
    """
    champ = champ_signe(grille)
    somme = np.zeros_like(champ)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            somme += np.roll(np.roll(champ, di, axis=0), dj, axis=1)
    return somme < 0


def caracteristique_euler(noir):
    """
    Caractéristique d'Euler d'un ensemble de cases : nombre de morceaux
    moins nombre de trous.

    Calculée sans chercher les morceaux un par un, avec la formule
    sommets - arêtes + faces sur les carrés des cases (les cases qui se
    touchent par un coin sont comptées dans le même morceau), bords
    périodiques.

    Paramètres
    ----------
    noir : np.ndarray de booléens, forme (L, L)
        Cases de l'ensemble.

    Renvoie
    -------
    chi : int
        Positif pour des taches isolées, négatif pour un fond troué,
        proche de 0 pour des rayures.
    """
    haut = np.roll(noir, 1, axis=0)
    gauche = np.roll(noir, 1, axis=1)
    coin = np.roll(haut, 1, axis=1)

    faces = int(noir.sum())
    aretes = int((noir | haut).sum()) + int((noir | gauche).sum())
    sommets = int((noir | haut | gauche | coin).sum())
    return sommets - aretes + faces


def classer_motif(grille):
    """
    Donne le type de motif de la grille.

    D'abord les fractions : s'il y a moins de FRACTION_MIN de noir (ou de
    jaune), la peau est « tout jaune » (ou « tout noir »). Ensuite, si le
    spectre a une direction privilégiée (anisotropie supérieure à
    SEUIL_ANISOTROPIE), ce sont des « rayures », même interrompues. Sinon
    on compte, sur la grille lissée, la caractéristique d'Euler du noir par
    maille du motif, chi·λ²/L² :
      - supérieure à SEUIL_EULER : le noir forme des morceaux isolés,
        « taches noires » ;
      - inférieure à -SEUIL_EULER : le noir est un fond troué de jaune,
        « taches jaunes » ;
      - entre les deux : « rayures » en labyrinthe.

    C'est un classement automatique simple : un motif intermédiaire (des
    rayures très morcelées, par exemple) peut tomber d'un côté ou de
    l'autre.

    Paramètres
    ----------
    grille : np.ndarray, forme (L, L)
        État de la peau.

    Renvoie
    -------
    motif : str
        "tout jaune", "taches noires", "rayures", "taches jaunes" ou
        "tout noir" ; "vide" pour une grille sans cellules (run à peine
        commencé).
    """
    vide, jaune, noir = fractions(grille)
    if jaune < FRACTION_MIN and noir < FRACTION_MIN:
        return "vide"
    if noir < FRACTION_MIN:
        return "tout jaune"
    if jaune < FRACTION_MIN:
        return "tout noir"

    if anisotropie(grille) > SEUIL_ANISOTROPIE:
        return "rayures"

    L = grille.shape[0]
    chi_par_maille = caracteristique_euler(lisser(grille)) * longueur_onde(grille) ** 2 / L**2
    if chi_par_maille > SEUIL_EULER:
        return "taches noires"
    if chi_par_maille < -SEUIL_EULER:
        return "taches jaunes"
    return "rayures"


def analyser_runs(runs):
    """
    Mesure chaque run et calcule les longueurs d'onde théoriques.

    Paramètres
    ----------
    runs : liste de dict
        Sortie de charger_runs.

    Renvoie
    -------
    lignes : liste de dict
        Une entrée par run, avec les clés de COLONNES. Les longueurs d'onde
        théoriques valent nan si lX <= bX (pas d'état mélangé) ou, pour le
        noyau exact, si aucun mode n'est amplifié à ce h.
    """
    lignes = []
    for run in runs:
        grille = run["grille"]
        vide, jaune, noir = fractions(grille)
        lam = longueur_onde(grille)

        # Théorie : éq. (20) au seuil, et mode le plus instable du noyau exact
        lam_taylor = lam_exact_1d = lam_exact_2d = float("nan")
        if run["lX"] > BX:
            lam_taylor = longueur_onde_turing(BX, run["lX"], S)
            lam_exact_1d = longueur_onde_instable(BX, run["lX"], S, run["h"], noyau="exact", dim=1)
            lam_exact_2d = longueur_onde_instable(BX, run["lX"], S, run["h"], noyau="exact", dim=2)

        lignes.append({
            "L": run["L"], "h": run["h"], "lX": run["lX"], "bande": run["bande"],
            "graine": run["graine"], "t": run["t"],
            "vide": round(vide, 3), "jaune": round(jaune, 3), "noir": round(noir, 3),
            "motif": classer_motif(grille),
            "lambda_mc": round(lam, 2), "lambda_sur_h": round(lam / run["h"], 2),
            "lambda_taylor_20": round(lam_taylor, 2),
            "lambda_exact_1d": round(lam_exact_1d, 2),
            "lambda_exact_2d": round(lam_exact_2d, 2),
        })
    return lignes


def ecrire_csv(lignes, chemin):
    """
    Écrit le tableau d'analyse dans un CSV et l'affiche.

    Paramètres
    ----------
    lignes : liste de dict
        Sortie de analyser_runs.
    chemin : str
        Fichier CSV à écrire.
    """
    with open(chemin, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLONNES)
        writer.writeheader()
        writer.writerows(lignes)

    for ligne in lignes:
        print(f"h = {ligne['h']:3d}  lX = {ligne['lX']:4}  graine = {ligne['graine']}  "
              f"{ligne['motif']:14s} λ = {ligne['lambda_mc']:6}  (exact 2D : {ligne['lambda_exact_2d']})")
    print(f"Tableau sauvegardé dans : {chemin}")


def tracer_motifs(runs, lignes, chemin):
    """
    Planche des grilles finales, rangées comme la Fig. 2 : lX en abscisse,
    h en ordonnée (croissant vers le haut), avec le type de motif et la
    longueur d'onde sous chaque grille.

    Seuls les runs sans iridophores sont affichés, et un seul par couple
    (h, lX) : celui de plus petite graine.

    Paramètres
    ----------
    runs : liste de dict
        Sortie de charger_runs.
    lignes : liste de dict
        Sortie de analyser_runs (même ordre que runs).
    chemin : str
        Fichier PNG à écrire.
    """
    # 1. Un run par couple (h, lX)
    cases = {}
    for run, ligne in zip(runs, lignes):
        if run["bande"] == 0 and (run["h"], run["lX"]) not in cases:
            cases[(run["h"], run["lX"])] = (run["grille"], ligne)
    if not cases:
        return
    h_valeurs = sorted({h for h, lX in cases}, reverse=True)
    lX_valeurs = sorted({lX for h, lX in cases})

    # 2. Une grille par case du tableau
    fig, axes = plt.subplots(len(h_valeurs), len(lX_valeurs), squeeze=False,
                             figsize=(2.6 * len(lX_valeurs), 2.9 * len(h_valeurs)))
    for i, h in enumerate(h_valeurs):
        for j, lX in enumerate(lX_valeurs):
            ax = axes[i, j]
            ax.set_xticks([])
            ax.set_yticks([])
            if j == 0:
                ax.set_ylabel(f"h = {h}")
            if i == 0:
                ax.set_title(f"lX = {lX}")
            if (h, lX) not in cases:
                ax.axis("off")
                continue
            grille, ligne = cases[(h, lX)]
            ax.imshow(grille, cmap=ListedColormap(COULEURS), vmin=0, vmax=2, interpolation="nearest")
            legende = ligne["motif"]
            if not math.isnan(ligne["lambda_mc"]):
                legende = legende + f", λ = {ligne['lambda_mc']:.1f}"
            ax.set_xlabel(legende)

    plt.savefig(chemin, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Graphique sauvegardé avec succès dans : {chemin}")


def tracer_lambda(lignes, chemin):
    """
    Trace la longueur d'onde en fonction de h : points mesurés sur les
    grilles Monte Carlo et courbes du noyau exact en dimension 2, une
    couleur par valeur de lX, plus la droite λ = 2h.

    Paramètres
    ----------
    lignes : liste de dict
        Sortie de analyser_runs.
    chemin : str
        Fichier PNG à écrire.
    """
    avec_motif = [l for l in lignes if not math.isnan(l["lambda_mc"])]
    if not avec_motif:
        return
    h_max = max(l["h"] for l in avec_motif)
    h_courbe = np.arange(2, h_max + 2)

    plt.figure(figsize=(8, 5))
    couleurs = ["blue", "orange", "green", "red", "purple", "brown"]
    for lX, couleur in zip(sorted({l["lX"] for l in avec_motif}), couleurs):
        points = [l for l in avec_motif if l["lX"] == lX]
        plt.scatter([l["h"] for l in points], [l["lambda_mc"] for l in points],
                    color=couleur, label=f"Monte Carlo, lX = {lX}")
        theorie = [longueur_onde_instable(BX, lX, S, h, noyau="exact", dim=2) for h in h_courbe]
        plt.plot(h_courbe, theorie, color=couleur, linestyle="--", label=f"noyau exact 2D, lX = {lX}")
    plt.plot(h_courbe, 2 * h_courbe, color="black", linewidth=0.8, label="λ = 2h")

    plt.xlabel("Portée h (cases)")
    plt.ylabel("Longueur d'onde λ (cases)")
    plt.title("Longueur d'onde des motifs : Monte Carlo et théorie linéaire")
    plt.legend(fontsize=8)
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.savefig(chemin, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Graphique sauvegardé avec succès dans : {chemin}")


if __name__ == "__main__":
    # 1. Grilles finales de tous les runs
    args = lire_arguments()
    runs = charger_runs(args.out_dir)
    if not runs:
        raise SystemExit(f"Aucun run dans {args.out_dir} : lancer d'abord simulation.py")
    print(f"{len(runs)} run(s) trouvé(s) dans {args.out_dir}")
    os.makedirs(args.dossier_resultats, exist_ok=True)

    # 2. Mesures et comparaison avec la théorie
    lignes = analyser_runs(runs)
    ecrire_csv(lignes, os.path.join(args.dossier_resultats, "analyse.csv"))

    # 3. Figures
    tracer_motifs(runs, lignes, os.path.join(args.dossier_resultats, "analyse_motifs.png"))
    tracer_lambda(lignes, os.path.join(args.dossier_resultats, "analyse_lambda.png"))
