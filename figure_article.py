"""
Figure 2 de l'article (diagramme de bifurcation) refaite avec nos runs.

Dans le plan (lX, h), chaque run est dessiné à sa position (lX, h) par sa
grille finale, par-dessus les deux courbes du champ moyen de l'article :
  - en bleu, h_T(lX) (éq. 19) : au-dessus, motif de Turing (fond bleu clair),
    en dessous la peau reste toute jaune ;
  - en rouge, la demi-longueur d'onde critique λ/2 (éq. 20).

Aucune valeur à choisir : le script lit tous les out/*/params.json et
trouve tout seul les valeurs de h, lX, L, bande et graine présentes.
  - Une figure par couple (L, bande) : results/figures/fig2_L{L}_b{bande}.png.
  - Plusieurs graines pour le même (L, bande, h, lX) : on montre le run le
    plus avancé (terminé de préférence), à égalité la plus petite graine.
  - Un run pas encore terminé (par exemple prolongé de 10⁷ à 10⁹ pas et en
    cours) est montré avec sa dernière sauvegarde et un cadre orange.
  - Un dossier sans params.json ou sans grille.npy, ou un params.json
    illisible, est ignoré (avec un message) : on peut ajouter n'importe quoi
    dans out/ sans casser la figure.
  - La taille des vignettes suit l'écart le plus petit entre deux valeurs de
    lX et entre deux valeurs de h : le balayage de l'article (4 x 3 runs)
    comme celui de 01_prepare_chunks.sh (30 x 19 runs) donnent une figure
    lisible.

Dépendances : numpy et matplotlib. États, couleurs, taux (BX, S) et masque
des iridophores viennent de simulation.py.

Utilisation (depuis la racine du projet) :
    python figure_article.py
    python figure_article.py --out-dir out --sortie results/figures
"""

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D

from simulation import COULEURS, BX, S, masque_iridophores

# Clés de params.json nécessaires pour placer un run sur la figure
CLES = ["L", "h", "lX", "bande", "graine", "n_pas", "pas_faits", "t"]

# Écart entre deux vignettes quand une seule valeur de lX ou de h est
# présente (pas d'écart mesurable) : celui de la Fig. 2 de l'article
ECART_LX_DEFAUT = 1.0
ECART_H_DEFAUT = 5.0


def lire_arguments():
    """
    Lit les options --out-dir (défaut out/) et --sortie
    (défaut results/figures/).

    Renvoie
    -------
    args : argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="Fig. 2 de l'article avec tous les runs de out/")
    parser.add_argument("--out-dir", default="out", help="dossier contenant les dossiers de run")
    parser.add_argument("--sortie", default="results/figures", help="dossier des figures")
    return parser.parse_args()


def lister_runs(out_dir):
    """
    Lit le params.json de chaque dossier de run de out_dir (un niveau :
    out_dir/<run>/params.json), sans charger les grilles.

    Un dossier sans grille.npy, ou un params.json illisible ou incomplet
    (clés CLES), est signalé et ignoré.

    Paramètres
    ----------
    out_dir : Path
        Dossier contenant les dossiers de run.

    Renvoie
    -------
    runs : list de dict
        Une entrée par run : les clés CLES, "termine" (pas_faits >= n_pas)
        et "dossier" (Path du run).
    """
    runs = []
    for fichier in sorted(out_dir.glob("*/params.json")):
        dossier = fichier.parent
        if not (dossier / "grille.npy").exists():
            print(f"  ignoré (pas de grille.npy) : {dossier}")
            continue
        try:
            with open(fichier) as f:
                infos = json.load(f)
            run = {c: infos[c] for c in CLES}
        except (OSError, ValueError, KeyError) as e:
            print(f"  ignoré (params.json illisible : {e}) : {dossier}")
            continue
        run["termine"] = run["pas_faits"] >= run["n_pas"]
        run["dossier"] = dossier
        runs.append(run)
    return runs


def choisir_representants(runs):
    """
    Garde un seul run par point (L, bande, h, lX) de la figure.

    Ordre de préférence : run terminé, puis plus de pas faits, puis plus
    petite graine.

    Paramètres
    ----------
    runs : list de dict
        Sortie de lister_runs.

    Renvoie
    -------
    groupes : dict
        (L, bande) -> list de dict, un run par couple (h, lX).
    """
    meilleurs = {}
    for run in runs:
        cle = (run["L"], run["bande"], run["h"], run["lX"])
        score = (run["termine"], run["pas_faits"], -run["graine"])
        if cle not in meilleurs or score > meilleurs[cle][0]:
            meilleurs[cle] = (score, run)

    groupes = {}
    for (L, bande, h, lX), (score, run) in sorted(meilleurs.items()):
        groupes.setdefault((L, bande), []).append(run)
    return groupes


def h_critique(lX, bX=BX, s=S):
    """
    Portée critique h_T(lX) du champ moyen (éq. 19 de l'article), pour
    lX > bX (NaN sinon : pas d'état mélangé, donc pas de Turing).

    Paramètres
    ----------
    lX : np.ndarray
        Valeurs de lX.
    bX, s : float
        Taux de naissance du jaune et d'élimination (ceux de simulation.py).

    Renvoie
    -------
    hT : np.ndarray, même forme que lX
    """
    lX = np.asarray(lX, dtype=float)
    hT = np.full_like(lX, np.nan)
    ok = lX > bX
    l = lX[ok]
    num = 2 * s * (l - bX) + l * (s + l) + 2 * np.sqrt(l * (s + bX) * s * (l - bX))
    hT[ok] = np.sqrt(num / ((l - bX) * (s + l)))
    return hT


def longueur_onde_critique(lX, bX=BX, s=S):
    """
    Longueur d'onde critique λ(lX) du champ moyen (éq. 20 de l'article, a = 1),
    pour lX > bX (NaN sinon).

    Paramètres
    ----------
    lX : np.ndarray
        Valeurs de lX.
    bX, s : float
        Taux de naissance du jaune et d'élimination.

    Renvoie
    -------
    lam : np.ndarray, même forme que lX
    """
    lX = np.asarray(lX, dtype=float)
    lam = np.full_like(lX, np.nan)
    ok = lX > bX
    l = lX[ok]
    q = (np.sqrt(l * (s + bX) * s * (l - bX)) - s * (l - bX)) / (s * bX)
    lam[ok] = np.sqrt(2) * np.pi / np.sqrt(q)
    return lam


def ecart_min(valeurs, defaut):
    """
    Plus petit écart entre deux valeurs distinctes (arrondi à 1e-6 pour
    absorber les 0.30000000000000004 de seq), ou defaut s'il n'y a qu'une
    valeur.

    Paramètres
    ----------
    valeurs : list de float
    defaut : float

    Renvoie
    -------
    ecart : float
    """
    v = np.unique(np.round(valeurs, 6))
    if len(v) < 2:
        return defaut
    return float(np.diff(v).min())


def grille_en_image(grille, bande):
    """
    Image RGB d'une grille (couleurs COULEURS de simulation.py), avec la
    bande d'iridophores éclaircie en gris argenté comme dans grille.png.

    Paramètres
    ----------
    grille : np.ndarray uint8, forme (L, L)
    bande : int

    Renvoie
    -------
    image : np.ndarray float, forme (L, L, 3), valeurs dans [0, 1]
    """
    palette = np.array([to_rgb(c) for c in COULEURS])
    image = palette[grille]
    masque = masque_iridophores(grille.shape[0], bande)
    image[masque] = 0.5 * image[masque] + 0.5 * np.array(to_rgb("silver"))
    return image


def tracer_figure2(runs, L, bande, chemin):
    """
    Trace la Fig. 2 pour un groupe de runs de même L et même bande.

    Mise en page : axes en unités de données (lX en abscisse, h en
    ordonnée), rapport d'aspect dx/dy pour que les vignettes, de taille
    0.9 dx x 0.9 dy, soient carrées à l'écran ; dx et dy sont les plus
    petits écarts entre valeurs de lX et de h. La taille de la figure
    suit le nombre de vignettes.

    Paramètres
    ----------
    runs : list de dict
        Un run par couple (h, lX) (sortie de choisir_representants).
    L, bande : int
        Taille de grille et bande d'iridophores communes aux runs.
    chemin : Path
        Fichier PNG à écrire.
    """
    lXs = [r["lX"] for r in runs]
    hs = [r["h"] for r in runs]
    dx = ecart_min(lXs, ECART_LX_DEFAUT)
    dy = ecart_min(hs, ECART_H_DEFAUT)

    # 1. Cadre : de 0 (comme l'article) jusqu'à la dernière vignette
    x_min = min(0.0, min(lXs) - 0.6 * dx)
    x_max = max(lXs) + 0.6 * dx
    y_min = min(0.0, min(hs) - 0.6 * dy)
    y_max = max(hs) + 0.6 * dy

    # 2. Taille de la figure : ~0.6 pouce par vignette, entre 8 et 30 pouces
    nx, ny = (x_max - x_min) / dx, (y_max - y_min) / dy
    largeur, hauteur = 0.6 * nx, 0.6 * ny
    facteur = min(max(8 / max(largeur, hauteur), 1.0), 30 / max(largeur, hauteur))
    fig, ax = plt.subplots(figsize=(largeur * facteur + 1.5, hauteur * facteur + 1.5))

    # 3. Courbes du champ moyen et zone de Turing (h > h_T) en bleu clair
    lx_fin = np.linspace(max(x_min, BX) + 1e-4, x_max, 2000)
    hT = h_critique(lx_fin)
    demi_lambda = longueur_onde_critique(lx_fin) / 2
    ax.fill_between(lx_fin, np.minimum(hT, y_max), y_max, color="#c9cdf2", zorder=0)
    ax.plot(lx_fin, hT, color="tab:blue", lw=2, zorder=5)          # par-dessus les vignettes :
    ax.plot(lx_fin, demi_lambda, color="tab:red", lw=2, zorder=5)   # visibles même si le balayage est dense

    # 4. Une vignette par run, centrée sur (lX, h)
    w, hh = 0.9 * dx, 0.9 * dy
    for run in runs:
        grille = np.load(run["dossier"] / "grille.npy")
        x, y = run["lX"], run["h"]
        ax.imshow(grille_en_image(grille, bande), interpolation="nearest", zorder=3,
                  extent=(x - w / 2, x + w / 2, y - hh / 2, y + hh / 2))
        couleur_cadre = "black" if run["termine"] else "darkorange"
        ax.add_patch(Rectangle((x - w / 2, y - hh / 2), w, hh, fill=False, zorder=4,
                               lw=0.5 if run["termine"] else 1.5, edgecolor=couleur_cadre))

    # 5. Axes, légende, titre
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect(dx / dy)
    ax.set_xlabel(r"$l_X$")
    ax.set_ylabel(r"$h$")

    legende = [Line2D([], [], color="tab:blue", lw=2, label=r"$h_T$ (éq. 19)"),
               Line2D([], [], color="tab:red", lw=2, label=r"$\lambda/2$ (éq. 20)"),
               Patch(facecolor="#c9cdf2", label="Turing (champ moyen)")]
    n_en_cours = sum(1 for r in runs if not r["termine"])
    if n_en_cours:
        legende.append(Patch(facecolor="none", edgecolor="darkorange", lw=1.5,
                             label=f"run non terminé ({n_en_cours})"))
    ax.legend(handles=legende, loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)

    temps = [r["t"] for r in runs]
    t_txt = f"t = {min(temps):.0f}" if min(temps) == max(temps) else f"t de {min(temps):.0f} à {max(temps):.0f}"
    ax.set_title(f"L = {L}, bande = {bande}, {len(runs)} runs, {t_txt} "
                 f"($b_X$ = {BX:g}, s = {S:g})")

    chemin.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(chemin, dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    args = lire_arguments()

    # 1. Runs présents dans out/
    runs = lister_runs(Path(args.out_dir))
    if not runs:
        raise SystemExit(f"Aucun run lisible dans {args.out_dir}")
    print(f"{len(runs)} run(s) trouvé(s) dans {args.out_dir}")

    # 2. Une figure par (L, bande), un run par (h, lX)
    for (L, bande), groupe in choisir_representants(runs).items():
        chemin = Path(args.sortie) / f"fig2_L{L}_b{bande}.png"
        tracer_figure2(groupe, L, bande, chemin)
        print(f"Écrit : {chemin} ({len(groupe)} runs)")
