"""
Timeline d'un run : GIF et planche d'images à partir de son historique.

Indépendant de Slurm : se lance en local (après avoir rapatrié out/) ou sur
la frontale, sur le dossier d'un run lancé avec --historique.

Entrée : out/<run>/historique/grille_<pas>.npz, une grille tous les
PAS_ENTRE_SAUVEGARDES pas, stockée en deux matrices binaires compactées et
compressées "jaune" et "noir" (voir sauvegarder_run dans simulation.py).

Sorties, écrites dans le dossier du run :
  - timeline.gif : animation de la grille au cours du temps ;
  - timeline.png : planche de quelques instants régulièrement espacés, avec
                   le temps t = pas / L² sous chaque image.

Utilisation :
    python timeline.py out/L100_h16_lX2.5_b0_g0
    python timeline.py out/L100_h16_lX2.5_b0_g0 --fps 5 --n_images 6
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import to_rgb

from simulation import VIDE, JAUNE, NOIR, COULEURS


def lire_arguments():
    """
    Lit la ligne de commande.

    Arguments :
      dossier     dossier du run (positionnel), par exemple
                  out/L100_h16_lX2.5_b0_g0 ;
      --fps       images par seconde du GIF (défaut 10) ;
      --n_images  nombre d'instants sur la planche (défaut 8).

    Renvoie
    -------
    args : argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="GIF et planche d'un run lancé avec --historique")
    parser.add_argument("dossier", help="dossier du run, par exemple out/L100_h16_lX2.5_b0_g0")
    parser.add_argument("--fps", type=int, default=10, help="images par seconde du GIF")
    parser.add_argument("--n_images", type=int, default=8, help="nombre d'instants sur la planche")
    return parser.parse_args()


def lister_instantanes(dossier):
    """
    Liste les grilles sauvegardées dans dossier/historique/, dans l'ordre du
    temps.

    Paramètres
    ----------
    dossier : Path
        Dossier du run.

    Renvoie
    -------
    instantanes : list de (int, Path)
        Couples (pas, fichier), triés par pas croissant. Le pas est lu dans
        le nom du fichier grille_<pas>.npz.

    Lève
    ----
    FileNotFoundError
        Si dossier/historique/ est absent ou vide (run lancé sans
        --historique).
    """
    fichiers = sorted((dossier / "historique").glob("grille_*.npz"))
    if not fichiers:
        raise FileNotFoundError(f"Aucune grille dans {dossier / 'historique'} : run lancé sans --historique ?")

    # Le pas est le nombre écrit dans le nom : grille_000010000000.npz -> 10000000
    instantanes = [(int(f.stem.split("_")[1]), f) for f in fichiers]
    instantanes.sort()
    return instantanes


def charger_instantanes(instantanes):
    """
    Charge les grilles listées par lister_instantanes.

    Chaque grille est reconstruite à partir de ses deux matrices binaires :
    L = jaune.shape[0], puis
    grille = JAUNE * np.unpackbits(jaune, axis=-1, count=L)
             + NOIR * np.unpackbits(noir, axis=-1, count=L),
    les cases ni jaunes ni noires valant VIDE.

    Paramètres
    ----------
    instantanes : list de (int, Path)
        Sortie de lister_instantanes.

    Renvoie
    -------
    pas : np.ndarray d'entiers, forme (n,)
        Pas de chaque grille.
    grilles : np.ndarray uint8, forme (n, L, L)
        Les grilles, dans le même ordre.
    """
    pas = []
    grilles = []
    for p, fichier in instantanes:
        donnees = np.load(fichier)
        jaune = donnees["jaune"]
        noir = donnees["noir"]
        L = jaune.shape[0]

        # Reconstruction : JAUNE où jaune vaut 1, NOIR où noir vaut 1, VIDE ailleurs
        grille = JAUNE * np.unpackbits(jaune, axis=-1, count=L) + NOIR * np.unpackbits(noir, axis=-1, count=L)
        pas.append(p)
        grilles.append(grille.astype(np.uint8))

    return np.array(pas), np.array(grilles)


def en_couleurs(grilles):
    """
    Convertit les grilles en images : VIDE en blanc, JAUNE en jaune, NOIR
    en noir (mêmes couleurs que grille.png de simulation.py).

    Paramètres
    ----------
    grilles : np.ndarray, forme (n, L, L)
        Grilles de VIDE, JAUNE et NOIR.

    Renvoie
    -------
    images : np.ndarray uint8, forme (n, L, L, 3)
        Images RGB.
    """
    # Une couleur RGB par état, dans l'ordre VIDE, JAUNE, NOIR
    palette = (np.array([to_rgb(c) for c in COULEURS]) * 255).astype(np.uint8)
    return palette[grilles]


def creer_gif(images, temps, chemin, fps):
    """
    Écrit l'animation des images, avec le temps t en titre de chaque image.

    Paramètres
    ----------
    images : np.ndarray uint8, forme (n, L, L, 3)
        Images RGB, dans l'ordre du temps.
    temps : np.ndarray, forme (n,)
        Temps t = pas / L² de chaque image.
    chemin : Path
        Fichier GIF à écrire.
    fps : int
        Images par seconde.
    """
    fig = plt.figure(figsize=(5, 5))
    image = plt.imshow(images[0], interpolation="nearest")
    titre = plt.title(f"t = {temps[0]:.0f}")
    plt.axis("off")

    # Image k de l'animation : on change la grille et le titre
    def mise_a_jour(k):
        image.set_data(images[k])
        titre.set_text(f"t = {temps[k]:.0f}")
        return image, titre

    anim = animation.FuncAnimation(fig, mise_a_jour, frames=len(images))
    anim.save(chemin, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)


def creer_planche(images, temps, chemin, n_images):
    """
    Écrit une planche de n_images instants régulièrement espacés, du premier
    au dernier, côte à côte, avec le temps t sous chaque image.

    Paramètres
    ----------
    images : np.ndarray uint8, forme (n, L, L, 3)
        Images RGB, dans l'ordre du temps.
    temps : np.ndarray, forme (n,)
        Temps t = pas / L² de chaque image.
    chemin : Path
        Fichier PNG à écrire.
    n_images : int
        Nombre d'instants affichés (au plus n).
    """
    # Indices de n_images instants régulièrement espacés, premier et dernier inclus
    n = min(n_images, len(images))
    indices = np.unique(np.linspace(0, len(images) - 1, n).round().astype(int))

    fig, axes = plt.subplots(1, len(indices), figsize=(2.5 * len(indices), 3), squeeze=False)
    for ax, k in zip(axes[0], indices):
        ax.imshow(images[k], interpolation="nearest")
        ax.set_xlabel(f"t = {temps[k]:.0f}")
        ax.set_xticks([])
        ax.set_yticks([])

    plt.savefig(chemin, dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    # 1. Grilles sauvegardées pendant le run
    args = lire_arguments()
    dossier = Path(args.dossier)
    pas, grilles = charger_instantanes(lister_instantanes(dossier))
    temps = pas / grilles.shape[1] ** 2
    print(f"{len(pas)} instantanés, de t = {temps[0]} à t = {temps[-1]}")

    # 2. Images en couleurs
    images = en_couleurs(grilles)

    # 3. GIF et planche, dans le dossier du run
    creer_gif(images, temps, dossier / "timeline.gif", args.fps)
    creer_planche(images, temps, dossier / "timeline.png", args.n_images)
    print(f"Écrits : {dossier / 'timeline.gif'} et {dossier / 'timeline.png'}")
