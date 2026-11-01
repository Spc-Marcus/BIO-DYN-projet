"""
Timeline d'un run : GIF et planche d'images à partir de son historique.

Indépendant de Slurm : se lance en local (après avoir rapatrié out/) ou sur
la frontale, sur le dossier d'un run lancé avec --historique.

Entrée : out/<run>/historique/grille_<pas>.npy, une grille tous les
PAS_ENTRE_SAUVEGARDES pas (voir simulation.py).

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

from simulation import VIDE, JAUNE, NOIR


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
    raise NotImplementedError


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
        le nom du fichier grille_<pas>.npy.

    Lève
    ----
    FileNotFoundError
        Si dossier/historique/ est absent ou vide (run lancé sans
        --historique).
    """
    raise NotImplementedError


def charger_instantanes(instantanes):
    """
    Charge les grilles listées par lister_instantanes.

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
    raise NotImplementedError


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
    raise NotImplementedError


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
    raise NotImplementedError


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
    raise NotImplementedError


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
