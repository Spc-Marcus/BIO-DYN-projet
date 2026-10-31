"""
04_merge_results.py — Récupère les paramètres et l'état de tous les runs dans un seul CSV.

Usage (depuis la racine du projet) :
    python slurm_pipeline/04_merge_results.py [--out-dir out] [--output results/final_runs.csv]

Entrées : les params.json écrits par simulation.py, un par dossier de run
(out/L{L}_h{h}_lX{lX}_b{bande}_g{graine}/).

Sortie : une ligne par run, avec les colonnes
    L, h, lX, bande, graine, n_pas, pas_faits, t, termine, dossier
où termine vaut True si pas_faits >= n_pas. Les runs non terminés sont à
relancer avec 03_submit_all.sh.
"""

import argparse
import csv
import json
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("merge_results")

SCRIPT_DIR   = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUT_DIR      = PROJECT_ROOT / "out"
OUTPUT_PATH  = PROJECT_ROOT / "results" / "final_runs.csv"
COLUMNS = ["L", "h", "lX", "bande", "graine", "n_pas", "pas_faits", "t", "termine", "dossier"]


def lire_arguments():
    """
    Lit les options --out-dir (défaut out/) et --output
    (défaut results/final_runs.csv).

    Renvoie
    -------
    args : argparse.Namespace
    """
    raise NotImplementedError


def charger_runs(out_dir):
    """
    Lit le params.json de chaque dossier de run de out_dir.

    Un fichier illisible est signalé dans le log et ignoré.

    Paramètres
    ----------
    out_dir : Path
        Dossier contenant les dossiers de run.

    Renvoie
    -------
    runs : list de dict
        Une entrée par run, avec les clés de COLUMNS (termine calculé ici,
        dossier = nom du dossier du run), triée par h, lX, bande puis graine.
    """
    raise NotImplementedError


def ecrire_csv(runs, output):
    """
    Écrit les runs dans un CSV (colonnes COLUMNS, dossier parent créé si
    besoin) et logue le nombre de runs terminés et à relancer.

    Paramètres
    ----------
    runs : list de dict
        Sortie de charger_runs.
    output : Path
        Fichier CSV de sortie.
    """
    raise NotImplementedError


def main():
    """Point d'entrée : chargement des params.json et écriture du CSV."""
    args = lire_arguments()
    runs = charger_runs(Path(args.out_dir))
    ecrire_csv(runs, Path(args.output))


if __name__ == "__main__":
    main()
