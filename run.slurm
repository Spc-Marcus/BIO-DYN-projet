#!/usr/bin/env bash
#SBATCH --job-name=biodyn_run
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm_%j.out

# Usage (depuis la racine du projet, logs/ doit exister : mkdir -p logs) :
#   Lancer   : sbatch run.slurm --h 16 --lX 2.5
#   Relancer : sbatch run.slurm --reprise out/L100_h16_lX2.5_g0
#
# Les arguments sont passés tels quels à simulation.py.

# ── Environnement Conda ───────────────────────────────────────────
# Création une seule fois : conda create -n biodyn python numpy matplotlib
. /local/env/envconda.sh
conda activate biodyn

mkdir -p out

echo "[LOG: Démarrage] Job=${SLURM_JOB_ID} args=$*" >&2

python simulation.py "$@"

EXIT_CODE=$?
echo "[LOG: Fin] Job=${SLURM_JOB_ID} code_sortie=${EXIT_CODE}" >&2
exit ${EXIT_CODE}
