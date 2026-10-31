#!/usr/bin/env bash
#SBATCH --job-name=biodyn_run
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm_%A_%a.out

# Une tâche du job array, soumise par lancer.sh (ne pas lancer directement) :
#   sbatch --array=1-N run.slurm <liste_L> <liste_h> <liste_lX> <liste_graine> <n_pas>
#
# Chaque liste est de la forme "v1,v2,...". La tâche SLURM_ARRAY_TASK_ID
# prend une combinaison des listes (toutes les combinaisons sont couvertes
# par les tâches 1 à N).

IFS=',' read -ra LISTE_L <<< "$1"
IFS=',' read -ra LISTE_H <<< "$2"
IFS=',' read -ra LISTE_LX <<< "$3"
IFS=',' read -ra LISTE_GRAINE <<< "$4"
N_PAS=$5

# ── Combinaison de la tâche ────────────────────────────────────────
# I = TASK_ID - 1 est décomposé en indices : la graine varie le plus vite,
# puis lX, puis h, puis L.
I=$(( SLURM_ARRAY_TASK_ID - 1 ))
GRAINE=${LISTE_GRAINE[$(( I % ${#LISTE_GRAINE[@]} ))]}; I=$(( I / ${#LISTE_GRAINE[@]} ))
LX=${LISTE_LX[$(( I % ${#LISTE_LX[@]} ))]};             I=$(( I / ${#LISTE_LX[@]} ))
H=${LISTE_H[$(( I % ${#LISTE_H[@]} ))]};                I=$(( I / ${#LISTE_H[@]} ))
L=${LISTE_L[$I]}

if [[ -z "$L" || -z "$H" || -z "$LX" || -z "$GRAINE" || -z "$N_PAS" ]]; then
    echo "[LOG: Erreur] Paramètres invalides (TASK_ID=${SLURM_ARRAY_TASK_ID:-?})" >&2
    exit 1
fi

# ── Environnement Conda ───────────────────────────────────────────
# Création une seule fois : conda create -n biodyn python numpy matplotlib
. /local/env/envconda.sh
conda activate biodyn

echo "[LOG: Démarrage] Job=${SLURM_ARRAY_JOB_ID} Tâche=${SLURM_ARRAY_TASK_ID} L=${L} h=${H} lX=${LX} graine=${GRAINE} n_pas=${N_PAS}" >&2

python simulation.py --L "$L" --h "$H" --lX "$LX" --graine "$GRAINE" --n_pas "$N_PAS"

EXIT_CODE=$?
echo "[LOG: Fin] Tâche=${SLURM_ARRAY_TASK_ID} code_sortie=${EXIT_CODE}" >&2
exit ${EXIT_CODE}
