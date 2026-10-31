#!/usr/bin/env bash
#SBATCH --job-name=biodyn_worker
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm_%A_%a.out

# Usage : sbatch --chdir=<PROJECT_ROOT> --array=1-NUM_CHUNKS%LIMIT 02_run_worker.slurm
#
# La tâche N traite params_chunks/chunk_N.txt : une ligne
# "L h lX bande graine n_pas" par run, lancés l'un après l'autre avec
# simulation.py. Un run déjà terminé est sauté, un run interrompu reprend
# depuis sa dernière sauvegarde.
#
# Les chemins ci-dessous sont relatifs au répertoire racine du projet
# (positionné via --chdir lors de la soumission par 03_submit_all.sh).

set +u

# ── Environnement Conda ───────────────────────────────────────────
# Création une seule fois : conda create -n biodyn python numpy matplotlib
. /local/env/envconda.sh
conda activate biodyn

mkdir -p out logs

# ── Chunk de la tâche ─────────────────────────────────────────────
CHUNK_ID=${SLURM_ARRAY_TASK_ID}
CHUNK_FILE="params_chunks/chunk_${CHUNK_ID}.txt"

if [[ -z "$CHUNK_ID" || ! -f "$CHUNK_FILE" ]]; then
    echo "[LOG: Erreur] Paramètres invalides (TASK_ID=${SLURM_ARRAY_TASK_ID:-?}, fichier=${CHUNK_FILE})" >&2
    exit 1
fi

echo "[LOG: Démarrage] Job=${SLURM_JOB_ID} Tâche=${SLURM_ARRAY_TASK_ID} chunk=${CHUNK_ID} cwd=$(pwd)" >&2

EXIT_CODE=0
while read -r L H LX BANDE GRAINE N_PAS; do
    echo "[LOG: Run] L=${L} h=${H} lX=${LX} bande=${BANDE} graine=${GRAINE} n_pas=${N_PAS}" >&2

    python simulation.py \
        --L "${L}" \
        --h "${H}" \
        --lX "${LX}" \
        --bande "${BANDE}" \
        --graine "${GRAINE}" \
        --n_pas "${N_PAS}" \
        < /dev/null || EXIT_CODE=1
done < "${CHUNK_FILE}"

echo "[LOG: Fin] Tâche=${SLURM_ARRAY_TASK_ID} chunk=${CHUNK_ID} code_sortie=${EXIT_CODE}" >&2
exit ${EXIT_CODE}
