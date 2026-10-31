#!/usr/bin/env bash
# Usage : bash slurm_pipeline/03_submit_all.sh
#
#   Soumet un seul Job Array : une tâche par chunk de params_chunks/
#   (créés par 01_prepare_chunks.sh).
#
#   Relancer : refaire la même commande. Les runs terminés sont sautés, les
#   runs interrompus (limite de temps, nœud perdu) reprennent depuis leur
#   dernière sauvegarde.

# ── Limite de concurrence ──────────────────────────────────────────
MAX_CONCURRENT_CPU=40

if [[ $# -gt 0 ]]; then
    echo "[LOG: Erreur] Option inconnue : $1" >&2
    echo "Usage : bash 03_submit_all.sh (les plages se règlent dans 01_prepare_chunks.sh)" >&2
    exit 1
fi

# ── Chemins ────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CHUNKS_DIR="${PROJECT_ROOT}/params_chunks"
WORKER_SCRIPT="${SCRIPT_DIR}/02_run_worker.slurm"

# ── Vérification préalable ────────────────────────────────────────
if [[ ! -f "${WORKER_SCRIPT}" ]]; then
    echo "[LOG: Erreur] Script worker introuvable : ${WORKER_SCRIPT}" >&2
    exit 1
fi

# ── Détection des chunks ───────────────────────────────────────────
NUM_CHUNKS=$(ls "${CHUNKS_DIR}"/chunk_*.txt 2>/dev/null | wc -l)
if [[ $NUM_CHUNKS -eq 0 ]]; then
    echo "[LOG: Erreur] Aucun chunk trouvé dans ${CHUNKS_DIR}." >&2
    echo "             Exécutez d'abord : bash slurm_pipeline/01_prepare_chunks.sh" >&2
    exit 1
fi

mkdir -p "${PROJECT_ROOT}/out"
mkdir -p "${PROJECT_ROOT}/logs"

echo "[LOG: Config] ${NUM_CHUNKS} chunk(s) dans ${CHUNKS_DIR}" >&2

# ── Pool CPU — un seul Job Array (un chunk par tâche) ─────────────
JOB_OUT=$(sbatch \
    --array="1-${NUM_CHUNKS}%${MAX_CONCURRENT_CPU}" \
    --job-name="biodyn_cpu" \
    --chdir="${PROJECT_ROOT}" \
    "${WORKER_SCRIPT}")
echo "[LOG: Soumis] Pool CPU — ${NUM_CHUNKS} tâches → ${JOB_OUT}" >&2

echo "[LOG: Terminé] Job soumis. Utilisez 'squeue -u \$USER' pour suivre l'avancement." >&2
