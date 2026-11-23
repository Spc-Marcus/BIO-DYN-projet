#!/usr/bin/env bash
# Usage : bash slurm_pipeline/05_tests_extensions.sh
#
#   Teste les extensions du modèle aux points de la Fig. 2 : seuls h (5, 10, 15) et
#   lX (0.5, 1.5, 2.5, 3.5) varient, le reste est celui de l'article (L = 100,
#   10⁹ pas, graine 0). Une série de 12 runs par variante (VARIANTES), tous avec
#   --historique, écrits dans params_chunks_tests/ et lancés par 02_run_worker.slurm
#   (Job Array). Quand toutes les tâches sont finies, 06_timelines.slurm fait le GIF
#   et la frise de chaque run. Valeurs choisies et raisons : doc/tests_extensions.md.
#
#   Relancer : refaire la même commande. Les runs terminés et les GIF à jour sont sautés.

# ── Paramètres configurables ───────────────────────────────────────
CHUNK_SIZE=10
MAX_CONCURRENT_CPU=40
L=100
N_PAS=1000000000
H_VALUES=(5 10 15)
LX_VALUES=(0.5 1.5 2.5 3.5)

# Taux de blessure : 4 à 7 blessures par run selon lX (voir doc/tests_extensions.md)
LB=2.5e-08

# Témoin des blessures : morts naturelles (4) et (5) au taux dX = dM = LB × |D_R|,
# |D_R| = 81, 317, 709 cases pour R = h = 5, 10, 15
declare -A MORTS=([5]=2.025e-06 [10]=7.925e-06 [15]=1.7725e-05)

# Variantes : options de simulation.py ajoutées à chaque run (TAUX_MORT est remplacé par MORTS[h])
VARIANTES=(
    "--reseau hex"                         # réseau hexagonal
    "--rho0 0.01"                          # départ avec 1 % de jaunes, témoin du contact
    "--contact --rho0 0.01"                # naissances au contact
    "--lB ${LB}"                           # blessures au hasard
    "--lB ${LB} --contact --rho0 0.01"     # blessures au hasard et naissances au contact
    "--tB 50000"                           # une seule blessure au centre, à t = 50 000
    "--dX TAUX_MORT --dM TAUX_MORT"        # morts dispersées au même taux que les blessures
)

# ── Chemins ────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CHUNKS_DIR="${PROJECT_ROOT}/params_chunks_tests"

# ── Une ligne "L h lX bande graine n_pas options" par run ─────────
RUNS=()
for OPTIONS in "${VARIANTES[@]}"; do
    for H in "${H_VALUES[@]}"; do
        for LX in "${LX_VALUES[@]}"; do
            RUNS+=("${L} ${H} ${LX} 0 0 ${N_PAS} ${OPTIONS//TAUX_MORT/${MORTS[$H]}}")
        done
    done
done

# ── Chunks ─────────────────────────────────────────────────────────
rm -rf "${CHUNKS_DIR}"
mkdir -p "${CHUNKS_DIR}" "${PROJECT_ROOT}/out" "${PROJECT_ROOT}/logs"
NUM_CHUNKS=$(( (${#RUNS[@]} + CHUNK_SIZE - 1) / CHUNK_SIZE ))
for (( i=0; i<NUM_CHUNKS; i++ )); do
    printf '%s\n' "${RUNS[@]:i*CHUNK_SIZE:CHUNK_SIZE}" > "${CHUNKS_DIR}/chunk_$(( i + 1 )).txt"
done
echo "[LOG: Préparation] ${#VARIANTES[@]} variantes × ${#H_VALUES[@]} h × ${#LX_VALUES[@]} lX = ${#RUNS[@]} run(s) → ${NUM_CHUNKS} chunk(s) dans ${CHUNKS_DIR}" >&2

# ── Job Array des runs (avec historique), puis GIF et frises une fois tout fini ──
JOB_ID=$(sbatch --parsable \
    --array="1-${NUM_CHUNKS}%${MAX_CONCURRENT_CPU}" \
    --job-name="biodyn_tests" \
    --chdir="${PROJECT_ROOT}" \
    --export=ALL,CHUNKS_DIR=params_chunks_tests \
    "${SCRIPT_DIR}/02_run_worker.slurm" 1)
JOB_ID="${JOB_ID%%;*}"
if [[ -z "${JOB_ID}" ]]; then
    echo "[LOG: Erreur] Soumission des runs de test impossible" >&2
    exit 1
fi
echo "[LOG: Soumis] Runs de test — ${NUM_CHUNKS} tâches → job ${JOB_ID}" >&2

GIF_ID=$(sbatch --parsable \
    --dependency="afterany:${JOB_ID}" \
    --chdir="${PROJECT_ROOT}" \
    "${SCRIPT_DIR}/06_timelines.slurm")
echo "[LOG: Soumis] GIF et frises après le job ${JOB_ID} → job ${GIF_ID%%;*}" >&2
