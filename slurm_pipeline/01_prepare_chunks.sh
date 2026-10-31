#!/usr/bin/env bash
# Usage : bash slurm_pipeline/01_prepare_chunks.sh [--L N] \
#             [--h-start A] [--h-end B] [--h-step S] \
#             [--lx-start A] [--lx-end B] [--lx-step S] \
#             [--graines N] [--n-pas N]
#
#   Génère une ligne "L h lX graine n_pas" par run (toutes les combinaisons
#   de h, lX et graine) et les répartit dans params_chunks/chunk_1.txt,
#   chunk_2.txt, … (CHUNK_SIZE lignes max). Une tâche Slurm = un chunk.
#
#   --L         Taille de la grille. Défaut : 200
#   --h-start   Première valeur de h. Défaut : 2
#   --h-end     Dernière valeur de h, incluse. Défaut : 20
#   --h-step    Pas du balayage de h. Défaut : 1
#   --lx-start  Première valeur de lX. Défaut : 1.1
#   --lx-end    Dernière valeur de lX, incluse. Défaut : 4.0
#   --lx-step   Pas du balayage de lX. Défaut : 0.1
#   --graines   Nombre de répétitions (graines 0 à N-1). Défaut : 3
#   --n-pas     Nombre de pas par run. Défaut : 100000 × L²

# ── Paramètre configurable ─────────────────────────────────────────
CHUNK_SIZE=10

# ── Valeurs par défaut (balayage de pseudocode.md) ─────────────────
DEFAULT_L=200
DEFAULT_H_START=2
DEFAULT_H_END=20
DEFAULT_H_STEP=1
DEFAULT_LX_START=1.1
DEFAULT_LX_END=4.0
DEFAULT_LX_STEP=0.1
DEFAULT_GRAINES=3

# ── Parsing des arguments ──────────────────────────────────────────
L="${DEFAULT_L}"
H_START="${DEFAULT_H_START}"
H_END="${DEFAULT_H_END}"
H_STEP="${DEFAULT_H_STEP}"
LX_START="${DEFAULT_LX_START}"
LX_END="${DEFAULT_LX_END}"
LX_STEP="${DEFAULT_LX_STEP}"
GRAINES="${DEFAULT_GRAINES}"
N_PAS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --L)        L="$2";        shift 2 ;;
        --h-start)  H_START="$2";  shift 2 ;;
        --h-end)    H_END="$2";    shift 2 ;;
        --h-step)   H_STEP="$2";   shift 2 ;;
        --lx-start) LX_START="$2"; shift 2 ;;
        --lx-end)   LX_END="$2";   shift 2 ;;
        --lx-step)  LX_STEP="$2";  shift 2 ;;
        --graines)  GRAINES="$2";  shift 2 ;;
        --n-pas)    N_PAS="$2";    shift 2 ;;
        *)
            echo "[LOG: Erreur] Option inconnue : $1" >&2
            echo "Usage : bash 01_prepare_chunks.sh [--L N] [--h-start A] [--h-end B] [--h-step S] [--lx-start A] [--lx-end B] [--lx-step S] [--graines N] [--n-pas N]" >&2
            exit 1
            ;;
    esac
done

N_PAS="${N_PAS:-$(( 100000 * L * L ))}"

# ── Chemins ────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CHUNKS_DIR="${PROJECT_ROOT}/params_chunks"

# ── Nettoyage et recréation ────────────────────────────────────────
echo "[LOG: Nettoyage] Suppression de ${CHUNKS_DIR}" >&2
rm -rf "${CHUNKS_DIR}"
mkdir -p "${CHUNKS_DIR}"

# ── Grilles de h et lX (LC_ALL=C : point décimal) ─────────────────
mapfile -t H_VALUES  < <(LC_ALL=C seq "${H_START}"  "${H_STEP}"  "${H_END}")
mapfile -t LX_VALUES < <(LC_ALL=C seq "${LX_START}" "${LX_STEP}" "${LX_END}")

# ── Une ligne par run ─────────────────────────────────────────────
RUNS=()
for H in "${H_VALUES[@]}"; do
    for LX in "${LX_VALUES[@]}"; do
        for (( G=0; G<GRAINES; G++ )); do
            RUNS+=("${L} ${H} ${LX} ${G} ${N_PAS}")
        done
    done
done

TOTAL=${#RUNS[@]}
if [[ $TOTAL -eq 0 ]]; then
    echo "[LOG: Erreur] Aucun run : vérifier les plages (start <= end ?) et --graines" >&2
    exit 1
fi

NUM_CHUNKS=$(( (TOTAL + CHUNK_SIZE - 1) / CHUNK_SIZE ))
echo "[LOG: Préparation] L=${L} h=[${H_START}..${H_END}/${H_STEP}] lX=[${LX_START}..${LX_END}/${LX_STEP}] graines=${GRAINES} n_pas=${N_PAS}" >&2
echo "[LOG: Préparation] ${TOTAL} run(s) → ${NUM_CHUNKS} chunk(s) de ${CHUNK_SIZE} max" >&2

# ── Écriture des chunks ───────────────────────────────────────────
for (( i=0; i<NUM_CHUNKS; i++ )); do
    CHUNK_ID=$(( i + 1 ))
    CHUNK_FILE="${CHUNKS_DIR}/chunk_${CHUNK_ID}.txt"

    START=$(( i * CHUNK_SIZE ))
    printf '%s\n' "${RUNS[@]:START:CHUNK_SIZE}" > "${CHUNK_FILE}"

    echo "[LOG: Chunk ${CHUNK_ID}] $(wc -l < "${CHUNK_FILE}") run(s) → ${CHUNK_FILE}" >&2
done

echo "[LOG: Terminé] ${NUM_CHUNKS} chunk(s) créé(s) dans ${CHUNKS_DIR}" >&2
echo "Utilisez --array=1-${NUM_CHUNKS} dans vos jobs Slurm." >&2
echo "${NUM_CHUNKS}"
