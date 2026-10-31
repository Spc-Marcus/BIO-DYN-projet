#!/usr/bin/env bash
# Lance (ou relance) les runs sur GenOuest : un job array Slurm avec une
# tâche par combinaison de paramètres.
#
# Usage : bash lancer.sh [--L V] [--h V] [--lX V] [--graine V] [--n_pas N]
#
#   Chaque V est une valeur seule (16) ou une plage MIN:MAX:PAS, bornes
#   incluses (2:20:2 donne 2 4 6 ... 20). Toutes les combinaisons des
#   plages sont lancées. L, h, graine et n_pas sont des entiers.
#
# Exemples :
#   bash lancer.sh --h 16 --lX 2.5                                # un seul run
#   bash lancer.sh --h 2:20:2 --lX 1.1:4.0:0.1 --graine 0:2:1     # 10 x 30 x 3 runs
#
# Relancer : refaire la même commande. Les runs terminés sont sautés, les
# runs interrompus reprennent depuis leur dernière sauvegarde.

# ── Valeurs par défaut (Fig. 1) ────────────────────────────────────
L="100"
H="16"
LX="2.5"
GRAINE="0"
N_PAS="1000000000"

# ── Nombre maximal de tâches simultanées ──────────────────────────
MAX_SIMULTANES=40

# ── Parsing des arguments ──────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --L)      L="$2";      shift 2 ;;
        --h)      H="$2";      shift 2 ;;
        --lX)     LX="$2";     shift 2 ;;
        --graine) GRAINE="$2"; shift 2 ;;
        --n_pas)  N_PAS="$2";  shift 2 ;;
        *)
            echo "[LOG: Erreur] Option inconnue : $1" >&2
            echo "Usage : bash lancer.sh [--L V] [--h V] [--lX V] [--graine V] [--n_pas N]" >&2
            exit 1
            ;;
    esac
done

# ── Valeur seule ou plage MIN:MAX:PAS → liste "v1,v2,..." ─────────
# LC_ALL=C : point décimal même si la session est en français.
plage() {
    local MIN MAX PAS
    IFS=':' read -r MIN MAX PAS <<< "$1"
    if [[ -z "$MAX" ]]; then
        echo "$MIN"
    else
        LC_ALL=C seq "$MIN" "${PAS:-1}" "$MAX" | paste -sd, -
    fi
}

taille() {
    awk -F, '{ print NF }' <<< "$1"
}

LISTE_L=$(plage "$L")
LISTE_H=$(plage "$H")
LISTE_LX=$(plage "$LX")
LISTE_GRAINE=$(plage "$GRAINE")

N=$(( $(taille "$LISTE_L") * $(taille "$LISTE_H") * $(taille "$LISTE_LX") * $(taille "$LISTE_GRAINE") ))

if [[ $N -eq 0 ]]; then
    echo "[LOG: Erreur] Aucune combinaison : vérifier les plages (MIN <= MAX ?)" >&2
    exit 1
fi

# ── Limite de taille des job arrays sur le cluster ────────────────
MAX_ARRAY=$(scontrol show config 2>/dev/null | awk '/^MaxArraySize/ { print $3 }')
if [[ -n "$MAX_ARRAY" && $N -ge $MAX_ARRAY ]]; then
    echo "[LOG: Erreur] ${N} combinaisons, limite du cluster : $(( MAX_ARRAY - 1 )). Réduire les plages." >&2
    exit 1
fi

echo "[LOG: Config] L=[${LISTE_L}] h=[${LISTE_H}] lX=[${LISTE_LX}] graine=[${LISTE_GRAINE}] n_pas=${N_PAS}" >&2

# ── Soumission ─────────────────────────────────────────────────────
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/out"

JOB_OUT=$(sbatch \
    --array="1-${N}%${MAX_SIMULTANES}" \
    --chdir="${PROJECT_ROOT}" \
    "${PROJECT_ROOT}/run.slurm" "${LISTE_L}" "${LISTE_H}" "${LISTE_LX}" "${LISTE_GRAINE}" "${N_PAS}")

echo "[LOG: Soumis] ${N} tâche(s) → ${JOB_OUT}" >&2
echo "Suivre l'avancement : squeue -u \$USER" >&2
