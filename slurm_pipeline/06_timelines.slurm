#!/usr/bin/env bash
#SBATCH --job-name=biodyn_gif
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm_gif_%j.out

# Usage : sbatch --chdir=<PROJECT_ROOT> 06_timelines.slurm
#
# GIF et frise (timeline.py) de chaque run de out/ lancé avec --historique,
# sauf si son timeline.gif est déjà plus récent que son historique. Soumis
# par 05_tests_extensions.sh une fois toutes les tâches des runs finies.
# Mémoire : environ 5 Mo par grille de l'historique (jusqu'à ~900 grilles
# pour un run avec des blessures).

set +u
. /local/env/envconda.sh
conda activate biodyn

EXIT_CODE=0
for DOSSIER in out/*/; do
    DOSSIER="${DOSSIER%/}"
    [[ -d "${DOSSIER}/historique" ]] || continue
    [[ "${DOSSIER}/timeline.gif" -nt "${DOSSIER}/historique" ]] && continue
    echo "[LOG: Timeline] ${DOSSIER}" >&2
    python timeline.py "${DOSSIER}" < /dev/null || EXIT_CODE=1
done

echo "[LOG: Fin] GIF et frises, code_sortie=${EXIT_CODE}" >&2
exit ${EXIT_CODE}
