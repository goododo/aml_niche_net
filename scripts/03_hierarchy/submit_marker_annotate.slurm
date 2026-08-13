#!/bin/bash
# submit_marker_annotate.slurm ----
# Marker-based cell typing (05_marker_annotate.R) for the whole cohort, one array task per sample.
#
# WHY THE WHOLE COHORT and not just GSE253355: the marker call is kept ALONGSIDE the BMM projection
# for every sample, so the two can be compared per cell. On GSE253355 that comparison already paid
# for itself -- the two methods agree 97% on Stromal and 99% on B_Plasma but only 3% on LMPP_GMP,
# and the bins they disagree on are exactly the bins whose healthy-donor FPR was worst. A dataset
# with no healthy donor has no FPR to check, so cross-method agreement is the only projection-QC
# signal available for it.
#
# STEP 1 build the task list (dataset<TAB>sample, one per line):
#   Rscript -e 'source("scripts/config/config_paths.R"); source("scripts/config/utils.R");
#               data.table::fwrite(qc_rds_roster(on_extra="ignore")[,.(dataset,sample)],
#                                  "marker_anno_tasks.tsv", sep="\t", col.names=FALSE)'
#
# STEP 2 submit from the project root:
#   N=$(wc -l < marker_anno_tasks.tsv)
#   sbatch --array=1-${N}%12 scripts/03_hierarchy/submit_marker_annotate.slurm marker_anno_tasks.tsv
#
#SBATCH -J markanno
#SBATCH -p gr10634c
#SBATCH -t 04:00:00
#SBATCH --rsc p=1:t=4:c=4:m=32G
#SBATCH -o logs/%x.%A_%a.out
#SBATCH -e logs/%x.%A_%a.err
#SBATCH --no-requeue

set -eux -o pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
TASKS="${1:?pass tasks tsv: dataset<TAB>sample per line}"
export TMPDIR="/LARGE1/gr10634/gaozy/tmp/markanno_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}"; trap 'rm -rf "${TMPDIR}"' EXIT
module load miniforge3

# Same .here anchor walk as 45_infercnv_submit.slurm: robust to submitting from any subdir.
_find_here() {
  local d; d="$(cd "${1:-$PWD}" 2>/dev/null && pwd)" || return 1
  while [[ "$d" != "/" ]]; do
    [[ -e "$d/.here" ]] && { printf '%s\n' "$d"; return 0; }
    d="$(dirname "$d")"
  done; return 1; }
PROJECT_ROOT="$(_find_here "${SLURM_SUBMIT_DIR:-$PWD}")" || true
[[ -n "${PROJECT_ROOT}" ]] || { echo "ERROR: no .here anchor found"; exit 1; }
source "${PROJECT_ROOT}/scripts/config/config_paths.sh"
WORKER="${PROJECT_ROOT}/scripts/03_hierarchy/05_marker_annotate.R"
[[ -f "${WORKER}" ]] || { echo "ERROR: ${WORKER} not found"; exit 1; }

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${TASKS}")
DS=$(echo "${LINE}" | cut -f1); SID=$(echo "${LINE}" | cut -f2)
echo "TASK ${SLURM_ARRAY_TASK_ID}: ${DS} :: ${SID}"
[[ -n "${DS}" && -n "${SID}" ]] || { echo "ERROR: empty dataset/sample"; exit 1; }

srun --export=ALL -n 1 conda run -p "${ENV_PREFIX}" \
  Rscript "${WORKER}" --dataset "${DS}" --sample "${SID}" --force
