#!/bin/bash
# 45_infercnv_submit.slurm ----
# Parallel inferCNV via SLURM array; one task per sample, calling 44_infercnv_run_one.R.
#
# STEP 0 (ONCE) build the external-reference cache so parallel tasks don't race:
#   module load miniforge3
#   conda run -p /FAST/gr10634/gaozy/general_env \
#     Rscript scripts/02_malignancy/44_infercnv_run_one.R --prebuild_ext_ref
#
# STEP 1 build the task list (dataset<TAB>sample) for the datasets you want:
#   awk -F, '$2=="GSE227903"||$2=="GSE239721"{print $2"\t"$1}' \
#       /FAST/gr10634/gaozy/aml_niche_net/results/tables/02_malignancy/ref_norm_summary.csv \
#       > infercnv_tasks.tsv
#   (the header row won't match GSE..., so it's skipped automatically)
#
# STEP 2 submit (from the project root or any subdir inside the tree):
#   mkdir -p logs
#   N=$(wc -l < infercnv_tasks.tsv)
#   sbatch --array=1-${N}%6 scripts/02_malignancy/45_infercnv_submit.slurm infercnv_tasks.tsv
#
# 3a-2 rewire (no logic change): the worker is located via the .here anchor (walked up from
#   ${SLURM_SUBMIT_DIR}, robust to submitting from any subdir + SLURM spooling) and is now
#   44_infercnv_run_one.R under scripts/02_malignancy (was: c46 assumed a sibling of the submit
#   dir); the conda env is ${ENV_PREFIX} from config_paths.sh (was hard-coded). SBATCH directives,
#   resources, and task logic byte-preserved.
#
#SBATCH -J infercnv
#SBATCH -p gr10634c
#SBATCH -t 24:00:00
#SBATCH --rsc p=1:t=8:c=8:m=64G
#SBATCH -o logs/%x.%A_%a.out
#SBATCH -e logs/%x.%A_%a.err
#SBATCH --no-requeue

set -eux -o pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export INFERCNV_THREADS=8
TASKS="${1:?pass tasks tsv: dataset<TAB>sample per line}"
# TODO(config-completion): this /LARGE1/gr10634/gaozy/tmp base is hard-coded (not under
#   LARGE1_DIR, so not directly derivable); lift to a SCRATCH_TMP_BASE config var later.
export TMPDIR="/LARGE1/gr10634/gaozy/tmp/icnv_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}"; trap 'rm -rf "${TMPDIR}"' EXIT
module load miniforge3

# Locate the project root via the .here anchor (walked up from the submit dir), then source
# config for ENV_PREFIX. Robust to submitting from the project root or any subdir, and to spooling.
_find_here() {
  local d; d="$(cd "${1:-$PWD}" 2>/dev/null && pwd)" || return 1
  while [[ "$d" != "/" ]]; do
    [[ -e "$d/.here" ]] && { printf '%s\n' "$d"; return 0; }
    d="$(dirname "$d")"
  done; return 1; }
_start="${SLURM_SUBMIT_DIR:-$PWD}"
PROJECT_ROOT="$(_find_here "${_start}")" || true
[[ -n "${PROJECT_ROOT}" ]] || { echo "ERROR: no .here anchor found walking up from ${_start} (sbatch from inside the project tree)"; exit 1; }
source "${PROJECT_ROOT}/scripts/config/config_paths.sh"
MAL_DIR="${PROJECT_ROOT}/scripts/02_malignancy"
[[ -f "${MAL_DIR}/44_infercnv_run_one.R" ]] || {
  echo "ERROR: 44_infercnv_run_one.R not found in ${MAL_DIR}"; exit 1; }

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${TASKS}")
DS=$(echo "${LINE}" | cut -f1); SID=$(echo "${LINE}" | cut -f2)
echo "TASK ${SLURM_ARRAY_TASK_ID}: ${DS} :: ${SID}"
[[ -n "${DS}" && -n "${SID}" ]] || { echo "ERROR: empty dataset/sample"; exit 1; }

srun --export=ALL -n 1 conda run -p "${ENV_PREFIX}" \
  Rscript "${MAL_DIR}/44_infercnv_run_one.R" --dataset "${DS}" --sample "${SID}"
