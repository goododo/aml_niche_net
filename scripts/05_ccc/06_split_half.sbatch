#!/bin/bash
# 06_split_half.sbatch ----
# SLURM job array for the split-half reliability run (08_scoring D5). One task per HALF, so
# 2 tasks per paired sample. Each task runs
#   Rscript scripts/05_ccc/02_run_cellchat.R --dataset=<ds> --sample=<smp> \
#           --cell_subset=<half csv> --out_suffix=split_<half>
# reading row $SLURM_ARRAY_TASK_ID of DIR_CCC/split_half/tasks.csv (from 06_split_half_reliability.R).
#
# THE PRODUCTION TENSORS ARE NEVER TOUCHED. --out_suffix diverts every write to
#   <CCC_TENSOR_DIR>__split_A / __split_B  and  <CCC_GRAPH_DIR>__split_A / __split_B.
# The runner's staleness guard also compares against files in the suffixed directory, so a
# resubmit fills gaps rather than redoing finished halves.
#
# Submit (count comes from the manifest; nothing hard-coded that can go stale):
#   cd <project root>
#   Rscript scripts/05_ccc/06_split_half_reliability.R          # writes tasks.csv first
#   mkdir -p logs
#   N=$(( $(wc -l < results/tables/05_ccc/split_half/tasks.csv) - 1 ))
#   echo "half-tasks = ${N}"      # expect 74
#   sbatch --array=1-${N}%8 scripts/05_ccc/06_split_half.sbatch
#
# Resources copied from 02_run_cellchat.sbatch. A half has at most ~4.7k cells, i.e. strictly
# smaller than the production run this budget was sized for, so the ceiling is unchanged.
#
#SBATCH -J ccc_splithalf
#SBATCH -p gr10634c
#SBATCH -t 02:00:00
#SBATCH --rsc p=1:t=2:c=2:m=36G
#SBATCH -o logs/%x.%A_%a.out
#SBATCH -e logs/%x.%A_%a.err
#SBATCH --requeue

set -eux -o pipefail
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
       NUMEXPR_NUM_THREADS=2 TBB_NUM_THREADS=2

export TMPDIR="/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}"; trap 'rm -rf "${TMPDIR}"' EXIT

_find_here() {
  local d; d="$(cd "${1:-$PWD}" 2>/dev/null && pwd)" || return 1
  while [[ "$d" != "/" ]]; do
    [[ -e "$d/.here" ]] && { printf '%s\n' "$d"; return 0; }
    d="$(dirname "$d")"
  done; return 1; }
_start="${SLURM_SUBMIT_DIR:-$PWD}"
PROJECT_ROOT="$(_find_here "${_start}")" || true
[[ -n "${PROJECT_ROOT}" ]] || { echo "ERROR: no .here anchor found walking up from ${_start}"; exit 1; }

module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
source "${PROJECT_ROOT}/scripts/config/config_paths.sh"
conda activate "${ENV_PREFIX}"

CCC_R="${PROJECT_ROOT}/scripts/05_ccc/02_run_cellchat.R"
TASKS="${PROJECT_ROOT}/results/tables/05_ccc/split_half/tasks.csv"
[[ -f "${CCC_R}" ]] || { echo "ERROR: 02_run_cellchat.R not found at ${CCC_R}"; exit 1; }
[[ -f "${TASKS}"  ]] || { echo "ERROR: tasks.csv not found -- run 06_split_half_reliability.R first"; exit 1; }

# row N of tasks.csv (header is row 1). Columns: row,dataset,sample,half,n_cells,subset_file
LINE=$(awk -F, -v r="${SLURM_ARRAY_TASK_ID}" 'NR>1 && $1==r {print; exit}' "${TASKS}")
[[ -n "${LINE}" ]] || { echo "ERROR: no task row ${SLURM_ARRAY_TASK_ID} in ${TASKS}"; exit 1; }
DS=$(   cut -d, -f2 <<< "${LINE}")
SMP=$(  cut -d, -f3 <<< "${LINE}")
HALF=$( cut -d, -f4 <<< "${LINE}")
SUBSET=$(cut -d, -f6 <<< "${LINE}")

echo "TASK ${SLURM_ARRAY_TASK_ID}: ${DS}/${SMP} half ${HALF} <- ${SUBSET}"
srun --export=ALL -n 1 Rscript "${CCC_R}" \
     --dataset="${DS}" --sample="${SMP}" \
     --cell_subset="${SUBSET}" --out_suffix="split_${HALF}" \
     --nboot=100 --workers=1
