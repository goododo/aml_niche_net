#!/bin/bash
# 02_run_cellchat.sbatch ----
# SLURM job array: one task per ELIGIBLE CCC sample. Each task runs
#   Rscript scripts/05_ccc/02_run_cellchat.R --row=$SLURM_ARRAY_TASK_ID
# --row indexes into the ccc_eligible==TRUE rows of DIR_CCC/ccc_sample_manifest.csv (from 01),
# sorted (dataset, sample). Stable as long as the manifest is fixed. The runner is resume-safe
# (skips a sample whose tensor already exists), so re-submitting only fills gaps and --requeue is safe.
#
# Submit (manifest-driven count; no hard-coded 159 that goes stale):
#   mkdir -p logs
#   MAN=results/tables/05_ccc/ccc_sample_manifest.csv
#   N=$(awk -F, 'NR>1 && $3=="TRUE"{c++} END{print c}' "${MAN}")     # ccc_eligible is column 3
#   echo "eligible = ${N}"
#   sbatch --array=1-${N}%8 scripts/05_ccc/02_run_cellchat.sbatch
#   # %8 = at most 8 samples at once (8 x t2/36G is light on the shared node; CCC tasks are ~2-4 min).
#
# Modeled byte-for-byte on 02_malignancy/11_submit_array.slurm: same #SBATCH --rsc resource syntax,
# same .here-anchor walk from ${SLURM_SUBMIT_DIR}, same module load + conda activate pattern.
#
#SBATCH -J ccc_cellchat
#SBATCH -p gr10634c
#SBATCH -t 02:00:00
#SBATCH --rsc p=1:t=2:c=2:m=36G
# NOTE on resources:
#   -t 02:00:00 is PER TASK (one sample). A sample (load + normalize + CellChat nboot=100) runs in
#      ~2-4 min; 02:00:00 is generous headroom. Lower to the QOS MaxWall cap if needed
#      (check: sacctmgr show qos format=Name,MaxWall -p).
#   m=36G (2 cores @ ~18G/core): an 8k-14k-cell CellChat peaks well under this; 36G is safety margin
#      against the largest samples and avoids the scheduler bumping cores.
#SBATCH -o logs/%x.%A_%a.out
#SBATCH -e logs/%x.%A_%a.err
#SBATCH --requeue

set -eux -o pipefail
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
# R BLAS threads -> match allocated cores (2).
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
       NUMEXPR_NUM_THREADS=2 TBB_NUM_THREADS=2

# Per-task scratch for R tempfiles (cleaned on exit). Mirrors the BAM-chain template's tmp hygiene.
export TMPDIR="/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}"; trap 'rm -rf "${TMPDIR}"' EXIT

# SLURM spools THIS script to a temp dir, so $0/dirname point there. Locate the project root via the
# .here anchor, walked up from the submit dir (works from the project root OR any subdir inside).
_find_here() {
  local d; d="$(cd "${1:-$PWD}" 2>/dev/null && pwd)" || return 1
  while [[ "$d" != "/" ]]; do
    [[ -e "$d/.here" ]] && { printf '%s\n' "$d"; return 0; }
    d="$(dirname "$d")"
  done; return 1; }
_start="${SLURM_SUBMIT_DIR:-$PWD}"
PROJECT_ROOT="$(_find_here "${_start}")" || true
[[ -n "${PROJECT_ROOT}" ]] || { echo "ERROR: no .here anchor found walking up from ${_start} (sbatch from inside the project tree)"; exit 1; }

# Activate the R env: source config_paths.sh for ENV_PREFIX, then conda activate (vet-script pattern).
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
source "${PROJECT_ROOT}/scripts/config/config_paths.sh"
conda activate "${ENV_PREFIX}"

CCC_R="${PROJECT_ROOT}/scripts/05_ccc/02_run_cellchat.R"
[[ -f "${CCC_R}" ]] || { echo "ERROR: 02_run_cellchat.R not found at ${CCC_R}"; exit 1; }

echo "TASK ${SLURM_ARRAY_TASK_ID}: running CellChat for eligible sample row ${SLURM_ARRAY_TASK_ID}"
srun --export=ALL -n 1 Rscript "${CCC_R}" --row="${SLURM_ARRAY_TASK_ID}" --nboot=100 --workers=1
