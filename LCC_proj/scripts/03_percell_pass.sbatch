#!/bin/bash
# 03_percell_pass.sbatch ----
# SLURM job array: one task per row of LCC_proj/results/tables/03_sample_manifest.csv (220 samples).
# Each task runs
#   Rscript LCC_proj/scripts/03_percell_pass.R --row=$SLURM_ARRAY_TASK_ID --workers=2
# The runner is resume-safe (skips a sample whose four outputs all exist), so re-submitting only
# fills gaps and --requeue is safe.
#
# Submit (manifest-driven count; no hard-coded 220 that goes stale):
#   mkdir -p LCC_proj/logs
#   N=$(( $(wc -l < LCC_proj/results/tables/03_sample_manifest.csv) - 1 ))
#   sbatch --array=1-${N}%8 LCC_proj/scripts/03_percell_pass.sbatch
#
# Modeled on scripts/05_ccc/02_run_cellchat.sbatch: same --rsc syntax, same .here-anchor walk from
# ${SLURM_SUBMIT_DIR}, same module load + conda activate pattern.
#
#SBATCH -J lcc_scan
#SBATCH -p gr10634c
#SBATCH -t 02:00:00
#SBATCH --rsc p=1:t=2:c=2:m=36G
# Timing measured on Chen2023/AML320_NICHE_IMMUNE (7100 cells): 7.5 min with xz compression and
# single-core UCell; ~2-3 min after switching to gzip and --workers=2. 02:00:00 covers the largest
# GSE185381 libraries with wide margin.
#SBATCH -o LCC_proj/logs/%x.%A_%a.out
#SBATCH -e LCC_proj/logs/%x.%A_%a.err
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

SCAN_R="${PROJECT_ROOT}/LCC_proj/scripts/03_percell_pass.R"
[[ -f "${SCAN_R}" ]] || { echo "ERROR: 03_percell_pass.R not found at ${SCAN_R}"; exit 1; }

echo "TASK ${SLURM_ARRAY_TASK_ID}: LCC per-cell pass for manifest row ${SLURM_ARRAY_TASK_ID}"
srun --export=ALL -n 1 Rscript "${SCAN_R}" --row="${SLURM_ARRAY_TASK_ID}" --workers=2
