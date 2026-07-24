#!/bin/bash
# 11_submit_array.slurm ----
# SLURM job array: one task per sample, runs the full BAM chain via 10_prealign_one_sample.sh.
# Submit with the array range + concurrency cap, passing DATASET and the sample sheet:
#
#   mkdir -p logs
#   N=$(grep -vcE '^\s*#|^\s*$' samples_GSE227903.tsv)     # number of samples
#   sbatch --array=1-${N}%3 scripts/02_malignancy/11_submit_array.slurm GSE227903 samples_GSE227903.tsv
#
# %3 = at most 3 samples run at once (3 x t16/300G is polite on the shared node).
#
# 3a-1 rewire (no logic change): the worker script is located via the .here anchor walked up
#   from ${SLURM_SUBMIT_DIR} (robust to submitting from the project root OR any subdir, and to
#   SLURM spooling $0), and is now 10_prealign_one_sample.sh under scripts/02_malignancy (was:
#   c10_run_one_sample.sh assumed a sibling of the submit dir). This wrapper uses no config vars,
#   so it does the walk but does NOT source config (the worker sources it). SBATCH directives,
#   resources, and task logic are byte-preserved.
#
#SBATCH -J bamchain
#SBATCH -p gr10634c
#SBATCH -t 72:00:00
#SBATCH --rsc p=1:t=80:c=106:m=1200G
# NOTE on resources:
#   -t 72:00:00 is PER TASK (one sample), not the whole array. If the QOS MaxWall is
#      lower (check: sacctmgr show qos format=Name,MaxWall -p), drop this to the cap.
#      A single sample (download + STARsolo + cellsnp + phasing + numbat) typically
#      finishes in a few hours, so even -t 24:00:00 is usually plenty.
#   m=80G: STAR loads the hg38 index (~30G) + sort buffers; Numbat peaks ~4G. 80G is
#      ample and avoids the scheduler bumping cores (m=300G forced c=16 -> 17).
#SBATCH -o logs/%x.%A_%a.out
#SBATCH -e logs/%x.%A_%a.err
#SBATCH --no-requeue

set -eux -o pipefail
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
       NUMEXPR_NUM_THREADS=16 TBB_NUM_THREADS=16

DATASET="${1:?pass DATASET, e.g. GSE227903}"
SHEET="${2:?pass sample sheet, e.g. samples_GSE227903.tsv}"

# Per-task general scratch on the big disk (cleaned on exit). NOTE: STAR's own FIFO
# tmp is set separately to node-local /tmp in config (Lustre may not support FIFOs).
# TODO(config-completion): this /LARGE1/gr10634/gaozy/tmp base is hard-coded (not under
#   LARGE1_DIR, so not directly derivable); lift to a SCRATCH_TMP_BASE config var later.
export TMPDIR="/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TMPDIR}"; trap 'rm -rf "${TMPDIR}"' EXIT

module load miniforge3

# SLURM copies THIS batch script to a spool dir, so $0/dirname point there, NOT to where the
# worker lives. Locate the project root via the .here anchor, walked up from the submit dir
# (works whether you sbatch from the project root or any subdir inside the tree).
_find_here() {
  local d; d="$(cd "${1:-$PWD}" 2>/dev/null && pwd)" || return 1
  while [[ "$d" != "/" ]]; do
    [[ -e "$d/.here" ]] && { printf '%s\n' "$d"; return 0; }
    d="$(dirname "$d")"
  done; return 1; }
_start="${SLURM_SUBMIT_DIR:-$PWD}"
PROJECT_ROOT="$(_find_here "${_start}")" || true
[[ -n "${PROJECT_ROOT}" ]] || { echo "ERROR: no .here anchor found walking up from ${_start} (sbatch from inside the project tree)"; exit 1; }
MAL_DIR="${PROJECT_ROOT}/scripts/02_malignancy"
[[ -f "${MAL_DIR}/10_prealign_one_sample.sh" ]] || {
  echo "ERROR: 10_prealign_one_sample.sh not found in ${MAL_DIR}"; exit 1; }

# This task's sample = the Nth non-comment, non-blank line of the sheet.
LINE="$(grep -vE '^\s*#|^\s*$' "${SHEET}" | sed -n "${SLURM_ARRAY_TASK_ID}p")"
SID="$(echo "${LINE}" | cut -f1)"
SRR="$(echo "${LINE}" | cut -f2)"
echo "TASK ${SLURM_ARRAY_TASK_ID}: dataset=${DATASET} sample=${SID} acc=${SRR}"
[[ -n "${SID}" && -n "${SRR}" ]] || { echo "ERROR: empty sample/SRR for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

srun --export=ALL -n 1 bash "${MAL_DIR}/10_prealign_one_sample.sh" "${DATASET}" "${SID}" "${SRR}"
