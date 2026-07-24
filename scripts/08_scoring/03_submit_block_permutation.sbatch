#!/bin/bash
# 03_submit_block_permutation.sbatch ----
# Block-level permutation run for 08_scoring/03_block_permutation.py (~83 min at n_perm=10000, same
# cost as 02: each permutation rebuilds two barycenters; only the statistic changed, 49 edges -> 5 blocks).
#
# Submit:
#   mkdir -p logs
#   sbatch scripts/08_scoring/03_submit_block_permutation.sbatch          # defaults 10000 300
#   sbatch scripts/08_scoring/03_submit_block_permutation.sbatch 2000 300 # quicker check
#
#SBATCH -J ccc_block
#SBATCH -p gr10634c
#SBATCH -t 24:00:00
#SBATCH --rsc p=1:t=4:c=4:m=72G
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err
#SBATCH --no-requeue
# --no-requeue: the permutation loop has no mid-run checkpoint; a requeue would silently restart at 0.

set -eux -o pipefail
# Do NOT export PYTHONNOUSERSITE here (see 02_submit_permutation.sbatch): dependencies live partly in
# the conda env and partly in ~/.local, and both must stay visible.
export PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 \
       NUMEXPR_NUM_THREADS=4 TBB_NUM_THREADS=4

N_PERM="${1:-10000}"
MAX_ITER="${2:-300}"

export TMPDIR="/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}"
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
conda activate "${ENV_PREFIX}" || echo "WARN: conda activate returned nonzero; using absolute interpreter"

# Resolve the interpreter explicitly (conda activate does not reliably win the PATH race in batch).
PYBIN="${ENV_PREFIX}/bin/python"
[[ -x "${PYBIN}" ]] || { echo "ERROR: interpreter not executable: ${PYBIN}"; exit 1; }
PY="${PROJECT_ROOT}/scripts/08_scoring/03_block_permutation.py"
[[ -f "${PY}" ]] || { echo "ERROR: 03_block_permutation.py not found at ${PY}"; exit 1; }

echo "--- env preflight ---"
echo "PYBIN = ${PYBIN}"
"${PYBIN}" -c "import numpy,pandas,ot; print('numpy :',numpy.__version__,numpy.__file__); print('pandas:',pandas.__version__); print('POT   :',ot.__version__,ot.__file__)"
echo "---------------------"

echo "BLOCK PERMUTATION RUN: n_perm=${N_PERM} max_iter=${MAX_ITER}"
srun --export=ALL -n 1 "${PYBIN}" "${PY}" --n_perm "${N_PERM}" --max_iter "${MAX_ITER}"
