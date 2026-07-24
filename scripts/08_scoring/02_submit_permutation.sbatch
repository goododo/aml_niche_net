#!/bin/bash
# 02_submit_permutation.sbatch ----
# Full permutation run for 08_scoring/02_permutation_emergent.py. The n_perm=20 smoke test showed the
# p-value FLOOR is 1/(n_perm+1), so a 20-perm run can only ever reach p=0.0476 -- useless against a
# 49-edge BH-FDR correction (rank-1 needs p < 0.05/49 ~ 0.001). n_perm=2000 gives a floor of 0.0005,
# enough headroom for FDR. Each permutation rebuilds TWO barycenters (~20 + ~93 graphs), so this runs
# well past the 09:00-17:00 login window -> must be a batch job.
#
# Submit (from anywhere inside the project tree):
#   mkdir -p logs
#   sbatch scripts/08_scoring/02_submit_permutation.sbatch
#   # optional override: sbatch scripts/08_scoring/02_submit_permutation.sbatch 1000 300
#
# Modeled on the proven 02_malignancy/11_submit_array.slurm: same --rsc syntax, same .here-anchor walk
# from ${SLURM_SUBMIT_DIR}, same module load + conda activate pattern.
#
#SBATCH -J ccc_perm
#SBATCH -p gr10634c
#SBATCH -t 24:00:00
#SBATCH --rsc p=1:t=4:c=4:m=72G
# NOTE on resources:
#   Single-threaded numpy/POT work on 7x7 matrices -- tiny memory. c=4/m=72G is generous headroom
#   (memory couples to cores at ~18G/core on this cluster). -t 24:00:00 covers a slow 2000-perm run;
#   drop to the QOS MaxWall cap if lower (sacctmgr show qos format=Name,MaxWall -p).
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err
#SBATCH --no-requeue
# NOTE: --no-requeue because the permutation loop has no mid-run checkpoint; a requeue would restart
#   from perm 0 and silently waste the elapsed work. If the job dies, resubmit deliberately.

set -eux -o pipefail
# NOTE: the bamchain template exports PYTHONNOUSERSITE=1 (it runs STAR/samtools/R and wants a hermetic
# conda env). DO NOT do that here: numpy / pandas / POT are installed in the USER site-packages
# (~/.local/lib/python3.x/site-packages), and PYTHONNOUSERSITE=1 hides exactly that directory, which is
# why the batch job failed with ModuleNotFoundError while the identical interactive command worked.
export PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 \
       NUMEXPR_NUM_THREADS=4 TBB_NUM_THREADS=4

N_PERM="${1:-2000}"
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
conda activate "${ENV_PREFIX}" || echo "WARN: conda activate returned nonzero; continuing with absolute interpreter"

# Resolve the interpreter EXPLICITLY rather than trusting `conda activate` to win the PATH race.
# Root cause of the previous failure: in the batch shell `python` resolved to an interpreter that saw
# ~/.local (numpy, scipy, pandas) but NOT ${ENV_PREFIX}/lib/pythonX/site-packages (POT) -- so numpy
# imported and `ot` did not. The env's own interpreter sees BOTH (user-site stays enabled because we do
# not set PYTHONNOUSERSITE), which is exactly the interactive situation that already works.
PYBIN="${ENV_PREFIX}/bin/python"
[[ -x "${PYBIN}" ]] || { echo "ERROR: interpreter not executable: ${PYBIN}"; exit 1; }

PY="${PROJECT_ROOT}/scripts/08_scoring/02_permutation_emergent.py"
[[ -f "${PY}" ]] || { echo "ERROR: 02_permutation_emergent.py not found at ${PY}"; exit 1; }

# Pre-flight: prove THIS interpreter resolves every dependency before burning ~80 min of compute.
echo "--- env preflight ---"
echo "PYBIN = ${PYBIN}"
"${PYBIN}" -c "import sys; print('exe :', sys.executable); print('ver :', sys.version.split()[0]); print('path:', *sys.path, sep='\n      ')"
"${PYBIN}" -c "import numpy,pandas,ot; print('numpy :',numpy.__version__,numpy.__file__); print('pandas:',pandas.__version__); print('POT   :',ot.__version__,ot.__file__)"
echo "---------------------"

echo "PERMUTATION RUN: n_perm=${N_PERM} max_iter=${MAX_ITER}"
srun --export=ALL -n 1 "${PYBIN}" "${PY}" --n_perm "${N_PERM}" --max_iter "${MAX_ITER}"
