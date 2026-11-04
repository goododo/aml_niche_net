#!/bin/bash
#SBATCH -J ws_a1_g04
#SBATCH -p gr10634c
#SBATCH -t 12:00:00
#SBATCH --rsc p=1:t=1:c=4:m=64G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
#SBATCH --no-requeue
# WP-A1 G04 (PREREGISTRATION_A1_G04.md, tag prereg-a1-g04-v1): one mode of scripts/ws/41 per job.
#   sbatch scripts/ws/submit_ws_a1_g04.sbatch --build | --known | --primary | --sens | --readouts | --children | --combine [--smoke]
#   sbatch scripts/ws/submit_ws_a1_g04.sbatch smoke-all      (every written mode on synthetic data, in order)
# Each job runs a dated copy of the script, so an edit on disk cannot land in a running job; the copy's md5 is in every output.
set -eu -o pipefail
export R_FUTURE_GLOBALS_MAXSIZE=68719476736 PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1   # parallel work is forked by mclapply instead
export OMP_PROC_BIND=false   # the job environment sets OMP_PROC_BIND=spread, which pins R (and every forked worker) to one core
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-4}
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
CODE=/LARGE1/gr10634/gaozy/aml_niche_net/ws_a1_g04/code; mkdir -p "$CODE"
run() {   # run [args]: copy the script with a date and job stamp, then run the copy
  local c="$CODE/41_a1_g04_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID}.R"; cp scripts/ws/41_a1_g04.R "$c"; chmod a-w "$c"
  echo "running $c (md5 $(md5sum "$c" | cut -d' ' -f1)) $*; $(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc) CPUs allowed to this shell"; Rscript "$c" "$@"; }
if [ "${1:-}" = "smoke-all" ]; then
  for m in --build --known --primary --sens --readouts --children --combine; do run "$m" --smoke; done
else
  run "$@"
fi
