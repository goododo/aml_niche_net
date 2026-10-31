#!/bin/bash
#SBATCH -J ws_a1_batch
#SBATCH -p gr10634c
#SBATCH -t 12:00:00
#SBATCH --rsc p=1:t=1:c=2:m=48G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
#SBATCH --no-requeue
# WP-A1 batch amendment (PREREGISTRATION_A1_batch.md, tag prereg-a1-batch-v1), scripts/ws/37: one mode per job.
#   sbatch scripts/ws/submit_ws_a1_batch.sbatch --prep | --arms | --combine        (after scripts/ws/38 has run)
#   smoke: sbatch scripts/ws/submit_ws_a1_batch.sbatch --smoke-all
set -eu -o pipefail
export R_FUTURE_GLOBALS_MAXSIZE=68719476736 PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
if [ "${1:-}" = "--smoke-all" ]; then
  Rscript scripts/ws/37_a1_batch.R --prep --smoke
  Rscript scripts/ws/37_a1_batch.R --arms --smoke
  for t in $(awk -F, 'NR>1{print $NF}' /LARGE1/gr10634/gaozy/aml_niche_net/ws_a1_batch/jobs_smoke.csv | sort -un); do Rscript scripts/ws/37_a1_batch.R --task=$t --smoke; done
  Rscript scripts/ws/37_a1_batch.R --combine --smoke
else
  Rscript scripts/ws/37_a1_batch.R "$@"
fi
