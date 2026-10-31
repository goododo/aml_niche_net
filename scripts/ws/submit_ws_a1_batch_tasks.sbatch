#!/bin/bash
#SBATCH -J ws_a1_batch_task
#SBATCH -p gr10634c
#SBATCH -t 04:00:00
#SBATCH --rsc p=1:t=1:c=2:m=24G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.err
#SBATCH --no-requeue
# WP-A1 batch amendment, scripts/ws/37 --task=N: one slice of the reference fits (jobs in /LARGE1/.../ws_a1_batch/jobs.csv).
#   sbatch --array=1-100 scripts/ws/submit_ws_a1_batch_tasks.sbatch        (after --prep and --arms)
set -eu -o pipefail
export R_FUTURE_GLOBALS_MAXSIZE=68719476736 PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/ws/37_a1_batch.R --task="${SLURM_ARRAY_TASK_ID:-1}" "$@"
