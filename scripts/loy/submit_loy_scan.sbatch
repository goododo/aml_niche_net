#!/bin/bash
#SBATCH -J loy_scan
#SBATCH -p gr10634c
#SBATCH -t 03:00:00
#SBATCH --rsc p=1:t=1:c=2:m=48G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.err
#SBATCH --no-requeue
# LOY step 0: one dataset per array task; the task index is a line of results/tables/loy/01_scan/datasets.txt.
#   sbatch --array=1-$(wc -l < results/tables/loy/01_scan/datasets.txt) scripts/loy/submit_loy_scan.sbatch
set -eu -o pipefail
export R_FUTURE_GLOBALS_MAXSIZE=68719476736 PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 OMP_PROC_BIND=false
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/loy/01_scan_y_counts.R --task="${SLURM_ARRAY_TASK_ID:-1}" "$@"
