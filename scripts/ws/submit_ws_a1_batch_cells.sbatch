#!/bin/bash
#SBATCH -J ws_a1_batch_cells
#SBATCH -p gr10634c
#SBATCH -t 08:00:00
#SBATCH --rsc p=1:t=1:c=4:m=128G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
#SBATCH --no-requeue
# WP-A1 batch amendment (PREREGISTRATION_A1_batch.md, tag prereg-a1-batch-v1): W4 doublet scores and D9 soup profile (scripts/ws/38).
#   smoke:  sbatch scripts/ws/submit_ws_a1_batch_cells.sbatch --smoke        real: sbatch scripts/ws/submit_ws_a1_batch_cells.sbatch
set -eu -o pipefail
export R_FUTURE_GLOBALS_MAXSIZE=68719476736 PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/ws/38_a1_batch_cells.R "$@"
