#!/bin/bash
#SBATCH -J ws_anno
#SBATCH -p gr10634c
#SBATCH -t 12:00:00
#SBATCH --rsc p=1:t=1:c=8:m=1024G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
#SBATCH --no-requeue
# Marker-score annotation at res 1.0, following 03.2_Anno_pipeline_260206.r. Deletes nothing.
set -eu -o pipefail
export R_FUTURE_GLOBALS_MAXSIZE=68719476736   # 64 GiB -- see submit_ws_rpca.sbatch
export PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/ws/09_annotate_clusters.R
