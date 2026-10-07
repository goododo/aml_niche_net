#!/bin/bash
#SBATCH -J ws_cluster
#SBATCH -p gr10634c
#SBATCH -t 24:00:00
#SBATCH --rsc p=1:t=1:c=8:m=512G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
#SBATCH --no-requeue
# Resolution scan on the re-integrated object, following 03.1_Diff_reso.r.
# Reads the 250 MB PCA sidecar, not the 42 GB object. Deletes nothing.
set -eu -o pipefail
# future's globals guard defaults to 500 MiB and killed two jobs at this cell count
# (3641248 after 6h45m in IntegrateData, 3641585 after 1m in FindNeighbors). The env var
# sets the option for every script run from here, so a new script cannot forget it.
export R_FUTURE_GLOBALS_MAXSIZE=68719476736   # 64 GiB
export PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/ws/07_cluster_resolutions.R
