#!/bin/bash
#SBATCH -J ws_mixdiag
#SBATCH -p gr10634c
#SBATCH -t 6:00:00
#SBATCH --rsc p=1:t=1:c=8:m=512G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
#SBATCH --no-requeue
# Reads the integrated object once, writes a small PCA+metadata sidecar, and recomputes mixing
# within each SampleType. Deletes nothing; does not modify the integrated object.
set -eu -o pipefail
# future's globals guard defaults to 500 MiB and killed two jobs at this cell count
# (3641248 after 6h45m in IntegrateData, 3641585 after 1m in FindNeighbors). The env var
# sets the option for every script run from here, so a new script cannot forget it.
export R_FUTURE_GLOBALS_MAXSIZE=68719476736   # 64 GiB
export PYTHONNOUSERSITE=1 R_LIBS_USER=
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/ws/06_mixing_diagnostic.R
