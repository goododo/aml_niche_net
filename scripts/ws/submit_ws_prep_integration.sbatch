#!/bin/bash
#SBATCH --job-name=ws_prep_inte
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=8:m=1024G
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/ws_prep_inte_%j.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/ws_prep_inte_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 R_LIBS_USER=
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
echo "=== prepare re-integration input  $(date -Is)"
Rscript scripts/ws/03_prepare_integration_input.R
echo "=== finished $(date -Is)"
