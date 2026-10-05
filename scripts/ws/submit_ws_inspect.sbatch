#!/bin/bash
#SBATCH --job-name=ws_inspect
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=4:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/ws_inspect_%j.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/ws_inspect_%j.err
set -euo pipefail
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/ws/90_inspect_atlas_object.R
