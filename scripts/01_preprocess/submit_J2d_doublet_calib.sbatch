#!/bin/bash
# =============================================================================
# submit_J2d_doublet_calib.sbatch
#
# Measures scDblFinder vs DoubletFinder SEPARATELY across the library-size range,
# under two settings of scDblFinder's dbr.sd:
#   task 1  package default (40% of dbr)  -- dbr treated as a soft prior
#   task 2  dbr.sd = 0                    -- dbr treated as the expected rate
# so union-vs-intersection is decided on measurement rather than on preference.
#
# Each task writes its OWN filenames via --out_suffix. The previous version had
# both tasks write 07_doublet_calibration.csv and then mv it, which raced: one mv
# failed and the surviving file carried the other condition's label. A job must
# not depend on winning a rename.
#
# n_per_bin=7 over 4 size bins -> 28 samples, drawn round-robin across datasets
# within each bin so one platform cannot speak for a whole size regime.
# =============================================================================
#SBATCH --job-name=J2d_dblcalib
#SBATCH -p gr10634c
#SBATCH --array=1-2
#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=24:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J2d_dblcalib_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J2d_dblcalib_%A_%a.err
set -euo pipefail
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
cd /FAST/gr10634/gaozy/aml_niche_net
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then SD="-1"; TAG="default_n28"; else SD="0"; TAG="sd0_n28"; fi
echo "=== doublet calibration  dbr.sd=$TAG  $(date -Is)"
Rscript scripts/01_preprocess/07_doublet_calibration.R \
        --n_per_bin 7 --dbr_sd "$SD" --out_suffix "$TAG"
echo "=== finished $(date -Is)"
