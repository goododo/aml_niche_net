#!/bin/bash
# =============================================================================
# submit_J2d_doublet_calib.sbatch
# Measures scDblFinder vs DoubletFinder separately across the library-size range,
# with the package-default dbr.sd (task 1) and with dbr.sd=0 (task 2), which
# forces scDblFinder to honour the expected rate instead of treating it as a soft
# prior. Decides union-vs-intersection on measurement rather than on preference.
# =============================================================================
#SBATCH --job-name=J2d_dblcalib
#SBATCH -p gr10634c
#SBATCH --array=1-2
#SBATCH --rsc p=1:t=1:c=4:m=256G
#SBATCH --time=12:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J2d_dblcalib_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J2d_dblcalib_%A_%a.err
set -euo pipefail
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
cd /FAST/gr10634/gaozy/aml_niche_net
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then SD="-1"; TAG="default"; else SD="0"; TAG="sd0"; fi
echo "=== doublet calibration, dbr.sd=$TAG  $(date -Is)"
Rscript scripts/01_preprocess/07_doublet_calibration.R --n_per_bin 3 --dbr_sd "$SD"
mv results/tables/01_preprocess/07_doublet_calibration.csv \
   "results/tables/01_preprocess/07_doublet_calibration__${TAG}_predbl.csv"
mv results/figures/01_preprocess/07_doublet_calibration.png \
   "results/figures/01_preprocess/07_doublet_calibration__${TAG}_predbl.png"
echo "=== finished $(date -Is)"
