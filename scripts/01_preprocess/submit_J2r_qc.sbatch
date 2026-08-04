#!/bin/bash
# =============================================================================
# submit_J2r_qc.sbatch   (J2r -- targeted per-sample QC re-run)
# Re-runs 03_per_sample_qc.R for the three re-ingested datasets only.
# GSE185991 is here because its QC FAILED outright in the first pass (demux
# prefilter emptied it); the other two need their per-sample objects refreshed so
# they carry the corrected Timepoint / Tissue.
# =============================================================================
#SBATCH --job-name=J2r_qc
#SBATCH -p gr10634c
#SBATCH --array=1-3
#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=72:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J2r_qc_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J2r_qc_%A_%a.err
set -euo pipefail
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
cd /FAST/gr10634/gaozy/aml_niche_net
DATASETS=(GSE116256 GSE185991 GSE239721)
DS="${DATASETS[$((SLURM_ARRAY_TASK_ID - 1))]}"
echo "=== per-sample QC $DS  $(date -Is)"
Rscript scripts/01_preprocess/03_per_sample_qc.R --dataset "$DS"
echo "=== finished $(date -Is)"
