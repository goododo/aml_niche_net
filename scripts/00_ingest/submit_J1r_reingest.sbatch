#!/bin/bash
# =============================================================================
# submit_J1r_reingest.sbatch   (J1r -- targeted re-ingest)
#
# Re-ingests ONLY the three datasets affected by the defects fixed in 1604df5:
#   GSE116256   19 samples carried the retired "Post_treatment" value
#   GSE185991   PB draw mislabelled BM (TISSUE_OVERRIDE keyed on the curated id)
#   GSE239721   Tissue blank for all 20 (ingest bypasses make_seurat)
# All three defects live in metadata written AT INGEST, so the RDS must be
# rebuilt; the other 10 datasets are unaffected and are not touched.
#
# Depends on afterANY (not afterok) of the J2 QC array: one of its tasks failed
# by design of the bug above, and the only thing that matters here is that no QC
# task is still READING these RDS files when we overwrite them.
# =============================================================================
#SBATCH --job-name=J1r_reingest
#SBATCH -p gr10634c
#SBATCH --array=1-3
#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=12:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J1r_reingest_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J1r_reingest_%A_%a.err
set -euo pipefail
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
cd /FAST/gr10634/gaozy/aml_niche_net
DATASETS=(GSE116256 GSE185991 GSE239721)
DS="${DATASETS[$((SLURM_ARRAY_TASK_ID - 1))]}"
echo "=== re-ingest $DS  $(date -Is)"
Rscript "scripts/00_ingest/ingest_${DS}.R"
echo "=== finished $(date -Is)"
