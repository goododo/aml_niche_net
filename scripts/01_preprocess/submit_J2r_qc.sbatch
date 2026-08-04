#!/bin/bash
# =============================================================================
# submit_J2r_qc.sbatch   (J2r -- targeted per-sample QC re-run)
#
# Re-runs 03_per_sample_qc.R for specific datasets, named in DATASETS below.
#
# REQUIRES the --dataset precedence fix: 03_per_sample_qc.R used to let
# SLURM_ARRAY_TASK_ID override --dataset, so an array doing its OWN task-id ->
# dataset mapping silently ran all_datasets[i] instead. The first attempt at this
# job re-ran Chen2023 / E-MTAB-11536 / GSE116256 and never touched the two
# datasets it was submitted for -- with no error, because those names are valid.
# =============================================================================
#SBATCH --job-name=J2r_qc
#SBATCH -p gr10634c
#SBATCH --array=1-1
#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=72:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J2r_qc_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J2r_qc_%A_%a.err
set -euo pipefail
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
cd /FAST/gr10634/gaozy/aml_niche_net
DATASETS=(GSE185991)
DS="${DATASETS[$((SLURM_ARRAY_TASK_ID - 1))]}"
echo "=== per-sample QC $DS  (array id ${SLURM_ARRAY_TASK_ID})  $(date -Is)"
Rscript scripts/01_preprocess/03_per_sample_qc.R --dataset "$DS"
echo "=== finished $DS  $(date -Is)"
