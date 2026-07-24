#!/bin/bash
# =============================================================================
# submit_03_per_sample_qc.sbatch
# SLURM array: one task per dataset (13 datasets -> array 1-13).
# Each task picks its dataset via SLURM_ARRAY_TASK_ID inside 03_per_sample_qc.R.
#
# Recommended workflow:
#   1) First CHECK serially on one dataset (login/interactive or a single job):
#        Rscript scripts/01_preprocess/03_per_sample_qc.R --dataset Chen2023
#   2) Then submit the array for throughput:
#        sbatch submit_03_per_sample_qc.sbatch
#   3) After all tasks finish, rebuild the combined report:
#        Rscript -e 'library(here); source(here::here("scripts","config","config_paths.R"));
#          library(data.table);
#          f <- list.files(DIR_PREPROCESS, "^03_qc_report__(?!ALL).*csv$", full.names=TRUE, perl=TRUE);
#          fwrite(rbindlist(lapply(f, fread), fill=TRUE),
#                 file.path(DIR_PREPROCESS, "03_qc_report__ALL.csv"))'
# =============================================================================
#SBATCH --job-name=01_aml_qc
#SBATCH -p gr10634c
#SBATCH --array=1-13

#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=72:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/01_qc_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/01_qc_%A_%a.err

set -euo pipefail
#mkdir -p logs

# --- adjust module / conda environment to your cluster ---
module load miniforge3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /FAST/gr10634/gaozy/general_env

# DoubletFinder's pK sweep is the slow part; --cpus-per-task lets you set
# cfg_parallel_workers > 1 inside 03_per_sample_qc.R to parallelise samples.
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/01_preprocess/03_per_sample_qc.R
