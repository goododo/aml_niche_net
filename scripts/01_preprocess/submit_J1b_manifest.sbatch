#!/bin/bash
# =============================================================================
# submit_J1b_manifest.sbatch   (job J1b of the v2 re-run)
#
# Runs AFTER J1 (all 13 ingests). Three serial steps, all cheap but all needing
# every ingest RDS present, which is why they cannot go inside the J1 array:
#   1. ingest_MASTER_summary.R  stack the per-dataset QC summaries
#   2. 01_dataset_roles.R       build 01_sample_role_manifest.csv -- now also
#                               carrying Tissue / n_donors_in_library /
#                               patient_resolved, which downstream stages filter on
#   3. 02_study_split.R         study-level 70/30 discovery/validation split
#
# NOTE 01_dataset_roles.R is run WITHOUT --cache on purpose: the cached
# 01_sample_metadata_raw.csv is from the v1 cohort and would silently reinstate
# the sample set this whole re-run exists to change.
#
# Usage: sbatch --dependency=afterok:<J1_JOBID> scripts/01_preprocess/submit_J1b_manifest.sbatch
# =============================================================================
#SBATCH --job-name=J1b_manifest
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=8:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J1b_manifest_%A.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J1b_manifest_%A.err

set -euo pipefail

module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env

cd /FAST/gr10634/gaozy/aml_niche_net

echo "=== [1/3] ingest_MASTER_summary  $(date -Is)"
Rscript scripts/00_ingest/ingest_MASTER_summary.R

echo "=== [2/3] 01_dataset_roles  $(date -Is)"
Rscript scripts/01_preprocess/01_dataset_roles.R

echo "=== [3/3] 02_study_split  $(date -Is)"
Rscript scripts/01_preprocess/02_study_split.R

# Self-verifying: print what the new manifest actually contains, so the log
# answers "did the cohort change the way it was supposed to" without a second run.
Rscript -e '
suppressPackageStartupMessages({library(here);library(data.table)})
source(here::here("scripts","config","config_paths.R"))
m <- fread(file.path(DIR_PREPROCESS,"01_sample_role_manifest.csv"))
cat("\n=== manifest:", nrow(m), "samples /", uniqueN(m$dataset), "datasets ===\n")
print(m[, .N, by=dataset][order(dataset)])
cat("\n--- timepoint ---\n"); print(m[, .N, by=Timepoint][order(-N)])
if ("Tissue" %in% names(m)) { cat("\n--- tissue ---\n"); print(m[, .N, by=Tissue]) }
if ("patient_resolved" %in% names(m)) {
  cat("\n--- patient_resolved=FALSE ---\n"); print(m[patient_resolved==FALSE, .(dataset,Sample)]) }
cat("\n--- HEALTHY ARM ---\n"); print(m[Timepoint=="Healthy", .N, by=dataset][order(-N)])
cat("healthy total:", m[Timepoint=="Healthy", .N], "\n")'

echo "=== finished $(date -Is)"
