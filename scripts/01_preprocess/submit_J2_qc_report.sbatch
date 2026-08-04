#!/bin/bash
# =============================================================================
# submit_J2_qc_report.sbatch   (job J2c of the v2 re-run)
#
# Runs AFTER the 03_per_sample_qc array. Combines the per-dataset QC reports,
# verifies the objects landed, draws the cohort QC figures, and then runs the
# metadata audit chain ending in the BLOCKING roster gate.
#
#   1. combine 03_qc_report__<ds>.csv -> 03_qc_report__ALL.csv
#   2. 05_qc_landing_check.R    objects on LARGE1 vs report rows on FAST
#   3. 06_qc_figures.R          Q1-Q5 cohort QC panels
#   4. 99_admin/10 -> 11 -> 12  skeleton, join+label audit, roster diff
#   5. 99_admin/12 --post       EXITS NON-ZERO if the roster does not match
#                               curated with every difference diagnosed
#
# Step 5 is the point of the whole re-run: it makes "the sample set is correct"
# a checked property rather than a claim.
#
# Usage: sbatch --dependency=afterok:<J2_ARRAY_JOBID> scripts/01_preprocess/submit_J2_qc_report.sbatch
# =============================================================================
#SBATCH --job-name=J2c_qcreport
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=4:m=256G
#SBATCH --time=4:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J2c_qcreport_%A.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J2c_qcreport_%A.err

set -euo pipefail

module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env

cd /FAST/gr10634/gaozy/aml_niche_net

echo "=== [1/5] combine per-dataset QC reports  $(date -Is)"
Rscript -e '
suppressPackageStartupMessages({library(here);library(data.table)})
source(here::here("scripts","config","config_paths.R"))
source(here::here("scripts","config","config_qc.R"))
source(here::here("scripts","config","utils.R"))
f <- list.files(DIR_PREPROCESS, "^03_qc_report__(?!ALL).*csv$", full.names=TRUE, perl=TRUE)
if (!length(f)) stop("no per-dataset 03_qc_report__*.csv found")
a <- rbindlist(lapply(f, fread), fill=TRUE)
fwrite_safe(a, file.path(DIR_PREPROCESS,"03_qc_report__ALL.csv"))
cat(sprintf("combined %d dataset reports -> %d samples\n", length(f), nrow(a)))
print(a[, .(samples=.N, pass=sum(status=="PASS", na.rm=TRUE),
            cells_final=sum(n_final, na.rm=TRUE)), by=dataset][order(dataset)])
if (length(f) != 13) stop(sprintf("expected 13 dataset reports, got %d -- an array task failed", length(f)))'

echo "=== [2/5] 05_qc_landing_check  $(date -Is)"
Rscript scripts/01_preprocess/05_qc_landing_check.R

echo "=== [3/5] 06_qc_figures  $(date -Is)"
Rscript scripts/01_preprocess/06_qc_figures.R

echo "=== [4/5] metadata audit chain  $(date -Is)"
Rscript scripts/99_admin/10_build_metadata_skeleton.R
Rscript scripts/99_admin/11_audit_labels_vs_curated.R
Rscript scripts/99_admin/12_cohort_roster_check.R

echo "=== [5/5] ROSTER GATE (blocking)  $(date -Is)"
Rscript scripts/99_admin/12_cohort_roster_check.R --post

echo "=== finished $(date -Is)"
