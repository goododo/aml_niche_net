#!/bin/bash
# =============================================================================
# submit_J1_ingest.sbatch   (job J1 of the v2 re-run)
#
# Re-ingests all 13 datasets. One SLURM array task per dataset.
#
# WHY THIS RE-RUN EXISTS: results/tables/01_preprocess/ is dated 2026-07-15 while
# the config fixes land after it, so every downstream stage rests on pre-fix
# data. On top of that the COHORT ITSELF changed:
#   E-MTAB-11536  12 -> 8   (bone marrow only; 4 blood libraries dropped)
#   Petti2019      9 -> 5   (4 uncurated healthy donors dropped)
#   GSE185381     63 -> 47 libraries -> 52 donor-level samples (5' arm dropped;
#                 AML001/AML2123/PAWWEE recovered from v1's library-level gate)
#   Chen2023      20 libraries -> 10 donor-level samples (pools A+B merged)
#   GSE289435     13 -> 12  (PDX moved into INGEST_EXCLUDE)
# Healthy arm 53 -> 40.
#
# Order (do not reorder):
#   J1   this           per-dataset ingest  -> LARGE1/01_processed_counts/rds/
#   J1b  MASTER_summary combine per-dataset QC summaries
#   J2   01_dataset_roles + 02_study_split + 03_per_sample_qc + 05_landing_check
#   then 99_admin 10 -> 11 -> 12, and 12 --post as the blocking roster gate
#
# Usage:  sbatch scripts/00_ingest/submit_J1_ingest.sbatch
# Watch:  squeue -u $USER ; tail -f 00_project/logs/J1_ingest_*_<task>.err
# =============================================================================
#SBATCH --job-name=J1_ingest
#SBATCH -p gr10634c
#SBATCH --array=1-13
#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=12:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J1_ingest_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J1_ingest_%A_%a.err

set -euo pipefail

module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env

cd /FAST/gr10634/gaozy/aml_niche_net

# Sorted so the task-id -> dataset mapping is stable across submissions; a task
# that fails can then be resubmitted by id alone (--array=5).
mapfile -t SCRIPTS < <(ls scripts/00_ingest/ingest_*.R \
  | grep -v -E 'ingest_(MASTER_summary)\.R$' | sort)

N=${#SCRIPTS[@]}
if [ "$N" -ne 13 ]; then
  echo "FATAL: expected 13 ingest scripts, found $N" >&2
  printf '  %s\n' "${SCRIPTS[@]}" >&2
  exit 1
fi

SCRIPT="${SCRIPTS[$((SLURM_ARRAY_TASK_ID - 1))]}"
echo "=== task ${SLURM_ARRAY_TASK_ID}/${N} : ${SCRIPT}"
echo "=== started $(date -Is)"

Rscript "$SCRIPT"

echo "=== finished $(date -Is)"
