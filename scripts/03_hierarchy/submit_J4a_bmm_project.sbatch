#!/bin/bash
# =============================================================================
# submit_J4a_bmm_project.sbatch   (job J4a of the v2 re-run)
#
# Project every QC'd sample onto the BoneMarrowMap scaffold. Runs in ONE job
# rather than an array because the 362MB reference is loaded once and reused
# across samples; an array would pay that cost per task.
#
# WHY IT IS BEING RE-RUN. Every projection on disk dated 2026-07-16 and described
# a QC object that the v2 re-run has since replaced. 77 of them belonged to
# samples that left the cohort entirely (quarantined by 99_admin/14); the other
# 143 were computed from superseded cells. The resume guard tested existence, so
# without the freshness test added alongside this script the job would have
# reported "already projected" for all 143 and written nothing.
#
# IT GATES TWO THINGS:
#   * the per-bin false-positive rate in 96_malignancy_fpr_healthy -- the most
#     diagnostic piece of the malignancy calibration, since the v1 failure showed
#     up as FPR rising monotonically with distance from the lymphoid reference
#     (T_NK 0.036 -> HSC_MPP 0.398 -> Stromal 0.903);
#   * PERBIN_DIR, which 03_hierarchy/00_node_status.R needs to split `present`
#     from `absent_biological`. Until it exists every unmasked node stays
#     `present` by default.
#
# Set BMM_FORCE=1 to recompute projections that are already current.
#
# Usage: sbatch scripts/03_hierarchy/submit_J4a_bmm_project.sbatch [datasets]
# =============================================================================
#SBATCH --job-name=J4a_bmm
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=4:m=256G
#SBATCH --time=48:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J4a_bmm_%A.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J4a_bmm_%A.err
set -euo pipefail
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
cd /FAST/gr10634/gaozy/aml_niche_net
echo "=== 01_bmm_project  $(date -Is)"
Rscript scripts/03_hierarchy/01_bmm_project.R ${1:+--datasets "$1"}
echo "=== finished $(date -Is)"
