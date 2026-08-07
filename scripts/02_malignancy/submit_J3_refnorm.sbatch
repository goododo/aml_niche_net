#!/bin/bash
# =============================================================================
# submit_J3_refnorm.sbatch   (job J3 of the v2 re-run)
#
# Phase 1.1 of malignancy: find each sample's AUTOLOGOUS normal reference (high-
# confidence mature T/B cells, blast-gated) via SingleR against BoneMarrowMap.
# This anchors the inferCNV baseline; it must be redone because the per-sample QC
# objects it reads were all rewritten by the v2 re-run (new cohort, absolute
# bounds now wired, dbr.sd=0 doublet setting).
#
# It also now skips samples with patient_resolved=FALSE: inferCNV's reference
# subtraction assumes ONE genome, and on a 2-donor chimera it returns a smooth,
# plausible average of two people rather than failing.
#
# DOWNSTREAM IS GATED. 40/44_infercnv_run.R must NOT start until the malignancy
# calibration is rechecked: on the v1 calls healthy median frac_malignant (0.135)
# was HIGHER than diagnosis (0.115), AUC 0.546, because REF_KEEP_LABELS is
# lymphoid-only and the false-positive rate rises monotonically with
# transcriptional distance from lymphocytes (T_NK 0.036 -> HSC_MPP 0.398 ->
# Stromal 0.903). J0 rebuilt the lineage-matched external reference; whether that
# fixes the calibration is an empirical question this job's successor must answer.
#
# Usage: sbatch scripts/02_malignancy/submit_J3_refnorm.sbatch
# =============================================================================
#SBATCH --job-name=J3_refnorm
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=72:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J3_refnorm_%A.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J3_refnorm_%A.err
set -euo pipefail
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
cd /FAST/gr10634/gaozy/aml_niche_net
echo "=== 20_refnorm_identify  $(date -Is)"
# Optional dataset filter, forwarded to --datasets. A subset run UPDATES the summary
# in place (rows for the samples it processed are replaced, the rest carried through),
# so it cannot silently shrink inferCNV's task list.
#   sbatch submit_J3_refnorm.sbatch E-MTAB-11536,GSE253355
Rscript scripts/02_malignancy/20_refnorm_identify.R ${1:+--datasets "$1"}
echo "=== finished $(date -Is)"
