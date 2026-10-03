#!/bin/bash
# =============================================================================
# submit_ws_T3a_units.sbatch
#
# Three independent tasks of HANDOFF_within_sample_v1.5, submitted together so
# they survive the interactive session:
#   1  T3a   per-library doublet scores, both callers (the long one: DoubletFinder's
#            pK sweep runs once per library, ~95 libraries)
#   2        GATE_W1.2b recounted with bmm_fine (29 projected states) as the state unit
#   3        GATE_W1.2b recounted with top_MP (10 data-driven malignant meta-programs)
# Tasks 2 and 3 answer the objection that a 7-bin receiver is not one state; they
# apply the SAME frozen counting rule to a finer and to a data-driven definition.
# Neither can overwrite the production tables: a non-default --unit requires --suffix.
# =============================================================================
#SBATCH --job-name=ws_T3a_units
#SBATCH -p gr10634c
#SBATCH --array=1-3
#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/ws_T3a_units_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/ws_T3a_units_%A_%a.err
set -euo pipefail
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
cd /FAST/gr10634/gaozy/aml_niche_net
echo "=== task $SLURM_ARRAY_TASK_ID start $(date -Is)"
case "$SLURM_ARRAY_TASK_ID" in
  1) Rscript scripts/ws/02_doublet_scores.R ;;
  2) Rscript scripts/ws/01_samebin_counts.R --unit=bmm_fine --suffix=bmmfine ;;
  3) Rscript scripts/ws/01_samebin_counts.R --unit=top_MP   --suffix=topMP ;;
esac
echo "=== task $SLURM_ARRAY_TASK_ID finished $(date -Is)"
