#!/bin/bash
# submit_loy_h2.sbatch ----
# input : step number 04 | 05 | 06 as $1, remaining args passed to the R script
# output: results/tables/loy/04_h2/* (per step); log 00_project/logs/loy_h2_<jobid>.{out,err}
# what it does: runs one H2 step of the LOY line (scripts/loy/PREREGISTRATION_H2.md) on a compute node.
# Usage : sbatch scripts/loy/submit_loy_h2.sbatch 04            # Stage 0: counts, matching, K3, gene-set freeze
#         sbatch scripts/loy/submit_loy_h2.sbatch 05 --world null   # etc.
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=2:m=48G
#SBATCH -t 04:00:00
#SBATCH -J loy_h2
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
source ~/.bashrc; conda activate general_env; export OMP_PROC_BIND=false
cd /FAST/gr10634/gaozy/aml_niche_net
STEP=${1:?usage: sbatch scripts/loy/submit_loy_h2.sbatch 04|05|06 [args]}; shift
SCRIPT=$(ls scripts/loy/${STEP}_h2_*.R | head -1)
[ -n "$SCRIPT" ] || { echo "no script for step $STEP"; exit 1; }
echo "[job] $SLURM_JOB_ID $SCRIPT $* $(date -Iseconds)"
Rscript "$SCRIPT" "$@"
