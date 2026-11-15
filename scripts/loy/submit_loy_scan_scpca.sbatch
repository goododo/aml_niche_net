#!/bin/bash
# submit_loy_scan_scpca.sbatch ----
# input : the ScPCA SCPCP000007 download (config section 7f); output: per-cell / summary / checks tables of 03c
# what it does: runs 03c_scan_scpca.R on a compute node (30 samples, BoneMarrowMap projection each).
# Usage : sbatch scripts/loy/submit_loy_scan_scpca.sbatch [--force]
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=2:m=48G
#SBATCH -t 06:00:00
#SBATCH -J loy_scan_scpca
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
source ~/.bashrc; conda activate general_env; export OMP_PROC_BIND=false
cd /FAST/gr10634/gaozy/aml_niche_net
echo "[job] $SLURM_JOB_ID $* $(date -Iseconds)"
Rscript scripts/loy/03c_scan_scpca.R "$@"
