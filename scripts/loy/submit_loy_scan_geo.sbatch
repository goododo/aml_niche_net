#!/bin/bash
# submit_loy_scan_geo.sbatch ----
# input : the LOY_GEO_MTX rows of config_loy.R section 7e (array index = row); output: per-cell / summary / checks tables of 03b
# what it does: runs 03b_scan_geo_mtx.R for one raw-matrix GEO series per array task on a compute node (BoneMarrowMap projection
#   per library: GSE223844 has ~125 libraries, hence the long wall time).
# Usage : sbatch scripts/loy/submit_loy_scan_geo.sbatch          |  sbatch --array=1 scripts/loy/submit_loy_scan_geo.sbatch --force
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=2:m=48G
#SBATCH -t 12:00:00
#SBATCH -J loy_scan_geo
#SBATCH --array=0-3
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.err
source ~/.bashrc; conda activate general_env; export OMP_PROC_BIND=false
cd /FAST/gr10634/gaozy/aml_niche_net
DS=(GSE315856 GSE223844 GSE343017 GSE198681)
echo "[job] $SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID ${DS[$SLURM_ARRAY_TASK_ID]} $* $(date -Iseconds)"
Rscript scripts/loy/03b_scan_geo_mtx.R --dataset=${DS[$SLURM_ARRAY_TASK_ID]} "$@"
