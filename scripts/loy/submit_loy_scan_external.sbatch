#!/bin/bash
# submit_loy_scan_external.sbatch ----
# input : results/tables/loy/01_scan/external_datasets.txt (one Zeng-pack dataset id per line; array index = line)
# output: per-cell Y/XIST counts for the external pack via scripts/loy/03_scan_external.R; logs in 00_project/logs
# what it does: array wrapper, one task per study; Lasry (7.2 GB RDS) and Lambo (5.3 GB) need the 48 GB slot.
# Usage : sbatch --array=1-22 scripts/loy/submit_loy_scan_external.sbatch      (or a subset, e.g. --array=21)
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=2:m=48G
#SBATCH -t 03:00:00
#SBATCH -J loy_scan_ext
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.err
source ~/.bashrc; conda activate general_env; export OMP_PROC_BIND=false
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/loy/03_scan_external.R --task=${SLURM_ARRAY_TASK_ID}
