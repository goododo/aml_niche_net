#!/bin/bash
# submit_loy_call.sbatch ----
# input : all per-cell files from 01 (+ whatever external ones from 03 exist); output: results/tables/loy/02_calls/*
# what it does: runs 02_call_loy_status.R on a compute node (800k+ cells in memory). Pass --partial through $1 if needed.
# Usage : sbatch scripts/loy/submit_loy_call.sbatch            |  sbatch scripts/loy/submit_loy_call.sbatch --partial
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=2:m=48G
#SBATCH -t 02:00:00
#SBATCH -J loy_call
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
source ~/.bashrc; conda activate general_env; export OMP_PROC_BIND=false
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/loy/02_call_loy_status.R "$@"
