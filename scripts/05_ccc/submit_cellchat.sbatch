#!/bin/bash
# =============================================================================
# submit_cellchat.sbatch -- CellChat over the CCC-eligible cohort, one array task per sample.
#
# WHY SLURM. This is the heaviest stage in the pipeline and the login node is shared by 137
# users; the CytoTRACE2 run was SIGKILLed 32 times there before being moved here.
#
# WHY IT IS BEING RE-RUN AT ALL. The previous tensor set was pinned by an existence-only resume
# guard (fixed: 02_run_cellchat.R now calls is_stale()). The manifest it was built from listed
# 220 samples of which 77 no longer exist, and omitted 71 that do. --row indexes the REBUILT
# manifest, so the array bound below must match its eligible count.
#
#   n=$(Rscript -e 'suppressMessages({library(data.table);library(here)});
#        source(here::here("scripts","config","config_paths.R"));
#        source(here::here("scripts","config","config_ccc.R"));
#        cat(fread(file.path(DIR_CCC,"ccc_sample_manifest.csv"))[ccc_eligible==TRUE,.N])')
#   sbatch --array=1-${n} scripts/05_ccc/submit_cellchat.sbatch
# =============================================================================
#SBATCH --job-name=cellchat
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=4:m=64G
#SBATCH --time=12:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/ccc_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/ccc_%A_%a.err

set -uo pipefail
cd /FAST/gr10634/gaozy/aml_niche_net
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

/FAST/gr10634/gaozy/general_env/bin/Rscript scripts/05_ccc/02_run_cellchat.R \
  --row "${SLURM_ARRAY_TASK_ID}"
