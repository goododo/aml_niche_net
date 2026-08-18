#!/bin/bash
# =============================================================================
# submit_cytotrace2.sbatch
#
# CytoTRACE2 over the cohort, one array task per slice of the roster.
#
# WHY SLURM AND NOT A BACKGROUND SHELL. Run on the shared login node this job was
# SIGKILLed 32 times in 8 retry rounds and stalled at 94/214 -- not because of its own
# footprint (~4 GB per worker, and the node showed 399 GB free) but because a 137-user
# interactive node reaps long-running processes under other people's peaks. The retries
# made progress precisely because 11_cytotrace2_score.R skips samples whose output
# already exists, but "make progress until something kills you" is not a run.
#
# SINGLE-THREADED ON PURPOSE. cytotrace2 forks its model and smoothing workers, and the
# conda R links a threaded BLAS; forking from inside a threaded region deadlocks. Five
# concurrent slices x 12 forked workers each sat at 0.1% CPU with 385 GB free -- hung,
# and hung in a way the log renders as "still running". One core per task cannot deadlock,
# and the OMP/OPENBLAS/MKL exports below keep the BLAS out of it.
#
#   sbatch scripts/03_hierarchy/submit_cytotrace2.sbatch
# =============================================================================
#SBATCH --job-name=ct2
#SBATCH -p gr10634c
#SBATCH --array=0-7
#SBATCH --rsc p=1:t=1:c=2:m=48G
#SBATCH --time=24:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/ct2_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/ct2_%A_%a.err

set -uo pipefail
cd /FAST/gr10634/gaozy/aml_niche_net
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

/FAST/gr10634/gaozy/general_env/bin/Rscript scripts/03_hierarchy/11_cytotrace2_score.R \
  --offset "${SLURM_ARRAY_TASK_ID}" --stride 8
