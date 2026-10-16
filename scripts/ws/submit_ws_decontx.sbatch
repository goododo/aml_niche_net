#!/bin/bash
#SBATCH -J ws_decontx
#SBATCH -p gr10634c
#SBATCH -t 04:00:00
#SBATCH --rsc p=1:t=1:c=4:m=64G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.err
#SBATCH --no-requeue
# v2.0 S1: one DecontX unit (10x channel) per array task.
#   build the index first:  Rscript scripts/ws/22_decontx_library.R --dataset=GSE185381 --index
#   then:                   sbatch --array=1-$(($(wc -l < results/tables/ws/decontx/GSE185381/ws_decontx_units.csv)-1)) \
#                               scripts/ws/submit_ws_decontx.sbatch --dataset=GSE185381
set -eu -o pipefail
export R_FUTURE_GLOBALS_MAXSIZE=68719476736 PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/ws/22_decontx_library.R --task="${SLURM_ARRAY_TASK_ID:-1}" "$@"
