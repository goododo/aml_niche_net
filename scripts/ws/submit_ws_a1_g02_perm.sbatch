#!/bin/bash
#SBATCH -J ws_a1_g02_perm
#SBATCH -p gr10634c
#SBATCH -t 04:00:00
#SBATCH --rsc p=1:t=1:c=2:m=24G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%A_%a.err
#SBATCH --no-requeue
# WP-A1 gap G02, D2 (PREREGISTRATION_A1_G02_reporting.md section 4): one slice of 20 relabellings per array task (scripts/ws/36).
#   first, once:  Rscript scripts/ws/36_a1_g02_perm_crossfit.R --draw        (in the conda env; writes ws_a1_g02_perm_labels.csv)
#   then:         sbatch --array=1-50 scripts/ws/submit_ws_a1_g02_perm.sbatch
#   finally:      Rscript scripts/ws/36_a1_g02_perm_crossfit.R --combine
set -eu -o pipefail
export R_FUTURE_GLOBALS_MAXSIZE=68719476736 PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/ws/36_a1_g02_perm_crossfit.R --task="${SLURM_ARRAY_TASK_ID:-1}" "$@"
