#!/bin/bash
#SBATCH -J ws_rpca
#SBATCH -p gr10634c
#SBATCH -t 72:00:00
#SBATCH --rsc p=1:t=58:c=58:m=1024G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
#SBATCH --no-requeue
# Re-integration step 2: Seurat RPCA on the new 9-dataset manifest (895,228 cells).
# Resource profile copied from 01.1_Unify_gene_names.sh, the original pipeline's heaviest step.
# Deletes nothing except the TMPDIR this job itself creates.
set -eu -o pipefail
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg R_LIBS_USER=
export OMP_NUM_THREADS=58 OPENBLAS_NUM_THREADS=58 MKL_NUM_THREADS=58 \
       NUMEXPR_NUM_THREADS=58 TBB_NUM_THREADS=58 VECLIB_MAXIMUM_THREADS=58 BLIS_NUM_THREADS=58
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
echo "=== ENV ==="; hostname; df -h "$TMPDIR" | tail -1
cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/ws/05_rpca_integrate.R "$@"
