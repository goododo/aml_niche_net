#!/bin/bash
# e10_array.sbatch ----
# SLURM array driver for per-sample cNMF (e10). One array task per sample.
# Build the sample list first (dataset<TAB>sample), then submit with the right array size.
#
# Build list (run once, from scripts/06_hierarchy/):
#   ls -d /LARGE1/gr10634/gaozy/aml_niche_net/06_cnmf/input/*/*/ \
#     | awk -F/ '{print $(NF-2)"\t"$(NF-1)}' > cnmf_sample_list.tsv
#   wc -l cnmf_sample_list.tsv          # -> set --array=1-N below
#
# Submit:
#   sbatch --array=1-41 e10_array.sbatch
#
#SBATCH -p gr10634c
#SBATCH --job-name=cnmf
#SBATCH --chdir=/FAST/gr10634/gaozy/aml_niche_net/scripts/06_hierarchy
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/scripts/06_hierarchy/logs/cnmf_%A_%a.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/scripts/06_hierarchy/logs/cnmf_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --rsc p=1:t=8:c=8:m=32G

set -euo pipefail
ENV=/FAST/gr10634/gaozy/general_env
SCRIPTS=/FAST/gr10634/gaozy/aml_niche_net/scripts/06_hierarchy
LIST=$SCRIPTS/cnmf_sample_list.tsv

cd "$SCRIPTS"
mkdir -p logs

# pick the line for this array index
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
ds=$(echo "$line" | cut -f1)
sm=$(echo "$line" | cut -f2)
echo "[array ${SLURM_ARRAY_TASK_ID}] dataset=$ds sample=$sm  $(date)"

# limit BLAS/OMP threads so cNMF's own parallelism doesn't oversubscribe the 8 cpus
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

conda run -p "$ENV" python e10_run_cnmf.py --dataset="$ds" --sample="$sm"
echo "[array ${SLURM_ARRAY_TASK_ID}] done  $(date)"
