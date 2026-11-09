#!/bin/bash
#SBATCH -J dbl_smoke
#SBATCH -p gr10634c
#SBATCH -t 12:00:00
#SBATCH --rsc p=1:t=1:c=8:m=96G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
#SBATCH --no-requeue
# Doublet-pairing SMOKE (PROPOSAL_doublet_pairing_external_v1.3; decisions of 2026-11-09): one stage of scripts/dbl/ per job.
#   sbatch scripts/dbl/submit_dbl.sbatch 01 [--lanes=single,pooled] [--force]
#   sbatch scripts/dbl/submit_dbl.sbatch 02 | 03 [--scdbl=devel] | 04 | 05     (stage-specific flags are passed through)
# Each job runs a dated read-only copy of the stage script, so an edit on disk cannot land in a running job; md5 in the log.
set -eu -o pipefail
export R_FUTURE_GLOBALS_MAXSIZE=68719476736 PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_PROC_BIND=false   # the job environment sets OMP_PROC_BIND=spread, which pins R (and forked workers) to one core
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-8}
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
CODE=/LARGE1/gr10634/gaozy/aml_niche_net/dbl/smoke/code; mkdir -p "$CODE"
stage="${1:?give a stage: 01 02 03 04 05}"; shift
f=$(ls scripts/dbl/${stage}_*.R | head -1)
[ -f "$f" ] || { echo "no script for stage $stage" >&2; exit 2; }
c="$CODE/$(basename "$f" .R)_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID}.R"; cp "$f" "$c"; chmod a-w "$c"
echo "running $c (md5 $(md5sum "$c" | cut -d' ' -f1)) $*; $(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc) CPUs allowed"
# the copy lives outside the repo; here::here() still finds the project root because the working directory is the repo root
Rscript "$c" "$@"
