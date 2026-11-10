#!/bin/bash
#SBATCH -J ws_a1_g09
#SBATCH -p gr10634c
#SBATCH -t 12:00:00
#SBATCH --rsc p=1:t=1:c=16:m=64G
#SBATCH -o /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.out
#SBATCH -e /FAST/gr10634/gaozy/aml_niche_net/00_project/logs/%x_%j.err
#SBATCH --no-requeue
# WP-A1 G09 stage 1 (PREREGISTRATION_A1_G09.md, tag prereg-a1-g09-v1): one mode of scripts/ws/45 per job.
#   sbatch scripts/ws/submit_ws_a1_g09.sbatch --known | --run | --combine [--smoke]
# Each job runs a dated read-only copy of the script, so an edit on disk cannot land in a running job; the copy's md5 is in every output.
set -eu -o pipefail
ulimit -c 0   # no core files in the repo root (the 23:30 cron runs git add -A)
export R_FUTURE_GLOBALS_MAXSIZE=68719476736 PYTHONNOUSERSITE=1 R_LIBS_USER=
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1   # parallel work is forked by mclapply instead
export OMP_PROC_BIND=false   # the job environment sets OMP_PROC_BIND=spread, which pins R (and every forked worker) to one core
export LC_ALL=en_US.UTF-8    # K0b: S8's locale; table() orders S8's pasted bin labels by collation, which sets N0's head/tail split
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-16}
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
JOBTAG=${SLURM_JOB_ID:-local_$$}   # two local runs never share a tmp dir
export TMPDIR=/LARGE1/gr10634/gaozy/tmp/${JOBTAG}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd /FAST/gr10634/gaozy/aml_niche_net
CODE=/LARGE1/gr10634/gaozy/aml_niche_net/ws_a1_g09/code
SRC=scripts/ws/45_a1_g09.R
case " $* " in
  *" --smoke "*) CODE="$CODE/smoke" ;;   # smoke copies sit apart, so the 8.5 clean-up finds them
  *) git ls-files --error-unmatch "$SRC" >/dev/null 2>&1 && git diff --quiet HEAD -- "$SRC" || { echo "commit $SRC before a real run"; exit 2; }
     echo "real run of $SRC at commit $(git rev-parse HEAD)" ;;
esac
mkdir -p "$CODE"
c="$CODE/$(basename "$SRC" .R)_$(date +%Y%m%d_%H%M%S)_${JOBTAG}.R"; cp "$SRC" "$c"; chmod a-w "$c"
echo "running $c (md5 $(md5sum "$c" | cut -d' ' -f1)) $*; $(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc) CPUs allowed to this shell"
Rscript "$c" "$@"
