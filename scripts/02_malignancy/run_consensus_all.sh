#!/usr/bin/env bash
# run_consensus_all.sh ----
# PRODUCTION full-cohort consensus (union default). Iterates every inferCNV burden CSV across all
# datasets, converts to percell (41) if needed, then runs the consensus (50). Numbat is added when
# a percell file exists (mostly GSE227903); everything else is inferCNV-only -> tier C_single.
#
# RESUME-SAFE: a sample whose <sample>__consensus_summary.csv already exists is skipped, so you can
# re-run after a session cutoff and it continues where it stopped.
#
# Run (assumes the general_env is active; see the sbatch wrapper below):
#   bash scripts/02_malignancy/run_consensus_all.sh

set -uo pipefail

## ---- locate project root via the .here anchor ----
_find_here() {
  local d; d="$(cd "${1:-$PWD}" 2>/dev/null && pwd)" || return 1
  while [[ "$d" != "/" ]]; do
    [[ -e "$d/.here" ]] && { printf '%s\n' "$d"; return 0; }
    d="$(dirname "$d")"
  done; return 1; }
_self="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(_find_here "${_self}")" || true
[[ -n "${PROJECT_ROOT}" ]] || { echo "ERROR: no .here anchor found from ${_self}"; exit 1; }
source "${PROJECT_ROOT}/scripts/config/config_paths.sh"

MAL_DIR="${PROJECT_ROOT}/scripts/02_malignancy"
BURDEN="${LARGE1_DIR}/05_cnv_snv/infercnv_burden"
[[ -d "${BURDEN}" ]] || { echo "ERROR: no burden root: ${BURDEN}"; exit 1; }

n_run=0; n_done=0; n_skip=0
for f in "${BURDEN}"/*/*_infercnv_burden.csv; do
  [[ -e "$f" ]] || continue
  ds=$(basename "$(dirname "$f")")
  s=$(basename "$f" _infercnv_burden.csv)
  out="${MAL_TAB}/${ds}"
  summ="${out}/${s}__consensus_summary.csv"
  if [[ -f "$summ" ]]; then n_done=$((n_done+1)); continue; fi          # resume-skip

  qc="${QC_RDS_DIR}/${ds}/${s}.rds"
  ic="${INFERCNV_ROOT}/${ds}/${s}/${s}__infercnv_percell.csv"
  nb="${CNV_ROOT}/numbat/${ds}/${s}/numbat/${s}__numbat_percell.csv"
  if [[ ! -f "$qc" ]]; then echo "[skip] no QC rds: ${ds}/${s}"; n_skip=$((n_skip+1)); continue; fi

  # 41: burden -> percell (skip if already made)
  [[ -f "$ic" ]] || Rscript "${MAL_DIR}/41_infercnv_to_percell.R" \
      --sample "$s" --burden_csv "$f" --outdir "${INFERCNV_ROOT}/${ds}/${s}" >/dev/null

  args=(--sample "$s" --dataset "$ds" --qc_rds "$qc" --infercnv "$ic" --outdir "$out")
  [[ -f "$nb" ]] && args+=(--numbat "$nb")
  echo "[run] ${ds}/${s}$([[ -f $nb ]] && echo ' (+numbat)')"
  Rscript "${MAL_DIR}/50_consensus_malignancy.R" "${args[@]}" >/dev/null && n_run=$((n_run+1))
done

echo "[consensus] ran=${n_run}, already-done(skipped)=${n_done}, no-qc(skipped)=${n_skip}"
echo "[consensus] outputs under ${MAL_TAB}/<dataset>/"
