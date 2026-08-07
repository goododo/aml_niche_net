#!/usr/bin/env bash
# run_consensus_all.sh ----
# PRODUCTION full-cohort consensus (union default). Iterates every inferCNV burden CSV across all
# datasets, converts to percell (41) if needed, then runs the consensus (50). Numbat is added when
# a percell file exists (mostly GSE227903); everything else is inferCNV-only -> tier C_single.
#
# RESUME-SAFE, WITH A FRESHNESS TEST. An output is reusable only if it POSTDATES the input it was
# computed from -- same invariant as the QC report combiner and 44. Testing existence alone made
# this script unrunnable after any upstream change, in the direction that matters: 130 consensus
# summaries dated 2026-07-14 sit on disk, 68 of them for samples whose burden CSV is being
# recomputed right now, so a plain -f guard would skip every one and leave ALL_consensus_summary
# carrying v1 calls while reporting "already-done". The percell step (41) had the same guard.
# Set CONSENSUS_FORCE=1 to rebuild regardless.
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

FORCE="${CONSENSUS_FORCE:-}"
n_run=0; n_done=0; n_skip=0; n_stale=0
for f in "${BURDEN}"/*/*_infercnv_burden.csv; do
  [[ -e "$f" ]] || continue
  ds=$(basename "$(dirname "$f")")
  s=$(basename "$f" _infercnv_burden.csv)
  out="${MAL_TAB}/${ds}"
  summ="${out}/${s}__consensus_summary.csv"
  qc="${QC_RDS_DIR}/${ds}/${s}.rds"
  ic="${INFERCNV_ROOT}/${ds}/${s}/${s}__infercnv_percell.csv"
  nb="${CNV_ROOT}/numbat/${ds}/${s}/numbat/${s}__numbat_percell.csv"

  # The roster check comes FIRST. A sample that left the cohort has no QC object
  # (they are quarantined by 99_admin/14), and it must not be reported as
  # "already-done" just because a burden CSV from the previous ingest survives.
  if [[ ! -f "$qc" ]]; then echo "[skip] not in cohort / no QC rds: ${ds}/${s}"; n_skip=$((n_skip+1)); continue; fi

  if [[ -z "$FORCE" && -f "$summ" && "$summ" -nt "$f" && "$summ" -nt "$qc" ]]; then
    n_done=$((n_done+1)); continue                                      # resume-skip: output is current
  fi
  [[ -f "$summ" ]] && { echo "[recompute] stale summary: ${ds}/${s}"; n_stale=$((n_stale+1)); }

  # 41: burden -> percell. Same freshness test: a percell file older than the
  # burden it derives from is not a reusable output, it is a stale one.
  if [[ -n "$FORCE" || ! -f "$ic" || ! "$ic" -nt "$f" ]]; then
    Rscript "${MAL_DIR}/41_infercnv_to_percell.R" \
      --sample "$s" --burden_csv "$f" --outdir "${INFERCNV_ROOT}/${ds}/${s}" >/dev/null
  fi

  args=(--sample "$s" --dataset "$ds" --qc_rds "$qc" --infercnv "$ic" --outdir "$out")
  [[ -f "$nb" ]] && args+=(--numbat "$nb")
  echo "[run] ${ds}/${s}$([[ -f $nb ]] && echo ' (+numbat)')"
  Rscript "${MAL_DIR}/50_consensus_malignancy.R" "${args[@]}" >/dev/null && n_run=$((n_run+1))
done

echo "[consensus] ran=${n_run} (of which ${n_stale} replaced a stale summary), current(skipped)=${n_done}, not-in-cohort(skipped)=${n_skip}"
echo "[consensus] outputs under ${MAL_TAB}/<dataset>/"
