#!/bin/bash
# 19_verify_all.sh ----
# Run the whole pipeline audit in dependency order and stop at the first batch that fails.
#
# The order is not cosmetic: batch 0 validates the cells every later stage reads, a wrong bin map
# (batch 1) makes batch 5's disagreement statistics meaningless, and a scrambled reference (batch 2)
# makes every malignancy number downstream meaningless. Advancing past a failure produces numbers
# that look fine and mean nothing, which is the exact failure mode this suite exists to prevent.
#
# Each batch is independently runnable; this is just the gate.
#   bash scripts/99_admin/19_verify_all.sh

set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
RS=/FAST/gr10634/gaozy/general_env/bin/Rscript

BATCHES=(
  "0|per-sample QC objects (thresholds recomputed from counts)|scripts/99_admin/18_verify_qc_objects.R"
  "0b|BoneMarrowMap projection (now the sole source of every bin)|scripts/99_admin/26_verify_bmm_projection.R"
  "1|config layer (bin vocabulary, mapping lint)|scripts/99_admin/20_verify_annotation_config.R"
  "2|inferCNV reference construction (leave-one-out, block sizes)|scripts/99_admin/21_verify_infercnv_reference.R"
  "3|per-cell calls and consensus (arms vs evidence on disk)|scripts/99_admin/22_verify_malignancy_percell.R"
  "4|marker scoring (independent re-implementation)|scripts/99_admin/23_verify_marker_scoring.R"
  "5|annotation reconciliation (correction rule)|scripts/99_admin/24_verify_reconciliation.R"
  "6|annotated objects (id join, QC objects untouched)|scripts/99_admin/25_verify_annotated_objects.R"
)

pass=0; failed=""
for b in "${BATCHES[@]}"; do
  IFS='|' read -r num desc script <<< "$b"
  printf '\n############ BATCH %s: %s ############\n' "$num" "$desc"
  if [[ ! -f "$script" ]]; then echo "  [SKIP] $script not found"; continue; fi
  # PIPESTATUS, not a second run: batch 4 re-reads Seurat objects and running it twice to recover
  # an exit code the pipe swallowed would double the cost of the slowest check in the suite.
  "$RS" "$script" 2>&1 | grep -E "^\s*\[(PASS|FAIL)\]|^batch .* PASS|=========== BATCH"
  if [[ ${PIPESTATUS[0]} -eq 0 ]]; then
    pass=$((pass+1))
  else
    failed="$num"
    printf '\n>>> BATCH %s FAILED -- stopping. Later batches read what this one validates.\n' "$num"
    break
  fi
done

printf '\n############ SUITE: %d batch(es) passed ############\n' "$pass"
if [[ -n "$failed" ]]; then printf 'first failing batch: %s\n' "$failed"; exit 1; fi
echo "all batches pass"
