#!/bin/bash
# THE PANEL SCREEN AT A PERMUTATION COUNT THAT CAN RESOLVE ITS OWN THRESHOLD.
#
# WHY 1e6 AND NOT 10000. The 2026-08-26 screen reported 7 of 114 features passing BH-FDR; a rerun
# differing only in the permutation stream gave 4; six realizations gave {4,6,7,7,8,8}, and of the 8
# features passing at least once only 2 passed every time. That is not instability in the biology --
# beta was bit-identical (0.0e+00) across all six. The `cs` family clears BH only through its k=3 step,
# with a margin of 0.59 Monte Carlo SE at n_perm=10000, so ~28% of streams lose all three at once;
# `pg` sits at 0.39 SE. At n_perm=1e6 the MC SE shrinks 10x and those margins become 5.9 and 3.9 SE.
# Fixing the shared RNG (done, tests T1-T3 pass) makes the answer REPRODUCIBLE; only this makes it
# RESOLVED. They are different problems and both had to be fixed.
#
# BOTH SCALINGS, because the config's choice of within_sample_rank was made by comparing the p-value
# of the hypothesis under test -- no valid provenance -- and because the two scalings' discovery hit
# sets at n_perm=10000 were disjoint (4 vs 1, zero overlap). Inputs are matched: identical edge tables
# (md5 d2d3138f), identical cohort, n_cells-truncation fixed in both.
set -uo pipefail
cd /FAST/gr10634/gaozy/aml_niche_net
PY=/FAST/gr10634/gaozy/general_env/bin/python
S=/LARGE1/gr10634/gaozy/tmp/screen1e6
NPERM=1000000

one () {
  local sc="$1"
  echo "[$sc] start $(date -Is)" >> "$S/PROGRESS.log"
  $PY scripts/08_scoring/07_feature_decomposition.py \
      --root "$S/$sc" --split discovery --split_csv "$S/$sc/01_preprocess/02_sample_split.csv" \
      --subsets panels --alphas 0 --n_perm $NPERM > "$S/$sc.log" 2>&1
  echo "[$sc] exit=$? $(date -Is)" >> "$S/PROGRESS.log"
}

one within_sample_rank &
P1=$!
one global_z &
P2=$!
wait $P1 $P2
echo "SCREEN COMPLETE $(date -Is)" >> "$S/PROGRESS.log"
