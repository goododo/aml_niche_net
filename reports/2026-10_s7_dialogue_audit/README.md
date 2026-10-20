# S7 DIALOGUE audit (2026-10-19)

Audit of the first S7 evaluation (job 3650058, 2026-10-17) and of the corrected rebuild (jobs 3650849, 3650854).
Scripts in `code/` are the exact diagnostics run on the saved objects under `/LARGE1/.../ws_dialogue/`; they read
only. Pre-registration record: `scripts/ws/PREREGISTRATION_S7_dialogue.md` section 8.

## Findings

1. **RNG reset** (`insp_draws.R`). `DIALOGUE1.PMD` calls `set.seed(seed1)`. In the first evaluation, stability draws
   2-20 and all 20 shuffles were one draw each (identical sample sets and p matrices). Fixed in 30 by drawing every
   subset and permutation before the first DIALOGUE call (distinct by construction).
2. **emp.p semantics** (`insp_p.R`). DIALOGUE's empirical p is a one-sided Wilcoxon rank-sum p of the observed
   canonical correlation against 99 permutations: rank 1 of 100 -> 0.0448, rank 2 -> 0.0482, rank 3 -> 0.0517,
   rank 4 -> 0.0555. The lookup reproduces every stored value exactly; the exact permutation p (rank/100) is now
   reported and, by the user's decision of 2026-10-19, used for criterion 1.
3. **Compartment assignment bug in DIALOGUE 1.0** (`ict_bug.R`). `DIALOGUE.identify.cell.types` indexes
   `R$cell.types` (T, NK, B, Mono) with a logical vector computed in `generic.vector2mat`'s alphabetical row order
   (B, Mono, NK, T). R1/R2 MCP1-2 are T+NK+B, not T+B; the planted run's MCP2-3 were T+NK, not B+Mono. Patched in
   `config_ws.R` section 13 (`ws_patch_dialogue()`), applied before every `DIALOGUE.run`.
4. **Planted control, first version** (`vis_null.R`). The factor had been added to the expression matrix after the
   PCA, so step I never saw it (planted emp.p identical to R1). Rebuilt with planting before the PCA: the factor is
   visible along the planted direction in PC space (cor 0.96-0.98 vs a chance P95 of about 0.5).
5. **DIALOGUE step I does not find the planted program at n = 13** (`planted_diag.R`, `stepI_diag.R`,
   `stepI_diag2.R`, `output/stepI_diag2.txt`). In DIALOGUE's own processed sample space (13 common samples, sample
   medians of 28-30 per-PC-standardised PCs), the planted direction correlates 0.95-0.99 with the factor in every
   compartment and 0.96-1.00 across compartments (mean pairwise 0.974). The three canonical variates DIALOGUE
   returned have |cor| <= 0.36 with the factor and mean pairwise canonical correlations of 0.80-0.87, which is the
   level its own permutation null reaches on noise (30 permuted data sets: median 0.86, P95 0.91). Sparsity is not
   the limit (9-14 non-zero weights per compartment); the sparse CCA converges to noise-level local optima at this
   sample size and feature count. Step II then has nothing to find: the planted MCP had no two-compartment gene
   signature, and planted genes were not enriched among the step-I gene correlations.

## Implication

The pre-registered positive control fails before the gene level. A null result from S7 in GSE185381 adults
(13 patients with all four compartments) is therefore uninformative about H-A1b: DIALOGUE, as configured
(30 PCs, k = 3, abn.c = 15), cannot recover a 100-gene, 0.5-SD shared program here. On the real data no MCP passes
criterion 1 under either reading once compartments are assigned correctly (the weak third pair enters the maximum),
so the corrected evaluation (job 3650855) will report "none found" together with the failed known answer.
Whether any configuration has power here (fewer PCs; the 23-sample R2 input) is a tool-calibration question to be
answered on planted data before any further real-data run.

## Tool checks on planted data (2026-10-20; prereg section 8, item 7)

| configuration | step I links | largest exact p | planted genes in top 200 (T/NK/B/Mono) | section 4.1 rule |
|---|---|---|---|---|
| 13 patients, 10 PCs | all four | 0.05 | 0.94 / 0.92 / 0.88 / 0.88 | PASS |
| 23 samples, 30 PCs | T+NK+B | 0.01 | 0.34 / 0.84 / 0 / 0 | FAIL |
| 23 samples, 10 PCs | all four | 0.01 | 1.00 / 0.99 / 0.99 / 0.96 | PASS |

Reading: the failure in finding 5 is driven by the 30 per-PC-standardised features per compartment, not by n alone.
With 10 PCs the sparse CCA recovers a 100-gene, 0.5-SD shared program at n = 13 and n = 23. Specificity under
10 PCs (shuffled-sample false-positive rate) has not been measured yet and is the gate of the v2 protocol.
Evaluator: `scripts/ws/31_s7_planted_toolcheck.R`; builds: `scripts/ws/29_s7_dialogue_build_run.R --npcs --tag`.

## v2 specificity gate (2026-10-20; prereg section 9)

Built R1 and R2 cell.type objects with 10 PCs (no DIALOGUE on the real pairing), then 20 sample relabellings within
each compartment per run. R1: 10 of 20 shuffles gave an MCP with largest-pair exact p <= 0.05 (5 of 20 with DIALOGUE's
p), identical to the 30-PC rate; R2: 6 of the first 11 (final count in `results/tables/ws/ws_s7_shuffle_R2_v2.csv`).
Gate ceiling 1 of 20 -> FAIL -> S7 closed as untestable, no real-data 10-PC run. Reading: DIALOGUE's empirical p comes
from permuting each PC independently across samples, which is not the null for a sample-level pairing question; under
sample relabelling its step I finds "significant" cross-compartment structure in about half of the data sets at
n = 13-23. Sensitivity (planted recovery) was fixed by the PC count; specificity was not.
