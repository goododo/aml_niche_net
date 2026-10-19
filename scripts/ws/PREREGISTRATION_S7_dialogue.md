# PRE-REGISTRATION — S7: unsupervised multicellular programs (DIALOGUE), the last test of H-A1b

Written 2026-10-17, after S5 (prereg-s5-v1, section 11) and before any DIALOGUE run on project data. Decision
that motivates it (handoff v2.0 §16.1 item 15): the supervised imprint does not co-vary across compartments
(GATE_A2 failed). If an unsupervised method that looks for exactly such co-variation also finds none, H-A1b
("the imprint is coordinated across compartments within a patient") is closed.

Status: **BINDING** once committed with tag `prereg-s7-v1`; later changes are appended sections.

## 1. Question

Is there any multicellular program (MCP) in GSE185381 adults, i.e. a set of genes in each of several
compartments whose sample-level expression co-varies across patients, that survives the same confounders the
imprint had to survive and is stable under patient resampling?

## 2. Method and inputs (fixed)

- Tool: DIALOGUE (Jerby-Arnon & Regev 2022), R package from github.com/livnatje/DIALOGUE, version 1.0 as
  installed 2026-10-17; `DIALOGUE.run` with `DLG.get.param`. No re-implementation.
- Cells: the S4 kept cells (concordant, not high_error, four compartments T, NK, B, Mono; config_ws.R
  section 12). Per compartment: Seurat LogNormalize, 2,000 variable genes, 30 PCs computed within the
  compartment on the run's cells. `cellQ` = scaled log total UMI.
- Samples: a compartment contributes a sample when it has >= 15 cells there (DIALOGUE `abn.c = 15`); DIALOGUE
  uses the samples present in every compartment of the run. Cells per sample and compartment are capped at 500
  (random, seed 491638) so that the gene-level mixed models stay tractable; the cap is applied before the PCA.
- Two runs:
  - R1 adult AML only (expected n about 11 patients present in all four compartments);
  - R2 adult AML + healthy donors (expected n about 20), with `status` as a covariate.
- Confounders handed to DIALOGUE: `conf = "cellQ"`; covariates in its gene-level mixed models:
  `covar = c("cellQ", "hb_share_nonery", "gran_share_lymph", "prog_frac")` plus `status` in R2.
- k = 3 MCPs per run; `n.genes = 200`; all other parameters at their defaults; `seed1 = 491638`.

## 3. What counts as "an MCP was found"

All three, for at least one MCP in at least one run:
1. DIALOGUE's own empirical p for the MCP <= 0.05 (its permutation of the sample pairing across compartments).
2. After the MCP's sample-level scores per compartment are residualised on prog_frac, the two composition PCs,
   hb_share_nonery and gran_share_lymph (the GATE_A2 covariates), the mean pairwise Spearman across the
   compartments that carry the MCP is > 0.30, with a patient-bootstrap 95% CI excluding 0.
3. Stability: in 20 random 70% patient subsets, DIALOGUE is rerun with the same parameters and the MCP is
   matched to the full-run MCP by its compartment score correlation; the median Jaccard of the top-50 genes
   per compartment is >= 0.5 in at least two compartments.

If no MCP meets all three in either run, H-A1b is closed and recorded in GATES.md. If one does, it is reported
as an unsupervised candidate and is NOT a result: it gets its own pre-registered confirmatory test
(replication study, ambient and clonal checks) before anything is claimed.

## 4. Known answers, run before the real result is read

1. Planted MCP: in a copy of R1's inputs, a random sample-level factor is added to 100 random genes in each
   compartment (same sign per sample across compartments, effect 0.5 SD of the gene's log expression);
   DIALOGUE must recover it (empirical p <= 0.05; >= 50% of the planted genes among the MCP's top 200 in each
   compartment).
2. Shuffled samples: sample labels are permuted independently within each compartment (destroying any
   cross-compartment pairing); no MCP may pass criterion 1 at the nominal rate above 5% over 20 shuffles.

## 5. Not done

No parameter is changed after seeing the MCPs; no compartment is dropped after seeing the MCPs (Mono stays in,
and its clonal admixture is a reason for criterion 3 and for the confirmatory test, not for exclusion now); no
gene list from S5 enters this analysis.

## 6. Outputs

`results/tables/ws/ws_s7_mcps.csv` (MCP x compartment: genes, empirical p), `ws_s7_scores.csv` (sample x
compartment x MCP), `ws_s7_criteria.csv` (criteria 1-3 per MCP and run), `ws_s7_known_answers.csv`,
`ws_s7_checks.csv`; DIALOGUE's own result objects under `/LARGE1/.../ws_dialogue/`; GATES.md block "A1b-S7".

## 7. Build details fixed before the first successful run (appended 2026-10-17)

- Genes handed to DIALOGUE per compartment: those detected in >= 1% of that compartment's cells (DIALOGUE
  holds a dense expression matrix; the full 36,601-gene matrix exceeded the login node's memory).
- Dependencies installed for DIALOGUE 1.0: nnls, plotrix (CRAN), in addition to its declared imports.
- The first two attempts failed before producing any result (missing package; an error under investigation),
  so nothing has been seen.

## 8. Audit of the first evaluation (job 3650058, 2026-10-17) and corrections (appended 2026-10-19)

What had been seen before anything below was changed: the criterion-1 and criterion-2 values of every MCP in R1
and R2 (DIALOGUE p 0.052-0.912; the two T+B programs at 0.056 (R1 MCP1) and 0.052 (R2 MCP1)), their compartments,
and the fact that DIALOGUE returned no two-compartment gene signature for any MCP except R2 MCP3. The corrections
are therefore not blind. They change no parameter of section 2 and no threshold of sections 3 and 4.

Defects found (each verified on the saved objects under /LARGE1/.../ws_dialogue):

1. The planted control (section 4.1) never reached DIALOGUE step I. 29 added the factor to the expression matrix
   after the PCA; step I works on the sample-averaged PCs, so the planted run's empirical p matrix was identical
   to R1's and the control tested nothing. Correction: the factor is added to the log-expression of the 100 genes
   before ScaleData/RunPCA, the planted genes are added to the variable features, and a build check requires the
   factor to be visible in the sample-mean PCs of every compartment (max |cor| >= 0.8). Only the planted run is
   rebuilt; the R1 and R2 objects of 2026-10-17 are reused unchanged. The invalid outputs are kept as
   planted_v1_tpm_only/, R_planted_v1_tpm_only.rds, rA_planted_v1_tpm_only.rds, planted_truth_v1_tpm_only.rds.
2. DIALOGUE1.PMD calls set.seed(seed1). Subsets and permutations drawn inside the loops after a DIALOGUE call
   therefore repeated: stability draws 2-20 and all 20 shuffles were one draw each (identical sample sets and
   p matrices on disk). Correction: the 20 subsets per run and the 20 permutations are drawn, distinct, before
   the first DIALOGUE call, and a check records it. The invalid outputs are kept as stability_v1_rngbug/ and
   shuffle_v1_rngbug/.
3. Criterion 1's statistic. DIALOGUE's emp.p is a one-sided Wilcoxon rank-sum p (normal approximation) of the
   observed canonical correlation against 99 permutations (n1 = 100, fixed inside DIALOGUE1), i.e. a deterministic
   function of the observed value's rank: rank 1 of 100 gives 0.045, rank 2 0.048, rank 3 0.052, rank 4 0.056,
   rank 6 0.064, rank 11 0.088. The pre-registered "emp.p <= 0.05" therefore meant "top 2 of 100" (exact p <= 0.02),
   stricter than intended. Decision (user, 2026-10-19, after seeing the values above): criterion 1 uses the exact
   permutation p = rank / 100, recovered by lookup from DIALOGUE's p (the lookup reproduces every stored emp.p
   entry to < 1e-8); DIALOGUE's p and the verdict under the original reading are reported alongside (columns
   crit1_prereg and found_prereg). Effect on the first run's values: R1 MCP1 (exact p 0.04) and R2 MCP1 (0.03)
   pass criterion 1 under the new reading; both lack a two-compartment gene signature, so criterion 3 is undefined
   for them and the verdict cannot change through this decision alone. The planted and shuffled known answers
   use the same rule.
4. Minor: R2's sample set for the resampling is DIALOGUE's own (23 samples, from its full output) instead of the
   intersection of the score tables (24); 29's check table is renamed ws_s7_build_checks.csv so that it no longer
   overwrites 30's ws_s7_checks.csv; 30 prints lme4/DIALOGUE warnings in-line with the run label (9 rank-deficiency
   warnings in the main runs, 3 per stability draw, 12 per shuffle; not yet diagnosed); an --out-dir option
   routes a smoke test's outputs away from results/tables/ws and GATES.md.

GATES.md: the 2026-10-17 block is marked withdrawn by a note below it; the rerun regenerates the block.
Not changed: k, abn.c, covariates, the 500-cell cap, feature and PC counts, seeds, 20 draws at 70%, and the
thresholds of criteria 2 and 3.
