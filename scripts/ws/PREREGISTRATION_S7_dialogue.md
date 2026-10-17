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
