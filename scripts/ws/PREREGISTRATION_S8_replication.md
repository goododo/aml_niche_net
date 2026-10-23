# PRE-REGISTRATION — S8: structural replication of the T-cell imprint in independent adult AML cohorts

Written 2026-10-23; reviewed and approved by the user the same day. **BINDING** with tag `prereg-s8-v1`; later changes
are appended, dated sections. Nothing below has been run on the replication sets. Context: handoff v2.0 §16.1 item 15 (narrowed claim) and §16.2 item 13 (S7 closed).

## 1. Question

The T-cell imprint of S5 (prereg-s5-v1 section 11) separates adult AML at diagnosis from young healthy donors in
GSE185381. Is it the same axis of T-cell variation in adult AML cohorts that were never used to define it? This is
a *structural* replication: the cohorts have no healthy donors, so neither the AML-vs-healthy direction nor the age
confound can be tested here. Three structural properties are tested, each against matched random gene sets:
(T2) the up- and down-halves of the imprint move oppositely across patients; (T3) a patient's score is the same
whether it is read in CD4 or in CD8 T cells; (T1) within a patient's T cells the imprint genes co-vary more than
matched random genes.

## 2. Inputs (fixed)

- Gene sets. T imprint = `results/tables/ws/ws_s5_imprint_genes.csv`, compartment T, imprint_called TRUE:
  402 genes (171 up / 231 down). The 32 genes in the extended ambient ranking (GATE_A1 (2); list in
  `results/tables/ws/ws_s8_t_imprint_ambient_prone.csv`, 15 up / 17 down) are removed for the primary analysis:
  **primary set 370 genes (156 up / 214 down)**; the full 402 is a sensitivity analysis.
- Replication set (decisive): GSE289435, adult AML_Dx samples with >= 20 concordant T cells under the S4 rules:
  8 patients (54-543 T cells). Second set, same tests, reported alongside but not decisive: GSE227903 diagnosis
  samples, 7 patients (36-386 T cells); its MRD and relapse samples are not touched (WP-A2). Descriptive only:
  GSE116256 (Seq-Well, different platform; 6 adults with >= 20 T cells, 2 healthy donors with >= 20).
- Cells: S4 kept cells (concordant, not high_error), T group = CD4_T + CD8_T (`scripts/ws/s4_cell_groups.tsv`);
  pseudobulk per sample and group by `26_s4_pseudobulk.R`, unchanged, run for GSE227903 and GSE116256
  (GSE289435's object exists). Dataset roles stay as frozen in S4: GSE289435 replicates, the others describe.
- Scoring (as S5): per-sample T pseudobulk log CPM (edgeR `cpm`, prior count 1; genes with CPM > 1 in at least
  half of the set's samples), z-scored per gene across the set's samples; score = mean z(up) - mean z(down).
- Null (as S5): 1,000 random gene sets with the same number of up and down genes, drawn from the tested genes
  within mean-log-CPM deciles matched to the imprint, run through the identical computation.
- Covariates: `hb_share_nonery` and `gran_share_lymph` from the S4 meta tables, reported as partial correlations
  with the score (informational, as GATE_A1 (1)); no adjustment of the tests themselves.

## 3. Tests and the gate

All on the primary set in GSE289435; the 402-gene set and the other cohorts are reported with the same statistics.

- T2, opposite movement: Spearman rho across patients between mean z(up) and mean z(down). Pass if rho < 0 and
  below the random-set P5.
- T3, subtype concordance: scores recomputed on CD4_T and on CD8_T pseudobulks separately (samples with >= 20 cells
  in each); Spearman rho across patients between the two. Pass if rho > 0.5 and above the random-set P95. Requires
  >= 6 patients with both subtypes; otherwise "not testable".
- T1, within-patient coherence: per patient, T cells capped at 500, genes detected in >= 10% of that patient's
  T cells (each half must keep >= 20 genes, else the patient is skipped); coherence = mean Pearson r among the
  up-genes' log-normalised expression minus the mean r between up- and down-genes; compared with 200 random matched
  sets per patient. Pass if the imprint exceeds the random P95 in >= 6 of 8 patients.
- **GATE_S8 (structural replication): T2 AND T3 pass in GSE289435**; T1 is supporting. If T3 is not testable, the
  gate is T2 AND T1.

## 4. Known answers, run before the replication result is read

1. Positive: the same three tests on the discovery cohort (GSE185381 adult AML, 20 patients) must pass; this
   checks the machinery, not the hypothesis (the imprint was defined there).
2. Negative: 100 sign-scrambled versions of the imprint (the same 370 genes with up/down labels permuted) must pass
   T2 in at most 5 of 100 runs on GSE289435.
3. The random-set null is calibrated by construction; its P5/P95 are reported.

## 5. Reading the result

- Gate passes: the narrowed claim gains its replication half as *structure*: "a reproducible T-cell state axis in
  adult AML marrow at diagnosis". Direction (AML vs healthy) and age still require older healthy donors
  (handoff §16.2 item 6; separate decision), and the 32 ambient-prone genes stay excluded from any claim.
- Gate fails: WP-A has no replicated positive result. Per the stop rule agreed on 2026-10-23, no further method is
  tried on this question; the thesis line is re-planned around what the data support (descriptive atlas plus the
  methodological audits), recorded as a new §16.1 decision.
- No threshold is changed after seeing a replication number; later changes are appended, dated sections.

## 6. Outputs and code

`results/tables/ws/ws_s8_scores_<set>.csv`, `ws_s8_tests.csv` (set x gene-set x test: statistic, random P5/P95,
pass), `ws_s8_known_answers.csv`, `ws_s8_checks.csv`; GATES.md block "A1-S8". Code after approval: `26_s4_pseudobulk.R`
for GSE227903 and GSE116256 (no change to the script); new `scripts/ws/32_s8_replication.R` (+ sbatch), about 130
lines, no new dependencies.

## 7. Deviations found at the first run (appended 2026-10-23, before the corrected run)

- The gene filter in `32_s8_replication.R` added edgeR's `filterByExpr` to the pre-registered rule (CPM > 1 in at least
  half of the set's samples). With 8 ungrouped samples that function requires expression in every sample, and only
  32% of the primary genes were in GSE289435's tested universe; the coverage check (section 6 of the script) caught
  it. The first full run (job 3658344) was cancelled minutes after starting and the filter was corrected to the
  pre-registered rule. Seen before the correction, from the smoke run on GSE289435 (100 random sets): T2 -0.88
  (random P5 -0.74, pass), T1 8 of 8 patients (pass), T3 0.94 with a random-set P95 of 0.94 (fail), and the scrambled
  imprint passing T2 in 0 of 10 scrambles. No threshold is changed.
- T3 as pre-registered is degenerate at n = 6: random matched gene sets scored in CD4 and CD8 T cells of the same six
  patients reach a Spearman concordance of 0.94, the largest value below 1 that six ranks allow, so no gene set can
  exceed the random P95. The corrected run reports the share of random sets at or above the observed concordance.
  How T3 is read is a decision for the user, recorded as a further appended section before any new statistic is
  computed; T2 and T1 are unaffected.

## 8. Outcome (appended 2026-10-23; job 3658345, corrected filter; tables results/tables/ws/ws_s8_*.csv; GATES.md block A1-S8)

| set | n | T2 rho (random P5) | T3 rho (random P95; share of random sets >= observed) | T1 patients above P95 |
|---|---|---|---|---|
| disc GSE185381 (positive control) | 15 | -0.29 (0.41) pass | 0.78 (0.84; 0.12) FAIL | 13 of 15 pass |
| rep GSE289435 (decisive) | 8 | -0.55 (0.38) pass | 1.00 (0.94; 0.03) pass, 6 patients | 8 of 8 pass |
| rep2 GSE227903 diagnosis | 7 | -0.43 (0.21) pass | not testable (4 patients) | 7 of 7 pass |
| desc GSE116256 (Seq-Well) | 6 | +0.60 (0.54) fail | not testable (2) | 3 of 5 fail |

Full 402-gene set: same pattern in every set (rep T2 -0.57, T3 1.00, T1 8/8; disc T3 0.74 fail; desc T2 +0.60 fail).
Sign-scrambled imprint: 0 of 100 scrambles pass T2 on GSE289435. Score vs ambient shares (primary): disc hb 0.07 /
gran -0.35, rep 0.14 / -0.24, rep2 0.43 / 0.04, desc -0.68 / 0.23 (informational).

- **GATE_S8 as pre-registered (T2 and T3 on GSE289435): PASS.** Under the fallback reading (T3 not testable -> T2 and
  T1): also PASS. The structural replication therefore does not depend on how T3 is read.
- Known answer 1 (positive control) FAILS through T3 only: in the discovery cohort 12% of matched random gene sets
  are at least as concordant between CD4 and CD8 T cells as the imprint. T3's null is dominated by patient-level
  factors shared by both subtypes of the same patient (its median is 0.59-0.77), so T3 cannot separate a real axis
  from any large gene set; its pass in GSE289435 sits at the ceiling (rho 1.00; 3% of random sets also reach it).
  T3 is therefore recorded as uninformative in this design; the evidence for replication rests on T2 and T1, both
  of which pass in the discovery cohort and in both 10x replication cohorts and fail for scrambled imprints.
  How T3 is formally recorded (fail of a known answer vs declared uninformative) is the user's decision (§9 to append).
- GSE116256 (Seq-Well, 6 patients, 21 T cells median) does not show the structure (T2 positive, T1 3 of 5). It was
  descriptive by design; the structural replication is shown on 10x Chromium data only, and is reported so.
- Open item: 26's "HB share equals the S1 table" check fails for GSE227903 (passes for the other three datasets);
  HB share is only an informational covariate here; to be looked at before any use of GSE227903 beyond description.
