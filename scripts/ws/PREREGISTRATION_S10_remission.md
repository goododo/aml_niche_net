# PRE-REGISTRATION — S10: does the T-cell imprint of adult AML persist into remission (MRD time point)?

Written 2026-10-26. **DRAFT — not binding until the user approves it and the tag `prereg-s10-v1` is set.** No MRD or
post-treatment T pseudobulk has been scored or compared with its diagnosis sample before the tag. Later changes are appended,
dated sections. Context: handoff v2.0 §6 (WP-A2), §16.1 items 2, 8 and 15 (T imprint only; rank preservation + self-
recognition; written down before any MRD expression is read), appendix A (template), §17 S10 row. Upstream results this
document relies on: S5 (imprint definition, tag prereg-s5-v1), S8 (structural replication in GSE227903 diagnosis samples,
GATE_S8 PASS), S9 (age share bounded, GATE_S9 PASS). The relapse pairs (S11) are a separate pre-registration, not part of
this one.

## 1. Question and hypotheses

The 370-gene T imprint (S5 primary set: 402 called genes minus the 32 ambient-prone genes of S8) separates adult AML T cells
at diagnosis from healthy donors' T cells. S10 asks whether a patient's imprint level and imprint profile at diagnosis are
still recognisable in the same patient's T cells at the MRD time point after induction.

- **H-A2a rank preservation (test 1, gating):** across patients, the diagnosis imprint score and the MRD imprint score are
  positively rank-correlated (a patient with a stronger imprint at diagnosis has a stronger imprint at MRD).
- **H-A2b self-recognition (test 2, gating):** each MRD T profile over the imprint genes is more similar to its own diagnosis
  profile than to the other patients' diagnosis profiles, and this self-similarity exceeds what matched random gene sets
  give (random genes also carry individuality: genotype, sex, batch).
- **H-A2c (reference, reported):** MRD scores relative to GSE185381 healthy donors, in the joint space of S9.
- **H-A2d (relapse)** is S11, not here.

Both gating tests are one cohort, six adults, one batch. Section 8 states the power: adults give a direction; the formal
conclusion waits for the paediatric pairs (Mumme GSE235923, Lambo GSE235063; handoff §7).

## 2. Data and eligibility (frozen before the tag)

**2.1 Cohort.** GSE227903 (Ennis 2023, 10x 3'), S4 pseudobulk `ws_pb_GSE227903.rds` (26_s4_pseudobulk.R: kept = projected
T-group cells, concordant with the marker typing, not high_error). Diagnosis (`_Dg`) and MRD (`_MRD`) samples of the same
patient form a pair. Eligibility, read from `ws_s4_pseudobulk_meta_GSE227903.csv` (cell counts only; no expression):
both samples of the pair have >= 20 kept T cells (the S4/S9 minimum, not appendix A's 30), and both are in the same
library batch (`library_batch` of `meta_GSE227903_v2.3.csv`).

| patient | Dg kept T | MRD kept T | batch Dg / MRD | MRD status (rule B) | eligible |
|---|---|---|---|---|---|
| 1216 | 338 | 413 | nuig_batch2 / nuig_batch2 | unknown -> treated as remission | yes |
| 1886 | 57 | 273 | batch2 / batch2 | positive | yes |
| 3904 | 86 | 110 | batch2 / batch2 | unknown | yes |
| 40389 | 386 | 1385 | batch2 / batch2 | unknown | yes |
| 4978 | 176 | 278 | batch2 / batch2 | positive | yes |
| 5750 | 36 | 111 | batch2 / batch2 | unknown | yes |
| 6323 | 12 | 636 | nuig_batch1 / nuig_batch2 | unknown | **no** (12 kept T at diagnosis; cross-batch) |

**n = 6 pairs.** 5143 and 3853 have no MRD sample. These six are the registered set; no patient is added or dropped after a
number is seen. 1886 and 4978 are MRD-positive: they stay in the main test (the hypotheses are about the T cells, not the
blasts); the four MRD-negative/unknown patients are a reported sensitivity set, never a gate.

**2.2 Descriptive set.** GSE116256 (van Galen 2019, Seq-Well): patients with a diagnosis sample of >= 20 kept T cells and a
post-treatment sample with clinical blasts < 5 %: AML329 (D0 80 T; D20, 1 %, 88 T), AML420B (D0 149 T; D35, 1 %, 124 T;
D14 at 5 % is not < 5 %), AML556 (D0 24 T; D15, 0 %, 547 T; D31, 4 %, 64 T). AML707B-D0 has 8 kept T cells and is out;
AML328's post-treatment samples carry 15-30 % blasts and are out. Three patients, four post-treatment samples: scores and
per-patient differences are reported; no test, no wording (Seq-Well; S8 showed the imprint's structure does not replicate
on this platform).

**2.3 Gene sets.** Primary set = the 370 T-imprint genes (`ws_s5_imprint_genes.csv`, compartment T, minus
`ws_s8_t_imprint_ambient_prone.csv`), 156 up / 214 down. Test 2 additionally removes sex-chromosome and HLA genes (appendix
A item 6): the 10 chrX/chrY genes of the primary set (ARAF, IDS, HMGB3, GTPBP6, LINC00685, MSL3, LINC01278, ARHGEF9, CHM,
SEPTIN6; `gencode_GRCh38_gene_order.txt`) and HLA-A, HLA-B -> **358 genes (151 up / 207 down)**, before the expression
filter of 3.1. The full 402 set is a reported sensitivity for test 1 only. Blacklist for every pseudobulk: `black_for()` of
S8/S9 (HB genes, MPO/ELANE/AZU1/PRTN3, the 79-gene S1 ambient list, XIST, chrY, MT-, RPL/RPS).

## 3. Space and score

**3.1 Remission space R (GSE227903 only).** Raw T sums of the 12 samples (6 Dg + 6 MRD); blacklist removed; U_R = genes with
CPM > 1 in at least 3 of the 6 Dg samples AND in at least 3 of the 6 MRD samples; TMM on the 12 samples together (edgeR
calcNormFactors on U_R); log-CPM, prior count 1; per-gene mean m_g and SD s_g over the **6 diagnosis samples only**, applied
to all 12 (z_D for diagnosis, z_M for MRD). Diagnosis samples therefore have mean 0 / SD 1 per gene; MRD samples are read
in the diagnosis coordinates. Score = mean z(up) - mean z(down) over primary ∩ U_R (S5/S8/S9 formula). Coverage
|primary ∩ U_R| / 370 >= 0.80 is a known answer (K3).

**3.2 Joint space J (reference only, H-A2c).** As S9 section 2.4 with G = the 25 eligible GSE185381 T pseudobulks and E = the
12 GSE227903 samples: z with G's per-gene mean/SD; positions of the 6 Dg and 6 MRD samples relative to G's healthy mean,
AML mean and their midpoint are reported. No study offset can be estimated (GSE227903 has no healthy donors), so nothing
here enters a gate or a sentence.

## 4. Gating tests

**4.1 Test 1 — rank preservation.** Statistic: Spearman rho between the 6 diagnosis scores and the 6 MRD scores (primary
set, space R). Null: the 720 pairings of the 6 MRD samples to the 6 diagnosis samples (exact). One-sided (rho > 0).
**Pass: exact p < 0.05, which at n = 6 is rho >= 0.8286** (sum of squared rank differences <= 6; p = 12/720 = 0.0167; the
next value, rho = 0.7714, has p = 0.051 and fails). Reported beside it: Kendall tau, Pearson r of the scores, the 6 pairs
of scores, and the same rho on the full 402 set.

**4.2 Test 2 — self-recognition.** For each MRD sample i and each diagnosis sample j, Pearson correlation c_ij over the test-2
genes (358 ∩ U_R) between z_M,i and z_D,j. Rank r_i of j = i among the 6 (1 = own diagnosis sample is the most similar;
ties by average rank). Statistic M = mean of the 6 ranks (3.5 expected under no self-recognition; 1 is perfect).
Two null distributions, both must be passed:
(a) **pairing null:** M recomputed under the 720 pairings; exact one-sided p_a = share of pairings with M <= M_observed;
    pass p_a < 0.05;
(b) **random-gene null:** 1,000 random gene sets of the same size as 358 ∩ U_R, drawn from U_R minus chrX/chrY/HLA genes
    with S5's bin matching (deciles of mean log-CPM x detection bins over the 6 diagnosis samples; `matched()` of S8/S9,
    SEED), each giving M_rand with the TRUE pairing; pass M_observed < the 5th percentile of M_rand (the imprint identifies
    its owner better than 95 % of matched random sets).
Passing (a) but not (b) means "individuals resemble themselves" and is recorded as not passing (appendix A item 6). Spearman
instead of Pearson for c_ij is a reported sensitivity.

**4.3 Gate.** GATE_S10 = K0, K1, K3, K4 pass (section 6) AND test 1 passes AND test 2 passes (a and b). A test that fails
stops the line for discussion (P5); the wording table of section 7 is fixed now.

## 5. Reported analyses (never in claim wording)

1. **Per-patient change:** delta_i = MRD score - diagnosis score in space R; the six deltas, their mean and range; the share
   of patients whose MRD score stays above the diagnosis cohort's minimum. A common shift of all six does not touch test 1.
2. **MRD-negative/unknown subset (1216, 3904, 40389, 5750):** test-1 rho and the test-2 ranks recomputed; descriptive (at
   n = 4 the exact one-sided p of rho = 1 is 1/24; nothing is read as a pass or a fail).
3. **Composition:** per sample the kept-T shares of CD4 Naive, CD8 Naive and CD8 Effector Memory (1 + 2); the diagnosis ->
   MRD change; Spearman of the naive share between the two time points (does composition alone preserve rank?); a J-space
   age-share style covariate is not fitted (n = 6).
4. **Ambient and doublets (S1/S4 covariates):** per sample hb_share_nonery, gran_share_lymph, decontX median contamination
   (`results/tables/ws/decontx/GSE227903/<sample>__checks.csv`) and T3a doublet mean where it exists (NA for MRD samples).
   Scores are on raw counts; decontX-corrected counts are not used (handoff §16.1 item 1).
5. **H-A2c reference (3.2):** J-space positions of the 12 samples against GSE185381 healthy/AML means.
6. **Subtype pseudobulks (W-type sensitivity):** test-1 rho on CD4 Naive and on CD8 EM (1 + 2) pseudobulks for the pairs with
   >= 20 kept cells of the subtype at both time points (computed only if >= 4 such pairs); reported.
7. **GSE116256 descriptive (2.2):** scores of the diagnosis and post-treatment samples in a Seq-Well space built like R on
   the six S8 `desc` diagnosis samples (z parameters from them), per-patient differences; no test.

## 6. Known answers (phase 1; no MRD score is paired with its own diagnosis score here)

| id | what | pass | consequence of a fail |
|---|---|---|---|
| **K2** metadata | the 12 samples' patient_id, timepoint_class (Diagnosis / Post_treatment_unspecified), library_batch and mrd_status in `meta_GSE227903_v2.3.csv` equal the table of 2.1 | all equal | table corrected and recorded before anything else |
| **K0** eligibility | 6 pairs with >= 20 kept T at both time points and one batch, from the S4 meta | exactly the six of 2.1 | stop: the registered set is wrong |
| **K3** coverage | \|primary ∩ U_R\| / 370 and \|test-2 set ∩ U_R\| / 358 | both >= 0.80 | stop |
| **K1** regression | the exact S8 recipe (`logz` over the 7 diagnosis samples 1216, 1886, 3904, 40389, 4978, 5143, 5750; z over the 7) reproduces `ws_s8_scores_rep2.csv` score_primary | Spearman > 0.99 and max \|diff\| < 0.05 | stop: scoring machinery changed |
| **K6** covariates | the table of 5.4 is written; no MRD sample has hb_share_nonery > 0.05 or gran_share_lymph > 0.01 | reported; a breach demotes that pair to "flagged" in every table (not excluded) | — |
| **K4** planted preservation (machinery) | the 6 MRD samples are assigned to the 6 diagnosis samples by one random permutation pi (SEED; redrawn while pi fixes any patient, i.e. a derangement); the real pairing is never used in phase 1. Plant: on primary ∩ U_R genes, each MRD sample's counts are multiplied by 2^(c x (L_D,pi(i),g - m_g)), i.e. the assigned diagnosis sample's log-CPM deviation is added to the MRD counts (c = 1 full, c = 1/3 third); log-CPM, TMM and z are re-derived with U_R, m_g, s_g fixed. Recovery under the ASSIGNED pairing: full plant -> test-1 rho >= 0.8286 and test-2 M <= 1.5; third plant -> rho and M reported. Unplanted scrambled pairing -> rho and M reported as the negative control (an accidental pass has probability ~0.05 and is not a failure) | full plant passes both | stop: the machinery cannot see a planted preservation |

Reading real MRD counts under the scrambled pairing (K4) exposes the MRD score distribution but not any MRD-diagnosis pairing;
this is the same compromise as S9's K5 (permuted ages). Phase 1 writes `ws_s10_known_answers.csv`, `ws_s10_checks.csv`
and `ws_s10_space.csv` (U_R size, coverage, K1 numbers) before phase 2 starts.

## 7. Evaluation order, gate and wording

Order: K2, K0, K3, K1, K6, K4; then test 1 (rho, exact p); then test 2 (a), then (b); then section 5 in its order. Nothing in
section 5 changes a verdict.

| outcome | reading (wording fixed here; bracketed values from the S10 tables) |
|---|---|
| tests 1 and 2 pass | "In six adults (GSE227903, one 10x batch), the between-patient ordering of the T imprint at diagnosis is preserved at the MRD time point (Spearman [rho], exact p [p]), and each patient's MRD T profile identifies its own diagnosis profile better than [P5-rank] of 1,000 matched random gene sets (mean self-rank [M] of 6). The imprint persists into remission in these adults; with six pairs this is a direction, and the formal conclusion waits for the paediatric pairs." Plus, if 1886 and 4978 (MRD-positive) change the sign of rho when removed: "the preservation depends on the two MRD-positive patients" |
| test 1 passes, test 2 fails only (b) | "the ordering is preserved ([rho]); the profile's self-similarity does not exceed that of random matched genes ([M] vs random P5 [..]) — individuality, not the imprint, identifies the owner"; GATE_S10 FAIL |
| test 1 passes, test 2 fails (a) | "the ordering is preserved; the profile does not identify its owner"; GATE_S10 FAIL |
| test 1 fails | "no evidence that the imprint's ordering persists at n = 6 ([rho], p [p]; a pass needed rho >= 0.83)"; test 2 reported; GATE_S10 FAIL; the claim stays "diagnosis-phase imprint" and the remission half of WP-A2 waits for the paediatric cohorts |
| K0/K1/K3/K4 fail | no pairing is read; correct and record, or stop (section 6) |

GATE_S10 FAIL stops the line for discussion (P5). A pass does not change the S9 wording (direction versus healthy stays in
GSE185381 only; the age share stays a bound).

## 8. Power and limits (stated before the run)

- Test 1 at n = 6 can pass only for rho >= 0.8286 (two adjacent swaps at most). Simulated power (bivariate normal, 20,000
  draws): true rho 0.5 -> 0.10, 0.7 -> 0.21, 0.8 -> 0.32, 0.9 -> 0.52, 0.95 -> 0.69. A fail is therefore weak evidence
  against persistence; a pass is strong evidence for it (exact false-pass 1.7 %).
- Test 2's (a) null is exact; (b) is calibrated by the random sets' own self-recognition, which is expected to be well above
  chance (individuality), so (b) is the demanding half.
- Chemotherapy and haematopoietic recovery change T composition: the CD8 EM share falls and the CD4 Naive share rises at MRD
  in these samples (S4 cell counts). A shift common to all patients leaves test 1 untouched; patient-specific shifts are a
  real limitation, reported through section 5.3, not corrected (n = 6).
- MRD status is known for two patients only (positive); regimens are not in the metadata (`treatment_agents_detail`
  unknown), so no treatment covariate exists.
- One cohort, one batch: a pass says "persists in these six adults"; no healthy-donor or second-cohort comparison is made.

## 9. Outputs, code, smoke run, order of operations

- Script `scripts/ws/34_s10_remission.R` (helpers copied verbatim from 32/33: `black_for`, `logz`, `score`, `matched`,
  `split_random`); SLURM wrapper `scripts/ws/submit_ws_s10_remission.sbatch`. No new dependency; no pipeline stage is rerun
  (GSE227903 and GSE116256 are already through S4); no configuration edit.
- Inputs: `ws_pb_GSE227903.rds`, `ws_pb_GSE116256.rds`, `ws_pb_GSE185381.rds` (3.2), `ws_s4_pseudobulk_meta_*.csv`,
  `ws_s5_imprint_genes.csv`, `ws_s8_t_imprint_ambient_prone.csv`, `ws_s8_scores_rep2.csv`, `ws_s8_scores_desc.csv`,
  `ws_s5_scores.csv` (the 25), `decontx/GSE185381/ws_s1_ambient_blacklist.csv`, `decontx/GSE227903/*__checks.csv`,
  `meta_GSE227903_v2.3.csv`, `gencode_GRCh38_gene_order.txt`.
- Outputs: `results/tables/ws/ws_s10_known_answers.csv`, `ws_s10_checks.csv`, `ws_s10_space.csv` (phase 1);
  `ws_s10_tests.csv`, `ws_s10_scores.csv` (12 + GSE116256 samples: scores, cell counts, composition, covariates),
  `ws_s10_self_recognition.csv` (the 6 x 6 correlations and ranks), `ws_s10_positions_J.csv`,
  `/LARGE1/.../ws_pseudobulk/ws_s10_null_draws.rds`; GATES.md block `A2-S10`.
- `--smoke`: the MRD-to-diagnosis pairing is PERMUTED (a derangement, SEED), 50 random sets, outputs suffixed `_smoke`,
  GATES.md not written. No `_smoke` table may contain a real patient id on an MRD row next to a score (the S9 lesson: every
  fwrite is audited for this before the first smoke).
- Order: this file approved and tagged `prereg-s10-v1` -> script written, smoke-tested, reviewed by independent agents
  against this file -> one real run -> outcome appended as a dated section -> GATES.md, handoff §17 S10 row.

## 10. Rules

No threshold, gene set, patient list, correlation measure or wording is changed after a number is seen; later changes are
appended, dated sections. S11 (relapse pairs: GSE227903 diagnosis -> relapse, GSE201966) is a separate pre-registration
with its own tag; nothing from it is read here. This document deletes nothing and adds no dependency.
