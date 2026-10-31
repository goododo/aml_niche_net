# PRE-REGISTRATION A1-BATCH — is the T imprint separable from the processing factors measurable inside GSE185381?

WP-A1 gaps G05 (pools, capture period, depth, cell quality), the core of G06 (doublets, decontX, sex), lane clustering of
the headline statistic (the part of G03 that cannot wait), and two credibility checks (a bound on T-ALL soup; B and NK
cells as negative-control compartments).

**Status: DRAFT v3, 2026-10-31 (v1 → v2 after review round 1; v2 → v3 after review round 2; section 10). Not binding
until the user approves it and the tag `prereg-a1-batch-v1` is set. No quantity declared below has been computed. Later
changes are appended as dated sections.** Order (user decision 2026-10-28): A1-G02 (PASS) → this amendment → G01 → G04.
Handoff v2.0 §17 S15.

## 0. Question, scope, and what each outcome lets the paper say

S5 found a T-cell imprint in 15 adult AML vs 10 healthy donors of GSE185381. A1-G02 showed that it exists out of sample:
cross-fitted gap 1.378 (1.98 out-of-fold healthy-SD), full-pipeline permutation p 0.008, leave-pair-out AUC 0.93. The
two arms were obtained and processed differently. Healthy marrow was bought (STEMCELL or StemExpress; Lasry et al. 2023
Methods, recorded outside the repository in `reports/2026-10_next_directions/handoff_v2_review/review_results.json`);
AML came from a tissue bank. 7/10 healthy vs 5/15 AML libraries carry 2019 file-name dates. Every healthy donor was
hashed into a pool with leukaemia samples.

What this amendment can and cannot do:
- **It cannot separate disease from sample source** (bought vs tissue bank): the source differs for every healthy donor.
  The only on-disk contrast whose design separates them is G04(a) (Chen2023: one centre, one protocol, allograft donors;
  5 vs 4, power about 0.5 at the observed effect, so a fail there is uninformative).
- **From the design alone** (metadata, least squares): a gene whose only signal is the capture period, the UMI depth of
  T cells or "late AML" reaches at most |t| 1.69, 1.49 or 2.60 for status, below the smallest P370 |z| of 2.78. None of
  these factors can carry the whole signal of an imprint gene; they can only weaken it. The one within-study pattern
  that can carry it (cap 4.18) is the **Control1-5 group**: the five donors with 2019-08 file dates, both
  metadata-female healthy donors, the deepest healthy T cells, the least foreign leukaemic soup in their lanes, and the
  only T-lineage-leukaemia pool partners.
- Adjusting for a factor cannot tell confounding from over-adjustment when the factor can be a consequence of disease
  (depth, doublet score, mitochondrial fraction; period nests Lasry's inflammation label). Wording follows this.

So two readouts can stop the line (is the imprint carried by the Control1-5 group? does the out-of-sample test survive
pool-level clustering?), four covariates can only weaken it, and everything else is descriptive.

| outcome | the manuscript may say | next |
|---|---|---|
| PASS | "None of the six listed within-study processing factors, each tested one at a time, explained or weakened the frozen direction (R a-b; block-respecting p c). Sample source, storage, sorting method and chemistry version are constant within arms or unrecorded, and were not tested." | G01 |
| PARTIAL | "The direction is partly sensitive to <factor> (R x; share carried f̂ y)." | G01, carrying <factor> (the same covariate or the same removal) in numerator and denominator |
| FAIL | "Inside GSE185381 the direction is not separable from <factor>." | stop for discussion: methods/negative paper; a composition paper only if G01 shows the composition difference survives <factor> |

## 1. Already known before this draft (disclosed)

Every readout below needs fits that have not been run. Public tables give proxies for most of them, so no arm is
presented as a blind test.
- **Out-of-fold scores** (`ws_a1_g02_oof_scores.csv`) were read in full. Healthy: Control4 -2.20, Control0005 -1.07,
  Control0082 -1.01, Control4003 -1.00, Control5 -0.98, Control3 -0.93, Control2 -0.44, Control0004 -0.05, Control1
  -0.03, Control0058 +0.17. AML from -0.94 (AML1371) to 2.12 (AML0310). Computed from them: the Control1-5 donors
  average -0.92 and the other five -0.59 (a crude proxy of G1); within-AML means by calendar year and Lasry label are
  2019 low 0.77, late low 0.85, late high 0.48 (this **is** the Lasry-label comparison, so it is disclosed here and not
  run). Other proxies computable from this table with public covariates were not computed (within-day contrasts, period
  and depth regressions, leave-pool-out); the round-1 review quoted AML2451 (-0.41) and Control0004 (-0.05) side by side.
- **In-sample scores**: `ws_s5_scores.csv` (an input) holds each sample's S5 score over the 402 genes next to every MAIN
  covariate; its score column was not read in full for this draft. `GATES.md:237` prints four AML scores
  (1.07/0.48/0.83/0.53) and the healthy range (-2.27 to -0.43). With the covariates it gives close proxies of every
  covariate and subset arm.
- **Leave-pair-out pairs** (`ws_a1_g02_lpo_pairs.csv`): first rows seen (AML0048 vs Control0004, 0005, 0058, 0082, all
  correctly ordered; same-day pair AML0048 vs Control0082 1.54 vs -0.88); misordered pairs counted by sample (AML1371 4,
  AML2451 4, AML2123 3; Control0004 3, Control0058 3, Control1 3, Control2 1, Control3 1).
- **Fold imprint sizes**: 29 (leaving out AML4340), 75 (Control4), 86 (Control3), 127 (AML0114), others up to 1,075.
- **D2 null** (all 1,000 relabellings): the gap is unrelated to the 2019 share (r 0.008) and to the side of Control1-5
  (r -0.031). The S5 permutation null's split by 2019 capture is in `ws_a1_g02_s5_perm_draws.csv`.
- **S5 tables**: limma-voom arm 678 hits, 383 in the 402 (`ws_s5_sensitivity.csv`). "Multiplexed" arm 405 hits, 308
  overlap: its coding (`27_s5_imprint.R:61`) equals "hashed" in 24 of 25 samples (AML4897 is hashed but coded 0), so it
  is nearly a hashing check already. Same-channel pair AML2451 vs Control0004: concordance 0.639 vs random P95 0.525.
  GATE_A1(4) T shift 0.28 (registered 3 vs 7 donors, run 2 vs 8). The multiplexed covariate removed the NK imprint (91 →
  2 genes). NK and Mono imprints share 22 genes each with P370 (same sign 22/22 and 20/22). `ws_s5_de_B.csv` and
  `ws_s5_de_NK.csv` test 174 and 179 of the P370 genes: a public proxy of D10, not computed.
- **Gene-level reads** (`ws_s5_de_T.csv`): JUN -1.36, NFKBIA -0.87, IER2 -1.27, DUSP4 +2.26, GNA15 +2.27, HMGB3 +1.59;
  ATF4, MTHFD2, BAG3, NDRG1 and PLIN2 were named as up-arm genes linked to processing delay.
- **September T_NK** (2026-09-16): 1 baseline hit, reported as null after its permutation test: AC007952.4 (in P370,
  down, selected in all 25 folds). No other P370 gene was read in the September tables; this amendment reads none.
- **Hashtag and metadata facts** (no expression): section 2; the T-ALL hashtag background of the healthy donors' kept T
  cells (Control3 0.085-0.104 of HTO UMI; others 0.002-0.039); T3a per-cell doublet scores of the 15 AML
  (`ws_doublet_scores_percell.csv.gz`) and their `dbl_mean`.
- **Design-only numbers computed for this draft**: the t caps above; the variance factors of the drop references
  (section 3).

## 2. Inputs and frozen definitions (fixed before any expression is read)

### 2.1 Samples, design, genes, engine

- The 25 S5 T samples in S5 order (`ws_s5_scores.csv`, compartment T). `ws_pb_GSE185381.rds` is subset to these 25
  columns at load (checked: 25 columns, no paediatric sample). `$raw$T`; rounded `$dx$T` for D4; `$raw$B` and `$raw$NK`
  for D10.
- `MAIN = ~ status + hb + gran + sex + subtype_share` (subtype_share = CD8 share of kept T cells), metadata sex
  (Control4 male) unless an arm says otherwise; the sex term is dropped inside any fit whose samples have one sex.
  Engine: 27's edgeR QL branch (`fit_de`, copied as in 35/36): blacklist, filterByExpr(group = status) and CPM > 1 in
  half, TMM; hit rule FDR < 0.10 and |log2FC| > 0.5; imprint = hits minus the 32 ambient-prone genes.
- Frozen gene sets: P370 (156 up / 214 down) and F402; main logFC = the S5 fit. `library_id` from
  `ws_s4_pseudobulk_meta_GSE185381.csv` (blocks holding >= 2 of the 25: P01 4 AML, P02 4 healthy, P15 2 AML, P05 1 + 1).
- SEED = 491638; each reference has its own seed (SEED + a fixed offset per arm, listed in the script header).

### 2.2 File-name dates, period and hashing

Dates are those in GEO file names. They belong to libraries, not donors (5' and 3' libraries of the same donor carry
different dates), and are not documented capture dates.

| date | AML (15) | healthy (10, all hashed) |
|---|---|---|
| 2019-07-01 | AML0310, AML0612, AML3133, AML3762 (pool P01, hashed) | |
| 2019-08-01 | | Control1-4 (pool P02) |
| 2019-08-29 | | Control5 (pool P03) |
| 2019-10-25 | AML2451 (hashed; same channel GSM5613756) | Control0004 (same channel) |
| 2019-11-04 | | Control0058 (pool P06; two hashtags, Healthy_0058 and Healthy_0059) |
| 2020-03-18 | AML0048 (hashed with its undeposited PR, R1, R2 samples); AML3266 (unhashed, no ADT) | Control0082 (pool P12) |
| 2020-07-29 | AML2123 (hashed with 3 unnamed) | |
| 2020-08-26 | | Control4003 |
| 2020-09-11 | AML1371, AML4340 (pool P15, hashed with 4 unnamed) | |
| 2020-10-06 | AML4897 (hashed with 2 unnamed) | Control0005, channel GSM5613781 |
| 2021-03-09 | AML0114 (Z-1 libraries; hashed with 5 unnamed; 163-feature ADT, 5P(BC) pattern); AML0693, AML3730, AML3948 (unhashed) | Control0005, channel GSM5613787: **date not established** |

- Control0005's second channel: GEO titles Control0005_AML009a/b, identical feature references, count1..count4 numbering
  across 2020-10-06 and 2021-03-09, and its 5' library GSM5613799 is dated 2020-10-06. Settling the date needs read
  headers (first reads of SRR16208958 and SRR16208964); SRA run tables carry no sequencing date.
- **Matched days** (D1): 2019-10-25, 2020-03-18, 2020-10-06 (Control0005's GSM5613781 only): 4 AML vs 3 healthy units.
- **Period** (W1): early = libraries dated up to 2020-03-18 (AML3266 by its date; ADT panels of 69 or 84 features);
  late = from 2020-07-29 (277 features; AML0114 163). Early 7 AML vs 8 healthy (as pool units 4 vs 5, since P01 and
  P02 are single-arm pools); late 8 vs 2. **Calendar year** (D6): 2019 5 vs 7; 2020-21 10 vs 3.
- **Hashed** (D11): all 10 healthy and 11 of 15 AML; unhashed AML0693, AML3266, AML3730, AML3948.

### 2.3 Lane composition of the healthy donors (top-hashtag shares of called barcodes; no expression)

| healthy donor | own share | T-ALL/ETP share | other samples | partners |
|---|---|---|---|---|
| Control1-4 (P02) | 78-91 % | 9-22 % | 0 % | T-ALL 3783, T-ALL 3862, ETP-ALL 9153 (undeposited) |
| Control5 (P03) | 21 % | 8.5-9.0 % | about 70 % | adult AML0361, AML3050, AML2910, one unnamed; T-ALL 4954 |
| Control0004, 0058, 0082, 4003, 0005 | | 0 % | 54-79 % | adult and paediatric AML, MDS, treated, end-of-treatment and relapse samples, unnamed samples |

So the Control1-5 group is not "the T-ALL soup group" alone: four of its five donors have the least foreign soup of all.

### 2.4 Control0005 by channel (D1 only)

Its kept T cells are re-summed per channel from the 04_annotated counts: GSM5613781 (523 cells) and GSM5613787 (411).
The 25 samples' log-CPM is computed exactly as in A1-G02 R2 (TMM over the 25); the two channel units are normalised
with `calcNormFactors(refColumn = R2's reference sample)` and R2's mean library size for the prior count.

### 2.5 Per-sample covariates (computed before any fit; written to a table with each one's correlation with status and its t cap)

- period: early/late (2.2).
- depth: log10 median UMI per kept T cell (`median_umi`).
- mito: median `percent.mt` over the sample's kept T cells (04_annotated meta.data).
- t_dbl: mean doublet score over the sample's kept T cells (W4).

## 3. Readouts and references

- **R (retention of the frozen direction).** R = Σ s_g · lfc_arm,g / Σ |lfc_main,g| over the P370 genes tested in both
  fits; s_g = sign of the main logFC. Reported with R_up, R_down, sign concordance and coverage; F402 alongside. Computed
  for an arm or a reference fit only if >= 0.8 of P370 is tested (D10 waives this; section 4); excluded fits are counted.
- **Selection share.** R does not fall to 0 even when a factor carries the whole gap: the arm inherits the selection
  noise of the frozen logFC, so R ≈ 1 − f(1 − φ), with f the share of the gap the factor carries and φ the selection
  share. A synthetic calibration matched to A1-G02's shrinkage (review round 2) gives φ̂ = 0.40 (q10-q90 0.29-0.51).
  Every R is also reported as **f̂ = (1 − R)/(1 − φ̂)** with its range over φ̂ 0.29-0.51. A narrowed effect is quoted
  as (1 − f̂) × 1.378, never as R × 1.378.
- **Percentile** = (1 + #{R_ref <= R_arm})/(N + 1) for drawn references; for enumerated references that contain the arm,
  the share at or below it. Mid-ranks for ties. A failed or under-covered reference fit is excluded, counted and replaced
  by a declared reserve draw. Every reference is drawn once in a `--draw` step and written before use; the identity draw
  is excluded.
- **References.**
  - G1 (drop Control1-5): the drops of 5 of the 10 healthy donors whose status-coefficient variance factor
    ((XᵀX)⁻¹ under MAIN, design only) is within ±15 % of the real drop's: 40 drops including the real one. The real
    drop's factor is 2.12 times the full sample's (rank 7 of 252), because it removes both female healthy donors.
    Secondary: z = (R − 1)/√(Δv) over all 252 drops.
  - W arms: 1,000 draws that keep each covariate's relation to MAIN: regress the covariate on MAIN, permute its
    residual, add back the fitted part. (A free permutation breaks the covariate's correlation with status and flagged
    3-5 times too often in the round-2 simulation.)
  - G2: block-respecting relabellings (section 4).
- **CF (cross-fitted gap)** only where an arm changes the counts or the samples (D3, D4, D5); never for covariate arms,
  whose effect A1-G02's label-free scoring cannot see. A subset CF is read against the median CF of its reference drops.

## 4. Arms

### Gating: can make the gate FAIL

**G1: the Control1-5 group.** MAIN refit on the 20 samples without Control1-5 → R and f̂; percentile among the 40
matched drops. Any result is worded "the Control1-5 donors carry <share> of the effect", never "T-lineage soup explains
it": the drop also removes 2019-08 dates, both female healthy donors, the deepest healthy T cells and the
least-souped lanes. D9 bounds the soup mechanism.

**G2: lane clustering of the headline statistic.** A1-G02's D2 (the full-pipeline permutation of the cross-fitted gap;
36's fold loop, copied into 37) re-run with relabellings that keep P01, P02 and P15 whole and P05 split. All 22,284 such
labellings with 15 AML are enumerated (checked in K0); 1,000 are drawn uniformly without replacement, the true labelling
excluded (seed SEED + 6). p_block = (1 + #{gap_b >= 1.378})/(1,000 + 1); free D2 p = 0.008 for comparison. 3.7 % of
block labellings recreate at least 2/3 of the true contrast (0.24 % of free ones), so p_block exceeds the free p even
without clustering; in synthetic worlds with a real effect and no clustering it stayed <= 0.10.

### Weakens-only: can make the gate PARTIAL, never FAIL (their t caps are below 2.78)

One covariate at a time is added to MAIN; R and f̂ are read against 1,000 residual-permutation draws.
- **W1 period** (early/late).
- **W2 depth**: log10 median UMI per kept T cell; read as "UMI content per T cell (sequencing depth or cell state)".
- **W3 cell quality**: median `percent.mt` of kept T cells (the per-sample QC ceilings differed by arm: median 15.5 % in
  AML vs 10.9 % in healthy).
- **W4 residual doublets**: t_dbl. Score: scDblFinder 1.20.2 as it runs in this environment; its xgboost classifier
  does not train with xgboost 3.2, so the score is the pre-training heuristic, as at QC and in T3a, and it does not depend
  on the expected doublet rate. Input: the author-labelled barcodes (`Demuxed == TRUE` in
  `01_processed_counts/rds/GSE185381.rds`) of the 35 libraries that hold the 25 samples' cells: 153,935 barcodes, all
  14,609 kept T cells included, no undeposited, post-treatment or relapse barcodes. Per library, `set.seed(SEED)` before
  each. Scores are written only for the 25 samples' 04_annotated cells, with a kept-T flag (other labelled cells,
  including children's, are scoring context only). Union doublet removal already ran at QC, so W4 tests residual
  doublets. It replaces the registered `dbl_mean` arm, which could not run (controls had no scores); that definition (mean
  over all the sample's 04_annotated cells) is reported in D11's table.

### Descriptive: never change the gate

- **D1 capture-day contrast.** v2's B1 machinery on the 3 matched days: training = the 25 minus every sample with kept T
  cells on day d (Control0005 leaves training whole); MAIN fit; imprint_d; held-out units scored with z from the
  training columns; δ_a = s_a − s_h(d) for the 4 AML; W = mean δ_a; W / 1.378; per-day values; the same-channel pair
  alone; p_W over the 12 within-day arrangements (minimum 1/12; reported, never a test). The 2021-03-09 AML against
  Control0005's GSM5613787 channel are reported separately as "date not established". Same date is not same processing:
  only AML2451 and Control0004 share a channel; AML3266 is unhashed without ADT; AML0048 shares its lanes with its own
  undeposited later samples; AML4897 is pooled with 2 unnamed samples. Wording, declared now: W/1.378 >= 2/3 →
  "consistent within capture days"; < 1/3 → "not reproduced within capture days"; otherwise "inconclusive".
- **D2 within-healthy contrast.** The 10 healthy donors only: edgeR QL `~ group + sex`, group = Control1-5 (1) vs the
  other five (0), passed to filterByExpr; readout = share of P370 genes whose group logFC has the sign opposite to the
  imprint (above 0.5: Control1-5 look more healthy-like); exact reference over the 252 splits. Not attributable
  (confounded as G1).
- **D3 companions of G1.** CF without Control1-5, read against the median CF of the 39 matched drops (each a 20-fold
  CF); R without P02 against the 13 drops of 4 healthy donors within ±15 % of its variance factor (P02's drop has the
  largest factor of all 210).
- **D4 decontX.** MAIN on rounded `$dx$T`: R, R_up, R_down, f̂; dx/raw logFC ratio for GNA15 and HMGB3 (named before
  the run); CF / 1.378. decontX saw only the cells in our objects (not the T-ALL cells) and recovers about 60 % of
  spiked ambient RNA.
- **D5 sex.** Control4 coded female: R, CF / 1.378 and the leave-pair-out AUC over 150 pairs (declared in A1-G02 §7;
  reported next to D1/D2 there, never replacing them). Males only by expression (7 AML vs 7 healthy), MAIN without sex:
  R, no percentile.
- **D6 strata.** Early stratum (7 AML vs 8 healthy; 4 vs 5 pool units): R against 500 subsets of 7 AML and 8 healthy in
  which P01 and P02 each appear whole or not at all; the same stratum without P01 and P02 (3 vs 4): R only. MAIN +
  calendar year: R against 1,000 residual-permutation draws. The late stratum (2 healthy donors) is not analysed.
- **D7 leave one pool out.** R without P01 (against the drops of 4 of the 15 AML within ±15 % of its variance factor;
  count printed), P15, P05 and Control5.
- **D8 handling-artefact score.** Per sample, the mean label-free z (over the 25) of `WS_B0_IEG` (26 genes; IER2 is in
  F402 and NR4A1 in the ambient blacklist, so at most 24 enter), with and without KLF2, KLF4 and KLF6 (quiescence
  factors that track naive share): AML vs healthy difference, Spearman with the out-of-fold score within each arm, and
  MAIN + score → R. A drop in R cannot separate handling from biology. If the user approves (section 9), the same for
  the Savage et al. 2021 processing-delay genes.
- **D9 T-ALL soup bound.** In the six P02/P03 lanes, barcodes whose top hashtag is a T-ALL/ETP tag with top share >= 0.8
  (HTO rows of the GEO matrices) are summed into one soup profile (CPM). (i) T-ALL/normal-T ratio of the 214 down and the
  156 up genes against 1,000 random gene sets matched on expression (normal T = the 10 healthy donors' pooled kept-T
  pseudobulk). (ii) Bound: the profile is added to each Control1-5 T pseudobulk at dose = lane T-ALL share × the donor's
  decontX contamination of kept T cells (S1 per-cell estimates; median, and 95th percentile / 0.62); the 25 samples are
  rescored with the frozen genes and the z reference unchanged; the shift of the healthy mean is reported as a fraction
  of the in-sample gap 1.971 (same scoring). Reading: < 0.1 at the high dose → "T-ALL soup cannot carry a material share
  of the gap"; otherwise the share is reported.
- **D10 negative-control compartments.** R_B and R_NK of the frozen P370 direction, refitting S5's B and NK designs (13
  AML + 10 healthy; 12 + 9) on every P370 gene that passes that compartment's filter (about half; the 0.8 coverage rule is
  waived and coverage reported). Reference: 1,000 random gene sets matched on expression in T and in the compartment, R
  computed from the same compartment fit. Reading: R above the random P95 and >= 1/2 → "a sample-level component
  (handling, soup or systemic disease) is present"; inside the random range → "not a generic sample-level artefact".
- **D11 hashing.** MAIN + hashed (2.2): R. A table per sample: hashed, lane soup shares, author kept fraction per lane,
  own-hashtag share, QC mito ceiling, registered-definition `dbl_mean` from W4's scores.

Moved out: the September/S5 hit-count ladder (v2's B3) goes to the G03 step after G04's tag; its urgent question
(lane clustering) is G2. Lasry's label: disclosed in section 1, not run.

## 5. Known answers (run first; no arm result is printed until all pass)

- **K0** Tables built by the script from file names, feature references and HTO rows equal section 2 (dates, period by
  the date rule, hashing, pool partners; lane shares). The block label space has 22,284 labellings (25,716 for 24
  samples, not used). The G1 reference has 40 drops including the real one; the P02 reference 13. The covariate table is
  printed with correlations and t caps.
- **K1** The MAIN refit reproduces the 402 genes; R(MAIN) = 1 within 1e-6. `ws_a1_g02_perm_labels.csv` has md5
  b2f6e70b7e16f4ea374dce14ef6684e7, 1,000 rows of 25-character label strings (read as character) with 15 ones.
- **K2** The two Control0005 channel units sum exactly to its `$raw$T` column.
- **K3** Code paths, with no real arm evaluated: (a) the removal function with nothing removed gives R = 1 (1e-6) and
  CF 1.37816 (± 0.001); (b) the covariate function with no covariate gives R = 1; (c) the covariate function with S5's
  `multiplexed` column reproduces S5's plus_multiplexed arm (405 hits, 308 in the 402); (d) G2's fold loop with the true
  labels reproduces gap 1.378 (± 0.001); (e) every drawn covariate keeps its fitted part (max |difference| < 1e-9).
- **K4** D1: (a) adding a constant to every held-out unit of one day changes W by < 1e-9; (b) +0.5 on the 4 AML units
  raises W by 0.5 (± 1e-9); (c) the 12 arrangement values of W sum to 0; (d) one-hot scores give the hand-computed
  per-day deltas; (e) scoring one sample alone as a "day" with D1's normalisation reproduces its out-of-fold score within
  0.001.
- **K5** W4: barcodes removed as union doublets at QC (in the ingest object, absent from 04_annotated) score higher than
  kept cells of the same library: per-library AUC reported, median AUC > 0.5.
- **K6** D9: per lane, the number of barcodes with a T-ALL/ETP top hashtag equals the 2026-10-28 tally (403, 2,757,
  2,408, 2,558, 1,911, 1,924) before the 0.8 filter; CD3E or CD7 is in the top 1 % of the soup profile by CPM.
- **Calibration report** (after the arms, not a gate): a healthy-ward shift proportional to each gene's main logFC (0.25,
  0.5 and 1.0 × logFC) is planted on Control1-5's P370 counts; the main fit, G1 and its 39 references are refit with P370
  fixed; reported: the smallest plant at which G1's percentile is <= 0.05.

## 6. Decision rules (bracketed values are proposals; the user fixes them before the tag)

- **G1**: explains = R < [0.6] and percentile <= [0.05]; weakens = percentile <= [0.05] and R < [0.85]; otherwise no
  effect.
- **G2**: p_block <= [0.05] → survives lane clustering; in ([0.05], [0.10]] → weakens; > [0.10] → explains ("not
  separable from lane clustering"). Monte Carlo SE with 1,000 draws: 0.007 at p 0.05, 0.009 at p 0.10.
- **W1-W4**: weakens = percentile <= [0.05] and R < [0.85]; they cannot explain.
- **GATE_A1B**: FAIL if G1 or G2 explains. PARTIAL if neither explains and at least one of G1, G2 and W1-W4 weakens.
  PASS otherwise.
- **Operating characteristics** (synthetic, review round 2: G1 plus covariates like W1, W2 and W4, references and
  thresholds as above): no processing effect → PASS 0.86-0.89, PARTIAL 0.11-0.14, FAIL 0; the Control1-5 group carries
  20 % of the gap → PASS 0.78; half → PASS 0.47, PARTIAL 0.33, FAIL 0.19; all → FAIL 0.68-0.92; period or depth carries
  half → PARTIAL 0.94-1.00. G2 and W3 were not simulated; each adds about 0.05 to chance PARTIAL, so with no processing
  effect expect PASS about 0.8. The chance of at least one weakens flag among the six gating readouts is at most 0.30
  (Bonferroni).

## 7. What this amendment cannot do

- Separate the source (bought vs tissue bank), or the vendor of each healthy donor (two vendors; not on disk; the
  Control1-5 naming series may mark vendor or lot).
- Test storage time, thaw-to-capture time, per-sample sorting method (dead-cell depletion or DAPI-negative sort) or 10x
  v3 vs v3.1: unrecorded per sample.
- Match processing: the same file-name date is not the same capture, library type, hashing or ADT staining.
- Establish capture dates (file-name dates; Control0005's second channel not established).
- Attribute an effect of the Control1-5 group to T-lineage soup (confounded with date, sex, depth, soup amount and naming
  series); D9 bounds the soup mechanism only.
- Test doublets that union removal already deleted, or ambient RNA from cells outside our objects (except via D9). The
  authors' SoupX matrices cannot serve as a sensitivity (not a plain correction; cell-called).
- Test factors jointly: each is tested one at a time.
- Tell confounding from over-adjustment for factors that can be downstream of disease (depth, doublets, mito, period).
- Test blood admixture in bought aspirates, or age (bounded by S9).
- Narrow S9's study offsets (0.275, -0.120, 0.838 against about 1.97, biased up by about 0.3): they remain the only
  on-disk measure of how far source and protocol move healthy scores.
- Open the GSE185381 5' pan-T arm (sorted T cells, no blast soup): it would test pooling and soup directly but is fully
  date-confounded (adult AML 5' libraries 2021-04-29; healthy 2019-10-21 to 2020-10-06).
- Change A1-G02's D1/D2.

## 8. Code, outputs, compute (plan; scripts are written after the tag)

- `scripts/ws/37_a1_batch.R`: tables and known answers; G1, G2, W1-W3, D1-D8, D10, D11. Modes `--draw` (every reference,
  per-arm seeds, written before use), `--task=N` (arrays), `--combine`. Helpers copied from 35/36 with explicit inputs
  (counts, log-CPM, design), never globals.
- `scripts/ws/38_a1_batch_cells.R`: loads the ingest object once on a >= 64 GB node (layers checked; no JoinLayers over
  918,864 barcodes); W4 scores; D9's soup profile from the six lanes' GEO matrices. Writes only the 25 samples' per-cell
  scores and one profile vector under `/LARGE1/gr10634/gaozy/aml_niche_net/ws_a1_batch/`.
- Allowlist (exact paths): `ws_pb_GSE185381.rds` (25 columns), `ws_s5_*`, `ws_a1_g02_*`, `ws_s4_*_GSE185381.csv`,
  decontX per-cell tables of the 25, the 04_annotated objects of the 25 samples, the ingest object (38 only), GEO file
  names, feature references and the six P02/P03 lanes' matrices. Never opened: anything under `12_pseudobulk_de/`, any
  Chen2023 file, any paediatric sample's expression, any post-treatment or relapse sample.
- One write helper accepts only `results/tables/ws/ws_a1_batch_*` and `/LARGE1/.../ws_a1_batch/*` and refuses to
  overwrite a non-smoke file. GATES markers `<!-- A1-BATCH start/end -->`. Scripts 02, 22, 26, 27, 33, 35, 36 and
  `12_pseudobulk_de/*` are never run or sourced. No new dependency.
- Compute (one edgeR QL fit ≈ 7 s): G2 1,000 × 25 folds ≈ 50 CPU-h (array of 50); W1-W4 and D6's year arm 5 × 1,000 ≈
  10 CPU-h; D3 CFs over 39 drops ≈ 1.5 h; D6/D7 references ≈ 3 h; D5's 150 pairs ≈ 0.3 h; D10 minutes; 38 ≈ 2 h. About
  65 CPU-h; a few hours of wall time as arrays. The known-answer job runs before any array is submitted.
- Smoke runs: real metadata; synthetic negative-binomial counts drawn around one profile pooled over all samples with
  each sample's real column sum, and S5's main logFC planted on P370 in the AML columns (so fold imprints are non-empty);
  a check stops the run if a real count column enters a smoke fit; no per-sample label-free quantity is written; every
  smoke file ends in `_smoke` and is deleted once the real run has reproduced.

## 9. For the user to decide before the tag

1. **Scope.** v3 is the slim version: two gating arms, four weakens-only covariates, the rest descriptive; the hit-count
   ladder moves to G03. The alternative is v2's full set with the round-2 fixes: about 2-3 more days, 30-40 more CPU-h
   and one more review round, and it can conclude nothing more about disease versus processing. Recommended: v3.
2. **Thresholds** (section 6): G1 explains at R < 0.6 and weakens at R < 0.85 (both with percentile <= 0.05); G2 at 0.05
   and 0.10; W arms weaken at R < 0.85.
3. **D8**: also download the Savage et al. 2021 processing-delay gene list (a few MB; DUSP4, ATF4 and MTHFD2 are on the
   up arm). Recommended: yes.
4. **W4 score**: the heuristic score as at QC and in T3a (recommended), or approve an older xgboost so the classifier
   trains (a dependency change).

## 10. Review history

- **Round 1, 2026-10-28** (two reviewers, v1 → v2): B1 demoted to descriptive (public proxy; 3 healthy donors; 2021
  date not established); calibration rationale corrected; Monte Carlo for infeasible enumerations; K4/K5 rewritten; B3
  references fixed; depth arm added; facts corrected; section 1 expanded.
- **Round 2, 2026-10-29 to 10-31** (five lenses: statistics, facts and blinding, implementation and compute, scientific
  adequacy, synthetic operating characteristics; two skeptics on each of 11 blocking findings, all confirmed real;
  workflow wf_72a14baf-4eb, interrupted twice by session ends, merged by hand). v2's gate did not discriminate (no effect:
  PASS 0.53, PARTIAL 0.36, FAIL 0.12; the Control1-5 group carrying everything: FAIL 0.16-0.24). Changes, v2 → v3:
  covariate references keep the covariate's relation to MAIN; "explains" on the scale complete explanation reaches (R
  inherits the selection share φ ≈ 0.40), with f̂ reported; a size floor for "weakens"; the Control1-5 drop read against
  design-matched drops; B2a and the ladder no longer gate; lane clustering tested on the headline statistic (G2);
  cell-quality covariate (W3) and doublet input and score fixed (W4); T-ALL soup bound (D9) and negative-control
  compartments (D10) added; section 1 completed (all arms have public proxies); facts corrected (AML0048 pooled with its
  later samples; ADT panels 69/84/277/163; lane composition; Control0005's date evidence; AML-side hashing); smoke runs on
  synthetic counts with real metadata; write guard, allowlist, per-arm seeds, sampler by enumeration; compute restated.
  Not adopted: SRA run tables for Control0005's date (they hold no dates); a gating hashing arm (S5's arm already covers
  24 of 25 samples).
