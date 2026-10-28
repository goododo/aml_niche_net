# PRE-REGISTRATION A1-BATCH — does any processing factor measurable inside GSE185381 explain the T imprint?

WP-A1 gaps G05 (capture day, pools, period, depth), the core of G06 (decontX, doublets, sex), the blind-safe part of G03
(library blocking on S5's own samples) and G08(a) (Lasry's label).

**Status: DRAFT v2, 2026-10-28 (v1 rewritten after review round 1, section 10). Not binding until the user approves it
and the tag `prereg-a1-batch-v1` is set. No quantity declared below has been computed. Later changes are appended as
dated sections.** Order (user decision 2026-10-28): A1-G02 (PASS) → this amendment → G01 → G04. Handoff v2.0 §17 S15.

## 0. Question and scope

S5 found a T-cell imprint in 15 adult AML vs 10 healthy donors of GSE185381; A1-G02 showed that it exists out of sample
(cross-fitted gap 1.378 = 1.98 out-of-fold healthy-SD, full-pipeline permutation p 0.008, leave-pair-out AUC 0.93).
The two arms were obtained and processed differently. Healthy marrow was bought (STEMCELL / StemExpress; Lasry et al.
2023 Methods, via `reports/2026-10_next_directions/handoff_v2_review/review_results.json`; the curated metadata does not
record it); AML came from a tissue bank. 7/10 healthy vs 5/15 AML were captured in 2019 (GEO file-name dates). Every
healthy donor was hashed into a pool with leukaemia samples, five of them with T-lineage leukaemias.

**This amendment cannot decide disease versus processing.** Inside GSE185381 the source differs for every healthy
donor, and only three healthy donors share a capture day with an AML sample. The decisive processing test is G04(a)
(Chen2023: one centre, one protocol, allograft donors). What this amendment can do is ask, factor by factor, whether a
processing factor that varies inside the study **explains** the imprint: does adjusting for it, or removing the samples
it marks, remove most of the frozen direction, beyond what a random adjustment or removal of the same size does?

Outcomes: PASS (no factor explains or weakens it) or PARTIAL (some factor weakens it; the effect is narrowed) → G01 next.
FAIL (a factor explains it) → stop for discussion (methods/negative or composition paper).

## 1. Already known before this draft (disclosed)

Public tables from A1-G02 (commit a7e9858) and S5, and what was read from them:
- **Out-of-fold scores** (`ws_a1_g02_oof_scores.csv`) were read in full. Healthy: Control4 -2.20, Control0005 -1.07,
  Control0082 -1.01, Control4003 -1.00, Control5 -0.98, Control3 -0.93, Control2 -0.44, Control0004 -0.05, Control1
  -0.03, Control0058 +0.17. AML from -0.94 (AML1371) to 2.12 (AML0310). From these, the five T-lineage-pool donors
  average -0.92 and the other five -0.59 (computed: this is a crude proxy of B2), and the within-AML means by calendar
  year and Lasry label are 2019 low 0.77, late low 0.85, late high 0.48 (B10). The within-day out-of-fold contrast (a
  proxy of B1) is computable from this table and was **not** computed.
- **Leave-pair-out pairs** (`ws_a1_g02_lpo_pairs.csv`): the first rows were seen (AML0048 against Control0004, 0005,
  0058, 0082: all correctly ordered; AML0048 vs Control0082 is a same-day pair, 1.54 vs -0.88), and the misordered pairs
  were counted by sample (AML1371 4, AML2451 4, AML2123 3; Control0004 3, Control0058 3, Control1 3, Control2 1, Control3
  1). So no misordered pair involves Control0005 or Control0082, and at most one same-day pair (AML2451 vs Control0004)
  is misordered. The table holds a leave-pair-out score for all 8 same-day pairs.
- **Fold imprint sizes**: leaving out AML4340, Control4, Control3 or AML0114 leaves 29, 75, 86 or 127 genes (others
  up to 1,075).
- **D2 null tail**: unrelated to the 2019 split (r 0.01) and to the side of the T-lineage-pool donors (r -0.03).
- **S5 results already on disk**: limma-voom arm 678 hits, 383 of the 402 (`ws_s5_sensitivity.csv`; = rung 3a below,
  used as a known answer); "multiplexed" covariate arm 405 hits, 308 overlap (coded "more than one channel OR pooled",
  `27_s5_imprint.R:61`, a disclosed deviation); same-channel T pair AML2451 vs Control0004, sign concordance 0.639 vs
  random P95 0.525 (`ws_s5_same_channel.csv`); GATE_A1(4) T shift 0.28 (registered as 3 vs 7 donors, run as 2 vs 8).
- **Gene-level reads** (`ws_s5_de_T.csv`): JUN -1.36, NFKBIA -0.87, IER2 -1.27, DUSP4 +2.26 (2026-10-28 addendum), so
  the AML-vs-healthy direction of the IEG list (B10) is partly known.
- **September T_NK** (2026-09-16): its single baseline hit, reported as null after the permutation test, is
  AC007952.4, which is in P370 (down) and selected in all 25 folds. No other P370 gene has been read in the September
  tables, and this amendment reads none (that read waits for G04's result).
- **Metadata** (no expression): pool partners and hashtag counts (section 2.3); Control4 is female (A1-G02 §7); union
  doublet removal at QC for 25/25; `$dx$T` exists for all 25 and is non-integer.

Every readout below needs fits that have not been run, but B1 and B2 have public proxies. They are therefore
descriptive (B1) or read against calibrated references (B2), not presented as blind tests.

## 2. Inputs and frozen definitions (fixed before any expression is read)

### 2.1 Samples, design, genes, engine

- The 25 S5 T samples in S5 order (`ws_s5_scores.csv`, compartment T). Pseudobulks `ws_pb_GSE185381.rds` (`$raw$T`;
  `$dx$T` rounded to integers for B4; `$raw$NK` for rung 4). `MAIN = ~ status + hb + gran + sex + subtype_share`
  (subtype_share = CD8 share of kept T cells), metadata sex (Control4 male) unless an arm says otherwise. Engine: 27's
  edgeR QL branch (`fit_de`), copied verbatim as in 35/36; blacklist as 27; hit rule FDR < 0.10 and |log2FC| > 0.5;
  imprint = hits minus the 32 ambient-prone genes.
- Frozen gene sets: P370 (156 up / 214 down) and F402; main logFC = the S5 fit. `library_id` from
  `ws_s4_pseudobulk_meta_GSE185381.csv` (blocks with >= 2 of the 25: P01 4 AML, P02 4 healthy, P15 2 AML, P05 1 + 1).
- Depth: median UMI per kept T cell (`median_umi`, same table). SEED = 491638.

### 2.2 Capture day and period (dates are the dates in GEO file names; they are not documented capture dates)

| file-name date | AML (15) | healthy (10) |
|---|---|---|
| 2019-07-01 | AML0310, AML0612, AML3133, AML3762 (pool P01) | |
| 2019-08-01 | | Control1-4 (pool P02) |
| 2019-08-29 | | Control5 (pool P03) |
| 2019-10-25 | AML2451 | Control0004 (same channel GSM5613756) |
| 2019-11-04 | | Control0058 (pool P06) |
| 2020-03-18 | AML0048, AML3266 (no hashtag, no ADT) | Control0082 (pool P12) |
| 2020-07-29 | AML2123 | |
| 2020-08-26 | | Control4003 |
| 2020-09-11 | AML1371, AML4340 (pool P15) | |
| 2020-10-06 | AML4897 | Control0005, channel GSM5613781 ("count2") |
| 2021-03-09 | AML0114 (Z-1 libraries, different ADT chemistry), AML0693, AML3730, AML3948 | Control0005, channel GSM5613787 ("count1_0.95") |

- **The 2021-03-09 date of Control0005's second channel is doubtful.** Its two channels are titled Control0005_AML009a/b,
  have identical feature references, and are numbered count1 (2021-03-09) and count2 (2020-10-06) next to count3
  (2021-03-09) and count4 (2020-10-06): probably one four-channel capture with files written on two dates. The AML
  libraries of 2021-03-09 use other formats. So **2021-03-09 is not treated as a shared capture day.**
- **Matched days** (used by B1): 2019-10-25, 2020-03-18, 2020-10-06 (Control0005 channel GSM5613781 only): 4 AML vs 3
  healthy units.
- **Period** (primary): early = libraries dated up to 2020-03-18 (AML3266, which has no ADT, by its date); late = from
  2020-07-29. Early 7 AML vs 8 healthy; late 8 vs 2. This is also the ADT panel switch (small panels of 67-84 features
  before, 276-277 after; AML0114 a third chemistry). **Calendar year** (sensitivity): 2019 5 vs 7; 2020-21 10 vs 3.

### 2.3 Pool partners of the healthy donors (feature references; hashtag counts 2026-10-28; no expression)

| healthy donor | captured with | soup type |
|---|---|---|
| Control1-4 (P02) | T-ALL 3783, T-ALL 3862, ETP-ALL 9153 (not deposited; 8.8-22.2 % of called barcodes) | T-lineage |
| Control5 (P03) | T-ALL 4954 (not deposited; 8.5-9.0 %), AML0361, AML3050, AML2910, one unnamed | T-lineage |
| Control0004 | AML2451, two undeposited AML | myeloid |
| Control0058 (P06) | paediatric PAUMTZ, PAXMIJ and their undeposited end-of-treatment and relapse samples; its kept cells carry two hashtags (Healthy_0058, Healthy_0059), both mapped to Control0058 | myeloid |
| Control0082 (P12) | paediatric AML022, AML052 and two treated samples of them | myeloid |
| Control0005 | paediatric AML009, two unnamed | myeloid |
| Control4003 | two undeposited MDS, two undeposited AML | myeloid |

### 2.4 Control0005 by channel

For B1 only, Control0005's kept T cells are re-summed per channel from the 04_annotated counts: GSM5613781 (523 cells)
and GSM5613787 (411 cells).

## 3. Readouts and references

- **R (retention of the frozen direction).** R = Σ s_g · lfc_arm,g / Σ |lfc_main,g| over the P370 genes tested in both
  fits; s_g = sign of the main logFC. Reported with R_up, R_down, sign concordance, coverage; F402 alongside. R is
  computed only if >= 0.8 of P370 is tested in the arm. R is compared only between fits of the same engine.
- **Why references are needed.** If a factor has no effect, R stays near 1 in expectation (a subset or an added null
  covariate inherits the selection bias of the full fit), but its spread grows as samples are removed or covariates
  added. Each gating R is therefore read against a reference of the same size:
  - covariate arms: 200 fits with the covariate permuted across the 25 samples (seeded);
  - sample-removal arms: every subset with the same AML/healthy counts when there are <= 1,365, otherwise 500 seeded
    Monte Carlo subsets;
  - cell-removal arms (B5 iii): 200 removals of the same number of kept T cells per sample, at random.
  The arm's **percentile** = share of reference fits with R at or below the arm's R.
- **CF (cross-fitted gap within the arm)**: A1-G02 R2's procedure on the arm's samples and design. CF falls when
  training samples are removed, so subset arms read **CF_ratio = CF_arm / median CF over the reference subsets**; arms
  on the same 25 samples read CF / 1.378. A narrowed effect is quoted as CF_ratio × 1.378.

## 4. Arms (G = gating)

**B1 (descriptive): capture-day-matched contrast, cross-fitted.** For each matched day d: training = the 25 samples
minus every sample with kept T cells on d (Control0005 leaves training whole). MAIN fit on training; imprint_d = hits
minus the 32. Label-free log-CPM over the 25 samples plus the two Control0005 channel units (gene filter over the 25 as
in R2; TMM over 27 columns). Held-out units scored with z from training columns. δ_a = s_a − s_h(d) for the 4 AML
samples; W = mean δ_a; W / 1.378; per-day values; the same-channel pair alone. p_W over the 2 × 3 × 2 = 12 within-day
arrangements (minimum 1/12, compared with a 1e-9 tolerance), reported, never interpreted as a test. The 2021-03-09 AML
samples against Control0005's GSM5613787 channel are reported separately, labelled "date not established". Same day is
not same processing: only AML2451 and Control0004 share a channel; the other healthy units sit in hashed pools while
their AML partners are unhashed (AML3266, also without ADT) or singly hashed (AML0048) or pooled with undeposited
samples (AML4897). The public leave-pair-out proxy (section 1) is disclosed next to B1; disagreement between them is
reported, not resolved.

**B2 (G): T-lineage pools.**
- B2a: CF on the 20 samples without Control1-5 (15 AML vs the 5 myeloid-soup donors); CF_ratio and percentile against
  all 252 drops of 5 of the 10 healthy donors, each with its own 20-fold CF.
- B2b: R without Control1-5 against the same 252 drops (refit only).
- B2c: R without P02 (Control1-4) against all 210 drops of 4 healthy donors.
- B2d (descriptive): the 10 healthy donors only, edgeR QL `~ soup + sex`, soup coded T-lineage = 1, myeloid = 0. Readout:
  share of P370 genes whose soup logFC has the sign opposite to the imprint (above 0.5 = T-lineage-soup donors look more
  healthy-like than myeloid-soup donors); exact permutation over the 252 splits. Confounded with capture date (all five
  in August 2019), sex (every female healthy donor is in P02) and depth (P02 median UMI 2,434-3,948 vs 1,068-1,855).

**B3 (for G03; read separately, not gating): library blocking and the September design, S5's samples only.**
- rung 1: S5 as run (edgeR QL, MAIN). rung 2: edgeR QL, `~ status`. rung 3a: limma-voom, MAIN, no block (known answer,
  section 5). rung 3b: limma-voom, MAIN, `duplicateCorrelation(block = library_id)`, two voom passes as in
  `12_pseudobulk_de/02_de_limma.R`. rung 4: rung 3b on T + NK sums (subtype_share unchanged). rung 5: rung 3b on the
  September projection-only T_NK bin (24 samples: 14 AML, 10 healthy; AML0310 has no September file). rung 6: rung 5 with
  the September design `~ arm` (blocked) and hit rule (adj. p < 0.05 and |logFC| >= 1) on the September gene universe
  (gene names shared by Chen2023 and GSE185381, minus chrY and XIST, filterByExpr on the design); the September depth term
  is not used (stated difference).
- Per rung: hit count under the rung's rule, and two permutation p values from 500 relabellings each, with the full
  pipeline (including duplicateCorrelation) re-run per relabelling: (i) free relabellings (rungs 1-4: the first 500 D2
  relabellings; rungs 5-6: 500 new relabellings of the 24 with 14 AML, drawn once with SEED + 5 and written before use);
  (ii) block-respecting relabellings (P01, P02 and P15 move as units, P05's two samples separately; 500 drawn with SEED + 6
  from the labellings that keep 15 AML, or 14 for rungs 5-6). A rung **keeps** the contrast if both p <= 0.05, **loses**
  it if either p > 0.10, otherwise ambiguous.
- Reading: first rung that loses it. Rung 2 → the covariates are needed (reported). Rung 3b → lane-level clustering is a
  validity problem for S5: the imprint is reported as exploratory, and this is carried to the GATE as FAIL. Rungs 4-6 → a
  design difference (cell set, bin, hit rule, universe), reported with the September null alongside; the attribution of
  the September null to pooling with children and Chen2023 waits for G04.

**B4 (G): decontX-corrected counts.** MAIN on rounded `$dx$T`: R, R_up, R_down; dx/raw logFC ratio for GNA15 and HMGB3
(named before the run); CF / 1.378. No reference (same samples, same design). decontX recovers about 60 % of spiked
ambient RNA, so a pass is partial evidence.

**B5 (G for i and iii): doublets.** scDblFinder 1.20.2 per library on all deposited barcodes of the 35 libraries that
hold kept T cells of the 25 (`01_processed_counts/rds/GSE185381.rds`, split by GSM prefix), seed SEED, expected doublet
rate by the project's T3a rule (`02_doublet_scores.R`), run on a node with >= 64 GB. Other cells in those libraries
(children, undemultiplexed and undeposited samples) serve only as context; scores are written for the 25 samples' kept
cells only, with the call rate per library. (i) MAIN + t_dbl (mean score over the sample's kept T cells): R vs 200
permuted-covariate fits. (ii) MAIN + the mean new score over all the sample's 04_annotated cells (the definition of the
registered `dbl_mean`, `26_s4_pseudobulk.R`): R, descriptive; it replaces the registered arm, which could not run.
(iii) kept T cells called doublet dropped, raw T sums rebuilt: R vs 200 random removals of the same number of cells per
sample; CF / 1.378.

**B6: sex.** (i) (G) Control4 coded female: R and CF / 1.378, no reference; plus the leave-pair-out AUC over 150 pairs
(declared in A1-G02 §7; reported next to D1/D2, never replacing them). (ii) (descriptive) males only by expression (7 AML
vs 7 healthy), MAIN without sex: R against 500 Monte Carlo subsets of 7 AML and 7 healthy.

**B7 (G for i and ii): period.** (i) MAIN + period: R vs 200 permuted-covariate fits. (ii) early stratum only (7 AML vs 8
healthy): R vs 500 Monte Carlo subsets of 7 and 8. (iii) (descriptive) MAIN + calendar year, R vs permuted; late stratum
(8 vs 2), R only.

**B8 (G): depth.** MAIN + log10(median UMI per kept T cell): R vs 200 permuted-covariate fits.

**B9 (descriptive): leave one pool out.** R without P01 (vs all 1,365 drops of 4 of the 15 AML), P02 (= B2c), P15, P05
and Control5 (P03).

**B10 (descriptive): handling-artefact score.** Per sample, the mean label-free z (over the 25) of the project IEG and
heat-shock list (`WS_B0_IEG`, 26 genes, minus IER2, the one F402 member): AML vs healthy difference, Spearman with the
out-of-fold score within each arm, and MAIN + IEG score: R vs 200 permuted. A processing-delay set (Savage et al. 2021)
runs the same way only if the user approves downloading it.

**B11 (descriptive): Lasry's label.** Out-of-fold scores in three groups by calendar year (2019 low: AML0310, AML0612,
AML2451, AML3133, AML3762; later low: AML0048, AML3948; later high: 8), already readable (section 1). No p value.

## 5. Known answers (run first; no arm result is printed until all pass)

- K0: capture-day, period and pool tables built by the script from file names and feature references equal sections
  2.2-2.3 (AML3266 by date).
- K1: rung 1 reproduces the 402 genes; R(rung 1) = 1; the D2 relabelling file has its committed md5.
- K2: the two Control0005 channel units sum exactly to Control0005's `$raw$T` column.
- K3: rung 3a reproduces S5's limma-voom arm (678 hits, 383 in the 402).
- K4 (B1 pairing logic, with each day's fold genes and z references frozen): (a) adding a constant to every held-out
  unit of one day changes W by < 1e-9; (b) shifting the 4 AML units' scores by +0.5 raises W by 0.5 (± 1e-9).
- K5: the no-drop / unpermuted reference reproduces R = 1 and CF = 1.378.
- Calibration report (after the arms, not a gate): planting a healthy-ward shift of 0.1, 0.2 and 0.3 log2 on the P370
  counts of Control1-5, the smallest plant that moves B2b's percentile to <= 0.05.

## 6. Decision rules (values in brackets are proposals; the user fixes them before the tag)

- **Gating readouts:** B2a, B2b, B4, B5(i), B5(iii), B6(i), B7(i), B7(ii), B8.
  - Referenced arms: **explains** = R (or B2a's CF_ratio) < [1/3] and percentile <= [0.05]; **weakens** = percentile <=
    [0.05] and R >= [1/3]; otherwise **no effect**.
  - Unreferenced arms (B4, B6 i): **explains** = R < [1/3] or CF / 1.378 < [1/3]; **weakens** = R or CF / 1.378 in
    [[1/3], [2/3]); otherwise no effect. In B4, R_up < [1/3] with R_down >= [2/3] also counts as weakens: the claim
    narrows to the down arm and GNA15/HMGB3 leave the interpretation.
- **GATE_A1B:** FAIL if any gating readout explains, or B3 loses the contrast at rung 3b. PARTIAL if none explains and at
  least one weakens: the narrowed effect = the smallest weakening R or CF_ratio × 1.378, and every later step uses it.
  PASS otherwise.
- **Expected rates with no processing effect:** nine gating readouts, seven referenced at 0.05; they are correlated, so
  the chance of at least one "weakens" flag is below 1 − 0.95^7 = 0.30 but not small. PARTIAL is therefore a common
  outcome without any real processing effect; the report gives the reference distributions so the reader can judge.
  "Explains" needs R < 1/3 as well; its chance rate is reported from the references.
- B1, B2d, B3 (except rung 3b), B5(ii), B6(ii), B7(iii), B9, B10, B11 are descriptive and never change the gate.

## 7. What this amendment cannot do

- Separate the bought (healthy) from the tissue-bank (AML) source: constant in every day and pool.
- Match processing: same file-name date is not the same capture or the same library type (section 4, B1).
- Establish capture dates: they come from file names; Control0005's second channel is probably mis-dated.
- Test T-lineage soup apart from August-2019 capture, female sex and higher depth (B2d).
- Test doublets that union removal already deleted, or ambient RNA that decontX cannot model.
- Account for age (Control0005 is 19-26 years old; healthy donors are 19-53 vs AML mostly >= 60; bounded by S9).
- Change D1/D2 (A1-G02 §7).

## 8. Code, outputs, compute (plan; scripts are written after the tag)

- `scripts/ws/37_a1_batch.R`: tables and known answers; B1-B4 and B6-B11; `--task=N` for the arrays (B2a's 252 CF runs;
  B3's permutations) and `--combine`. Helpers copied from 35/36, not sourced. File allowlist: only the 24 adult September
  files named in `ws_s5_scores.csv`; never `results/tables/12_pseudobulk_de/de/`, `sensitivity_pooled_bh.csv`,
  `hitlist_frozen.csv` or `validation/`.
- `scripts/ws/38_a1_batch_doublets.R`: scDblFinder per library → per-cell scores for the 25 samples under
  `/LARGE1/gr10634/gaozy/aml_niche_net/ws_a1_batch/`.
- SLURM wrappers. Outputs prefix `results/tables/ws/ws_a1_batch_`; GATES.md block `A1-BATCH`; outcome appended here.
  Scripts 27, 33, 35 and 36 are never re-run. No new dependency.
- Compute (one edgeR QL fit ≈ 7 s): B2a 5,040 fits ≈ 10 CPU-h; B3 rungs 3b-6 with two references ≈ 4 × 1,000 blocked voom
  fits at ≈ 15-20 s ≈ 17-22 CPU-h, plus edgeR rungs ≈ 2 h; Monte Carlo and permuted references ≈ 3,500 fits ≈ 7 h;
  scDblFinder 2-4 h at >= 64 GB; the rest under 1 h. About 40 CPU-h, run as arrays (a few hours of wall time).
- Smoke runs use permuted status labels and write `_smoke` files only; no real covariate is written next to a smoke
  score.

## 9. For the user to decide before the tag

1. Threshold values in section 6 (proposed 1/3, 2/3 and 0.05, the project's existing lines).
2. B1 as descriptive (recommended, because of its public proxy, its 3 healthy donors and the uncertain 2021 date).
3. B10 with the IEG list only, or also download the Savage 2021 processing-delay list (a few MB).
4. Optional: fetch SRA run metadata (flowcell and run dates) for Control0005's two channels to settle the 2021-03-09
   date (a network read; needs approval). The draft does not depend on it.

## 10. Review history

- **Round 1, 2026-10-28** (two independent reviewers, read-only). Changes from v1: B1 demoted from primary to descriptive
  (public leave-pair-out proxy; 3 healthy donors, w SD ≈ 0.39 without any batch effect; the 2021-03-09 pairing not
  established, leaving 12 arrangements); the calibration rationale corrected (R does not fall by itself; CF does, so CF
  is read relative to its reference); infeasible exact enumerations replaced by Monte Carlo; K5 replaced (it computed the
  real B2b and passed by construction) and K4 made exact; B3 references fixed (24-sample rungs need their own
  relabellings; block-respecting references added; "pseudoreplication" replaced by lane-level clustering); depth arm
  added; facts corrected (pool partners of Control0005 and Control0058, ADT panels, AML3266, fold sizes, September
  wording); section 1 expanded with everything already readable; GATE logic made consistent and its chance rates stated.
