# PRE-REGISTRATION — S9: is the AML-vs-healthy direction of the T-cell imprint an age effect? (older healthy donors)

Written 2026-10-23; reviewed and approved by the user the same day. **BINDING with tag `prereg-s9-v1`**; later changes are
appended, dated sections, never edits of the text above them. Nothing below has been run on the three new cohorts; no age
relation is read on them before the tag. Context: handoff v2.0 §16.2 item 14(b) (older healthy donors before WP-A2), §16.2 item 6, S5 §3.4 ("AML-or-age"),
S8 §9 (claim carried forward: "direction versus healthy and age not yet replicated"). Label note: handoff §17's "S9" is the
remission pre-registration (appendix A); this document keeps the name S9 for the age/direction control and the remission
row becomes S10 (S10–S13 shift) in the §17 row appended for this work (user decision 2026-10-23).

Design reviewed before writing by three independent adversarial reads, and the written draft by two more; the changes they
forced are recorded in section 9 so that the reasoning is auditable.

## 1. Question and hypotheses

The T imprint (S5) separates adult AML at diagnosis from healthy donors in GSE185381 by 3.87 healthy-SD on the primary 370-gene
set (section 2). The healthy donors are 19–53 years old and the eligible AML patients mostly ≥ 60, so every S5 statement is
"AML-or-age". S8 replicated the imprint's *structure* in two independent AML cohorts but could not touch direction or age.
S9 asks whether the *direction* is an age effect, using healthy marrow across the adult age range.

- **H_AML**: the direction is disease. In healthy marrow the imprint score does not rise with age; healthy donors ≥ 60 do not
  score like AML patients.
- **H_age**: the direction is ageing. The score rises with age in healthy donors by about the full AML–healthy gap between the
  median ages of the eligible healthy (32.5) and AML (64) samples of GSE185381, i.e. about 1.2 healthy-SD per decade; healthy
  donors ≥ 60 score like AML patients.

The yardstick is circular by construction (the 370 genes were chosen on the same 25 samples, AUC 1.0); section 3.5 estimates its
inflation and the gate uses the deflated gap. The gate bounds the age share of the direction; it does not exclude an age
component (section 5 wording).

## 2. Inputs (fixed)

### 2.1 Gene sets and the defining cohort

| item | value |
|---|---|
| imprint | `results/tables/ws/ws_s5_imprint_genes.csv` (md5 977f58b0…), compartment T, imprint_called TRUE: 402 genes (171 up / 231 down) |
| ambient-prone | `results/tables/ws/ws_s8_t_imprint_ambient_prone.csv` (md5 acb2b043…): 32 genes, removed as in S8 |
| **primary set** | **370 genes (156 up / 214 down)**; the 402 set is a sensitivity arm |
| GSE185381 T pseudobulks | `/LARGE1/…/aml_niche_net/ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds` `$raw[["T"]]`, the 25 eligible samples of `ws_s5_scores.csv` compartment T (md5 58374fc6…): 15 adult AML, 10 healthy; unchanged |
| AML log2FC per gene | `results/tables/ws/ws_s5_de_T.csv` (md5 b13de32b…) |
| S5 covariates | `hb`, `gran`, `sex`, `subtype_share` per sample from `ws_s5_scores.csv` (for the cross-fitted yardstick, 3.5) |

**Yardstick (computed 2026-10-23, before any new-cohort score; S5 4.1 space = `logz()` of `32_s8_replication.R` on the 25
samples: genes with CPM > 1 in at least half of the 25, TMM, log-CPM prior count 1, z per gene across the 25; no filterByExpr):**
primary set AML mean 0.788 (n 15), healthy mean −1.182 (n 10), healthy SD 0.509 → **gap_in = 1.970 = 3.87 healthy-SD**, AUC 1.00
(402 set: 1.966 = 3.92). The S9 code must reproduce 1.970 within 0.01 (K6); it is not recomputed. Its own relative standard
error is about 10 % (0.509 × √(1/15 + 1/10) / 1.970), which every CI below is conditional on.

**Age gap for the yardstick: Δage = 31.5 years** = median age of the 15 eligible AML (64; range 32–80; 11 ≥ 60) minus the
median of the 10 eligible healthy (32.5; 19–53: five donors at 19–26 and five at 39–53; S5 §3.4's "6 and 4" is a miscount of
the same ten). Not 35: the 20-adult median 67.5 includes 5 patients without an eligible T pseudobulk.

### 2.2 New cohorts and their roles

| cohort | role | donors used | ages (y) | not used |
|---|---|---|---|---|
| GSE120221 Oetjen 2018 | **DECISIVE** for the age test: 20 donors, one protocol (10x 3' v2, Ficoll MNC, cryopreserved, unsorted, Cell Ranger 2.0), 25 libraries | 20 | young ≤ 53: T 24, W 28, E 30, R 31, F 41, J 43, U 46, B 47, H 50, O 50 (n 10, mean 39.0); middle 54–59: S 56, L 57, G 58, P 58, A 59 (n 5); old ≥ 60: C 60, M 60, Q 66, N 67, K 84 (n 5, mean 67.4) | — |
| GSE227690 Griffin 2023 | descriptive direction check; author labels reported under K1 | BM1–BM5 (13 10x 3' v3/v3.1 libraries, 2–3 per donor) | young: BM1 31, BM2 37, BM3 45; old: BM4 75, BM5 74 | BM6.1 (Seq-Well, 52 y) |
| GSE276360 Rodriguez-Sevilla 2025 | descriptive direction check | YHD658 (matrix prefix HD-36), YHD636 (HD-95), HD2, HD3 | young: YHD658 23, YHD636 24; old: HD2 61, HD3 70 | CCUS-50, CCUS-53, CCUS-PR (clonal cytopenia, not healthy); NK/T-sorted libraries (not downloaded) |

Raw files and provenance: `/LARGE1/…/aml_niche_net/00_raw/public/GEO/<GSE>_RAW/PROVENANCE.txt` (GSE120221_RAW.tar md5 fd94debc…;
GSE227690_BM_counts c30520b6…, BM_metadata 7f360952…; GSE276360_bm_scrna_expr.rds.gz 831d5d3d…). **Donor age and sex
(authoritative for every test):** `/LARGE1/…/aml_niche_net/00_raw/metadata/<GSE>/donor_ages_<GSE>_v1.0.csv` (md5 178ced1c…,
ccf1b395…, 07a23af0…), joined on donor id. Oetjen donor ages for the replicate donors are the age at the first-listed library:
**C = 60** (GEO lists Ck as 59), **S = 56** (Sk1/Sk2 55); second aspirates are pooled into the donor. The `meta_<GSE>_v1.0.csv`
rows written for the pipeline must agree with these tables (K7).

Facts recorded before the run: (i) **Oetjen is one of the six source studies of the BoneMarrowMap reference** (64,813 of
263,159 reference cells; all 22 donor entries; the reference metadata's ExactAge and DonorSex agree with the donor table for every
donor). Its projection is therefore in-sample: `bmm_fine`, `high_error` and marker concordance are optimistic relative to
GSE185381, and K1 is an age known-answer only. Griffin's author labels are compared with ours per cell under K1 and reported;
**no statement about annotation transfer is made from that number.** (ii) During feasibility checks the Griffin *author* labels
showed CD8Naive shares of T cells of 24 % (BM1, 31 y), 10 % (BM2, 37), 6.5 % (BM3, 45), 3.3 % (BM4, 75), 5 % (BM5, 74); disclosed
here, not our labels. (iii) Oetjen library depth: 1,052–7,247 cells, median genes/cell 575–1,390, median percent.mt 0.5–4.2,
CD3E+CD3D+ cells 236 (Q) to 3,601 (K); above every QC gate of `config_qc.R`. (iv) In `ws_s5_de_T.csv` the SD of log2FC over the
primary genes is about 0.97 against about 0.39 for expression-matched random genes (2.5×); this is why A2b uses a slope and not a
correlation (3.2).

### 2.3 Units, eligibility, pipeline conventions (pre-registered procedures)

- **Sample = library** at ingest and QC (25 + 13 + 4 samples; ids GSM-prefixed, e.g. `GSM3396183_T`, `GSM7105550_BM1_BioRep1`,
  `GSM8497031_YHD658`), `Patient_ID` = donor, `GSM` and `Library` per cell (ws/22 runs decontX per GSM channel). Reason: a
  donor-level Sample would triple (C, 3 libraries) or quadruple (S, 4) the doublet prior (`03_per_sample_qc.R` computes it from
  the Sample's n_raw), and C is one of the five donors ≥ 60.
- **Donor pooling inside the S9 script**: a donor's T pseudobulk is the sum of its libraries' kept-cell T pseudobulks from
  `ws_pb_<GSE>.rds` (`$raw[["T"]]`, S4 kept cells: group4 T, `anno_source == "concordant"`, not `high_error`). **Eligibility is per
  DONOR: ≥ 20 kept T cells.** This is a documented deviation from S4's per-sample flag; every other S4 rule (`config_ws.R`
  section 12) is unchanged and the file is not edited.
- Age groups: **young = ≤ 53** (the GSE185381 healthy range), **old = ≥ 60**; Oetjen's 54–59 donors enter the slope only.
- GSE120221 gene axis: Cell Ranger 2.0 symbols are mapped by Ensembl ID to the project's 2020-A symbol set (`REF_FEATURES_2020A`
  in `config_paths.R`, the GSE289435 features table; 32,351 of 33,694 genes map, 10,784 symbols change); unmapped IDs keep the
  Cell Ranger 2.0 symbol unless it collides with a mapped symbol (then suffixed with the Ensembl ID). Primary-gene presence by
  symbol 326/370; after the remap 359/370 (the 11 missing are lncRNAs/pseudogenes, listed in `ws_s9_checks.csv`). The remap table
  is written to `results/tables/00_ingest/GSE120221_gene_remap.csv`. No shared reader is changed.
- Exclusions live in `INGEST_EXCLUDE` (`config_qc.R`): GSE276360 `^CCUS-` (not healthy; required because the dataset enters
  `TP_HEALTHY_DATASETS`, whose rule labels every cell Healthy), GSE227690 `^BM6` (Seq-Well). The excluded libraries keep curated
  rows in the meta tables (as the GSE289435 PDX row) with DISPOSITION rows in `99_admin/12_cohort_roster_check.R`.
- Covariates: **sex** from the donor table (the sex-adjusted fit decides; the unadjusted fit is reported); `hb_share_nonery` and
  `gran_share_lymph` from the S4 meta tables as partial correlations with the score (informational, as S8). Doublet scores are
  reported where T3a has run (NA otherwise).
- Random-set null (as S5/S8): 1,000 matched random gene sets (same number of up and down genes; drawn within bins of GSE185381
  mean-log-CPM decile × detection class (−∞, 0.5, 0.8, 0.95, ∞), `matched()` of `32_s8_replication.R`), SEED 491638
  (`config_paths.R`). 100 sign-scrambled imprints are reported (K5).

### 2.4 Scoring spaces (one row per test; a 1/3 line is applied only within a row)

| test | space | genes | z parameters | numerator | denominator |
|---|---|---|---|---|---|
| A1 gating, A3, K5 | **J-E** (joint with GSE185381; E = Oetjen, Griffin, R-S): the S5 §3.3 blacklist (`black_for()` of 32: HB, MPO/ELANE/AZU1/PRTN3, the 79-gene ambient blacklist, XIST, chrY, MT-, RPL/RPS) removed in both cohorts before the CPM filter and TMM, as `logz()`; universe U_J = remaining genes with CPM > 1 in ≥ half of the 25 GSE185381 samples AND in ≥ half of E's eligible donors; TMM **within each cohort separately** on U_J; log-CPM prior count 1 | primary 370 ∩ U_J (402 ∩ U_J as sensitivity) | per-gene mean m_g and SD s_g over the **25 GSE185381 samples**, applied to both cohorts; score = mean z(up) − mean z(down) | slope_J × 31.5 y; contrast_J = mean(old) − mean(young) | **D_res_J** = mean(AML) − mean(GSE185381 healthy), same space and genes, deflated by r (3.5) where section 3 says so |
| K2 | J-Oetjen at the **library** level: the Oetjen library T pseudobulks form their own DGEList on U_J, TMM, log-CPM prior 1; U_J, m_g, s_g unchanged | same | the 25 GSE185381 samples | library score differences | (D_res_J / r) / 3 |
| A1 reported | **W-Oetjen**: `logz()` on Oetjen's eligible donors alone | same | over the eligible donors | Spearman ρ(score, age); slope per decade | Oetjen young-donor SD (descriptive) |
| A2a reported | **S5 space** (GSE185381 alone, the yardstick space) | 370 minus the removed genes | over the 25 | gap_robust | gap_in = 1.970 |
| 3.4 informational | S5 space | 370 | over the 25 | log2FC sign concordance (4 AML < 60 vs 5 healthy 39–53); ρ(score, age) within AML and within healthy; positions of the four | n/a |
| K6 regression | `logz()` on the 15 GSE185381 AML alone (32's disc set) | primary | over the 15 | scores vs `ws_s8_scores_disc.csv` | Spearman > 0.99, max \|diff\| < 0.05 |
| 3.5 cross-fitted yardstick | S5 space, leave-one-sample-out | fold-selected imprint (FDR < 0.10, \|log2FC\| > 0.5, the 32 ambient-prone removed) | fold mean/SD over the 24 | gap_oof | gap_in = 1.970 |
| A2b modifier | Oetjen edgeR QL (`~ age + sex`) on the eligible donor pseudobulks; GSE185381 log2FC from `ws_s5_de_T.csv` | primary genes in U_J (all carry a log2FC_AML) | n/a | OLS slope b of log2FC_age(31.5 y) on log2FC_AML | random-set P95 of the same slope |

Why J space gates: an additive platform offset per gene is a constant within Oetjen and cancels in a within-Oetjen slope or
contrast, and numerator and denominator then share the per-gene weights 1/s_g. A ratio of an Oetjen-SD quantity to a
GSE185381-SD quantity (the first draft) is not unit-consistent (the same imprint z-scored within different sets has score SD
0.45–0.80) and is biased toward passing. J space amplifies Oetjen's noise on genes with small s_g; that widens the CI without
biasing the point estimate and is what the SD-ratio check in A3 and the random-set ranks measure.

## 3. Tests

### 3.1 A1 — age gradient within Oetjen (decisive)

In J-Oetjen, on the eligible donors (K0):
- `lm(score ~ age + sex)`, age in years. **f_slope = β_age × 31.5 / (D_res_J / r)**.
- **f_contrast = (mean score of donors ≥ 60 − mean score of donors ≤ 53) / (D_res_J / r)**.
- **f_age = max(f_slope, f_contrast)** (conservative for the claim). **The point estimate decides: f_age < 1/3 passes.** The 1/3
  line is GATE_A1 (4)'s "less than one third of the AML–healthy difference".
- Donor bootstrap: 1,000 draws with replacement over the eligible donors' (score, age, sex) rows; U_J, Oetjen's TMM factors,
  m_g, s_g, D_res_J and r fixed; draws with fewer than 2 old or 2 young donors are redrawn; f_age recomputed as the max within
  each draw. 95 % percentile CI, **reported** as a donor-resampling CI conditional on the yardstick (relative SE about 10 %);
  a pass whose upper bound exceeds 1/3 is worded as a bound (section 5).
- Leave-one-donor-out (applies to a pass only): the fit and the group means are recomputed without each eligible donor in turn,
  U_J, Oetjen's TMM factors, m_g, s_g, D_res_J and r fixed; **if dropping any single donor moves a passing f_age to ≥ 1/3 the
  result is "not decisive"** and goes to discussion (K 84 y sits 17 years beyond the next oldest donor and has leverage ≈ 0.3 in
  the slope fit).
- Negative f_age (older donors score *less* AML-like) is read as "age opposes the direction", reported with its CI; it passes the
  gate like any value below 1/3.
- Random sets: the 1,000 matched sets scored in the same space (each set's slope and contrast divided by the *imprint's*
  D_res_J / r, since a random set's own gap is near 0) give the P5/P95 of f_slope and f_contrast and the imprint's rank (a
  *relative* statement, since ageing shifts many genes' expression). Not gating; a passing f_slope above the random P95 adds
  one sentence to the wording (section 5).

Reported alongside (never gating, never in claim wording): the 402 set; the unadjusted fit (`score ~ age`); `score ~ age + sex
+ naive_share` (share of CD4 Naive + CD8 Naive among the donor's kept T cells) to separate age-via-composition from age-via-state;
W-Oetjen ρ(score, age) with its CI and its rank among the random sets; the same ρ on the CD4 Naive and CD8 EM (EM1 + EM2) subtype
pseudobulks (W space; donors with ≥ 20 kept cells of the subtype); partial Spearman of the score with the two ambient shares.

Power (fixed now; assumes a within-age residual SD of Oetjen scores of about 1 healthy-SD and uses the observed age SD of 15.0 y):
the 1/3 line (1.3 SD over 31.5 y) corresponds to ρ ≈ 0.5 (one-sided Spearman power at n = 20 ≈ 0.75); full H_age (3.9 SD over
31.5 y) to ρ ≈ 0.9 (power > 0.95). The gate is on r × f_age, not on the ρ test. SE(f_slope) ≈ 0.13 r and SE(f_contrast) ≈ 0.14 r
(5 old vs 10 young donors), and taking the max adds about +0.04 r under H_AML. The false-fail rate under H_AML is about 1 % at
r = 1, 5 % at r = 1.3 and 9 % at r = 1.6; an age effect that truly explains 40 % (50 %) of the cross-fitted gap passes about 25 %
(7 %) of the time at r = 1 and about 30 % (16 %) at r = 1.6. r is computed and reported before any Oetjen score is read, and these
rates are restated for the observed r in `ws_s9_yardstick.csv`. This is why the CI is reported and why a pass is worded as a bound.

### 3.2 A2 — age-sensitive genes and gene-level age–AML alignment

- **A2a (reported only).** Age-sensitive genes = genes with raw p < 0.05 for age in edgeR QL (`glmQLFit`, robust = TRUE) on the
  eligible Oetjen donor T pseudobulks, design `~ age + sex`, genes with CPM > 1 in at least half of the donors after the same
  blacklist removal as 2.4 (the FDR < 0.10 set is reported too). Remove them from the 370; rescore GSE185381 in the S5 space; ratio = gap_robust / gap_in. Reference:
  1,000 random removals of the same number of genes from the 370. Also reported: the count removed and, for the removed genes,
  the sign concordance between their age direction (Oetjen) and their AML direction (GSE185381). Stated limitation: at n = 20
  this screen has little power (a gene whose age effect is a third of its AML effect has t ≈ 1.4), so the ratio can only fall
  when A1 has already failed; it is the user's requested analysis and is kept as a report.
- **A2b (wording modifier; role to be confirmed at approval, section 8 D6).** Over the primary genes in U_J (all carry a
  log2FC_AML): **the OLS slope b of each gene's Oetjen age effect, log2FC_age × 31.5 years (the A2a model's age coefficient),
  on its GSE185381 AML log2FC, divided by the gene set's reliability λ = 1 − mean(SE²_AML) / var(log2FC_AML)** (SE_AML per gene
  = |log2FC_AML| / √F from the S5 model refit in 33 with 27's settings). b/λ is the gene-level analogue of f_age: under a shared
  per-gene age share a the raw OLS slope estimates λ·a (regression dilution by the noise in log2FC_AML), and λ differs between
  sets (about 0.90 for the imprint, whose log2FC are large and precisely estimated, against about 0.47 for expression-matched
  random genes), so without the correction the imprint would exceed the random P95 under a *shared* rate (simulated 44 % of
  runs at a = 0.15, 63 % at 0.20). Spearman ρ and sign concordance are reported alongside; b/λ read as an age share is reported,
  never gating. Null: b/λ for each of the 1,000 matched random sets with its *own* log2FC_AML and its own λ, drawn with
  `matched()` from U_J ∩ tested(`ws_s5_de_T.csv`); the genome-wide median is reported as a fact in its own right (if the
  AML-vs-healthy contrast of GSE185381 is partly age, genes genome-wide show it). **Reading: b_imprint/λ below the random-set P95
  → "the imprint genes' age effects are not more aligned with their AML effects than the 95th percentile of matched random
  sets"; at or above → the pass wording must add "the imprint genes show an age alignment in healthy marrow beyond the
  genome-wide rate".** A slope and not a correlation because the imprint's log2FC spread is 2.5× that of matched random genes
  (2.2 (iv)): under a shared alignment rate any correlation-type statistic of the imprint exceeds the random P95 (simulated FAIL
  rates 41–99 %), whereas the disattenuated slope estimates the rate itself. It detects an imprint-specific per-gene age share
  of about 10 % of the AML effect, far below the 1/3 line, which is why it modifies wording rather than gates (D6).

### 3.3 A3 — direction against each external healthy cohort (gates only the sentence it adds)

In J-E for E ∈ {Oetjen, Griffin, R-S}: D_res_J; **S_E** = mean(E young) − mean(GSE185381 healthy) (the study effect: same age
range and biology, different study and platform); **A_E** = mean(E old) − mean(E young) (the age effect within E; for Oetjen
A_E / (D_res_J / r) is f_contrast); D_dir = mean(AML) − mean(E old) = D_res_J − S_E − A_E (reported).
- **Interpretable** if |S_E| < D_res_J / 2 **and** SD(E young) / SD(GSE185381 healthy) ∈ [0.5, 2] (both undeflated: they compare
  the study offset and spread with the *observed* gap, not with an effect size). S_E estimated from ≤ 3 young donors (Griffin 3,
  R-S 2) is **descriptive**; the sentence for such a cohort is given with its n and never enters the gate. In practice only
  Oetjen can be interpretable.
- **Direction holds against E** if A_E / (D_res_J / r) < 1/3 **and** at most one of E's old donors has a study-corrected score
  (score − S_E) above the midpoint mean(GSE185381 healthy) + D_res_J / 2 (the observed classification boundary, undeflated).
  The position of every old and young donor relative to that midpoint is reported for all three cohorts. (Under H_AML with the
  old donors' corrected scores spread like the GSE185381 healthy scores, "every donor below" would fail about 13 % of the time
  at SD ratio 1 and 59 % at ratio 2; "at most one above" fails about 1 %, 8 % and 20 % at ratios 1, 1.5 and 2.)
- If Oetjen fails this after GATE_S9 passes, the pass wording adds "k of n older Oetjen donors score above the AML–healthy
  midpoint after study correction" and the contradiction is recorded for discussion.

### 3.4 Informational analyses inside GSE185381 (no new data)

As S5 §3.4(a) promised: sign concordance of the 370 genes' log2FC between the 4 AML patients < 60 (32, 55, 57, 58 y) and the
5 healthy donors aged 39–53 (39, 42, 47, 50, 53; S5 §3.4(a) wrote 4 from the miscount recorded in 2.1) (no test). Extending it: Spearman ρ(score, age) within the 15 eligible AML and within the 10 healthy
donors (S5 space, primary set), and the position of the four AML patients < 60 relative to the healthy range. Confounded by
disease heterogeneity; reported only.

### 3.5 Cross-fitted yardstick (user decision D1, 2026-10-23: adopted)

The in-sample gap is inflated by gene selection (winner's curse), so a fixed 1/3 of it has no fixed meaning. Leave-one-sample-out
over the 25 GSE185381 T pseudobulks: for each left-out sample, refit the S5 differential model on the other 24 (edgeR QL robust,
design `~ status + hb + gran + sex + subtype_share`, blacklists removed before filtering, `filterByExpr` + CPM > 1 in half, as in
`27_s5_imprint.R`), select the fold's imprint (FDR < 0.10, |log2FC| > 0.5, up and down; the 32 ambient-prone genes removed),
take the fold's per-gene mean/SD from the 24, and score the left-out sample (log-CPM computed once on the 25, label-free: TMM and
the CPM filter do not use the status label; the in-fold filterByExpr and gene selection do). A fold whose up or down half selects
no gene is recorded and that half-score is set to 0. gap_oof = mean(out-of-fold AML scores) − mean(out-of-fold healthy scores).
**r = gap_in / gap_oof, floored at 1 (never a loosened yardstick); it requires gap_oof > 0. If gap_oof ≤ 0 or r > 2 the yardstick
does not cross-validate or is under-resolved (false-fail under H_AML ≥ 16 % at r = 2): GATE_S9 is not evaluated and the line
stops for discussion (P5).** Reported in `ws_s9_yardstick.csv`: gap_oof, gap_oof / SD(out-of-fold healthy scores) beside 3.87,
r, D_res_J beside gap_in, and the 3.1 error rates restated for the observed r. **Every 1/3 comparison (3.1, the A_E rule of 3.3,
K2) uses the deflated denominator D_res_J / r.**

## 4. Known answers (run before the age relation A1–A3 is read on a new cohort; one statistic, one threshold, one consequence each)

| id | statistic | threshold | consequence if it fails |
|---|---|---|---|
| **K0** eligibility | eligible Oetjen donors (≥ 20 kept T cells after pooling) | ≥ 15, including ≥ 4 of the 5 donors ≥ 60 | no decisive cohort; GATE_S9 is not evaluated; stop for discussion (P5). No fallback rule. |
| **K1** age positive control | per Oetjen donor, share of `CD8 Naive` among **all projected T cells** (T-group fine labels, `high_error` dropped, concordance not required); Spearman with age | ρ < 0, one-sided p < 0.05 | S9 cannot see a known age effect: uninformative, stop. `CD4 Naive` share reported. In Griffin and R-S: descriptive (both old donors below every young donor); if not, that cohort's A_E is labelled "age effect not visible in labels" and not read as an age estimate. Griffin author labels (CD8Naive, CD4Naive) vs ours: per-cell agreement reported, no statement made from it. |
| K1b gene-level age control (reported; D3) | mean age coefficient over CCR7, SELL, LEF1, TCF7, IL7R; over GZMH, GZMB, PRF1, NKG7, KLRG1, B3GAT1 (A2a model) | first < 0, second > 0 | informational for the K1 discussion only |
| **K2** resolution | technical replicates scored as libraries (2.4, K2 row; no age is read): \|score(C1) − score(C2)\|, \|score(S1) − score(S2)\| | both < (D_res_J / r) / 3 | measurement noise reaches the effect to be resolved: stop before A1 is read. Second aspirates (Ck vs C1/C2; Sk1, Sk2 vs S1/S2) reported. |
| **K3** identity | per library, XIST vs RPS4Y1 UMI sums over kept cells (S4 rule: female if XIST > 3·max(RPS4Y1, 1), male if RPS4Y1 > 3·max(XIST, 1); neither → "unclear", which counts as a disagreement for its donor) vs the donor table; 29 donors | all agree (HD-36 male → YHD658, HD-95 female → YHD636 asserted in the ingest as well) | a disagreeing donor's identity is unresolved → excluded before any score, recorded; two or more → stop |
| **K4** coverage | \|primary set ∩ U_J\| / 370 per cohort | ≥ 0.80 (S8's coverage check) | Oetjen below → no decisive cohort, stop; Griffin or R-S below → that cohort descriptive only |
| **K5** planted H_age (machinery; the real age relation is never read) | donor ages permuted once (SEED; redrawn while \|Spearman ρ(age_perm, age)\| > 0.3, so the real association is destroyed). On the eligible donors' T pseudobulk **counts**, each primary gene's counts are multiplied by 2^(d_g × (age_perm − 32.5) / 31.5) and rounded, d_g = that gene's GSE185381 AML-minus-healthy mean log-CPM difference in the J space of 2.4 (log2 units; the plant therefore adds D_res_J per 31.5 years); J log-CPM and the A2a model recomputed on the planted counts with the permuted ages, **U_J, m_g, s_g and D_res_J fixed at their unplanted values**. Recovery, undeflated: f_slope(planted) − f_slope(unplanted), f_contrast(planted) − f_contrast(unplanted), both under the permuted ages; second plant with d_g / 3; planted A2b. Expected: slope recovery about 0.98 (the prior count and integer rounding attenuate low-count genes), contrast about 0.885; TMM re-estimation is immaterial (the plant changes a library by < 2 %) | full plant: slope recovery ∈ [0.95, 1.03], contrast recovery within ± 0.05 of (mean age of eligible old − mean age of eligible young) / 31.5 (0.902 if all 20 are eligible); third plant: ∈ [0.30, 0.37] and within ± 0.03 of 0.301; planted b/λ above the random-set P95 | code defect: corrected and recorded as a dated deviation before the real age relation is read (S8 §7 precedent). 100 sign-scrambled imprints: P5/P95 of f_slope reported (calibration report, not a test). |
| **K6** regression | S9 code vs S8: `logz()` on the 15 GSE185381 AML T pseudobulks alone (z over the 15, as 32's disc set) vs `ws_s8_scores_disc.csv` score_primary; the yardstick on the 25 | Spearman > 0.99 and max \|diff\| < 0.05; \|gap − 1.970\| < 0.01 | as K5 |
| **K7** metadata | per library, donor-table age vs `meta_<GSE>_v1.0.csv` age_years and sex; per donor, the first-listed library's age and sex vs the BoneMarrowMap reference ExactAge/DonorSex (Oetjen); kept donors only | all agree, with Ck 59 vs C 60 and Sk1/Sk2 55 vs S 56 the only tolerated per-library differences | the table is corrected and the correction recorded before any score |

## 5. Gate and readings

**Evaluation order:** K7, K3, K0, K4 (no score involved; K3 before K0 because an identity exclusion changes the counts K0
tests); r (3.5; GSE185381 only); K6; K2 (library scores, no age); K5 (planted, permuted ages); K1 (labels, no score). Then the
f_age point estimate (≥ 1/3 is a fail whatever the leave-one-donor-out shows); then A2b (wording); then leave-one-donor-out (a
pass only; a flip ends the evaluation: no claim wording is issued, the CI and A3 are reported in the tables only); then the CI
wording; then A3.

**GATE_S9 (age bound on the direction) = K0, K1, K2 and K4(Oetjen) pass AND f_age < 1/3.** A2b and the J-space random rank
modify the wording; A3 gates only the sentence it adds. Reported quantities (402 set, subtype scores, W-space ρ, partial
correlations, A2a ratio, b/λ read as an age share — the wording quotes it only beside its P95 — and D_dir) never enter claim
wording.

| outcome | reading (wording fixed here; bracketed values filled from the S9 tables) |
|---|---|
| gate passes | "In one healthy marrow cohort (Oetjen 2018: [n_elig] of 20 donors, [min]–[max] y, 10x 3' v2, one protocol, part of the BoneMarrowMap reference), the T imprint's age shift over the GSE185381 age gap is [f_age] [95 % CI] of the cross-fitted AML–healthy gap, below one third. The age alternative to the direction is **bounded, not excluded**. Direction versus healthy remains shown in the defining cohort only (no independent AML + healthy cohort exists in this work); the 32 ambient-prone genes stay excluded from every statement. The S5 label 'AML-or-age' is narrowed to 'AML, with an age share bounded below one third in healthy marrow (one protocol)'." Plus the A2b clause: below P95 → "the imprint genes' age effects in healthy marrow are not more aligned with their AML effects than the 95th percentile of matched random sets (b = [b], P95 = [P95])"; at or above → "the imprint genes show an age alignment in healthy marrow beyond the genome-wide rate (b = [b], P95 = [P95])". Per interpretable A3 cohort: "older healthy donors in [E] do not score like AML (A_E / (D_res/r) = [..], n = [..])"; descriptive cohorts are given with their n. |
| pass, CI upper bound > 1/3 | the same wording as a bound: "age share below one third at the point estimate; upper bound [X]" |
| pass, f_slope above the J-space random-set P95 | the wording adds "a weak age association exists (f_slope = [..], random P95 = [..])" |
| pass, Oetjen fails the A3 at-most-one rule | the wording adds "[k] of [n] older Oetjen donors score above the AML–healthy midpoint after study correction"; contradiction recorded for discussion |
| pass, Oetjen not interpretable in A3 | the wording adds "Oetjen's study offset or spread is outside the interpretable range (S_E = [..], SD ratio = [..])"; its at-most-one rule is reported, not worded |
| f_age ≥ 1/3 | **GATE_S9 FAIL**: the claim stays "AML-or-age" with the age share quantified; stop for discussion (P5). Stage 2 (more older donors: Bergiers GSE232988 needs per-donor ages; Wu GSE196052) is a separate decision, not a continuation. |
| LOO-donor flip (a pass only) | "not decisive"; discussion |
| K0, K1, K2 or K4(Oetjen) fails; gap_oof ≤ 0 or r > 2 | as in sections 4 and 3.5; no age relation is read on the new cohorts |
| K3, K5, K6, K7 fail | as in section 4 (identity exclusion, or a stop when two or more donors disagree in K3, or a dated code/table correction before the real run) |

No threshold, gene set, donor list or wording is changed after a number is seen; later changes are appended, dated sections.

## 6. Outputs, code, order of operations, configuration edits

Outputs: `results/tables/ws/ws_s9_scores_<cohort>.csv` (donor × gene set: J and W scores, library scores for K2, cell counts,
ambient shares), `ws_s9_tests.csv` (test × gene set × statistic × threshold × pass), `ws_s9_age_genes.csv` (A2a/A2b per gene),
`ws_s9_direction.csv` (A3 per cohort and donor positions), `ws_s9_yardstick.csv` (3.5), `ws_s9_known_answers.csv`,
`ws_s9_checks.csv`; GATES.md block "A1-S9". The age-robust set is written as `ws_s9_robust_imprint.csv` for the record only
(section 7).
Code after the tag: `scripts/ws/33_s9_age_direction.R` (+ sbatch), about 300 lines; edgeR, limma, data.table, Matrix,
SeuratObject as in 26/27/32; **no new dependency**.

Order of operations (each step gets its own plan, approval and commit; no age relation is read on a new cohort before the tag):
1. `meta_GSE120221_v1.0.csv` (25 rows), `meta_GSE227690_v1.0.csv` (13 + BM6.1), `meta_GSE276360_v1.0.csv` (4 + 3 CCUS), 72
   columns (meta_v2.3), template = the GSE185381 Control rows (disease healthy_donor, control_type healthy_donor, timepoint_class
   Healthy, age_years, age_group, sex) — not the E-MTAB-11536 rows (disease non_heme_condition would make every donor
   "other_disease" in `26_s4_pseudobulk.R`). Platform 10x_3p_v2 / 10x_3p_v3 or 10x_3p_v3.1 per library / unknown (Chip A, kit
   not stated); sorting none / viability_only / viability_only; cell_prep mononuclear_ficoll / mononuclear_unspecified ×2;
   sample_state cryopreserved / unknown ×2; raw_data_type fastq; matrix_level filtered_umi / processed_object ×2.
2. `config_qc.R`: ROLE_TABLE rows (study_role exactly `Healthy-control`; is_longitudinal FALSE; discovery_candidate FALSE;
   l2_capable **FALSE** (D2) with l2_reason "S9 age control; not a CCC-cohort member"); PLATFORM_TABLE rows (10x|3prime,
   10x|3prime, 10x|unknown); `TP_HEALTHY_DATASETS` += the three; `INGEST_EXCLUDE` rows GSE276360 `^CCUS-`, GSE227690 `^BM6`;
   `DATASETS_UPSTREAM_FILTERED` += GSE227690, GSE276360. Dataset counts 13 → 16 in `submit_J1_ingest.sbatch`,
   `submit_J2_qc_report.sbatch` (the `03_qc_report__ALL.csv` combiner) and `submit_03_per_sample_qc.sbatch`. DISPOSITION rows for
   BM6.1 and CCUS-50/53/PR in `99_admin/12_cohort_roster_check.R` (D4). `config_hierarchy.R` HEALTHY_DATASETS is used only by
   archived scripts and is left alone.
3. `ingest_GSE120221.R` (explicit `ReadMtx` paths for the `<GSM>_{barcodes,genes,matrix}_<lib>` layout; Ensembl remap as in 2.3;
   asserts: ≥ 0.95 of the 370 present, no primary gene among duplicates or collisions, max percent.mt > 0 per library),
   `ingest_GSE227690.R` (dense tsv → one object per replicate library; `ingest_keep` drops BM6; author CellType kept as
   `anno_CellType`; replicate → GSM join asserted complete), `ingest_GSE276360.R` (`readRDS` of the on-disk `.rds`; `ingest_keep`
   drops CCUS; HD-36 → YHD658 and HD-95 → YHD636 asserted by XIST/RPS4Y1; `stopifnot` every kept prefix matches `^HD-`). All
   three: Sample = GSM-prefixed library, `Patient_ID` = donor, `GSM`/`Library` set, `extra_meta` with age_years, age_group, sex
   (→ `anno_*`, read by ws/26), timepoint "Baseline", disease "Healthy". Run by hand, not through the J1 array.
4. `00_curated_manifest.R` and `01_dataset_roles.R` (not `02_study_split.R`). The manifest build must not run before step 1: a
   metadata directory without `meta_*.csv` breaks it; each directory holds exactly one `meta_*.csv` (the donor tables keep their name).
5. `03_per_sample_qc.R --dataset` × 3; rebuild `03_qc_report__ALL.csv` **before** any `03_hierarchy` or `99_admin` call
   (`qc_rds_roster(on_extra = "error")` stops otherwise; `99_admin/14` is not run).
6. `01_bmm_project.R --datasets GSE120221,GSE227690,GSE276360`; `05_marker_annotate.R` array over the new samples only; `06` and
   `07` `--dataset` × 3.
7. `ws/22_decontx_library.R --index` and per unit × 3; `ws/26_s4_pseudobulk.R --dataset` × 3.
8. Tag `prereg-s9-v1` (set at approval, before step 3 at the latest) → `33_s9_age_direction.R`, known answers first, in the
   order of section 5.

## 7. Not done, by decision

The age-robust set (A2a) never replaces the 370 primary set anywhere. No imprint is re-defined with age in the design (that would
be a new pre-registration). The three cohorts are never pooled into one age gradient. Only the T compartment is scored. The
paediatric cohort, GSE227903 remission samples, CCUS, BM6 and the NK/T-sorted libraries are not touched. No AML sample from any
new cohort exists or is sought; "direction versus healthy" is not replicated here. `02_study_split.R` is not rerun. `config_ws.R`
sections 0–14 are not edited. No file is deleted; no new dependency is installed.

## 8. Decisions taken at approval (user, 2026-10-23)

D1 cross-fitted yardstick deflation r applied to the gate: yes. D2 l2_capable FALSE: yes. D3 K1b reported: yes. D4 curated rows
plus DISPOSITION rows for the excluded libraries: yes. D5 keep "S9" here, renumber the remission pre-registration S10: yes.
**D6 (decided at approval, 2026-10-23): A2b is a wording modifier, not a gate**; the gate is the 1/3 magnitude line alone.
The pass label of section 5 ("AML, with an age share bounded below one third in healthy marrow (one protocol)", stated as a
bound, not an exclusion) was confirmed at the same time. The paragraph below records the alternative that was not taken.
This draft makes it a wording modifier (3.2): the gate is the 1/3 magnitude line alone. The alternative is to make "b_imprint/λ below the random P95" a second gate condition; that would fail
the direction claim for an imprint whose genes carry an age share of about 10 % of their AML effect, far inside the 1/3 the
design tolerates, and it rests on the reliability correction λ (3.2): without it the imprint exceeds the random P95 under a
shared genome-wide alignment of 15 % in about 44 % of runs. The first-round plan had A2b gating with a Spearman statistic and a
"5 % false-fail by construction" claim that the second review showed to be wrong; the statistic is corrected either way.

## 9. Record of the pre-approval reviews (2026-10-23)

First round (three independent adversarial reads of the plan: statistical design; data and pipeline feasibility;
pre-registration discipline; 42 findings). Changes they forced: (1) the gating age effect and the yardstick are measured in one
space (2.4), because the first draft divided an Oetjen-SD quantity by a GSE185381-SD quantity, which is not unit-consistent and
biased toward passing; (2) A2's removal ratio is demoted to a report, because at n = 20 it can only fail when A1 has already
failed; (3) A3's rule is study-corrected and symmetric in the sign of the platform offset; (4) Sample = library at QC, donor
pooling in the S9 code, because a donor-level Sample triples or quadruples the doublet prior of donors C and S; (5) Δage = 31.5 y
from the eligible samples, not 35; (6) f_age takes the larger of slope and contrast, with a leave-one-donor-out rule, because a
step-type age effect at 60 is under-estimated by a slope; (7) every known answer carries a consequence and every outcome pattern a
reading; (8) the sign-scramble "known answer" (which cannot fail) is replaced by a planted H_age; (9) the pass wording is a bound,
not an exclusion, and keeps the 32-gene acknowledgement and the "direction versus healthy only in the defining cohort" limitation;
(10) the Oetjen-in-BoneMarrowMap fact and the Ensembl remap are registered as inputs; (11) the cross-fitted yardstick (3.5)
answers the winner's-curse objection to the 1/3 line.

Second round (two reads of the written draft: internal consistency; statistics of the new elements; 26 findings). Changes:
(12) A2b's statistic is the slope of log2FC_age on log2FC_AML, not a Spearman correlation, because the imprint's log2FC spread is
2.5× that of matched random genes and any correlation statistic then fails under a shared genome-wide alignment; its role is put
to the user (D6); (13) K5 plants on counts under permuted ages, so the real age relation is never read before the known answers
pass, and its windows are tight enough to catch a 5 % scale error; (14) A3 uses an at-most-one-donor rule (the every-donor rule
fails 13–59 % of the time under H_AML); (15) r requires gap_oof > 0 and is bounded at 2, with a stop otherwise; (16) an evaluation
order and the "LOO applies to a pass only" rule remove the precedence ambiguity between outcome rows; (17) the W-space ρ no longer
adds wording (a J-space quantity does); (18) the pass wording takes its donor count and age range from the tables; (19) K2, K3, K6
and K7 are fully specified (library-level DGEList; "unclear" counts as disagreement; z over the 15; per-library vs per-donor
comparison with the two tolerated differences); (20) S5 §3.4(a)'s promised concordance analysis is included and S5's 6/4 healthy
age split is recorded as a miscount of the same ten donors.

Third round (one read of the rewritten draft against the 26 second-round findings and the new elements; all 26 confirmed
resolved; 11 new findings, one major). Changes: (21) A2b's slope is divided by each gene set's reliability λ, because regression
dilution differs between the imprint (λ ≈ 0.90) and matched random genes (λ ≈ 0.47) and the uncorrected comparison would flag
"beyond the genome-wide rate" under a shared rate; (22) the J-E recipe states the blacklist removal that `logz()` performs, and
A2a's filter follows it; (23) K5 fixes U_J, m_g, s_g and D_res_J at unplanted values, redraws permutations correlated with the
true ages, and states the expected recoveries (≈ 0.98 slope, ≈ 0.885 contrast) with a window of [0.95, 1.03]; (24) K3 precedes
K0 in the evaluation order; a leave-one-donor-out flip ends the evaluation; a non-interpretable Oetjen adds its own sentence;
(25) the random-set f_slope/f_contrast use the imprint's denominator; (26) 2.4 carries one row per test (3.4 and K6 added) and
3.4 names five, not four, healthy donors aged 39–53.

## 10. Execution notes and approved amendments (appended 2026-10-23, before any score on the new cohorts)

10.1 **Ingest (section 6 step 3), commit 4cf755e.** Run as three named SLURM jobs (3658414 GSE276360, 3658415 GSE120221,
3658418 GSE227690), not through the J1 array, so that no array index could re-ingest an existing dataset; same scripts and
environment otherwise.
- GSE120221: 25 libraries, 90,653 cells. Ensembl remap as in 2.3: 32,351 genes mapped (mapped genes take 2020-A's
  make.unique'd names, so duplicate symbols point to the same Ensembl ids as in GSE185381), 80 Cell Ranger 2.0 names suffixed
  for collisions; table `results/tables/00_ingest/GSE120221_gene_remap.csv`. Primary genes on the axis: 359/370, as
  registered; POLR2J3 is present under its Cell Ranger 2.0 name (its old Ensembl id ENSG00000168255 is not in 2020-A). Sex
  genes agree with the donor table in 25/25 libraries.
- GSE227690: two deposit defects found at ingest. (i) The `BM_counts.tsv.gz` line for CCDC60 carries 23,679 values for 20,411
  cells; it is dropped (not an imprint gene). The first run (job 3658416) read only the 2,449 genes before that line because
  fread stops there with a warning; the ingest's marker check caught it after the outputs were written, and the rerun
  overwrote them; the reader now keeps well-formed lines only and checks gene count and per-cell totals before writing (median
  1 UMI per cell below the authors' nCount). (ii) The deposit has no MT- genes, XIST, MALAT1 or NEAT1: percent.mt is 0 for
  every cell (the authors' own cell QC applies) and sex can only be read from Y-linked genes (10.2). 13 libraries, 17,937 cells.
- GSE276360: 4 libraries, 10,213 cells; HD-36 -> YHD658 (male) and HD-95 -> YHD636 (female) confirmed by XIST/RPS4Y1 at ingest.

10.2 **K3 amendment for a deposit without XIST (approved by the user, 2026-10-23).** For GSE227690 only, a library's sex is
called from the summed UMI of the Y-linked genes present (RPS4Y1, DDX3Y, KDM5D, UTY, EIF1AY, ZFY, USP9Y) per kept cell:
male if >= 1.0, female if <= 0.2, otherwise "unclear", which counts as a disagreement as in K3. GSE120221 and GSE276360 keep the
S4 rule. The rule was written while the GSE227690 rerun was in progress and before its Y-gene numbers were read; the
thresholds come from the other two cohorts (GSE276360: RPS4Y1 alone 3.5 per cell in the male donor, 0.013-0.063 in the three
female donors; GSE120221: 1.5-3.2 in male libraries, <= 0.005 in female libraries). Seen afterwards, at ingest, on all cells:
male GSE227690 libraries 6.2-10.3 Y-UMI per cell, the female donor's (BM2) <= 0.0006. Without the amendment K3 would have
counted all five GSE227690 donors as disagreements and stopped S9 on a property of the deposit.

10.3 **Step 4 of section 6 amended (approved by the user, 2026-10-23): `01_dataset_roles.R` is not run.** No script on the S9
path reads its output (`01_sample_role_manifest.csv`; its consumers are the CCC line, `02_study_split.R` and the malignancy
stage), and a full run re-reads all 16 ingest objects and rewrites a table the CCC line uses. `00_curated_manifest.R` is run;
its rows for the 13 existing datasets must equal the committed version, otherwise the line stops for discussion.

10.4 **Code review and smoke runs before the real run (2026-10-26; appended before the real run; nothing in sections 1-9 changed).**

`scripts/ws/33_s9_age_direction.R` implements sections 2-5 and 10. Before the real run it was smoke-tested five times with
`--smoke` (donor ages of the three new cohorts permuted; 50 random sets, 100 bootstrap draws, 5 of the 25 cross-fit folds;
outputs suffixed `_smoke`; GATES.md not written) and reviewed by three independent agents (statistics/definitions, data flow and
R correctness, discipline/outputs) against this document. Smoke history: job 3669045 stopped at K1 (expected: under permuted
ages the CD8 Naive share has no age relation; smoke now reports phase-1 failures and continues), job 3669047 crashed in the
subtype pseudobulk (26's `$fine` matrices hold only the libraries where the subtype occurs; aligned on all libraries), job
3669048 completed end to end, jobs 3669078 and 3669088 completed after the review fixes below. Known answers in the smoke runs: K7 42/42
and 20/20, K3 29/29, K0 20 eligible incl. 5 old, K4 358/370, 329/370, 355/370, K6 Spearman 1.0000 / max diff 0.00000 and
|gap_in - 1.970| = 0.0014, K2 0.053 / 0.033 against 0.356, K5 recoveries 0.975, 0.880 vs 0.902, 0.322, 0.294 vs 0.301 and
planted A2b 1.006 vs P95 0.064 (5-fold r = 1.84; the 25-fold r of the real run is read there).

Disclosure on the yardstick: the registered 1.970 (section 2.1) was computed on 10-23 without the blacklist step of S8's
scoring; with the exact S8 `logz` (blacklist removed before the CPM filter) the primary-set gap on the 25 is 1.9714 =
3.84 healthy-SD (healthy SD 0.513), not 3.87. K6's tolerance (0.01) covers it; the number the real run uses is 1.9714.

Review findings fixed before the real run (none changes a threshold, a gene set, a donor list or a registered sentence):
1. (major) The wording clause `f_age [CI]; r` was built with an `if` without `else` inside `sprintf`; when the CI upper bound
   is <= 1/3 the clause would have vanished from the log and GATES.md. Fixed (`else ""`).
2. (major) In smoke mode the Oetjen scores table wrote the 25 library rows with the REAL library ages and sexes next to the
   library-level J scores, so the real age-score relation of the decisive cohort was on disk after job 3669048. Nobody opened
   those rows (the session only printed row counts and column names); the file was overwritten by job 3669078, which writes NA
   there under `--smoke`. Recorded here because the point of the permutation is that no one can read the relation early.
3. K7 against the BoneMarrowMap reference passed vacuously on an empty merge; now every donor must be found. K7 per library is
   evaluated on the kept donors' libraries only (BM6 dropped), as section 4 says.
4. NA results of gating or defect known answers were printed as "reported" and let phase 2 proceed; every gating pass is now
   `isTRUE(...)`, so NA counts as a failure.
5. The K3 consequence for exactly one disagreeing donor (excluded before any score, recorded; K0/K4 on the reduced set; its
   libraries dropped from K2) was not implemented; only the two-or-more stop was. Implemented; K3 is 29/29 in every smoke run.
6. A2b's imprint set is now primary ∩ U_J at both call sites (planted K5 and the real A2b), the same universe as the random
   sets (section 3.2); the smoke value was already computed on 358 genes.
7. Reported quantities aligned with their registered definitions: the 402-set sensitivity uses its own D_res_J / r (section
   2.4, same space and genes); the naive-share covariate is the share among the donor's KEPT T cells (section 3.1); the 100
   sign-scrambled imprints are reported over f_slope's denominator D_res_J / r.
8. Wording: the donor count and age range in the pass sentence and the GATES block come from the tables; one sentence per A3
   cohort (interpretable: "older healthy donors in [E] do not score like AML"; descriptive: given with n young / n old; not
   interpretable or K4-demoted: stated as such); no claim clause is issued on a LOO flip or a FAIL; the A3 rows of descriptive
   cohorts carry pass = NA in ws_s9_tests.csv (direction_holds stays in ws_s9_direction.csv). Pipes inside GATES table cells
   are escaped.
9. Subtype pseudobulk: donors restricted to the eligible set; empty donor set guarded.
10. (third review, major) `ws_s9_scores_<cohort>.csv` lacked the ambient shares and the doublet score that section 6 lists; each
    donor row now carries the cell-weighted hb_share_nonery, gran_share_lymph and dbl_mean (NA: T3a has not run on these
    cohorts) from 26's library rows, and the ambient partial correlation uses the same per-donor table.
11. (third review) Yardstick disclosure in the outputs: `ws_s9_yardstick.csv` carries the registered 1.970 / 3.87 beside the
    recomputed values with a note, and the GATES yardstick line prints both. Deviation recorded: r (section 3.5) and the A2a
    ratio use the recomputed 1.9714, not the registered 1.970, so that gap_in and gap_oof share the blacklisted gene universe;
    the effect on r and f_age is 0.07 %. Every "3.87 healthy-SD" in sections 1, 2.1 and 3.5 reads 3.84 under the exact recipe;
    the section-3.1 error rates are in r units and unaffected.
12. (third review) The f_age / CI / r clause is reported under every outcome (section 3.1 "reported"), including a FAIL, which
    now quotes f_age ("age share quantified"); the pass verdict is the registered sentence with the values filled, in the
    "below one third at the point estimate; upper bound [X]" form when the CI upper bound is above 1/3.
13. (third review) K1 drops a K3-excluded donor; the per-library K7 requires the full join (42 libraries); the W-space rho
    carries a Fisher-z CI; folds with an empty up or down half are counted in the yardstick table and the 3.5 detail; the
    6-donor minimum of the subtype rho is stated in its detail; the positions file is listed in the script header.
The real run is `sbatch scripts/ws/submit_ws_s9_age_direction.sbatch` on the file with md5 7ceeeb721912d15dcc3261e8f6dce6f8 (commit
c222ebc).

## 11. Outcome (appended 2026-10-26 after the real run; sections 1-9 unchanged, tag prereg-s9-v1)

Run: `sbatch scripts/ws/submit_ws_s9_age_direction.sbatch`, SLURM job 3669091, 2026-10-26 11:15-11:20 JST, script md5
7ceeeb721912d15dcc3261e8f6dce6f8 (commit c222ebc). Outputs: `results/tables/ws/ws_s9_{known_answers,checks,yardstick,tests,
direction,direction_positions,age_genes,robust_imprint}.csv`, `ws_s9_scores_{GSE120221,GSE227690,GSE276360}.csv`,
`/LARGE1/.../ws_pseudobulk/ws_s9_null_draws.rds`, GATES.md block A1-S9. 8 checks and 23 known-answer rows, none failed. The
age relation of the new cohorts was read once, in this run.

**11.1 Phase 1 (known answers, section-5 order; all PASS).** K7 42/42 libraries, 20/20 Oetjen donors vs the BoneMarrowMap
reference. K3 29/29 donors (S4 rule; Y-gene rule for GSE227690). K0 20 eligible Oetjen donors (kept T 104-3,721 per donor),
5 of them >= 60 y. K4 primary genes in U_J 358/370 (0.968), 329/370 (0.889), 355/370 (0.959); U_J 11,417 / 10,542 / 11,522
genes. 3.5: gap_oof 1.378 (1.98 out-of-fold healthy SD), r = 1.9714 / 1.378 = 1.430; fold imprints 29-1,075 genes, no
empty half; error rates restated for r = 1.43: false-fail under H_AML ~7 %, false pass at a true age share of 40 % ~29 %,
of 50 % ~13 %. K6 Spearman 1.0000, max |diff| 0.00000 vs S8; |gap_in - 1.970| = 0.0014. K2 C1-C2 0.053, S1-S2 0.033
against (D_res_J/r)/3 = 0.459; second aspirates Ck 0.269, Sk1 0.074, Sk2 0.093. K5 full plant slope 0.978, contrast 0.885
vs 0.902 expected; third plant 0.327, 0.299 vs 0.301; planted A2b b/lambda 0.976 vs random P95 0.034; sign-scrambled
f_slope P5/P95 -0.036 / 0.030. K1 CD8 Naive share vs age rho = -0.767, one-sided p < 1e-4 (n 20); CD4 Naive rho = -0.245
(p 0.149, reported); Griffin 75 y 0.031 and 74 y 0.039 below every young donor (0.234, 0.101, 0.060); Rodriguez-Sevilla 61 y
0.048 and 70 y 0.047 below 0.089 and 0.149; Griffin author labels vs ours: CD8Naive 0.791, CD4Naive 0.910 agreement (author
CD8Naive share 0.100, ours 0.092); K1b naive genes -0.0131, cytotoxic +0.0212 log2 per year (expected signs). D_res_J
1.968 (Oetjen), 1.970, 1.962; gate denominator D_res_J/r = 1.376 (Oetjen).

**11.2 Phase 2.** A1: f_slope 0.102 (sex-adjusted; unadjusted 0.164), f_contrast 0.073 (old n 5, mean 67.4 y, vs young
n 10, mean 39.0 y) -> **f_age = 0.102 < 1/3, PASS.** 95 % donor-bootstrap CI -0.085 to 0.463 (conditional on the yardstick):
the upper bound is above 1/3, so the pass is worded as a bound. Leave-one-donor-out range 0.045-0.185: no flip. Matched
random sets (1,000): f_slope P5/P95 -0.032 / 0.046 with the imprint at rank 1.000, f_contrast P5/P95 -0.018 / 0.065 -> the
"weak age association exists" clause is issued. A2b: b 0.046, lambda 0.92, b/lambda 0.050 vs random P95 0.206 (random median
0.074 = genome-wide rate), Spearman 0.167, sign concordance 0.601, n 358 -> below P95: the imprint genes' age effects are not
more aligned with their AML effects than random genes'. A2a (reported): 29/370 genes at raw p < 0.05, 0 at FDR < 0.10;
removal ratio 1.007 (random removals 0.992-1.008); removed genes' age sign equals the AML sign in 0.62; robust set of 341
genes written as a record only. A3: Oetjen S_E 0.275 (< D_res/2 = 0.984), SD ratio 0.586 -> interpretable; A_E 0.100,
A_E/(D_res/r) 0.073; 0 of 5 older donors above the midpoint after study correction -> direction holds; D_dir 1.593. Donor N
(67 y) scores 0.031 raw, above the midpoint (-0.197) before study correction and 0.047 below it after (-0.244); the
at-most-one rule would have tolerated one such donor in any case. GSE227690 (3 young, 2 old: descriptive): S_E -0.120,
A_E -0.030 (-0.022 of D_res/r), 0 of 2 above. GSE276360 (2 / 2: descriptive): S_E 0.838 (0.43 D_res, inside the
interpretable band), SD ratio 1.53, A_E 0.407 (0.296), 0 of 2 above after correction; uncorrected, all four donors sit at
-0.91 to 0.23, i.e. the whole cohort is shifted toward AML by 0.84 units (D_dir 0.718).

**11.3 Reported analyses.** 402 set: f_slope 0.121, f_contrast 0.097, f_age 0.121 over its own D_res_J/r 1.373.
`score ~ age + sex + naive share (kept T cells)`: f_slope -0.047 -- with the naive-T share in the model the age slope is
gone, so the weak age association of the J score runs through T composition (the CD8 Naive share, K1), not through a
within-cell age programme, consistent with A2b. W space (Oetjen's own z): rho(score, age) 0.277, p 0.237, Fisher-z CI
-0.189 to 0.641, 0.122 W units per decade, young-donor SD 0.568, |rho| rank among random sets 0.847 (inside the random
range). Subtype pseudobulks: CD4 Naive rho 0.321 (p 0.18, n 19), CD8 EM rho 0.043 (p 0.86, n 18). Ambient partials given
age: hb_share_nonery -0.571, gran_share_lymph -0.020 (informational; HB genes are blacklisted; donor G has hb share 0.35,
104 kept T cells and scores -1.26). Section 3.4 inside GSE185381: sign concordance 0.903 of the 370 genes between the 4 AML
< 60 y and the 5 healthy 39-53 y; rho(score, age) within AML -0.074 (p 0.79), within healthy 0.261 (p 0.47); the four AML
< 60 score 0.48-1.07 against a healthy range of -2.27 to -0.43. Oetjen J scores: young -1.31 to -0.25 (mean -0.906),
middle mean -0.918, old -1.42 to 0.03 (mean -0.806); GSE185381 healthy -1.181, AML 0.787.

**11.4 Reading (section 5, filled from the tables; the sentence written to GATES.md).** GATE_S9 PASS. "In one healthy marrow
cohort (Oetjen 2018: 20 of 20 donors, 24-84 y, 10x 3' v2, one protocol, part of the BoneMarrowMap reference), the T imprint's
age shift over the GSE185381 age gap is 0.102 (95 % CI -0.085 to 0.463) of the cross-fitted AML-healthy gap, below one third
at the point estimate; upper bound 0.463 (worded as a bound). The age alternative to the direction is bounded, not excluded.
Direction versus healthy remains shown in the defining cohort only (no independent AML + healthy cohort exists in this work);
the 32 ambient-prone genes stay excluded from every statement. The S5 label 'AML-or-age' is narrowed to 'AML, with an age
share bounded below one third in healthy marrow (one protocol)'." Clauses: A2b below the random P95 (b/lambda 0.050 vs
0.206); a weak age association exists (f_slope 0.102 above the random P95 0.046); older healthy donors in GSE120221 do not
score like AML (A_E/(D_res/r) 0.073, n young 10, n old 5); GSE227690 (3 / 2) and GSE276360 (2 / 2) are descriptive, with
0 of 2 older donors above the midpoint in each.

What this does and does not establish. (i) Within healthy marrow, 31.5 years of age move the T-imprint score by about one
tenth of the AML-healthy gap, and that tenth is carried by the naive-T share (f_slope -0.047 once the share is in the model;
A2b shows no gene-level alignment beyond the genome-wide rate). (ii) The bound is loose: the bootstrap upper bound (0.463)
and the false-pass rate at r = 1.43 (~29 % at a true share of 40 %) leave an age share near 40 % compatible with the data;
the pre-registered line is the point estimate, and the wording says so. (iii) Nine healthy donors >= 60 y in three studies
(5 + 2 + 2) stay on the healthy side of the AML-healthy midpoint after study correction; one (Oetjen N, 67 y) is above it
before correction. (iv) No independent AML cohort enters, so direction versus healthy is still shown in GSE185381 only.
(v) One protocol (10x 3' v2), and Oetjen is inside the BoneMarrowMap reference: its projection is in-sample and K1 is an
age check, not an annotation check (Griffin's author labels agree at 0.79 / 0.91 on the naive subsets).

**11.5 Disclosures and deviations (complete list).** (1) Yardstick: the registered 1.970 = 3.87 healthy-SD omitted the
blacklist step; the exact logz recipe gives 1.9714 = 3.84 (K6 passes, 0.0014 < 0.01); r and the A2a ratio use 1.9714
(section 10.4). (2) The real run was made once and read the age relation once; before it, five smoke runs under permuted
ages and three independent reviews (section 10.4), including the smoke-table rows that carried real library ages next to
library scores and were overwritten unread (10.4 item 2). (3) The subtype rho used an unregistered >= 6-donor minimum
(reported only). (4) No threshold, gene set, donor list or wording was changed after a number was seen; sections 1-9 are as
tagged. Stage 2 (more older donors) is not triggered by a pass and is a separate decision; the next pre-registration in
the handoff is S10 (remission).

## 12. Erratum (appended 2026-10-28; sections 1-11 unchanged)

12.1 **Section 11.4 (iii) misstates the uncorrected positions.** It reads "one (Oetjen N, 67 y) is above it before
correction". From `results/tables/ws/ws_s9_direction_positions.csv` (job 3669091), three of the nine healthy donors aged
>= 60 y score above the AML-healthy midpoint before study correction: Oetjen N (67 y; 0.031 against a midpoint of -0.197),
Rodriguez-Sevilla HD2 (61 y; -0.068 against -0.196) and HD3 (70 y; 0.202 against -0.196). After study correction none of
the nine is above the midpoint (corrected scores -0.244, -0.906 and -0.636). Over all 29 donors of the three cohorts, 4
are above the midpoint before correction (the three above plus Rodriguez-Sevilla YHD636, 24 y) and 0 after. The corrected
statement of 11.4 (iii) is therefore: "Nine healthy donors >= 60 y in three studies (5 + 2 + 2) stay on the healthy side of
the AML-healthy midpoint after study correction; three of them (Oetjen N; Rodriguez-Sevilla HD2 and HD3) are above it
before correction. The Griffin and Rodriguez-Sevilla corrections rest on 3 and 2 young donors." Section 11.2 already gave
the per-cohort numbers correctly. No gate, test or wording of GATES.md depends on this sentence. The error was found by an
independent review of the WP-A1 evidence on 2026-10-28 (reports/2026-10_wpa1_paper_gap/).
