# Pre-registration v1: within-patient -X vs X-retained leukaemic program in a female t(8;21) AML (LOX-1) (2026-11-26)

Status: v1. Design and file split approved 2026-11-26 ("预注册草稿和文件划分可以"). Frozen copy, tagged `prereg-lox-1-v1`
before any script reads an expression value by X status. Sections 1-9 are the approved draft; the only edit is the order of
the two smoke worlds in section 9 (see section 10, item 1). Section 10 was added at the freeze.

What has been read: per-cell XIST UMIs, chrY UMIs, depth and labels of every female sample (census scripts s10-s13 under
`reports/2026-11_monthly_PI/pi_questions/DATA_D/`); the donor's diagnostic karyotype (curated manifest); and every result of
the -Y line (H2 and its post-hoc checks), which is where the two hypotheses come from.
What has not been read: any gene other than XIST and the chrY genes, split by XIST status, in any female sample.

Everything not restated here is as `scripts/loy/PREREGISTRATION_H2.md` (v2): scoring, model, permutation scheme, constants.

## 1. Question

In the -Y branches of the two male t(8;21) cases (AML15, AML13), HLA class II genes were higher and the HALLMARK
TNF-alpha/NF-kB set was lower than in the Y-retaining leukaemic cells of the same patient. Is either observation repeated in
the -X branch of a female t(8;21) leukaemia?

- H-X1: the HLA-II score is higher in -X than in X-retained leukaemic cells.
- H-X2: the TNF-alpha/NF-kB score is lower in -X than in X-retained leukaemic cells.

Confirmatory and directional, one patient. It is the only independent test the public data allow for the HLA-II result
(the -Y cases are capped at three).

Reference effects from the -Y line (score-level d, same model; 05e, post-hoc):

| | AML15 | AML13 | PT43 (TP53, not t(8;21)) |
|---|---|---|---|
| HLA-II (11 structural genes) | +0.32 | +0.17 | +0.90 |
| F7 TNF-alpha/NF-kB | -0.30 | -0.21 | +0.20 |
| F7 without the genes shared with MP9 | -0.14 | -0.22 | +0.18 |

## 2. Case (chosen on counts, disclosed as such)

GSE185381, pool P06 (three lanes, 10x 3'). The cells labelled PAUMTZ are donor PAXMIJ (label swap, `config_loy.R` section 6).
Diagnostic karyotype: t(8;21) in all 14 abnormal metaphases; `idem,-X` in 7; `idem,del(9q)` in 5; 6 normal metaphases.
scRNA census: XIST-zero excess 0.685 in the leukaemic compartment and 0.67 in Mono_DC; 0.000 in T cells (n = 40).

Why only this case: 4 readable female t(8;21) patients in the cohorts with fusion data. AML005 is clonal (0.91, 18 X-retained
cells: no comparison arm), AML055 has no X loss, Lambo AML16 is at 0.06. Lambo AML14 is clonal by karyotype.

Counts under the rule of section 3 (XIST and depth only):

| fine label | -X | X-retained | median depth -X / X-retained |
|---|---|---|---|
| Early GMP | 534 | 165 | 6,209 / 6,630 |
| GMP-Neut | 291 | 123 | 5,774 / 5,954 |
| GMP-Mono | 258 | 75 | 5,279 / 5,706 |
| CD14 Mono | 226 | 53 | 4,339 / 4,269 |
| total (eligible labels) | 1,309 | 416 | |

Of these 1,725 cells, 1,719 carry the depositors' malignant call. Eligibility as H2: a fine label enters with >= 50 cells on
each side after matching.

## 3. Labels and misclassification

- Cells: the depositors' malignant call, depth 3,000-11,999 UMIs, no chrY UMI.
  - The depth window was chosen after seeing that the XIST-zero fraction is flat inside it (0.73-0.75) and falls above 12,000.
  - The other two donors of the pool are male; 222 chrY-positive cells under this label are removed.
- -X (group A): 0 XIST UMIs. X-retained (group B): >= 2 XIST UMIs. One UMI: dropped.
- X-retained cells called -X (XIST dropout): 6-15 % of X-intact leukaemic cells have 0 XIST in this window (AML3948 0.056,
  AML055 0.06, AML012 0.155), which puts about 2-7 % of group A on the wrong side.
- -X cells called X-retained (ambient XIST): 12 of 1,040 cells in the clonal -X case AML005 (1.2 %), about 4 % of group B.
- Both errors pull a true difference toward zero. Their sum c is reported; d_corr = d / (1 - c).
- XIST is removed from every feature. All chrX and chrY genes are removed from the feature scores and flagged in the
  genome-wide table, because one X fewer changes their dose directly.

## 4. Matching and model

- Within fine label, each X-retained cell (the smaller side here) is matched to up to 3 -X cells by nearest log10 depth,
  caliper 0.1, without replacement (`match_depth()` of 04, sides swapped).
- Model on the matched set: score ~ label + fine label + lane + log10 depth. d = label coefficient / residual SD; d > 0 means
  higher in -X.
- p: 1,000 label permutations within fine label x depth bin. Holm over the two primary features.

## 5. Features

Primary (two, direction fixed in section 1):
- P1 HLA-II: HLA-DRA, -DRB1, -DRB5, -DPA1, -DPB1, -DQA1, -DQB1, -DMA, -DMB, -DOA, CD74 (the list of 05c).
- P2 F7: HALLMARK_TNFA_SIGNALING_VIA_NFKB from the pinned GMT.

Descriptive (no claim rests on them): F7 without the MP9-shared genes; the NF-kB core list of 05e; MP9 (handling program);
MP5; interferon alpha / gamma; MHC class I (HLA-A, -B, -C, -E, -F, B2M); the rest of the H2 v2 panel; genome-wide table S1
with chrX genes flagged. If MP9 moves with |d| >= 0.2, H-X2 is read on "F7 without MP9-shared genes" and said so.

Dose readout for the UTX/UTY question (descriptive): per-gene -X / X-retained ratio for KDM6A, KDM5C, DDX3X, EIF1AX, ZFX and
the PAR1 genes (CSF2RA, IL3RA, CD99, SLC25A6), next to the -Y ratios of the male cases.

## 6. Gates (in order; a failed gate is a stop for discussion)

- GATE_LX_0 (counts): >= 2 eligible fine labels and >= 300 matched X-retained cells. Met on the counts above.
- GATE_N1 (null world): labels shuffled within fine label x depth bin, 200 times; pooled rejection rate at 0.05 inside
  0.03-0.07 and every feature inside the binomial band 0.020-0.080 (the band H2 amendment 14 showed 200 shuffles can resolve).
- GATE_P (planted world): d = 0.3 planted into a random 60-gene set is detected at Holm p < 0.05 in >= 80 % of 200 plants.
  If it fails, a null result can only be called "not readable".
- GATE_LX_K2 (is the label a lost chromosome, not a silenced XIST): the PAR-gene score (the 29 PAR genes of
  `gene_sets_h2_v2.tsv`, those detected in >= 5 % of cells) is lower in group A, permutation p < 0.001. PAR genes are expressed
  from both X chromosomes, so one X fewer means a lower dose; in -Y cells CD99 and SLC25A6 fell to 0.43-0.59 of the
  Y-retained level (05e). If it fails, the contrast is reported as "XIST-negative vs XIST-positive" and no statement about
  X loss is made.
- GATE_LX_K1 (dropout control): the identical labels, matching and model on X-intact female AML samples of the same dataset,
  where XIST-zero is dropout.
  - Samples by rule: XIST-zero excess < 0.03 in both LEUK and Mono_DC, and >= 50 XIST-zero cells in the eligible labels and
    window. This gives AML3948 (179 zero / 3,043) and AML012 (313 / 1,708).
  - Sample x fine label is the stratum.
  - A primary feature whose control d has the same sign as its main estimate with p < 0.05 and |d| >= 0.2 is "removed by K1".
  - The control is conservative: there both sides are selected on XIST capture, in the case only the X-retained side is.

## 7. Wording (fixed now), per primary feature

- "Repeats": Holm p < 0.05, expected direction, |d| >= 0.2, same sign in every eligible fine label, not removed by K1.
- "Same direction, small": Holm p < 0.05, expected direction, |d| < 0.2.
- "Opposite": Holm p < 0.05 in the other direction.
- "No difference": 95 % CI inside (-0.2, 0.2) and GATE_P passed.
- Otherwise "not resolved".

Sensitivity runs, reported next to the primary and never instead of it: (a) without the CD14 Mono label (progenitor labels
only, as in the male t(8;21) cases); (b) depth below 8,000; (c) with the S/G2M score as covariate.

What "repeats" would and would not mean: the same direction in a third t(8;21) patient, with a different sex chromosome lost.
It would not show that losing the chromosome causes it. The X-retained side carries del(9q) in 5 metaphases, so this is a
branch contrast, as in AML15.

## 8. Known weaknesses (stated now)

- One patient.
- Pooled lanes at 4,000-6,500 UMIs per cell, against 9,000-12,000 in the Lambo libraries.
- The label rests on one gene (XIST) where the -Y label rests on about twelve, so dropout misclassification is larger.
- Doublet asymmetry: a blast paired with one of the donor's own normal XIST-positive cells lands in the X-retained group.
  That pushes HLA-II toward "lower in -X", against H-X1. Depositor doublet removal and the depth window limit it; it is not
  measured.
- Leukaemic cells are the depositors' call; the T-cell control of H2 is unavailable (40 T cells).
- The dataset has antibody tags. Whether HLA-DR is in the panel is not verified; protein is outside this document.

## 9. File split (for confirmation before any code is written)

| file | role | writes |
|---|---|---|
| `scripts/loy/PREREGISTRATION_LOX.md` | frozen copy of this document; tag `prereg-lox-1-v1` | - |
| `scripts/loy/config_loy.R` | one appended section 9: sample ids, depth window, XIST thresholds, the gene lists of section 5 | - |
| `scripts/loy/07_h2_lox_contrast.R` | one script: labels and matching from the per-cell table; counts from the `04_annotated` objects of PAUMTZ, AML3948, AML012; scores, model, permutations, K1 / K2 / N1 / P; `--world real\|null\|planted`; S1 | `results/tables/loy/05_lox/` |

Run with the existing `submit_loy_h2.sbatch 07`. No new dependency. Order: tag, planted world, null world, real run, report.
Smoke outputs carry `_smoke` names and are deleted after the real run, on request.

## 10. Implementation details fixed at the freeze (2026-11-26, before any expression value is read by X status)

None of these changes a rule of sections 1-7; they say how each rule is computed.

1. Smoke order: the planted world runs before the null world, so that the null world can reach the "no difference"
   branch (which needs GATE_P). The draft had the two the other way round.
2. Both smoke worlds shuffle the labels of both arms within stratum x depth bin before any expression is loaded, and
   write no S1, dose or sensitivity table (H2 amendment 9).
3. Gene filter: genes detected in >= 5 % of the donor's malignant-called, chrY-negative cells (any depth, any XIST count).
   The same gene list is used in the control arm.
4. Matching direction: the smaller side is the anchor, with up to 3 partners. Case: X-retained cells are anchors. Control:
   XIST-zero cells are anchors. A control stratum (sample x fine label) enters with >= 10 matched XIST-zero cells.
5. Control samples are named in `config_loy.R` (AML3948, AML012, the outcome of the section 6 rule); the script stops if
   either has fewer than 50 XIST-zero cells under the section 3 filter.
6. Model: lane is a fixed effect; the design matrix is reduced to its rank (lane is nested in sample in the control arm).
   Permutation strata: fine label (control: sample x fine label) x 10 depth bins.
7. "Same sign in every eligible fine label": per-label fits (label + lane + log depth) in labels with >= 50 / 50 matched cells.
8. GATE_N1 is evaluated on the two primary features, 199 permutations per shuffle as in H2: each rate inside 0.020-0.080 and
   the pooled rate inside 0.03-0.07.
9. GATE_P: Holm over the two primary features plus the planted set, 1,000 permutations per plant.
10. The MP9 clause of section 5: if |d| of MP9 is >= 0.2, the second primary feature is "F7 without MP9-shared genes"
    (expected lower), Holm is taken over it and HLA-II, and the all-gene F7 row is still printed.
11. GATE_LX_K1 uses the unadjusted permutation p of the control arm.
12. Misclassification, counts only. e and r2 = XIST-zero and XIST >= 2 rates among the control samples' cells; a = XIST >= 2
    rate in AML005; n0 and n2 = the case's XIST-zero and XIST >= 2 counts, all under the section 3 filter in the eligible
    labels. c_A = (n2 / r2) * e / n0, c_B = n0 * a / n2, c = c_A + c_B.
13. Sensitivity (a) and (b) are subsets of the matched set with the model unchanged; (c) adds the S and G2M scores as two
    covariates. Scores are not recomputed for a subset.
14. Dose readout: raw detection and UMIs per 10,000 for the listed genes, in the case arm and in the control arm.
15. Seed: LOY_H2_SEED for the matching, so the matched set is identical in every world; then LOY_H2_SEED + the length of
    the world name (as 05) for shuffles and permutations.
