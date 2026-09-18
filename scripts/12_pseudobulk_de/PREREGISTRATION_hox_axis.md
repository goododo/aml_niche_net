# PRE-REGISTRATION — the HOXA/HOXB–MEIS1 axis: is it blast burden, and does it track NPM1?

Written 2026-09-18. Committed before `05_hox_axis.R` exists.

---

## 0. This hypothesis is post-hoc in origin, and that is the reason for this document

The Mono_DC result in `PREREGISTRATION_pseudobulk_de.md` was pre-registered and is
closed. **The hypothesis below was formed on 2026-09-18 by looking at which genes
came out of it.** That is legitimate hypothesis generation and illegitimate
evidence. Everything downstream of this file is pre-registered, or it is nothing.

Nothing here amends the pseudobulk pre-registration. Its hit lists, its universes,
its Discovery/Validation boundary and its conclusions are **frozen inputs** to this
document and are not recomputed, re-thresholded, or re-interpreted by it.

---

## 1. The observation that triggered this (stated once, not re-litigated)

`04_validation.R`, run 2026-09-18, on the held-out Validation arm (GSE116256,
Seq-Well; the Discovery arm is Chen2023 + GSE185381, 10x):

| | |
|---|---|
| Mono_DC frozen hits | 385 (of 537 baseline; 152 struck by the depth arm) |
| permutation p | 0.002 (`sample_within_dataset`), 0.001 (`library_arm_pure`) |
| Validation sign concordance | 163 / 211 = 77.3%, one-sided binomial p = 3.8e-16 |

The genes at the top of the up-in-AML side are the HOXA and HOXB clusters plus
MEIS1 (`HOXB-AS3` +6.6, `HOXA3` +4.6, `HOXB5` +4.4, `HOXB3` +4.3, `HOXB6` +4.3,
`MEIS1` +4.2, `HOXA9` +4.0, `HOXA10` +3.4, `HOXB4` +3.3, `HOXB7` +3.0, plus `MYB`
+1.9 and `FLT3` +1.7), every one replicating with the same sign across platform.
The down side is mature monocyte / cDC2 identity (`FCER1A` −3.4, `CLEC10A` −3.3,
`CLEC4D` −3.4, `CLEC4E` −3.3, `FCGR3A` −2.8, `FPR1`, `FPR2`, `G0S2`, `RETN`).

The pseudobulk pre-registration §10 states the relevant design decision plainly:
*"Malignant vs non-malignant. Not a split here. inferCNV under-calls ~9x; every
§4.2-surviving cell in a bin enters its pseudobulk regardless of CNV call."* So
the leading candidate explanation is that the AML `Mono_DC` bin carries monocytic
blasts and the contrast is reading **cell identity, not niche state**.

**This document tests that explanation rather than assuming it.**

---

## 2. The two questions

**Q1 — Is the axis blast burden?**
Score each sample on the axis, correlate with the clinical blast percentage, the
one ruler in this project that the expression data never saw.

**Q2 — Does the axis track NPM1 status?**
HOXA/HOXB–MEIS1 is high in NPM1-mutant and MLL-rearranged AML and low in
TP53-mutant and core-binding-factor AML. GSE116256 carries a mutation call for
all 9 of its Validation AML samples. The gene set is derived **entirely from
Discovery**; GSE116256 is the held-out arm. This test is therefore non-circular
by construction.

### 2.1 What each outcome buys — fixed in advance, including the bad case

| Q1 | Q2 | Reading |
|---|---|---|
| BLAST-PROXY | — | The Mono_DC DE is re-labelled cell-identity contamination. In exchange the project gets its **first working blast proxy**: inferCNV `malignant_frac` gives rho ~ 0 against this same ruler and the van Galen signature caller gave rho = −0.165 (2026-09-16). A Discovery-derived, externally-validated, annotation-free proxy is a deliverable. |
| NOT-BLAST | SUPPORT | The axis is not proportional to blast fraction but separates a genotype. That is a patient-stratification axis and the strongest available seed for an article. |
| NOT-BLAST | NULL | An axis that is neither blast burden nor NPM1. Reported, not built on, until something external anchors it. |
| BLAST-PROXY | SUPPORT | Expected if the proxy works: NPM1-mut AML is both blast-rich and HOX-high. Both readings stand; neither is upgraded by the other. |
| INDETERMINATE | any | No claim. n is reported and the question is left open. |

The fourth and fifth rows exist so that "怎么答都赢" is a statement about the
design, not a licence to narrate whatever comes out.

---

## 3. Frozen inputs

| Input | Path | Frozen because |
|---|---|---|
| per-sample pseudobulk counts | `PB_RDS_DIR/<ds>/<sample>__pseudobulk.rds` (67 files) | written 2026-09-16 by `01_aggregate.R` |
| roster | `results/tables/12_pseudobulk_de/pb_sample_manifest.csv` | 67 `in_analysis`; Discovery 42 AML + 14 healthy, Validation 9 AML + 2 healthy |
| frozen hit list | `results/tables/12_pseudobulk_de/hitlist_frozen.csv` | Discovery-only (`02_de_limma.R:82`) |
| per-bin tested universe | `results/tables/12_pseudobulk_de/de/*.csv` | defines which score genes are available |
| clinical blast % | `results/tables/05_ccc/undercall_contamination.csv` | 59 samples; 30 overlap the pseudobulk roster (21 Discovery + all 9 Validation) |
| mutation calls | `00_project/metadata/GSE116256__samples.tsv` | all 9 Validation AML annotated |

No input is regenerated. If any file's mtime differs from the value recorded at
first run, the script stops.

---

## 4. The score — every degree of freedom is fixed here

**4.1 Two gene sets, both specified before any outcome is touched.**

- **Set L (primary, locus-defined, external).** Every gene in the bin's tested
  universe whose symbol begins `HOXA` or `HOXB`, plus `MEIS1` and `PBX3`.
  Membership is a genomic-locus rule, not a p-value rule, so it cannot be accused
  of cherry-picking from our own result. In Mono_DC this is **19 genes**
  (`HOXA-AS3 HOXA10 HOXA10-AS HOXA3 HOXA4 HOXA5 HOXA6 HOXA7 HOXA9 HOXB-AS1
  HOXB-AS3 HOXB2 HOXB3 HOXB4 HOXB5 HOXB6 HOXB7 MEIS1 PBX3`); in LMPP_GMP 17; in
  T_NK 7; in B_Plasma 2.
- **Set D (secondary, data-derived).** The same rule applied to the Mono_DC
  **frozen hit list** — **12 genes**, dropping `HOXA4 HOXA5 HOXA6 HOXA-AS3
  HOXB-AS1 HOXB2 PBX3` and adding none.

Set L is primary precisely because it includes genes our DE did *not* pick.

**4.2 Normalisation.** log2 CPM computed within the bin, using that bin's own
column sum as the library size, prior count 1. No cross-bin pooling.

**4.3 Standardisation.** Each gene is z-scored **within dataset** across all
samples of that dataset scored in that bin, healthy included — healthy anchors
the low end. Dataset is the dominant nuisance everywhere in this project, and
both tests below are within-dataset, so within-dataset z costs nothing.

**4.4 Missing genes.** A score gene absent from *any* sample's pseudobulk in a
dataset is dropped for that **whole dataset**, so every sample in a dataset is
scored on an identical gene list. The list actually used is written out per
(dataset, bin, set).

**4.5 Score.** Unweighted mean of the z-scores. No PCA, no eigengene, no
weighting by logFC — a weighted score would import the Discovery effect sizes and
re-open the circularity Set L was chosen to close.

**4.6 Bins scored.** `Mono_DC` (primary), `LMPP_GMP` (positive control — blasts
legitimately live there), `T_NK` (negative control). `B_Plasma` is **not** scored:
2 available genes is not an axis.

---

## 5. Q1 — blast burden

**Test.** Spearman between the Mono_DC Set-L score and `blast_pct_clinical`,
**within dataset**, AML samples only (healthy carry no blast %).

Testable strata: **GSE116256 n=9** (primary — held out, gene set never saw it)
and **GSE185381 n=18** (secondary — selection-adjacent: it contributed to the
AML-vs-healthy contrast, though never to any blast-% comparison). Chen2023 n=3 is
reported and not tested.

**Decision rule, fixed now:**

- **BLAST-PROXY** — rho ≥ +0.60 in **both** strata.
- **NOT-BLAST** — |rho| < 0.30 in **both** strata.
- **INDETERMINATE** — anything else, including discordant signs, one stratum
  passing and the other not, and everything between 0.30 and 0.60. Declared
  indeterminate at the stated n; **no directional claim is made in this case.**

**Reported beside it on the identical samples**, so the comparison is like for
like: Spearman of `malignant_frac` (inferCNV) with the same blast %, and the
recorded van Galen figure (rho = −0.165, 2026-09-16).

**Minimum detectable effect.** By exact permutation, n=9 reaches two-sided 0.05
only at **|rho| >= 0.700**; at n=18 the boundary is **|rho| ~ 0.47**. Both are
printed next to the result. Note that the BLAST-PROXY threshold of 0.60 in the
GSE116256 stratum is therefore *below* that stratum's own significance boundary:
the rule is deliberately an effect-size rule, not a significance rule, because at
n=9 a significance rule would only ever fire on near-perfect monotonicity.

---

## 6. Q2 — NPM1

**Samples.** All 9 GSE116256 Validation AML, mutation calls verbatim:

| sample | NPM1 | other |
|---|---|---|
| AML210A-D0 | **W288fs** | DNMT3A.R882C |
| AML329-D0 | **W288fs** | **FLT3.ITD**, SMC3 |
| AML419A-D0 | **W288fs** | **FLT3.ITD**, FLT3.A680V, FLT3.N841K, DNMT3A.R882C |
| AML556-D0 | **W288fs** | DNMT3A.R882C, NRAS.Q61H, TET2 ×2, ATM |
| AML1012-D0 | wt | KRAS.G13D, NRAS.G13D, SF3A1 |
| AML328-D0 | wt | TP53.P152R, TP53.Q144P |
| AML707B-D0 | wt | KIT.Y823C, RAD21, BRCC3 |
| AML916-D0 | wt | TP53.C238Y |
| AML921A-D0 | wt | DNMT3A.R882H, RUNX1, SETD2 |

**4 mutant vs 5 wild-type.**

**Test.** One-sided Wilcoxon rank-sum, NPM1-mut > NPM1-wt. The direction is
pre-specified from external literature, not from our data, which is what makes
one-sided legitimate here.

**Minimum attainable p.** C(9,4) = 126 orderings, so the smallest one-sided p is
**1/126 = 0.0079**. Significance at 0.05 therefore requires near-complete
separation. This is stated before the test, not after.

**Decision rule:**

- **SUPPORT** — one-sided p < 0.05 **and** all 4 mutants above the wild-type median.
- **NULL** — otherwise. **A null at 4 vs 5 is not evidence of absence** and will
  not be written as one.

**Confound declared in advance.** 2 of the 4 NPM1-mutants (AML329-D0, AML419A-D0)
are also FLT3-ITD. NPM1 and FLT3-ITD cannot be separated at this n. Any positive
result is reported as **"NPM1-mutant and/or FLT3-ITD"**, never NPM1 alone.

**Direction check, not a test.** The 2 TP53-mutants (AML328-D0, AML916-D0) are
expected HOX-low. Their ranks are reported; n=2 supports no inference.

**The table is the result.** At n=9 all nine per-sample scores are printed in
full, so the reader can see the separation rather than take a p-value for it.

---

## 7. Controls — each one can kill the axis on its own

1. **Bin specificity.** If the signal is blasts sitting in the Mono_DC bin,
   LMPP_GMP should score high and track blast %, and **T_NK should not**. If T_NK
   tracks blast % about as well as Mono_DC, the axis is ambient RNA or doublets
   and both bins' readings are artifacts. **KILL.**
2. **Depth.** Spearman(score, log10 bin library size) and Spearman(score, cells
   in bin). **|rho| > 0.50 on either → the score is a depth readout. KILL.**
3. **Healthy floor.** Mono_DC score, healthy vs AML, reported as AUC. If healthy
   are not clearly lower, the score does not measure what the name says.
4. **Set agreement.** Spearman(Set L score, Set D score) across samples must be
   **≥ 0.80**. Below that the two definitions are not describing one axis and
   neither is carried forward. **KILL.**
5. **Random-set null.** 1000 random gene sets drawn from the Mono_DC tested
   universe, matched to Set L on size and mean-expression decile. The observed
   Q1 Spearman must exceed the **95th percentile** of that null. Otherwise any 19
   genes would have done as well and the HOX identity is decoration.

Controls run and are reported **whatever Q1 and Q2 return**, including when both
are null.

---

## 8. What will not be claimed, whatever comes out

- **Not a cell-level malignancy caller.** This scores samples. Nothing here
  assigns a cell.
- **Not a niche or microenvironment finding.** If Q1 returns BLAST-PROXY, the
  Mono_DC DE is explicitly re-labelled cell-identity contamination in
  `FINDINGS_*`, and the 385 genes are not presented as monocyte biology.
- **Not risk group, treatment response, or survival.** None of these exist in
  this project — karyotype is 0/244 and mutations 21/244, GSE116256 only.
- **Not NPM1 alone**, per §6's FLT3-ITD confound.
- **No p-value appears without its n beside it**, and 9 and 4-vs-5 are small.

---

## 9. Files

| File | Role |
|---|---|
| `scripts/12_pseudobulk_de/PREREGISTRATION_hox_axis.md` | this document |
| `scripts/12_pseudobulk_de/05_hox_axis.R` | the only new script; one file, no helpers |
| `results/tables/12_pseudobulk_de/hox_axis_scores.csv` | one row per (sample, bin, set): score, genes used, bin depth, cells |
| `results/tables/12_pseudobulk_de/hox_axis_genes.csv` | the gene list actually used per (dataset, bin, set) |
| `results/tables/12_pseudobulk_de/hox_axis_tests.csv` | every test in §5–§7 with its n, rho/p, and decision label |
| `results/tables/12_pseudobulk_de/hox_axis_null.csv` | the §7.5 random-set null |

It lives in `12_pseudobulk_de/` because it consumes that stage's frozen objects
and adds no new stage to the pipeline. **No new dependencies**: base R,
`data.table`, `Matrix` — all already used by `01`–`04`.

---

## 10. Order of execution

1. This document committed.
2. `05_hox_axis.R` written, printing §7's controls **before** §5 and §6, so a
   killed axis is never scored against an outcome.
3. One run. Results recorded whatever they are, in a new section of
   `FINDINGS_topology_null.md` or a new findings file, with §2.1's table as the
   frame.
4. If a control kills the axis, §5 and §6 are **not** run and that is the result.
