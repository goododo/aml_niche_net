# PRE-REGISTRATION — is the Mono_DC HOX signal reducible to the LMPP_GMP one?

Written 2026-09-24. Committed before `06_hox_compartment.R` exists.

---

## 0. Exactly what I have and have not already seen

This matters more than usual, because the previous step's results are already on the table.

**Already seen** (so nothing below may be presented as a first look):
- Mono_DC HOX score vs NPM1: GSE116256 4v3 perfect separation (pre-registered, p at the
  floor); GSE185381 7v25 AUC 0.891, p = 0.0004 (exploratory — the ranking was seen first);
  pooled 16v28 AUC 0.893.
- Mono_DC / LMPP_GMP / T_NK HOX score vs clinical blast %: `PREREGISTRATION_hox_axis.md` §7.1,
  rho −0.434 / −0.248 / −0.265 (GSE185381) and +0.500 / +0.100 / 0.000 (GSE116256).

**Not computed, not seen, and the whole subject of this document:**
- the correlation between the Mono_DC score and the LMPP_GMP score,
- the LMPP_GMP score vs NPM1,
- any adjusted or partial association between the two compartments.

---

## 1. The question, and why it is the one that matters

`hox_axis` left Q1 at **INDETERMINATE**: the Mono_DC HOX score does not track clinical blast %
(rho −0.434 at n=17 and +0.500 at n=7 — opposite signs). So the simple story, "the AML Mono_DC
bin is blasts and the score counts them", is unconfirmed. It is also unrefuted, because this
project has no working malignancy call (inferCNV rho +0.14/+0.04, van Galen −0.165 on that same
ruler).

The next decidable question does not need a malignancy call at all:

> **Does the HOX signal measured in the Mono_DC compartment carry information that the
> LMPP_GMP compartment — where blasts legitimately live — does not already carry?**

- **No** → Mono_DC is a diluted copy of the blast compartment. The Mono_DC line closes honestly,
  and what remains is a working NPM1 classifier built without any malignancy annotation.
- **Yes** → the Mono_DC compartment is not reducible to the blast compartment. That is a claim
  about compartments that this data can actually support, and it is the first thing in this
  project that would be a finding rather than a validation.

---

## 2. Frozen inputs — nothing is recomputed

| Input | Path |
|---|---|
| per-(sample,bin,set) HOX scores + `n_cells` + `lib_size` | `results/tables/12_pseudobulk_de/hox_axis_scores.csv` (written 2026-09-18) |
| NPM1 status | `results/tables/01_preprocess/00_curated_manifest.csv`, `driver_mutations`, `grepl("NPM1", .)` |

Scores are taken verbatim from the frozen file. Set **L** (locus-defined, 19 genes in Mono_DC)
is used throughout, as in `hox_axis` §4.1. No new pseudobulk, no new normalisation, no new gene
set. If either file's mtime differs from the value recorded at first run, the script stops.

### 2.1 The roster, counted before any test

Samples scored in **both** Mono_DC and LMPP_GMP **and** carrying a mutation call — AML only:

| dataset | n | NPM1-mut | in the analysis? |
|---|---|---|---|
| GSE185381 | 29 | 5 | yes |
| GSE116256 | 7 | 4 | yes |
| Chen2023 | 5 | **5** | **no — see below** |
| **total used** | **36** | **9** | |

**Chen2023 is excluded, decided here.** All 5 of its samples are NPM1-mutant, so it contributes
zero within-dataset contrast and would separate the outcome perfectly on the dataset term alone.
Its exclusion is a consequence of its composition, not of its values, and is therefore not an
outcome-dependent choice.

Adding T_NK (the §5.1 control) reduces the roster to **30 samples, 7 NPM1-mut** (GSE185381 26/5,
GSE116256 4/2). The control is run on that smaller roster and its n is reported beside it.

**9 events over 36 samples.** This is too few for an asymptotic two-predictor model, which is why
§4 uses a rank statistic with a permutation null rather than a logistic-regression p-value.

**The roster is not the one already seen.** The AUC 0.891 in §0 was measured on 32 GSE185381
samples with 7 NPM1-mutants. Requiring an LMPP_GMP score as well drops that to 29 with **5**
mutants: two of the seven mutants do not clear the 30-cell gate in LMPP_GMP. So the already-seen
number does not transfer, and the test below is weaker than it, not a restatement of it.

---

## 3. Precondition: can the question be answered at all?

Spearman(Mono_DC score, LMPP_GMP score) across the 36 samples, within dataset.

**If |rho| > 0.90 the two compartments are collinear and no statistic can separate their
contributions at n=36. The analysis then stops and reports UNANSWERABLE.** Deciding this in
advance is the point: a collinear pair would otherwise yield an arbitrary-looking split of
credit that reads like a result.

---

## 4. The statistic

Logistic regression is deliberately avoided: at 9 events with 2 correlated predictors, separation
and unstable coefficients are likely, and the p-value would be the fragile part of the answer.

**Adjusted score.** Within each dataset, rank-transform both scores, regress the Mono_DC ranks on
the LMPP_GMP ranks by ordinary least squares, and keep the residual. That residual is the part of
the Mono_DC signal that the LMPP_GMP signal does not explain.

**Primary statistic.** AUC of the LMPP-adjusted Mono_DC residual for NPM1 status, computed on the
pooled 36 samples (residuals are within-dataset, so pooling them is on a common scale).

**Symmetric comparator, given equal standing.** AUC of the Mono_DC-adjusted LMPP_GMP residual for
NPM1 status. The comparison of the two is the result; neither is privileged.

**Null.** NPM1 labels permuted **within dataset**, 10,000 times, recomputing the full adjustment
each time. One-sided (adjusted score higher in NPM1-mut), because the direction is fixed by the
external literature and by the already-seen unadjusted result, not chosen here.

**Reported beside them for reference, not tested:** the unadjusted Mono_DC AUC and unadjusted
LMPP_GMP AUC on the same 36 samples.

---

## 5. Decision rule, fixed now

Let `p_mono` and `p_lmpp` be the two permutation p-values.

| outcome | condition | what it means |
|---|---|---|
| **MONO_ADDS** | `p_mono < 0.05` and `p_lmpp >= 0.05` | Mono_DC carries non-redundant information. The strongest outcome available. |
| **LMPP_ONLY** | `p_lmpp < 0.05` and `p_mono >= 0.05` | Mono_DC is the diluted copy. The Mono_DC line closes. |
| **BOTH_ADD** | both `< 0.05` | Two partially independent readouts of the same axis. Reported as such; neither compartment is called primary. |
| **REDUNDANT** | neither `< 0.05` | No separation demonstrable at n=36 / 9 events. Not evidence that Mono_DC is redundant — see §7. |
| **UNANSWERABLE** | §3 precondition fails | Collinear; the question is not decidable here. |

### 5.1 Controls, each able to change the reading

1. **T_NK negative control.** The identical adjusted test with T_NK in place of Mono_DC, on the
   30-sample roster. **If T_NK also "adds", the adding is generic** — ambient RNA, doublets, or
   the adjustment itself manufacturing residual signal — and MONO_ADDS is downgraded to
   BOTH_ADD-with-a-failed-control, which supports nothing.
2. **Cell-count confound.** Spearman(Mono_DC `n_cells`, NPM1) and the same for LMPP_GMP. A bin
   with more cells is measured more precisely and can "add" for that reason alone. If either
   |rho| > 0.40 the primary test is repeated with log10 `n_cells` of both bins added to the
   adjustment regression, and **both versions are reported**.
3. **Set D repeat.** The whole §4 procedure on Set D (the 12 Discovery-frozen genes) instead of
   Set L. Set L and Set D agreed at rho 0.984 in `hox_axis` §7.4, so a disagreement here would
   mean the adjustment, not the axis, is driving the answer.

---

## 6. What will not be claimed, whatever comes out

**This is the section that bounds the paper.** MONO_ADDS would show that the Mono_DC
compartment's HOX signal is not reducible to the LMPP_GMP one. It would **not** show that the
cells carrying it are non-malignant. Three explanations survive a positive result:

1. non-malignant monocytes in AML marrow carry a HOX-linked state change;
2. the blasts sitting in the Mono_DC bin are a **different** blast population from the LMPP-like
   blasts, and therefore carry different information — still "it is blasts";
3. the two bins differ in measurement precision (addressed only partly by §5.1 control 2).

**Explanation 2 cannot be excluded without a working cell-level malignancy call, which this
project does not have.** Therefore the claim ceiling is a statement about **compartments**, not
about malignancy, and the wording is fixed here:

> *"The HOX-cluster signal measured in the monocyte/DC compartment is not reducible to the signal
> in the LMPP/GMP compartment."*

Also not claimed: anything about NPM1 alone (the FLT3-ITD and DNMT3A co-occurrence confounds from
`hox_axis` §6 carry over unchanged); risk group, treatment response or survival; and no
per-gene statement of any kind — this step tests two scalars per sample.

---

## 7. Power, stated before the result

9 NPM1-mutants over 36 samples, two correlated predictors. **REDUNDANT is the outcome this design
is most likely to return even if Mono_DC genuinely adds a modest amount**, and it will be reported
as "not demonstrable at n=36", never as "shown to be redundant". The permutation null's minimum
attainable one-sided p is 1/10,001 = 9.999e-5, which is not the binding constraint; the binding
constraint is 9 events.

---

## 8. Files

| File | Role |
|---|---|
| `scripts/12_pseudobulk_de/PREREGISTRATION_hox_compartment.md` | this document |
| `scripts/12_pseudobulk_de/06_hox_compartment.R` | the only new script; one file, no helpers |
| `results/tables/12_pseudobulk_de/hox_compartment_tests.csv` | every test in §3–§5 with its n and decision label |
| `results/tables/12_pseudobulk_de/hox_compartment_null.csv` | the §4 permutation distributions |

No new dependencies: base R + `data.table`, both already used by `01`–`05`.

## 9. Order of execution

1. This document committed.
2. `06_hox_compartment.R` written; §3's precondition printed **first**, and the script exits
   there if it fails.
3. One run. §5.1's controls printed before the §5 decision label is assigned.
4. Result recorded whatever it is, with §6's wording ceiling and §7's power statement attached.
