# Pre-registration 2: the genotype-conditioned design (NPM1-like within FLT3-ITD-negative AML)

**Status: DRAFT — awaiting user approval. No number in §3–§6 has been computed.**
Written 2026-09-29, immediately after `PREREGISTRATION_bulk_transfer.md` §7–§8 closed the
unconditioned design, and BEFORE the third cohort is joined to any outcome. This document exists
because the previous design's failure mode was diagnosed, not because a new hypothesis appeared.

## 0. Honesty ledger — everything already seen

This is the third pass over the same idea, so the ledger matters more than usual.

**Seen, and it is why this design exists:**
- BeatAML, unconditioned (§7): instrument AUC **0.923**, beating 100% of matched random sets.
  NPM1-like n = **17** of 330 WT. Venetoclax raw median AUC diff **−90.8** (p 0.0014) →
  **p 0.060** with blasts adjusted → **p 0.073** with FLT3-ITD adjusted. Survival null
  (HR 0.87, p 0.73). Q3: 22/112 drugs q < 0.10, all same direction, four of top five FLT3-active.
- The confound that motivates conditioning: NPM1-like are **58.8% FLT3-ITD+ vs 14.4%** in WT-low
  (Fisher p = 5.9e-5).
- Stratified, underpowered, seen: ITD-negative NPM1-like n = 6, venetoclax p = 0.30 (same
  direction); ITD-positive n = 9, p = 0.065.
- FIMM external gate (§8): instrument AUC 0.794; both leads same sign but **+0.16 and +0.43 sDSS**
  at p 0.47 / 0.96 against a **+8.37** NPM1-mut positive control. Neither graduated.
- Set L behaves identically on the shared 14 genes (BeatAML AUC 0.923 either way).

**Consequence for this design, stated plainly:** the venetoclax and palbociclib leads have already
been tested twice and are NOT the primary here. Carrying them forward as primaries would be
testing the same hypothesis a third time on overlapping data. They appear in §4 only as a
pre-specified *replication check with no claim attached*.

## 1. The question this design can actually answer

Restricting to **FLT3-ITD-negative AML**, is the NPM1-like HOX state associated with
**overall survival**?

Why survival and not drug: survival is available in every third-party AML cohort at n in the
hundreds (BeatAML 458, TCGA-LAML ~151, AMLCG GSE146173 246, Leucegene 691), whereas ex vivo drug
data exists only in BeatAML and FIMM — both already used. Conditioning on ITD-negativity costs
~40–60% of the NPM1-like group in any single cohort, so only a design that can pool cohorts has a
chance. Survival is the only endpoint that pools.

**This flips the primary endpoint relative to pre-registration 1, where survival was Q2 and null
(HR 0.87 [0.41–1.86], p 0.73, unconditioned, n = 17 NPM1-like).** That null is the honest prior
for this design. It is registered here as such: **the expected outcome is no effect**, and the
purpose of the design is to bound it in the conditioned stratum rather than to find it.

## 2. Cohorts, and the rule for adding them

| cohort | n (RNA) | survival | NPM1 | FLT3-ITD | role |
|---|---|---|---|---|---|
| BeatAML2 | 458 (1pp, frozen) | yes | consensus | consensus + allelic ratio | discovery stratum |
| TCGA-LAML | ~151 | yes | to be verified | **to be verified — see §2.1** | second stratum |
| AMLCG GSE146173 | 246 | yes | yes | yes (ITD/TKD) | third stratum, resistance-enriched |
| Leucegene GSE232130 | 691 | none in GEO | — | — | NOT usable (no outcome) |

Cohorts enter as **independent strata in a fixed-effect meta-analysis**, never pooled at the
patient level (platform and population differ). A cohort is admissible only if it supplies, from
an open source: gene-level expression, NPM1 status, FLT3-ITD status, OS time and vital status.

### 2.1 Hard precondition: FLT3-ITD status

**If FLT3-ITD status cannot be obtained openly for a cohort, that cohort is excluded — it is not
an ITD-negative design without an ITD call.** ITDs are missed by standard MAF pipelines, so this
is a real risk for TCGA-LAML and is being verified before this document is finalised. If TCGA-LAML
fails this test, the design proceeds with BeatAML + AMLCG only (two strata) and that reduction is
recorded here rather than worked around.

## 3. Frozen definitions

- **Score**: Set L, **the shared 14 protein-coding genes** fixed in §8 (HOXA3/4/5/6/7/9/10,
  HOXB2/3/4/5/6, MEIS1, PBX3) — chosen because it is the set that exists on every platform,
  and it scores identically to the 17-gene version in BeatAML (AUC 0.923 both ways). Per-sample
  score = mean within-sample percentile rank, computed **within each cohort separately**.
- **Stratum**: FLT3-ITD-negative patients only. NPM1-mutant patients are used **only** to set the
  threshold and to verify the instrument; they are excluded from the outcome comparison.
- **NPM1-like**: ITD-negative, NPM1-WT patient whose score ≥ the 25th percentile of that cohort's
  **NPM1-mutant** score distribution (identical rule to pre-registration 1, computed within
  cohort). Comparator: all remaining ITD-negative NPM1-WT.
- **Primary model**: Cox `Surv(OS, dead) ~ NPM1like + age`, fitted **per cohort**, combined by
  inverse-variance fixed-effect meta-analysis on log HR. Heterogeneity reported (I², Q).
- **Primary result is the meta-analytic HR with its CI.** A null is reported as an interval, not
  as "p > 0.05".

## 4. Pre-specified secondary analyses, each with its claim status fixed now

| analysis | claim status |
|---|---|
| Continuous score (no dichotomy), same Cox model | co-primary sensitivity — if it disagrees with the dichotomised result, the dichotomy is reported as an artifact |
| ELN2017-adjusted (BeatAML/AMLCG where available) | sensitivity only; ELN partly encodes genotype |
| Blast-% adjusted | sensitivity; the composition analogue |
| ITD-**positive** stratum, same model | descriptive contrast, no claim |
| Venetoclax / palbociclib in ITD-negative BeatAML | **already-tested leads; reported for completeness, cannot be called a finding here whatever the result** |

## 5. Decision rules, fixed before running

| gate | rule | consequence |
|---|---|---|
| **instrument, per cohort** | Set-14 score separates NPM1-mut from WT at AUC ≥ 0.75 in that cohort | below 0.75 → cohort contributes to no outcome analysis, reported as instrument failure |
| **group size** | ≥ 10 NPM1-like in a cohort | below 10 → cohort reported but excluded from the meta-analysis |
| **minimum evidence** | ≥ 2 admissible cohorts | 1 cohort → the whole design is reported as not evaluable; no HR is promoted |
| **confound table** | ITD-negative NPM1-like vs comparator on age, blast %, TP53, RUNX1, ASXL1, sex, and secondary/de-novo status, printed BEFORE any survival model runs | imbalances are a mandatory caveat, not a stop |
| **equivalence reporting** | if the meta HR CI includes 1, report the CI bounds as the exclusion statement | prevents "absence of evidence" being written as "evidence of absence" |
| **no redefinition** | the threshold, the score, the model and the strata may not be changed after any survival output is seen | any deviation must be appended to this file, dated, and labelled post hoc |

## 6. Wording ceiling

Best case licenses: *"in FLT3-ITD-negative NPM1-WT AML, an NPM1-like HOX expression state is
associated with survival comparable to / different from NPM1-WT HOX-low disease, HR x.xx
[CI], meta-analysed across k cohorts"* — an association, in a defined stratum, with its CI.
It does NOT license: mechanism; the word "phenocopy" without the CI beside it; any statement
that the state is independent of unmeasured lesions (the DeltaNp73 objection from pre-registration
1 §5 applies unchanged); any drug claim.

**Most likely honest outcome, stated in advance:** a null with an interval — "an HR larger than
X is excluded in ITD-negative disease" — which, given §7's unconditioned null, would close the
discordance line on public bulk data and should be written as a bounded negative rather than
padded. This document exists so that outcome is publishable as a result rather than as a failure.

## 7. Outputs

`results/tables/13_bulk_anchor/`: `itdneg_cohort_table.csv` (per cohort: n, instrument AUC, group
sizes), `itdneg_confounds.csv`, `itdneg_survival.csv` (per-cohort and meta HR/CI), every row
carrying `kind_of_evidence`. Implementing scripts: `13.5_tcga_intake.R`, `13.6_amlcg_intake.R`,
`13.7_itdneg_meta.R` — written only after this document is approved and committed.
