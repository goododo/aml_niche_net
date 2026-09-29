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

---

## 8. Results (2026-09-29, `13.5_tcga_intake.R`, `13.7_itdneg_meta.R`, `13.8_reducibility_interrogation.R`)

### 8.1 §2.1 precondition: PASSED, by derivation rather than by a curated column

No open source carries a curated FLT3-ITD boolean. GDC's open masked somatic MAFs are useless for
it (6 FLT3 rows in all 153 files, every one a missense SNP; zero insertions; NPM1 sensitivity ~20%).
The working channel is cBioPortal's `laml_tcga_pub` MAF (the WashU curated calls), where exon-14
juxtamembrane in-frame insertions — many annotated literally as `p.S584_D600dup`, `p.Y597_K602dup` —
are the ITDs. Derivation rule fixed in `13.5`: ITD = FLT3 ∧ In_Frame; TKD = FLT3 missense at protein
position 820–860. **All four known answers reproduced independently: NPM1 54, ITD 30, TKD 17,
ITD∩TKD 0**, with 45/46 concordance against the orthogonal RT-PCR NPMc channel. TCGA-LAML is
admissible; AMLCG was not needed.

### 8.2 The dichotomised design is NOT EVALUABLE — and that failure is the result

| cohort | instrument AUC (gate ≥0.75) | ITD-neg mut / WT | ITD-neg internal AUC | NPM1-like | gate ≥10 |
|---|---|---|---|---|---|
| TCGA-LAML | **0.976** PASS | 33 / 108 | 0.982 | **1** (0.9%) | FAIL |
| BeatAML | **0.932** PASS | 68 / 260 | 0.941 | **7** (2.7%) | FAIL |

Zero cohorts pass both gates, so per §5 no HR from the dichotomy is promoted. What the failing gate
measured is worth stating as a finding in its own right: **of 368 FLT3-ITD-negative NPM1-WT
patients across two cohorts, 8 (2.2%) reach the NPM1-mutant 25th percentile and 0 (0.0%, exact 95%
CI 0–1.00%) reach the mutant median HOX score.** In ITD-negative AML the HOX/MEIS1 program is
almost perfectly determined by NPM1 genotype (within-stratum AUC 0.98 / 0.94). **The discordance
hypothesis has no substrate there** — which also explains pre-registration 1 §7: 58.8% of that
cohort's "NPM1-like" cases were FLT3-ITD+, and in BeatAML ITD+ NPM1-WT patients do carry elevated
HOX (median 0.679 vs 0.453, p = 5e-7). The phenocopy group was largely an ITD group.

### 8.3 The co-primary continuous analysis: a real, replicated association that is fully reducible

Registered co-primary (Cox OS ~ z(HOX) + age, ITD-negative NPM1-WT only, inverse-variance fixed
effect):

| arm | META HR per 1 SD [95% CI] | p | I² | verdict |
|---|---|---|---|---|
| **registered** | **1.253 [1.091–1.439]** | 0.0014 | 0% | positive, both cohorts same direction |
| exclude HOX-fusion patients (KMT2A/MLLT3/NUP98/MECOM) | 1.271 [1.096–1.473] | 0.0015 | 0% | survives — not a fusion detector |
| adjust for blast % | 1.272 [1.093–1.480] | 0.0019 | 0% | survives |
| initial-diagnosis specimens only | 1.211 [1.037–1.414] | 0.015 | 0% | survives, weaker |
| **adjust for cytogenetic/molecular risk class** | **0.917 [0.763–1.101]** | 0.35 | 0% | **ABOLISHED** |
| exclude fusions AND adjust risk | 0.906 [0.747–1.097] | 0.31 | 0% | abolished |

The abolition is due to the adjustment, not to the subset it requires: on the identical
risk-classified subsets the unadjusted HRs are 1.288 (TCGA, n=107) and 1.171 (BeatAML, n=185),
versus 0.974 and 0.887 once risk enters. The mechanism is visible in the score itself — it tracks
risk class monotonically in both cohorts (TCGA Good 0.283 / Intermediate 0.515 / Poor 0.542;
BeatAML Favorable 0.261 / Intermediate 0.511 / Adverse 0.522), i.e. **HOX-low marks
favorable-risk core-binding-factor AML, which is textbook.** Tertile medians are monotone
(TCGA median OS 30.6 / 14.6 / 8.6 months), and 5 within-cohort score permutations give meta HRs
0.94–1.09 (p 0.22–0.93), so the machinery is not manufacturing the effect.

### 8.4 Conclusion at the registered wording ceiling

**Licensed:** in FLT3-ITD-negative NPM1-WT AML, a higher single-cell-derived HOX/MEIS1 score is
associated with shorter overall survival (HR 1.25 per SD [1.09–1.44], two cohorts, I² = 0%), but
**none of that information is independent of established cytogenetic/molecular risk
classification** (risk-adjusted HR 0.92 [0.76–1.10]); the score behaves as a favorable-risk
(CBF-AML) proxy.
**Not licensed:** any claim of new prognostic information, any mechanism, any drug statement.

**What survives the whole L2 line as a positive:** the instrument. The HOX program derived from
36 single-cell samples reads NPM1 genotype in **three independent bulk cohorts** — BeatAML AUC
**0.923**, TCGA-LAML **0.976**, FIMM **0.794** — with a matched-random-gene-set specificity control
(100th percentile in BeatAML). Combined with §8.2's rarity bound, the coherent unit of work is:
*a validated single-cell→bulk transfer instrument, plus the quantitative demonstration that the
genotype it reads has essentially no transcriptional phenocopy in ITD-negative AML, and that its
prognostic signal is not independent of risk class.* That is a bounded negative with a working
positive control, which §6 anticipated and which is reportable as a result.

**Not attempted further on public bulk data.** AMLCG GSE146173 would add a third stratum to an
already homogeneous (I² = 0%) null-after-adjustment; it cannot change the reducibility verdict.
