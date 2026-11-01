# PRE-REGISTRATION A1-G01 — how much of the T-cell AML–healthy difference is T-subset composition, and how far does it reach beyond T cells?

WP-A1 gap G01 (composition versus within-subtype state), with the four items that A1-BATCH section 12.5 carried forward:
the processing-delay sensitivity, the comparison with B, NK and monocytes, the Control4 sex coding, and a sensitivity
without the three expression-female healthy donors. G09 (power of the S8 structural replication) is in section 5.

**Status: DRAFT v1, 2026-11-01 23:30. Not approved, not tagged. Nothing declared below has been computed.**
Order (user decisions 2026-10-28 and 2026-11-01): A1-G02 (PASS) → A1-BATCH (PASS with qualifications) → **this** →
G04(a). Handoff v2.0 §17 S15. Decision lines in section 6 are proposals; the user fixes them before the tag
`prereg-a1-g01-v1`.

## 0. Question, scope, and what each outcome lets the paper say

S5 found a T-cell imprint in 15 adult AML vs 10 healthy donors of GSE185381. A1-G02 showed it out of sample
(cross-fitted gap 1.378, i.e. 1.98 out-of-fold healthy-SD; full-pipeline permutation p 0.008). A1-BATCH found no
processing factor measurable inside GSE185381 that removes it, with three qualifications: the same direction appears in
B and NK cells of the same donors (and partly in monocytes), so it is not T-specific; it is aligned with a
processing-delay program; and the part that goes with the Control1-5 pool is tied to the three expression-female healthy
donors.

S5 adjusted for one composition covariate (the CD8 share of kept T cells). AML marrow T cells are known to carry fewer
naive and more effector-memory cells. **G01 asks how much of the gap is the mix of T subtypes (composition) and how much
is a shift inside the subtypes (state).**

Two estimates answer it from opposite sides:
- **S (primary)**: the gap measured inside BoneMarrowMap subtype labels, as a fraction of the whole-T gap on the same
  samples. A pure state shift gives S ≈ 1; a pure composition shift gives S ≈ 0.
- **C (co-reported)**: the gap that the observed subtype mix alone would produce if every subtype kept the pooled
  healthy expression profile, as a fraction of the whole-T gap. Pure composition gives C ≈ 1; pure state gives C ≈ 0.

What each outcome lets the paper say (wording fixed in section 6):
- **"within-state"**: the difference is mostly a shift inside T subtypes. Skeletons A/B of the assessment stay open.
- **"mixed"**: part composition, part state; both reported with their sizes. Skeleton A with that wording.
- **"composition"**: the difference is mainly the subtype mix. The paper becomes the composition-first skeleton C.

What G01 cannot do, under any outcome:
- **It cannot separate disease from sample source or handling.** A handling effect on the whole sample would also
  appear inside every subtype, so "within-state" does not mean "caused by AML". Only G04(a) (Chen2023, one protocol for
  both arms; 5 vs 4, so a pass is informative and a fail is not) addresses that.
- "State" means "inside a BoneMarrowMap fine label". A composition shift at a finer level than the labels (for example
  TEMRA-like cells inside "CD8 Effector Memory 2") is counted as state here.
- Label purity in AML is not verified: if AML cells labelled "CD4 Naive" are less naive than healthy ones, S is biased
  upwards and C downwards. The protein (ADT) purity check is deferred (assessment G01 option b).
- Adjusting for composition also removes any AML effect that acts through composition, so "beyond composition" is
  conservative and "composition" does not mean "unrelated to AML".

## 1. Already known before this draft (disclosed)

Seen, from earlier steps:
- Whole-T: gap_oof 1.378164 (A1-G02 R2), per-sample out-of-fold scores (`ws_a1_g02_oof_scores.csv`) and per-fold gene
  lists (`ws_a1_g02_fold_genes.csv`, 25 leave-one-out folds, 32 ambient-prone genes removed).
- Composition-related, all in-sample or descriptive: S5 adjusted for the CD8 share; the T score vs the CD8 share is 0.17
  (AML) / 0.22 (healthy) (not in a table); GSE227903 diagnosis score vs CD8 share 0.86 (n 7) / 0.83 (n 6) (memory only);
  in Oetjen the S9 age slope goes from 0.102 to -0.047 with the naive share in the model; refits inside CD4 Naive (12 AML
  + 8 healthy, 85 hits, sign concordance 1.00) and CD8 EM2 (12 + 9, 7 hits, 0.996), not cross-fitted; canonical subset
  markers (CCR7, TCF7, LEF1, GZMB, GZMH, NKG7, PRF1, GZMK) are not hits; the down arm holds naive/resting-associated
  genes (LINC02446, TMIGD2, RGS10, LDHB, NOSIP, SATB1-AS1, TXK, CAMK4) and TRDC (log2FC -1.86); the ambient filter
  removed SELL, IL7R, FHIT, RCAN3, RUNX3, CD5 and MAF from the 402.
- The composition table was read for eligibility counts only (section 2.2). The G01 verifier (2026-10-28) noted that
  naive share and status "barely overlap"; no share value is recorded.
- A1-BATCH section 12.4 (post-hoc, in-sample): P370 retention in B 0.69 and NK 0.93, sign concordance 0.97 / 0.95, above
  all 60 random relabellings; monocytes R 0.57, sign concordance 0.87 (11 AML + 10 healthy); genome-wide T-B logFC
  correlation 0.60 lies inside donor-grouping nulls, T-NK 0.70 at their top. Delay: the imprint overlaps the Savage 2021
  delay genes in the same direction; the cross-fitted gap without delay genes is 1.318.
- A1-BATCH: the 9 of 40 matched drops that remove Control1, Control3 and Control4 have R 0.75-0.92, the other 31
  0.89-1.18.

Not computed (blind at this draft): any within-label score or gap; any composition counterfactual; any held-out B, NK
or monocyte score; any naive-share coefficient; S8's T1/T2 among healthy donors or against composition-matched random
sets; gamma-delta/MAIT flags; any GSE227903 post-treatment or relapse expression (S10 stays blind); Chen2023 and the
GSE185381 children (G04).

## 2. Inputs and frozen definitions (fixed before any expression is read)

### 2.1 Samples, folds, engine
- Samples: the 25 of the S5 T design (`ws_s5_scores.csv`, compartment T; 15 AML, 10 healthy), file order.
- Folds: the 25 leave-one-out folds of A1-G02 R2. Fold i's gene lists U_i (up) and D_i (down) are read from
  `ws_a1_g02_fold_genes.csv` and **not reselected** for the primary analysis (selected on the other 24 whole-T
  pseudobulks; design MAIN = status + hb + gran + sex + subtype_share, metadata sex; edgeR QL; FDR < 0.10 and
  |log2FC| > 0.5; 32 ambient-prone genes removed). The 402-analogue (ambient-prone genes kept) refits the 25 folds with
  the same code; known answer K0b checks the refit against the stored lists.
- Log scale: `logz()` of 37_a1_batch.R (blacklist; CPM > 1 in at least half of the columns; TMM; log-CPM, prior 1),
  applied to the matrix being scored.

### 2.2 Subtype labels and eligibility
- Labels: the nine BoneMarrowMap fine labels mapped to T in `s4_cell_groups.tsv`. Counts: `PB$fine[[label]]` of
  `ws_pb_GSE185381.rds` (raw counts summed over kept cells); kept-cell numbers from `PB$cells`.
- Eligible: at least 20 kept cells of the label in the sample (WS_S4_MIN_CELLS).
- Primary labels: the four with at least 17 eligible samples (verifier's counts, rechecked as K1): CD4 Naive (12 AML +
  8 healthy), CD8 Effector Memory 2 (12 + 9), CD8 Effector Memory 1 (11 + 7), CD4 Central Memory (10 + 7). The other
  five labels enter C (section 3.2) but not S.

### 2.3 Score units
- **u-score (primary unit)**: for sample i scored with fold i's genes in a log-CPM matrix X,
  u_i = mean over g in U_i of (X[g, i] − m_g) − mean over g in D_i of (X[g, i] − m_g), where m_g is the mean of X[g, ]
  over the other samples of the same matrix. Genes outside X's universe are skipped (coverage reported).
- **z-score (A1-G02's unit, reported)**: the same with (X[g, i] − m_g) / s_g, s_g the SD over the other samples.
- Why the primary unit is u, not z: a label pseudobulk has fewer cells than the whole-T pseudobulk, so its
  between-sample SD carries more counting noise; dividing by it shrinks the within-label z-gap even when the per-cell
  shift is identical. The u-score is a plain log-fold difference and has no such scaling. Its remaining bias (the prior
  count weighs more in small libraries) is measured by the planted known answers K3.

## 3. Primary and co-reported estimates

### 3.1 S — the within-label fraction (primary)
For each primary label k, with E_k its eligible samples:
1. X_k = logz(label-k counts of E_k); X_T,k = logz(whole-T counts of the same E_k).
2. For every i in E_k: s_ik = u-score of i in X_k; t_ik = u-score of i in X_T,k (both centred on E_k without i).
3. gap_k = mean of s_ik over AML − mean over healthy; gapT_k the same with t_ik. S_k = gap_k / gapT_k.

Pooled: **S = (sum over the four labels of gap_k) / (sum of gapT_k)**. CI: 2,000 bootstrap draws resampling AML and
healthy samples separately (a drawn sample enters every label where it is eligible); per-sample held-out scores and fold
genes held fixed, so the CI is conditional on gene selection. 95% percentile interval. Also reported: S_k per label, the
z-score versions, gene coverage per label, and the number of folds with an empty half.

### 3.2 C — the composition counterfactual (co-reported)
For each fold i:
1. Reference profile per label k (all nine T labels): r_k = the summed label-k counts of the healthy samples other than
   i, divided by their summed kept label-k cells (a per-cell profile pooled over healthy donors, so labels that are
   rare in single donors still have one).
2. Synthetic whole-T counts for every sample s: y_s = sum over k of n_ks · r_k, where n_ks is s's kept cells of label k.
3. Y = logz of the 25 synthetic columns; c_i = u-score of i in Y (centred on the other 24).

gap_C = mean of c_i over AML − mean over healthy; gapT = the u-score whole-T gap over all 25 (fold genes as above).
**C = gap_C / gapT.** Same bootstrap. C uses only each sample's label counts, never its expression, so it answers: how
much of the gap would the observed subtype mix produce at healthy states? Also reported: S + C (near 1 if the two parts
add up).

<!-- part 2 (sections 4-9: sensitivities and carried-forward items, gamma-delta/MAIT, G09, known answers, decision
rules, code plan, user decisions) follows in the next draft step -->
