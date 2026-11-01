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

## 4. Sensitivities (each reports S and C; none changes the wording by itself, section 7)

- **4.1 The 402-analogue**: the 25 folds refitted without removing the 32 ambient-prone genes (they include SELL, IL7R
  and other subset markers, so composition should weigh more here).
- **4.2 Processing delay (carried forward from A1-BATCH 12.5)**: fold lists minus the Savage 2021 18-h delay genes
  (295 up / 497 down; `savage2021_T_delay18h_consistent.csv`, md5 5ee444a9284ae52b7afab1c8aeeced2c). Known answer: the
  whole-T z gap with these genes removed is 1.3184 (A1-BATCH 12.4, same construction).
- **4.3 Control4 sex coding**: the 25 folds refitted with Control4 coded female (its cells are female: XIST 7,797 vs
  RPS4Y1 56 UMI); A1-G02 section 7 made this the G06 sex sensitivity.
- **4.4 Without the three expression-female healthy donors** (Control1, Control3, Control4): 22 folds refitted on 15
  AML + 7 healthy (each fold trains on 21); S and C on the 22. Eligible counts per label are reported.
- **4.5 Leave one label out**: S recomputed four times, each time without one of the four primary labels.

## 5. Further declared readouts

### 5.1 Beyond T cells, cross-fitted (new; carried forward from A1-BATCH 12.5)
A1-BATCH 12.4 measured the sharing in-sample: the 370 genes were chosen on the same donors whose B and NK cells were then
read, so donor-level noise picked up by the selection can also appear in B and NK. Here each held-out donor's B, NK and
monocyte pseudobulk (`PB$raw`, S5's compartment designs; only donors among the 25, since only they have a fold) is
scored with the T genes of the fold that left that donor out (u and z, centred on the other eligible donors of that
compartment). Reported per compartment: n, the held-out gap, and its ratio to the whole-T held-out gap on the same
donors. Reading: a ratio near 1 means the shift is as large in that compartment as in T cells; near 0, confined to T.
Caveat stated with the result: monocytes in AML samples are often clonal, so a monocyte shift may be leukaemic biology
rather than a shared environment or handling effect. The in-sample numbers of A1-BATCH 12.4 are re-derived in the same
script and written to a durable table (known answers: R_B 0.692, R_NK 0.928, monocytes 0.568).

### 5.2 Gamma-delta and MAIT cells (verifier P3)
BoneMarrowMap has no gamma-delta or MAIT label, and Lasry 2023 reported MAIT loss in this dataset. From the raw counts of
the kept T cells (04_annotated objects): gamma-delta flag = TRDC ≥ 1 and TRGC1 + TRGC2 ≥ 1 UMI; MAIT flag = TRAV1-2 ≥ 1
and SLC4A10 ≥ 1 UMI. Detection is incomplete at this depth, so the flags are lower bounds. Reported: flagged shares per
sample by status (sanity range: healthy median gamma-delta share between 0.5% and 15% of T cells); the whole-T z gap
after removing flagged cells from the pseudobulks (fold genes fixed); the gap with TRDC removed from the fold lists.

### 5.3 G09: what S8's structural replication can show (verifier P4)
S8's T2 (up and down halves move oppositely across patients) and T1 (imprint genes co-vary within a patient's T cells)
would pass for any gene set aligned with a T-differentiation axis, because subtype heterogeneity alone produces both.
Rerun, with S8's code and thresholds:
- (a) **among healthy donors**: GSE185381 healthy (10) and Oetjen (20). If T1 or T2 passes there as often as in AML,
  the manuscript may not cite it as evidence of an AML-specific program.
- (b) **against random sets matched also on a composition loading**: each gene's log2 difference between pooled naive
  (CD4 + CD8 Naive) and effector-memory (CD8 EM1 + EM2, CD4 EM) T pseudobulks of the Oetjen donors, in terciles, on top
  of S8's expression and detection bins. Cohorts GSE289435 and GSE227903 (manifest role AML_Dx only: the 9 diagnosis
  samples are selected by metadata before any expression is loaded, so S10 stays blind).
Section 9 item 3 asks whether G09 stays here or becomes its own short amendment after G01.

### 5.4 Status plus composition regression (verifier P5; secondary)
The A1-G02 out-of-fold z-scores regressed on status + naive share (CD4 + CD8 Naive / kept T) + CD8 EM share (EM1 + EM2 /
kept T), **only if an overlap check passes**: at least 3 AML and 3 healthy samples inside the common range of naive
share (from the larger of the two group minima to the smaller of the two group maxima). Otherwise reported as "not
estimable: no overlap", and nothing else is fitted.

### 5.5 Reported for every estimate
Both gene sets (370-analogue primary, 402-analogue in 4.1); n per label and compartment; gene coverage; the number of
folds with an empty half.

## 6. Known answers (run first; no readout is printed until all pass)

- K0: the whole-T z gap from the stored fold lists reproduces 1.378164 (|difference| < 1e-6).
- K0b: the refit code (used by 4.1, 4.3, 4.4) run with the primary settings reproduces the 25 stored fold lists exactly.
- K1: eligible counts per primary label equal 12 + 8, 12 + 9, 11 + 7, 10 + 7 (section 2.2).
- K2: for every sample, the nine T-label count vectors sum exactly to the whole-T pseudobulk.
- K3, planted on the 10 healthy donors only (5 drawn as pseudo-AML, 20 draws, fixed seed), every T label's counts
  rebuilt into whole-T as their sum:
  - K3a pure state: every T label of the pseudo-AML donors gets the S5 main-fit log2FC on the P370 genes (counts times
    2^log2FC), scored with P370 as the gene list. Expected S near 1 and C near 0. Pass: median S ≥ 0.8 and median
    |C| ≤ 0.2. The median S here is also the "pure state" value that real S is read against (it measures the shrinkage
    from small label libraries).
  - K3b pure composition: no expression change; the pseudo-AML donors' naive labels (CD4 Naive, CD8 Naive) are halved and
    their effector-memory labels (CD8 EM1, EM2, CD4 EM) doubled, counts and cell numbers alike; scored with the 100 genes
    most higher and the 100 most lower in pooled healthy naive versus effector-memory profiles. Expected S near 0 and C
    near 1. Pass: median |S| ≤ 0.2 and median C ≥ 0.8.
- K4: Control4's whole-T counts are female by XIST/RPS4Y1 (as recorded); Control1 and Control3 likewise.
- K5: the Savage delay file has md5 5ee444a9284ae52b7afab1c8aeeced2c; 4.2's whole-T z gap is 1.3184 (± 1e-4).
- K6: the GSE227903 sample list has exactly the 9 AML_Dx samples, and no other column is loaded.
- K7: the 5.1 in-sample re-derivation gives R_B 0.692, R_NK 0.928, monocytes 0.568 (± 0.001).

## 7. Decision rules (proposed; the user fixes the values before the tag)

- **Not decisive**: any known answer fails; or gapT_k ≤ 0 for two or more of the four labels.
- **"within-state"**: S ≥ 1/2, the lower bound of S's 95% CI > 0, and C ≤ 1/2.
- **"composition"**: S < 1/3, or C ≥ 2/3.
- **"mixed"**: anything else.
- The lines echo the assessment (1/2, 1/3) and the 1/3 line of GATE_A1(4) and S9. C acts as a veto so that a high S
  that the subtype mix could also produce is not called "within-state".
- Sensitivities (section 4) and readouts (section 5) never change the wording. If one of them would put S or C in a
  different category, the manuscript reports it next to the primary.
- Under every outcome the manuscript states the three limits of section 0: no separation of disease from source or
  handling; "state" means inside a BoneMarrowMap label; label purity in AML unverified.

## 8. Code, outputs, compute (plan; scripts are written after the tag)

- One script, `scripts/ws/39_a1_g01.R`, modes: `--known` (K0-K7), `--primary` (S, C, bootstrap), `--sens` (4.1-4.5;
  refits 25 + 25 + 22 folds), `--share` (5.1), `--regress` (5.4), `--combine` (tables, GATES block A1-G01). A second
  script, `scripts/ws/40_a1_g01_cells.R`, for the per-cell steps: 5.2 (reads the 25 annotated objects) and, if kept, G09
  (5.3: per-cell T1 in healthy donors and the two AML cohorts). sbatch wrappers as for A1-BATCH.
- Outputs: `results/tables/ws/ws_a1_g01_{checks,labels,primary,boot,sens,share,gdmait,regress,g09}.csv`.
- Smoke runs first, on synthetic counts with real metadata; smoke tables deleted after the real run reproduces.
- Compute: every mode except 40 runs in under an hour on one node; 40 needs about 64 GB for the annotated objects.
- No new package. No download.

## 9. For the user to decide before the tag

1. Decision lines of section 7 (recommended as written: 1/2 and 1/3 for S; C ≤ 1/2 and ≥ 2/3).
2. C as a veto on "within-state" (recommended) or S alone.
3. G09 (5.3) inside this amendment, as the assessment ordered, or as its own short amendment right after G01
   (recommended: separate, so this one stays small and finishes in one or two sessions).
4. The gamma-delta/MAIT arm (5.2) now (recommended: yes; one job of about an hour).

## 10. Review history

- 2026-11-01: DRAFT v1 (sections 0-3 at ff53d2b; 4-10 in the next commit). Short reviews follow.
