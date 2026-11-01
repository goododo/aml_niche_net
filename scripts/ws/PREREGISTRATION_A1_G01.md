# PRE-REGISTRATION A1-G01 — how much of the T-cell AML–healthy difference is T-subset composition, and how far does it reach beyond T cells?

WP-A1 gap G01 (composition versus within-subtype state), with the four items that A1-BATCH section 12.5 carried forward:
the processing-delay sensitivity, the comparison with B, NK and monocytes, the Control4 sex coding, and a sensitivity
without the three expression-female healthy donors. G09 (what S8's structural replication can show) is in section 5.4.

**Status: DRAFT v2, 2026-11-01 23:40, after short review 1 (section 10; v1 at bd29ce6). Not approved, not tagged. Nothing
declared below has been computed.** Order (user decisions 2026-10-28 and 2026-11-01): A1-G02 (PASS) → A1-BATCH (PASS
with qualifications) → **this** → G04(a). Handoff v2.0 §17 S15. Decision lines in section 7 are proposals; the user
fixes them before the tag `prereg-a1-g01-v1`.

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

- **C (primary): the composition counterfactual.** The gap that the observed subtype mix alone would produce if every
  subtype kept the pooled healthy expression profile, as a fraction of the whole-T gap. It uses all nine T labels and
  all 25 samples. Pure composition gives C ≈ 1; a pure state shift gives C ≈ 0 wherever it sits. 1 − C is the state
  share (composition is evaluated at healthy states, the state part at the AML mix).
- **S_w (required check): the shift inside the four large labels.** The gap measured inside each label, as a fraction of
  the whole-T gap on the same samples and genes, weighted by the label's share of T-cell UMIs and summed. It counts only
  the state shift inside these four labels, so it is a lower bound on the state share: a shift confined to the other
  five labels gives S_w ≈ 0, and small label libraries shrink it further (measured by K3a). It shows directly, without
  a model, that a shift exists inside subtypes.

What each outcome lets the paper say (rules in section 7):
- **"within-state"**: mostly a shift inside T subtypes, seen inside the main subtypes. Skeletons A/B stay open.
- **"mixed"**: part composition, part state; both sizes reported. Skeleton A with that wording.
- **"composition"**: the difference is mainly the subtype mix. The paper becomes the composition-first skeleton C.
- **"not composition; state not located"**: composition explains little, but the shift was not resolved inside the
  four main labels. Skeleton A with that wording.

What G01 cannot do, under any outcome:
- **It cannot separate disease from sample source or handling.** A handling effect on the whole sample would also
  appear inside every subtype, so "within-state" does not mean "caused by AML". Only G04(a) (Chen2023, one protocol for
  both arms; 5 vs 4, so a pass is informative and a fail is not) addresses that.
- "State" means "inside a BoneMarrowMap fine label". A composition shift at a finer level than the labels (for example
  TEMRA-like cells inside "CD8 Effector Memory 2") is counted as state here; section 5.5 asks whether the
  within-label shift runs along the naive-to-effector axis.
- Label purity in AML is not verified: if AML cells labelled "CD4 Naive" are less naive than healthy ones, S_w is
  biased upwards and C downwards. The protein (ADT) purity check is deferred (assessment G01 option b).
- Adjusting for composition also removes any AML effect that acts through composition, so "beyond composition" is
  conservative and "composition" does not mean "unrelated to AML".

## 1. Already known before this draft (disclosed)

Seen, from earlier steps:
- Whole-T: gap_oof 1.378164 (A1-G02 R2), per-sample out-of-fold scores (`ws_a1_g02_oof_scores.csv`) and per-fold gene
  lists (`ws_a1_g02_fold_genes.csv`, 25 leave-one-out folds, 32 ambient-prone genes removed; the smallest lists are
  11 up / 18 down with AML4340 left out and 21 up / 65 down with Control3 left out).
- Composition-related, all in-sample or descriptive: S5 adjusted for the CD8 share; the T score vs the CD8 share is 0.168
  (AML) / 0.224 (healthy), and in GSE227903 diagnosis samples 0.857 (n 7) / 0.829 (n 6) (`ws_a1_g02_descriptive.csv`);
  in Oetjen the S9 age slope goes from 0.102 to -0.047 with the naive share in the model; refits inside CD4 Naive (12 AML
  + 8 healthy, 85 hits, sign concordance 1.00) and CD8 EM2 (12 + 9, 7 hits, 0.996), not cross-fitted; canonical subset
  markers (CCR7, TCF7, LEF1, GZMB, GZMH, NKG7, PRF1, GZMK) are not hits; the down arm holds naive/resting-associated
  genes (LINC02446, TMIGD2, RGS10, LDHB, NOSIP, SATB1-AS1, TXK, CAMK4) and TRDC (log2FC -1.86); the ambient filter
  removed SELL, IL7R, FHIT, RCAN3, RUNX3, CD5 and MAF from the 402.
- Composition table: read for eligibility counts (section 2.2). The G01 verifier (2026-10-28) also recorded pooled
  shares of the 14,609 kept T cells (both arms together: CD4 Naive about 17%, CD8 EM2 about 39%) and that CD8 Naive
  reaches 20 kept cells in 1/15 AML vs 5/10 healthy samples (CD8 naive cells are scarcer in AML), and noted that naive
  share and status "barely overlap". No per-arm share value is recorded.
- A1-BATCH section 12.4 (post-hoc, in-sample): P370 retention in B 0.69 and NK 0.93, sign concordance 0.97 / 0.95, above
  all 60 random relabellings; monocytes R 0.57, sign concordance 0.87 (11 AML + 10 healthy); genome-wide T-B logFC
  correlation 0.60 lies inside donor-grouping nulls, T-NK 0.70 at their top. Delay: the imprint overlaps the Savage 2021
  delay genes in the same direction; the cross-fitted gap without delay genes is 1.3184.
- A1-BATCH: the 9 of 40 matched drops that remove Control1, Control3 and Control4 have R 0.75-0.92, the other 31
  0.89-1.18.

Not computed (blind at this draft): any within-label score, gap or gene-level shift; any composition counterfactual;
any held-out B, NK or monocyte score; any per-arm subtype share or naive-share coefficient; S8's T1/T2 among healthy
donors or against composition-matched random sets; gamma-delta/MAIT flags; any GSE227903 post-treatment or relapse
expression (S10 stays blind); Chen2023 and the GSE185381 children (G04).

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

### 2.2 Subtype labels, eligibility, weights
- Labels: the nine BoneMarrowMap fine labels mapped to T in `s4_cell_groups.tsv`. Counts: `PB$fine[[label]]` of
  `ws_pb_GSE185381.rds` (raw counts summed over kept cells); kept-cell numbers from `PB$cells`.
- Eligible: at least 20 kept cells of the label in the sample (WS_S4_MIN_CELLS).
- Primary labels: the four with at least 17 eligible samples (verifier's counts, rechecked as K1): CD4 Naive (12 AML +
  8 healthy), CD8 Effector Memory 2 (12 + 9), CD8 Effector Memory 1 (11 + 7), CD4 Central Memory (10 + 7). The other
  five labels enter C but not S_w.
- Weights: π_k = label k's share of whole-T UMIs, summed over the 25 samples (both arms together).

### 2.3 Score units
- **u-score (the unit of every ratio)**: for sample i scored with gene lists U and D in a log-CPM matrix X,
  u_i = mean over g in U of (X[g, i] − m_g) − mean over g in D of (X[g, i] − m_g), where m_g is the mean of X[g, ]
  over the other samples of the same matrix.
- **z-score (A1-G02's unit, reported)**: the same with (X[g, i] − m_g) / s_g, s_g the SD over the other samples.
- Why ratios use u, not z: a label pseudobulk has fewer cells than the whole-T pseudobulk, so its between-sample SD
  carries more counting noise; dividing by it shrinks the within-label z-gap even when the per-cell shift is identical.
  The u-score is a plain log-fold difference. Its remaining bias (the prior count weighs more in small libraries) is
  measured by the planted known answers K3.

### 2.4 One gene list per comparison
Every ratio scores its numerator and its denominator on the same genes: fold i's genes that lie in the universes of
both matrices. A sample (or sample-label pair) whose matched list has fewer than 5 genes in either half is dropped from
that comparison, and the number dropped is reported. Reason: an empty half scores as 0 (`score_col`, 37_a1_batch.R),
and the smallest stored lists have 11 up / 18 down genes.

## 3. Primary estimate and required check

### 3.1 C — the composition counterfactual (primary)
For each fold i:
1. Reference profile per label k (all nine T labels): r_k = the summed label-k counts of the healthy samples other than
   i, divided by their summed kept label-k cells (a per-cell profile pooled over healthy donors, so labels that are
   rare in single donors still have one). A label with fewer than 50 pooled healthy kept cells in a fold takes the
   pooled profile of its group (CD4_T or CD8_T in `s4_cell_groups.tsv`); how often this happens is reported.
2. Synthetic whole-T counts for every sample s: y_s = sum over k of n_ks · r_k, where n_ks is s's kept cells of label k.
3. Y = logz of the 25 synthetic columns; X_T = logz of the 25 real whole-T columns. With fold i's matched list
   (section 2.4): c_i = u-score of i in Y, t_i = u-score of i in X_T (both centred on the other 24).

gap_C = mean of c_i over AML − mean over healthy; gapT = the same with t_i. **C = gap_C / gapT.** CI: 2,000 bootstrap
draws resampling AML and healthy samples separately, per-sample held-out scores and fold genes held fixed (so the CI
is conditional on gene selection); 95% percentile interval. C uses only each sample's label counts, never its
expression: it answers how much of the gap the observed subtype mix would produce at healthy states.

### 3.2 S_w — the shift inside the four large labels (required check)
For each primary label k, with E_k its eligible samples:
1. X_k = logz(label-k counts of E_k); X_T,k = logz(whole-T counts of the same E_k).
2. For every i in E_k, with fold i's list matched to X_k and X_T,k: s_ik = u-score of i in X_k; t_ik = u-score of i in
   X_T,k (both centred on E_k without i).
3. gap_k = mean of s_ik over AML − mean over healthy; gapT_k the same with t_ik; S_k = gap_k / gapT_k.

**S_w = sum over the four labels of π_k · S_k.** Same bootstrap (a drawn sample enters every label where it is
eligible). Also reported: each S_k, the z-score versions, gene coverage per label, and K3a's median S_w as the
pure-state benchmark (it combines the four labels' UMI share and the small-library shrinkage).

## 4. Sensitivities (each reports C and S_w; none changes the wording, section 7)

- **4.1 The 402-analogue**: the 25 folds refitted without removing the 32 ambient-prone genes (they include SELL, IL7R
  and other subset markers, so composition should weigh more here).
- **4.2 Processing delay (carried forward from A1-BATCH 12.5)**: fold lists minus the Savage 2021 18-h delay genes
  (295 up / 497 down; `savage2021_T_delay18h_consistent.csv`, md5 5ee444a9284ae52b7afab1c8aeeced2c).
- **4.3 Control4 sex coding**: the 25 folds refitted with Control4 coded female; A1-G02 section 7 made this the G06
  sex sensitivity.
- **4.4 Without the three expression-female healthy donors** (Control1, Control3, Control4): 22 folds refitted on 15
  AML + 7 healthy (each fold trains on 21); C and S_w on the 22. Eligible counts per label are reported.
- **4.5 Leave one label out**: S_w recomputed four times, each time without one of the four primary labels.

## 5. Further declared readouts

### 5.1 Beyond T cells, cross-fitted (new; carried forward from A1-BATCH 12.5)
A1-BATCH 12.4 measured the sharing in-sample: the 370 genes were chosen on the same donors whose B and NK cells were then
read, so donor-level noise picked up by the selection can also appear in B and NK. Here each held-out donor's B, NK and
monocyte pseudobulk (`PB$raw`, S5's compartment designs; only donors among the 25, since only they have a fold) is
scored with the T genes of the fold that left that donor out, matched (section 2.4) to that compartment's universe and
to the whole-T universe of the same donors. Reported per compartment: n, gene coverage, the held-out gap, and its ratio
to the whole-T held-out gap on the same donors and genes, with a bootstrap CI.
- Reading is one-sided: a ratio whose CI lower bound is above 0 shows the shift in that compartment out of sample. A
  ratio near 0 is **not** evidence that the shift is confined to T cells, because smaller libraries and the loss of
  T-specific genes both shrink it.
- Monocytes in AML samples are often clonal, so a monocyte shift may be leukaemic biology rather than a shared
  environment or handling effect; stated with the result.
- The in-sample numbers of A1-BATCH 12.4 are re-derived in the same script and written to a durable table (K7).

### 5.2 Gamma-delta and MAIT cells (verifier P3)
BoneMarrowMap has no gamma-delta or MAIT label, and Lasry 2023 reported MAIT loss in this dataset. From the raw counts of
the kept T cells (04_annotated objects): gamma-delta flag = TRDC ≥ 1 and TRGC1 + TRGC2 ≥ 1 UMI; MAIT flag = TRAV1-2 ≥ 1
and SLC4A10 ≥ 1 UMI. Detection is incomplete at this depth, so the flags are lower bounds. Reported: flagged shares per
sample by status (sanity range: healthy median gamma-delta share between 0.5% and 15% of T cells); the whole-T z gap
after removing flagged cells from the pseudobulks (fold genes fixed); the gap with TRDC removed from the fold lists.

### 5.3 Status plus composition regression (verifier P5; secondary)
The A1-G02 out-of-fold z-scores regressed on status + naive share (CD4 + CD8 Naive / kept T) + CD8 EM share (EM1 + EM2 /
kept T), **only if an overlap check passes**: at least 3 AML and 3 healthy samples inside the common range of naive
share (from the larger of the two group minima to the smaller of the two group maxima). Otherwise reported as "not
estimable: no overlap", and nothing else is fitted.

### 5.4 G09: what S8's structural replication can show (verifier P4)
S8's T2 (up and down halves move oppositely across patients) and T1 (imprint genes co-vary within a patient's T cells)
would pass for any gene set aligned with a T-differentiation axis, because subtype heterogeneity alone produces both.
Rerun with S8's code and thresholds:
- (a) **Among healthy donors.** Only Oetjen (20 donors) counts for the rule: GSE185381's healthy donors were in the gene
  selection, so their result is reported but not used. T1 "passes in healthy" if the imprint exceeds its random P95
  in at least 75% of evaluable Oetjen donors (S8's 6 of 8); T2 "passes in healthy" if rho < 0 and below the random P5
  across Oetjen donors. If either passes in healthy marrow, the manuscript may not cite it as evidence of an
  AML-specific program.
- (b) **Against random sets matched also on a composition loading**: each gene's log2 difference between pooled naive
  (CD4 + CD8 Naive) and effector-memory (CD8 EM1 + EM2, CD4 EM) T pseudobulks of the Oetjen donors, in terciles, on top
  of S8's expression and detection bins. Cohorts GSE289435 and GSE227903 (AML_Dx only, K6).
Section 9 item 3 asks whether G09 stays here or becomes its own short amendment after G01.

### 5.5 Does the within-label shift run along the naive-to-effector axis? (descriptive)
Within CD4 Naive and CD8 EM2 (the two labels with most eligible samples), the AML − healthy difference in mean log-CPM
for each gene of the union of the 25 fold lists, correlated (Spearman) with the Oetjen naive-vs-effector-memory loading
of 5.4(b). A clear positive correlation (the AML shift points towards effector memory) means part of the "state" may be
composition below the label resolution; a correlation near 0 means the within-label shift is not that axis. In-sample
and descriptive; no threshold.

### 5.6 Reported for every estimate
Both gene sets (370-analogue primary, 402-analogue in 4.1); n per label and compartment; gene coverage; samples dropped
by the 5-gene rule; the share of kept T cells that the four primary labels cover, by status.

## 6. Known answers (run first; no readout is printed until all pass)

- K0: the whole-T z gap from the stored fold lists reproduces 1.378164 (|difference| < 1e-6).
- K0b: the refit code (used by 4.1, 4.3, 4.4) run with the primary settings reproduces the 25 stored fold lists exactly.
- K1: eligible counts per primary label equal 12 + 8, 12 + 9, 11 + 7, 10 + 7 (section 2.2).
- K2: for every sample, the nine T-label count vectors sum exactly to the whole-T pseudobulk.
- K3, planted on the 10 healthy donors only: 5 drawn as pseudo-AML, 20 draws, fixed seed; the other 5 are the healthy
  arm and give C's reference; eligibility recomputed after each plant; whole-T rebuilt as the sum of the labels. Each
  plant is scored twice: with its own gene list (pass/fail below) and with the real fold list of the fold that left each
  donor out (realistic; reported).
  - K3a uniform state: every T label of the pseudo-AML donors gets the S5 main-fit log2FC on the P370 genes (counts
    times 2^log2FC); list P370. Pass: median |C| ≤ 0.2 and median S_w > 0. Its median S_w is the pure-state benchmark.
  - K3b composition: no expression change; the pseudo-AML naive labels (CD4 Naive, CD8 Naive) are halved and the
    effector-memory labels (CD8 EM1, EM2, CD4 EM) doubled, counts and cell numbers alike; list = the 100 genes most higher
    and the 100 most lower in pooled healthy-arm naive vs effector-memory profiles. Pass: median C ≥ 0.8 and median
    |S_w| ≤ 0.2.
  - K3c state in one label: as K3a, but only in CD4 Naive. Pass: median |C| ≤ 0.2 and median S_w > 0.
  - K3d state outside the four labels: as K3a, but only in the five other T labels. Pass: median |C| ≤ 0.2; the median
    S_w is reported (expected near 0: the "state not located" case).
- K4: Control1, Control3 and Control4 are female by expression: XIST / (XIST + RPS4Y1) > 0.9, sums over all cells of the
  sample as in 26_s4_pseudobulk.R:81 (Control4: 7,797 vs 56).
- K5: the Savage delay file has md5 5ee444a9284ae52b7afab1c8aeeced2c; 4.2's whole-T z gap is 1.3184 (± 1e-4; A1-BATCH
  12.4, same construction).
- K6: the GSE227903 sample list is exactly the 9 AML_Dx samples of the manifest; non-diagnosis columns are dropped right
  after loading, before any use; per-cell objects are loaded for AML_Dx samples only.
- K7: the 5.1 in-sample re-derivation gives R_B 0.692, R_NK 0.928, monocytes 0.568 (± 0.001).

## 7. Decision rules (proposed; the user fixes the values before the tag)

Decisions use the point estimates of C and the CI of S_w; every CI is reported.
- **Not decisive**: any known answer fails; gapT_k ≤ 0 for two or more of the four labels; C's 95% CI wider than 1;
  or S_w + C > 1.5 (the two estimates together claim far more than the whole gap).
- **"within-state"**: C ≤ 1/2 and the lower bound of S_w's 95% CI > 0.
- **"composition"**: C ≥ 2/3.
- **"mixed"**: 1/2 < C < 2/3.
- **"not composition; state not located"**: C ≤ 1/2 and S_w's 95% CI includes 0.
- The lines are the assessment's (state share ≥ 1/2 for "state", ≤ 1/3 for "composition"), restated for C = 1 − state
  share; they echo the 1/3 line of GATE_A1(4) and S9.
- Sensitivities (section 4) and readouts (section 5) never change the wording. If one of them would put C in a
  different category, the manuscript reports it next to the primary.
- Under every outcome the manuscript states the limits of section 0: no separation of disease from source or handling;
  "state" means inside a BoneMarrowMap label; label purity in AML unverified.

## 8. Code, outputs, compute (plan; scripts are written after the tag)

- One script, `scripts/ws/39_a1_g01.R`, modes: `--known` (K0-K7), `--primary` (C, S_w, bootstrap), `--sens` (4.1-4.5;
  refits 25 + 25 + 22 folds), `--share` (5.1), `--regress` (5.3), `--combine` (tables, GATES block A1-G01). A second
  script, `scripts/ws/40_a1_g01_cells.R`, for the per-cell steps: 5.2 (reads the 25 annotated objects) and, if kept, G09
  (5.4: per-cell T1 in Oetjen and the two AML cohorts). sbatch wrappers as for A1-BATCH.
- Outputs: `results/tables/ws/ws_a1_g01_{checks,labels,primary,boot,sens,share,gdmait,regress,axis,g09}.csv`.
- Smoke runs first, on synthetic counts with real metadata; smoke tables deleted after the real run reproduces.
- Compute: every mode except 40 runs in under an hour on one node; 40 needs about 64 GB for the annotated objects.
- No new package. No download.

## 9. For the user to decide before the tag

1. **Primary = C, with S_w as the required check** (recommended after review 1), or the assessment's original
   (the unweighted within-label fraction as primary). Reason: the within-label fraction depends on how the shift is
   spread over subtypes (it exceeds 1 when the shift sits in one small label and is 0 when it sits outside the four);
   C uses all nine labels and all 25 samples.
2. Decision lines of section 7 (recommended as written: C ≤ 1/2 and C ≥ 2/3).
3. G09 (5.4) inside this amendment, as the assessment ordered, or as its own short amendment right after G01
   (recommended: separate, so this one stays small and finishes in one or two sessions).
4. The gamma-delta/MAIT arm (5.2) now (recommended: yes; one job of about an hour).

## 10. Review history

- 2026-11-01: DRAFT v1 (ff53d2b, 0ef4ae0, bd29ce6).
- 2026-11-01: short review 1 (one agent; no blind quantity computed): 1 blocking, 5 should-fix, 4 minor; all accepted.
  Blocking: v1's rules could return "composition" when C said otherwise (a shift outside the four labels gave S ≈ 0
  and C ≈ 0) and had no category for that case. Should-fix: the unweighted within-label fraction equals about 1/(4w)
  for a shift confined to one label of UMI share w, so C became primary and the within-label part a weighted check
  (S_w); numerator and denominator now share one gene list, with a 5-gene floor; K3 now includes one-label and
  outside-labels plants and realistic fold lists, and the calibration step was dropped; 5.1 is read one-sided; v1 had
  not disclosed the verifier's recorded pooled shares and CD8 Naive eligibility. Minor: two values are in
  `ws_a1_g02_descriptive.csv`, not only in memory; K4 states its rule; K6 reworded for S8's loader; G09(a) given a
  proportional rule with Oetjen as the only clean healthy cohort. Added: 5.5 (does the within-label shift follow the
  naive-to-effector axis).
- 2026-11-01 23:40: DRAFT v2.
