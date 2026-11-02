# PRE-REGISTRATION A1-G01 — how much of the T-cell AML–healthy difference is T-subset composition, and how far does it reach beyond T cells?

WP-A1 gap G01 (composition versus within-subtype state), with the four items that A1-BATCH section 12.5 carried
forward: (a) "carry the delay score as a declared sensitivity" (4.2); (b) "ask whether any within-subtype residual is
also present in B and NK cells of the same donors": 5.1 tests the whole-T direction there, because B and NK contain no T
subsets, so the residual itself is not separable in them; (c) word the result "beyond T-subset composition", not
"T-cell-intrinsic" (section 7); (d) declare the Control4 sex coding and a sensitivity without the three
expression-female healthy donors (4.3, 4.4). G09 (what S8's structural replication can show) becomes its own short
amendment after G01 (user decision 2026-11-02).

**Status: v1, APPROVED by the user 2026-11-02 (all section-9 items as recommended), after reviews 3 and 4 (section 10);
tag `prereg-a1-g01-v1`. Nothing declared below was computed before the tag.** Order (user decisions 2026-10-28,
2026-11-01 and 2026-11-02): A1-G02 (PASS) → A1-BATCH (PASS with qualifications) → **this** → G04(a)+(b). Handoff v2.0
§17 S15.

## 0. Question, scope, and what each outcome lets the paper say

S5 found a T-cell imprint in 15 adult AML vs 10 healthy donors of GSE185381. A1-G02 showed it out of sample
(cross-fitted gap 1.378, i.e. 1.98 out-of-fold healthy-SD; full-pipeline permutation p 0.008). A1-BATCH 12.5 allows:
"none of the six within-study processing factors listed, each tested one at a time, explained or weakened it", and
requires: "Sample source, storage, sorting method and chemistry version are constant within arms or unrecorded, and were
not tested." It also found that the same direction appears in B and NK cells of the same donors (and partly in
monocytes), so it is not T-specific; that it is aligned with a processing-delay program; and that the part that goes
with the Control1-5 pool is tied to the three expression-female healthy donors.

S5 adjusted for one composition covariate (the CD8 share of kept T cells). Blood studies report fewer naive and more
effector-memory T cells in AML; marrow evidence is weaker (Knaus 2018: the marrow naive and TEMRA shifts were not
significant; Lasry 2023, this dataset: MAIT and NK cells diminished, cytotoxic and Treg populations unchanged). **G01
asks how much of the gap is the mix of T subtypes (composition) and how much is not.**

- **C (primary): the composition share, in both orders.** The gap that the observed subtype mix alone would produce if
  every subtype kept one reference expression profile per cell, as a fraction of the whole-T gap. All nine T labels,
  all 25 samples. C_H uses pooled healthy profiles, C_A pooled AML profiles, and C is their mean. Pure composition
  gives about 1 and a pure shift inside subtypes about 0, in both orders. The orders differ when a shift inside a label
  interacts with that label's change in share (for example a shift inside CD4 Naive, which AML depletes); C splits that
  interaction equally (the two-factor Shapley split). 1 − C is the **non-composition residual**: shifts inside subtypes,
  but also ambient RNA beyond the 32 removed genes, per-label differences in RNA content, handling, label impurity and
  composition below the label level.
- **S_w (required check): is there a shift inside the four S_w labels?** The four S_w labels ("the main labels") are the
  four largest T labels, fixed by name in 2.2. S_w is the gap measured inside each of them against the whole-T gap of
  the same samples and genes, combined with UMI-share weights. It is a sign check: its size depends on where the shift
  sits and on that label's share of each gene's expression (about 0.3-3 for a shift in one label; 0.5-1.1 with
  realistic lists; review 3), it can exceed 1 − C, and it is not a lower bound on anything. A shift confined to the
  other five labels gives S_w about 0.

What each outcome lets the paper say (rules and handover in section 7):
- **"composition"**: the difference is mainly the subtype mix, which in GSE185381 is confounded with age. The paper
  becomes the composition-first skeleton C.
- **"beyond T-subset composition"**, worded either "detected inside the main labels" or "not resolved inside the four
  main labels". Skeletons A/B stay open.
- **"mixed"**: part composition, part not; both sizes reported, with "order-dependent" or "composition not established"
  where it applies. Skeletons A/B stay open, with part-composition wording.
- **"not decisive"**: a scope-P known answer failed, or C_H or C_A could not be determined. Skeleton A, with composition
  stated as unresolved.

What G01 cannot do, under any outcome:
- **Separate disease from sample source or handling.** A handling effect on the whole sample also appears inside every
  subtype, so "beyond composition" does not mean "caused by AML". Only G04(a) (Chen2023, one protocol for both arms;
  5 vs 4, so a pass is informative and a fail is not) addresses that.
- **Separate a cell-intrinsic shift from ambient RNA.** Soup from blasts is in every AML sample and appears inside every
  label, so it lands in 1 − C and in S_w. The 32-gene filter removes only known ambient-prone genes, S5's ambient gate
  failed, and A1-BATCH D9 (not interpretable) recorded the up genes enriched and the down genes depleted in soup.
  G04(a)'s FACS-separated T cells reduce this, but no step tests it.
- **Separate composition from age.** The healthy donors are about 31 years younger than the patients (S9: 31.5 y), and
  the CD8 Naive share falls with age in healthy marrow (S9 K1, Oetjen: rho −0.77). So "composition" here is partly age,
  and because the healthy profiles come from younger donors, 1 − C can also carry age-related drift inside labels.
  Chen2023 (G04(a)) is age-matched (AML 39-59 y, donors 47-53 y).
- **Resolve below the label level.** "Inside a subtype" means inside a BoneMarrowMap fine label. A composition shift at
  a finer level (TEMRA-like cells inside "CD8 Effector Memory 2"; gamma-delta and MAIT cells, which have no label)
  counts as non-composition here; 5.2 and 5.5 measure the two cases they can.
- **Verify label purity in AML.** If AML cells labelled "CD4 Naive" are less naive than healthy ones, S_w is biased
  upwards and C downwards; if shifted cells are relabelled as memory, C is biased upwards. The protein (ADT) check is
  deferred (assessment G01 option b).
- **Cover all T cells.** C and S_w use kept cells only (BoneMarrowMap label concordant with marker typing, not high
  error; "T Proliferating" is not in T). If AML T cells fail these filters more often, the result describes a selected
  subset; kept fractions are reported (5.6).
- **Capture composition along the CD4/CD8 axis.** The fold lists were selected with the CD8 share in the model and
  |log2FC| > 0.5, so they hold few composition-driven genes, and adjusting for composition also removes any AML effect
  that acts through composition. "Beyond composition" is conditional on these lists; "composition" does not mean
  "unrelated to AML".

## 1. Already known before this draft (disclosed)

Seen, from earlier steps:
- Whole-T: gap_oof 1.378164 (A1-G02 R2), per-sample out-of-fold scores (`ws_a1_g02_oof_scores.csv`) and per-fold gene
  lists (`ws_a1_g02_fold_genes.csv`, 25 leave-one-out folds, 32 ambient-prone genes removed; median 329 genes; the
  smallest lists are 11 up / 18 down with AML4340 left out, 27 / 48 with Control4 and 21 / 65 with Control3; every other
  fold has at least 51 / 76).
- Composition-related, all in-sample or descriptive: S5 adjusted for the CD8 share. Spearman of the S5 T score with the
  CD8 share is 0.168 (15 AML) / 0.224 (10 healthy), and of the S8 score with the CD8 share in GSE227903 diagnosis
  samples 0.857 (n 7) / 0.829 (n 6, the S10 subset) (`ws_a1_g02_descriptive.csv`). In Oetjen, S9's age share f_slope
  (a fraction of the yardstick gap) goes from 0.102 to −0.047 with the naive share in the model, and the CD8 Naive share
  falls with age (rho −0.767). Refits inside CD4 Naive (12 AML + 8 healthy, 85 hits; sign concordance 1.00 over the
  imprint genes tested there) and CD8 EM2 (12 + 9, 7 hits, 0.996), in-sample, not cross-fitted. Canonical subset
  markers (CCR7, TCF7, LEF1, GZMB, GZMH, NKG7, PRF1, GZMK) are not hits. The up arm holds CXCR4 (log2FC 1.30), S1PR5
  (1.18) and ZEB2 (0.80); the down arm holds naive/resting-associated genes (LINC02446, TMIGD2, RGS10, LDHB, NOSIP,
  SATB1-AS1, TXK, CAMK4), PRKCQ-AS1 (−1.17), NUCB2 (−0.82) and TRDC (−1.86). The ambient filter removed SELL, IL7R,
  FHIT, RCAN3, RUNX3, CD5 and MAF from the 402.
- The per-sample CD8 share (S5's covariate) is on disk in `ws_s5_scores.csv`, and the A1-G02 permutation table holds the
  mean CD8 share of each permuted "AML" group; no arm mean under the real labels appears in any outcome section. Review
  3 saw two per-sample CD8 shares while checking column layout.
- Composition table: read for eligibility counts (section 2.2). The G01 verifier (2026-10-28) also recorded pooled
  shares of the 14,609 kept T cells (both arms together: CD4 Naive about 17%, CD8 EM2 about 39%) and that CD8 Naive
  reaches 20 kept cells in 1/15 AML vs 5/10 healthy samples, and wrote that these counts suggest naive share and status
  barely overlap. No per-arm share value is recorded.
- A1-BATCH 12.3-12.4 (post-hoc parts in-sample). P370 retention in B 0.69 and NK 0.93, sign concordance 0.97 / 0.95,
  above all 60 random relabellings; monocytes R 0.57, sign concordance 0.87 (11 AML + 10 healthy); genome-wide T-B logFC
  correlation 0.60 lies inside donor-grouping nulls, T-NK 0.70 at their top. Delay: the imprint overlaps the Savage 2021
  delay genes in the same direction; the cross-fitted gap without delay genes is 1.3184; the delay score still tracks
  the out-of-fold score after that removal (Spearman 0.67 within healthy, 0.55 within AML); MAIN + delay score gives
  R 0.823. D5 (Control4 coded female): R 1.019, CF/1.378 1.03, leave-pair-out AUC 0.913; males only by expression
  R 0.878. D3: CF without Control1-5 is 0.123, with 14 of its 20 fold imprints empty (CF tracks the hit count, Spearman
  0.82). D9 (not interpreted): P370 up genes log2 soup/normal-T +0.31 (random −0.13), down genes −0.72 (random −0.11).
  The 9 of 40 matched drops that remove Control1, Control3 and Control4 have R 0.75-0.92, the other 31 0.89-1.18.
- Review 3 (2026-11-02, section 10) read, for its checks: the XIST shares of Control1, 3 and 4 (0.994, 0.994, 0.993;
  Control4 7,797 vs 56 UMI) and their stored sex calls; pseudobulk names and dimensions; the D10 and monocyte tables
  above. Review 4 read documents, code and headers only (one incidental `head -5` of
  `ws_s4_pseudobulk_meta_GSE227903.csv` showed one diagnosis sample's per-compartment cell numbers, not a G01
  quantity). All simulations of both reviews were synthetic.

Not computed (blind at this draft): any within-label score, gap or gene-level shift; any composition counterfactual;
any held-out B, NK or monocyte score; any per-arm subtype share, kept fraction or naive-share coefficient; S8's T1/T2
among healthy donors or against composition-matched random sets (G09); gamma-delta/MAIT flags; any GSE227903
post-treatment or relapse expression (S10 stays blind); Chen2023 and the GSE185381 children (G04).

## 2. Inputs and frozen definitions (fixed before any expression is read)

### 2.1 Samples, folds, engine
- Samples: the 25 of the S5 T design (`ws_s5_scores.csv`, compartment T; 15 AML, 10 healthy), file order. Every
  pseudobulk object is subset to these 25 columns right after loading (the GSE185381 file also holds the paediatric
  samples, which G04 keeps blind).
- Folds: the 25 leave-one-out folds of A1-G02 R2. Fold i's lists U_i (up) and D_i (down) are read from
  `ws_a1_g02_fold_genes.csv` and **not reselected** for the primary analysis (selected on the other 24 whole-T
  pseudobulks; design MAIN = status + hb + gran + sex + subtype_share, metadata sex; edgeR QL; FDR < 0.10 and
  |log2FC| > 0.5; 32 ambient-prone genes removed). Refits (4.1, 4.2a, 4.3, 4.4) use the same code; K0b checks it.
- Log scale: `logz()` of 37_a1_batch.R (blacklist; CPM > 1 in at least half of the columns; TMM; log-CPM, prior 1),
  applied to the matrix being scored.

### 2.2 Subtype labels, eligibility, weights
- Labels: the nine BoneMarrowMap fine labels mapped to T in `s4_cell_groups.tsv`. Counts: `PB$fine[[label]]` of
  `ws_pb_GSE185381.rds` (raw counts summed over kept cells). A sample without kept cells of a label has no column there,
  so absent columns are zero-filled (n_ks = 0). Kept-cell numbers from `PB$cells`.
- Eligible: at least 20 kept cells of the label in the sample (WS_S4_MIN_CELLS).
- **The four S_w labels ("the main labels"), fixed by name**: CD4 Naive (12 AML + 8 healthy eligible), CD8 Effector
  Memory 2 (12 + 9), CD8 Effector Memory 1 (11 + 7), CD4 Central Memory (10 + 7) (the verifier's counts,
  verify_G01.json items 7 and 10; rechecked as K1). The other five labels enter C but not S_w.
- Weights: π_k = label k's share of whole-T UMIs over the genes of the whole-T logz universe, summed over the 25 samples
  (both arms together); fixed in the bootstrap and in K3.

### 2.3 Score units
- **u-score (the unit of every ratio)**: for sample i scored with lists U and D in a log-CPM matrix X,
  u_i = mean over g in U of (X[g, i] − m_g) − mean over g in D of (X[g, i] − m_g), where m_g is the mean of X[g, ] over
  the other samples of the same matrix. A per-sample constant (library size, TMM factor) cancels between the halves.
- **z-score**: the same with (X[g, i] − m_g) / s_g, s_g the SD over the other samples. Used only to reproduce A1-G02
  (K0, K5). No z-version of a ratio is reported: a label pseudobulk's SD carries more counting noise than the whole-T
  one, and the synthetic matrix of 3.1 carries none, so z-ratios would mix units.
- Remaining bias of u: the prior count weighs more in small libraries, which shrinks shifts measured in label matrices
  (K3a measures this); and when the U and D genes sit at different expression levels, a difference in library size
  between the arms moves a gap with no shift at all (the thinned K3b measures this offset).

### 2.4 One gene list per comparison
Every ratio scores its numerator and its denominator on the same genes: fold i's genes that lie in the universes of all
matrices of that comparison (for C: both synthetic matrices and X_T). A sample (or sample-label pair) whose matched list
has fewer than 5 genes in either half is dropped from that comparison, and the dropped IDs are reported. Reason: an
empty half scores as 0 (`score_col`, 37_a1_batch.R), and the smallest stored lists have 11 up / 18 down genes.

## 3. Primary estimate and required check

### 3.1 C — the composition share, in both orders (primary)
For each fold i and each order o (H = healthy reference, A = AML reference):
1. Reference profile per label k (all nine T labels): r_k^o = the summed label-k counts of the arm-o samples other than
   i, divided by their summed kept label-k cells (a per-cell profile pooled over donors). If those samples hold fewer
   than 50 kept label-k cells, r_k^o is the label's profile pooled over all samples other than i (both arms), so the
   label keeps its identity; if those hold none, r_k^o = 0 (the pair is reported). The fallback puts the other arm's
   cells into the reference: it acts like an added within-label shift in that label and can move C_H or C_A in either
   direction (review 4: by up to about 0.05 at 15 vs 10). C_H and C_A are therefore also reported with the arm-only
   profile (any number of cells > 0; both arms only if the arm has none), together with the fallback labels per fold
   and the share of bootstrap draws in which each label fell back. In K3 the both-arm pool contains planted donors
   (section 6).
2. Synthetic whole-T counts for every sample s: y_s^o = sum over k of n_ks · r_k^o, where n_ks is s's kept cells of
   label k.
3. Y_o = logz of the 25 synthetic columns; X_T = logz of the 25 real whole-T columns. With fold i's matched list
   (section 2.4): c_i^o = u-score of i in Y_o; t_i = u-score of i in X_T (both centred on the other 24).

gap_C^o = mean of c_i^o over AML − mean over healthy; gapT = the same with t_i.
**C_H = gap_C^H / gapT; C_A = gap_C^A / gapT; C = (C_H + C_A) / 2.** C_H and C_A are undefined when gapT ≤ 0. C uses
only each sample's label cell numbers and the pooled reference profiles, never the sample's own expression.

- Per-label RNA content: C weights labels by cell numbers times the reference UMI per cell, so a label whose cells carry
  more or fewer UMIs in AML counts as non-composition (4.6 weights by UMI instead).
- Count noise: the synthetic columns have none and X_T has some. Under pure composition this biases C slightly
  downwards (about 3-6% in review 3's simulations at realistic list sizes); 4.7 adds matched count noise.
- C outside [0, 1] is reported as such: above 1, composition would produce more than the whole gap (an opposing shift
  inside subtypes); below 0, composition works against the gap.
- CI: 2,000 bootstrap draws (seed offset 32) resampling AML and healthy samples separately. In each draw, steps 1-3 are
  repeated on the drawn samples (references re-pooled from them; references and centring exclude every copy of i), with
  the fold lists held fixed; 95% percentile interval. A draw in which C_H or C_A cannot be computed is undefined;
  undefined draws are counted and reported, and enter the interval as −∞ for its lower bound and +∞ for its upper bound,
  so they can only widen it. Holding the lists fixed understates the spread (1.2- to 1.6-fold in review 3's
  simulations; re-pooling did not change the width in review 4's), so the interval is descriptive and section 7 uses it
  only with a margin.

### 3.2 S_w — is there a shift inside the four S_w labels? (required check)
For each S_w label k, with E_k its eligible samples:
1. X_k = logz(label-k counts of E_k); X_T,k = logz(whole-T counts of the same E_k).
2. For every i in E_k, with fold i's list matched to X_k and X_T,k: s_ik = u-score of i in X_k; t_ik = u-score of i in
   X_T,k (both centred on E_k without i).
3. gap_k = mean of s_ik over AML − mean over healthy; gapT_k = the same with t_ik.

**S_w = (sum of π_k) × (sum over k of π_k · gap_k) / (sum over k of π_k · gapT_k)**, over the four labels. To first
order this is the π-weighted sum of the per-label ratios S_k = gap_k / gapT_k, but one label with a small gapT_k cannot
dominate it. Same bootstrap as C (a drawn sample enters every label where it is eligible). A label with no eligible
sample in an arm (in the estimate or in a draw) leaves all three sums. S_w is undefined if the denominator is ≤ 0;
undefined draws are handled as in 3.1. S_w is a sign check: its size depends on where the shift sits and in which
genes, and it exceeds 1 − C for example when the shift sits inside a label that AML depletes, or in genes that label
expresses below its UMI share. Also reported: each S_k with its gapT_k (flagged when gapT_k < gapT / 4); S_w and each
S_k for the up and down halves separately; gene coverage per label; the dropped sample IDs.

## 4. Sensitivities (each reports C_H, C_A, C and S_w with 500 bootstrap draws, seed offset 33; none changes the wording)

- **4.1 The 402-analogue**: the 25 folds refitted without removing the 32 ambient-prone genes (they include SELL, IL7R
  and other subset markers, so composition should weigh more here).
- **4.2 Processing delay (A1-BATCH 12.5(a))**: (a) the 25 folds refitted with the delay score added to MAIN (the D8
  score of 37_a1_batch.R: Savage 2021 up minus down genes, z-scored; recomputed in each fold over its 24 training
  samples, excluding that fold's own 402-analogue hits (4.1) instead of F402, so the held-out sample cannot shape the
  covariate; K8 checks the 25-sample construction); (b) the stored fold lists minus all 792 Savage 18-h delay genes
  (295 up / 497 down; `savage2021_T_delay18h_consistent.csv`, md5 5ee444a9284ae52b7afab1c8aeeced2c), from both halves
  regardless of direction.
- **4.3 Control4 sex coding**: the 25 folds refitted with Control4 coded female (A1-G02 section 7's G06 sensitivity; the
  whole-T result is already known, section 1).
- **4.4 Without the three expression-female healthy donors** (Control1, Control3, Control4): 22 folds refitted on 15 AML
  + 7 healthy (each fold trains on 21). Reported as "not estimable" if more than a third of the 22 folds fail the
  5-gene rule, or if fewer than 5 of the 7 healthy or fewer than 10 of the 15 AML samples remain in C (A1-BATCH D3
  found most 20-sample refits empty). An S_w label with fewer than 3 eligible samples in an arm leaves S_w here.
  Eligible counts per label are reported.
- **4.5 Leave one S_w label out**: S_w recomputed four times, each time without one label.
- **4.6 UMI-weighted C**: y_s^o = sum over k of U_ks · r̃_k^o, where U_ks is the sample's own label-k UMI total and
  r̃_k^o the reference profile scaled to sum 1, so per-label RNA content counts as composition. UMI per kept cell per
  label is reported by arm.
- **4.7 Count-noise-matched C**: each y_s^o is scaled to the sample's real whole-T UMI total and drawn as Poisson counts;
  u-scores are averaged over 20 draws (seed offset 34 plus the bootstrap draw index; index 0 is the estimate).

A sensitivity is reported as moving C out of the primary band (≤ 1/2, between, ≥ 2/3) only if its CI for C lies
entirely outside that band; otherwise as "consistent within uncertainty". The manuscript reports every sensitivity next
to the primary.

## 5. Further declared readouts

### 5.1 Beyond T cells, cross-fitted (A1-BATCH 12.5(b))
A1-BATCH 12.4 measured the sharing in-sample: the 370 genes were chosen on the same donors whose B and NK cells were then
read, so donor-level noise picked up by the selection can also appear in B and NK. Here each held-out donor's B, NK and
monocyte pseudobulk (`PB$raw`, S5's compartment designs restricted to the 25; only donors among the 25 have a fold) is
scored with the T genes of the fold that left that donor out, matched (section 2.4) to that compartment's universe and
to the whole-T universe of the same donors. Reported per compartment: n, gene coverage, the held-out u-gap, and its
ratio to the whole-T held-out u-gap on the same donors and genes, with a bootstrap CI (2,000 draws, seed offset 35).
- This tests the whole-T direction in each compartment. B, NK and monocyte pseudobulks contain no T subsets, so a shift
  there cannot come from T-subset composition; it can come from those compartments' own subset composition, which 5.1
  does not control.
- Reading is one-sided: a ratio whose CI lower bound is above 0 shows the direction in that compartment out of sample.
  A ratio near 0 is **not** evidence that the shift is confined to T cells, because smaller libraries and the loss of
  T-specific genes both shrink it.
- Monocytes in AML samples are often clonal, so a monocyte shift may be leukaemic biology rather than a shared
  environment or handling effect; stated with the result.
- The in-sample numbers of A1-BATCH 12.4 are re-derived in the same script and written to a durable table (K7).

### 5.2 Gamma-delta and MAIT cells (verifier P3)
BoneMarrowMap has no gamma-delta or MAIT label, and Lasry 2023 reported MAIT loss in this dataset. From the raw counts of
the kept T cells (04_annotated objects): gamma-delta flag = TRDC ≥ 1 and TRGC1 + TRGC2 ≥ 1 UMI; MAIT flag = TRAV1-2 ≥ 1
and SLC4A10 ≥ 1 UMI; a cell with both flags counts as gamma-delta. Detection is incomplete at this depth, so the flags
are lower bounds. Reported: flagged shares per sample by status (sanity range: healthy median gamma-delta share between
0.5% and 15% of T cells); C_H, C_A and C with the flagged cells moved out of their labels into two extra labels; the
whole-T u-gap after removing flagged cells (fold genes fixed); the gap with TRDC removed from the fold lists. Reading,
fixed now: under every outcome the manuscript states, next to C, "with gamma-delta and MAIT cells as separate labels the
composition share is x (the flags are lower bounds)"; the category does not change.

### 5.3 Status plus composition regression (verifier P5; secondary)
The A1-G02 out-of-fold z-scores regressed on status + naive share (CD4 + CD8 Naive / kept T) + CD8 EM share (EM1 + EM2 /
kept T), **only if an overlap check passes**: at least 3 AML and 3 healthy samples inside the common range of naive
share (from the larger of the two group minima to the smaller of the two group maxima). Otherwise reported as "not
estimable: no overlap", and nothing else is fitted.

### 5.4 G09 — moved out
By the user's decision (2026-11-02, section 9 item 3), G09 becomes its own short amendment, written after G01's outcome.
The v3 design notes carry over there: T1/T2 among Oetjen donors with a proportional rule (GSE185381's healthy donors
reported, not used); random sets matched also on the 5.5 loading; GSE227903 AML_Dx only, subset right after loading
(former K6). Nothing of G09 runs under this tag.

### 5.5 Does the within-label shift run along the naive-to-effector axis? (descriptive)
Loading: for each gene, log2((CPM_EM + 1) / (CPM_naive + 1)) between the pooled CD8 EM1 + EM2 + CD4 Effector Memory and
the pooled CD4 + CD8 Naive kept-cell pseudobulks of all Oetjen libraries (`ws_pb_GSE120221.rds`; the 20 donors of S9,
none excluded), inside its logz universe; e_g = the log-CPM of the pooled naive + EM Oetjen pseudobulk (status-free).
Within CD4 Naive and CD8 EM2 (the two S_w labels with most eligible samples), for every gene of the label's logz
universe that the loading covers: d_g = the AML − healthy difference in mean log-CPM. Statistic: the Spearman
correlation of d_g with the loading inside each decile of e_g, averaged over the ten deciles. Null: 1,000 status
relabellings within the label's eligible samples, keeping its arm sizes (seed offset 36), the statistic recomputed each
time. Reading, fixed now: if the statistic is > 0 and its relabelling p < 0.05, the manuscript adds "the
within-label shift partly follows the naive-to-effector axis, so composition below the label level may contribute";
otherwise nothing. The same correlation over the union of the 25 fold lists is reported as a number only (selection
makes it positive whenever there is whole-T composition; review 4). In-sample and descriptive; the category does not
change.

### 5.6 Reported for every estimate
Both gene sets (370-analogue primary, 402-analogue in 4.1); n per label and compartment; gene coverage; dropped sample
IDs (5-gene rule); the share of kept T cells in the four S_w labels, by status; the kept fraction of T-lineage cells per
label and the "T Proliferating" share, by status (section 0, kept cells only).

## 6. Known answers (run first)

Each known answer has a scope. P covers the estimator code that every mode uses, so it blocks every mode; the others
block, in addition, only the sensitivity or readout named. The checks table has three states: PASS, FAIL
(deterministic) and FAIL (statistical, final).
- Deterministic checks must pass exactly. If one fails, the cause is a code or input error: it is fixed, the fix is
  recorded in a dated note in this file before any readout of that scope is printed, and all known answers rerun.
- The statistical bars of K3 run once, with the declared seeds, after the deterministic checks pass. If a bar fails,
  the result is "not decisive" with no change to methods, thresholds or seeds; the modes of its scope still run, their
  outputs are labelled "not decisive, descriptive", and the GATES block names the failed bar. A bar failure traced to a
  demonstrable code error (shown on synthetic data) is handled as a deterministic failure.
- A mode runs only when every deterministic check of its scope passes (non-smoke checks table).

- K0 (P): the whole-T z gap from the stored fold lists reproduces 1.378164 (|difference| < 1e-6).
- K0b (4.1, 4.2a, 4.3, 4.4): the refit code run with the primary settings reproduces the 25 stored fold lists exactly.
- K1 (P): eligible counts of the four S_w labels equal 12 + 8, 12 + 9, 11 + 7, 10 + 7. The pooled eligible totals of all
  nine labels (not split by status) are printed for the record; nothing is asserted about the other five.
- K2 (P; construction check): for each of the 25 samples, the nine zero-filled T-label count vectors sum exactly to the
  whole-T pseudobulk.
- K2b (P; leave-one-out check): for every fold i and both orders, the reference profiles computed with sample i's label
  counts set to zero are bit-identical to those of fold i.
- K3 (P), planted on the 10 healthy donors only. 20 splits into 5 pseudo-AML and 5 pseudo-healthy, drawn without
  repetition from the 252 (seed offset 31). In each split, sections 2-3 run within the 10 donors (C_H's reference = the
  pseudo-healthy donors other than i, C_A's = the pseudo-AML donors other than i, with 3.1's fallback); eligibility is
  recomputed after each plant; whole-T is rebuilt as the sum of the labels; π_k stays fixed from the 25; S_k needs at
  least 2 eligible donors per arm (otherwise the label is dropped, leaves all three sums of S_w, and is reported). Each
  plant is scored with its own list (the bars) and with the real fold list of the fold that left each donor out
  (reported). Paired value of a ratio: computed from (gap with plant − gap without plant) over (gapT with plant − gapT
  without plant), in the same split; for S_k and S_w, every gap and gapT entering them is paired this way.
  - K3a uniform state: every T label of the pseudo-AML donors gets the S5 main-fit log2FC on the P370 genes (counts
    × 2^log2FC); list P370.
  - K3b composition: no expression change. Binomial thinning (rbinom, p 0.5; cell numbers halved, floored; seed offset
    38, one pass per split in split order) of CD4 Naive and CD8 Naive in the pseudo-AML donors and of CD8 EM1, CD8 EM2
    and CD4 Effector Memory in the pseudo-healthy donors. List: from the 5 unplanted pseudo-healthy donors, pooled CD4 +
    CD8 Naive against pooled CD8 EM1 + EM2 + CD4 Effector Memory, ranked by log2((CPM_EM + 1) / (CPM_naive + 1))
    inside the whole-T logz universe, ties broken by gene name; the top 100 are "up", the bottom 100 "down".
  - K3c state in one label: as K3a, only in CD4 Naive.
  - K3d state outside the four S_w labels: as K3a, only in the five other T labels.
  - K3e composition with a state inside a depleted label: K3b's thinning plus K3c's CD4 Naive state; list = the K3b list
    together with P370. Reported only (C_H > C_A is expected).
  - Deterministic (reference tests). For K3a, K3c and K3d, in every split, for every held-out i and every label k whose
    reference did not fall back: r_k^H is bit-identical with and without the plant; r_k^A differs with the plant when
    k is planted, and is bit-identical with and without the same plant applied to the pseudo-healthy donors instead
    (one extra plant run; references only). Labels that fell back are counted per split, order and i, and paired C_H
    is reported for those splits (not a bar).
  - Statistical bars (medians over the 20 splits, own lists, paired). K3a: |C_A| ≤ 0.10, every kept label's S_k in
    [0.3, 1.3], and S_w ≥ 0.3 × (sum of π_k over the kept labels). K3b: C_H ≥ 0.8 and C_A ≥ 0.8. K3c: S_w ≥ 0.2. The
    reference tests catch reference construction errors (a C_A built from the healthy arm passes every bar; review 4);
    the bars catch gross estimator errors and are loose because the prior count shrinks shifts in small healthy label
    libraries.
  - Reported only (operating characteristics; 200 bootstrap draws per split for K3b, K3d and K3e, seed offset 37):
    unpaired C_H, C_A and C; S_w (K3b's S_w is the no-state offset from library-size differences between the arms);
    each S_k; the section 7 category of each split; for K3d, the share of splits whose S_w lower bound exceeds 0.10.
- K4 (4.4): Control1, Control3 and Control4 have sex_inferred "female" (26_s4_pseudobulk.R:91: XIST > 3 × max(RPS4Y1, 1),
  sums over all cells of the sample). Their XIST shares were seen in review 3 (section 1), so this is a construction
  check.
- K5 (4.2b): the Savage file has md5 5ee444a9284ae52b7afab1c8aeeced2c; the stored fold lists minus all 792 Savage genes
  (both halves, regardless of direction), with z-scores as in K0 and no list matching, give 1.3184 (± 1e-4; source
  `/FAST/gr10634/gaozy/reports/2026-11_a1_batch_posthoc/code_audit/a1_posthoc_d8_shared_genes.R` and `.out`).
- K6: moved with G09 (5.4).
- K7 (5.1): R = R_of(fit_de(PB$raw[[cmp]], S5 compartment design restricted to the 25, MAIN), min_cov = 0) gives
  0.692 (B, 174 genes), 0.928 (NK, 179 genes) and 0.568 (monocytes, 329 genes) (± 0.001; `ws_a1_batch_d10.csv`,
  `p5_mono.csv` of the post-hoc folder).
- K8 (4.2a): MAIN + the D8 delay score over the 25, with 37's code path (F402 excluded), reproduces A1-BATCH's R 0.823
  (± 0.001; `ws_a1_batch_d8.csv`, row Savage, 0.823025).

## 7. Decision rules (fixed by the user 2026-11-02)

Evaluated top to bottom; the first rule that matches decides.
1. **"Not decisive"**: a known answer of scope P failed under section 6's rule (deterministic, or a statistical bar);
   or C_H or C_A is undefined.
2. **"Composition"**: C_H ≥ 2/3 and C_A ≥ 2/3, and the lower bound of C's 95% CI > 1/2.
3. **"Beyond T-subset composition, detected inside the main labels"**: C_H ≤ 1/2 and C_A ≤ 1/2, and the lower bound of
   S_w's 95% CI > 0.10. An undefined S_w (estimate or bound) does not meet this rule.
4. **"Beyond T-subset composition, not resolved inside the four main labels"**: C_H ≤ 1/2 and C_A ≤ 1/2.
5. **"Mixed"**: every other case. Worded "order-dependent" when C_H and C_A fall in different bands (≤ 1/2, between,
   ≥ 2/3) and differ by at least 0.10, and "composition not established" when both are ≥ 2/3 but C's lower bound is
   ≤ 1/2.

Why these lines:
- 1/2 and 2/3 are the assessment's lines (state share ≥ 1/2 for "state", ≤ 1/3 for "composition"), **moved** from the
  within-label fraction to C (1 − C, the non-composition residual, takes the place of the assessment's state share).
  They echo the 1/3 line of GATE_A1(4) and S9.
- Both orders must agree, because C_H and C_A can differ by more than a band when a shift sits inside a label whose
  share changes (review 3: 13 of 24 synthetic configurations changed category between the orders). The 0.10 floor on
  "order-dependent" is above every non-interacting difference simulated in review 4 (at most 0.043) and below every
  interacting one (at least 0.12).
- Only "composition" moves the paper to skeleton C, so only it carries a CI condition. The lower bound must clear 1/2, a
  margin of 1/6 below the line, because the CI is conditional on gene selection.
- The 0.10 margin on S_w: with lists held fixed, the bootstrap understated S_w's spread 1.6- to 2-fold in review 3's
  simulations, and a shift outside the four labels then passed "> 0" in 5-14% of runs instead of 2.5%. A margin of 0.10
  brought this to at most 2.5% (review 4: 1/30), at a cost of about 0.1 in detection when the shift sits inside one
  depleted label.

Wording under every outcome: C_H, C_A, C and S_w with their CIs, and the limits of section 0. Under "composition",
also S_w's sign and CI, S9's f_age (0.102; 95% CI −0.085 to 0.463), and "in GSE185381 the subtype mix is confounded with
age". Sensitivities and readouts never change the category (section 4 says how they are reported).

What each outcome hands on (fixed now, so that the next steps' estimands are not chosen after this result). "The
counterfactual" handed to G04 is computed in both orders and their mean, as in 3.1; S10 uses a within-study reference
(GSE227903 has no healthy donors) that S10's pre-registration fixes before any of its expression is read.
- **"Composition"**: G04(a)'s primary is the counterfactual and the naive share, AML vs donors (Chen2023 is age-matched,
  so this is also the age check). G01 does not propose G04(c). S10 reports the counterfactual next to the whole-T
  score, because T reconstitution after chemotherapy is itself a composition change.
- **"Beyond T-subset composition"** (either wording) **and "mixed"**: G04(a)'s primary is the whole-T direction; the
  residual (whole-T minus the counterfactual) is reported. Its expected size is about (1 − C) × 1.98 healthy-SD, and
  5 vs 4 has power of about 0.27 at 1 SD and 0.74 at 2 SD, so G04(a) can support the direction but rarely the
  residual's size. G01 therefore proposes G04(c) (E-MTAB-16321; the download needs approval). S10 inherits the whole-T
  score, plus the counterfactual under "mixed".
- **"Not decisive"**: G04(a) tests the whole-T direction only; G01 does not propose G04(c); S10 reports the
  counterfactual next to the whole-T score.
- G04(b) (the GSE185381 children) shares G04(a)'s primary under every outcome and is read one-sided: an AML-side result
  counts, a null cannot be read (age and naive composition both push children towards healthy).

## 8. Code, outputs, compute (plan; scripts are written after the tag)

- `scripts/ws/39_a1_g01.R`, modes: `--known` (K0, K0b, K1, K2, K2b, K3, K4, K5, K7, K8), `--primary` (C, S_w, bootstrap, 5.6), `--sens`
  (4.1-4.7; refits 25 + 25 + 25 + 22 folds), `--share` (5.1), `--regress` (5.3), `--axis` (5.5), `--gdmait` (C with
  the 5.2 labels; after 40), `--combine` (tables, rules, GATES block A1-G01).
- `scripts/ws/40_a1_g01_cells.R`, the per-cell step of 5.2: reads the 25 annotated objects; writes the flagged shares
  and the label counts with the two extra labels.
- Both copy the functions they need (logz, score_col, fit_de and R_of from 37), each with a comment naming its source
  line. They never run or source 26, 27, 32, 33, 35, 36 or 37, which write frozen tables and GATES blocks.
- Outputs: `results/tables/ws/ws_a1_g01_{checks,labels,primary,boot,sens,share,gdmait,regress,axis}.csv`; each records
  the md5 of the script copy that wrote it. Each sbatch job runs a dated copy of its script (an edit on disk lands in a
  running Rscript). The checks table holds the three states of section 6.
- Smoke runs first, on synthetic counts **and synthetic label cell numbers** (real sample IDs and status only), so that
  no per-arm composition is printed before the known answers pass; every fwrite is audited; smoke tables carry `_smoke`
  and are deleted after the real run reproduces.
- Compute (16 cores): `--known` about 30 min (K3: 6 states × 20 splits, plus 200 bootstrap draws for K3b, K3d and
  K3e); `--primary` about 10 min; `--sens` about 1 h; 40 needs about 64 GB for the annotated objects.
- Seeds are offsets added to SEED (491638, config_paths.R): K3 splits 31; primary bootstrap 32; sensitivities 33; 4.7
  Poisson draws 34 plus the draw index; 5.1 bootstrap 35; 5.5 relabellings 36; K3 bootstraps 37; K3b thinning 38.
- No new package. No download.

## 9. Decisions (the user, 2026-11-02: all four as recommended)

1. Primary = the composition share in both orders (C_H, C_A and their mean C), with S_w as the sign check.
2. The decision lines and the handover of section 7 (both orders ≤ 1/2 or ≥ 2/3; "composition" also needs C's lower
   bound > 1/2; "detected inside the main labels" needs S_w's lower bound > 0.10). Review 4's corrections after the
   decision change no line; they add G04(b) to the handover (same primary as G04(a), read one-sided), name the order of
   the counterfactual handed on, and add a 0.10 floor to the "order-dependent" wording.
3. G09 becomes its own short amendment after G01 (5.4; K6 moves with it).
4. The gamma-delta/MAIT arm (5.2) runs now.

## 10. Review history

- 2026-11-01: DRAFT v1 (ff53d2b, 0ef4ae0, bd29ce6).
- 2026-11-01: short review 1 (one agent; no blind quantity computed): 1 blocking, 5 should-fix, 4 minor; all accepted.
  Blocking: v1's rules could return "composition" when C said otherwise (a shift outside the four labels gave S ≈ 0
  and C ≈ 0) and had no category for that case. Should-fix: the unweighted within-label fraction exceeds 1 for a shift
  confined to one small label, so C became primary and the within-label part a weighted check (S_w); numerator and
  denominator now share one gene list, with a 5-gene floor; K3 now includes one-label and outside-labels plants and
  realistic fold lists, and the calibration step was dropped; 5.1 is read one-sided; v1 had not disclosed the
  verifier's recorded pooled shares and CD8 Naive eligibility. Minor: two values are in `ws_a1_g02_descriptive.csv`,
  not only in memory; K4 states its rule; K6 reworded for S8's loader; G09(a) given a proportional rule with Oetjen as
  the only clean healthy cohort. Added: 5.5 (does the within-label shift follow the naive-to-effector axis).
- 2026-11-01 23:40: DRAFT v2 (a239a0e, a0a0626). Review 2 (one agent) ended with the session, without output.
- 2026-11-02: review 3, a workflow: four lenses (estimator simulation; operating characteristics of the decision rules;
  known answers against the code; consistency, disclosure and blindness), each followed by an adversarial verifier,
  then a completeness critic. All simulations were synthetic; no declared quantity was computed. Record:
  `/FAST/gr10634/gaozy/reports/2026-11_a1_g01_prereg_review/` (outside the repository). After verification: 2 blocking,
  15 should-fix, 28 minor; the critic added 3 should-fix and 4 minor.
  - Blocking, both accepted. (1) v2's C used healthy profiles only; the AML-profile order put C in a different category
    in 13 of 24 configurations, and "composition" could be returned with S_w clearly above 0 → both orders, and
    decisions need them to agree (3.1, 7). (2) The rules were not exhaustive (C ≤ 1/2 with S_w's CI wholly below 0
    matched nothing; undefined values and precedence were unstated) → first match wins, and "mixed" is the remainder
    (7).
  - Should-fix, accepted. S_w is not a lower bound on the state share → worded as a sign check, and the "S_w + C > 1.5"
    guard is dropped (0, 3.2, 7). S_w's percentile CI is too narrow with per-fold lists → the 0.10 margin (7). "C's CI
    wider than 1" fired mostly under composition → dropped; "composition" needs C's lower bound > 1/2 instead (7).
    K3's |C| bars measured 5-vs-5 noise, and a correct pipeline failed them in up to a third of simulated cohorts → a
    leak test, paired ratios, and a failure rule that separates code errors from statistical failures (6). K3b scaled
    counts, which hides library-size offsets → binomial thinning (6). C counts per-label RNA content as non-composition
    → stated, plus 4.6 (3.1). Ambient RNA, age and kept-cell selection were missing from the limits, and 1 − C is the
    non-composition residual, not a cell-intrinsic state (0). The categories now read "beyond T-subset composition",
    and the lines are said to be moved, not restated (0, 7). Section 1 lacked A1-BATCH D3, D5, D9 and 12.4.3 (1).
    12.5(a) asked for the delay score, not gene removal → 4.2(a) (4); 12.5(b) is answered by the whole-T direction in
    B and NK (intro, 5.1). The marrow literature was overstated (0). Nothing handed the result on to G04(a) and S10
    (7). Gamma-delta/MAIT and sub-label composition had no reading rule (5.2, 5.5).
  - Minor, applied: S_w as a ratio of sums; the bootstrap re-pools references; C outside [0, 1] worded; K3b's list,
    π_k and the per-arm minimum specified; K1, K2, K4, K5, K6 and K7 constructions stated; the label fallback keeps the
    label's identity; Control4's list size and dropped IDs; 5.5's loading defined in 5.5 and given a mode; the CD8 share
    on disk disclosed; section 1 wording; A1-BATCH 12.5's sentences quoted; smoke runs with synthetic label numbers;
    seeds; execution safeguards; no z-version of a ratio; count-noise-matched C (4.7).
  - Changed from the reviewers' proposals: K3's bars are looser than proposed (S_k in [0.3, 1.3] rather than
    [0.7, 1.1]; K3c S_w ≥ 0.2 rather than 0.3), because prior-count shrinkage in small healthy label libraries is real
    and the bars must catch construction errors without measuring it; "composition in at least 70% of K3b splits" is
    reported, not a bar (it would test 5-vs-5 noise); a size threshold on S_w for "detected inside the main labels" was
    replaced by the 0.10 margin, because S_w's size depends on where the shift sits; 4.4's per-arm minimum was first
    replaced by a fold-count rule (v3) and is now both (review 4).
  - Not adopted: direct standardisation (critic; a third estimator whose fallback has a known downward bias, and no
    decision needs it); an age-predicted composition share (critic, optional; Chen2023 in G04(a) is the age check); an
    equal-depth thinning sensitivity (the library-size offset is small and works against "detected").
- 2026-11-02: DRAFT v3 (52b0561). The user decided section 9 (all four as recommended).
- 2026-11-02: review 4, a short check of v3 (workflow: a fix audit and a synthetic simulation of the v3 rules and K3,
  each followed by an adversarial verifier; record in `review4/` of the same folder). After verification: 1 blocking
  (found by both checkers), 7 should-fix and 12 minor, partly overlapping; all accepted with the verifiers' corrections.
  - Blocking: v3's both-arm fallback lets planted donors into C_H's reference inside K3 (4-5 donors per arm), so the
    deterministic leak test would fail with correct code (10/10 synthetic cohorts; probability 0.69-1.00 from the
    disclosed eligibility counts alone) → per-label reference tests (6); the fallback's cross-arm effect is stated and
    C is also reported with arm-only profiles (3.1).
  - Should-fix: a C_A built from the healthy arm passed every K3 bar → a mirror reference test for C_A, K2b, and K3e
    (composition with a CD4 Naive state, reported) (6); a failed statistical bar would have stopped the descriptive
    readouts that section 6 promises → three check states (6, 8); 4.4 lacked the per-arm minimum asked for in review 3
    (4); "undefined" was not operational → definitions, and undefined draws can only widen an interval (3.1, 3.2, 7);
    the handover lacked G04(b) and did not name the order of the counterfactual (7); 5.5's union-versus-random rule
    fires on selection alone (synthetic: rho 0.15-0.31 above a random P95 of about 0.04 with no sub-label composition)
    → an expression-stratified all-gene statistic with a status-relabelling null (0-3% false positives and full
    detection at 5% admixture in the verifier's simulation) (5.5).
  - Minor: seeds for K3b's thinning and 4.7; paired S_k and S_w in the K3 bars and the handling of dropped labels; S_w's
    size range and the margin's cost in detection restated; "1 − C, the non-composition residual"; 12.5(b) quoted;
    scopes and K8's source; one name for the four S_w labels; 4.2(a) excludes the fold's own 402-analogue hits; "only
    composition moves the paper to skeleton C"; "order-dependent" needs a difference of at least 0.10.
  - Not adopted: a |paired C_H| ≤ 0.05 bar for fallback splits and a within-arm fallback that includes i (the first
    fails with correct code in 2 of 10 synthetic cohorts, the second puts i's own expression into its reference); the
    sentence "the fallback can only move C_A towards C_H" (false in the verifier's runs); a 1/6 floor on
    "order-dependent" (it would drop real interactions of 0.12-0.17).
  - Operating characteristics in review 4's synthetic world (15 vs 10, gap about 2 healthy-SD): pure composition gave
    "composition" in 12/12; a uniform state "detected" in 12/12; a shift outside the four labels was falsely "detected"
    in 1/30; a CD4 Naive shift with naive depletion scattered over all categories (order-dependent 3/12, "composition"
    3/12, two of them with a true C_A of about 0.6). Every K3 bar passed in 10/10 synthetic cohorts with correct code.
- 2026-11-02: v1, approved by the user (section 9); tag `prereg-a1-g01-v1`.
