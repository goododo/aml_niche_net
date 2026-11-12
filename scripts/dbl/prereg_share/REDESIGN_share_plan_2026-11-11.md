# Lineage-level T-LEUK test: redesign (workflow wf_f6844a8c-01b, 2026-11-11)

Status: PROPOSAL, not approved. Nothing coded or run. Saved here because the session scratchpad is not durable.

# Recommended plan: a dbr-free within-patient share, checked on real cross-donor doublets

**The idea.** Stop counting T–LEUK doublets in absolute terms. Inside each patient, test whether T–LEUK calls make up a larger share of LEUK-partner doublet calls than random co-encapsulation of that patient's own cells predicts. This is design 2, from v1.1 Q1's share family. Every doublet type carries the same doublet count, so the count cancels. Hashtag cross-donor doublets stay, but only as a real-doublet check (from designs 1, 3 and 4), used on the side where leakage cannot bias them. A one-day untyped Stage 0, with a close rule written in advance, comes first.

## 1. Estimand
- **Per patient i:** ψ_i = R_i(T–LEUK) / R_i(ref).
  - R(u) = observed / randomly expected same-donor doublets of type u.
  - ref = the D-weighted mix of LEUK|X pairs, X ∈ {Tamb, NK, B, Ery, MonoDC, Other}.
  - Lineage level: strict T; LEUK = 5 states + LEUK_other.
- **Cohort:** log ψ_i = log ψ + δ_dataset + u_i, with u_i ~ N(0, τ²).
  - δ_dataset is sum-to-zero over the datasets, and is dropped if any dataset has fewer than 2 eligible patients.
  - One common ψ for adults and children.
  - H0 ψ = 1 vs ψ > 1, one-sided (v1.1 Q1); band [2/3, 1.5].
- **What a positive result says:** T cells are over-represented among leukemic cells' doublet partners, beyond random co-encapsulation. Never "contact".
  - Cancels: the doublet rate, generic stickiness, and proportional typing loss.
  - Pulls ψ down: an excess of NK–, Tamb– or MonoDC–LEUK pairs.
  - Not excluded: T-generic stickiness (see the ψ_T qualifier, §5) and aggregation during thawing.

## 2. Normalisation: why the doublet rate drops out
- **Flags.**
  - Y = M1, unchanged.
  - D = M1 OR [called AND (unambiguous origin {LEUK-lineage, X ∈ P} OR NNLS top two {LEUK, X ∈ P}, minor ≥ 0.15, gap ≥ 0.05)], with P = {T, Tamb, NK, B, Ery, MonoDC, Other}.
  - M1 ⊂ D is asserted for every barcode.
- **Expected share:** π_il(ψ) = (ψ·S^Y + A) / (ψ·S^D + B). Every doublet term is Nd_il × composition weight × synthetic rate, so Nd cancels exactly, whether assumed or true.
- **Where the rate still enters.** Only in two second-order terms:
  - the singlet weight ρ = Ns/Nd: d log ψ̂ / d log ρ ≤ 0.2 for the worst template, ≈ 0.03 information-weighted;
  - the within-donor triplet weight λp/3: 0.04–0.14 across lanes, modelled with synthetic triplets.

  Pools take λ from hashtag-counted multiplets. Other lanes use the 10x rule on in-QC barcodes, not 8e-6 × deposited, which is several-fold wrong in debris-heavy lanes such as AML0114. Both terms get a per-lane LN(0, 0.5) factor in the reference distribution.
- **scDblFinder's own dbr input** stays at M1's rule (8e-6 × deposited, cap 0.4) in both calibration and real runs, so rates are measured at the real operating point. Gate: |threshold(calibration run) − threshold(masked real run)| ≤ 0.05.
- **Known answers.**
  - With singlet rates set to 0, expected shares are identical to 1e-9 under dbr ×0.5 and ×2.
  - Null worlds with true dbr ×0.5 and ×2 must pass the size bar.

## 3. Data scope and honest count
- **GSE185381.**
  - Adults: AML0612 and AML3762 (P01, 3 hashed lanes each); AML4340 (777/78, unnamed tags); AML3948 (785; a second lane only if its barcodes are disjoint); AML0114 (3 debris-heavy lanes); AML0693.
  - Children: AML005 and AML055 (764/65); AML073 and AML025 (762/63); AML009, which sits in a genotype-demultiplexed pool, so only its author-assigned A/B barcodes are used.
- **GSE289435:** 6 patients. GSM8791437 (AEL) is excluded; GSM8791439 is MDS.
- **Petti 2019:** 5 patients.
- **Count.** 22 candidates (17 adults).
  - Realistic after strict T, lane typing gates (5′ chemistry is untested) and floors: 11–18, central ≈ 14 (9–14 adults). Rough P(n ≥ 10) ≈ 0.7–0.8.
  - Likely losses: Ery-rich patients (Petti 809653, GSM8791433, GSM8791444, where expected T–Ery ≥ T–LEUK) on the precision gate; AML0114 and AML0693, whose dbr input lies outside M1's validated range.
- **Pools.** A donor's analysis set is every line-QC-passing barcode that our hashtag rule calls single(i), in any class.
  - Two-tag barcodes are random doublets formed after pooling, so they are excluded; AML3762's two-tag share is reported as a null control.
  - Escapes (cross-donor doublets that lost a tag) are modelled from synthetic pairs with summed hashtag rows, per compartment pair.
- **K1 lanes:** the cross-donor barcodes of the analysis pools, plus GSM5613756 and 760/61.
- **Excluded:** the T-ALL co-captured pools 747–755, and GSE201966.

## 4. Calibration (stage 12; only synthetic or untyped data before freeze 2)
- **Runs.** Per lane, 20 held-out runs: 5 folds × 4 batches, removing one fifth of class A per run as in the smoke. Devel 1.27.12 with proc_320 and kNN labels. Injections are knownDoublets with knownUse "discard", at most 2.5 % of the lane per run.
- **Synthetic set per qualifying donor-lane** (design 2's allocation); halved for donors spanning 3 lanes:

  | Type | Count |
  |---|---|
  | T\|LEUK | 150 per LEUK label |
  | X\|LEUK (incl. LEUK_other partners) | 1,500 |
  | LEUK\|LEUK | 800 |
  | T\|X | 700 |
  | non-T non-LEUK (MonoDC\|MonoDC, B\|B, NK\|MonoDC …) | 700 |
  | triplets | 300 |
  | cross-donor pairs with summed HTO rows (pools only) | 900 |
  | T+ambient and T+platelet pseudo-singlets (window check) | 100 each |
- **Rate tables.** k/n per lane × pair type for Y, D, the T-partner flag and the window variant Y_w, plus 2,000 parent-cluster bootstrap replicates of every rate.
- **One masked real run per lane.** It returns only its threshold, the untyped called fraction and the cluster-28 called fraction.
- **Compositions.**
  - f_A: class A.
  - f_Q: class A plus author-assigned barcodes that failed only the lower QC tails, with lineage-level labels.
  - π0 is computed under both.
- **Hashtags (untyped).** Tag maps (author majority ≥ 0.80 for unnamed tags), per-compartment HTO sensitivity, donor-pair propensities, and λ̂ per pool.
- **Replica known answer** (v1.1 §4.7.2), on the most T-rich eligible lane of each dataset.
  - Lane = class-A cells only, plus synthetic random same-donor doublets added as ordinary kNN-labelled barcodes at the lane's Nd, plus planted T|LEUK at ψ ∈ {1, 1.5, 2}.
  - Only synthetic barcodes are counted.
  - It also measures how the share transfers from discard mode to ordinary mode.
- **Singlets are not calibrated here (masking).** See U0 in §7.

## 5. Analysis model (stage 14, after freeze 2)
- **Likelihood.** Y_il | D_il ~ Binomial(D_il, π_il(ψ_i)).
  - Expectations are computed per lane, and the likelihood is multiplied over a patient's lanes (A.3).
  - u_i is integrated by 20-node Gauss–Hermite, checked against brute force.
  - ML by L-BFGS-B with |log ψ| ≤ 3. A fit that ends on a bound counts as not identified, giving D.
- **Test.** One-sided LRT, with 1,999 parametric draws under the fitted null. Each draw redraws:
  - u_i, with τ ≥ 0.15;
  - a shared relative bias: a common N(0, 0.08²) plus a per-dataset N(0, σ_d²), with σ_d = 0.08 for GSE185381 and 0.12 for GSE289435 and Petti;
  - lane rates, from the bootstrap;
  - singlet rates, from U0's posterior × a per-lane LN(0, 0.5);
  - the ρ and λ factors.

  Each draw is then refitted with the point rates.
- **Interval.** 999 draws under the fitted alternative, in stage 11's studentised form.
- **Binding variants.**
  - Y_w: Y restricted to the lane's LEUK-share window, set at 2.5–97.5 % of synthetic T|LEUK, with all rates recomputed. This restores v1.1 §4.5.6 as a guard.
  - Conservative composition: per patient-lane, whichever of f_A and f_Q gives the larger π0, fixed before freeze 2.
- **Hidden doublets in class A** inflate the singlet rate at roughly the doublet mix. This barely moves π0, and is conservative under H1.
- **K1 = b̂0**, the cross-donor share's O/E at ψ_x = 1. It uses the donor-pair composition with propensities, the cross-donor synthetic rates, and terms for false cross-donor calls on singlets and for triplets.
  - Pre-staining T–LEUK units leak into the cross-donor set: in P01 at ψ 1.5, ψ_x ≈ 1.15–1.35. So b̂0 is clean only under H0.
  - K1 therefore gates the null side only, and cannot certify A or B in any dataset.
  - A/B protection comes from the propagated bias priors, the binding variants, the replica and the 2-of-3 dataset rule.
- **Qualifier ψ_T:** the same model with the T-partner flag as denominator.
- **Reported, never deciding:**
  - per-patient ψ̂_i and per-dataset ψ̂;
  - adults and children separately;
  - leave-one-dataset-out;
  - without MonoDC; without Tamb + NK; without class C;
  - v1.1's paediatric-covariate model;
  - b̂0 against its leakage prediction at ψ̂;
  - AML3762's two-tag share;
  - LEUK-share distribution of real M1 barcodes against synthetic ones;
  - Petti: mutant reads among M1 barcodes against the predicted precision.

## 6. Simulation worlds (stage 13, realised templates)
- **Singlet grid.** a_s ∈ {0.001, 0.002, 0.004, 0.008} and b_s − a_s ∈ {0.0015, 0.003}.
- **Gating nulls:** the base world; true dbr ×0.5 and ×2; τ 0 and 0.5; n = 10; singlet rate at the grid maximum.
- **Gating robustness** (false A/B under the full rule):
  - T under-represented by 15 % in class A, with f_Q recovering half;
  - non-class-A singlets at 3×, the extras T–LEUK-typed, their LEUK share drawn from the T+ambient pseudo-singlets;
  - shared transfer ratio T|LEUK ×1.10 and ×0.90;
  - NK|LEUK and Tamb|LEUK ×1.5, including their leak into Y; Tamb alone ×1.5;
  - T|MonoDC ×2;
  - every pair type ×1.5 (ψ̂ must stay ≈ 1);
  - escape ×1.5 for T-containing pairs;
  - triplet weight ×1.5;
  - ρ ×3 in debris-heavy lanes;
  - calibration effective n halved.
- **Reported stress (not gating):** transfer ×1.15 and ×1.20; singlets ×6; shared bias 0.15.
- **Power:**
  - ψ 1.25, 1.5 and 2 at τ 0.3;
  - ψ 1.5 at τ 0.5;
  - P(C | ψ 1.5) in the base world and under transfer ×0.90 and ×0.85.

  Power always means P(outcome A or B) under the full rule.
- **Edge cases:** n = 9 gives all D without a fit; a lane with D = 0; an absent LEUK state; a compartment with fewer than 5 cells; a patient without reference partners. The null world must reach the no-result branch.
- **Function checks:**
  - Gauss–Hermite against brute force, to 1e-5;
  - 400 Poisson patients at ψ 2 give |ψ̂ − 2| < 0.05;
  - the LRT is 0 when the data sit below expectation;
  - decide() unit cases;
  - doublet-rate cancellation to 1e-9.

## 7. Reading rule (frozen at freeze 1)
- **Stage 0 → close with BN-1** if either:
  - planning power at ψ 2 is below 0.80 at n = 12; or
  - K1 is unusable on P01 count-1, meaning any of: pair-matrix χ²/df > 2 after donor propensities; median nUMI(cross)/nUMI(single) < 1.4; or two-tag observed/predicted (doublet-like nUMI only, after propensities) outside [0.67, 1.5]. The authors' metadata on what the second tag marks is checked too.

  v1.1 §4.8 stops the line on a failed known answer.
- **Stage 1 (after freeze 1):** strict T ≥ 30 and LEUK ≥ 50 in class A, summed over a patient's lanes. Fewer than 10 patients → BN-1.
- **Stage 2 lane gates:**
  - doublet-level natural precision of M1 ≥ 0.70 and natural recall ≥ 0.50;
  - D purity ≥ 0.50;
  - SE(logit π0) ≤ 0.10;
  - threshold match ≤ 0.05;
  - the window removes ≥ 90 % of flagged T+ambient and T+platelet pseudo-singlets and keeps ≥ 95 % of synthetic T|LEUK.

  Final eligibility adds predicted Y ≥ 5, predicted LEUK|X ≥ 5, and slope |d logit π / d log ψ| ≥ 0.2. Fewer than 10 patients → BN-1.

  Replica: ψ̂ must be monotone and within ×/÷ 1.2 of each planted ψ, pooled over replicas. A failure stops the line.
- **Stage 3 machinery.** Otherwise the run is not read → BN-2. One rerun is allowed, only for a demonstrated code bug and with the same seeds. Settings: 1,000 cohorts per world, 199 reference + 199 interval draws each.
  - Gating null size ≤ 0.071 (0.05 + 3 SE, SE 0.0069). P(false fail) is 0.002 per world, ≈ 0.02 over 8 worlds.
  - Robustness false A/B ≤ 0.128 (0.10 + 3 × 0.0095).
  - 95 % coverage ≥ 0.90 at ψ 1 and 1.5. SE 0.007, so a true 0.95 cannot fail.
  - Mean log ψ̂ at ψ 1.5 within log 1.5 ± 0.10 (Monte Carlo SE ≈ 0.005).
  - P(C | ψ 1.5) ≤ 0.071.
  - Empty world all D; no errors; non-converged fits count as D and stay ≤ 1 %.
- **Power verdict at freeze 2** (at a_s 0.002):
  - power(ψ 1.5) ≥ 0.80 → "adequately powered";
  - else power(ψ 2) ≥ 0.80 → "only a doubling is detectable" (decision 12);
  - else BN-2, without unblinding.
- **Stage 4.**
  - **U0.** Singlet rates are calibrated inside the worker: 10 discard-mode runs per lane, each re-entering a disjoint tenth of class A. A lane whose rate exceeds 0.008 has its patients made descriptive, and n ≥ 10 is rechecked.
  - **U1.** If K1's 95 % CI lies entirely below 0.80 → D (BN-3).
  - **U2 outcomes (A.5):**
    - A: p < 0.05; lower(ψ) > 1.5; both binding variants reject; ψ̂_d > 1 in at least 2 of 3 datasets (v1.1 §4.8).
    - B-small / B-undetermined: as A, without the 1.5 bound. B-small if the interval lies inside the band.
    - C: p ≥ 0.05; interval ⊂ [2/3, 1.5] under the primary, under Y_w and under both compositions; and K1's CI ⊂ [0.67, 1.5].
    - D: everything else, including a rejection that a binding variant does not reproduce. If fewer than 2 datasets have ≥ 2 eligible patients, A/B is capped at B-undetermined.
    - "T–LEUK-specific" may be written only if ψ_T also rejects.
- **Pre-written bounded negatives:**
  - BN-1: "lineage-level T–LEUK pairing has no validated null in these data (n < 10, or the cross-donor check failed)".
  - BN-2: "at n = k the share test detects ψ = 1.5 with power x only".
  - BN-3: "synthetic calibration does not transfer to real random doublets (b̂0 = x)".

## 8. Masking and freezes
- **Before freeze 2**, no typed call of a real barcode is written or printed. That covers within-donor, cross-donor, two-tag, hashtag-negative and re-entered singlet barcodes, in any split or in total.
- **Allowed before freeze 2:**
  - barcode classes;
  - compositions (after freeze 1);
  - hashtag calls, loads and pair matrices;
  - synthetic barcodes' flags;
  - masked real-run thresholds and untyped called fractions;
  - replica ψ̂ (synthetic barcodes only);
  - simulations.
- **Freeze 1**, committed before Stage 1 touches a new lane.
  - Contents: this design; config §11; the headers of stages 12–14 with every rule; the world list and bars; the lane table.
  - Amendments: the share replaces the absolute O/E; one common ψ for children; recall floor 0.50; the window restored as a binding variant.
  - Exposures: v1.1 App. A and §2.1; stage 11's templates; the 2 typed AML3948 singlets of 11-09; the 10-28 hashtag tally; P01 count-1 hashtag counts and called fractions; all Stage 0 outputs.
- **Freeze 2**, after stage 13.
  - Contents: the eligible list, gates, rates, templates and verdict.
  - Record: git rev-parse HEAD:scripts/dbl (a tree hash, so the 23:30 cron does not break it), plus the md5 of the parameter table and of stage 14.
  - Stage 14 refuses to run otherwise, and runs U0 → U1 → U2 in that order.
- **No third exception.**

## 9. Code plan (scripts/dbl/; each step is planned, approved, then coded)

| File | Change | ≈ lines |
|---|---|---|
| config_dbl.R | §11 constants: P, flags, allocation, gates, σ, τ floor, singlet grid, worlds, bars, seeds SEED + 921/931/941/951, freeze slots | +70 |
| dbl_lanes_share.tsv | lane table: dataset, GSM, donors, hashed, role | ~30 rows |
| 01 / 02 | read the lane table; GSE289435/Petti ingest layouts; second-lane overlap check | +40 / +15 |
| 03_call_hashtags.R | loop over hashed lanes; unnamed tags; `--desk` mode (Stage 0); λ̂ and propensities. The caller itself is unchanged | +90 |
| 12_share_calibrate.R | `--reuse-smoke` (Stage 0: re-score stored smoke synthetics, no scDblFinder call) and `--calibrate` (runs, flags, window, replica, masked real run). call_scdbl, proc_320 and build_pairs (04) and make_basis/nnls (09/10) copied verbatim; nnls_lin generalised to nnls_top2 | ~480 |
| 13_share_model.R | `--planning` / `--audit`: shares, Gauss–Hermite GLMM, LRT, reference, worlds, bars, verdict, freeze-2 table. decide() and chk() copied from 11 | ~450 |
| 14_share_unblind.R | freeze guard, U0/U1/U2, outputs; 12's worker and 13's functions copied verbatim | ~350 |
| submit_dbl.sbatch | lane-array variant, c = 16, OMP_PROC_BIND=false | +10 |
| smoke/code/check_masking_share.R | stub-scDblFinder harness: no real typed column before freeze 2; 14 refuses without the record | ~60 |

- About 1,600 lines, roughly 40 % copied.
- Per-lane outputs are merged in a separate step, so no array task writes a cohort table.
- Jobs run dated read-only copies of the scripts.
- Outputs go to the smoke tables folder with a `_share` suffix.

## 10. Compute

| Stage | CPU-hours | Notes |
|---|---|---|
| 0 | < 3 | ≈ 1 analyst day |
| 1 | ≈ 5 | one load per ingest object, ≤ 40 GB |
| 2 | ≈ 40–60 | ~28 lanes × ~23 devel runs × 2–5 min, as a 16-core array |
| 3 | ≈ 30–60 | ~25 worlds × 1,000 cohorts; 2–4 h wall |
| 4 | ≈ 20–25 | includes U0 at ≈ 12–20 CPU-h |

Total ≈ 100–150 CPU-hours, 4–6 days wall-clock with queues, and ≈ 2–3 weeks of side-line analyst time over five approved steps.

## 11. Planning numbers
These are summary-level, using design 2's assumed D rates; Stage 0 replaces them.
- **Power, with the bias priors, τ 0.3 and every binding condition:**
  - n 12: 0.67 at ψ 1.5 (0.78 for the primary test alone) and 0.98 at ψ 2.
  - n 16: 0.75 at ψ 1.5 and 0.99 at ψ 2.
- **Other outcomes:**
  - Size 0.05 for the primary test, 0.024 under the full rule.
  - P(C | ψ 1) is 0.37 at n 12 and 0.52 at n 16.
  - P(A | ψ 2) is 0.40–0.47.
- **Bias beyond the prior.** A shared bias of 0.10 in every dataset gives false A/B 0.07; 0.15 gives 0.14. At τ 0.5, power at ψ 1.5 falls to 0.50–0.58.
- **Expected verdict:** "only a doubling is detectable", unless the measured templates beat the planning ones. Unlike every absolute design, outcome C is reachable.
- **Scratch scripts (not durable):**
  - /tmp/claude-39170/-FAST-gr10634-gaozy/24ae6318-b751-4c66-856d-38903ce24d07/scratchpad/synth/power_share.R
  - /tmp/claude-39170/-FAST-gr10634-gaozy/24ae6318-b751-4c66-856d-38903ce24d07/scratchpad/synth/arith1.R

## Decisions for the user (synthesis defaults)

- 1. Now: run Stage 0 only, then decide again. Stage 0 = 12 --reuse-smoke, 03 --desk and 13 --planning: about 1 analyst day, under 3 CPU-h, no real typed count. Default: yes.
- 2. Stage 0 close rule: close with BN-1 if planning power at ψ 2 is below 0.80 at n = 12, or if the P01 cross-donor set fails the K1 usability checks. v1.1 §4.8 stops the line on a failed known answer. Default: adopt as written.
- 3. Primary estimand: the within-patient LEUK-partner share ψ (v1.1 Q1 family) replaces the absolute O/E of v1.3 A.7, recorded as an amendment at freeze 1. Default: yes. Not recommended: a hashed-only absolute O/E with a cross-donor baseline (about 7 patients, diluted by leakage).
- 4. Children: one common ψ, which amends v1.1 §4.8's paediatric covariate; adults-only, children-only and the covariate model are reported; n ≥ 10 is counted over all eligible patients. Default: yes. Alternative: v1.1's covariate model with n counted on adults; 9–14 adults are expected, so D comes up more often.
- 5. Keep n ≥ 10 and the no-third-exception policy: a machinery rerun only for a demonstrated code bug, with the same seeds. Default: keep.
- 6. Per-lane typing gate: doublet-level M1 precision ≥ 0.70 and recall ≥ 0.50. Once calibrated, low recall only costs power. Default: yes. Alternative: 0.70/0.70, which loses more lanes.
- 7. Shared relative-bias prior in the reference distribution: common sd 0.08 plus a per-dataset sd of 0.08 (GSE185381) and 0.12 (GSE289435, Petti). Default: yes. Alternative: 0.10 plus 0.10/0.15, safer against transfer bias but about 0.1 lower power at ψ 1.5.
- 8. Singlet calibration masking: strict. Singlet rates are computed only inside the unblind worker (U0, about 12–20 CPU-h extra), and stage 13 uses a pre-registered grid. Default: strict. Alternative: write sealed counts pooled per dataset before freeze 2, as a registered exposure.
- 9. Binding variants: the per-lane LEUK-share window (v1.1 §4.5.6) and the conservative composition must both reproduce any A, B or C. Default: yes; costs about 0.1 power.
- 10. Partner set P = {T, Tamb, NK, B, Ery, MonoDC, Other}, with 'without MonoDC' and 'without Tamb+NK' reported. Default: yes. Alternative: drop MonoDC from the primary, as v1.1 kept T–MonoDC out.
- 11. AML009 (genotype-demultiplexed pool): include it using author-assigned A/B barcodes only, with an exclusion sensitivity. Default: include.
- 12. At freeze 2, if the verdict is 'only a doubling is detectable': proceed to unblinding (default), or stop with BN-2.

## Rejected and why

- hto-dbr as primary (absolute O/E, hashtag-measured dbr, cross-donor anchor): pre-staining units leak into the anchor at P01 loading. Power is 0.26–0.36 at R 1.5, reject-and-gate is 0.38–0.59 at R 2, n ≤ 7, and P01 holds about 86 % of the signal.
- skeptic-minimal as primary (absolute O/E, hashtag-counted N_dbl, keeps n ≥ 10): its own rules leave 4–7 patients, so every path ends in D, and P01's two-tag/cross ratio of 1.75–1.82× already fails its G2.
- empirical-baseline as primary (cross-donor rows as the baseline): cancels the doublet rate exactly, but has about 7 patients in 4 pools with AML0612 at about 57 % of the weight. Its κ rows also absorb leaked signal, and per-donor singlet pre-selection gives 0.21–0.29 false rejection. Kept only as the K1 idea.
- Stage 11's absolute O/E with an assumed dbr: power is capped near 0.63 at any n, and the rule is already about 1.3× off in P01.
- Transporting the hashtag-measured dbr to unhashed lanes, or estimating it from called totals: size reaches 0.18 at a 1.3× structural error, and called totals bring back scDblFinder's dbr input.
- Cross-donor clearance for A/B (lower ψ > upper ψ_x, as design 2's critic asked): leakage puts ψ_x at about 1.15–1.35 when ψ = 1.5, so true signals would fail it. A leakage-corrected joint fit identifies ψ only through a contrast of about 0.15 log units.
- Full 2×2 interaction (log odds ratio including non-T non-LEUK doublets): cancels composition, but needs X–X typing that is sparse and unvalidated.
- T-partner share as primary: its small denominator makes it more sensitive to Y-skewed singlets, and T–Ery false positives (0.45–0.62) weaken the contrast. Kept as a wording qualifier.
- Re-tuning M1 or its dbr input, or a fixed 0.25–0.80 window as primary: would be a third exception. The per-lane window is kept only as a binding guard.
- Writing pooled re-entered singlet counts before freeze 2 (the default in designs 2 and 3): in single-donor lanes this is a typed per-patient count. Replaced by calibration inside the unblind worker.
- Lowering the 10-patient minimum: would be a third exception, informed by eligibility numbers that are already exposed.
- Closing now: stays as the Stage 0 fallback (BN-1). One untyped day decides it with a pre-written rule.

## Biggest risks

- Eligible n < 10 → BN-1 (rough P 0.2–0.3). Causes: strict T, the precision gate on Ery-rich patients, the debris-heavy AML0114/AML0693 lanes, and M1 untested on 5′ chemistry, where TRAC/TRBC capture is weak and the gate rests on CD3D/E/G.
- Shared relative bias beyond the prior, from the synthetic-to-real transfer ratio of T|LEUK vs X|LEUK or from class-A vs loaded composition. A bias of 0.10 in every dataset gives false A/B 0.07; 0.15 gives 0.14. Only GSE185381 has a real-doublet check, and only on the null side.
- Y-skewed singlet false positives: ambient-heavy T singlets, and T–platelet complexes typed T|LEUK_other because Mk maps to LEUK_other. The window variant is the guard, but singlet rates outside class A are never measured directly.
- K1 may be unusable. The unexplained two-tag excess (1.75×) may reflect donor stickiness or aliquot aggregation; if so, the line closes at Stage 0.
- Relative estimand. NK/Tamb/MonoDC–LEUK biology moves ψ. T-generic stickiness is excluded only if ψ_T agrees. Aggregation during thawing and FACS doublet exclusion push ψ toward 1. A positive result is never 'contact'.
- Likely end state: 'only a doubling is detectable' (planning power 0.67–0.75 at ψ 1.5), after about 100–150 CPU-h and 2–3 weeks of side-line time.
- A machinery failure ends the line (BN-2); stage 11 failed twice, once by 0.01. The bars here sit at 3 SE, so Monte Carlo noise alone cannot fail them, but a small real excess can.
- M1 passed through a post-hoc exception with thin margins (worst 0.71–0.73) on 3 patients. New lanes are calibrated on synthetic doublets only, and T-poor lanes (48–200 strict T) force parent reuse.

## Scores (1-5)

| design | dbr robustness | power | coverage | confound control | cost | integrity |
|---|---|---|---|---|---|---|
| hto-dbr | 3 | 1 | 1 | 2 | 3 | 2 |
| conditional-share | 5 | 4 | 4 | 3 | 2 | 3 |
| empirical-baseline | 5 | 2 | 2 | 3 | 3 | 3 |
| skeptic-minimal | 3 | 1 | 1 | 2 | 4 | 4 |
