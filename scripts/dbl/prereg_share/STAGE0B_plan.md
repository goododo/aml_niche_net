# Stage 0b: re-score the smoke synthetics for the share design, then planning power (plan v1, 2026-11-12)

Status: PLAN for the user's approval. Nothing is coded or run yet.

0a finished on 2026-11-12 with "K1 usable" (SMOKE_REPORT section 6h). By the redesign (REDESIGN_share_plan.md section 7), Stage 0 still holds one close test:

> **close with BN-1 if planning power at ψ 2 is below 0.80 at n = 12.**

The planning numbers in redesign section 11 used per-patient information assumed from design 2's templates (SE of log ψ_i 0.194 / 0.348 / 0.535). 0b replaces those assumptions with rates measured on the smoke's own synthetic doublets. **It changes exactly one input; everything else in the planning calculation is kept as it was.**

## 1. What exists (checked 2026-11-12)
- **Stage 04's stored per-barcode tables** for the run tags tune2_knn and f3_knn: dbl_typing_injected_devel_<tag>_smoke.csv.gz, held-out runs.
  - Synthetic doublets per smoke lane: T|LEUK 900 (180 per LEUK state); Tamb|LEUK*, NK|LEUK*, B|LEUK*, Ery|LEUK*, MonoDC|LEUK* and LEUK*|LEUK* 60 each; T|Tamb, T|NK, T|B, T|Ery, T|MonoDC, T|LEUK_other and T|T 60 each.
  - Columns used: called, origin, ambiguous, share.
  - The same tables also hold re-entered **real singlets** (type "singlet:*"). 0b drops those rows at read and never flags them (masking, redesign section 8).
- **Stage 10 rebuilds every held-out injection deterministically** (seeds DBL_SEED_SYN + 10f + b). From the rebuild it computes the NNLS T+LEUK call, the T-identity UMI gate and M1. Its job took 4-6 min.
- **There is no per-barcode M1 file.** Only k/n tables exist (dbl_nbs_rates_{tune,confirm}_smoke.csv). D needs M1 and the other flags jointly per barcode, so the rebuild has to be redone. No scDblFinder call is needed.
- **Not injected in the smoke:** Other|LEUK, non-T non-LEUK pairs (MonoDC|MonoDC, B|B, NK|MonoDC ...) and triplets. These enter the planning only as stated assumptions (section 4).
- **Templates are the three smoke patients** from stage 11: AML3948 (Nd 1263), AML0612 (Nd 393) and AML3762 (Nd 175).
  - Class-A lineage compositions come from the calib files.
  - For the P01 donors only count-1 is used. The real analysis pools three lanes, so this is conservative.

## 2. Step A: re-score (new scripts/dbl/12_share_calibrate.R, mode --reuse-smoke)
- **Inputs.** The stage 04 tables above, synthetic rows only; the lane objects and kNN labels (class-A cells only), as in stage 10.
- **Rebuild.** Stage 10's rebuild block is copied verbatim: build_injection, make_basis, nnls_lin, leuk_share and the T-identity gate. The classifier is dropped, because M1 does not use it.
- **Flags per synthetic doublet.**
  - **Y = M1**, stage 10 verbatim: (V1a OR V3) AND T-identity UMI ≥ 1.
  - **D = M1 OR [called AND (an unambiguous origin naming the LEUK lineage plus one X ∈ P, OR the NNLS top two groups = {LEUK, X ∈ P} with minor ≥ 0.15 and gap ≥ 0.05)].**
    - P = {T, Tamb, NK, B, Ery, MonoDC, Other}; LEUK = the 5 states plus LEUK_other.
    - nnls_top2 is stage 10's nnls_lin generalised to return the top two groups (redesign section 9).
  - **Y_w:** Y inside the lane's LEUK-share window, set at 2.5-97.5 % of the lane's synthetic T|LEUK. It is stored for later and does not enter 0b.
- **Validity checks (the output is not used unless all pass).**
  - The rebuilt injection equals stage 04's rows, types and LEUK share, as in stage 10.
  - **M1 k/n reproduces dbl_nbs_rates_tune_smoke.csv and dbl_nbs_rates_confirm_smoke.csv exactly**, row by row for the pair types. This is the known answer for the rebuild.
  - M1 ⊂ D holds for every barcode.
  - No row of type singlet:* enters any flag.
- **Output.** DBL_TAB_DIR/dbl_share_rates_smoke.csv: lane x pair type x truth x flag, as n and k. Plus dbl_checks_12_smoke.csv. No per-barcode file.

## 3. Step B: planning power (new scripts/dbl/13_share_model.R, mode --planning)
1. **Expected counts per template t at ψ** (redesign section 2). Every doublet term is Nd x composition weight x synthetic rate.
   - The weight of pair type a|b is 2·f_a·f_b, or f_a² for a homotypic pair.
   - Y side: S^Y = Σ_k (T|LEUK_k terms) × r^Y; A = Σ (all other pair types) × r^Y + Ns × a_s.
   - D side: S^D and B in the same way, with r^D and b_s.
   - π_t(ψ) = (ψ S^Y + A) / (ψ S^D + B).
2. **Per-patient information.** Y ~ Binomial(D, π). SE(log ψ̂_t) = 1 / sqrt(E[D] · π(1-π) · slope²), where slope = d logit π / d log ψ at ψ = 1.
   - Calibration uncertainty is added by 2,000 draws of every rate from its Beta(k + ½, n − k + ½) posterior, with the variance of log ψ̂ added by the delta method.
3. **Cohort power.** redesign/power_share.R's run() is copied verbatim; the per-template variances replace v_t. Everything else is unchanged:
   - n = 12 cycling the three templates; dataset shares 45/30/25 %; τ 0.3;
   - shared bias: common sd 0.08 plus dataset sd 0.08 / 0.12 / 0.12; equal dataset weights;
   - the binding window and composition approximations and the 2-of-3 dataset rule;
   - 20,000 simulated cohorts.

## 4. Frozen before Step A runs
**Close rule.** BN-1 iff power_binding(n = 12, ψ = 2) < 0.80 under the base assumptions:
- singlet rates a_s = 0.002 and b_s = a_s + 0.003 (the redesign's power-verdict grid point);
- triplet weight 0.10 per doublet;
- the missing-type rules below.
The text if it triggers: BN-1, "lineage-level T-LEUK pairing has no validated null in these data (planning power x at ψ 2, n = 12)".

**Missing-type rules (base).**
- Other|LEUK takes the lane's composition-weighted mean of its measured X|LEUK rates, for each flag.
- T|Other takes the mean of the measured T|X rates.
- Non-T non-LEUK pairs take the mean D rate of the measured T|X pairs, as a proxy for "a pair without LEUK typed as LEUK + X", and Y rate 0.
- Triplets take the lane's pooled rates over all measured doublet types.

**Reported, never deciding.**
- Power at ψ 1.25, 1.5 and 2 for n = 10, 12 and 16; P(C | ψ 1); size at ψ 1.
- Sensitivity rows:
  - the singlet grid (a_s 0.001 to 0.008);
  - missing-type rates ×0.5 and ×2;
  - triplet weight 0.04 and 0.14;
  - P01 templates with three-lane Nd;
  - τ 0.5;
  - the redesign's alternative bias prior (0.10 + 0.10 / 0.15).
- Per template: E[D], E[Y], π0, slope, SE(log ψ̂), and the SE share that comes from calibration.
- A binomial simulation of ψ̂_t per template, checking the normal approximation for small counts.

**Expected outcome (redesign section 11, assumed rates).** Power at ψ 2 is 0.98, so 0b most likely passes. At ψ 1.5 it is 0.67, which gives the verdict "only a doubling is detectable". Measured rates can move either number.

## 5. Known answers (13 refuses to report unless all pass)
- **Regression.** run() with design 2's v_t reproduces the redesign's numbers at n = 12: binding power 0.67 at ψ 1.5 and 0.98 at ψ 2, within ±0.02 (Monte Carlo SE about 0.003).
- **Doublet-rate cancellation.** With a_s = b_s = 0, π_t(ψ) is identical, to 1e-9, under Nd ×0.5 and ×2.
- **Fisher SE.** For one template at E[D] ×100, the binomial-simulation SD of log ψ̂ is within 5 % of the Fisher SE.
- **Size.** At ψ = 1 the primary rejection rate is 0.05 ± 0.01.

## 6. Code, runs, freeze
| File | Change | ~ lines |
|---|---|---|
| config_dbl.R | section 11b: P, the D thresholds 0.15 / 0.05, planning constants (τ, priors, grid point, triplet weight, n, ψ grid, bar 0.80, nsim), seeds SEED + 1071 / 1081 | +15 |
| 12_share_calibrate.R | new; only the --reuse-smoke mode; about 60 lines copied verbatim from 10 | ~160 |
| 13_share_model.R | new; only the --planning mode; run() copied verbatim from power_share.R | ~150 |

- **Files.** 12 and 13 are the files named in the redesign's code plan (section 9). Later stages add their other modes; nothing is restructured.
- **Order.**
  1. Code.
  2. Unit and known-answer tests.
  3. An independent review of both scripts.
  4. One commit by explicit path, tagged prereg-dbl-0b-v1.
  5. `sbatch submit_dbl.sbatch 12 --reuse-smoke`, about 10 min.
  6. `sbatch submit_dbl.sbatch 13 --planning`, minutes.
- **No real barcode is typed.** Only synthetic doublets are flagged; singlet rows are dropped at read.
- **Note on 03b.** The config edit changes config_dbl.R's md5, so 03b's guard will refuse to rerun. 0a is complete, so this is intended.

## 7. After the run
- **Pass:** SMOKE_REPORT section 6i and memory are updated. The next step is the design of freeze 1, i.e. Stage 1 (eligibility, the lane table, the headers of stages 12-14), as its own plan for approval.
  - 0a's lesson goes into freeze 1 as an amendment: K1 uses doublet-like cross calls only, or models the singlet share. Its basis is f_hat 0.43 and G2 at 1.433.
- **Fail:** BN-1, and the close-out in SMOKE_REPORT section 11.

## 8. Review amendments (2026-11-12, before the tag; no stage-12 rate has been computed)
One independent reviewer ran only `13 --planning --mock`; 12 was never run. The mock's D is fake: D = 1.5 x M1 + 1, built from stage 10's M1 counts.

**Changes adopted.**
- **F1, triplets.** A triplet is now flagged at the donor's composition-weighted rate per doublet (S + O).
  - The literal "lane's pooled rate over all measured synthetic doublets" imported the smoke's allocation, which is 54 % T|LEUK. Triplets then got roughly the T|LEUK recall, inflating A and B.
  - In the mock this one change moved power at ψ 2 from 0.71 to 0.81.
- **F2, calibration draws.** These are a joint posterior: r_D ~ Beta(k_D + ½, n − k_D + ½) and r_Y = r_D × Beta(k_Y + ½, k_D − k_Y + ½), because Y ⊂ D on the same doublets. Section 3 said "every rate from its Beta", which is corrected here.
- **F4, md5 guard.** 12 and 13 now stop if config_dbl.R or config_paths.R differ from the literals taken at the tag. Both record the md5 of their running copy.
- **F5, reporting.** The templates and power tables are written and printed only if every known answer passes.
  - Rerun rule: one rerun, only for a demonstrated code bug, with the same seeds.
- **F6, ψ̂ at the upper end.** A share at or above S^Y / S^D now maps to the upper clip; it used to map to the lower one. D = 0 gives NA.
- **F7, Other|LEUK.** It is weighted by **donor** composition, as the code does; section 4's "lane" is corrected.
- **F8, missing types.** Their rates are capped at 1. The ×0.5 / ×2 rows scale only the three imputed types, not triplets.
- **F9, reported row.** A non-deciding row now re-solves ψ̂ per draw (brute force). The delta method understates calibration variance where the slope is small.
- **F10, housekeeping.**
  - 12 checks that every flagged barcode is one of stage 04's synthetic doublets.
  - A killed fork stops 12.
  - The template table carries the total SE and the calibration share of the variance.

**F3, the calibration basis of the base rule: the user decides (config DBL_0B_CAL_BASIS).**
- "design": k/n is held at the smoke's rates, with n at stage 12's planned sizes per donor-lane, halved for P01. Those sizes are T|LEUK 150 per label, X|LEUK 250 per type, LEUK|LEUK 800 and T|X 100 per type.
- "smoke": the smoke's own counts, about 60 per non-T|LEUK type per lane.
- Whichever basis is not chosen is reported. In the mock, calibration variance alone moved power at ψ 2 across the bar (0.71 vs 0.85 with sampling variance only).
