# Preregistration, freeze 1: the within-patient T-LEUK share test (doublet-pairing redesign)

Status: FREEZE 1, the binding layer. Tag `prereg-dbl-share-f1`.
- The tag is placed before stage 02 runs on any non-smoke lane, i.e. before any class-A composition of a new lane is tabulated in this line.
- Stage 01 extracted the lane objects first. It tabulates barcode classes only, never cell types.

**Precedence.** This file comes first. Next is `REDESIGN_share_plan_2026-11-11.md` (a verbatim copy in this folder), then the stage plans (`STAGE0A_plan_v2.md`, `STAGE0B_plan.md`, `STAGE1_plan.md`, also copied here). Anything not restated below holds as written in the redesign copy.

**Decisions.** The user adopted the redesign's synthesis defaults, decisions 3-12, on 2026-11-12 ("全部采用").

## 1. Estimand and reading (redesign sections 1, 7)
- **Per patient:** ψ_i = R_i(T–LEUK) / R_i(ref). R(u) is the observed / random-expected count of same-donor doublets of lineage pair type u. ref is the D-weighted mix of LEUK|X, X ∈ P = {T, Tamb, NK, B, Ery, MonoDC, Other}.
- **Cohort:** log ψ_i = log ψ + δ_dataset + u_i.
  - δ is sum-to-zero; a dataset with fewer than 2 eligible patients loses its δ.
  - u_i ~ N(0, τ²) with τ ≥ 0.15.
  - One common ψ for adults and children (decision 4).
- **Test:** H0 ψ = 1 vs ψ > 1, one-sided; band [2/3, 1.5].
- **Outcomes:** A / B-small / B-undetermined / C / D, and the bounded negatives BN-1 / BN-2 / BN-3, exactly as in redesign section 7.
- **Wording.** A positive result is never "contact". "T–LEUK-specific" may be written only if ψ_T also rejects.

## 2. Candidates and lanes (Stage 1 plan section 2, S1a)
- **The lane table** is `scripts/dbl/dbl_lanes_share.tsv` (commit 7f22eb3), built by metadata only.
- **Candidates:** 51 patients.
  - GSE185381: 36 AML diagnosis patients, 3′ (17 adults, 19 children);
  - GSE289435: 10, 5′;
  - Petti 2019: 5, 5′.
- **Excluded:**
  - the T-ALL co-capture pools GSM5613747-55;
  - the timepoint-mixed pool GSM5613757-59 (PAUMTZ, PAXMIJ);
  - GSE289435's AEL and MDS samples.
- **Amendment to redesign section 3.** Its 22-candidate list had been screened on compositions seen in earlier work, so it is replaced by the metadata rule. Composition screening is the Stage 1 gate only.
- **Analysis set per lane.**
  - Hashed pools: line-QC-passing barcodes that stage 03's rule calls single(i), in any class (redesign section 3).
  - Multiplexed lanes without hashtag rows: author-assigned barcodes only.
  - Single-donor lanes: every line-QC-passing barcode.
- **K1 lanes:** the cross-donor barcodes of hashed analysis pools (`k1_lane` in the table).

## 3. Labels and flags
- **Compartments.** Stage 02's label map (`dbl_label_map.tsv`) with the strict-T rule: bmm_prob ≥ 0.5, else LEUK_other / Tamb. LEUK = the five states plus LEUK_other. B/C barcodes take kNN labels from the same lane's class A.
- **Y = M1**, stage 10's confirmed typing (V1a OR V3, AND T-identity UMI ≥ 1), unchanged and never re-tuned.
- **D = M1 OR [called AND (an unambiguous origin naming the LEUK lineage plus one X ∈ P, OR the NNLS top two groups = {LEUK, X ∈ P}, minor ≥ 0.15, gap ≥ 0.05)].** This is implemented in `12_share_calibrate.R` (`lx_origin`, `nnls_lx`). M1 ⊂ D is asserted for every barcode.
- **Y_w** = Y inside the lane's LEUK-share window (2.5-97.5 % of synthetic T|LEUK). It is a binding variant.
- **scDblFinder:** devel 1.27.12 with proc_320, and the dbr input at M1's rule (8e-6 × deposited, cap 0.4) in calibration and real runs alike.

## 4. Stage gates
**Stage 1** (this freeze → stage 02 → eligibility):
- strict T ≥ 30 and LEUK ≥ 50 class-A cells, summed over a patient's **analysis** lanes in the lane table (excluded lanes never count);
- **LEUK here is the lineage-level definition: the five states plus LEUK_other.** The smoke's qualification rule (stages 08-11) counted the five states only; that rule is not used here.
- The count is made by `00_share_eligibility.R` from stage 02's class-A composition tables.
- fewer than 10 such patients → BN-1 "lineage-level T–LEUK pairing has no validated null in these data (n < 10)".

**Stage 2** (12 --calibrate; redesign sections 4 and 7). Lane gates, each computed on synthetic doublets:
- natural M1 precision ≥ 0.70 and natural recall ≥ 0.50 (decision 6);
- D purity ≥ 0.50;
- SE(logit π0) ≤ 0.10;
- threshold match ≤ 0.05;
- the window removes ≥ 90 % of the flagged T+ambient and T+platelet pseudo-singlets and keeps ≥ 95 % of synthetic T|LEUK.

Final eligibility adds predicted Y ≥ 5, predicted LEUK|X ≥ 5 and |d logit π / d log ψ| ≥ 0.2; fewer than 10 → BN-1. The replica known answer must give ψ̂ monotone and within ×/÷ 1.2 of each planted ψ; a failure stops the line.

**Stage 3** (13 --audit): the machinery bars and the power verdict at freeze 2, as in redesign section 7. 1,000 cohorts per world; 199 reference and 199 interval draws.

**Stage 4** (14_share_unblind.R): U0, then U1, then U2, as in redesign section 7, with the K1 amendment in section 5.

## 5. Amendments from Stage 0 (binding)
1. **K1 uses doublet-like cross calls only** (from 0a; SMOKE_REPORT section 6h).
   - Doublet-like means nCount ≥ the lane's q75 of class-A single calls, the 0a definition.
   - Reason: f_hat 0.43 and G2 at 1.433 indicate that many raw cross calls look like singlets.
   - b̂0 (U1) and every K1 statement are computed on that restricted set.
2. **Calibration variance in planning uses stage 12's design sizes** (from 0b; user decision 2026-11-12).
3. **The pre-registered 0b result stands:** planning power 0.985 at ψ 2 and 0.669 at ψ 1.5 (n = 12). The expected end state is "only a doubling is detectable".

## 6. Masking and freezes (redesign section 8)
- **Before freeze 2**, no typed call of a real barcode is written or printed. That covers within-donor, cross-donor, two-tag, hashtag-negative and re-entered singlet barcodes, in any split or in total.
- **Allowed after freeze 1:** barcode classes; compositions; hashtag calls, loads and pair matrices; synthetic flags; masked real-run thresholds and untyped called fractions; replica ψ̂ (synthetic only); simulations.
- **Freeze 2** follows stage 13 and records the eligible list, gates, rates, templates and verdict, plus `git rev-parse HEAD:scripts/dbl` and the md5s of the parameter table and of stage 14. Stage 14 refuses to run otherwise.
- **No third exception.** A machinery rerun is allowed only for a demonstrated code bug, with the same seeds.

## 7. Seeds (amendment: the redesign's SEED + 921/931/941/951 collide with stage 10's SEED + 901 + 10f)
| Stream | Seed |
|---|---|
| 12 --calibrate | SEED + 1201 (+ 100 × lane index + run) |
| 13 --audit | SEED + 1301 (+ 10000 × world + cohort) |
| 14 | SEED + 1401 |
| Stage 1 kNN transfer (02) | stage 02's DBL_SEED_TRANSFER, unchanged |

## 8. Exposures before this freeze
- v1.1 Appendix A and section 2.1; stage 11's templates; the 2 typed AML3948 singlets of 11-09; the 10-28 hashtag tally (GSM5613747-52).
- P01 count-1 hashtag counts and called fractions.
- All Stage 0 outputs (0a gates and reported items; 0b synthetic rates and power).
- **Class-A compositions of GSE185381, GSE289435 and Petti 2019 seen in other workstreams.** These are the ws atlas, S3-S9, smoke stages 02-11 on AML3948 / P01 count-1, and the redesign's composition-screened candidate list. The Stage 1 thresholds were fixed on 11-11 and are not changed here.
- Stage 01's class tables (A / B_mad / B_union / C counts per lane and donor) for every candidate lane.
