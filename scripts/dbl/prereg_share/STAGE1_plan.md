# Stage 1: freeze 1 and eligibility (plan v1, 2026-11-12)

Status: PLAN for the user's approval. Nothing is coded or run yet.

Stage 0 is complete:
- 0a: "K1 usable" (SMOKE_REPORT section 6h);
- 0b: "pass", with planning power 0.985 at ψ 2 and 0.669 at ψ 1.5 for n = 12 (section 6i).

By the redesign (REDESIGN_share_plan.md, sections 7 and 8), the next unit is **freeze 1**. Everything the share test needs is fixed before any new lane's composition is tabulated. The **Stage 1 gate** follows: at least 10 patients with strict T ≥ 30 and LEUK ≥ 50 class-A cells, summed over each patient's lanes; otherwise BN-1.

## 1. What changed since the redesign (found 2026-11-12)
- **The candidate list.** The redesign named 22 candidates: GSE185381 11, GSE289435 6, Petti 5. That list was evidently already screened on compositions seen in earlier work; for example AML3133 (LEUK 48) and AML0160 (T 3) are absent.
  - Freeze 1 therefore defines candidates **by metadata only**, and leaves all composition screening to the pre-registered Stage 1 gate.
  - By metadata: GSE185381 has 42 AML diagnosis samples (about 35 after the exclusions below); GSE289435 has 10 (AEL MLL_28824 and MDS MLL_28855 excluded); Petti has 5.
- **Exposure, declared in freeze 1.** Class-A compositions of all three datasets have been seen in other work: the ws atlas, S3-S9, and smoke stages 08-11. The Stage 1 gate (T ≥ 30, LEUK ≥ 50, n ≥ 10) was fixed in the redesign on 11-11 and is not changed.
- **Chemistry.** GSE289435 and Petti are 10x 5′ (curated manifest: 10x_5p_v1.1 and 10x_5p_v1). They have one lane per patient and no hashtags. M1 and the T-identity gate were developed only on 3′ data; Stage 2's per-lane typing gates decide whether 5′ lanes stay in.

## 2. Steps (each step is its own commit; the step after needs the one before)
**S1a. Lane table.**
- New scripts/dbl/00_share_lanes.R (about 80 lines) writes scripts/dbl/dbl_lanes_share.tsv, which is committed.
- Inputs: the curated manifest, qc_percell (Library, Demuxed) and the GEO feature files (are there HTO rows?). It reads no cell labels and no composition.
- Rules:
  - GSE185381: AML diagnosis, all of a patient's 3′ libraries. Excluded:
    - the T-ALL co-capture lanes GSM5613747-55 (memory: healthy pools held T-ALL);
    - genotype-only pools, except author-assigned barcodes (redesign decision 11).
  - GSE289435: AML diagnosis, excluding the AEL and MDS samples.
  - Petti 2019: all 5 samples.
- Columns: dataset, gsm, lane, patient, age group (adult / child), chemistry, donors in lane, hashed, role (analysis / K1 / excluded), reason.
- Check: the lanes of a patient that spans libraries do not share barcodes (the second-lane overlap check).

**S1b. Generalise 01 to the lane table and extract every analysis lane.**
- 01 gains about 40 lines: it reads the lane table, handles the GSE289435 and Petti ingest layouts (one sample per library, no HTO), and handles hashed pools with unnamed tags (author majority ≥ 0.80).
- Output: lane objects only.
- **02 does not run yet.** It tabulates compositions, which are allowed only after freeze 1.

**S1c. The freeze-1 package**, committed and tagged prereg-dbl-share-f1:
- STAGE1_freeze1.md holds every rule:
  - the estimand, the flags Y and D, the lane gates, the model and test, the binding variants, the worlds and bars, the reading rule, and the BN texts;
  - the amendments from Stage 0: K1 uses doublet-like cross calls only, or models the singlet share (0a); calibration variance at the design sizes (0b);
  - the full exposure list.
- config_dbl.R section 11: the remaining share constants. These are the allocation, gates, σ, τ floor, singlet grid, worlds, bars and seeds.
- Rule-text headers for 12 --calibrate, 13 --audit and 14_share_unblind.R. The code for these comes later, as its own approved steps.

**S1d. The Stage 1 run (after the tag).**
- 02 runs on every analysis lane (class-A labels and strict T; kNN labels for B/C).
- The eligibility table gives strict T and LEUK in class A per patient, summed over lanes.
- Verdict: n ≥ 10 → Stage 2 (calibration, about 40-60 CPU-h). Fewer → BN-1 "lineage-level T-LEUK pairing has no validated null in these data (n < 10)".

Effort: S1a + S1b about 1 analyst day; S1c about 1 day of rule text; S1d under 1 h of compute.

## 3. The one decision needed now: the content of freeze 1
Adopt the redesign's synthesis defaults (decisions 3-12) as freeze-1 rules, together with the Stage 0 amendments:
- the share estimand ψ replaces the absolute O/E (3);
- one common ψ for children and adults; adults-only and children-only are reported (4);
- n ≥ 10, and no third exception: a machinery rerun only for a demonstrated code bug (5);
- lane typing gate: M1 precision ≥ 0.70, recall ≥ 0.50 (6);
- shared-bias prior: common sd 0.08, plus per dataset 0.08 for GSE185381 and 0.12 for GSE289435 and Petti (7);
- strict singlet masking: singlet rates calibrated only inside the unblinding worker (8);
- binding variants: the LEUK-share window and the conservative composition (9);
- P = {T, Tamb, NK, B, Ery, MonoDC, Other} (10);
- AML009 included via author-assigned barcodes (11);
- if freeze 2 says "only a doubling is detectable", unblinding goes ahead (12).

The expected end state stays as 0b put it: only a doubling is detectable at n ≈ 12.
