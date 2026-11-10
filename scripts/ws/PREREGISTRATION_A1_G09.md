# PREREGISTRATION A1-G09, stage 1: is S8's replicated structure specific to AML?

**Status: DRAFT v2 (2026-11-10), after review 1. Not approved, not tagged. Nothing declared in sections 3-6 has been
computed.** Decisions 9.1-9.9 were fixed by the user on 2026-11-10. Items 9.10-9.12 are new and remain open.

**Sources.** All are kept outside the repo, under `/FAST/gr10634/gaozy/reports/2026-11_a1_g09_prereg/`:
- `inventory/`: a read-only inventory (lenses A obligations and blindness, B machinery and data, C design and power, plus
  a critic).
- `review1/`: R1 design and statistics; R2 conformance, blindness and wording.

Definition sources, mapped in 0.1:
- the paper-gap assessment, `reports/2026-10_wpa1_paper_gap/WPA1_PAPER_GAP_ASSESSMENT_2026-10-28.md` (:69, :97-98,
  :160-161);
- A1-G01 5.4 (:282-286);
- the G01 review-4 v2 rule, `reports/2026-11_a1_g01_prereg_review/review4/A/v2.md` (:189-200; K6 :237-238);
- `reports/2026-11_T_shift_progress/writable_inventory_2026-11-10.md`.

## 0. Question and scope

S8 showed two structural properties of the P370 axis in adult AML cohorts that played no part in defining it:
- **T2:** across patients, the up-half and down-half means move oppositely.
- **T1:** within a patient's T cells, the axis genes co-vary more than matched random genes.

Result: PASS in GSE289435, consistent in GSE227903, not reproduced on Seq-Well GSE116256. Those cohorts have no healthy
donors. A gene set aligned with T-cell differentiation would pass T2 and T1 wherever T-subset composition varies, because
subtype heterogeneity alone produces both (G01 review 4; assessment :97-98).

**Stage 1 asks two things:**
- **(a)** Do T2 and T1 also pass among healthy donors (Oetjen, 20 donors)?
- **(b)** In the AML cohorts, do they survive random sets matched also on the naive-to-effector loading?

The fixed sentences for each outcome are in section 7 only (decision 9.5).

**S8's null splits the halves by bin order (found in review 1).**
- `matched()` draws one random gene per imprint gene from the union of both halves and returns them ordered by stratum.
  `split_random` (`32_s8_replication.R`:60) and the T1 draw (:107) then take `head()` as the random "up" half and
  `tail()` as the random "down" half.
- So the random halves are not exchangeable. For T1 the random "up" half is always the least-detected genes; for T2 the
  halves follow the alphabetical order of the pasted bin labels.
- S8's published nulls (`ws_s8_tests.csv`) are therefore not per-half matched. G09 adds a per-half null (N1) and reports
  what it does to S8's AML passes (3.2, 7.2).

**What G09 is not.**
- It is not a test of level or direction. A1-G02 holds the level difference, and the direction stays "shown in one
  study" (A1-G04).
- It does not compare strengths between AML and healthy donors (7.3).
- It settles neither cause nor age.

**Stage 2 is not run here (decision 9.1).** It is a separate amendment covering two parts:
- the power of S8 at n 7-8, from planted structure;
- the Seq-Well boundary, at the two-dimensional target of about 21 T cells and the measured GSE116256 median UMI per T
  cell (decision 9.4).

Until stage 2 runs:
- S8's Seq-Well sentence stays as written in GATES.md A1-S8.
- No power statement is made for S8.
- The AML side of the Seq-Well result is already public (T2 +0.60, T1 3 of 5), so stage 2 cannot be blind on that side.
- Its reading bands are left to the amendment.

### 0.1 Obligations, where they are handled, and declared departures

| source | obligation | here |
|---|---|---|
| G01 5.4; review-4 v2 :189-197 | T1/T2 among Oetjen donors; GSE185381 healthy donors reported, not used; the T1 75 % rule and the T2 rho < 0 / P5 rule | 3.1, 7.1 |
| G01 5.4; v2 :198-200 | random sets matched also on the 5.5 loading; GSE289435 and GSE227903 | 3.2, 7.2 |
| v2 :237-238 (former K6) | GSE227903 AML_Dx subset | K6 |
| assessment :69, :160-161 | power at n 7-8; Seq-Well downsampling | stage 2 (declared deferral) |
| inventory C | equalise cells and depth between AML and healthy donors | **not adopted**: stage 1 reads (a) as existence among healthy donors, not as a comparison of strengths (7.3) |

**Declared departures from the carried notes:**
1. The loading-matched null uses L deciles or quintiles, not terciles (review R1-3).
2. The loading for Oetjen T1 is leave-one-donor-out (9.3). For Oetjen T2 it comes from GSE185381 healthy donors (9.7,
   9.10).
3. A per-half null (N1) is added, because S8's null is not per-half matched (above).

## 1. Already known before this draft (disclosed)

1. **S8's AML statistics are public** (`ws_s8_tests.csv`; GATES.md A1-S8):

   | cohort | T2 rho (random P5) | T1 |
   |---|---|---|
   | GSE289435 | −0.55 (0.38) | 8 of 8 |
   | GSE227903 | −0.43 (0.21) | 7 of 7 |
   | discovery | −0.29 | 13 of 15 |
   | Seq-Well | +0.60 | 3 of 5 |

   Only the N1 and N2 nulls are new for the AML cohorts.
2. **No S8-style T1 or T2 has been computed among healthy donors.** Analogous quantities exist for the 10 GSE185381
   healthy donors: G07's `ws_a1_g07_s5b_patients.csv` holds per-donor per-cell `mean_z_up` and `mean_z_down` (depth-thinned,
   per-cell z, not S8's pseudobulk z), and `/LARGE1/.../ws_a1_g07/percell_scores_T1000_P370.rds` holds per-cell scores.
   A T2-like correlation among those donors is therefore derivable from a public table. Those donors are reported, not
   used.
3. **What is known about Oetjen.**
   - S9 published a per-donor combined score against age (W_rho_score_age 0.277, GATES.md:232; f_slope 0.102 above the
     random P95 0.046) and per-gene age logFC (`ws_s9_age_genes.csv`; A2b sign concordance 0.601).
   - The CD8 Naive share falls with age (rho −0.767, S9 K1).
   - Per-donor combined scores were read in S9's erratum (`ws_s9_direction_positions.csv`).
   - S9's table has no up-half or down-half means, so T2 cannot be derived from it; its score columns were not opened
     for this draft.
4. **P370's placement along the loading is partly known.**
   - G01 5.5 printed union correlations of 0.355 (CD4 Naive) and 0.313 (CD8 EM2) between the within-label AML shift
     and the loading (`ws_a1_g01_axis.csv`).
   - G07's Oetjen naive/EM tags (`tags_per_gene.csv.gz`; public `ws_a1_g07_tags.csv`, n_tag_O) and its K4 marker
     values were read. They are not readable, and G07 7.4's CMV sentence travels with any citation of them.
   - G07's 4.4-only sets were CD4_NAIVE_2 and CD8_NAIVE_2 down, CD4_CYTOTOXIC up.
5. **Only counts were read for this draft:** cell numbers, rosters, library maps and median UMI.
   - Oetjen: 20 donors, 104 to 3,721 kept T cells; C pools 3 libraries and S pools 4.
   - GSE185381 healthy: 10 donors, 123 to 1,917 kept T cells.
   - Griffin GSE227690: 5 donors over 13 libraries (1250, 820, 349, 2080, 533 kept T cells).
   - GSE276360: 4 donors.
   - GSE116256 healthy: 2 samples (113, 727 kept T cells).
   - GSE227903: 9 AML_Dx samples, 7 eligible (3853_Dg has 0 T cells, 6323_Dg has 12).

## 2. Data and frozen definitions

**Gene sets.**
- Primary: P370, S8's "primary" (`32_s8_replication.R`:38-39: 370 / 156 / 214).
- Sensitivity: F402, in tables only and never in running text (9.9).

**Cohorts.**

| cohort | samples | role |
|---|---|---|
| Oetjen GSE120221 | 20 donors; libraries pooled per donor as in S9 (`33_s9_age_direction.R`:110-126) | decides (a) |
| GSE185381 healthy | 10 donors, and without Control1-5 (5) | reported, not used: these donors were in the gene selection |
| GSE227690 (Griffin) | 5 donors, pooled per donor | T1 only, supporting (9.6) |
| GSE276360 | 4 donors | T1 only, supporting (9.6) |
| GSE116256 healthy | 2 Seq-Well samples | T1 only, descriptive |
| GSE289435 AML_Dx | 8 patients | (b) |
| GSE227903 AML_Dx | 7 eligible of 9 | (b) |

Why the small cohorts get T1 only: T1-only is a choice made on power, not a code limit. At n = 4 the smallest Spearman
p is 1/24. Under `min(6, n_ok)` (`32_s8_replication.R`:109), a supporting cohort must have every donor above its P95.

**Cells and genes.** These are S8's own (`32_s8_replication.R`:23, 47-61, 70-82): kept T cells (`group4 == "T"`), ≥ 20
per sample, the blacklist, CPM > 1 in at least half of the run's samples, and z-scores per gene across the run's samples.
For T1: a cap of 500 cells, genes detected in ≥ 10 % of those cells, ≥ 20 genes per half. Three changes are needed:
- the sample filter accepts `role == "control"`;
- a T1-only run bypasses S8's "≥ 4 samples" skip (:73);
- for a donor with several libraries, pseudobulks are summed and cells are pooled before the cap.

One cap subsample per donor (seed SEED + 901) is shared by every null.

**The loading L.** G01 5.5's code is copied verbatim (`39_a1_g01.R`:536-540):
- naive = pooled CD4 Naive + CD8 Naive; EM = pooled CD8 EM1 + CD8 EM2 + CD4 EM. Both are kept-cell fine pseudobulks.
- The two-column matrix goes through G01's `logz`, which sets the universe.
- `LOAD = log2((CPM_em + 1) / (CPM_naive + 1))` on TMM CPM.

Sources of L by run:
- AML runs: all 20 Oetjen donors, G01's exact vector.
- Oetjen T1 for donor d: the other 19 donors, dropping all of d's libraries.
- Oetjen T2: GSE185381 healthy donors without Control1-5 (9.7, 9.10); all 10 donors and the all-Oetjen L are
  sensitivities.

Oetjen is 25 % of BoneMarrowMap, so its labels are in-sample; this is disclosed. Genes without L are outside the matched
universe; coverage is reported.

**Cohort naive-share spread.** The per-donor (CD4 + CD8 Naive) / kept T share is reported for every cohort from cell
counts. T2 depends on this spread.

## 3. Analyses

Three nulls are used. Each draws the up half and the down half separately: `gu <- matched(U, up)`, then
`gd <- matched(U[!gene %in% c(up, dn, gu)], dn)`.

| null | strata, T2 | strata, T1 | role |
|---|---|---|---|
| N0 | S8's published null for the AML cohorts, taken from `ws_s8_tests.csv` as registered; for healthy cohorts, S8's code unchanged (head/tail split) | as T2 | context |
| N1 | per-half, S8's bins: expression deciles × 4 detection bins | per-half, detection deciles | primary for (a) (9.11) |
| N2 | per-half, S8's bins × L deciles | per-half, detection quintiles × L quintiles | (b) |

- **Pools.** Sampling is without replacement. If a stratum's pool is smaller than 3 × need, adjacent detection bins are
  merged within the same L stratum until it is not. Fallbacks, merges and the mean Jaccard overlap between draws are
  logged to the checks table.
- **Draws.** 1,000 per test for T2; 200 per donor for T1.

### 3.1 (a) Healthy donors
- Run T2 and T1 on Oetjen with N1 and S8's pass rules (`32_s8_replication.R`:92 for T2, :109 for T1).
- N0 and N2 are reported.
- The supporting cohorts get T1 only.
- Descriptive only: 500 draws of 8 of the 20 Oetjen donors (seed SEED + 902). Each draw is a full recompute (logz,
  universe, bins) with 200 N1 sets. The pass share is reported with no threshold and is never compared with the AML
  passes (7.3).

### 3.2 (b) AML cohorts
- On GSE289435 and GSE227903 AML_Dx, run T2 and T1 with N1 and N2.
- The published N0 result is quoted from `ws_s8_tests.csv`.
- The re-run of S8's own null is a check only (K0b).
- Reading: 7.2.

### 3.3 Reported for every run
- n and P370 coverage (S8 requires ≥ 0.8);
- the median T1 coherence and the median random P95 per null;
- per-donor values in a work table outside the repo;
- the naive-share spread.
- Oetjen only: the up-half and down-half means against donor age, separately (new; S9 published only the combined score
  against age).

## 4. Sensitivities (reported; none changes a reading)
1. F402 in place of P370 (tables only).
2. Oetjen T2 with L from all 10 GSE185381 healthy donors, and with L from all 20 Oetjen donors (circular; stated).
3. GSE185381 healthy donors: all 10, and without Control1-5.
4. T1 on Oetjen donors C and S, and on Griffin donors, using the largest library only (library effects in pooled cells).

## 5. Multiplicity
Two tests (T2, T1) per cohort and null; each is read as pass or fail against its own null, as in S8. No p values are
combined. The supporting cohorts never change a reading.

## 6. Known answers (run first; no readout is printed until each has a state)

States are as in G07: PASS / FAIL deterministic (a code bug, fixed with a dated note) / FAIL statistical (the named reading
becomes "not readable").

| id | check | bar | type |
|---|---|---|---|
| K0 | S8's code path reproduces the observed T2 rho on GSE289435 and GSE227903, and the observed T1 coherence (P370) of every patient with ≤ 500 T cells | exact (< 1e-12) | deterministic |
| K0b | Re-run of S8's N0 null on the AML cohorts: pass flags equal those in `ws_s8_tests.csv` | equal | reported; if unequal, fixed sentence 7.2(iv) |
| K1 | Rosters. Oetjen: 20 donors, S9's library map, per-donor kept T cells equal to `n_kept_T`, read with `fread(select = c("level", "id", "n_libraries", "n_kept_T"))` so that no score column is loaded. Griffin: 5 donors with 1250 / 820 / 349 / 2080 / 533 kept T cells | exact | deterministic |
| K2 | L is G01's: recomputed all-Oetjen L, with G01's GSE185381 within-label shifts, reproduces `ws_a1_g01_axis.csv` `union_rho_number_only` and `n_union` | < 1e-12; integers equal | deterministic |
| K2b | Leave-one-donor-out with no donor left out equals all-Oetjen L; each leave-one-out run drops all of that donor's libraries (C 3, S 4) | exact | deterministic |
| K3 | Null world for T2: Oetjen raw counts permuted across donors independently per gene, then logz re-run; 100 permutations; N1 | T2 passes in ≤ 10 of 100 (false alarm 0.011) | statistical ((a) T2 "not readable") |
| K4 | Machinery check, not power: on a K3 permutation, a donor factor f ~ N(0, 1) planted on raw counts before logz, up × 2^(0.3 f) and down × 2^(−0.3 f) on P370; 20 plants. Plants at 0.05 and 0.1 are reported as numbers | T2 passes in ≥ 16 of 20 | statistical ((a) T2 "not readable") |
| K5 | Composition world for T2. Start from a K3 permutation (real marginals and noise kept); multiply each gene by (p_d·EM + (1 − p_d)·naive) / (p̄·EM + (1 − p̄)·naive), with p_d drawn uniformly over the observed Oetjen naive-share range; n = 20. Test set: genes drawn to match P370's L-decile profile (up = EM side, down = naive side); 20 worlds. The 100-extreme set is a reported stress case | N1: T2 passes in ≥ 16 of 20 (the confound is real). N2: T2 passes in ≤ 4 of 20 (matching works; false alarm 0.003 at 5 %) | statistical ((b) "not readable") |
| K5b | The same composition world for T1 within one Oetjen donor (K, 500 cells): cells reweighted to mix naive-like and EM-like profiles, with L leave-one-donor-out; 20 worlds | N1 passes in ≥ 16 of 20; N2 in ≤ 4 of 20 | statistical ((b) T1 "not readable") |
| K6 | GSE227903 (verbatim v2 :237-238): "the sample list is exactly the 9 AML_Dx samples of the manifest; non-diagnosis columns are dropped right after loading, before any use; per-cell objects are loaded for AML_Dx samples only". MRD and relapse samples are not touched (WP-A2). z-scores and the universe come from the 7 eligible samples | exact | deterministic |
| K7 | Null world for T1: log-normalised rows permuted across cells per gene, within Oetjen donor K (500 cells; library sizes kept); 20 permutations; N1 | coherence above P95 in ≤ 3 of 20 (false alarm 0.016) | statistical ((a) T1 "not readable") |
| K8 | Balance: the per-stratum counts of every random half equal those of its imprint half; for N2, the mean over draws of \|mean L(random half) − mean L(imprint half)\| / SD(L) < 0.1 per half | exact; < 0.1 | deterministic (counts); statistical, (b) "not readable" (balance) |

Run-time blindness: no healthy T2 or T1 value and no N1 or N2 AML value is printed until K0-K8 have states.

## 7. Reading rules and fixed sentences (the only verbatim texts)

**7.1 (a): Oetjen, P370, N1.** Name only the test that passed.
- **(i) T2 or T1 passes:** "Among the T cells of one healthy marrow cohort (Oetjen 2018, 20 donors, 10x 3' v2, one
  protocol; T labels in-sample to BoneMarrowMap), the axis genes show the structure that S8 tested ([T2: the up-half and
  down-half means move oppositely across donors, rho r, random P5 p] / [T1: they co-vary within donors, k of n donors
  above their random P95]). S8 therefore shows that the axis's structure is reproducible on 10x Chromium cohorts, not
  that it is specific to AML."
- **(ii) Neither passes:** "Among 20 healthy Oetjen donors (one cohort, one protocol), this structure was not detected
  (T2 rho r, random P5 p; T1 k of n donors). This is not evidence that the structure is specific to AML: depth, cell
  numbers, protocol, age and the spread of T-subset composition differ from the AML cohorts, and the power to detect T1
  among healthy donors was not established."
- Every (a) sentence adds: "Structure is covariance, not level: this does not change the A1-G02 numbers. It is not
  evidence for or against composition below the label level (A1-G01 11.5). The direction stays shown in one study
  (A1-G04)."

**7.2 (b): GSE289435 and GSE227903, P370.** One sentence per test and cohort; the two cohorts are named together only if
they agree.
- **(i) Passes N1 and N2:** "S8's [T2 / T1] in [cohort] exceeds random gene sets matched separately for each half on
  expression, detection and the naive-to-effector loading."
- **(ii) Passes N1, fails N2:** "S8's [T2 / T1] in [cohort] is consistent with the naive-to-effector differentiation
  axis."
- **(iii) Fails N1:** "Against random sets matched separately for each half, S8's [T2 / T1] in [cohort] does not pass;
  S8's published pass used a null that split the random halves by bin order." This case is also recorded as a dated
  erratum note after GATES.md block A1-S8 (decision 9.12).
- **(iv) K0b unequal:** "Re-running S8's own null changed its pass flag for [test, cohort] (published p, re-run p);
  Monte-Carlo error of the null." This sentence is reported only.
- Every (b) sentence adds: "This concerns the covariance of S8's structure, not the within-label shift of A1-G01 (5.5 did
  not fire)."
- **Joint case, (a)(ii) with (b)(ii):** add "The AML structure is consistent with the naive-to-effector axis, so its
  absence among healthy donors may reflect a narrower spread of T-subset composition."

**7.3 Never written, under any outcome:**
- "AML-specific" for S8's structure;
- "stronger", "weaker" or "more coherent in AML";
- any comparison of the 8-of-20 draw pass share with the AML passes;
- the draws as S8's power;
- "not explained by differentiation" (the ceiling is 7.2(i));
- "coherent T-cell program", or "program" for the axis;
- "S8 replicates the direction";
- "replicated on Seq-Well";
- any causal or age-excluding wording.

Every never-written item of S8, S9, A1-BATCH 12.5, A1-G01 11.5, A1-G04 11 and A1-G07 7.8 still applies.

**7.4 Context that travels with any G09 sentence:**
- the 32 ambient-prone genes "stay excluded from every statement";
- F402 appears in tables only, and no P370/F402 difference is attributed to the 32;
- S8 is "shown on 10x Chromium cohorts; not reproduced on Seq-Well GSE116256", and its T3 is uninformative;
- Oetjen is one cohort and one protocol (10x 3' v2), with labels in-sample to BoneMarrowMap;
- GSE185381's Control1-5 are T-ALL lanes, and Control4 is female by expression;
- with any (b) sentence: the A1-G01 category and its 11.5 sentence;
- with any healthy-versus-AML contrast: G07 7.7's acute-illness sentence, and that CMV is unrecorded (CMV expands
  TEMRA-like cells inside the CD8 EM labels and rises with age);
- the S9 age bound.

## 8. Code, outputs, compute

- **Script:** `scripts/ws/45_a1_g09.R`, new. It copies S8's and G01's functions with their sources; `32` and `39` stay
  untouched.
  - Modes: `--known`, `--run`, `--combine`, each also with `--smoke`.
  - G07's guards apply: unknown arguments stop; real runs use only the dated read-only copy of a committed script;
    nothing is overwritten.
- **Smoke runs:** synthetic counts on the real IDs, including the K3, K5, K5b and K7 null worlds. Smoke outputs carry
  `_smoke` and are written only under `/LARGE1/.../ws_a1_g09/smoke/`.
- **Outputs:** `results/tables/ws/ws_a1_g09_{checks,tests,summary}.csv`; per-donor work tables under
  `/LARGE1/.../ws_a1_g09/`; a GATES.md block "A1-G09".
- **Seeds:** one stream per block, with `set.seed` called per (cohort, gene set, null) block so that a subset of cohorts
  reproduces.
  - SEED + 901: cell caps.
  - SEED + 902: 8-of-20 draws.
  - SEED + 903: N0 re-run.
  - SEED + 904: N1.
  - SEED + 905: N2.
  - SEED + 906: K3-K5b worlds.
  - SEED + 907: K7.
- **Compute:** one `c=16, m=64G` job, mclapply with `OMP_PROC_BIND=false`. Reading per-sample objects for
  25 + 13 + 4 + 10 + 2 healthy and 8 + 7 AML samples, plus the 500-draw subsampling, should take under 3 h. No new
  package.

## 9. Decisions

Fixed by the user on 2026-11-10 ("按推荐"):
- **9.1** Stage 1 = (a) + (b); power and Seq-Well become stage 2.
- **9.2** Oetjen n = 20 with S8's pass rules is primary; the 8-of-20 draws are descriptive.
- **9.3** The loading for Oetjen T1 is leave-one-donor-out; the matched null applies to T1 and T2.
- **9.4** The Seq-Well depth comes from measured counts (stage 2).
- **9.5** Wording fixed now (section 7).
- **9.6** GSE185381 healthy donors reported only, with and without Control1-5; Griffin and GSE276360 are T1 support; the
  GSE116256 healthy samples are descriptive.

Fixed on 2026-11-10 ("可以，先按推荐继续"):
- **9.7** The Oetjen T2 loading comes from GSE185381 healthy donors; the all-Oetjen loading is a sensitivity.
- **9.8** The known-answer bars as in section 6.
- **9.9** F402 is a sensitivity only.

New after review 1 (open; recommendation first):
- **9.10** The GSE185381 donors for the 9.7 loading. Use the 5 donors without Control1-5, because T-ALL ambient RNA in
  P02/P03 can distort T-lineage profiles; all 10 is a sensitivity. Alternative: all 10.
- **9.11** The primary null for (a) is N1 (per-half), because S8's own null splits halves by bin order. Report N0 and
  N2. Alternative: keep S8's null unchanged as primary, which reproduces S8 exactly but carries the split.
- **9.12** If an AML pass fails under N1 (7.2(iii)), record a dated erratum note after GATES.md A1-S8, and the S8
  claim narrows accordingly. Alternative: report it in G09 only.

## 10. Review history
- 2026-11-10: inventory (four agents, read-only). DRAFT v1 (3e070f8).
- 2026-11-10: decisions 9.7-9.9 as recommended.
- 2026-11-10: review 1, two read-only agents (R1 design and statistics, 15 findings, 4 blocking; R2 conformance,
  blindness and wording, 15 findings, 3 blocking).
  - Blocking items fixed in v2:
    - per-half matching; S8's head/tail split disclosed and an N1 null added;
    - finer L strata with a balance check (K8);
    - K5 rebuilt from a permuted real world, and K5b added for T1;
    - the sentence "coherent T-cell program" removed;
    - "AML-specific" never written under any outcome;
    - "not explained by differentiation" never written.
  - Should-level items fixed in v2:
    - K2 is exact against G01's table;
    - K0b and the published-null rule;
    - K6 is verbatim;
    - section 1 is completed (G07 S5b per-donor means, S9 score-versus-age, G07 tags, S9 erratum);
    - the 7.4 context is completed;
    - the stage-2 disclosure is completed;
    - an obligations table (0.1);
    - seeds are separated;
    - T1-only is justified by power;
    - the T1 null is described as detection deciles;
    - the empty sensitivity is dropped;
    - K1 reads count columns only.
