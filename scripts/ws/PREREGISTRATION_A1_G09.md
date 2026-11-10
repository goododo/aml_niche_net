# PREREGISTRATION A1-G09, stage 1: is S8's replicated structure specific to AML?

**Status: v1, APPROVED by the user 2026-11-10 ("可以，然后根据修改后的正确步骤再进行test": all section-9 decisions fixed,
proceed with the corrected steps), after reviews 1 and 2 (section 10). Tag `prereg-a1-g09-v1`. Nothing declared in
sections 3-6 was computed before the tag.** Decisions 9.1-9.9 were fixed by the user on 2026-11-10. Items 9.10-9.12 were fixed the same day.

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

Result: PASS in GSE289435, consistent in GSE227903, not reproduced on Seq-Well GSE116256. S8 tested AML patients only
(GSE116256's two healthy samples were not scored). A gene set aligned with T-cell differentiation can pass T2 and T1
wherever T-subset composition varies, because subtype heterogeneity alone can produce both (G01 review 4; assessment
:97-98).

**Stage 1 asks two things:**
- **(a)** Do T2 and T1 also pass among healthy donors (Oetjen, 20 donors)?
- **(b)** In the AML cohorts, do they survive random sets matched also on the naive-to-effector loading?

The fixed sentences for each outcome are in section 7 only (decision 9.5).

**S8's null splits the halves by bin order (found in review 1).**
- `matched()` draws one random gene per imprint gene from the union of both halves and returns them ordered by stratum.
  `split_random` (`32_s8_replication.R`:60) and the T1 draw (:107) then take `head()` as the random "up" half and
  `tail()` as the random "down" half.
- So the random halves are not exchangeable.
  - T1: the strata are a `cut()` factor (:105), so the random "up" half always takes the draws for the lowest detection
    deciles of the combined imprint profile and the "down" half the highest (synthetic check, review 2: random-half
    detection 0.17 vs 0.47, imprint halves 0.35 vs 0.33).
  - T2: the strata are pasted character labels (:52), sorted by locale collation. Under S8's R build and locale the
    random "up" half takes expression deciles 2 upward and the "down" half the top deciles plus decile 1 (synthetic: mean
    decile 3.8 vs 7.1). The order depends on the R build and locale, which K0b pins.
- S8's published T2 P5 values were positive (+0.38, +0.21), so its T2 passes rested on rho < 0. Whether the split changed
  any pass is unknown until N1 runs.
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

**Cohort naive-share spread.** The per-donor share (CD4 + CD8 Naive) / kept T cells, from cell counts, is reported for
every cohort as its SD and range. T2 depends on this spread.

## 3. Analyses

Four nulls are used. N1, N1-L and N2 draw the up half and the down half separately with `matched_x(U, target, excl)`:
strata come from `target`; the pool of stratum s is `U[stratum == s & !gene %in% excl]`. Then
`gu <- matched_x(U, up, c(up, dn))` and `gd <- matched_x(U, dn, c(up, dn, gu))`. There is no whole-universe fallback
and no sampling with replacement.

| null | strata, T2 | strata, T1 | role |
|---|---|---|---|
| N0 | S8's published null for the AML cohorts, taken from `ws_s8_tests.csv` as registered; for healthy cohorts, S8's code unchanged (head/tail split) | as T2 | context |
| N1 | per-half, S8's bins on S8's universe: expression deciles × 4 detection bins | per-half, detection deciles on S8's per-cell universe (no blacklist) | primary for (a) (9.11); 7.2(iii) |
| N1-L | as N1, on N2's universe | per-half, detection quintiles, on N2's universe | comparator for 7.2(i)/(ii) |
| N2 | per-half, S8's bins × L deciles, on N2's universe | per-half, detection quintiles × L quintiles, on N2's universe | (b) |

N2's universe = the run's (T2) or donor's (T1) universe minus S8's `black_for()` genes, intersected with the genes that
have L. N1-L differs from N2 in L only, so 7.2(i)/(ii) contrasts N1-L with N2; N1 on S8's own universe isolates the
per-half change from N0 for 7.2(iii).

- **Strata.** L strata are `cut(rank(L, ties.method = "first"), k, labels = FALSE)` (G01's idiom, `39_a1_g01.R`:544),
  computed over the run's or donor's matched universe including the imprint genes; stratum IDs are integers.
- **Pools and merges.** need = n_up + n_dn of the stratum; pool = its genes excluding all imprint genes. If pool < 3 ×
  need, the stratum is merged, in this order: (1) with the adjacent detection bin with the smaller pool, within the same
  expression decile and L stratum; (2) for T2, with the adjacent expression decile; (3) with the adjacent L stratum
  toward the median. Merges are fixed once per run (T2) or donor (T1) and shared by both halves. If a pool is still below
  need, that is a deterministic FAIL (never replacement). Merges and the mean Jaccard overlap between draws are logged.
- **Draws.** 1,000 per test for T2; 200 per donor for T1. Per-stratum pools are precomputed once per run with `split()`.

### 3.1 (a) Healthy donors
- Run T2 and T1 on Oetjen with N1 and S8's pass rules (`32_s8_replication.R`:92 for T2, :109 for T1).
- N0 and N2 are reported.
- The supporting cohorts get T1 only.
- Descriptive only: 500 draws of 8 of the 20 Oetjen donors (seed SEED + 902). For T2 each draw is a full recompute
  (logz, universe, bins, merges) with 200 per-half N1 sets. For T1, the per-donor N1 flags from the n = 20 run are reused
  over the same draws (rule 6 of 8), because T1 does not depend on the other donors. The pass shares are reported with
  no threshold and are never compared with the AML passes (7.3).

### 3.2 (b) AML cohorts
- On GSE289435 and GSE227903 AML_Dx, run T2 and T1 with N1, N1-L and N2.
- The published N0 result is quoted from `ws_s8_tests.csv`.
- The re-run of S8's own null is a check only (K0b).
- One GSE289435 patient has 543 T cells; its N1 result uses the G09 cap subsample, so its N0-versus-N1 comparison mixes a
  cap-subsample change with the null change. This is disclosed beside its row.
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
| K0 | S8's code path reproduces the observed T2 rho on GSE289435 and GSE227903, and the observed T1 coherence (P370) of every patient with ≤ 500 T cells (7 of 8 GSE289435 patients; 7 of 7 GSE227903) | exact (< 1e-12 against the stored 15-digit values) | deterministic |
| K0b | Re-run of S8's N0 null on the AML cohorts, with S8's R build (general_env R 4.4.3) and locale (en_US.UTF-8): pass flags equal those in `ws_s8_tests.csv` | equal | reported; if unequal, fixed sentence 7.2(iv) |
| K1 | Rosters. Oetjen: 20 donors, S9's library map, per-donor kept T cells equal to `n_kept_T`, read with `fread(select = c("level", "id", "n_libraries", "n_kept_T"))` so that no score column is loaded. Griffin: 5 donors with 1250 / 820 / 349 / 2080 / 533 kept T cells | exact | deterministic |
| K2 | L is G01's. Inputs: `ws_pb_GSE120221.rds` `$fine` (5 labels); `ws_pb_GSE185381.rds` `$fine[["CD4 Naive"]]`, `$fine[["CD8 Effector Memory 2"]]`, `$cells` (kept, for cell counts) and `$genes` (G01's BLACK, `39_a1_g01.R`:75-76); `ws_s5_scores.csv` T rows (25 samples and status); `ws_a1_g02_fold_genes.csv` (the union); WS_S4_MIN_CELLS = 20; G01's `logz` (:91-94) and L block (:536-552), copied verbatim. The recomputation reproduces `ws_a1_g01_axis.csv` `statistic` and `union_rho_number_only`, and `n_genes` (9724, 10939) and `n_union` | < 1e-12; integers equal | deterministic |
| K2b | Leave-one-donor-out with no donor left out equals all-Oetjen L; each leave-one-out run drops all of that donor's libraries (C 3, S 4) | exact | deterministic |
| K3 | Null world for T2: Oetjen CPM permuted across donors independently per gene and rescaled to each donor's own library size, then logz re-run; 100 permutations; N1 | T2 passes in ≤ 10 of 100 (false alarm 0.011) | statistical ((a) T2 "not readable") |
| K4 | Machinery check, not power: on a K3 permutation, a donor factor f ~ N(0, 1) planted on raw counts before logz, up × 2^(0.3 f) and down × 2^(−0.3 f) on P370; 20 plants. Plants at 0.05 and 0.1 are reported as numbers | T2 passes in ≥ 16 of 20 | statistical ((a) T2 "not readable") |
| K5 | Composition world for T2. Start from a K3 permutation; per donor draw s_d ~ U(observed min, max Oetjen naive share); multiply each gene by r = (s_d·naive + (1 − s_d)·EM) / (s̄·naive + (1 − s̄)·EM), using naive and EM CPM profiles from the GSE185381 healthy donors without Control1-5 (multiplier 1 for genes without a profile; non-integer counts kept); n = 20. Test set: a fresh set per world (SEED + 906), drawn from the universe minus P370 with P370's per-L-decile counts (up = EM side, down = naive side). N2 uses the all-Oetjen L; 1,000 null draws per world; 20 worlds. Reported, not gating: the 100-extreme set (top 50 and bottom 50 genes by L) and the same-source best case (profiles from Oetjen) | N1-L: T2 passes in ≥ 16 of 20 (the confound is real). N2: T2 passes in ≤ 4 of 20 (matching works; false alarm 0.003 at 5 %) | statistical (only the 7.2(i)/(ii) distinction for T2 becomes "not readable") |
| K5b | Composition world for T1 in Oetjen donor K (500 cells): its raw counts permuted per gene across cells; per-cell multipliers r_gc as in K5 with weight w_c = 1 for EM-labelled cells, 0 for naive-labelled cells and the donor mean for the rest, applied by binomial thinning (r / max r); profiles from donor K's own naive and EM cells; N2 with the leave-one-donor-out L; S8's log-normalisation; 200 null draws per world; 20 worlds | N1-L passes in ≥ 16 of 20; N2 in ≤ 4 of 20 | statistical (only the 7.2(i)/(ii) distinction for T1 becomes "not readable") |
| K6 | GSE227903 (verbatim v2 :237-238): "the sample list is exactly the 9 AML_Dx samples of the manifest; non-diagnosis columns are dropped right after loading, before any use; per-cell objects are loaded for AML_Dx samples only". MRD and relapse samples are not touched (WP-A2). z-scores and the universe come from the 7 eligible samples | exact | deterministic |
| K7 | Null world for T1: log-normalised rows permuted across cells per gene, within Oetjen donor K (500 cells; library sizes kept); 20 permutations; N1 | coherence above P95 in ≤ 3 of 20 (false alarm 0.016) | statistical ((a) T1 "not readable") |
| K8 | Balance, on merged strata: the per-stratum counts of every random half equal those of its imprint half, and no draw is NA; for N2, the mean over draws of \|mean L(random half) − mean L(imprint half)\| / SD(L) < 0.1 per half, and for T1 N2 the same for detection. T2: per run; T1: per donor | counts and NA: exact; balance < 0.1 | deterministic (counts, NA); statistical for balance: a failing run, or any failing donor for T1, makes that cohort's 7.2(i)/(ii) distinction for that test "not readable" |

Run-time blindness: no healthy T2 or T1 value and no N1 or N2 AML value is printed until K0-K8 have states.

## 7. Reading rules and fixed sentences (the only verbatim texts)

**7.1 (a): Oetjen, P370, N1.** Name only the test that passed.
- **(i) T2 or T1 passes:** "Among the T cells of one healthy marrow cohort (Oetjen 2018, 20 donors, 10x 3' v2, one
  protocol; T labels in-sample to BoneMarrowMap), the axis genes show the structure that S8 tested ([T2: the up-half and
  down-half means move oppositely across donors, rho r, random P5 p] / [T1: they co-vary within donors, k of n donors
  above their random P95]). S8 therefore shows that the axis's structure is reproducible on 10x Chromium cohorts, not
  that it is specific to AML." If 7.2(iii) holds for the same test in GSE289435, the last sentence is replaced by: "This
  structure is therefore not specific to AML; S8's AML pass for [test] did not survive the per-half null (7.2(iii))."
- **(ii) Neither passes:** "Among 20 healthy Oetjen donors (one cohort, one protocol), this structure was not detected
  (T2 rho r, random P5 p; T1 k of n donors). This is not evidence that the structure is specific to AML: depth, cell
  numbers, protocol and age differ from the AML cohorts[, and the spread of T-subset composition is narrower among them],
  and the power to detect T2 or T1 among healthy donors was not established." The bracketed clause is written only if
  the Oetjen naive-share SD is below that of both AML cohorts.
- Every (a) sentence adds: "Structure is covariance, not level: this does not change the A1-G02 numbers. It is not
  evidence for or against composition below the label level (A1-G01 11.5). The direction stays shown in one study
  (A1-G04)."

**7.2 (b): GSE289435 and GSE227903, P370.** One sentence per test and cohort; the two cohorts are named together only if
they agree.
- **(i) Passes N1-L and N2:** T2: "S8's T2 in [cohort] exceeds random gene sets matched separately for each half on
  expression, detection and the naive-to-effector loading." T1: "S8's T1 in [cohort] exceeds random gene sets matched
  separately for each half on detection and the naive-to-effector loading."
- **(ii) Passes N1-L, fails N2:** "S8's [T2 / T1] in [cohort] is consistent with the naive-to-effector differentiation
  axis."
- **Not readable (K5, K5b or K8 balance failed):** "Whether S8's [T2 / T1] in [cohort] goes beyond the naive-to-effector
  axis is not readable: the loading-matched null failed its known answer ([check])." 7.2(iii) is unaffected.
- **(iii) Fails N1:** "Against random sets matched separately for each half, S8's [T2 / T1] in [cohort] does not pass;
  S8's published pass used a null that split the random halves by bin order." This case is also recorded as a dated
  erratum note in GATES.md (decision 9.12), written by `--combine` after the "Reading" paragraph of block A1-S8 and
  outside its `<!-- A1-S8 -->` markers, verbatim: "Erratum (A1-G09, [date]): S8's random-set null drew one gene per
  imprint gene and split the list by bin order, so the random halves were not matched separately. Against a per-half
  null (A1-G09 N1): [test, cohort, statistic, P5 or k of n] [pass / fail]. GATE_S8 under the per-half null: [PASS /
  FAIL] (T2 and T1 on GSE289435)."
- **(iv) K0b unequal:** "Re-running S8's own null changed its pass flag for [test, cohort] (published p, re-run p);
  Monte-Carlo error of the null." This sentence is reported only.
- Every (b) sentence adds: "This concerns the covariance of S8's structure, not the within-label shift of A1-G01 (5.5 did
  not fire)."
- **Joint case, (a)(ii) with (b)(ii):** add "The AML structure is consistent with the naive-to-effector axis, so its
  absence among healthy donors may reflect a narrower spread of T-subset composition." Only if the Oetjen naive-share SD
  is below that of both AML cohorts.

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

Fixed on 2026-11-10 ("可以"), as recommended:
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
- 2026-11-10: decisions 9.10-9.12 as recommended. Review 2 (one read-only agent: fix audit and implementability; 12
  findings, 1 blocking). Fixed in v3: the per-half draw (`matched_x`, no fallback, no replacement, NA draws fail); the
  head/tail description made concrete with its locale dependence; N1-L added so that N1-L vs N2 differs in L only; strata
  and merge rules fixed; K3, K5 and K5b defaults; narrower consequences of K5/K8 failures and a not-readable sentence; the
  joint case and the S8 erratum verbatim; test-specific 7.2(i); the spread clause made conditional; K2 inputs and K0
  coverage; the 8-of-20 draws for T2 only.

## 11. Dated note, 2026-11-10, before any real-data run (code and smoke runs)

Written after the tag and before any real-data mode ran. No G09 statistic on real data had been computed. Sections 0-10
are unchanged. Script: `scripts/ws/45_a1_g09.R`, run through `scripts/ws/submit_ws_a1_g09.sbatch`.

**11.1 Changes found necessary in smoke runs** (synthetic counts on the real IDs; no real expression read):
1. **Finer L strata.** With L deciles (T2) and L quintiles (T1), the K8 balance measure was 0.07-0.10 for T2 and
   0.18-0.21 for T1 in every smoke run, failing or nearly failing the 0.1 bar for reasons of design: the imprint halves
   sit at the extremes of L. The strata are now:
   - T2 N2: S8's bins × L ventiles (20);
   - T1 N1-L and N2: detection deciles; N2 adds L deciles.

   The merge rule, the bar (0.1) and everything else are unchanged. The smoke balance became 0.03-0.06 (T2) and
   0.09-0.10 for L and 0.04-0.06 for detection (T1).
2. **Seeds.** N1-L draws use SEED + 908 and the smoke worlds SEED + 909. Every block seed is SEED + offset × 100000 plus
   a deterministic hash of the block label, so seeds do not depend on the order of forked workers.
3. **Outputs.** `--run` writes K8 to `ws_a1_g09_balance.csv` before any test table, plus `ws_a1_g09_draws8.csv`. Work
   files under `/LARGE1/.../ws_a1_g09/`: `donor_values.csv`, `naive_share_spread.csv` and `k0b.csv`.
4. **K5b test set.** Drawn like K5's: P370's per-L-quantile counts, from genes outside P370, fresh per world.

**11.2 Two fixed sentences added for cases section 7 did not cover:**
- (a) not readable: "Whether the healthy Oetjen donors show S8's [T2 / T1] is not readable: its known answer failed
  ([K3/K4 / K7])."
- (b) N1 passes but N1-L fails: "S8's [test] in [cohort] does not pass on the loading-matched universe even before
  matching on the loading (N1-L); no (b) reading is made."

**11.3 Smoke-run bugs fixed** (no design change):
- the loading matrix lost its synthetic mark, so the smoke guard stopped the run;
- an empty K0b table crashed `--combine` in smoke;
- the smoke world first had one library column per donor and cohort-specific naive and EM biology; it was rebuilt with
  one column per library and a shared biology, so that K2b and K5 test the code rather than the synthetic world.

**11.4 Changes after the code review** (2026-11-10, still before any real-data run; a re-run smoke chain covers them):
1. **Conformance fixes** (the code now does what sections 3, 6 and 7 say):
   - Merge rule: a stratum with no neighbour left is a deterministic FAIL only when pool < need (section 3). The code
     had failed it at pool < 3 × need.
   - The 7.2(iii) erratum goes after the "Reading" paragraph of block A1-S8, as registered. The code had anchored it on
     the `<!-- A1-S8 end -->` marker, above that paragraph. Smoke now applies the insert to a copy of GATES.md.
   - K0b: the job script sets `LC_ALL=en_US.UTF-8`. R 4.4.3 and the collation locale are a deterministic check in
     `--known`, and both are written to `k0b.csv`.
   - The 7.2(iv) sentence quotes the published and re-run P5 (T2) or k of n (T1).
   - 7.4 context now carries the A1-G01 category with its 11.5 sentence, G07 7.7's acute-illness sentence, the CMV
     sentence, and "no P370/F402 difference is attributed to the 32".
   - 3.3: the T1 rows carry the median coherence and median random P95 per null. T2 rows flag coverage < 0.8, which is
     reported in the summary.
   - K6: no AML per-label matrix is kept after loading.
2. **Known-answer worlds.** A world that could not be computed (a stratum FAIL, or fewer than 20 genes per half) is a
   deterministic FAIL. It is never counted as "did not pass", which would have helped the "at most" bars of K3, K5b
   and K7.
3. **Balance.** The SD in the K8 balance measure is taken over the matched universe the strata were built on, not over
   every gene with L. A balance that cannot be computed is "not readable", not readable-by-default.
4. **Reading precedence for (b)**, not stated in section 7:
   - If N1 fails, only the (iii) sentence is written; N1-L and N2 stay in the tables.
   - Next: N1-L not computable, then N1-L failing (11.2), then N2 not computable, then not readable, then (i)/(ii).
   - "Not computable" uses the not-readable sentence with "the loading-matched null could not be computed ([reason])".
   - The not-readable sentence names the check that failed (K5, K5b and/or K8 balance).
5. **(a) edge cases:**
   - The 7.1(i) replacement sentence uses only tests that passed and are readable.
   - If one test is not readable and the other is readable but did not pass, the 7.1(ii) sentence is written for the
     readable test only, with "the power to detect [that test]".
6. **The erratum lists all four N1 results** (T2 and T1, both AML cohorts) with pass or fail, so it shows why GATE_S8
   is PASS or FAIL under the per-half null.
7. **Engineering** (no effect on any statistic):
   - 16 workers, as in section 8 (the node had reported 32);
   - K5b drops all-zero genes before permuting (memory; this changes which random world is drawn, not its design);
   - K8's NA test no longer builds a 1e9-element vector (about 8 GB);
   - `--combine` checks the md5 of every input table;
   - the GATES block cites 7.3 instead of quoting the never-written phrases.
