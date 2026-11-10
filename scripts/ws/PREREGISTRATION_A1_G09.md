# PREREGISTRATION A1-G09, stage 1: is S8's replicated structure specific to AML?

**Status: DRAFT v1 (2026-11-10). Not approved, not tagged. Nothing declared in sections 3-6 has been computed.**
The user fixed scope and five design choices on 2026-11-10 ("按推荐", section 9 items 9.1-9.6). Section 9 lists the
decisions that remain open.

**Sources.** All are kept outside the repo, under `/FAST/gr10634/gaozy/reports/2026-11_a1_g09_prereg/`:
- `inventory/inventory.json`: a read-only inventory with three lenses and a critic. The lenses are A (obligations and
  blindness), B (machinery and data) and C (design and power).

**Definition sources:**
- the paper-gap assessment (`reports/2026-10_wpa1_paper_gap/WPA1_PAPER_GAP_ASSESSMENT_2026-10-28.md`, G09 row at :69
  and the text at :160-161);
- A1-G01 section 5.4, the carried note at :282-286;
- the G01 review-4 v2 rule (`reports/2026-11_a1_g01_prereg_review/review4/A/v2.md`:189-200 and K6 at :237-238);
- the open risk in `reports/2026-11_T_shift_progress/writable_inventory_2026-11-10.md`.

## 0. Question and what each outcome lets the paper say

S8 showed two properties of the P370 axis in adult AML cohorts that played no part in defining it:
- **T2:** across patients, the mean z of the up genes and the mean z of the down genes move oppositely.
- **T1:** inside one patient's T cells, the axis genes co-vary more than matched random genes.

The result: PASS in GSE289435, consistent in GSE227903, not reproduced on Seq-Well GSE116256.

Those cohorts have no healthy donors. A gene set lined up with T-cell differentiation would pass T2 and T1 wherever
T-cell subtypes vary, because subtype heterogeneity alone produces both (G01 review 4).

**G09 stage 1 asks two things:**
- **(a)** Do T2 and T1 also pass among healthy donors?
- **(b)** In the AML cohorts, do they survive random sets matched on the naive-to-effector loading as well?

| outcome | the paper may write |
|---|---|
| T2 or T1 passes among the 20 Oetjen donors | "The axis genes also form a coherent T-cell program in healthy marrow (Oetjen, 20 donors); S8 shows that this program is reproducible, not that it is specific to AML." S8 may not be cited as evidence of an AML-specific program (G01 review 4). |
| Neither passes among the Oetjen donors | "The structure seen in two AML cohorts was not detected among 20 healthy donors, a larger set than either AML cohort." It is still not direction, not causation and not T-specificity. |
| (b): the AML pass is lost under the loading-matched null | "S8's structure is consistent with the naive-to-effector differentiation axis." |
| (b): the AML pass survives the loading-matched null | "S8's structure exceeds that of random gene sets matched on the naive-to-effector loading." |

**What G09 is not.**
- It is not a test of direction. G02 holds the AML-versus-healthy level difference, and G04 the independent direction.
- A pass among healthy donors does not weaken the level difference. G02 measures a level, and S8 measures covariance.
- It settles neither cause nor age. The Oetjen donors are 24-84 y, and S9 bounded the age share of the level
  difference (0.102, upper bound 0.463).

**Stage 2, outside this pre-registration** (decision 9.1). Two parts follow as their own short amendment:
- the power of S8 at n 7-8, estimated with planted structure;
- the Seq-Well boundary: GSE185381 AML samples thinned to Seq-Well depth.

The AML side of the Seq-Well result is already public (T2 positive, T1 3 of 5; S8 section 8), so stage 2 cannot be blind
on the AML side. Its depth target comes from measured counts, not from the assessment's ~1,900 UMI (decision 9.4). It
will read: pass probability < 0.2 means "the Seq-Well failure is expected at that depth (uninformative)"; ≥ 0.8 means
"the failure is evidence against generality or of a platform effect".

## 1. Already known before this draft (disclosed)

1. **S8's AML statistics are public.** They are in `ws_s8_tests.csv` and the GATES.md block A1-S8:
   - GSE289435: T2 −0.55 (random P5 0.38), T1 8 of 8;
   - GSE227903: T2 −0.43 (0.21), T1 7 of 7;
   - the discovery cohort: T2 −0.29, T1 13 of 15.

   Part (b) on the AML cohorts is therefore blind only to the new, loading-matched null.
2. **No healthy-donor T1 or T2 exists anywhere.** Every S8 table holds AML_Dx samples only. S9's Oetjen table
   (`ws_s9_scores_GSE120221.csv`) has one combined score per donor and no up-half or down-half means. It cannot give T2,
   and its score columns were not opened for this draft (inventory A).
3. **Partial prior knowledge about Oetjen.** S9 found a weak age association in Oetjen that runs through the naive-T share
   (f_slope 0.102 above the random P95 0.046). Donor composition therefore partly varies with age. Composition is exactly
   what can drive T2 across donors, and what (b) controls for.
4. **The loading has been used before.** G01 5.5 used the Oetjen naive-to-effector loading and printed union
   correlations of 0.355 and 0.313 (A1-G01 11.4). The loading as a per-gene vector was never written to a table.
5. **Only counts were read for this draft:** cell numbers, donor rosters, library-to-donor maps and median UMI.
   - Oetjen: 20 donors, all eligible, 104 to 3,721 kept T cells; donors C and S pool 3 and 4 libraries.
   - GSE185381 healthy donors: 10, with 123 to 1,917 kept T cells.
   - GSE227690 (Griffin): 5 donors over 13 libraries.
   - GSE276360: 4 donors.
   - GSE116256: 2 healthy samples with ≥ 20 T cells (113, 727).

## 2. Data and frozen definitions

**Gene sets.**
- Primary: P370 = S8's "primary" (`32_s8_replication.R`:38-39 checks 370 / 156 / 214).
- Sensitivity: F402 = S8's "full". It is reported as numbers only.

**Cohorts and their roles.**

| cohort | samples | role |
|---|---|---|
| Oetjen GSE120221 | 20 healthy donors; libraries pooled per donor as in S9 (`33_s9_age_direction.R`:110-126, `donor_ages_GSE120221_v1.0.csv`) | **decides (a)** |
| GSE185381 healthy | 10 donors; also run without Control1-5 (T-ALL lanes; 5 donors) | reported, not used: these donors were in the gene selection |
| GSE227690 (Griffin) | 5 donors, libraries pooled per donor | T1 only, supporting |
| GSE276360 | 4 donors | T1 only, supporting |
| GSE116256 healthy | 2 Seq-Well samples | T1 only, descriptive |
| GSE289435 AML_Dx | 8 patients (S8 rep) | (b) |
| GSE227903 AML_Dx | 9 in the manifest, 7 eligible (3853_Dg 0 and 6323_Dg 12 T cells) | (b); subset right after loading (K6) |

T2 needs ≥ 4 samples (S8 rule), which is why the Griffin, GSE276360 and GSE116256 healthy runs get T1 only.

**Cells and genes.** These are S8's own, unchanged (`32_s8_replication.R`:23, 47-61, 70-82):
- kept T cells (`group4 == "T"`), with ≥ 20 per sample;
- the blacklist;
- CPM > 1 in at least half of the run's samples;
- z-scores per gene across the run's samples;
- T1: a cap of 500 cells per sample, genes detected in ≥ 10 % of those cells, and ≥ 20 genes per half.

**Pooling by donor.** For a donor with several libraries (Oetjen C and S; Griffin), pseudobulks are summed. For T1,
cells are pooled across libraries before the 500-cell cap.

**The naive-to-effector loading** (G01 5.5). For each gene:

L = log2((CPM_EM + 1) / (CPM_naive + 1))

- EM = the pooled kept-cell pseudobulks of CD8 Effector Memory 1, CD8 Effector Memory 2 and CD4 Effector Memory.
- naive = the pooled CD4 Naive and CD8 Naive pseudobulks.
- The source is Oetjen (`ws_pb_GSE120221.rds`, `$fine`), with BoneMarrowMap labels. Oetjen is 25 % of BoneMarrowMap, so
  its labels are in-sample. This is disclosed, not corrected.

Which donors the loading comes from (decision 9.3):
- AML runs: all 20 Oetjen donors.
- Oetjen T1 for donor d: the other 19 donors (leave-one-donor-out).
- Oetjen T2: the GSE185381 healthy donors. T2 is one statistic over all 20 donors, so a leave-one-donor-out loading is
  not defined for it (decision 9.7). The all-Oetjen loading is a sensitivity.

## 3. Analyses

### 3.1 (a) T2 and T1 among healthy donors, with S8's null
- Run S8's T2 and T1 unchanged on each healthy cohort.
- Null: 1,000 random sets for T2 and 200 per donor for T1, matched on S8's expression × detection bins
  (`32_s8_replication.R`:52-57, 105).
- **Rule** (G01 review 4, as coded in S8 at `32_s8_replication.R`:92 and :109):
  - T2 passes if rho < 0 and below the random P5.
  - T1 passes if n_above ≥ min(6, n_ok) and n_above / n_ok ≥ 0.75. With 20 donors only the 75 % part binds.
- Read on Oetjen, P370. Every other cohort, and F402, are reported.
- **Why n = 20 is used** (decision 9.2). The question is whether the structure exists in healthy marrow. A larger n makes
  a healthy pass easier, which is the conservative direction for any AML-specific claim. A healthy fail at n = 20 is
  also more informative than one at n = 8.
- Reported alongside, descriptive only: 1,000 random draws of 8 of the 20 Oetjen donors (seed SEED + 901). For each draw,
  T2 rho and pass, and the share of draws that pass. No margin threshold is applied (decision 9.2).

### 3.2 (b) The loading-matched null, for T2 and T1
- Random sets are matched on S8's bins × terciles of L. The terciles are computed within each run's tested universe.
- Run on GSE289435 and GSE227903 AML_Dx, which is the reading of (b), and on Oetjen, reported.
- Same pass rules as 3.1. T1 uses the same loading strata within each donor's detected genes (decision 9.3).

**Reading of (b)** (AML cohorts, P370), per test:
- pass under S8's null and fail under the loading-matched null → "consistent with the naive-to-effector axis";
- pass under both → "exceeds loading-matched random sets";
- fail under S8's null → nothing to read, because this does not happen for the AML statistics already known (1.1).

### 3.3 Reported for every run
- n per run, and the coverage of P370 in the tested universe (S8 requires ≥ 0.8);
- the median T1 coherence and the median random P95;
- per-donor T1 values in a work table outside the repo;
- donor age next to T2 for Oetjen: a Spearman correlation of the per-donor up-minus-down score with age, descriptive.

## 4. Sensitivities (reported; none changes a reading)

1. F402 in place of P370.
2. Oetjen T2 with the all-Oetjen loading in place of the GSE185381 loading.
3. GSE185381 healthy donors, all 10 and without Control1-5.
4. Oetjen without the donors that S9 excluded under its own K3.

## 5. Multiplicity

Two tests (T2, T1) × the readings above. No p values are combined. Each reading is pass or fail against its own null, as
in S8. The supporting cohorts never change a reading.

## 6. Known answers (run first; no readout is printed until each has a state)

States are as in G07: PASS / FAIL deterministic (a code bug, fixed with a dated note) / FAIL statistical (the named
reading becomes "not readable").

| id | check | bar | type |
|---|---|---|---|
| K0 | Re-running S8's code path on GSE289435 and GSE227903 reproduces the observed T2 rho and the per-patient T1 coherence of every patient with ≤ 500 T cells (no subsampling) | exact (< 1e-12) | deterministic |
| K1 | The Oetjen roster: 20 donors; library-to-donor map equal to S9's; per-donor kept T cells equal to `ws_s9_scores_GSE120221.csv` `n_kept_T` (a count column) | exact | deterministic |
| K2 | The loading: the per-gene L recomputed from the pooled fine pseudobulks reproduces G01 5.5's union correlations 0.355 and 0.313 | ± 0.001 | deterministic |
| K3 | Null world for T2: on Oetjen, donor labels permuted independently per gene (structure destroyed), 100 permutations | T2 passes in ≤ 10 of 100 | statistical ((a) T2 "not readable") |
| K4 | Positive world for T2: on a K3 permutation, a donor factor f ~ N(0, 1) planted as up × 2^(0.3 f), down × 2^(−0.3 f) on the P370 genes, 20 plants | T2 passes in ≥ 16 of 20 | statistical ((a) T2 "not readable") |
| K5 | Composition-only world: synthetic donors built only by mixing Oetjen's pooled naive and EM profiles at random proportions; the test set is the 100 genes most EM-high and the 100 most naive-high by L; 20 worlds | under S8's null T2 passes in ≥ 16 of 20 (the confound is real); under the loading-matched null it passes in ≤ 4 of 20 (the matching works) | statistical ((b) "not readable") |
| K6 | GSE227903: exactly the 9 AML_Dx samples of the manifest; non-diagnosis columns dropped right after loading; z-scores and the universe built from the 7 eligible samples only | exact | deterministic |
| K7 | Null world for T1: cells permuted per gene within one donor (Oetjen donor K, 500 cells), 20 permutations | coherence above its P95 in ≤ 3 of 20 | statistical ((a) T1 "not readable") |

K3-K5 and K7 exercise the no-result branches (the G07 lesson: smoke and known answers must include a null world).

Run-time blindness: no healthy T2 or T1 value, and no loading-matched AML value, is printed until K0-K7 have states.

## 7. Reading rules and fixed sentences

**7.1 (a)** (Oetjen, P370, S8's null):
- **T2 or T1 passes:** "The axis genes also form a coherent T-cell program in healthy marrow (Oetjen, 20 donors; [T2
  rho r, random P5 p] / [T1 k of n donors above their random P95]). S8 therefore shows that this program is
  reproducible across AML cohorts, not that it is specific to AML."
- **Neither passes:** "Among 20 healthy Oetjen donors, the structure was not detected (T2 rho r, random P5 p; T1 k of n
  donors). This is more donors than either AML cohort (8 and 7)."
- Every (a) sentence adds: "Structure is covariance, not level: this does not bear on the AML-versus-healthy difference
  (A1-G02) or its direction (A1-G04)."

**7.2 (b)** As in 3.2, per test and per AML cohort. Both cohorts must agree for the sentence to name S8's structure as a
whole. Otherwise each cohort is written separately.

**7.3 Never written.**
- "AML-specific program" from S8 if 7.1's first case holds.
- "S8 replicates the direction".
- "Replicated on Seq-Well".
- "Not explained by differentiation" unless both T2 and T1 survive (b) in both cohorts.
- Any causal or age-excluding wording.
- Every never-written item of S8, S9, A1-G01, A1-G04 and A1-G07 7.8 still applies.

**7.4 Context that travels:**
- S8's T3 is uninformative;
- Oetjen's labels are in-sample to BoneMarrowMap;
- GSE185381 Control1-5 are T-ALL lanes and Control4 is female by expression;
- the S9 age bound.

## 8. Code, outputs, compute

- **Script:** `scripts/ws/45_a1_g09.R`, a new script that copies S8's functions with their sources and leaves
  `32_s8_replication.R` untouched.
  - Modes: `--known`, `--run`, `--combine`, each also with `--smoke`.
  - The guards are G07's: unknown arguments stop; real runs only from the dated read-only copy of a committed script;
    no overwrite.
- **Smoke runs** use synthetic counts on the real IDs and include K3-K5 and K7 null worlds.
- **Outputs:**
  - `results/tables/ws/ws_a1_g09_{checks,tests,summary}.csv`;
  - per-donor work tables under `/LARGE1/.../ws_a1_g09/`;
  - a GATES.md block "A1-G09".
- **Seeds:** SEED + 901 (subsampling), 902 (S8 nulls), 903 (loading-matched nulls), 904 (K3-K5), 905 (K7). One
  stream per offset.
- **Compute:** under 2 h on one `c=16, m=64G` job. Reading the 25 + 13 + 4 + 10 + 2 per-sample objects dominates. No
  new package.

## 9. Decisions

Fixed by the user on 2026-11-10 ("按推荐"):
- **9.1** Stage 1 = (a) + (b); power and the Seq-Well boundary become stage 2, a separate amendment.
- **9.2** Oetjen n = 20 with S8's rule is primary; the 8-of-20 draws are descriptive with no margin threshold.
- **9.3** The loading for Oetjen is leave-one-donor-out; the matched null applies to T1 and T2.
- **9.4** The Seq-Well depth target comes from the measured GSE116256 median UMI (stage 2).
- **9.5** The paper wording for each outcome is fixed now (sections 0 and 7).
- **9.6** GSE185381 healthy donors are reported only, with and without Control1-5. Griffin and GSE276360 are T1 support;
  the GSE116256 healthy samples are descriptive.

Still open (recommendation first):
- **9.7** The Oetjen T2 loading. Use the GSE185381 healthy donors' loading, because leave-one-donor-out is undefined for a
  cohort-level statistic; the all-Oetjen loading is a sensitivity. Alternative: use the all-Oetjen loading and state the
  circularity.
- **9.8** The K3/K4/K5/K7 bars as above. Alternative: tighter bars (K3 ≤ 5 of 100).
- **9.9** F402 as a sensitivity only, as in S8.

## 10. Review history

- 2026-11-10: inventory (four agents, read-only). DRAFT v1 written. No review yet.
