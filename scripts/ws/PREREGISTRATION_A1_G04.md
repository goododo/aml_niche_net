# PRE-REGISTRATION A1-G04 — does the T-cell imprint direction hold where both arms share one protocol? (Chen2023; GSE185381 children as a secondary arm)

WP-A1 gap G04, options (a) and (b) of the gap assessment (2026-10-28, section 2.2), with two obligations handed on:
A1-BATCH 12.5 ("its pre-registration should read the delay score there too") and A1-G01 section 7 (G04(a)'s primary is
the whole-T direction and the residual is reported; G04(b) shares that primary and is read one-sided).

**Status: DRAFT v1 (2026-11-04). Nothing declared below has been computed.** Written from a read-only inventory of the
data and code (no expression value of Chen2023 or of the children was read; section 1.3 lists what was seen). Order
(user decisions 2026-10-28 to 2026-11-04): A1-G02 (PASS) → A1-BATCH (PASS with qualifications) → A1-G01 ("beyond T-subset
composition, detected inside the main labels") → **this** → G03. Handoff v2.0 §17 S15. The tag `prereg-a1-g04-v1` follows
the user's section-9 decisions.

## 0. Question, scope, and what each outcome lets the paper say

The imprint (S5: 370 primary genes, 156 up / 214 down) separates 15 adult AML from 10 healthy donors of GSE185381 out of
sample (A1-G02: cross-fitted gap 1.378 = 1.98 out-of-fold healthy-SD; full-pipeline permutation p 0.008), and most of the
gap lies inside T-subset labels (A1-G01: composition share C 0.045). In GSE185381 the healthy marrow is a commercial
cryopreserved product and the AML samples come from a hospital tissue bank; the healthy donors are about 31.5 years
younger; and the arms were partly captured on different dates. The within-study checks could not separate disease from
sample source and handling (A1-BATCH 12.5; A1-G01 section 0). **G04 asks whether the direction holds in a contrast where
both arms went through one protocol.**

- **(a) primary: Chen2023.** 5 NPM1-mutated AML at diagnosis (39-59 y) vs 4 allogeneic-transplant donors (47-53 y), one
  study. Every aspirate was cryopreserved, FACS-sorted into five fractions with one gating scheme and captured as two
  10x 3' v3 libraries per donor (one donor per library). T cells come from the lymphoid sort (pool B), which separated
  them from the blasts before capture. The ages are matched, so age cannot produce the direction here.
- **(b) secondary: GSE185381 children.** The 18 paediatric AML samples at diagnosis with at least 20 kept T cells
  (PAUMTZ excluded, 2.3), against the 10 healthy donors scored out of fold. Read one-sided.

What (a) changes compared with GSE185381, and what it does not:
- Shared by design in both arms: sorting scheme, cryopreservation, chemistry, library construction, and every step of
  this project's processing.
- Still different between arms: the sample source (diagnostic aspirates vs aspirates from transplant donors). No
  collection site, date or sort batch is recorded on disk. A1-BATCH 12.5 called Chen2023 "one centre"; no local source
  supports that, so this document says "one study, one protocol".
- Population: all 5 AML carry NPM1 mutations (4 of 5 also DNMT3A). A result speaks for NPM1-mutated AML only.
- Power: with 5 vs 4, a pass is informative and a fail is not (3.3). A significant reversal is informative.

What each outcome lets the paper say (rules in section 7):
- **(a) "direction replicated"**: "In an independent, age-matched cohort processed with one protocol in both arms
  (Chen2023, 5 NPM1-mutated AML vs 4 transplant donors), AML T cells score higher on the imprint than donor T cells
  (exact p = k/126)." This argues against a pure source-and-handling artefact of GSE185381. It does not separate
  disease from the remaining source difference (diagnostic vs donor aspirates). The qualifier of 3.2 is added where it
  applies.
- **(a) "not shown"**: "The direction was not shown in Chen2023; 5 vs 4 cannot exclude an effect of the size seen in
  GSE185381" (with 3.3's power). The direction stays "shown in one study".
- **(a) "reversed"**: "AML T cells score lower than donor T cells in Chen2023." The GSE185381 direction then does not
  reproduce where the protocol is shared, and the imprint must be treated as study-specific (section 7 handover).
- **(b) "AML-side"**: "T cells of children with AML sit on the AML side of out-of-fold healthy donors", which older age
  cannot explain (14 of the 18 children are under 18 y; the healthy donors are 19-53 y). **(b) "not readable"**: no
  statement.

What G04 cannot do, under any outcome:
- Separate disease from the remaining source difference in (a), or from storage time; Chen2023 has no handling records.
- Exclude blast-derived ambient RNA in (a): pool B still carries some myeloid or progenitor-labelled cells, more in AML
  (2.7-8.0 % vs 0.4-1.5 % of pool-B cells; seen 2026-11-04, 1.3).
- Generalise (a) beyond NPM1-mutated AML, or measure the residual's size: 5 vs 4 can support a direction, not a size
  (5.1 is descriptive).
- Test T-specificity: B and NK cells (5.3) are descriptive, as in A1-G01 5.1.
- Address source or handling with (b): the children share GSE185381's design and its healthy donors, and 12 of the 18
  were captured on dates with no healthy donor.

## 1. Already known before this draft (disclosed)

**1.1 The September pooled pseudobulk DE** (pre-registration a66cbd0, 2026-09-16; no tag). T_NK bin (projection-only,
NK included, both Chen2023 pools): 39 AML (5 Chen2023, 20 GSE185381 children, 14 GSE185381 adults) vs 14 healthy (4
Chen2023, 10 GSE185381); `~ dataset + arm` with duplicateCorrelation on library; hit = q < 0.05 and |logFC| ≥ 1. One
baseline hit (AC007952.4); Monte Carlo permutation p 0.102 and 0.320 (1,000 draws each, not exact); reported as 0
reportable hits; recovery 0.905 at |logFC| 1.5 in the planted curve. Its pre-registration promised not to reframe this
null as underpowered (`scripts/12_pseudobulk_de/PREREGISTRATION_pseudobulk_de.md`:41-44). G04 keeps the promise: it does
not re-test genes; it tests one pre-specified sample-level score, and G03 relates the two afterwards.
- Prediction on record (gap-assessment addendum 1.3, 2026-10-28): at S5's per-gene sizes the September fit should have
  called about 156 of the 370 genes if all 39 AML carried them, and about 4 if only the 14 adults did; it called 1.
  Either the imprint is absent or weaker in the children and the Chen2023 AML, or S5's per-gene sizes are inflated. **The
  first reading predicts that (a) and (b) fail.**
- The September gene universe was a symbol intersection. Five P370 genes carry different Ensembl IDs on the Chen2023 and
  GSE185381 axes, among them AC007952.4, so the single September T_NK hit pooled two loci (Chen2023 ENSG00000266311,
  GSE185381 ENSG00000262202). G04 matches genes by Ensembl ID (2.2).
- The September per-gene tables stay closed until G04 has a result (section 8); only AC007952.4's values have been seen
  (A1-BATCH section 1).

**1.2 Other Chen2023 results on record** (none used this imprint): the Chen-only LMPP_GMP DE (5 vs 4; 304 hits, exact p
2/126), of which NBM3 alone carried 101 by leave-one-NBM-out; HSC_MPP 0 hits; the leave-one-dataset-out transfer of the
CCC-era Discovery AML-vs-healthy node-feature axis had the opposite sign in Chen2023 (`scripts/DECISIONS_pending.md` P11);
FGW topology null; inferCNV malignant-call rates in presumed-normal NK/CD8 effector-memory cells (AML 0.159 vs healthy
0.256) and per-donor malignant fractions.

**1.3 Seen while preparing this draft** (2026-11-04, read-only inventory; notes outside the repo): kept T counts per
BoneMarrowMap fine label for all 10 Chen2023 donors (all pools and pool B) and for the 19 eligible children, so the
per-arm T-subset composition of Chen2023 is known to the drafting side; per-pool counts of every hierarchy bin (the
pool-B myeloid/progenitor shares in section 0); kept NK, B and monocyte counts; QC thresholds and retention per library;
per-donor inferCNV malignant fractions; the children's hashtag tallies (identity), author labels (age band, sex,
inflammation group); and three children's all-cell QC covariates (AML001, AML052, PAUMTZ). No expression value of an
imprint gene, no score and no AML-vs-healthy expression contrast was computed for Chen2023 or for the children.
Consequence: the composition readouts (5.1) are descriptive only.

**1.4 The children** were in the September fit (20 of them). GATE_B0 used their T_NK IFNG/TNF/TGFB1 values. The 79-gene
ambient blacklist and the 32 ambient-prone genes come from 47 decontX units, 19 of which hold children (label-free
gene lists). The out-of-fold scores of the three young pool controls (Control0005, Control0058, Control0082) were seen in
A1-BATCH section 1. Lasry 2023 (the GSE185381 paper): CD4/CD8 T cells depleted only in children; inflammation remodelled
the paediatric T compartment (Treg and GZMK+ CD8 expansion), not the adult one.

**1.5 GSE185381 facts that bear on the readings.** The four adult AML under 60 y score 1.07 / 0.48 / 0.83 / 0.53 against
a healthy range of −2.27 to −0.43 (in-sample; S9 3.4). The imprint score follows the CD8 share in GSE227903 (Spearman
0.86, n 7) but not in GSE185381 (0.17 AML / 0.22 healthy). The direction is shared by B and NK cells of the same donors
(A1-BATCH; A1-G01 5.1) and is aligned with a processing-delay program (A1-BATCH D8).

**1.6 Status in the handoff.** Chen2023 is a direction "reference only" (`HANDOFF_within_sample_v2.0.md`:102, 160, 357,
445) and carries the CCC-era Discovery label (`02_study_split.csv`:2). The children were planned as a separate second
cohort without healthy comparators (HO:247), paused on 2026-10-28 (HO:511). Section 9 asks for the decisions.

## 2. Inputs and frozen definitions (fixed before any expression is read)

### 2.1 Shared inputs

| item | file | md5 |
|---|---|---|
| imprint (P370 primary; 402 sensitivity) | `results/tables/ws/ws_s5_imprint_genes.csv` (compartment T, imprint_called) | 977f58b06da076cc39a561ab7445b69d |
| 32 ambient-prone genes (removed) | `results/tables/ws/ws_s8_t_imprint_ambient_prone.csv` | acb2b043049a4286e4052d04cea15feb |
| GSE185381 kept-cell pseudobulks | `/LARGE1/…/ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds` | 6c2fb685c49d4a23ef20a6842aa7aa67 |
| the 25 S5 samples and status | `results/tables/ws/ws_s5_scores.csv` (compartment T) | 58374fc60aabc53ba020c3149764fe8f |
| A1-G02 fold lists | `results/tables/ws/ws_a1_g02_fold_genes.csv` | ede3280052ab8dfcc547e60c341f334b |
| A1-G02 out-of-fold scores | `results/tables/ws/ws_a1_g02_oof_scores.csv` | b1da65a5ddd4d050a1b0849a9aa1206d |
| r = 1.430 (gap_in / gap_oof) | `results/tables/ws/ws_s9_yardstick.csv` | dcbea817a8eb6c68575486e3576dfce9 |
| ambient blacklist (n_units ≥ 5: 79 genes) | `results/tables/ws/decontx/GSE185381/ws_s1_ambient_blacklist.csv` | acebed4dbeed2ef6ecd5a60791e28fb9 |
| processing-delay set (295 up / 497 down) | `/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv` | 5ee444a9284ae52b7afab1c8aeeced2c |
| S4 label map | `scripts/ws/s4_cell_groups.tsv` | f5a5a0283b243e5efaf17b1ddfdda5f1 |

### 2.2 (a) Chen2023

- **Objects:** `/LARGE1/…/02_seurat_objects/04_annotated/Chen2023/<donor>.rds` (md5: AML103 5c45be9c…, AML162 44866549…,
  AML24 d4fe0773…, AML304 ff64ce4e…, AML320 e31fd77f…, AML469 4b08d3df…, NBM1 29bcb07c…, NBM2 5364a9c1…, NBM3 3876099a…,
  NBM4 de8c3cc4…; full values in the K0 table). Metadata `meta_Chen2023_v2.xlsx` (c2407360893ccf5261dff4e2144cd996);
  deposit `features.tsv.gz` (identical in all 20 libraries, 16952af81ef735f42952b64c8402f69a);
  `results/tables/00_ingest/GSE120221_gene_remap.csv` (d8d7c37eeb8cb85f5bd4b9d53d1f0acb).
- **Donor table** (metadata coding, declared now; written by `41 --build` from the xlsx and checked by K2):

| donor | arm | age | sex | drivers (VAF) | pool-B kept T |
|---|---|---|---|---|---|
| AML103 | AML | 59 | F | NPM1 0.42, DNMT3A, FLT3-TKD, PTPN11 | 211 |
| AML24 | AML | 59 | M | NPM1 0.41, DNMT3A, IDH2, NRAS | 1,047 |
| AML162 | AML | 51 | M | NPM1 0.34, DNMT3A, TET2 | 890 |
| AML304 | AML | 56 | M | NPM1 0.43, DNMT3A, TET2 | 1,445 |
| AML320 | AML | 39 | F | NPM1 0.45, NRAS | 808 |
| NBM1 | donor | 47 | M | — | 2,306 |
| NBM2 | donor | 52 | F | — | 1,094 |
| NBM3 | donor | 50 | M | — | 2,795 |
| NBM4 | donor | 53 | F | — | 1,709 |
| AML469 | AML (excluded) | 50 | M | NPM1, FLT3-ITD, STAG2, NRAS | 2 |

- **Cells:** pool B (meta.data `Pool == "poolB_stroma_lymphoid"`; library suffix Niche_Immune in four spellings), S4 kept
  rule (bmm_fine maps to group4 T in `s4_cell_groups.tsv`, `anno_source == "concordant"`, `high_error != 1`). Pool-A T
  cells (1-339 per donor) lie outside pool A's CD34+ and CD117+/CD33+ gates and were captured with the blasts; they are
  excluded (sensitivity 4.3 adds them).
- **Unit and eligibility:** the donor (one pool-B library each); at least 20 pool-B kept T cells. Eligible: 5 AML and 4
  donors (table). **AML469 is excluded:** its pool-B library lost 526 of 552 barcodes to a donor-level haemoglobin QC
  threshold of 0, leaving 2 kept T cells. It is not rescued by another QC or by the leukemia_ml atlas cells (a forking
  path closed now).
- **Gene axis:** row i of each object ↔ row i of `features.tsv` ↔ row i of `GSE120221_gene_remap.csv` (K1). Rows with rule
  `mapped_by_ensembl` (32,351) take their 2020-A name; every other row (1,343) is dropped, so no gene is matched by symbol
  alone. The primary set on this axis has **358 genes (154 up / 204 down)**. Missing: AC016831.7, HNRNPA1P48,
  AL021155.5, CCDC18-AS1, AC253572.2, AC239799.2, SATB1-AS1, AL499604.1, PRKCQ-AS1, RPARP-AS1, LINC01278 (no 2020-A ID on
  the Chen2023 axis), and POLR2J3 (Chen2023 ENSG00000168255 vs 2020-A ENSG00000285437; dropped by the ID rule). The 402
  set has 388.
- **Matrices:** per donor, raw count sums of pool-B kept T (primary), kept B and kept NK cells (5.3).

### 2.3 (b) GSE185381 children

- **Samples:** the 19 paediatric samples with at least 20 kept T cells (S4), minus PAUMTZ: AML001, AML003, AML005, AML006,
  AML009, AML012, AML022, AML025, AML026, AML028, AML038, AML043, AML048, AML051, AML052, AML055, AML056, AML073 (18).
  **PAUMTZ is excluded:** 17 of its 31 kept T cells carry the PAUMTZ_EoT hashtag as raw top hashtag, its expression sex
  contradicts its karyotype (46,XY), and reading it would read post-treatment cells that S10 keeps blind.
- **Matrices:** `PB$raw$T` columns of these 18 (author labels, demultiplexed by hashtag plus genotype). Never read: the
  UNDEMUX cells of the ingest object, the P06 and P12 lane matrices beyond their hashtag rows, and the PAUMTZ and PAXMIJ
  columns of any matrix.
- **Flags** (sensitivity 4.7; from metadata and hashtags only): age at least 18 y, or no Table S1 age and no author band
  below 18 (AML022 19.75, AML038 21.08, AML043 19.58, AML073 none); an age conflict (AML052: Table S1 4.5 y, author band
  "teens"); therapy-related AML (AML043); more than 5 % of kept T cells whose raw top hashtag belongs to another
  deposited sample (AML056 119/250, AML038 10/35, AML028 2/20, AML012 3/58); the P19 naming mismatch (AML026, AML051);
  fewer than 50 kept T cells (AML028 20, AML043 29, AML038 35, AML026 41). Flag-free: AML001, AML003, AML005, AML006,
  AML009, AML025, AML048, AML055 (8).

### 2.4 Scores (all label-free; S5 signs)

- **J space ((a) primary; as S9 2.4).** `black_for()` removed on each axis; U_J = names with CPM > 1 in at least half of
  the 25 GSE185381 samples and in at least half of the 9 Chen2023 donors; TMM within each cohort on U_J; log-CPM with prior
  count 1; z with the 25 samples' per-gene mean and SD. Score = mean z(up) − mean z(down) over P370 ∩ U_J; each half needs
  at least 5 genes.
- **W space (4.2).** `logz()` on the 9 Chen2023 donors alone (blacklist, CPM > 1 in at least half, TMM, prior 1), z over
  the 9.
- **Out-of-fold space ((b)).** L = `logz()` log-CPM of the 25 (as A1-G02 R2). The children are placed with `make_L27()`
  (`37_a1_batch.R`:108-118): each child is TMM-normalised against the 25's reference column, with edgeR's prior scaled to
  the 25's mean library size, so L on the 25 is unchanged. For healthy fold h (folds 16-25 of
  `ws_a1_g02_oof_scores.csv`): z = (L_x − mean of L over the 24 others) / SD over the 24 others, and the score uses fold h's
  lists. The children enter no filter, TMM or z parameter, so the stored healthy scores are reproduced exactly (K7).

## 3. Primary estimate for (a) and the declared readout for (b)

### 3.1 (a) D and its exact p

D = mean J score of the 5 AML − mean J score of the 4 donors (P370). Exact one-sided p = (number of the 126 labellings of
the 9 donors into 5 + 4 whose D is at least the observed D, within 1e-12) / 126; the observed labelling is one of the 126.
Reported with: the rank-sum count U (AML-donor pairs with the AML donor higher, of 20) and its exact p; the per-donor
scores; D in units of the donors' score SD (descriptive); **ρ = D / (D_res_J / r)**, the share of GSE185381's
selection-corrected gap, where D_res_J is the AML − healthy mean of the 25 GSE185381 J scores on the same U_J and r =
1.430.

### 3.2 (a) leave-one-donor-out and the qualifier

Exact p for each of the 4 leave-one-healthy-donor-out sets (5 vs 3; 56 labellings) and each of the 5 leave-one-AML-out
sets (4 vs 4; 70 labellings). **Qualifier "rests on one healthy donor"**: any leave-one-healthy-donor-out p above 0.10. In a generic
simulation (normal scores, no data) this flags about 22 % of chance passes and 4-10 % of true passes at 1-2 SD (NBM3
carried a third of the September LMPP_GMP hits, 1.2).

### 3.3 (a) power (generic, no data)

Normal scores with equal SD in both arms: the exact mean-difference test at p ≤ 6/126 passes with probability 0.05 /
0.37 / 0.62 / 0.84 at a true difference of 0 / 1 / 1.5 / 2 within-arm SD (rank-sum at most 2 inversions, p ≤ 4/126: 0.03
/ 0.28 / 0.50 / 0.75). GSE185381's out-of-fold gap is 1.98 out-of-fold healthy-SD; how that transfers to Chen2023's own
SD is unknown before the run, and the September null (1.1) suggests smaller. Hence: a pass is informative, a fail is not.

### 3.4 (b) A_b, the within-fold AUC

For each healthy fold h (10 folds), a_h = share of the 18 children whose fold-h score exceeds healthy donor h's
out-of-fold score (ties count 0.5); **A_b = mean of the 10 a_h**. Both sides are out of sample for fold h's lists. 95 % CI
by 2,000 two-sample bootstrap draws that resample the 18 children and the 10 healthy folds with replacement (seed SEED +
17). Reported with: each child's mean score over the 10 healthy folds; the placement φ = (the children's mean score over
all 25 folds − the healthy out-of-fold mean) / (the adult AML out-of-fold mean − the healthy out-of-fold mean); and the
same-pool differences (AML022 and AML052 vs Control0082 in fold 19; AML009 vs Control0005 in fold 17), without a test.

## 4. Sensitivities (each reports its statistic with an exact p or a CI; none changes a category)

- **4.1** (a) the 402 set (388 genes).
- **4.2** (a) W space.
- **4.3** (a) all pools: pool-A kept T cells added.
- **4.4** (a) the up half and the down half scored separately.
- **4.5** (b) children under 18 y: Table S1 age below 18, or no age and an author band below 18 (14: AML001, AML003,
  AML005, AML006, AML009, AML012, AML025, AML026, AML028, AML048, AML051, AML052, AML055, AML056).
- **4.6** (b) at least 50 kept T cells (14).
- **4.7** (b) flag-free children (8; 2.3).

## 5. Further declared readouts (descriptive; no category)

- **5.1 Composition and the residual** (A1-G01 handover). (i) Naive share ((CD4 Naive + CD8 Naive) / kept T) and CD8 share
  per arm. (ii) The counterfactual of A1-G01 3.1 inside Chen2023: `est_C` with Chen2023's own pool-B label counts and
  leave-one-out arm profiles (fallback below 50 cells), the P370 list for every donor, u-units in `logz()` over the 9;
  C_H, C_A and C, and the residual gapT − gapC in each order and for the mean. Composition was seen before this draft
  (1.3) and the counterfactual's references depend on the labels, so no test.
- **5.2 Processing-delay score** (A1-BATCH 12.5). The D8 score (mean z of delay-up minus delay-down genes, F402 genes
  removed) in the same J space: D_delay with its exact p (126 labellings), and Spearman with the imprint score within
  each arm.
- **5.3 B and NK cells of the same donors** (is the direction shared, as in GSE185381?). The J-space imprint score of pool-B kept B (5 vs 4) and NK
  cells (4 vs 4; AML320 has 18): D and exact p. Monocytes are not scored: they sit in pool A with the blasts.
- **5.4 Gene-level sign concordance.** The share of the 358 genes whose AML − donor difference in mean J log-CPM has S5's
  sign, against its distribution over the 126 labellings. Random gene sets are not used as a null (they pass 27-40 % of
  null runs at 5 vs 4; addendum 1.7).
- **5.5 (b) strata.** Children's scores by Lasry inflammation group (7 high, 11 low) and by author age band; median per
  group, no test.

## 6. Known answers (run first)

Scope P blocks (a)'s category, scope B blocks (b)'s. A failure leaves the modes running, labelled "not decisive,
descriptive", as in A1-G01.

| id | check | expected | scope |
|---|---|---|---|
| K0 | md5 of every declared input (2.1-2.3) | as declared | P, B |
| K1 | `features.tsv` (ID, symbol) equals the remap's (ens, sym_cr2) in all 33,694 rows in order; every object's row names equal `make.unique(gsub("_", "-", symbol))`; P370 on the remapped axis | identical; 358 (154 / 204); 402 set 388 | P |
| K2 | donor table equals the xlsx (arm, age, sex) for 10/10; `Pool` agrees with the library name for every kept T cell; pool-B kept T per donor equals 2.2; each pseudobulk column sum equals the sum of its cells' counts | exact | P |
| K3 | expression sex from XIST and RPS4Y1 UMIs of each donor's pool-B kept cells (ws/26's rule) | 9/9 equal to the metadata coding; a mismatch gets a dated note before the primary is read (the Control4 precedent), not a FAIL | reported |
| K4 | the J-space code reproduces S9: D_res_J for GSE120221 from S9's inputs, and gap_in on the 25 in S5 space | 1.968 and 1.9714 within 0.001 | P |
| K5 | enumeration: 126 distinct labellings, the observed once; the rank-sum exact p equals `pwilcox` | exact | P |
| K6 | synthetic plant (no Chen2023 expression): 200 synthetic 5 vs 4 data sets from the pooled GSE185381 healthy T profile with Chen2023's pool-B library sizes; (i) no plant; (ii) S5 log2FC planted at full size into the 5 synthetic AML | (i) pass rate ≤ 0.08; (ii) pass rate ≥ 0.90 | P |
| K7 | `make_L27` keeps L on the 25 within 1e-10; the stored out-of-fold scores of all 25 folds are reproduced within 1e-6; every child's column is finite; fold-gene coverage of L reported | exact | B |
| K8 | the 18 IDs and their kept T counts equal `ws_s4_pseudobulk_meta_GSE185381.csv` | exact | B |

## 7. Decision rules (first match wins)

**(a):**
1. **"not decisive"**: a scope-P known answer failed.
2. **"direction replicated"**: D > 0 and exact p ≤ 6/126. Add "rests on one healthy donor" if 3.2's qualifier applies.
3. **"reversed"**: D < 0 and the mirror exact p (labellings with D at most the observed) ≤ 6/126.
4. **"not shown"**: otherwise.

**(b):**
1. **"not decisive"**: a scope-B known answer failed.
2. **"AML-side"**: the lower 95 % bound of A_b above 0.5.
3. **"not readable"**: otherwise.

(a) and (b) are read separately; there is no combined statistic, and (b) never changes (a)'s category. Sensitivities and
readouts never change a category.

Handover (fixed now):
- G03 runs after G04's result. The September per-gene tables are opened there, not before.
- **(a) "direction replicated"** is offered as the direction half of GATE_A3 ("at least 1 other study significantly
  concordant", HO:169) if the user agrees (9.6), with the qualifier and the source sentence of section 0.
- **(a) "reversed"**: the direction claim stops; discuss before G07 or any writing.
- **(a) "not shown"**: the direction stays "one study"; G04(c) remains a separate decision (it cannot separate disease
  from source either, gap-assessment addendum 1.1).
- S10 is unchanged (it inherits the whole-T score from A1-G01).

## 8. Code, outputs, compute (plan; scripts are written after the tag)

- **`scripts/ws/41_a1_g04.R`**, modes: `--build` (the Chen2023 donor table from the xlsx; pool-B kept T, B and NK
  pseudobulks with the remap, saved to `/LARGE1/…/ws_a1_g04/`; no expression-derived table in the repo), `--known` (K0-K8),
  `--primary` (3.1-3.2), `--sens` (4.1-4.4), `--children` (3.4, 4.5-4.7, 5.5), `--readouts` (5.1-5.4), `--combine` (summary,
  GATES.md block A1-G04). `--smoke`: synthetic counts and synthetic labels for both arms with real IDs only; no real age,
  sex or pool fraction is written next to a score; `*_smoke` outputs are deleted after the real run.
- **`scripts/ws/submit_ws_a1_g04.sbatch`**: the A1-G01 pattern (dated read-only copy, md5 in every table,
  `OMP_PROC_BIND=false`); about 4 CPUs, 64 GB, 12 h.
- Functions are copied, not sourced: `logz`, `black_for`, `jspace`, `zE_from_counts` (33), `make_L27` (37), `est_C` and
  helpers (39).
- **Not run on Chen2023:** ws/22 (decontX; its unit would merge the two 10x runs, and it is sensitivity-only by config),
  ws/25 (rewrites cohort-wide S3 tables), ws/26 (no pool filter; writes per-arm composition and sex calls to the repo),
  ws/02 (sorted AML are ineligible). No re-ingest: the remap is applied to the pseudobulks. `config_ws.R` is unchanged. No
  new dependency.
- Seeds: SEED + 17 ((b) bootstrap), + 18 (K6), + 19 (smoke); none is used elsewhere.
- Outputs: `results/tables/ws/ws_a1_g04_{checks,primary,loo,sens,children,readouts,summary}.csv`; GATES.md block A1-G04.
- **Never opened before G04's result** (replaces A1-BATCH's "any Chen2023 file" rule): `results/tables/12_pseudobulk_de/`
  per-gene and sample-level tables (`de/`, `hitlist_frozen.csv`, `sensitivity_pooled_bh.csv`, `depth_struck_hits.csv`,
  `perm_null.csv`, `secondary_*.csv`, `validation/`, `confound_inferred_sex.csv`) and `/LARGE1/…/12_pseudobulk_de/`
  pseudobulks; LCC_proj pseudobulks; `03_bmm_projected/*/*__cellstate_percell.csv` and `*__progeny_percell.csv` (per-cell
  CXCR4 and state scores) of Chen2023 and the children; `05_ccc/ccc_node_features*.csv`, `07_fgw*/fgw_nodes_long.csv` and
  the `08_scoring` panel-screen tables; `ws_b0_inputs*.csv`, `ws_t3c_preflight_*`, the truth-audit stage-2 results; the
  children's xist/rps4y1/sex_inferred/dbl_mean/median_umi columns of the S4 meta table; and the UNDEMUX cells, P06/P12 lane
  matrices and PAUMTZ/PAXMIJ columns of 2.3.

## 9. Decisions for the user (recommendation first)

1. **Chen2023's role.** Promote it from direction "reference only" to the G04(a) test; its CCC-era Discovery label does
   not bind. *(Recommended.)*
2. **(a)'s statistic.** The exact mean-difference test, pass at p ≤ 6/126 (power 0.37 / 0.62 / 0.84) *(recommended)*; or
   the rank-sum, pass at most 2 inversions (0.28 / 0.50 / 0.75; more robust to one extreme donor, effective alpha 0.032).
3. **The children.** Keep (b) as a secondary one-sided arm with the out-of-fold comparator *(recommended)*; or drop (b)
   from G04.
4. **PAUMTZ** excluded from (b). *(Recommended.)*
5. **Scope.** G03's rungs that add Chen2023 or the children go to G03's own pre-registration after G04's result
   *(recommended; the addendum had placed them here)*; G10(a) and G14(a) are not bundled; G04(c) stays separate.
6. **GATE_A3.** An (a) "direction replicated" counts as the direction half ("at least 1 other study significantly
   concordant") *(recommended, with the qualifier and the source sentence)*.

## 10. Review history

- DRAFT v1, 2026-11-04: written from the read-only inventory (workflow wf_e07a2214-567; notes outside the repo). Not yet
  reviewed.
