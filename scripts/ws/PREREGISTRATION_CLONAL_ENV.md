# Pre-registration: clone-tracking of AML-associated monocyte programs (CLONAL_ENV, frozen 2026-11-13)

Status: FROZEN 2026-11-13 on the user's approval, before any expression value was read for this question.
What had been read: cell counts by genotype label and cell type, mutation names per patient, and the CloneTracer
model inputs (which loci, which CNVs). Gene lists: `scripts/ws/gene_sets_clonal_env_v1.tsv` (md5 in the tag
message), built by `scripts/ws/ws_clonal_env_gene_sets.R`. Changes after the tag are dated appended sections only.

## 1. Question

Within one AML patient at diagnosis and within one mature compartment (monocytes / dendritic cells), do the
programs attributed to "AML-associated" monocytes differ between genotype-confirmed clone members and
wild-type cells of the same patient?

- If a program is higher in mutant cells than in same-patient wild-type cells, it travels with the clone.
- If it is equal in both, the program is not a clone property: it is either shared by the marrow (an
  environmental or systemic effect) or not AML-related at all. This design cannot tell those two apart without
  a same-protocol healthy arm (secondary analysis S2, one study only).

The patient is the block, so source, handling, age, capture period and ambient RNA cancel in the primary
contrast.

## 2. Data and eligibility (fixed now; rule chosen after cell counts, disclosed as such)

| source | platform | patient | driver used for labels | known answer |
|---|---|---|---|---|
| Petti 2019 | 10x 5' v1 | 721214 | NPM1 | HOX (K1) |
| Petti 2019 | 10x 5' v1 | 548327 | NPM1 | HOX (K1) |
| van Galen 2019 GSE116256 | Seq-Well | AML556 | NPM1 W288fs (+DNMT3A, NRAS, TET2) | HOX (K1) |
| van Galen 2019 GSE116256 | Seq-Well | AML419A | NPM1 W288fs (+DNMT3A, FLT3) | HOX (K1) |
| CloneTracer 2023 | 10x 3' + targeted | A.13 | NPM1, PTPN11, LMO2, DNMT3A, mtDNA 12539G>A | HOX (K1) |
| CloneTracer 2023 | 10x 3' + targeted | B.1 | KRAS, IDH2, NRAS (+chr8, chr3-part CNV) | chr8 dosage (K2) |
| CloneTracer 2023 | 10x 3' + targeted | B.4 | mtDNA 7075G>A (+chr7 CNV) | chr7 dosage (K2) |

Eligibility rule (applied on counts only): at diagnosis, >= 20 mutant and >= 20 wild-type cells in the
Mono/DC compartment, and the patient's T-cell mutant fraction <= 0.10. Excluded by the rule: 782328 (T 0.105),
AML921A (0.47), A.11 (0.12). Sirenko GSE230559 is excluded as a source because its GoT calls fail the T-cell
control and no read-level calls are deposited. Lilljebjorn 2025 has no per-cell calls.

Definitions:
- "Wild type" means no read of the leukaemic driver. In patients with a DNMT3A founder (AML556, AML419A, A.13),
  wild-type cells may still carry DNMT3A (clonal haematopoiesis). The contrast is therefore
  "leukaemic clone versus not leukaemic clone", not "mutant versus germline".
- Compartment: CD14 monocytes, CD16 monocytes and conventional DC, from cell-type labels that do not use
  genotype. van Galen "Mono" and "Mono-like" are merged (the "-like" suffix uses genotype). Petti and van Galen
  use this project's BoneMarrowMap projection where the sample is in the cohort; CloneTracer uses the authors'
  ct / ct_simple labels. pDC and progenitors are excluded.

Stage 0 (before expression): recount the eligibility table inside this compartment definition (Petti counts
in PROGRESS.md are myeloid-lineage, not Mono/DC). Any patient below the floor drops. Fewer than 5 remaining:
STOP.

## 3. Genotype labels and misclassification

- Petti and van Galen: mutant = at least 1 mutant read at the driver locus. Wild type = 0 mutant reads and at
  least 2 reference reads (primary); >= 1 reference read is a sensitivity arm.
- CloneTracer, primary: the authors' status (leukemic / healthy; unsure dropped). Their model uses
  expression-derived CNV counts (B.1 chr8 and chr3, B.4 chr7), so for these patients genes on those
  chromosomes are removed from every program score.
- CloneTracer, for the known answer K2 only: relabel cells from SNV and mtDNA reads alone (the M and N
  matrices in figshare 21982496), never using CNV counts, so that chromosome dosage is an independent check.
- Wild-type contamination: at a heterozygous site a true mutant cell with k reference reads and no mutant
  read is missed with probability 0.5^k. Per cell, c = pi * 0.5^k / (pi * 0.5^k + (1 - pi)), with pi the
  patient's mutant fraction in the compartment. The corrected difference is (mean_mut - mean_wt) / (1 - c_bar).
  For CloneTracer, c_bar = mean(1 - P(healthy)) over cells called healthy.
- Mutant contamination (ambient or barcode swap) is bounded by the T-cell mutant fraction.
  Both errors shrink differences toward zero; neither can create one. The known answer measures how much.

## 4. Scores and statistic

- Programs, primary: this project's MP5 (HLA-II), MP6 (non-classical / IFN) and MP10 (inflammatory
  secretory), gene lists as in `results/tables/04_cnmf/malignant_metaprograms.tsv` (62, 81 and 87 genes).
- Secondary feature panel (user 2026-11-13: "include more features, e.g. stemness"). Gene lists frozen in
  `scripts/ws/gene_sets_clonal_env_v1.tsv` with sources, md5 in the tag message. A feature without a citable gene table at
  freeze time is dropped, not improvised.
  Stemness: F1 LSC17 (Ng 2016, weighted; bulk-derived, noisy per cell); F2 van Galen normal HSC/Prog (50 genes);
  F3 van Galen malignant HSC/Prog-like (50 genes); F4 HSPC_core (14 genes). F1-F4 are already curated in
  `scripts/以前06_hierarchy/stemness_signatures.tsv`. F2/F3 were derived from van Galen's own cells, so for AML556
  and AML419A they are in-sample; their pooled estimate is also reported without those two patients.
  Myeloid state: F5 HALLMARK interferon-alpha response; F6 HALLMARK interferon-gamma response; F7 HALLMARK TNF-alpha
  via NF-kB (overlaps the handling-stress program MP9; flagged); F8 TET2-clonal-haematopoiesis monocyte inflammatory
  genes (DROPPED at freeze: the source table, Heimlich/Bhat 2022 Supplementary Data S7-S8, was not obtainable; may
  join later as a dated amendment); F9 van Galen malignant Mono-like program; F10 inhibitory ligands (LILRB4, CD274, PDCD1LG2, LGALS9, CD47, VSIR).
  HALLMARK sets from the pinned local GMT (h.all.v2024.1), never fetched at run time.
  Cell-level features (not gene sets): F11 subtype composition (CD14 / CD16 / cDC shares, mutant versus wild type);
  F12 maturation position on the monocyte trajectory (BoneMarrowMap pseudotime for Petti / van Galen with the
  placeholder-zero cells removed; authors' pseudo_myel for CloneTracer); F13 cell-cycle score.
- Score: within patient and compartment, log-normalised expression, genes detected in >= 5 % of the
  compartment's cells (the filter is set without looking at genotype; rarely detected genes are mostly zeros and
  would let a depth difference between mutant and wild-type cells pass for a program difference), z-scored per gene, mean over the program's genes.
- Per patient and program: linear model score ~ genotype + subtype (CD14 / CD16 / cDC) + log(nUMI).
  Effect d = genotype coefficient divided by the residual SD, then corrected for wild-type contamination.
  Subtype is adjusted because MP6 is in part the CD16-monocyte identity and MP5 the DC identity; if mutant
  cells are skewed toward one subtype, the raw difference would be composition.
- Per-patient p: 1,000 permutations of genotype labels within subtype.
- Pooled: random-effects meta-analysis of d across patients (DerSimonian-Laird), 95 % CI.

## 5. Gates (in order; any failure stops the line for discussion)

- GATE_S0: at least 5 patients pass Stage 0.
- GATE_K1 (known answer, NPM1 patients): HOX score (HOXA5, HOXA7, HOXA9, HOXA10, HOXB3, HOXB5, MEIS1, PBX3)
  is higher in mutant than wild-type Mono/DC, one-sided permutation p < 0.05, in at least 4 of the 5 NPM1
  patients. If it fails, the labels cannot be told apart from the biology in this compartment: STOP.
  HOX in the HSPC compartment is reported alongside, as a diagnostic only.
- GATE_K2 (known answer, B.1 and B.4): with SNV/mtDNA-only labels, mean relative expression of chr8 genes is
  higher (B.1) and of chr7 genes lower (B.4) in mutant cells, permutation p < 0.05 each. A failure drops that
  patient from the primary analysis; it does not stop the line.
- GATE_N1 (shuffle): genotype labels shuffled within patient and subtype, 200 times; the per-patient test at
  alpha 0.05 rejects in 0.03-0.07 of shuffles, and the pooled CI covers 0 in >= 93 % of shuffles.
- GATE_P (power, from the smoke run): planting d = 0.3 into a random 60-gene set in mutant cells of the real
  data is detected by the pooled test in >= 80 % of 200 plants. If it fails, a null primary result is
  reported as "not readable", not as "not clone-tracking".

## 6. Primary result and its wording (fixed now)

Per program, Holm-corrected over the 3 MPs. The feature panel F1-F13 is Holm-corrected within itself, uses the same
wording rules, and is labelled exploratory: a panel hit cannot stand in for a primary result.

- "Travels with the clone": pooled d > 0, CI excludes 0, and at least 5 of the eligible patients with d > 0.
- "Present in wild-type cells to the same degree": pooled CI inside (-0.2, 0.2) and GATE_P passed.
  (0.2 = a fifth of the cell-to-cell SD, the conventional "small" effect; without a fixed bound a non-significant
  result could not be told apart from low power.)
- Otherwise "not resolved". A negative pooled d with CI excluding 0 is reported as such and needs discussion.

## 7. Secondary (descriptive, no gate)

- S0: new data that meets section 2 (per-cell genotype with coverage, mature cells retained, T-cell control passed)
  may join at any time (user 2026-11-13) as a dated amendment, with its patients added to the pooled estimate and
  the result before and after reported.

- S1: the same contrast in B and NK cells where the clone enters them (expected rare; reported per patient).
- S2: van Galen only, which has 5 same-protocol healthy marrows: ratio R = (mut - wt) / (mut - healthy) for each
  program, i.e. the share of the AML-versus-healthy elevation that is clone-intrinsic. Cross-sample, so
  handling is not cancelled; reported with that caveat.
- S3: per-patient table of d, CI, n_mut, n_wt, c_bar, T-cell mutant fraction, platform.
- S4: wild-type >= 1 reference read as the alternative label rule.

## 8. Smoke before the real run

Two worlds on the real data, as required after the G07 lesson: a null world (labels shuffled) that must run
the "not resolved" and "same degree" branches end to end, and a planted world (d = 0.3) that must run the
"travels with the clone" branch. Smoke outputs are written under `*_smoke` names and never next to real
covariates.

## 9. Known weaknesses (stated now)

- Seven patients on three platforms with three labelling methods; pooling is across patients only.
- A.13 has 32 and 548327 has 25 wild-type cells: most power comes from AML556, B.4 and 721214.
- Five of seven are NPM1-mutant; a "travels with the clone" result is mostly an NPM1 statement.
- The eligibility rule was chosen after the counts were seen.
- CloneTracer's model sets T cells as a separate prior class, so its T-cell control is weaker than Petti's.
- A null result cannot separate "environmental" from "not AML-related" except in the van Galen S2 arm.

## 2026-11-13 RESULT, Stage 0: GATE_S0 FAILED (2 of the 7 listed patients pass; >= 5 required)

Counts only; no expression was read. Compartment exactly as in section 2 (CD14 Mono, CD16 Mono, cDC1/2; projection
labels for Petti and van Galen, authors' ct labels for CloneTracer); driver locus only (NPM1 for Petti/van Galen).

| patient | Mono/DC mutant | Mono/DC wild type (primary rule) | note |
|---|---|---|---|
| A.13 | 1,986 | 26 | passes (T 0.09, authors' T labels) |
| B.4 | 1,660 | 995 | passes (T 0.035) |
| B.1 | 19 | 173 | fails by one mutant cell |
| AML556 | 37 | 2 (45 with >= 1 wt transcript) | fails |
| AML419A | 16 | 1 (33) | fails |
| 548327 | 2 | 4 (22) | fails |
| 721214 | 0 | 0 (2) | fails: almost no mature monocytes after QC |

Why the section 2 table was wrong: the feasibility counts behind it used a wider compartment (van Galen
"ProMono-like" cells and, for Petti, all myeloid-lineage cells), most of which are blasts. That error was in the
feasibility step, before this file was frozen; it is recorded here, not corrected in section 2.
Also found: the BoneMarrowMap T_NK bin is not a valid T-cell floor in Petti (782328: 44 of 55 genotyped
"T_NK" cells mutant versus 6 of 57 with the curated labels), consistent with blasts mis-binned as T_NK.
Per section 5 the line stops here for discussion. Counts: code/stage0_recount.R in the report folder.
