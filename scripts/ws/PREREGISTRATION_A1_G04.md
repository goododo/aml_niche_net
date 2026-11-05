# PRE-REGISTRATION A1-G04 — does the T-cell imprint direction hold in a second study whose two arms share one protocol (Chen2023)? With a secondary direction check in GSE185381 children

WP-A1 gap G04, options (a) and (b) of the gap assessment (2026-10-28, section 2.2). (a) asks the whole-T direction
question in Chen2023. (b) asks the same whole-T direction question, one-sided, in younger patients of GSE185381.
Because the 10 healthy donors helped select P370, (b) uses the A1-G02 out-of-fold gene lists instead of P370, so the
children and the healthy donors are both out of sample (2.4). The obligations handed on by A1-G01 section 7, A1-BATCH
12.5, the gap-assessment addendum (sections 1.3 and 1.7), the September pre-registration and the handoff are listed in
0.1, each with where it is met or why it is declined.

**Status: v1, APPROVED by the user 2026-11-04 (all section-9 items as recommended), after reviews 1 and 2 (section
10); tag `prereg-a1-g04-v1`. Nothing declared below was computed before the tag.** Written from a
read-only inventory of metadata, file names and code. No imprint-gene value, score or AML-vs-healthy expression
contrast of Chen2023 or of the children was read; 1.3 lists what was seen, including expression-derived QC covariates
of three children. Order (user decisions 2026-10-28 to 2026-11-04): A1-G02 (PASS) → A1-BATCH (PASS with
qualifications) → A1-G01 ("beyond T-subset composition, detected inside the main labels") → **this** → G03. Handoff
v2.0 §17 S15.

## 0. Question, scope, and what each outcome lets the paper say

S5's imprint procedure separates 15 adult AML from 10 healthy donors of GSE185381 out of sample (A1-G02 cross-fit: gap
1.378 = 1.98 out-of-fold healthy-SD; full-pipeline permutation p 0.008). T-subset composition accounts for about 5 % of
the cross-fitted gap (A1-G01: C 0.045), and a shift is detected inside the main T labels (a sign check). In GSE185381
the healthy marrow is a commercial cryopreserved product and the AML samples come from a hospital tissue bank; the
healthy donors are about 31.5 years younger; and the arms were partly captured on different dates. The within-study
checks could not separate disease from sample source and handling (A1-BATCH 12.5; A1-G01 section 0). **G04(a) asks
whether the direction holds in a second study whose two arms were sorted, cryopreserved and sequenced under one
published protocol. G04(b) asks whether it holds in younger patients of the same study.**

- **(a) primary: Chen2023.** 5 NPM1-mutated AML at diagnosis (39-59 y) vs 4 allogeneic-transplant donors (47-53 y).
  Every aspirate was cryopreserved, FACS-sorted into five fractions with one gating scheme and captured as two 10x 3' v3
  libraries per donor (one donor per library). T cells come from the lymphoid sort (pool B), which separated them from
  the blasts before capture. The ages are matched, so age cannot produce the direction here.
- **(b) secondary: GSE185381 children.** 18 patients of Lasry's paediatric cohort at diagnosis (14 under 18 y; 2.3)
  with at least 20 kept T cells (PAUMTZ excluded), against the 10 healthy donors scored out of fold. Only an AML-side
  result is interpreted; any other result gives no statement (A_b and its CI are still reported, 3.4).

What (a) changes compared with GSE185381, and what it does not:
- Shared by design: the FACS gating scheme, 10x 3' v3 chemistry, one donor per library, and every step of this
  project's processing. Both arms were cryopreserved; the freezing protocol, freezing site and storage time are not
  recorded.
- Still different between arms: the sample source (diagnostic aspirates vs aspirates from transplant donors). No
  collection site, date, capture session or sort batch is recorded on disk. Every library holds one donor, so arm may
  coincide with capture session (GSE185381's hashed pools at least gave within-capture pairs).
- A1-BATCH 12.5 called Chen2023 "one centre"; no local source supports that, so this document says "one study, one
  protocol". This also narrows A1-BATCH 12.5 and the handoff S15 note of 2026-11-01, which call G04(a) the test whose
  design separates disease from sample source. G04(a) replaces GSE185381's source contrast (commercial vs tissue bank)
  with another one (transplant donor vs diagnostic aspirate); it does not separate disease from sample source.
- Population: all 5 AML carry NPM1 mutations (4 of 5 also DNMT3A). A result speaks for NPM1-mutated AML only.
- Power (3.3): a pass is informative; a fail bounds the share that transfers (ρ_U, 3.1); a significant reversal is
  informative.

What each outcome lets the paper say (rules in section 7; pass level 6/126, see 9.2):
- **(a) "direction replicated"**: "In a second study whose arms were cryopreserved, FACS-sorted and sequenced under one
  published protocol and are age-matched (Chen2023: 5 NPM1-mutated AML at diagnosis vs 4 transplant donors; pool-B T
  cells), AML T cells score higher on the frozen imprint than donor T cells (exact p = k/126; ρ with its bounds)." This
  is evidence that the direction is not specific to GSE185381's procurement and processing (commercial vs tissue-bank
  marrow, its capture dates and co-capture pools, T cells captured together with blasts). It cannot exclude a difference
  between donor and diagnostic marrow that both studies share (for example blood admixture, A1-BATCH 7), pre-freeze
  handling and storage time, capture and sort sessions (all unrecorded in Chen2023), or the clinical correlates of
  diagnosis (handoff section 14). A1-BATCH 12.5's may-not list applies unchanged. Added where they apply: 3.2's "rests
  on one donor (name)", and 5.2's sentence on the processing-delay program.
- **(a) "not shown"**: "The direction was not shown in Chen2023 (exact p = k/126); a common shift larger than ρ_U of
  GSE185381's selection-corrected gap is excluded (one-sided 95.2 %, common-shift model)." Only if ρ_U ≥ 1: "5 vs 4
  cannot exclude an effect of the GSE185381 size." The direction stays "shown in one study".
- **(a) "reversed"**: "AML T cells score lower than donor T cells in Chen2023 (exact p = k/126)." The direction does not
  reproduce in Chen2023 (NPM1-mutated AML vs allogeneic-transplant donors, one sorting protocol); study-, genotype- and
  source-specific explanations cannot be told apart. Reported with the mirrored 3.2 qualifier and the depth-equalised
  result (4.6).
- **(a) "not decisive"**: a scope-P known answer could not be passed, a statistical bar or the coverage floor of
  section 6 failed (K6b(i) only for a would-be pass), or K3's swap rule applies. The results are reported as
  descriptive, and the direction stays "shown in one study".
- **(b) "AML-side"**: "Within GSE185381, T cells of 18 paediatric-cohort AML patients (14 under 18 y; 0.2-21 y, two
  without a recorded age) score higher on the cross-fitted imprint lists than out-of-fold adult healthy donors more
  often than not (A_b, 95 % CI; placement φ, 95 % CI)." The patients are younger than the donors, so an older-age
  effect cannot explain it. No child shares a pool with an adult AML, so it does not need the adult-AML lanes: 14 of the 18 were captured
  before 2020-07-29, like 8 of the 10 healthy donors, 12 of them in children-only pools. It does not separate disease
  from sample source (the comparators are the same commercial donors), from a shift shared by the children's lanes,
  from blast-derived ambient RNA (unsorted marrow, blasts 23-94 %), from T-subset composition (not controlled in (b)),
  or from paediatric T-cell biology. (b) is the same study and never counts toward GATE_A3.
- **(b) "not readable"**: no statement. A null is not evidence against the adult direction, because paediatric AML
  differs in drivers and in T-cell biology (1.4). A_b, φ and their CIs still go to G03.
- **(b) "not decisive"**: a scope-B known answer could not be passed (section 6). Descriptive only.

What G04 cannot do, under any outcome:
- In (a), separate disease from the remaining source difference, pre-freeze handling, storage time, or capture and
  sort session; Chen2023 has no handling or session records.
- In (a), exclude blast-derived ambient RNA. Pool B still carries some myeloid or progenitor-labelled cells, more in AML
  (2.7-8.0 % vs 0.4-1.5 % of pool-B cells; seen 2026-11-04, 1.3). Per-donor ambient indicators are reported (3.1), not
  tested; 4.7 is a decontX sensitivity.
- In (a), exclude T cells derived from a founder (pre-leukaemic) clone. 4 of the 5 Chen2023 AML carry DNMT3A
  mutations, and no per-cell genotype exists. In Petti 2019, founder variants were found in 2 of 171 T cells and later
  variants in 0 of 395, so this is possible but probably minor; the same applies to GSE185381.
- Separate disease from its clinical correlates at diagnosis (infection, medication before aspiration; handoff section
  14).
- Generalise (a) beyond NPM1-mutated AML, or estimate the residual's size. 5 vs 4 can support a direction and bound its
  size (3.1); the residual (5.1) is descriptive.
- Test T-specificity: B and NK cells (5.3) are descriptive, as in A1-G01 5.1.
- With (b), address source or handling: the children share GSE185381's design and its healthy donors, and 12 of the 18
  were captured on dates with no healthy donor. Nor can (b) address blast-derived ambient RNA or paediatric T-cell
  biology outside S9's 24-84 y.

### 0.1 Inherited obligations

| source | obligation | where met, or declined and why |
|---|---|---|
| A1-G01 7 | G04(a)'s primary is the whole-T direction; the residual is reported | 3.1; 5.1 (residuals in u-units under every outcome) |
| A1-G01 7 | G04(b) shares that primary and is read one-sided | same whole-T direction, scored with fold lists because P370 is in-sample for the comparators (2.4); 7(b) |
| A1-G01 7 | G04(c) (E-MTAB-16321) once G01 shows a residual | due; not bundled (decision 9.5) |
| A1-BATCH 12.5 | read the delay score in Chen2023 | 5.2, with a declared reading; 4.5 |
| addendum 1.7 | exact enumeration of the 126 labellings | 3.1 |
| addendum 1.7 | "a pass needs at most 2 rank inversions (p = 4/126)" | deviation: one statistic (mean difference); its pass level is decision 9.2 (6/126 recommended; 4/126 keeps the addendum's alpha) |
| addendum 1.7 | Ensembl remap to the 2020-A axis; pool-B T cells | 2.2 |
| addendum 1.7 | ws/22's unit = the library | ws/22's decontX call runs inside 41, one unit per pool-B library (4.7); the ws/22 script itself is not run (8) |
| addendum 1.7 | ws/26's roles mapped explicitly | moot: 41 builds its own pseudobulks and takes the arms from the xlsx (2.2) |
| addendum 1.7 | a leave-one-healthy-donor-out rule | 3.2 (a donor-influence rule; the p-based rule of v1 could fire only at p = 6/126, review 1) |
| addendum 1.7 | a recorded decision to promote Chen2023 | decision 9.1 |
| addendum 1.3 | declare the September null as prior knowledge | 1.1 |
| addendum 1.3 | the rungs that add the children or Chen2023 go inside G04's pre-registration | 5.7 (declared here, run in G03; decision 9.5) |
| addendum 1.3 | disclose that S5's hit rule differs from September's | G03 (the addendum lists it among G03's consequences) |
| addendum 1.3; order item 2 | the blind-safe G03 rungs run before G04 | moved after G04's tag by A1-BATCH (section 4, "Moved out"); their urgent question, lane clustering, was answered by A1-BATCH G2 (block-respecting p 0.017) |
| September pre-registration :41-44 | the null is not reframed as underpowered or as needing a different model | 1.1 (binding wording) |
| assessment 5, item 4 | G10(a) and G14(a) bundled into G04 | not bundled (decision 9.5) |
| handoff S15 (2026-11-01) | "G04(a) is the only test whose design separates disease from sample source" | narrowed in section 0 |

## 1. Already known before this draft (disclosed)

**1.1 The September pooled pseudobulk DE** (pre-registration a66cbd0, 2026-09-16; no tag). T_NK bin (projection-only,
NK included, both Chen2023 pools): 39 AML (5 Chen2023, 20 GSE185381 children, 14 GSE185381 adults) vs 14 healthy (4
Chen2023, 10 GSE185381); `~ dataset + arm` with duplicateCorrelation on library; hit = q < 0.05 and |logFC| ≥ 1. One
gene met the hit rule (AC007952.4), but a count of 1 is inside both Monte Carlo permutation nulls (p 0.102 and 0.320,
1,000 draws each), so the bin was reported as having no reportable hits. A shift of |logFC| 1.5 planted into every AML
sample was detected for 90.5 % of genes. Its pre-registration promised not to reframe this null as underpowered or as
needing a different model (`scripts/12_pseudobulk_de/PREREGISTRATION_pseudobulk_de.md`:41-44). G04 binds itself to that
promise. G04 tests a different pre-specified quantity: one sample-level score of a gene set frozen in another cohort.
Whatever its outcome, the September T_NK result stays reported as a null for |logFC| ≥ 1.5 per gene, and no G04 result
is described as showing that the September analysis was underpowered or used the wrong model. 5.4 is descriptive.
- Prediction on record (addendum 1.3, 2026-10-28). Applying the September recovery curve to S5's per-gene sizes, the
  fit should have called about 156 of the 370 genes if all 39 AML carried them (about 63 with S5's sizes shrunk by r,
  about 17 halved), and about 4 if only the 14 adults did (about 1 shrunk). These are rough upper estimates, because the
  curve comes from mid-tertile genes planted into every AML sample. It called 1, which is what the adults-only, shrunk
  row predicts. The null fits S5 only if the imprint is absent or much weaker in the children and the Chen2023 AML, or
  if S5's per-gene sizes are strongly inflated (beyond r = 1.43). **The first reading predicts that (a) and (b) fail.**
- The September gene universe was a symbol intersection. Five P370 symbols name different Ensembl genes on the Chen2023
  and GSE185381 axes, among them AC007952.4, so the single September T_NK hit pooled two loci (Chen2023
  ENSG00000266311, GSE185381 ENSG00000262202). G04 matches genes by Ensembl ID (2.2).
- The September per-gene tables are public and git-tracked; that they are unread is attested, not enforced. They stay
  closed until G04 has a result (section 8). Seen so far: AC007952.4's values (A1-BATCH section 1) and, incidentally on
  2026-11-04, its September q value (1.1e-4); no other September per-gene value.

**1.2 Other Chen2023 results on record** (none used this imprint): the Chen-only LMPP_GMP DE (5 vs 4; 304 hits, exact p
2/126), of which NBM3 alone carried 101 by leave-one-NBM-out; HSC_MPP 0 hits; the September pooled B_Plasma (28 hits)
and Mono_DC (385 hits) bins, which included Chen2023 (B_Plasma bears on 5.3); the opposite-sign transfer of the
CCC-era Discovery node-feature axis to Chen2023 (`scripts/DECISIONS_pending.md` P11); the leukemia_ml 7-vs-8
fraction-level CCC null and the panel screen; the FGW topology null; inferCNV malignant-call rates in presumed-normal
NK/CD8 effector-memory cells (AML 0.159 vs healthy 0.256) and per-donor malignant fractions. Chen et al. 2023's own
T-cell results were not reviewed for this draft.

**1.3 Seen while preparing this draft** (2026-11-04; inventory and review 1; notes and copied evidence in
`/FAST/gr10634/gaozy/reports/2026-11_a1_g04_prereg/`, outside the repo):
- Chen2023: kept T counts per BoneMarrowMap fine label for all 10 donors (all pools and pool B), so the per-arm T-subset
  composition is known to the drafting side. Pool-B naive share median 0.40 in AML vs 0.48 in donors, CD8 share 0.32 vs
  0.42, with overlapping ranges, so no confident forecast follows. Also seen: per-pool counts of every hierarchy bin (the
  pool-B myeloid/progenitor shares of section 0); kept NK, B and monocyte counts; QC thresholds and retention per
  library; per-library pre-QC median nCount; the leukemia_ml atlas lineage counts per library; per-donor inferCNV
  malignant fractions.
- Children: kept T counts per fine label (naive share median 0.295, range 0-0.93; per sample, healthy adults 0.200 and
  AML adults 0.152); hashtag tallies (identity); author labels (age band, sex, inflammation group); karyotypes, fusions,
  blast % and capture dates; three children's all-cell QC covariates (AML001, AML052, PAUMTZ: kept UMI and HB/granule
  UMI shares, which are expression-derived).
- Review 1 (metadata only): all 18,322 S4-kept T cells of GSE185381 match an author cell. The kept T cells that the
  authors call non-T are NK (median share: children 4.3 %, healthy 6.3 %, adult AML 10.0 %; AML028 16 of 20, Control0004
  116 of 250) and one MPP (AML2123); no other kept T cell is author-labelled myeloid or progenitor. The authors flag 278
  kept T cells as malignant (277 CNV+; AML2123's single CNV− cell is the MPP above), spread over T subtypes: AML073 210
  of 605, AML0693 47 of 336, AML022 9, AML001 7, AML3266 3, AML2123 1, PAUMTZ 1.
- Consequence: the composition readouts (5.1, and 3.4's per-child shares) are descriptive only.

Not computed (blind at this draft): any Chen2023 or child imprint, delay, or B/NK score; any Chen2023 or child expression
value of a P370 gene; the Chen2023 pool-B library sizes, UMI per cell by arm and sex calls; the September per-gene
tables beyond 1.1.

**1.4 The children** were in the September fit (20 of them). GATE_B0 used their T_NK IFNG/TNF/TGFB1 values, and
GATE_W1.1/T3c their t3c preflight ligand and receptor values. The CCC audit found blasts in their ws T_NK bins (for
example AML009, 92 %); the S4 kept-T set does not repeat this (1.3). The 79-gene ambient blacklist and the 32
ambient-prone genes come from 47 decontX units, 19 of which hold children (label-free gene lists). All 10 healthy
out-of-fold scores were read in A1-BATCH section 1 (Control4 −2.20 … Control0058 +0.17), so only the children's side of
(b) is blind; with these comparators, AML-side needs A_b of about 0.70-0.74 or more. The S4 sex check (2026-10-17)
showed expression-sex mismatches for PAUMTZ and PAXMIJ. Comparator facts: Control0058 (fold 18, the highest healthy
out-of-fold score) carries two hashtags (Healthy_0058 1,231 and Healthy_0059 897 kept T cells); Control4 (fold 23, the
lowest) is expression-female, coded male in A1-G02's fold lists, which (b) uses as stored. Lasry 2023 (the GSE185381
paper): CD4/CD8 T cells were depleted only in children, and inflammation remodelled the paediatric T compartment (Treg
and GZMK+ CD8 expansion), not the adult one.

**1.5 GSE185381 facts that bear on the readings.** The four adult AML under 60 y score 1.07 / 0.48 / 0.83 / 0.53 against
a healthy range of −2.27 to −0.43 (in-sample; S9 3.4). The imprint score follows the CD8 share in GSE227903 (Spearman
0.86, n 7) but not in GSE185381 (0.17 AML / 0.22 healthy). The direction is shared by B and NK cells of the same donors
(A1-BATCH; A1-G01 5.1) and is aligned with a processing-delay program (A1-BATCH D8).

**1.6 Status in the handoff.** Chen2023 is a direction "reference only" (`HANDOFF_within_sample_v2.0.md`:102, 160, 357,
445) and carries the CCC-era Discovery label (`02_study_split.csv`:2). The handoff's GATE_A3 direction measure is the
gene-level direction-concordance rate of imprint genes (HO:160); 9.6 replaces it with the sample-level exact test,
because gene-level concordance against random gene sets fails at 5 vs 4 (addendum 1.7), and 5.4 keeps the gene-level
measure as a descriptive readout. The children were planned as a separate second cohort without healthy comparators
(HO:247), paused on 2026-10-28 (HO:511); keeping (b) reverses that plan and restarts part of the paused item (9.3).

## 2. Inputs and frozen definitions (fixed before any expression is read)

### 2.1 Shared inputs

| item | file | md5 |
|---|---|---|
| imprint (P370 primary; 402 sensitivity) | `results/tables/ws/ws_s5_imprint_genes.csv` (compartment T, imprint_called) | 977f58b06da076cc39a561ab7445b69d |
| 32 ambient-prone genes (removed) | `results/tables/ws/ws_s8_t_imprint_ambient_prone.csv` | acb2b043049a4286e4052d04cea15feb |
| GSE185381 kept-cell pseudobulks | `/LARGE1/…/ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds` | 6c2fb685c49d4a23ef20a6842aa7aa67 |
| the 25 S5 samples and status; S5's B (23) and NK (21) designs | `results/tables/ws/ws_s5_scores.csv` | 58374fc60aabc53ba020c3149764fe8f |
| A1-G02 fold lists | `results/tables/ws/ws_a1_g02_fold_genes.csv` | ede3280052ab8dfcc547e60c341f334b |
| A1-G02 out-of-fold scores | `results/tables/ws/ws_a1_g02_oof_scores.csv` | b1da65a5ddd4d050a1b0849a9aa1206d |
| A1-G02 leave-pair-out scores (K9) | `results/tables/ws/ws_a1_g02_lpo_pairs.csv` | 166623a1461ae59f0c0095e226788da4 |
| r = 1.430 (gap_in / gap_oof) | `results/tables/ws/ws_s9_yardstick.csv` | dcbea817a8eb6c68575486e3576dfce9 |
| S9 result, D_res_J at full precision (K4) | `results/tables/ws/ws_s9_direction.csv` | a95d8caf3607b8a0b05083346b7301b6 |
| K4 inputs (Oetjen) | `/LARGE1/…/ws_pseudobulk/GSE120221/ws_pb_GSE120221.rds`; `results/tables/ws/ws_s4_pseudobulk_meta_GSE120221.csv`; `/LARGE1/…/00_raw/metadata/GSE120221/donor_ages_GSE120221_v1.0.csv`; `/LARGE1/…/00_raw/metadata/GSE120221/meta_GSE120221_v1.0.csv` | b4ad90bd5c587d48cc66e02a537c5079; 11a282f9586be7f5d4090d93bd047f7d; 178ced1c73e531b67172ed413b009c07; 222ff9a6e85b08e551cd1ad95395c076 |
| ambient blacklist (n_units ≥ 5: 79 genes) | `results/tables/ws/decontx/GSE185381/ws_s1_ambient_blacklist.csv` | acebed4dbeed2ef6ecd5a60791e28fb9 |
| GENCODE gene order (the chrY list inside `black_for`) | `/LARGE1/…/reference/gencode_GRCh38_gene_order.txt` | c0950d592084a8b713de34ab30a7ddf8 |
| processing-delay set (295 up / 497 down) | `/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv` | 5ee444a9284ae52b7afab1c8aeeced2c |
| S4 label map | `scripts/ws/s4_cell_groups.tsv` | f5a5a0283b243e5efaf17b1ddfdda5f1 |

`/LARGE1/…` is `/LARGE1/gr10634/gaozy/aml_niche_net`.

### 2.2 (a) Chen2023

- **Objects:** `/LARGE1/…/02_seurat_objects/04_annotated/Chen2023/<donor>.rds`, md5: AML103
  5c45be9c3e98789944c3592afd0c2ad8, AML162 44866549e35fb6bbf50a8d1898859574, AML24 d4fe077315c91cd6ca240c3077e0cf88,
  AML304 ff64ce4ef252a71f3d1a2dbebae80b74, AML320 e31fd77fbab7b33b3ef87c3addbfa7ad, AML469
  4b08d3dfb5aa1aece3043fa6c3e823c9, NBM1 29bcb07cc25e58869926ef7f0c4324b5, NBM2 5364a9c13432d342d36805c400e8ab93, NBM3
  3876099a82edb40faeca04b7d7f95636, NBM4 de8c3cc44f518bb59cc216ac40b4f085. Metadata
  `/LARGE1/…/00_raw/metadata/Chen2023/meta_Chen2023_v2.xlsx` (c2407360893ccf5261dff4e2144cd996). Deposit
  `/LARGE1/…/00_raw/public/Zenodo/10.17632_gwjh3w6ztm.2/filtered_feature_bc_matrix.7z_*/filtered_feature_bc_matrix/features.tsv.gz`:
  the md5 of the decompressed content (`zcat | md5sum`) is 16952af81ef735f42952b64c8402f69a in all 20 libraries (the
  .gz files themselves differ). Remap `results/tables/00_ingest/GSE120221_gene_remap.csv`
  (d8d7c37eeb8cb85f5bd4b9d53d1f0acb).
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
| AML469 | AML (excluded) | 50 | M | NPM1, FLT3-ITD, FLT3-TKD 0.025, STAG2, NRAS | 2 |

- **Cells:** pool B (meta.data `Pool == "poolB_stroma_lymphoid"`; library suffix Niche_Immune in four spellings) and the
  S4 kept rule (bmm_fine maps to group4 T in `s4_cell_groups.tsv`, `anno_source == "concordant"`, `high_error != 1`).
  T cells should not pass pool A's CD34+ and CD117+/CD33+ gates; the 1-339 kept T cells per donor in pool A are sort
  impurities, doublets or mislabels (inferred), captured with the blasts. They are excluded (sensitivity 4.3 adds them).
- **Unit and eligibility:** the donor (one pool-B library each), with at least 20 pool-B kept T cells. Eligible: 5 AML
  and 4 donors (table). **AML469 is excluded:** its pool-B library lost 526 of 552 barcodes to a donor-level haemoglobin
  QC threshold of 0, leaving 2 kept T cells. It is not rescued by another QC or by the leukemia_ml atlas cells (a forking
  path closed now; decision 9.4).
- **Gene axis:** row i of each object ↔ row i of `features.tsv` ↔ row i of `GSE120221_gene_remap.csv` (K1). Rows with rule
  `mapped_by_ensembl` (32,351) take their 2020-A name; every other row (1,343) is dropped, so no gene is matched by symbol
  alone. The primary set on this axis has **358 genes (154 up / 204 down)**. Missing: AC016831.7, HNRNPA1P48,
  AL021155.5, CCDC18-AS1, AC253572.2, AC239799.2, SATB1-AS1, AL499604.1, PRKCQ-AS1, RPARP-AS1, LINC01278 (no 2020-A ID on
  the Chen2023 axis), and POLR2J3 (Chen2023 ENSG00000168255 vs 2020-A ENSG00000285437; dropped by the ID rule). The 402
  set has 388. The scored set is P370 ∩ U_J (2.4), at most 358.
- **Matrices** (`41 --build`, per donor, all on the remapped axis): pool-B kept T (whole), with its nine zero-filled
  BoneMarrowMap-label count matrices and the label × donor kept-cell numbers (5.1); pool-A kept T (4.3); pool-B kept B
  and NK (5.3); and the decontX-corrected, rounded pool-B kept T (4.7).

### 2.3 (b) GSE185381 children

- **Samples:** the 19 paediatric samples with at least 20 kept T cells (S4), minus PAUMTZ: AML001, AML003, AML005, AML006,
  AML009, AML012, AML022, AML025, AML026, AML028, AML038, AML043, AML048, AML051, AML052, AML055, AML056, AML073 (18).
  **PAUMTZ is excluded:** 17 of its 31 kept T cells carry the PAUMTZ_EoT hashtag as raw top hashtag; its expression sex
  contradicts its karyotype (46,XY); and reading it would read post-treatment cells, which this document closes for any
  later remission or relapse analysis.
- **Identity defaults:** Control0058 is used as stored in S5 (one donor). AML003 has no Table S1 record; its age band
  (infant) and sex come from author metadata, and its St Jude source is inferred from its pool (P08) and ID. (b) uses the A1-G02 fold lists as stored, with Control4
  coded male.
- **Matrices:** the `PB$raw$T` columns of these 18 (author labels, demultiplexed by hashtag plus genotype), subset
  immediately after loading. Never read: the UNDEMUX cells of the ingest object, any GEO lane matrix beyond the hashtag
  rows already tallied (1.3), and the PAUMTZ and PAXMIJ columns of any matrix.
- **Metadata:** `/LARGE1/…/00_raw/metadata/GSE185381/meta_GSE185381_v2.3.csv` (cf0c1ed37db4dabbac440ce26494f586; ages);
  the author per-cell metadata `/LARGE1/…/00_raw/public/GEO/GSE185381_RAW/GSM*_metadata.csv.gz` (age band,
  inflammation group, author cell identity and malignancy flag; 63 files). Their md5 is checked by running, inside
  `GSE185381_RAW`, `md5sum GSM*_metadata.csv.gz | sort -k2 | md5sum` (base names, md5sum's two-space format), which
  gives c9bf8d1d53653aad8d65d1da7f9da4e0. Also `ws_s4_pseudobulk_meta_GSE185381.csv`
  (8533d147444e3b70983693dfd77c5beb; named columns only, K8), and the children's rows of `PB$cells` (cell, sample, fine,
  group4 and kept only; no counts) for 3.4's per-child shares.
- **Flags** (sensitivity 4.10; from metadata and hashtags only):
  - age at least 18 y, or no Table S1 age and no author band below 18 (AML022 19.75, AML038 21.08, AML043 19.58, AML073
    none);
  - no Table S1 record (AML003);
  - an age conflict (AML052: Table S1 4.5 y, author band "teens");
  - therapy-related AML (AML043);
  - more than 5 % of kept T cells whose raw top hashtag belongs to another deposited sample (AML056 122/250, of which
    119 carry AML043's hashtag, so some may be AML043's cells; AML038 10/35; AML028 2/20; AML012 3/58);
  - P19 identity unverified (AML026, AML051): 619 of AML051's 653 kept T cells carry a hashtag with no deposited sample,
    and 32 of AML026's 41 carry the hashtag named AML051;
  - fewer than 50 kept T cells (AML028 20, AML043 29, AML038 35, AML026 41);
  - more than half of the kept T cells author-labelled NK (AML028 16/20);
  - more than 20 % of the kept T cells author-flagged malignant (AML073 210/605).
  - Flag-free: AML001, AML005, AML006, AML009, AML025, AML048, AML055 (7).

### 2.4 Scores (all label-free; S5 signs)

- **J space ((a) primary; as S9 2.4).** `black_for()` is applied on each axis. U_J = names with CPM > 1 in at least half
  of the 25 GSE185381 samples and in at least half of the 9 Chen2023 donors. TMM within each cohort on U_J; log-CPM with
  prior count 1; z with the 25 samples' per-gene mean and SD. Score = mean z(up) − mean z(down) over P370 ∩ U_J; each
  half needs at least 5 genes, and the coverage floor of K4b applies.
  - Why J decides: ρ and its bounds need GSE185381's units (S9 2.4), and in J the category and the bounds agree (ρ_L > 0
    exactly when the test passes). In simulation W has the same power and is less sensitive to depth (review 1), so it
    is sensitivity 4.2; decision 9.7.
  - Depth: shallower pseudobulks bias J scores downward, toward "reversed". At Chen-like cell numbers this is about
    −0.02 J at 2,000 UMI per kept T cell and about −0.05 at 1,000; W is unaffected (review 1, GSE185381 healthy donors
    thinned). K6c reports the implied bias at Chen2023's realised depths, and 4.6 equalises depth.
- **W space (4.2).** `logzZ()` on the 9 Chen2023 donors alone (blacklist, CPM > 1 in at least half, TMM, prior 1), with z
  over the 9.
- **Out-of-fold space ((b)).** L = `logzL()` log-CPM of the 25 (as A1-G02 R2). The children are placed with `make_L27()`
  (`37_a1_batch.R`:108-118): each child is TMM-normalised against the 25's reference column, with edgeR's prior scaled
  to the 25's mean library size, so L on the 25 is unchanged. For healthy fold h (folds 16-25 of
  `ws_a1_g02_oof_scores.csv`), z = (L_x − mean of L over the 24 others) / SD over the 24 others, and the score uses fold
  h's lists. The children enter no filter, TMM or z parameter, so the stored healthy scores are reproduced exactly (K7).
  - (b) keeps (a)'s whole-T direction but uses the A1-G02 fold lists (75-942 genes in the healthy folds; 29-1,075 over
    all 25), because P370 was selected on all 25 and is in-sample for every healthy comparator (in-sample healthy SD
    0.513 vs out-of-fold 0.696).
  - Depth: shallow columns are expanded and noisier under `make_L27`'s fixed prior. In adults cut to 20 / 41 / 71 / 118
    kept T cells and scored in their own fold, the shift is α + β × score with β about 0.26-0.29 / 0.17-0.20 / 0.12-0.13
    / 0.07-0.09 and noise SD about 0.37 / 0.23-0.25 / 0.15-0.16 / 0.11 (review 1, adults only). It is a UMI-depth effect.
    Under H0 it is conservative; it inflates the per-child scores and φ of shallow AML-like children. 5.6 measures it
    per child.

## 3. Primary estimate for (a) and the declared readout for (b)

### 3.1 (a) D, its exact p, and the bounds on ρ

D = mean J score of the 5 AML − mean J score of the 4 donors, over P370 ∩ U_J. Exact one-sided p = (number of the 126
labellings of the 9 donors into 5 + 4 whose D is at least the observed D, within 1e-12) / 126; the observed labelling
is one of the 126. The null is that the 9 label-free scores are exchangeable between the arms (K6a checks that the
scores do not depend on the labels). Arm-correlated technical factors are part of what a pass or a reversal can
reflect: pool-B kept T cells number 211-1,445 in AML vs 1,094-2,795 in donors, and the pool-B myeloid share differs
(section 0).

**ρ = D / (D_res_J / r)**, the share of GSE185381's selection-corrected gap, where D_res_J is the AML − healthy mean of
the 25 GSE185381 J scores on the same U_J and r = 1.430. **Bounds**, under a common-shift model (every AML score shifted
by δ): for each labelling S other than the observed one, t_S = (D_obs − D_S) / (1 − c_S), with c_S = k_S/5 − (5 − k_S)/4
and k_S the number of AML donors among S's five. δ_L is the 6th smallest t_S and δ_U the 6th largest; ρ_L and ρ_U are
δ_L and δ_U divided by D_res_J / r, each one-sided 95.2 %. ρ_L > 0 exactly when p ≤ 6/126 (K5). Reported under every
outcome.

Reported with D:
- the per-donor scores, each next to the donor's pool-B kept-T count, pseudobulk UMI total and number of detected genes;
- per-donor ambient indicators, no test: the pool-B myeloid/progenitor-labelled share (metadata, already seen), and the
  LYZ + S100A8 + S100A9 CPM of the T pseudobulk, read from raw counts during the run only (LYZ is blacklisted and so
  outside the score; S100A8 and S100A9 are in neither P370 nor the blacklist);
- the rank-sum count U (pairs with the AML donor higher, of 20) and its exact p;
- D over the pooled within-arm SD (7 df; descriptive);
- ρ, ρ_L and ρ_U, and K6c's implied depth bias beside D.

Fixed caveats for ρ: ρ is descriptive. It is attenuated by shallower AML pseudobulks (by about 0.05 J, about 0.04 of
ρ, at 1,000 UMI per cell; more if shallower) and by protocol differences. The GSE185381 gap it divides by includes an age share that
Chen2023 lacks (S9: f_age 0.102, CI −0.085 to 0.463), so ρ < 1 is expected even without a source effect.

### 3.2 (a) Does the result rest on one donor?

One healthy donor can carry a result in this dataset: NBM3 alone carried 101 of the 304 September LMPP_GMP hits (1.2).
The 9 scores of 3.1 are reused; nothing is re-normalised (S9 precedent). D_(−i) is D recomputed without donor i, for
each of the 9 donors. **Qualifier "rests on one donor (name)"**: after a pass, some single drop leaves D_(−i) below 2/3
of the observed D; the donor whose drop gives the smallest D_(−i) is named. "Reversed" gets the same rule with D's sign
mirrored.

Generic rates (synthetic normal scores, 5 vs 4, pass at 6/126; review 1, at least 100,000 runs per scenario): the
qualifier flags about 0.17 of chance passes and 0.11 / 0.08 / 0.05 of true passes at a common shift of 1 / 1.5 / 2
within-arm SD. It flags about 0.51 of passes made by one healthy donor at −3 SD, and 0.40 / 0.83 of passes made by one
AML donor at +3 / +5 SD with the other four null.

Also reported: the exact p of each leave-one-out set (5 vs 3: 56 labellings; 4 vs 4: 70). After a pass with full count
k, every healthy-drop p is at most k/56 and every AML-drop p at most k/70 by construction, so these p values are
reported only. The v1 rule ("any 5 vs 3 p above 0.10") could fire only when p was exactly 6/126 and is withdrawn
(review 1, A1).

### 3.3 (a) Power (synthetic; no data)

Generic, normal scores with equal SD in both arms. The exact mean-difference test at p ≤ 6/126 passes with probability
0.05 / 0.37 / 0.63 / 0.84 at a true shift of 0 / 1 / 1.5 / 2 within-arm SD; at p ≤ 4/126, 0.03 / 0.29 / 0.53 / 0.76. The
rank-sum with at most 2 inversions (also alpha 4/126) gives 0.03 / 0.28 / 0.52 / 0.74. The mean difference matched or
beat the rank-sum at equal alpha in every scenario tried (heavy tails, unequal SD, one wrong-way donor), and one extreme
donor inflated both equally. An AML SD of 1.5-3 times the healthy SD kept the type I error at 0.045-0.048 (review 1).

On the ρ scale (effect = ρ × D_res_J / r = ρ × 1.376 J units; pass at 6/126; review 1):

| within-arm SD (J units) | ρ 0.25 | ρ 0.5 | ρ 0.75 | ρ 1 |
|---|---|---|---|---|
| 0.30 in both arms (S9's external anchor: Oetjen's young donors had 0.586 × GSE185381's healthy SD) | 0.45 | 0.91 | 1.00 | 1.00 |
| 0.30 healthy, 0.60 AML | 0.23 | 0.60 | 0.88 | 0.98 |
| 0.51 in both (GSE185381 in-sample healthy) | 0.22 | 0.54 | 0.84 | 0.97 |
| 0.70 in both (GSE185381 out-of-fold healthy) | 0.16 | 0.36 | 0.62 | 0.83 |
| 0.70 healthy, 0.82 AML (GSE185381 out-of-fold) | 0.14 | 0.32 | 0.55 | 0.77 |

With real donor biology (GSE185381 healthy pseudobulks at their own depth, S5 sizes planted), ρ of about 0.5 gave a
pass rate of 0.51; with the AML arm thinned to 0.35 of its depth, 0.43. Chen2023's SD is unknown before the run, and the September null (1.1) suggests a smaller
transfer. Hence a pass is informative, and a fail bounds the share that transfers: if nothing transfers, ρ_U excludes
ρ = 1 in 1.00 / 0.98 / 0.84 / 0.57 of runs at a within-arm SD of 0.30 / 0.51 / 0.70 / 1.0.

### 3.4 (b) A_b, the within-fold AUC

For each healthy fold h (10 folds), a_h = the share of the 18 children whose fold-h score exceeds healthy donor h's
out-of-fold score (ties count 0.5); **A_b = mean of the 10 a_h**. Both sides are out of sample for fold h's lists.

**CI**: 2,000 two-way bootstrap draws (seed SEED + 17). Each draw takes 18 children with replacement (each keeps its 10
fold scores) and 10 healthy folds with replacement (each keeps its comparator and its column of 18 child scores), and
computes A_b* on that 18 × 10 grid; percentile interval (R quantile type 7). The same procedure and seed serve 4.8-4.10
and φ's CI.

**φ** = (the children's mean score over all 25 folds − the healthy out-of-fold mean) / (the adult AML out-of-fold mean −
the healthy out-of-fold mean). The children are scored in all 25 folds because the two means of the denominator come
from all 25 folds' out-of-fold scores. φ's CI resamples the children and the 10 healthy folds as above, conditional on
the 15 adult AML out-of-fold scores.

Reported under every outcome, without a category:
- A_b and φ with their two-way 95 % CIs; φ also over the 8 children with at least 118 kept T cells;
- each child's mean score over the 10 healthy folds, next to its kept-T count, T library size, kept / T-mapped
  fraction, naive and CD8 shares, author-NK share and author-malignant share (all but the score are metadata or depth);
- the same-pool differences (AML022 and AML052 vs Control0082 in fold 19; AML009 vs Control0005 in fold 17), without a
  test;
- whether 180 · A_b > 131, the exact Wilcoxon 2.5 % point (seed-free).

### 3.5 (b) Operating characteristics (synthetic, no child data; review 1)

- Under H0 the rule calls AML-side in about 2-4 % of runs at full depth (1.7 % at the children's depths; up to 3.7 %
  with strong pool clustering).
- Power is 0.17-0.24 / 0.62-0.68 / 0.94-0.95 / 1.00 when the children carry 0.25 / 0.5 / 0.75 / 1 of the adult out-of-fold
  gap, and 0.24 if only the 7 inflammation-high children carry it.
- The CI covers the sampling of children and healthy donors only. Any shift the children share is read as AML-side
  (St Jude/COG source, children-only lanes, blast ambient RNA, paediatric T biology). A shared shift of 0.5 / 1
  healthy-SD is called AML-side in 18 / 63 % of runs; this is the same power curve, since a shift of d healthy-SD equals
  k = d / 1.98.
- For scale, healthy cohorts of other studies sit −0.2 to +1.2 out-of-fold healthy-SD from GSE185381's healthy donors
  (S9, J space; other labs and chemistries; an upper reference for a within-study shift).
- Real-data controls (adults only, nested refits):
  - adult AML placed through the child path at the children's depths give A_b 0.92 and are AML-side in 50 of 50 runs;
  - healthy pseudo-children give A_b 0.48 (0 of 50);
  - 18 pseudo-children cut from one held-out healthy donor are called AML-side for the three highest-scoring healthy
    donors (Control0004, Control1, Control0058).

## 4. Sensitivities (each reports its statistic with an exact p or a CI, except 4.11; none changes a category)

(a) uses 3.1's statistic and pass level; (b) uses 3.4's bootstrap.
- **4.1** (a) the 402 set (402 ∩ U_J; 388 on the Chen2023 axis).
- **4.2** (a) W space.
- **4.3** (a) all pools: pool-A kept T cells added.
- **4.4** (a) the up half and the down half, scored separately.
- **4.5** (a) P370 without the genes of the Savage delay set (295 up / 497 down); the analogue of A1-G01 4.2b.
- **4.6** (a) depth-equalised: all 9 pool-B T pseudobulks are binomially thinned to the smallest donor's UMI total, 20
  draws (seed SEED + 20). In each draw, U_J, TMM and D_res_J are recomputed (and W's filter and TMM). The scores are
  averaged per donor over the draws; D and the exact p come from the averaged scores, and ρ = D / (mean D_res_J / r).
  D, exact p and ρ in J; D and exact p in W.
- **4.7** (a) decontX: ws/22's call (`celda::decontX`, z = NULL, seed SEED) with one unit per pool-B library of the 9
  eligible donors, on all cells of that library in the annotated object (the QC-passed cells, as in ws/22). The decontX
  counts of each donor's pool-B kept T are summed and then rounded (as A1-BATCH D4, `37_a1_batch.R`:354). D, exact p and
  ρ in J.
- **4.8** (b) children under 18 y: Table S1 age below 18, or no age and an author band below 18 (14: AML001, AML003,
  AML005, AML006, AML009, AML012, AML025, AML026, AML028, AML048, AML051, AML052, AML055, AML056).
- **4.9** (b) at least 50 kept T cells (14); this is also the depth sensitivity.
- **4.10** (b) flag-free children (7; 2.3).
- **4.11** (b) comparators outside the T-ALL co-capture pools P02/P03: the folds of Control0004, Control0005,
  Control0058, Control0082 and Control4003. A_b over these 5 folds, with no CI (a percentile bootstrap over 5 folds is
  unreliable). It is the harder test, because these donors score higher (mean −0.594 vs −0.915 for Control1-5). Also
  A_b over the 7 folds without Control1, Control3 and Control4 (expression-female, pool P02), with no CI.

## 5. Further declared readouts (descriptive; no category)

- **5.1 Composition and the residual** (A1-G01 handover).
  - (i) Naive share ((CD4 Naive + CD8 Naive) / kept T) and CD8 share, per donor and arm (pool B).
  - (ii) The counterfactual of A1-G01 3.1 inside Chen2023: `est_C` (`39_a1_g01.R`:151-195) with Chen2023's pool-B label
    matrices, leave-one-out arm profiles (fallback below FB_CELLS = 50 cells), the P370 list for every donor, and
    u-units in `logzL()` over the 9 (`black_for()` on the Chen2023 axis). Reported in every case: gapT, gapC_H, gapC_A,
    and the three residuals gapT − gapC (each order and their mean), in u-units. These are not comparable with A1-G01's
    handover figure, (1 − C) × 1.98 ≈ 1.89, which is in out-of-fold healthy-SD units. C_H, C_A and C only under "direction replicated", and only when
    gapT > 0 (`est_C` sets them to NA otherwise).
  - Expected fallbacks (metadata): order H, CD4 Regulatory in all 4 donor folds (leave-one-out healthy pools 45 / 41 /
    38 / 38 cells); order A, CD4 Regulatory for AML162 (42) and CD8 Tissue Resident Memory for AML304 (47). AML304's
    CD4 Regulatory pool is exactly 50, so it does not fall back.
  - Composition was seen before this draft (1.3), and the counterfactual's references depend on the labels, so there is
    no test.
- **5.2 Processing-delay score** (A1-BATCH 12.5). The D8 score (mean z of delay-up minus delay-down genes of the Savage
  set, F402 genes removed) in the same J space. D_delay with its one-sided exact p in GSE185381's direction (AML higher;
  there the delay score was 1.12 z higher in AML) and its mirror p, each over the 126 labellings. Spearman of the
  imprint and delay scores over the 9 donors after centring each arm (descriptive; the within-arm Spearman of v1 is
  dropped, because at n = 4 its smallest one-sided p is 1/24).
  - Reading, fixed now: if D_delay > 0 with exact p ≤ 6/126, the imprint-delay alignment persists under one
    post-collection protocol, and the (a) "replicated" sentence adds: "the processing-delay program also differs
    between the arms in Chen2023, so pre-freeze handling is not excluded there either." Otherwise no sentence. Either
    way, it cannot separate handling from biology (A1-BATCH D8).
- **5.3 B and NK cells of the same donors** (is the direction shared, as in GSE185381?). J space as in 2.4, with each
  compartment's own GSE185381 reference: the 23 (B) and 21 (NK) samples of S5's B and NK designs (`ws_s5_scores.csv`)
  and their `PB$raw$B` / `PB$raw$NK` columns give U_J and the z parameters. P370 with S5 signs and the 5-gene rule; the
  coverage |P370 ∩ U_J| is reported (T-specific genes drop out). Pool-B kept B: 5 vs 4 (126 labellings); NK: 4 vs 4 (70
  labellings; AML320 has 18 NK cells, below 20). Reported: D and its one-sided exact p (AML higher), and D_B / D_T and
  D_NK / D_T, with D_T computed on the same donors. Reading: a positive D repeats GSE185381's finding that the
  direction is not T-specific. Monocytes are not scored: they sit in pool A with the blasts.
- **5.4 Gene-level sign concordance** (the handoff's GATE_A3 measure, HO:160, kept as a descriptive readout). The share
  of the genes of P370 ∩ U_J (n reported) whose AML − donor difference in mean J log-CPM has S5's sign; p = (labellings
  with a share at least the observed one) / 126. Random gene sets are not used as a null (they pass 27-40 % of null runs
  at 5 vs 4; addendum 1.7).
- **5.5 (b) strata.** Each child's mean score over the 10 healthy folds, summarised by Lasry inflammation group (7 high,
  11 low) and by author age band: the median per group, no test.
- **5.6 (b) depth readout** (adults only, after the tag). For each child, every adult whose kept-T library size is at
  least the child's is binomially thinned to the child's kept-T library size (20 draws, seed SEED + 26), placed through
  `make_L27` and scored in its own fold. α, β and the noise SD (shift = α + β × full-depth score) are reported per
  child, with the number of adults used; a child with fewer than 5 such adults gets no calibration. Thinning reproduces
  the depth effect that cell subsampling shows (review 1), so no per-cell adult object is read. Also reported, per
  child and fold: the share of the fold's genes with zero counts.
- **5.7 G03 extension rungs** (declared here, run in G03 after G04's result; addendum 1.3).
  - Rung "+children": G03's final blind-safe rung (the September specification on the 24 S5 samples of the September
    fit: 14 adult AML, 10 healthy) plus the 18 children of (b). These are the September fit's children without PAUMTZ
    and PAXMIJ, which stay closed (2.3); nothing else changes.
  - Rung "+Chen2023": that, plus the 9 Chen2023 donors. It is the September T_NK fit without PAUMTZ and PAXMIJ (37 AML
    vs 14 healthy instead of 39 vs 14), so its count is not known; the September fit itself (1 hit) is the reference.
  - Specification: the September pre-registration's design, gene universe, hit rule and permutation procedure
    (`PREREGISTRATION_pseudobulk_de.md`).
  - Readout per rung: the hit count under the rung's rule with its permutation p, and, among the P370 genes in the
    universe, the number of hits and the share with S5's sign.
  - G03 states its reading rules before it runs. The September per-gene tables are opened only there.

## 6. Known answers (run first)

Each known answer has a scope: P gates `--primary`, `--sens` and `--readouts` for (a); B gates `--children`. The checks
table has three states: PASS, FAIL (deterministic) and FAIL (statistical, final), as in A1-G01 section 6.
- Deterministic checks must pass exactly. If one fails, the cause is a code or input error: it is fixed, the fix is
  recorded in a dated note in this file before any readout of that scope is printed, and all known answers rerun.
- The statistical bars of K6b and the coverage floor of K4b run once, with the declared seeds, after the deterministic
  checks pass. A failure of K4b or K6b(ii) makes (a) "not decisive". A failure of K6b(i), which concerns only the
  "replicated" tail, makes a would-be "direction replicated" "not decisive" and leaves the other categories readable.
  Nothing changes in methods, thresholds or seeds. The modes of the scope still run, their outputs are labelled "not
  decisive, descriptive" where the failure applies, and the GATES block names the failed bar. A bar failure traced to a
  demonstrable code error (shown on synthetic data) is handled as a deterministic failure.
- A mode runs only when every deterministic check of its scope passes (non-smoke checks table). Rows marked "reported"
  never block, except K3's swap rule (below; 7(a)1).

| id | check | expected | scope |
|---|---|---|---|
| K0 | md5 of every declared input (2.1-2.3); for `features.tsv`, the md5 of the decompressed content in all 20 libraries | as declared | P, B |
| K1 | `features.tsv` (ID, symbol) equals the remap's (ens, sym_cr2) in all 33,694 rows, in order; every object's row names equal `make.unique(gsub("_", "-", symbol))`; P370 on the remapped axis | identical; 358 (154 / 204); 402 set 388 | P |
| K2 | the donor table equals the xlsx (arm, age, sex) for 10/10; `Pool` agrees with the library name for every kept cell; pool-B kept T per donor equals 2.2; the column sums of the saved (remapped) pseudobulks equal their cells' summed counts over the 32,351 mapped rows; the nine label matrices sum to whole-T | exact (max \|diff\| 0) | P |
| K2b | `est_C`'s reference profiles are bit-identical when donor i's label counts are set to zero (both orders) | identical | P (5.1 only) |
| K3 | expression sex: the ws/26 thresholds (`26_s4_pseudobulk.R`:91) applied to the summed XIST and RPS4Y1 UMIs of each donor's pool-B kept T, B and NK cells, and separately to its pool-A kept cells | see the rule below | reported; the swap rule gates P |
| K3b | children's identity: the same thresholds on each child's kept-T XIST and RPS4Y1 UMIs (`PB$raw$T` rows), against the metadata sex; the UMIs of the other-sex gene are printed for AML056, AML038, AML026 and AML051, where hashtags suggest mixing (2.3) | a mismatch or an "unclear" call gets a dated note before (b) is read, and the child joins the 2.3 flags | reported |
| K4 | the J-space code reproduces S9: D_res_J for GSE120221 from S9's inputs and donor pooling (`load_cohort`); gap_in on the 25 in S5 space | 1.96770393747843 within 1e-6; 1.97138 within 1e-5 | P |
| K4b | coverage: \|P370 ∩ U_J\| / 358, and the genes per half | ≥ 0.80 (S9's floor) and ≥ 5 per half; else (a) is "not decisive" | P (bar) |
| K5 | enumeration: 126 distinct labellings, the observed one once; the rank-sum exact p equals `pwilcox`; the closed-form bounds of 3.1 equal bisection; ρ_L > 0 exactly when p ≤ 6/126 | exact | P |
| K6a | order invariance: the J and W scores of the real pool-B matrix, recomputed with the 9 columns in 10 random orders (seed SEED + 28), equal the re-ordered scores (only max \|diff\| is printed). Plant signs, on synthetic data only: S5 log2FC / r planted into the 5 synthetic AML gives D > 0, and into the 4 synthetic donors D < 0 | < 1e-12; the signs as stated | P |
| K6b | machinery (synthetic): NB size 10 around the pooled kept-T profile of the 10 GSE185381 healthy adults, on the 32,351 Chen2023 mapped genes; library sizes = the `--build` pool-B column sums, each kept with its donor; U_J from the 9 synthetic columns and the 25; 1,000 data sets (seed SEED + 18); (i) no plant; (ii) S5 log2FC planted on P370 in the 5 synthetic AML; pass = 3.1's rule | (i) "replicated" rate ≤ 0.08; the synthetic donors have no between-donor variance, so the depth offset dominates the other tail, which is reported only; (ii) pass rate ≥ 0.90 | P (bar) |
| K6c | depth calibration (no bar): the 3 deepest GSE185381 healthy T pseudobulks, thinned to 3.0M UMI and then to 2.0 / 1.0 / 0.5 / 0.25 / 0.1M (20 draws each, seed SEED + 25), each scored in a J space whose external cohort is the probe plus the other 9 healthy donors unchanged (U_J and TMM recomputed at each step), give a curve of J-score shift against UMI total (zero shift at or above 3.0M; linear in log UMI between grid points; the 0.1M value below it). At Chen2023's realised per-donor pool-B kept-T UMI totals the curve gives the implied D bias: `--known` computes it in J units, and `--primary` prints it beside D, also in pooled within-arm SD units | reported | reported |
| K7 | `make_L27` keeps L on the 25 within 1e-10; the stored out-of-fold scores of all 25 folds are reproduced within 1e-6 | exact | B |
| K7b | each of the 25 adults, placed as a new column via `make_L27`, reproduces its stored out-of-fold score in its own fold | within 1e-6, 25/25 | B |
| K7c | each of the 10 healthy adults, binomially thinned to the smallest child's kept-T library size (seed SEED + 27), gives finite scores in all 10 healthy folds; every child's column is finite | finite; the shifts are reported | B |
| K8 | the 18 IDs and their kept T counts equal `ws_s4_pseudobulk_meta_GSE185381.csv`, read with the columns sample, group4, n, eligible, role and paediatric only | exact | B |
| K9 | the (b) code (child and comparator scores per child-fold cell) applied to A1-G02's leave-pair-out scores, with each adult AML in a child's place | A_b = 139/150 (0.926666666666667, the stored D1_AUC_LPO) within 1e-12, and "AML-side" (bootstrap seed SEED + 17) | B |

K3 rule, fixed now: sex is not a covariate, and the arms come from the deposit names. A mismatch never removes a donor
or changes an arm; it is recorded in a dated note before the primary is read. Opposite confident calls between a
donor's pool-A and pool-B libraries indicate a library swap. If the calls imply a swap between an AML library and a
donor library, the run stops before the primary for discussion, and if the swap cannot be resolved from the deposit,
(a) is "not decisive". (Loss of Y in AML blasts gives "unclear", not "female", under ws/26's rule.)

## 7. Decision rules (first match wins)

**(a):**
1. **"not decisive"**: a scope-P known answer could not be passed; K4b or K6b(ii) failed; K6b(i) failed and the result
   would otherwise be "direction replicated"; or K3's swap rule applies (section 6).
2. **"direction replicated"**: D > 0 and exact p ≤ 6/126. Add "rests on one donor (name)" if 3.2's rule applies.
3. **"reversed"**: D < 0 and the mirror exact p (labellings with D at most the observed) ≤ 6/126. Add the mirrored 3.2
   qualifier and 4.6's result.
4. **"not shown"**: otherwise, with ρ_U (3.1) in its sentence.

**(b):**
1. **"not decisive"**: a scope-B known answer could not be passed under section 6's rule.
2. **"AML-side"**: the 2.5th percentile of the 2,000 A_b* exceeds 0.5 (one-sided alpha about 0.025; simulated
   0.017-0.037).
3. **"not readable"**: otherwise.

(a) and (b) are read separately; there is no combined statistic, and (b) never changes (a)'s category. Sensitivities and
readouts never change a category.

Handover (fixed now):
- G03 runs after G04's result. It opens the September per-gene tables, runs 5.7's rungs, relates the result to 1.1's
  prediction, and takes (b)'s A_b, φ and their CIs under any outcome.
- **(a) "direction replicated"** counts toward the direction half of GATE_A3 ("at least 1 other study significantly
  concordant", HO:169) as the user decides in 9.6, with the source sentence of section 0 (and the qualifier, under
  option B only).
- Under every (a) outcome, a line in the A1-G04 GATES block supersedes the S9 verdict clause "no independent AML +
  healthy cohort exists in this work" with that outcome's wording, without editing the S9 block.
- **(a) "reversed"**: the direction claim stops; discuss before G07 or any writing.
- **(a) "not shown"**: the direction stays "shown in one study", with ρ_U as its bound. G04(c) stays a separate decision
  (9.5); it cannot separate disease from source either (addendum 1.1).
- **(a) "not decisive"**: as "not shown", without the bound.
- **(b) "AML-side"**: G13 (age) gains one sentence (section 0's wording); S9's bound stands; no gate changes. **(b) "not
  readable"** adds nothing. (b) never counts toward GATE_A3 (same study, same healthy donors).
- S10 is unchanged (it inherits the whole-T score from A1-G01).

## 8. Code, outputs, compute (plan; scripts are written after the tag)

- **`scripts/ws/41_a1_g04.R`**, modes:
  - `--build`: the Chen2023 donor table from the xlsx, and the matrices of 2.2 saved to `/LARGE1/…/ws_a1_g04/`. No
    expression-derived table goes to the repo.
  - `--known` (K0-K9): writes and prints no imprint or delay score; its checks table holds only check values.
  - `--primary` (3.1-3.2), `--sens` (4.1-4.7), `--children` (3.4, 4.8-4.11, 5.5, 5.6), `--readouts` (5.1-5.4), and
    `--combine` (summary; GATES.md block A1-G04).
  - `--smoke`: synthetic counts drawn only from the pooled kept-T profile of the 10 GSE185381 healthy adults, with
    synthetic labels for both arms and real IDs only. A check stops the run if any Chen2023 or child count column enters
    a smoke computation. Smoke tables carry no real per-donor quantity. They are written under
    `/LARGE1/…/ws_a1_g04/smoke/`, not under `results/`, so the 23:30 cron cannot publish them; they are deleted after the
    real run, and the deleted files are listed to the user.
- **`scripts/ws/submit_ws_a1_g04.sbatch`**: the A1-G01 pattern (dated read-only copy, md5 in every table,
  `OMP_PROC_BIND=false`); about 4 CPUs, 64 GB, 12 h.
- Functions are copied, not sourced:
  - `logzZ` = 33:64-70 (returns Z, U and L; W space);
  - `logzL` = 39:91-94 = 37:65-68 (returns L; `est_C`, `make_L27` and (b));
  - `black_for` (33), applied on each axis; on the Chen2023 mapped axis it removes the same 290 genes as BLACK
    (review 1);
  - `jspace(G_ref, st_ref, E_raw)` and `zE_from_counts` (33:171-185), with the reference matrix and status vector as
    arguments (T, B and NK references);
  - `load_cohort`, with `age_group_of` and the constants it uses (33:43-45 and 33:112-126), for K4 only;
  - `make_L27` (37:108-118);
  - `est_C` and its helpers (39:151-195);
  - ws/22's decontX call (22:122).
- **Not run on Chen2023 as scripts:** ws/22 (its index step would build GSM units; the call itself runs inside 41 per
  pool-B library, 4.7), ws/25 (rewrites cohort-wide S3 tables), ws/26 (no pool filter; writes per-arm composition and
  sex calls to the repo), ws/02 (sorted AML are ineligible). No re-ingest: the remap is applied to the pseudobulks.
  `config_ws.R` is unchanged. No new dependency (celda is installed).
- Seeds:
  - SEED: 4.7 (ws/22's own call);
  - SEED + 17: (b) bootstrap (A_b, φ, 4.8-4.10, K9);
  - + 18: K6b;
  - + 19: smoke;
  - + 20: 4.6;
  - + 25: K6c;
  - + 26: 5.6;
  - + 27: K7c;
  - + 28: K6a's random orders.
  - None of the offsets is used elsewhere in `scripts/ws`. Already in use: 37 uses + 6, + 11 to + 16, + 21 to + 24, + 99 and
    + 121; 39 uses + 31 to + 38, + 99 and + 35 to + 534; 40 uses + 98; 29 uses + 1 to + 10 and + 101 to + 110.
- Outputs: `results/tables/ws/ws_a1_g04_{checks,primary,loo,sens,children,readouts,summary}.csv` and the GATES.md block
  A1-G04. They are public the same night (cron); nothing G04-derived is written to the repo before the tag.
- **Allowlist** (replaces A1-BATCH's allowlist plus its "any Chen2023 file" rule). Opened before G04's result, for a
  Chen2023 donor or a GSE185381 child:
  - the inputs of 2.1-2.3 by exact path: the ten 04_annotated Chen2023 objects (counts, and the meta.data columns
    Pool, Library, bmm_fine, anno_source, high_error and hierarchy_bin), the xlsx, `features.tsv` and the remap;
  - the `PB$raw$T` columns of the 18 children, subset immediately after loading (K3b reads their XIST and RPS4Y1 rows);
  - the children's rows of `PB$cells` (cell, sample, fine, group4 and kept only);
  - the metadata named in 2.3, 4.8 and 5.5.

  No other file holding expression, or a value derived from expression, of a Chen2023 donor or a GSE185381 child is
  opened: not in the repo, on /LARGE1, in `/FAST/gr10634/gaozy/leukemia_ml`, or among the GEO deposit files. Nothing
  under `12_pseudobulk_de/` is opened (as in A1-BATCH). Closed examples, with the reason where a file could preview a
  result:
  - the public CCC tensors (`results/tables/05_ccc/tensors__bmm_broad/{Chen2023,GSE185381}`,
    `tensors_retired_2026-07-17/Chen2023`), which name P370 genes as ligands or receptors (CXCR4, HLA-A, HLA-B, ITGAV,
    PRNP, ITGA6, GDF11) with T-cell receivers; and the related `ccc_edge_distance`, `06_distance`, `graph_geometry*`,
    `07_fgw*/fgw_edges_long` and `05_ccc/cooption` tables;
  - `/LARGE1` `05_ccc_graphs{,__bmm_broad}/Chen2023` and the children's graphs;
  - the per-cell `04_cnmf/malignant/mp_usage_all_bins/{Chen2023,GSE185381}` usage (MP9 is dissociation stress and
    previews 5.2), and `12_pseudobulk_de/confound_composition.csv`, `hox_axis_scores.csv`, `marker_gate.csv` and
    `depth/`;
  - the inferCNV objects of the Chen2023 donors and the children, `infercnv_matched_ref/Chen2023`, LCC_proj,
    `03_cytotrace2`, `03_vg_malignancy`, the stemness tables, and the `03_bmm_projected` `*__cellstate_percell.csv` and
    `*__progeny_percell.csv` tables (per-cell CXCR4 and state scores);
  - `ws_t3c_cache`, `ws_b0_cache/GSE185381`, `ws_b0_inputs*`, `ws_t3c_preflight_*`, `ws_tb_coexpression_check.csv`,
    `ws_w1_samebin_ligands*.csv`, `ws_doublet_scores_percell.csv.gz`, and the truth-audit stage-2 results;
  - the ws re-integration atlas objects; `01_processed_counts/rds/Chen2023.rds`;
    `02_seurat_objects/01_per_sample_qc/Chen2023` and `03_marker_annotated`;
  - leukemia_ml `05_CCC_per_sample` (`01_per_sample_counts`, `02_cellchat*`, `03_liana*`, `07_differential`);
  - `GSE185381_RAW/*_RNA_soupx_processed.csv.gz` and `*_ADT_processed.csv.gz`;
  - the September per-gene and sample-level tables (`de/`, `hitlist_frozen.csv`, `sensitivity_pooled_bh.csv`,
    `depth_struck_hits.csv`, `perm_null.csv`, `secondary_*.csv`, `validation/`, `confound_inferred_sex.csv`) and the
    `/LARGE1` 12_pseudobulk_de pseudobulks; `05_ccc/ccc_node_features*.csv`, `07_fgw*/fgw_nodes_long.csv` and the
    `08_scoring` panel-screen tables;
  - the children's xist_umi, rps4y1_umi, sex_inferred, dbl_mean, median_umi, kept_umi, hb_share_nonery,
    gran_share_lymph and n_lymph_for_share columns of the S4 meta table (K8 reads named columns only); the UNDEMUX
    cells, any GEO lane matrix beyond the hashtag rows already tallied, and the PAUMTZ/PAXMIJ columns of any matrix.

  Anything else holding Chen2023 or paediatric expression is closed.

## 9. Decisions (the user, 2026-11-04: all seven as recommended)

1. **Chen2023's role.** Promote it from direction "reference only" to the G04(a) test; its CCC-era Discovery split label
   does not restrict this use. *(Recommended.)*
2. **(a)'s pass level** (the statistic is the mean difference in both options).
   - Pass at p ≤ 6/126 *(recommended)*: alpha 0.048; power 0.37 / 0.63 / 0.84 at 1 / 1.5 / 2 within-arm SD; at ρ 0.5,
     0.91 / 0.55 / 0.36, and at ρ 1, 1.00 / 0.97 / 0.83, for a within-arm SD of 0.30 / 0.51 / 0.70 J units (3.3).
   - Pass at p ≤ 4/126: alpha 0.032; 0.29 / 0.53 / 0.76; at ρ 0.5, 0.86 / 0.45 / 0.28, and at ρ 1, 1.00 / 0.94 / 0.75.
     This is the addendum's level, which it stated for the rank-sum (0.28 / 0.52 / 0.74 at the same alpha).
   - The choice is about alpha, not the statistic: at equal alpha the mean difference gains at most 0.02. Single-donor
     dependence is handled by 3.2's rule under either option.
   - Under 4/126, "6" becomes "4" in sections 0, 3.1-3.2, 5.2, 7, K5 and K6b, and the bounds use the 4th t_S (one-sided
     96.8 %, also in section 0's "not shown" sentence). 3.2's rule then flags 0.11 of chance passes and 0.07 / 0.05 /
     0.03 of true passes, and 9.6's rates become: (A) 0.27 / 0.51 / 0.74, 0.03 under no effect, 0.01-0.04 when one
     donor makes the pass; (B) 0.29 / 0.53 / 0.76, 0.03 and 0.06-0.07 (synthetic; `review1/main_session/`).
3. **The children.** Keep (b) as a secondary arm with the out-of-fold comparator *(recommended)*. 3.5: 2-4 % false
   AML-side under H0; power 0.62-0.68 when the children carry half of the adult gap; a shift of 1 healthy-SD shared by
   the children's lanes is called AML-side in 63 % of runs. Keeping (b) reverses the handoff plan to study the children
   without healthy comparators (HO:247) and restarts part of the paused item 3 (HO:511). The alternatives:
   - drop (b);
   - score P370 against in-sample healthy scores deflated by r. This is rejected, because the in-sample healthy scores
     sit lower (mean −1.17 to −1.36 vs out-of-fold −0.75), which favours AML-side.

   Population: Lasry's paediatric cohort (18, including three patients aged 19.6-21.1 y and AML073 without an age), with
   under-18 as sensitivity 4.8 *(recommended)*; or under-18 as the primary (14).
4. **Exclusions.** PAUMTZ out of (b); AML469 out of (a), with no rescue. *(Recommended.)*
5. **Scope.**
   - (i) The G03 rungs that add the children or Chen2023 are declared now (5.7) and run in G03 after G04's result
     *(recommended)*. The alternative is to specify them in G03, at the cost that they are designed with G04's result
     known.
   - (ii) G10(a) (descriptive comparators: ICU organ donors, AEL/MDS, CCUS) and G14(a) (S8's structural replication in
     Chen2023) are not bundled *(recommended)*. G10(a) concerns other cohorts and is descriptive. G14(a) would read
     Chen2023's T-cell structure, which is cleaner after G04's direction result. Both can join the late-November
     decision point.
   - (iii) Options (d), (e) and the GSE185381 5' T arm are declined. (d), the placement of the 3' AML cohorts against
     the band of healthy studies, is supplement-only with the study as the unit (minimum p 1/21). (e), Bailur 2020 (8
     paediatric AML vs 4 adult healthy), needs a download and is a last resort. The 5' arm is fully date-confounded
     (A1-BATCH 7).
   - (iv) G04(c) (E-MTAB-16321) became due when G01 showed a residual. It needs a download approval and cannot separate
     disease from source either, because its healthy donors are also commercial marrow (4 of 10 the same STEMCELL
     product). Decide at the late-November
     decision point (addendum section 2, item 5) *(recommended)*, or now.
6. **GATE_A3** (the direction half: "at least 1 other study significantly concordant").
   - (A) Only an unflagged "direction replicated" counts *(recommended; GATE_A3 is a stop/go gate, HO:171)*. In
     simulation this clears 0.33 / 0.58 / 0.80 at a true shift of 1 / 1.5 / 2 within-arm SD, 0.04 under no effect, and
     0.01-0.05 when one donor makes the pass.
   - (B) Every pass counts: 0.37 / 0.63 / 0.84, 0.048 under no effect, and 0.08-0.11 when one donor makes the pass. A
     flagged pass still carries a likelihood ratio of about 4-6 against no effect.
   - Either way, this replaces the handoff's gene-level measure (HO:160) with the sample-level test (5.4 keeps the
     gene-level measure, descriptive), and (b) never counts. A flagged pass under (A) is reported, and GATE_A3's
     direction half stays open.
7. **Deciding space for (a).** J *(recommended: ρ and its bounds then agree with the category; S9's precedent; depth is
   reported by K6c and equalised in 4.6)*; or W (the same power and no depth bias, with ρ still reported in J).

## 10. Review history

- **DRAFT v1**, 2026-11-04 (commit 8eceb09): written from the read-only inventory (workflow wf_e07a2214-567; notes outside
  the repo).
- **Review 1**, 2026-11-04 (workflow wf_8890f84c-931). Four lenses, each with an adversarial verifier: A, the statistics
  of (a); B, the children's design; C, feasibility against the code and files; D, consistency and disclosure.
  - 58 findings, none refuted: 46 real, 12 partly real. Three were blocking after verification: A1 (the leave-one-out
    qualifier could fire only at p = 6/126), and C1 = D6 (K0's `features.tsv` md5 was that of the decompressed file, and
    the known-answer failure rule departed from A1-G01's three states, so a typo could have voided the primary).
  - The completeness critic hit a usage limit before returning. Its working notes were used: the conflicts between lens
    fixes, the status of the blind-safe G03 rungs, two metadata reads (1.3) and a GATE_A3 simulation.
  - Review 1 read only metadata, file names, code, GSE185381 adult data and synthetic data; it read no Chen2023 or child
    expression, score or September per-gene value. Notes, scripts and outputs are in
    `/FAST/gr10634/gaozy/reports/2026-11_a1_g04_prereg/review1/` (outside the repo). The scratch evidence of the
    inventory was copied to `inventory/scratch_evidence/` there.
- **DRAFT v2**, 2026-11-04: all findings folded in. Conflicts between lens fixes were resolved as follows.
  1. The qualifier is A1's donor-influence rule (D8's corrected rate for the old rule is moot).
  2. 9.2 is a choice of alpha with the mean difference fixed (A4); the rank-sum figures are kept for reference (D7).
  3. K6 is split. K6a is deterministic. K6b keeps machinery bars, one tail, at 1,000 draws, so that an exactly
     calibrated run almost never fails the bar (with 200 draws it failed 1.6 % of the time, D21). K6c is a depth
     calibration with no bar. The reviewer's common-factor K6b was not adopted: its median factor of 0.12 exaggerates
     the depth bias about 4 times (A2 verifier). Plants never touch Chen2023's real columns.
  4. B/NK use J space with each compartment's own GSE185381 reference (C4), not W (A5), so that 5.3 matches the
     primary's space.
  5. The sensitivities are renumbered: (a) 4.1-4.7, (b) 4.8-4.11.
  6. K3: a mismatch never removes a donor, and a cross-arm swap pattern stops the run before the primary (A8 and the C
     verifier).
  7. The allowlist names K3b's per-cell XIST/RPS4Y1 rows explicitly (B verifier).
  8. The (a) licence sentence follows D1's verifier wording, plus 5.2's handling clause (A7, D12).
  9. The GATE_A3 numbers were recomputed for the new rule (critic; main-session check
     `review1/main_session/sim_alpha4_influence.out`, which also gives 9.2's 4/126 rates).
  10. Seeds + 20, + 25, + 26 and + 27 were added.
  11. φ's CI is conditional on the 15 adult AML out-of-fold scores.

  Added beyond the findings: 4.7 (decontX per pool-B library), the D5 verifier's alternative to declining the ws/22
  obligation.

  Not adopted:
  - leave-two-healthy-out comparators for (b): infeasible, because 6 of 45 distinct nested refits (12 of 90 rows)
    collapse, each removing Control3 or Control4 with another P02 donor, or Control0082 with Control4 (B10);
  - weakening (b)'s one-sided reading (the B4 verifier: paediatric AML differs in drivers and T biology).

- **Review 2** (short), 2026-11-04 (workflow wf_6477cfeb-0ac). Three read-only checkers: fold-in coverage; consistency
  and numbers, including every md5; feasibility and blindness. Notes in
  `/FAST/gr10634/gaozy/reports/2026-11_a1_g04_prereg/review2/` (outside the repo). No Chen2023 or child expression,
  score or count was read.
  - Blocking, found by all three: the md5 recipe for the 63 GSE185381 author-metadata files was not stated, and the
    natural reading gives a different value. Fixed: the exact command is in 2.3.
  - Should-fix, fixed:
    - K3 is reported, but its swap rule now gates P (section 6 says so).
    - 7(a)1 gains the route for a scope-P check that cannot be passed; section 0 matches.
    - 5.7's rungs leave out PAUMTZ and PAXMIJ, so rung "+Chen2023" is no longer the September fit and its count is not
      known.
    - 3.3's real-biology power is described at the donors' own depth (0.51), with 0.43 when the AML arm is thinned.
    - The allowlist names `Library`, `hierarchy_bin` and the children's `PB$cells` rows.
    - 5.6 and K7c use binomial thinning only, and only adults at least as deep as the child. Cutting adults to deeper
      children is impossible, and cell subsampling would need undeclared per-cell objects.
    - K3b's per-cell hashtag split is dropped. Its only source, the GEO lane matrices, also holds the children's
      expression, and the split cannot change a category; the pseudobulk-level call and the other-sex UMIs are printed
      instead.
  - Minor, fixed:
    - A K6b(i) failure blocks only a would-be pass.
    - 9.2 gives ρ-scale power and the 4/126 rates for 9.6 (`review1/main_session/sim_rho_power_alpha.out`).
    - K9's expected value is exact.
    - Declared: 4.6's per-draw recomputation; 4.7's rounding after summing; K6c's probe cohort, interpolation and
      printing; meta_GSE120221_v1.0.csv; seeds SEED (4.7) and + 28 (K6a).
    - 4.11 adds the 7-fold version.
    - Counts and wording corrected: 278 author-malignant kept T cells; 46 real and 12 partly real review-1 findings;
      AML003's source; E-MTAB-16321's healthy donors; the (b) licence names the cross-fitted lists.
- **DRAFT v3**, 2026-11-04: review 2 folded in.
- **v1 APPROVED**, 2026-11-04: the user decided all seven section-9 items as recommended. With the approval, a dated
  note in the handoff S15 row narrows its 2026-11-01 sentence (D2). Tag `prereg-a1-g04-v1`.
- **Code review 1 and dated pre-run note**, 2026-11-04 (before any real-data mode of 41 was run; no Chen2023 or child
  score had been computed). `scripts/ws/41_a1_g04.R` and `submit_ws_a1_g04.sbatch` were written in three stages, each
  run end to end on synthetic data only (smoke jobs 3683524-3683538; smoke outputs under `/LARGE1/…/ws_a1_g04/smoke/`).
  Three read-only reviewers checked conformance with this document, correctness on real-data edge cases, and blindness
  and output safety (notes in `/FAST/gr10634/gaozy/reports/2026-11_a1_g04_prereg/code_review1/`, outside the repo). Two
  of them found the same blocking bug (`--combine` stopped when no statistical bar failed); it is fixed and that path was
  rerun on synthetic data. Fixed as declared above: the "not decisive" routes of 7(a)1 and 7(b)1 when a gate stops a mode;
  K3b's mismatched children leave the 4.10 set; the statistical bars (K4b, K6b) run only after every scope-P
  deterministic check has passed; C_H, C_A and C stay off the public readouts table unless (a) is "direction replicated";
  the "not shown" sentence always carries ρ_U; φ's bootstrap uses the drawn folds' child scores; K9 and (b) share one
  bootstrap function. Implementation details declared now, none of which can change a category:
  - Seeds: K5's synthetic test vectors use SEED + 28 (as K6a); every smoke draw derives from SEED + 19.
  - K2 checks the saved matrices' column sums against sums computed separately over the cells (T, pool-A T, B, NK and
    the decontX matrix).
  - K3 prints the Chen2023 sex calls only (not the UMIs); K3b prints, for AML056, AML038, AML026 and AML051, the UMIs of
    the gene of the other sex only.
  - A K3 mismatch stops `--primary`, and a K3b mismatch stops `--children`, until a dated note is added here; the mode is
    then rerun with `--k3-noted` or `--k3b-noted`.
  - The author-metadata join of 3.4 stops the run if any child matches fewer than 95 % of its kept T cells, or if the
    two counts disclosed in 1.3 are not reproduced (AML028: 16 author-NK of 20; AML073: 210 author-malignant of 605).
  - The 5-gene rule is applied to every J-space score of sections 4 and 5 (D and p are NA below 5 genes in a non-empty
    half; n_up and n_down are reported).
  - Every mode checks all of its output files before any work and refuses to overwrite one; a later mode stops unless
    the checks table was written by the same script version. A rerun after a fix therefore needs a dated note here and the
    earlier output renamed (`*_superseded_<date>`), never deleted.
- **Run 1 and a gate fix**, 2026-11-05. Run 1 (2026-11-04, script md5 `5385eb29…`):
  - `--build` (3683539) and `--known` (3683540) completed. The checks table has 70 rows: 0 deterministic FAIL, 0
    statistical FAIL and 0 swap stops. K3: all ten Chen2023 donors match in both pools. K3b: all 18 children match.
  - `--primary` (3683542), `--sens` (3683543) and `--readouts` (3683544) completed.
  - `--children` (3683545) stopped and named all 18 children as K3b mismatches. `--combine` (3683546) could not start
    and was cancelled.

  Cause: the K3b gate in `--children` and the K3 gate in `--primary` read the column `state`, which is always
  "REPORTED" for a reported row, instead of `got`, which holds the match call. So the K3b gate flagged every child, and
  the K3 gate could never fire. The K3 gate changed nothing in run 1 because all ten K3 calls matched.

  The same bug had emptied the flag-free set in smoke job 3683534. It was then patched where it surfaced (a guard in the
  bootstrap), not at its source.

  Fix: both lines now read `got`. The K3 gate, like the K3b gate, is skipped in a smoke run (smoke calls are synthetic).
  Checked on the stored tables: the real table gives 0 K3 mismatches and 0 flagged children; the smoke table gives 10
  and 14, as written. The before-fix copy of the script is in `reports/2026-11_a1_g04_prereg/code_review1/` (outside
  the repo).

  No run-1 analysis table (primary, loo, sens, readouts) was opened before the fix. The script's md5 changes, so every
  mode from `--known` on is rerun:
  - Run-1 outputs are renamed `*_superseded_20261105`. They were already committed by the 2026-11-04 23:30 snapshot.
  - `chen_pb.rds` is kept, because the `--build` code is unchanged.
  - Run 2's primary, loo, sens and readouts tables should equal run 1's apart from `script_md5`. Any difference will be
    reported.

## 11. Outcome (2026-11-05)

Run 2: jobs 3683796-3683801, script md5 `9db3b2c6…`, commit 19813ba. It ran after the gate fix in section 10.
- Run 2's checks, primary, loo, sens and readouts tables equal run 1's, apart from `script_md5`.
- `readouts_C.rds` is identical in both runs.
- The smoke outputs were deleted after this check.

**Known answers.** The checks table has 70 rows: 0 deterministic FAIL, 0 statistical FAIL and 0 swap stops. All K3 and
K3b sex calls match, so no child left the 4.10 set.

**(a) Chen2023: "not shown".**

| quantity | value |
|---|---|
| D (J units, AML minus donors) | 0.712 |
| exact p, upper / mirror | 17/126 / 110/126 |
| rank-sum p | 18/126 |
| ρ [ρ_L, ρ_U] | 0.518 [−0.247, 1.361] |
| qualifier | none (3.2 applies only under "direction replicated") |

Wording (section 0): "The direction was not shown in Chen2023 (exact p = 17/126); a common shift larger than 1.36 of
GSE185381's selection-corrected gap is excluded (one-sided 95.2%, common-shift model)." "5 vs 4 cannot exclude an
effect of the GSE185381 size." The direction stays "shown in one study".

Recorded although it changes nothing:
- Under 3.2's rule, dropping AML162 leaves D = 0.225, below 2/3 of D. A "replicated" result would therefore have
  carried "rests on one donor (AML162)".
- Per-donor J scores: AML162 +1.42; the other four AML −1.41 to −0.45; the four donors −1.71 to −0.60.

Sensitivities 4.1-4.7 (no category change):
- D 0.68-0.94, exact p 16/126 to 20/126.
- The up half and the down half alone give D 0.43 (p 20/126) and 0.28 (p 13/126).
- 4.6, depth-equalised, J: D 0.680, p 18/126. K6c's implied depth bias is −0.004 J.

Readouts (descriptive):

| readout | result |
|---|---|
| 5.1 composition | gapT 0.681 u; gapC_H 0.041; gapC_A 0.062; mean residual 0.630. The fallbacks are exactly the declared ones. C_H, C_A and C are not reported. |
| 5.2 delay score | D_delay 0.583, p 30/126. Spearman with the imprint score over the 9 donors, after centring each arm: 0.83. AML162 has the highest delay score (+1.28). |
| 5.3 B (5 vs 4) | D 0.347, p 27/126, D/D_T 0.49 |
| 5.3 NK (4 vs 4) | D 0.462, p 22/70, D/D_T 0.55 |
| 5.4 sign concordance | 0.832 of 357 genes, p 18/126 |

**(b) children: "not readable".**
- A_b 0.372 [0.161, 0.594]; φ −0.244 [−0.827, 0.092].
- Subsets 4.8-4.11: A_b 0.30-0.42.
- Under section 7 there is no statement. A_b, φ and their CIs go to G03.

**Gates and handover.**
- GATE_A3, direction half: does not count (option A; (a) is not "direction replicated").
- `--combine` appended the A1-G04 block to `GATES.md`. Its last line supersedes S9's clause "no independent AML +
  healthy cohort exists in this work".
- Per section 7: the direction stays "shown in one study", with ρ_U = 1.36 as its bound.
- G04(c) remains a separate decision (9.5). G03 is next.
