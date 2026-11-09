# PREREGISTRATION A1-G07: what the T-cell axis is (a pre-registered description)

**Status: DRAFT v1, 2026-11-09. Not approved; no tag.** Nothing declared in sections 3-6 has been computed. The
counts in section 2 use gene names, set memberships and cell numbers only (1.4).

**Sources.** Read-only inventories kept outside the repo, in
`/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/`:
- `A_obligations.md` (A);
- `B_machinery.md` (B);
- `C_resources_prior_art.md` (C);
- `D1_geneset_design.md` (D1, synthetic gene-set study);
- `D2_singlecell_layer.md` (D2, synthetic single-cell study);
- `critic.md` (critic: 44 citations checked, no blindness breach found).

Membership counts for this draft are in `.../drafting/membership_counts.py` and `.out`.

**Place in the sequence.** A1-G02 (PASS) → A1-BATCH (PASS with qualifications) → A1-G01 ("beyond T-subset
composition, detected inside the main labels") → A1-G04 ((a) "not shown", (b) "not readable") → A1-G03 ("held through
R6") → **this**, the last item of handoff v2.0 §17 S15.

The late-November decision point (assessment §5 item 6) has not been held. G07 is written now for two reasons:
- every interpreted part uses resources already on disk;
- downloaded lists enter only a descriptive table (3.6(b), decision 9.5).

## 0. Question, scope, and what each outcome lets the paper say

The frozen T axis is a gene list today:
- **P370**: 156 up, 214 down;
- **F402**: the 402 S5 hits before the post-hoc removal of 32 ambient-prone genes.

**G07 asks: what kind of program is it?** It answers at four levels, plus two descriptive tables:
1. **Gene sets (3.1-3.2).** Which pathway sets and T-cell-program sets carry more of the S5 AML-versus-healthy T contrast
   than the average tested gene? For each such set: does it rest on genes that also change with blood-processing delay
   or immediate-early stress?
2. **Compartments (3.3).** For each such set: is the shift larger in T than in B or NK cells of the same donors, or shared?
3. **Gene identity in healthy marrow (3.4).** Which axis genes mark naive, central-memory or effector-memory T cells in
   healthy marrow?
4. **Cells (3.5, the S5b layer).** Is the up-set shift carried by most T cells or by a minority?
5. **Tables (3.6).**
   - (a) The cross-fitted score with gene classes removed: lncRNA, clock/glucocorticoid, handling.
   - (b) Overlaps with lists on disk and, if approved, with published lists.

**What G07 is not.**
- **Not a gate.** No G07 result changes the existence, direction, composition, batch or bridge results (A1-G02, A1-BATCH,
  A1-G01, A1-G04, A1-G03).
- **Not independent evidence.** Every set test re-describes the S5 contrast that selected the genes: same donors, same
  design, same noise. The text says so every time (7.1).
- **Not a test of cause.** It does not test AML causation, T-specificity, handling, age or ambient RNA. These stay as
  worded in GATES.md (7.7).
- **Not blind at the level of gene identity (1.2).**
  - Already seen: the HALLMARK membership counts of the 370, the Savage-overlap counts and most gene classes.
  - Not yet seen: any test on the full ranked statistic, any per-cell distribution, overlaps with full published lists,
    the healthy-subtype identity of the 370, and set tests in B, NK or monocytes.

**What each outcome lets the paper say.** The fixed sentences are in section 7.
- **Sets flagged (3.1).** "In the S5 T-cell contrast, set X had a higher (lower) average AML-versus-healthy statistic
  than the average tested gene." Each flagged set carries:
  - its handling reading (3.2);
  - its compartment reading (3.3);
  - the clause "(same contrast that defined the axis; descriptive)".
- **No set flagged.** "No set differed from the average tested gene at FDR 0.05 with a test that accounts for inter-gene
  correlation (n = 25)." In D1's simulation this outcome is likely even when real set structure exists (1.4), so it is
  never read as "the axis is not a pathway".
- **Handling-carried (3.2, 7.3).** If at least half of the flagged sets in one direction are carried by the listed
  handling genes, or cannot be separated from them, the paper states as a main limitation: the set-level description of
  that arm cannot be read as AML biology.
- **S5b (3.5).** The cohort reading is "mostly uniform", "mostly sub-population", "mixed" or "not readable", with
  per-patient calls. "Sub-population" means "carried by a minority of cells", never "a new cell state".
- **Gene identity (3.4) and the tables (3.6)** are descriptive under every outcome.

### 0.1 What G07 owes, where it is handled, and declared departures

| source | obligation | here |
|---|---|---|
| addendum 1.6; verify_G07 item 9 | disclose the 10-28 look at the HALLMARK membership counts | 1.2 |
| verify_G07 best option; objections B1-B2 | the description predicted from gene names | 1.3, stated as a prior expectation |
| addendum 1.6; verify_G07 item 10 | camera/fry primary on the S5 design, fgsea secondary | 3.1 (departure 1); fgsea and fry reported (4.4) |
| verify_G07 item 10 | universe = tested genes minus extended-ambient genes | 2.2; all 7,379 as a sensitivity (4.1) |
| verify_G07 item 10 | report both the 370 and the 402 | 2.1; every membership, identity, class and S5b table |
| addendum 1.6; verify_G07 item 7 | PROGENy/CytoSig only with a recorded user decision | not used (departure 4; decision 9.3) |
| verify_G07 items 8, 11 | a durable overlap table from full lists, anchored on G02's always-selected core | 3.6(b) |
| verify_G07 best option | per-gene specificity across healthy T subtypes | 3.4 |
| S5 4.2 and :135; verify_G07 item 4 | S5b was registered and never produced: run it or disclose it | 3.5, a dated amendment of S5 4.2 (departure 2; decision 9.6) |
| verify_G07 best option; verify_G06 | where the GNA15/HMGB3 signal sits | 3.5(h), description only |
| verify_G06 | a score without lncRNA-like genes; one clock/glucocorticoid line | 3.6(a) |
| memory notes on Control4 and on the T-ALL pools; G03 disclosure 12 | declare the sex coding and the T-ALL handling before any result | 2.1, 4.2, 4.3 |
| A1-BATCH 12.5; A1-G01 11.5 | not T-specific; delay overlap; "beyond composition" does not mean "caused by AML" | 3.3; 7.6-7.7 |
| A1-G03 :466 | per-gene list status "stands", handed to G07 | used as is; claims at set level (7.4) |
| A1-G02 :116-118 | the 25-gene always-selected core | 3.6(b) |
| assessment 2.2 | generic labels may reopen the composition reading | 3.4; 7.4 |
| verify_G07 item 13 | GSE227903 diagnosis rows only | not used |
| verify_G08 item 4 | overlap of the 370 with Lasry's lists | 3.6(b) if the file is on disk before the run (9.5); otherwise G08 |

**Declared departures.**
1. **fry is not co-primary.** The verifier asked for "camera/fry". In D1's synthetic study (1.4), once a genome-wide donor
   component was present:
   - fry's self-contained "Mixed" test flagged 71-89 % of sets that carried no effect (FDR < 0.05);
   - fry "directional" drifted to 0.108 at p < 0.05.

   Only camera with estimated inter-gene correlation held its error rate in all three scenarios. fry directional is
   reported (4.4).
2. **S5b gets a reading pipeline and a rule.** S5 4.2 registered two summaries: the patient mean UCell, and the fraction
   of cells above the healthy same-label P95. In D2's synthetic study:
   - they cannot separate a uniform shift from a sub-population at a matched pseudobulk effect;
   - AML's deeper cells move them by +0.45 healthy SD with no expression change at all.

   G07 reports them as registered (pipeline P1). It reads the question with a depth-equalised pipeline (P3) and a shape
   rule (3.5). The amendment is recorded here because S5 is tagged; S5's text is not edited.
3. **S5b runs after G01.** The verifier placed S5b before G01's result was interpreted. G01's verdict stands whatever S5b
   shows. A sub-population call is described with G01's limits (7.5).
4. **No footprints and no TF activity** (if decision 9.3 is as recommended). GATE_W1.1 and GATE_B0 failed in these data
   (GATES.md:24, 83), and the paper outline lists the footprint control as a negative result (handoff :593).
5. **No dissociation list.** None is on disk. The handling union is Savage 2021 plus the project IEG list (2.3); A1-BATCH
   D8 scored those two only.
6. **Over-representation on the 370 is not a test.** Membership counts are reported with their expected counts and
   without p values. They were already seen (1.2).
7. **G06 leftovers: two description lines only (3.6(a)).**
   - Disclosed omissions: de novo patients only, a 50-cell floor, the GATE_A1(4) re-read, and a T-cell doublet
     drop-and-refit.
   - Already ran: the median-UMI covariate (A1-BATCH W2).
8. **The per-gene identity reference is the GSE185381 healthy donors.** They are one arm of the contrast that selected
   the genes. Oetjen is a sensitivity (3.4).

## 1. Already known before this draft (disclosed)

**1.1 Results G07 builds on** (GATES.md blocks; numbers as recorded there).

| step | result |
|---|---|
| S5 | T: 402 genes (171 up / 231 down), FDR < 0.10 and \|log2FC\| > 0.5. Fit: edgeR QL, `~ status + hb + gran + sex + CD8 share`, 7,379 genes tested. GATE_A1(2) failed (32 genes in the extended ambient ranking), giving P370 (S8). B 9 genes (not called), NK 91, Mono 275. S5b not produced. |
| A1-G02 | cross-fitted gap 1.378 (1.98 out-of-fold healthy SD), p 0.008, AUC 0.93; 25 genes selected in all 25 folds (10 up, 15 down, all in P370) |
| A1-BATCH | PASS. Direction shared with B and NK (R 0.69 / 0.93), so not T-specific. Overlaps the Savage delay program (D8: up 12 vs random 6; down 23 vs 10); delay score in MAIN gives R 0.82. IEG score +0.04 z (p 0.81). Post hoc: gap 1.318 without Savage genes, 1.323 without Savage + IEG genes. |
| A1-G01 | C 0.045 (−0.100 to 0.191); S_w 0.664 (0.497 to 0.772). Out of sample, the direction appears in B 0.857, NK 0.893 and Mono 0.732. 5.5 naive-to-effector axis not fired. Gamma-delta share 0.6 % (AML) vs 2.6 % (healthy). UMI per kept cell 1.15-2.07 times higher in AML in every label. |
| A1-G04 | Chen2023 "not shown" (p 17/126; ρ_U 1.36); children "not readable"; the direction stays "shown in one study" |
| A1-G03 | held through R6; per-gene list "stands"; the count is volatile when one AML is left out (38-1,106); with Control4 coded female, 402 becomes 895 |
| S8 / S9 | structure replicated in two 10x cohorts (T2, T1); age share bounded (f_age 0.102, CI −0.085 to 0.463) |

**1.2 Prior looks at the 370/402 against annotation** (ledger L1-L17 in A section 5; every item is disclosed).
- **HALLMARK membership counts, read 2026-10-28** (pinned v2024.1; counts only, no test):
  - E2F targets 7 up / 3 down;
  - UPR 6 / 0 (ATF4, BAG3, EIF4E, MTHFD2, SHC1, ZBTB17);
  - mTORC1 7 / 2;
  - IFN-alpha response 0 / 4 (IFIT3, IFITM1, TRIM14, MVB12A);
  - TNFA via NF-kB 5 / 1 (up includes DUSP4 and RELA; down IER2);
  - OXPHOS 3 / 6.
- **17 named-example lists, by gene name** (the table is lost; its prose is in `inventory_prior_art.json` items[13]).
  - Zero overlap: exhaustion, senescence and cytotoxic examples; Knaus up genes; Lasry's named genes.
  - Overlaps on the down side: Azimuth naive/TCM top-10 (LINC02446, NOSIP, RGS10, LDHB), gamma-delta (TRDC), Knaus down
    (TMIGD2), Szabo IFN (IFIT3).
  - Savage: DUSP4, up.
- **Gene classes by name.**
  - About 35 lncRNA-like genes down vs 3 up; only IER2 is on the IEG list.
  - Clock: PER2 and ARNTL up, DBP down. FKBP5 down.
  - HLA-A and HLA-B up.
  - Flat: canonical ISGs, and canonical exhaustion and senescence genes.
  - A sub-threshold effector lean: GZMB, GZMH, FGFBP2.
  - The identities of the 32 ambient-prone genes: SELL, IL7R, FHIT down; CD5, RUNX3, MAF, SRGN, CXCL8, ANXA2, EZR up.
- **Gene-level readouts in other steps.**
  - S9: per-gene age effects of the 370 in Oetjen.
  - A1-BATCH: D8 overlap counts; D9 soup ratios (not interpreted); post-hoc genome-wide T-versus-B/NK correlations.
  - G01: the 5.5 loading statistic.
  - G04: Chen2023 gene-level concordance and delay score.
  - G03: per-rung P370 readouts and the stored September table.
  - S10 review: the 370 self-identify better than about 94 % of expression-matched random sets in healthy Griffin
    libraries.
- **Not seen:**
  - any camera, fry or fgsea run on the ranked statistic;
  - any per-cell score;
  - overlaps with full published lists;
  - the healthy-subtype identity of the 370;
  - set tests in B, NK or monocytes.

**1.3 Prior expectations** (stated now; they are not hypotheses).
- 2026-10-28 (verifier; objections inventory): "post-activation and biosynthetic-stress genes up; a quiescent,
  naive-like layer and lncRNAs down; a slight effector lean."
- 2026-11-09 (inventory C, from values already quoted elsewhere): the same, plus "no exhaustion, no senescence, no
  ISG-up". C names processing/handling as the strongest competing reading, because it shares the up arm's most
  interpretable genes.

These descriptions were guessed from gene names. A set-level result that matches them is not a confirmation. The text
quotes the set-level result, not the guess.

**1.4 Seen while writing this draft** (no AML-versus-healthy expression value was computed).
- **Cell QC** (D2 2; B 4.3).
  - Kept T cells: 14,609 in total. AML 6,911 in 15 patients (20-1,208 each); healthy 7,698 in 10 donors (123-1,917).
  - AML T cells are deeper: median UMI 2,952 vs 2,000; genes detected 1,100 vs 919.
  - W4 doublet score in kept T cells: median 0.095 (AML) vs 0.075 (healthy). Its Spearman correlation with log UMI is
    0.62 in AML and 0.25 in healthy.
- **Membership counts** (names only; `drafting/membership_counts.out`).
  - Universes 7,379 and 6,957 (422 extended-ambient genes removed). 45 HALLMARK sets have ≥ 15 genes.
  - The handling union H has 698 genes in the ambient-free universe (base rate 0.100).
  - 12 HALLMARK sets have ≥ 25 % of their genes in H.
  - **20 of the 22 Gavish T-cell programs have ≥ 25 % of their genes in H**, including CD8 MEMORY_NAIVE_1 (0.76) and
    CD4 NAIVE_1 (0.74).

  So the published processing-delay list shares most genes with T-cell differentiation programs. At the level of gene
  sets, a delay reading and a naive-loss reading cannot be separated with these lists (7.3).

  Also counted: 305 / 274 lncRNA genes in the two universes, and the clock/glucocorticoid coverage (2.3).
- **D1 synthetic gene-set study** (D1 section 2; 45 HALLMARK-sized sets, 15 vs 10, S5's engine, 30 replicates per
  scenario).

  False-positive rate on sets without an effect, at p < 0.05 / at FDR < 0.05:

  | method | no genome-wide component, ρ 0.05 | genome-wide, ρ 0.05 | genome-wide, ρ 0.20 |
  |---|---|---|---|
  | camera, estimated correlation | 0.041 / 0.004 | 0.024 / 0.001 | 0.011 / 0.000 |
  | camera, 0.01 | 0.173 / 0.084 | 0.093 / 0.020 | 0.291 / 0.189 |
  | fry, directional | 0.056 / 0.010 | 0.108 / 0.032 | 0.069 / 0.008 |
  | fry, Mixed | 0.116 / 0.044 | 0.917 / 0.892 | 0.771 / 0.714 |
  | fgsea, preranked | 0.235 / 0.175 | 0.167 / 0.091 | 0.377 / 0.310 |

  - **Power of camera with estimated correlation** (FDR < 0.05) for a whole-set shift of 0.30 log2: 0.91, 0.32 and 0.05
    in the three scenarios.
  - **Cross-compartment interaction** (synthetic T + B):
    - T-only sets called "larger in T" 0.80 of the time; shared sets 0.00; null sets 0.002.
    - The declared reading classified flagged sets correctly in 0.94 (T-only) and 1.00 (shared) of cases, under fully
      shared noise.
    - Power asymmetry (a weaker or noisier effect in B) was not simulated.
- **D2 synthetic single-cell study** (D2 sections 4-5).
  - The registered pipeline P1 shifted null AML patients by +0.45 healthy SD and roughly doubled the share of cells
    above the healthy P95, from depth alone.
  - Thinning cells to 1,500 UMI (P3) left a +0.15 SD shift.
  - With the rule in 3.5(e) on P3:
    - no false "sub-population" call in 384 uniform synthetic patients;
    - a 20 % sub-population was called at rates of 0.95-1.00 (moderate effect) and 0.78-1.00 (weak effect) with ≥ 50
      cells;
    - null AML patients passed the shift gate in ≤ 0.025 of cases.
  - Uniform worlds with heterogeneous healthy cells, or with low-expressed sets, often ended "not separable" (0.33,
    0.83).

**1.5 Disclosures carried from earlier steps** (A section 6; each enters the paper's limitations).
1. **The 32-gene removal was post hoc.**
   - The 32 carry most recognisable lineage markers. So naive or effector labels can differ between P370 and F402, and
     both are reported.
2. **The narrowing to T cells chose where the axis is measured.** It was not a finding of T-specificity. B's
   "no imprint" (9 genes) is not evidence that B cells are unaffected.
3. **T3 is recorded as uninformative.**
   - No statement compares CD4 with CD8.
   - G01's S_k (0.72-0.93 in four labels) is the within-label evidence.
4. **The hashing covariate was mis-coded.** Gene membership is a property of one design's selection; with a multiplexed
   covariate the list was 405 hits, 308 of them shared with F402.
5. **Control4's sex.**
   - Metadata male, cells female.
   - Sample-level results barely move with the coding; the per-gene list moves a lot (402 → 895).
6. **Doublets and decontX.**
   - The registered doublet arm never ran (replaced by A1-BATCH W4).
   - The decontX arm ran in A1-BATCH D4, where GNA15 left the up-arm interpretation.
7. **Healthy lanes P02 and P03** (Control1-5) were captured with T-ALL/ETP-ALL samples. The T-lineage soup was not
   removed, and A1-BATCH D9 could not bound it.
8. **The healthy donors are 31.5 years younger in median**; age is bounded, not excluded (S9). CMV serostatus, infection,
   steroids and time of aspiration are unrecorded in both arms.
9. **The pre-registrations are internal git tags** with no public time stamp.

## 2. Data and frozen definitions

**2.1 Samples, fits, codings.**
- **T**: S5's 25 samples (15 adult AML, 10 healthy).
  - Design table: the T rows of `ws_s5_scores.csv`. Counts: `PB$raw$T` of
    `/LARGE1/gr10634/gaozy/aml_niche_net/ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds`.
  - Fit:
    - S5's blacklist (307 genes) before filtering;
    - `filterByExpr(group = status)` plus CPM > 1 in at least half the samples;
    - TMM;
    - `estimateDisp(robust = TRUE)` and `glmQLFit(robust = TRUE)`.

    This is 27's `fit_de`, as G03's K1 reproduced it, with `return_y = TRUE`.
  - Design columns are `(Intercept), statusAML, hb, gran, sexmale, subtype_share`. The tested coefficient is passed by
    name, `"statusAML"`, never by position: edgeR's `camera.DGEList` and `fry.DGEList` default to the last column (B 1.3).
- **B** (13 + 10), **NK** (12 + 9), **Mono** (16 + 10; 21 of them in T's 25): each refitted with its own S5 design and
  `subtype_share` (immature-B, CD56high NK, CD16 Mono).
- **Sex**: metadata coding, as in S5, G02 and G03, primary. Control4 coded female is a sensitivity (4.2). Level names are
  `female` and `male`.
- **T-ALL pools**: Control1-5 (lanes P02/P03) are flagged in every table. Leaving them out is a sensitivity (4.3).
- **Lists quoted**: P370 is primary. F402 stands beside it in every membership, identity, class and S5b table.

**2.2 Gene universes.**
- **U_T** = the 7,379 genes tested in S5's T fit.
- **U_A** = U_T minus the extended ambient ranking: `ws_s1_ambient_blacklist.csv` rows with `n_units ≥ 2`, as at
  `28_s5_coordination.R:144`; 926 genes, 422 of them in U_T. U_A has 6,957 genes and is primary for 3.1-3.3.
- B, NK and Mono use their own tested genes minus the extended ranking. The T × X interaction (3.3) uses
  (U_T ∩ X's tested genes) minus the extended ranking.
- **Implementation.** Fit on the full filtered DGEList, exactly as S5. Then subset rows to the universe with
  `keep.lib.sizes = TRUE` before camera. NB z-scores are computed gene by gene with fixed dispersions, so subsetting
  changes only the background.

**2.3 Gene-set families and gene lists (frozen).**
- **P-H, HALLMARK v2024.1 (pinned).** `/LARGE1/gr10634/gaozy/reference/h.all.v2024.1.Hs.symbols.gmt`, md5
  `d3606c32…`, sha256 `ee2463…` (PROVENANCE beside it).
  - 45 sets have ≥ 15 genes in U_A.
  - Not testable (< 15): ANGIOGENESIS 5, PANCREAS_BETA_CELLS 9, HEDGEHOG 10, APICAL_SURFACE 13, NOTCH 13.
- **P-T, Gavish 2023 3CA T-cell metaprograms (22 sets: 10 CD4, 12 CD8).** Decision 9.4.
  - Source: the repo's MSigDB v2026.1 GMT, `LCC_proj/scripts/panels/msigdb.v2026.1.Hs.symbols.gmt`, md5
    `30893650e0c0735a39819f9f7691f0d6`. It was committed in aac665f (2026-08-04) for another sub-project and has no
    provenance file.
  - If approved, the 22 sets are extracted to `/LARGE1/gr10634/gaozy/reference/gavish2023_3ca_T.v2026.1.gmt`, with a
    PROVENANCE note and an asserted md5.
  - All 22 have ≥ 15 genes in U_A (17-46). They are tumour-infiltrating T-cell programs (C 2.2).
- **H, the handling union.**
  - Members: the Savage 2021 T delay-18-h consistent set (792 genes; `savage2021_T_delay18h_consistent.csv`, md5
    `5ee444a9…`, as asserted at `37_a1_batch.R:372`) plus the project IEG list (26 genes, `config_ws.R` `WS_B0_IEG`).
  - 698 genes in U_A (base rate 0.100).
  - **Handling-adjacent** = at least 25 % of a set's U_A genes are in H (2.5 times the base rate). This flag is fixed now
    and does not depend on the outcome.
  - HALLMARK, 12 sets: TNFA_SIGNALING_VIA_NFKB 0.47, ALLOGRAFT_REJECTION 0.36, HYPOXIA 0.35, IL6_JAK_STAT3 0.30,
    INTERFERON_GAMMA 0.30, INFLAMMATORY_RESPONSE 0.30, EMT 0.29, INTERFERON_ALPHA 0.28, APICAL_JUNCTION 0.27, APOPTOSIS
    0.27, MYOGENESIS 0.25, TGF_BETA 0.25.
  - Gavish, 20 of 22: all except CD4 CELL_CYCLE (0.24) and CD8 NAIVE_3 (0.23).
- **lncRNA** = `gene_type "lncRNA"` (and no other type) in `refdata-gex-GRCh38-2020-A/genes/genes.gtf`, the annotation of
  the data's own gene axis. 305 genes in U_T, 274 in U_A.
- **Clock and glucocorticoid list** (declared from canonical knowledge; the axis status of PER2, ARNTL, DBP, FKBP5 and
  CXCR4 was seen, 1.2).
  - Clock: ARNTL, CLOCK, NPAS2, PER1, PER2, PER3, CRY1, CRY2, NR1D1, NR1D2, DBP, TEF, HLF, CIART, BHLHE40, BHLHE41.
  - Glucocorticoid-induced in lymphocytes: FKBP5, TSC22D3, ZBTB16, KLF9, DDIT4, PER1, SGK1, DUSP1, TXNIP, CXCR4.
  - In U_A: ARNTL, CLOCK, PER1, PER2, CRY1, CRY2, NR1D2, DBP, FKBP5, ZBTB16, KLF9, DUSP1, TXNIP, CXCR4.

**2.4 Cells (S5b).**
- Kept T cells of the 25 samples: `PB$cells` with `kept` and `group4 == "T"`, 14,609 cells.
- Counts: the RNA `counts` layer of `02_seurat_objects/04_annotated/GSE185381/<sample>.rds`, on the same 36,601-gene axis.
  Checks as in G01's script 40: the kept-T sums equal `PB$raw$T` exactly.
- Fine labels: the nine BoneMarrowMap T labels.
- Doublet score: `ws_a1_batch/w4_doublet_scores.csv.gz`. It is scDblFinder's pre-training heuristic, computed per
  library and untrained.

## 3. Analyses

### 3.1 Gene-set description (the two interpreted families)

- **Test.**
  ```
  edgeR::camera.DGEList(y_A, index, design, contrast = "statusAML",
                        inter.gene.cor = NA, allow.neg.cor = FALSE, sort = FALSE)
  ```
  `y_A` is S5's T DGEList with its dispersions, rows restricted to U_A. `inter.gene.cor = NA` estimates each set's
  correlation from the residual space of the fit, which is orthogonal to status.
- **Families.**
  - P-H: the 45 HALLMARK sets.
  - P-T: the 22 Gavish sets (if 9.4).
  - BH within each family. **Flagged** = FDR < 0.05. Direction = camera's `Direction`: "Up" means a higher AML-versus-
    healthy statistic than the average U_A gene.
- **Reported for every set:**
  - NGenes, the estimated correlation, Direction, PValue, FDR;
  - the set shift = mean logFC of the set's genes minus mean logFC of the other U_A genes (log2, descriptive);
  - the handling fraction (2.3);
  - the set's P370 and F402 members, up and down, each with its expected count (members in U_A × arm size / 6,957).
- **Leading genes of flagged sets**: fgsea's leading edge (4.4), as a description.

### 3.2 Handling: adjacency and removal (for every flagged set)

1. **Adjacency.** The flag from 2.3, fixed before the run.
2. **Removal.** Rerun 3.1, with the same settings, after removing every H gene from U_A (6,259 genes left). H genes leave
   both the sets and the background. BH within each family runs over the sets that keep ≥ 15 genes. Reading for each
   flagged set:
   - **"not carried by the listed handling genes"**: still FDR < 0.05, same direction, and the set shift is at least
     half of the primary's;
   - **"not separable from the listed handling genes"**: fewer than 15 genes remain;
   - **"carried by the listed handling genes"**: otherwise.
3. **Direction check** (descriptive). Take the set's Savage members that are in U_A. Count the share whose AML-versus-
   healthy sign equals their delay direction:
   - ≥ 0.75: "runs with the delay direction";
   - ≤ 0.25: "runs against it";
   - otherwise: "mixed".

### 3.3 Compartments (read only for flagged sets)

**(a) Each compartment alone (reported, not read).** camera as in 3.1, on each compartment's refitted S5 fit (B, NK,
Mono). Universe: that compartment's tested genes minus the extended ranking. Sets: P-H and P-T with ≥ 15 genes there. BH
within each compartment and family. Monocytes are reported only, never read, because they are often clonal in AML
(A1-G01 5.1).

**(b) T versus X, for X = B and X = NK (read).**
- **Donors**: those with both T and X pseudobulks (B 23, NK 21).
- **Genes**: tested in both S5 fits, minus the extended ranking.
- **Model**: one edgeR model on the stacked raw counts, with TMM and `estimateDisp(robust = TRUE)` on the stacked design.
  - Columns:
    - one fixed effect per donor;
    - `cT`, the T-compartment indicator;
    - `cT:AML`;
    - `cT:hb`, `cT:gran`, `cT:sexmale`;
    - `cT:CD8share` and `cT:Xshare` (X's own `subtype_share`).
  - Donor effects absorb every donor-level main effect, status included. `cT:AML` estimates the AML-versus-healthy shift
    in T minus the same shift in X.
- **Test**: camera, estimated correlation, contrast `cT:AML`, BH within each X and family over the sets with ≥ 15 genes.
- **Reading for each flagged set:**
  - "larger in T than in X": the interaction has FDR < 0.05 in the set's direction;
  - "larger in X than in T": FDR < 0.05 in the opposite direction;
  - "shared with X": otherwise, if the ratio ρ_X = set shift in X / set shift in T (each from its own S5 fit, on the
    common genes) is ≥ 0.5;
  - "not resolved": otherwise;
  - "not testable in X": fewer than 15 genes in common.
- **Never written: "T-specific"** (A1-BATCH 12.5).
- **Calibration**: K5 (section 6) plants T-only, shared and half-shared effects in the real stacked counts. If K5 fails,
  every reading of (b) is "not readable".

### 3.4 Gene identity in healthy marrow T subtypes (descriptive)

- **Reference.** The 10 GSE185381 healthy donors. Labels are tagged if at least 5 donors have ≥ 20 kept cells in them:

  | label | donors |
  |---|---|
  | CD4 Naive | 8 |
  | CD4 Central Memory | 7 |
  | CD4 Effector Memory | 6 |
  | CD8 Naive | 5 |
  | CD8 Central Memory | 6 |
  | CD8 Effector Memory 1 | 7 |
  | CD8 Effector Memory 2 | 9 |

  CD4 Regulatory (1 donor) and CD8 Tissue Resident Memory (3) are not tagged.
- **Statistic.** For gene g in U_T, label k and eligible donor d:

  L(g, k, d) = log2((CPM of g in label k of donor d + 1) / (CPM of g in donor d's other kept T cells + 1)).

  CPM = 10^6 × count / column total. Pairing within the donor cancels donor-level shifts. Over donors: M(g, k) is the
  median of L, and c(g, k) is the share of donors with L > 0.
- **Tag.** "k-high" if M ≥ 1, c ≥ 0.8 and the median CPM in k ≥ 10. Classes:
  - **naive**: CD4 or CD8 Naive high;
  - **central memory**: CD4 or CD8 Central Memory high;
  - **effector memory**: CD4 Effector Memory, CD8 Effector Memory 1 or 2 high;
  - none.

  A gene can belong to more than one class; overlaps are reported.
- **Readout.**
  - For the P370 up and down arms, the F402 arms and G02's 25-gene core: the count and share per class.
  - Base rate: the share in U_T.
  - Matched reference: 1,000 random sets of the same size drawn from U_T within 10 bins of status-free pooled healthy T
    logCPM (seed SEED + 704); the percentile of the observed share is reported.
  - No p value is interpreted.
- **Sensitivities.**
  - Threshold 0.585 (1.5-fold).
  - Oetjen (GSE120221): 20 donors, libraries pooled per donor as in S9, genes mapped by symbol, the same label rule.
    These donors are out of sample for the contrast, but in-sample for the BoneMarrowMap labels.
  - Without Control1-5: 5 donors; labels with ≥ 4 eligible donors; c ≥ 0.75.

### 3.5 S5b, the single-cell layer (dated amendment of S5 4.2; descriptive)

**(a) Scores.**
- `UCell::ScoreSignatures_UCell` (2.10.1), maxRank 1,500, ties "average".
- Run on the full 36,601-gene count matrix of the kept T cells.
- The up set and the down set are scored separately: P370, with F402 as a sensitivity.

**(b) Pipelines.**
- **P1 (as registered)**: raw counts.
- **P3 (reading)**: every cell above 1,500 UMI is thinned binomially to 1,500 UMI (seed SEED + 706); cells at or below
  1,500 UMI are left unchanged. Then as P1.

**(c) Reference.**
- For each fine label: the mean, SD and P95 of the healthy cells' scores, and z = (score − mean_k) / SD_k.
- AML patients are compared with all 10 healthy donors. Each healthy donor is compared with the other nine
  (leave-one-donor-out).

**(d) Registered summaries, P1 and P3.** Per patient: the mean UCell, the mean z, and F95 = the fraction of cells above
their label's healthy P95. P1 is reported as registered, with the depth warning of 1.4.

**(e) Reading rule (P3, up set, per AML patient).**
1. **Eligibility**: ≥ 50 kept T cells. 13 of the 15 patients qualify; AML2451 (47) and AML0310 (20) are reported as
   "too few cells".
2. **Shift gate**: the patient passes if

   mean z − 1.96 × sqrt(var(z) / (v_H × n)) > h_max.

   v_H is the variance of the pooled healthy reference z. h_max is the largest of the 10 healthy leave-one-donor-out mean
   z values. A patient who does not pass gets "no detectable shift".
3. **Shape**.
   - f0 = the kernel density of the pooled healthy reference z: R `density()`, 2,048 points, from min − 3 to max + 8,
     floored at 1e-8, interpolated linearly.
   - Uniform model: z ~ f0(z − δ). δ is fitted by maximum likelihood: a grid from −3 to 20 in steps of 1, then
     `optimize` within ± 1.
   - Mixture model: z ~ (1 − π) f0(z) + π f0(z − Δ), with π ∈ {0.1, …, 1.0} and Δ ∈ {0, 0.5, …, 20}.
   - LR = max(0, ℓ_mixture − ℓ_uniform). π̂ is the maximising π, set to 1 if LR < 10^−6.
   - Calls:
     - **"sub-population"**: π̂ ≤ 0.5 and LR ≥ 10;
     - **"uniform"**: π̂ ≥ 0.9 or LR < 2;
     - **"not separable"**: otherwise.

   (D2's `13_sim.R` `pstats` and `17_final_rule.R` implement this; the md5s are in 8.4.)
4. **Cohort reading**:
   - "mostly uniform" or "mostly sub-population": at least 2/3 of the patients with a detectable shift share that call;
   - "mixed": otherwise;
   - "not readable": fewer than 6 patients pass the gate, or a K6 bar fails.
5. **Down set**: the same statistics, but no shape call. A 10-20 % sub-population cannot carry a down gene of
   |log2FC| > 0.5 (D2 3.2).

**(f) Per label (description).** CD4 Naive, CD4 Central Memory, CD8 Effector Memory 1 and 2 are the labels with ≥ 50
cells in at least 7 AML patients (D2 2.1: 7, 7, 7 and 11). For each (patient, label) pair with ≥ 50 cells: the mean z
and the shape call.

**(g) Sensitivities.** The healthy reference without Control1-5 (4.3); F402.

**(h) Where GNA15 and HMGB3 sit (description only).**
- Doublet-like stratum: kept T cells above their library's P90 of the W4 score, computed over all scored cells of that
  library.
- For each of the two genes: the share of its UMIs in AML kept T cells that come from the stratum. It is compared with
  1,000 random U_T genes matched on pooled mean and detection (10 × 10 bins; seed SEED + 708). "Concentrated in
  doublet-like cells" if the gene lies above the random P95.
- GNA15 stays out of the biological reading whatever this shows (A1-BATCH D4).
- The set-level doublet readout proposed in D2 is not run, because it is uncalibrated.

**(i) Unit rule.** No per-cell p value is computed. Cells enter only through per-patient summaries (S5 P3).

### 3.6 Gene classes and the overlap table (descriptive)

**(a) Gene classes removed from the cross-fitted score.**
- Method: G02's stored fold lists (`ws_a1_g02_fold_genes.csv`), scored as A1-BATCH's post-hoc script does: log-CPM over
  the 25, out-of-fold z, mean z(up) − mean z(down). gap_oof = mean AML score − mean healthy score.
- Rows:

  | genes removed | gap_oof |
  |---|---|
  | none | known: 1.3782 |
  | Savage | known: 1.3184 |
  | Savage + IEG | known: 1.3228 |
  | lncRNA (2.3) | |
  | clock + glucocorticoid list (2.3) | |
  | lncRNA + clock/glucocorticoid + H | |

- Also reported: the lncRNA share in the P370 and F402 arms against the U_T base rate (305 / 7,379), with expected
  counts.

**(b) Overlap table.**
- Lists on disk: Savage up and down; IEG; the 50 HALLMARK sets; the 22 Gavish sets; the cNMF normal MP10 (naive) and MP4
  (stress) programs, labelled "derived from this project's cells".
- Downloaded lists, only if 9.5 is approved and the files are on disk before the run: Lasry 2023 (246-gene adult
  signature; 38-gene iScore), Rutella 2022 (IED172, IED68, PI24), Mazziotta 2024 (top-100 and, if printed, the 48 genes).
- Columns:
  - list and size;
  - genes in U_T and in U_A;
  - overlap with the P370 up and down arms, the F402 arms and G02's core;
  - expected counts (list genes in U_T × arm size / 7,379);
  - the names of the overlapping genes.
- No p values.
- **Public repo**: gene names, counts and citations only. A copied full list is never committed.
- Lasry's lists are malignant-cell signatures; their row is positioning only (C 2.2).

## 4. Sensitivities (each reported; none changes a reading)

4.1 Universe U_T (all 7,379 genes) for 3.1-3.3.
4.2 Control4 coded female for 3.1-3.2. The level is recoded with `factor("female", levels(sex))`, asserting levels
    `female` and `male`.
4.3 Without Control1-5 (15 vs 5) for 3.1, for 3.4 (rule in 3.4) and for S5b's healthy reference.
4.4 Other views of the same contrast:
    - camera with `inter.gene.cor = 0.01`;
    - `fgseaMultilevel` preranked on sign(logFC) × sqrt(F), with F from S5's QL refit (minSize 15, maxSize 500, eps 0,
      seed SEED + 705);
    - fry directional.
4.5 F402 for 3.4-3.6.

Reading: for each flagged set, the main text names, in one parenthesis, any sensitivity under which the set loses
FDR < 0.05 or flips direction. A set that is significant only in 4.4 is "not detected by the primary".

## 5. Multiplicity and what is read

| family | analysis | tests | FDR | status |
|---|---|---|---|---|
| P-H | 3.1, HALLMARK | 45 | BH within | interpreted |
| P-T | 3.1, Gavish T programs (if 9.4) | 22 | BH within | interpreted |
| removal | 3.2 | ≤ 45 + ≤ 22 | BH within | read only for flagged sets |
| X-B, X-NK | 3.3(b) | ≤ 41 / ≤ 43 HALLMARK, plus Gavish | BH within | read only for flagged sets |
| per compartment | 3.3(a), B, NK, Mono | each | BH within | reported |
| sensitivities | 4.1-4.4 | each | BH within | reported |

**Rule.** Nothing outside P-H and P-T is interpreted on its own. No secondary family can create a flagged set; it can only
qualify one (handling, compartment, stability). BH runs within each family, never pooled. The families are views of one
contrast, so there is no correction across families.

## 6. Known answers (run first; no readout is printed until all pass)

There are three states, as in G04:
- **PASS**;
- **FAIL deterministic**: a code bug. It is fixed, a dated note is written, and the run is repeated;
- **FAIL statistical**: the reading named in the row becomes "not readable" or "not interpretable". The other readouts
  still run, labelled descriptive.

| id | check | bar | type |
|---|---|---|---|
| K0 | md5s of the frozen inputs: `ws_pb_GSE185381.rds`, `ws_s5_scores.csv`, `ws_s5_de_{T,B,NK,Mono}.csv`, `ws_s5_imprint_genes.csv`, `ws_s8_t_imprint_ambient_prone.csv`, the ambient table, the chrY list, the HALLMARK GMT, the pinned Gavish GMT, the Savage set, `ws_a1_g02_fold_genes.csv`, `ws_a1_g02_oof_scores.csv`, the W4 scores, the GTF; the IEG list equals `WS_B0_IEG` | equal to recorded values | deterministic |
| K1 | S5 T refit: genes tested, hits (up / down); the hit set equals the frozen list; max \|logFC − `ws_s5_de_T.csv`\|; `colnames(design)[2] == "statusAML"` | 7,379; 402 (171 / 231); TRUE; < 1e-9; TRUE | deterministic |
| K1b | B, NK, Mono refits: genes tested; max \|logFC diff\| against the stored tables | 4,292 / 4,151 / 9,252; < 1e-9 | deterministic |
| K2a | camera calibration on the real T fit: 500 random sets from U_A, sizes drawn from the 45 HALLMARK sizes, matched on logCPM (10 bins; seed SEED + 701); share with p < 0.05 | ≤ 0.10 | statistical (3.1 "not interpretable") |
| K2b | plants in real counts (seed SEED + 702): 10 random 50-gene sets with y′ = y + Poisson(y(2^0.6 − 1)) in AML samples; 10 with binomial thinning at 2^−0.6; refit with S5's pipeline, then camera over the 45 HALLMARK sets plus the plant | planted set flagged in its direction in ≥ 8 of 10 per direction; never in the wrong direction | statistical (stop and inspect before any family is read) |
| K3 | U_A without H has 6,259 genes; every HALLMARK set keeps ≥ 15 genes; the 12 + 20 adjacency flags equal 2.3 | as stated | deterministic |
| K4 | healthy-subtype tags, positive control: among CCR7, LEF1, TCF7, SELL present in U_T, and among GZMB, GNLY, NKG7, PRF1, GZMH present in U_T | ≥ 3 present per group; ≥ 75 % tagged naive (first group) and effector memory (second); none tagged in the opposite class | statistical (3.4 "not readable") |
| K5 | cross-compartment plants in the real stacked counts (seed SEED + 703); 20 random 50-gene sets per type from the common universe: (i) T-only, 0.6 log2 in AML T samples; (ii) shared, 0.6 in AML T and AML X; (iii) half-shared, 0.6 in T and 0.3 in X; the reading of 3.3(b) applied to every planted set | (i) "larger in T" in ≥ 14/20; (ii) "larger in T" in ≤ 2/20 and "shared" in ≥ 14/20; (iii) reported | statistical (3.3(b) "not readable"); if (iii) gives "larger in T" in more than 10/20, the text writes "larger in T's measurement" |
| K6a | UCell engine: one healthy donor's scores computed twice with the same thinning seed | identical; thinned cells sum to 1,500; cells at or below 1,500 untouched | deterministic |
| K6b | planted shape on real healthy cells (seed SEED + 707). Each healthy donor in turn is held out as a pseudo-patient. Plants on the real up-set genes, with synthetic per-gene log2FC β_g = 0.6 × min(2.5, 0.5 + Exp(mean 0.4)); S5's logFCs are never used. (a) uniform: y′ = y + Poisson(μ̂_g(2^β_g − 1)) in every cell; (b) 20 %: a random fifth of cells get y′ = y + Poisson(μ̂_g(2^β_g − 1)/0.2). μ̂_g = cell UMI × the pooled profile of the cell's label in the other nine donors. Down set thinned at 2^−β_g in (a) only. Then the exact 3.5 pipeline: reference and h_max from the other nine. | (a) "uniform" or "not separable" in ≥ 9/10, "uniform" in ≥ 7/10, never "sub-population"; (b) "sub-population" in ≥ 8/10; unplanted held-out donors "no detectable shift" in ≥ 9/10 | statistical (S5b cohort reading "not readable"; P1 still reported) |
| K6c | matching: in K6b, the median over donors of \|up-set pseudobulk log2FC (a) − (b)\| | ≤ 0.05 | statistical (as K6b) |
| K6d | AML depth: each patient's median UMI before and after thinning | reported; the AML-depth null cannot be planted on shallower healthy cells, so D2's synthetic null is cited (false shift ≤ 0.025) | reported |
| K7 | gene-class removal: rebuilt out-of-fold scores equal `ws_a1_g02_oof_scores.csv`; gap_oof with no removal, without Savage, without Savage + IEG | max \|diff\| < 1e-9; 1.3782; 1.3184; 1.3228 (4 decimals) | deterministic |
| K8 | per-cell inputs: kept-T sums equal `PB$raw$T`; label sums equal `PB$fine`; 14,609 cells; lncRNA count in U_T | exact; 305 | deterministic |

Run-time blindness: K2, K5 and K6 read real counts after the tag. K6 uses healthy cells only. No AML-versus-healthy set
or per-cell value is printed until every check in a mode has a state.

## 7. Reading rules and fixed sentences

**7.1 Sets (3.1).** Template:

> "In the S5 T-cell contrast (15 adult AML vs 10 healthy donors, GSE185381), set X had a higher [lower] average
> AML-versus-healthy statistic than the average tested gene (camera with estimated inter-gene correlation, FDR q; same
> contrast that defined the axis; descriptive)."

Then the handling reading (7.3) and the compartment reading (7.6).

**If no set is flagged in either family:**

> "No HALLMARK or T-cell-program set differed from the average tested gene at FDR 0.05 with a test that accounts for
> inter-gene correlation (n = 25). In simulation, with a genome-wide donor component present, this test detected
> whole-set shifts of 0.3 log2 in about a third of cases at low within-set correlation and in 5 % at higher correlation,
> so this is not evidence that the axis lacks set structure."

**7.2 Secondary views.** A set flagged only by camera 0.01, fgsea or fry is "not detected by the primary". fry Mixed
results never appear in the text.

**7.3 Handling.**
- Each flagged, handling-adjacent set carries its 3.2 reading and its direction check.
- "Not carried by the listed handling genes" is never shortened to "not a handling effect".
- **Consequence rule (fixed now).** Take each direction separately, over P-H and P-T together. If at least one set is
  flagged and at least half of the flagged sets are "carried by" or "not separable from" the listed handling genes, the
  paper states as a main limitation:

  > "The set-level description of the [up/down] arm rests on genes that also change with blood-processing delay or
  > immediate-early stress. Without those genes it is [not flagged / not testable]. Lists from processing studies share
  > most of their T-cell genes with differentiation programs (20 of 22 T-cell programs here), so a handling reading and a
  > T-cell-state reading of this arm cannot be separated with these data."

**7.4 Gene identity and composition (3.4).**
- Template:

  > "k of the [up/down] genes are [naive/central-memory/effector-memory]-high in healthy marrow (paired within donor;
  > base rate b; matched-random percentile p)."
- Always followed by:

  > "A1-G01 attributes about 5 % of the gap to T-subset composition (C 0.045), so their presence does not by itself
  > indicate a composition shift."
- Tags are never used to re-estimate composition, or to claim a shift "inside subtypes". G01's S_w is the only statement
  on that.
- An effector-like description adds:

  > "CMV serostatus, which drives effector expansion, is unrecorded in both arms."

**7.5 S5b.**
- Template:

  > "With depth-equalised per-cell scores, the up-set shift was [carried by most cells / concentrated in a minority of
  > cells / mixed across patients] in k of 13 patients with ≥ 50 kept T cells (patient-level summaries)."
- "Not readable" and "not separable" are reported as such.
- A sub-population call lists what it can be: composition below the label level, a clonal expansion, label error, cells
  rich in ambient RNA, or doublets. It is not "a new cell state", and it does not change G01's category.

**7.6 Compartments (3.3).**
- Allowed readings: "larger in T than in B [NK]", "shared with B [NK] (ratio r)", "larger in B [NK] than in T", "not
  resolved", "not testable in B [NK]".
- "Shared" is read as "not restricted to T". It is not read as "systemic biology", because ambient RNA and handling are
  shared by every compartment of a sample.

**7.7 Context sentences that travel with any G07 result** (verbatim from GATES.md and the outcome sections).
- Not T-specific: the direction is shared with B and NK of the same donors, and appears out of sample in B, NK and
  monocytes.
- Processing: "None of the six listed within-study processing factors, each tested one at a time, explained or weakened
  the frozen direction … Sample source, storage, sorting method and chemistry version are constant within arms or
  unrecorded, and were not tested."
- Delay: it overlaps a blood-processing-delay program in the same direction (D8).
- Composition: "beyond T-subset composition" (C 0.045), which does not mean caused by AML.
- Age: "AML, with an age share bounded below one third in healthy marrow (one protocol)."
- Direction: "shown in one study" (A1-G04).

**7.8 Never written, under any outcome.**
- "pathway X is activated, induced or suppressed in AML T cells"; "drives", "reflects signalling", "caused by AML";
- "independent evidence", "validates", "confirms the imprint";
- "T-cell-specific program" or "T-cell-specific state";
- "not a handling effect", "not an artefact", "robust to processing";
- "the axis is not a pathway" (from a null result);
- "exhausted" or "senescent", unless a flagged set is named, and then only in the 7.1 form;
- "matches no published signature" as a novelty claim;
- "opposite to published ISG-up reports";
- any per-cell p value; "a new cell state";
- any biological reading of GNA15 or of the 32 ambient-prone genes;
- "first to compare";
- any statement that uses PROGENy, CytoSig or ligand-receptor evidence;
- the in-sample 3.84 SD or AUC 1.0 as an effect size.

## 8. Code, outputs, compute (plan; the scripts are written after the tag)

**8.1 Scripts.**
- `scripts/ws/43_a1_g07.R`, modes:
  - `--known` (K0-K5, K7);
  - `--sets` (3.1, 3.2, 4);
  - `--cross` (3.3);
  - `--tags` (3.4);
  - `--classes` (3.6(a));
  - `--overlap` (3.6(b));
  - `--combine` (summary and GATES block).
- `scripts/ws/44_a1_g07_cells.R`, modes:
  - `--known` (K6, K8);
  - `--score` (3.5(a)-(g));
  - `--doublet` (3.5(h)).
- `scripts/ws/submit_ws_a1_g07.sbatch`: dated read-only copies of the script, `OMP_PROC_BIND=false`, smoke switch.
- Pattern as in G03: `script_md5` in every table; a write guard on real outputs; inputs re-checked against `--known`'s
  md5s.

**8.2 Outputs.**
- `results/tables/ws/ws_a1_g07_{checks, sets, sets_removal, sets_sens, compartments, cross, tags, classes, overlap,
  s5b_patients, s5b_labels, s5b_checks, doublet_genes, summary}.csv`.
- GATES.md block `A1-G07`, with a "Reading" line, not a "Verdict".
- `/LARGE1/gr10634/gaozy/aml_niche_net/ws_a1_g07/`: per-cell scores, code copies, smoke outputs.

**8.3 Seeds** (SEED = 491638; offsets +701 to +720 are free, B 8.2):

| offset | stream |
|---|---|
| +701 | K2a random sets |
| +702 | K2b plants |
| +703 | K5 plants |
| +704 | 3.4 matched reference |
| +705 | fgsea |
| +706 | S5b thinning |
| +707 | K6 plants |
| +708 | 3.5(h) random genes |
| +709 | smoke counts |
| +710 | smoke cells |
| +711-720 | reserved |

**8.4 Simulation provenance.** The rules in 3.1, 3.3 and 3.5 rest on the D1 and D2 studies. Their scripts, with md5s:

| script | md5 |
|---|---|
| D1 `02_sim_geneset.R` | `94730713…` |
| D1 `03_sim_crosscomp.R` | `30079a4d…` |
| D1 `04_summarise.R` | `a2159db9…` |
| D1 `05_crosscomp_readout.R` | `76028d73…` |
| D2 `10_engine.R` | `1157306311…` |
| D2 `13_sim.R` | `bf631cdb…` |
| D2 `14_summarise.R` | `28c46716…` |
| D2 `16_rule_table.R` | `61e976a8…` |
| D2 `17_final_rule.R` | `e52077fd…` |

Decision 9.8: copy them into the repo with the code commit.

**8.5 Smoke runs.**
- Synthetic counts and synthetic cells on real IDs, with planted effects.
- No real covariate or status is written next to a smoke statistic (memory note on smoke leaks).
- Outputs go under `/LARGE1/.../ws_a1_g07/smoke/` and are deleted after the real run reproduces.

**8.6 Compute.** Under 3 hours in one `c=16, m=64G` job (B 9.2). K2b and K5 need about 60 refits at about 10 s each; the
UCell passes take minutes.

**8.7 After the run.** Outcome section 11; GATES block; handoff S15 note; memory.

## 9. Decisions for the user (recommended option first)

1. **Scope and timing.** Write and tag now, with 3.1-3.6 as specified; downloads stay optional (9.5). Alternative: wait
   for the late-November decision point.
2. **Primary test.** camera with estimated inter-gene correlation, FDR 0.05 (3.1). Alternatives: camera 0.01 (more
   power, 8-19 % false set calls in simulation); fry as co-primary.
3. **Footprints and TF activity.** Do not use PROGENy, CytoSig or CollecTRI (departure 4). Alternative: record that a
   descriptive use on T pseudobulks lies outside GATE_B0's verdict, with F402 genes removed from every footprint.
4. **T-cell programs.** Pin the 22 Gavish 3CA T-cell sets from the repo's MSigDB v2026.1 file and use them as the second
   interpreted family (P-T). Alternatives: descriptive only; not used. Without them, only 3.4 speaks to T-cell state,
   because HALLMARK says nothing about naive, effector or dysfunction programs.
5. **Downloads, one batch.** Lasry 2023 supplementary tables (xlsx, 4.8 MB), Rutella 2022 supplementary tables (xlsx,
   191 KB, CC BY 4.0), Mazziotta 2024 supplement (PDF, 2.3 MB; text-and-data-mining reserved, kept outside the repo).
   - Each gets a PROVENANCE note under `/LARGE1/gr10634/gaozy/reference/`.
   - They enter 3.6(b) only. Lasry's lists also serve G08.
   - Options: approve now; at the decision point; not at all.
6. **S5b.** Run it as the dated amendment in 3.5 (P1 as registered plus P3 and the shape rule). Alternatives: P1 only;
   disclose it as never run.
7. **Codings.** Metadata sex primary; Control4 female and "without Control1-5" as sensitivities; T-ALL flag on every
   table (2.1, 4).
8. **Simulation code.** Copy the nine D1/D2 scripts into the repo (`scripts/ws/g07_sim/`) with the code commit.
   Alternative: md5s only (8.4).
9. **ADT (option c).** Defer: the panels coincide with the capture period, and the extended markers cover 2 healthy
   donors.
10. **G10(a) non-AML comparators.** A limitation now; revisit at the decision point.
11. **The name "imprint".** Decide when writing (G11). A new name containing "T-cell state" would conflict with "not
    T-specific".
12. **Time stamp.** An internal git tag, worded "internally pre-registered (git-tagged)", as for every earlier step.

The real choices are 9.3, 9.4, 9.5 and 9.6. The others confirm defaults.

## 10. Review history

- 2026-11-09: inventories A, B, C, D1 and D2 plus a completeness critic (workflow wf_309ce460-130; outputs outside the
  repo). The critic checked 44 citations: 37 correct, 3 off by one or two lines, 4 overstated or wrong. Those four are
  corrected here:
  - the median-UMI covariate did run (A1-BATCH W2);
  - UCell code exists outside `ws/`;
  - no dissociation list was scored;
  - the age label is S9's narrowed form.

  No blindness breach was found.
- DRAFT v1 (this file).
