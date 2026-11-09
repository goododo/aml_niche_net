# PREREGISTRATION A1-G07: what the T-cell axis is (a pre-registered description)

**Status: DRAFT v2, 2026-11-09, after review 1 (four lenses; section 10). Not approved; no tag.** Nothing declared in
sections 3-6 has been computed. Section 2 counts use only gene names, set memberships and cell numbers (1.4).

**Sources.** All are kept outside the repo, under `/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/`.
- Read-only inventories in `inventory/`: `A_obligations.md` (A), `B_machinery.md` (B), `C_resources_prior_art.md` (C),
  `D1_geneset_design.md` (D1, synthetic gene-set study), `D2_singlecell_layer.md` (D2, synthetic single-cell study),
  and `critic.md`.
- Review 1 in `review1/`: L1 gene-set statistics, L2 single-cell layer, L3 conformance, L4 interpretation.
- Gene-name counts behind this draft in `drafting/`: `membership_counts.*`, `removal_counts.*`, `detection_floor.*`,
  `savage_celltypes.*`.

**Place in the sequence.** This is the last item of handoff v2.0 §17 S15. The steps before it:
- A1-G02: PASS.
- A1-BATCH: PASS with qualifications.
- A1-G01: "beyond T-subset composition, detected inside the main labels".
- A1-G04: (a) "not shown"; (b) "not readable".
- A1-G03: "held through R6".

The late-November decision point (assessment §5 item 6) has not been held. G07 is written now for two reasons:
- every interpreted part uses resources already on disk;
- downloaded lists can enter only a descriptive table (3.6(b), decision 9.5).

## 0. Question, scope, and what each outcome lets the paper say

The frozen T axis is a gene list today:
- **P370**: 156 up, 214 down. The text quotes P370 (2.1).
- **F402**: the 402 S5 hits before the post-hoc removal of 32 ambient-prone genes. It appears as numbers in tables only.

**G07 asks: what kind of program is it?** It answers at four levels, plus tables:
1. **Gene sets (3.1-3.2).**
   - Which pathway sets and T-cell-program sets carry more of the S5 AML-versus-healthy T contrast than the average
     tested gene?
   - For each such set: does it still hold without the genes that change with blood-processing delay in the same
     direction, or with immediate-early stress?
2. **Compartments (3.3).** For each such set: is the shift larger in T than in B or NK cells of the same donors, or
   shared? Does each compartment follow its own delay list?
3. **Gene identity in healthy marrow (3.4).** Which axis genes mark naive, central-memory or effector-memory T cells in
   two healthy references? How much of that is expected from composition alone?
4. **Cells (3.5, the S5b layer, T only).** Is the up-set shift carried by most T cells or by a minority?
5. **Tables (3.6).** The cross-fitted score with gene classes removed; overlaps with other gene lists.

**What G07 is not.**
- **Not a gate.** No result changes the existence, direction, composition, batch or bridge results (A1-G02, A1-BATCH,
  A1-G01, A1-G04, A1-G03).
- **Not independent evidence.** Every set test re-describes the S5 contrast that selected the genes: same donors, same
  design, same noise. The text says so every time (7.1).
- **Not able to settle cause.** It describes T versus B/NK sharing and alignment with handling lists. It cannot settle
  T-specificity, handling, AML causation, age or ambient RNA (7.7).
- **Not blind at the level of gene identity (1.2).**
  - Already seen: the HALLMARK membership counts of the 370, the Savage-overlap counts and their direction agreement,
    and most gene classes.
  - Not yet seen: any test on the full ranked statistic, any per-cell score, overlaps with full published lists other
    than Savage 2021, the healthy-subtype identity of the 370, and set tests in B, NK or monocytes.

**What each outcome lets the paper say.** The fixed sentences are in section 7.
- **Sets flagged (3.1).** "In the S5 T-cell contrast, set X had a higher (lower) average AML-versus-healthy statistic
  than the average tested gene." Every flagged set carries three things:
  - its handling reading (3.2);
  - its compartment reading (3.3);
  - the clause "(same contrast that defined the axis; descriptive)".
- **No set flagged.** "No set differed from the average tested gene at FDR 0.05 in a competitive test that accounts for
  inter-gene correlation (n = 25)." It is quoted with the real-fit power (K2c) and is never read as "the axis is not a
  pathway".
- **Handling-bound (7.3).** If at least half of the flagged units in one direction are bound to the listed handling
  genes, the paper states a main limitation: the set-level description of that arm cannot be separated from those genes.
  An arm that is not bound is still only associated with AML status in one study.
- **S5b (3.5).** One of "mostly uniform", "mostly sub-population", "mixed", "mostly not separable" or "not readable", with
  the four per-patient counts. "Sub-population" means "carried by a minority of cells", never "a new cell state".
- **Gene identity (3.4) and the tables (3.6)** are descriptive under every outcome.

### 0.1 What G07 owes, where it is handled, and declared departures

| source | obligation | here |
|---|---|---|
| addendum 1.6; verify_G07 item 9 | disclose the 10-28 look at the HALLMARK membership counts | 1.2 |
| verify_G07 best option; objections B1-B2 | the description predicted from gene names | 1.3, as a prior expectation |
| addendum 1.6; verify_G07 item 10 | camera/fry primary on the S5 design; fgsea secondary | 3.1 (departure 1); fgsea and fry reported (4.4) |
| verify_G07 item 10 | universe = tested genes minus extended-ambient genes | 2.2; all 7,379 as a sensitivity (4.1) |
| verify_G07 item 10 | report the 370 and the 402 | 2.1: P370 in the text, F402 as numbers in every table |
| addendum 1.6; verify_G07 item 7 | PROGENy/CytoSig only with a recorded user decision | not used (departure 4; 9.3) |
| verify_G07 items 8, 11 | durable overlap table from full lists, anchored on G02's core | 3.6(b) |
| verify_G07 best option | per-gene specificity across healthy T subtypes | 3.4, with two references (departure 8) |
| S5 4.2, :135; verify_G07 item 4 | S5b registered and never produced | 3.5: T only, as a dated amendment of S5 4.2 (departure 2; 9.6) |
| verify_G07 best option; verify_G06 | where the GNA15/HMGB3 signal sits | not run (departure 7) |
| verify_G06 | a score without lncRNA-like genes; one clock/glucocorticoid line | 3.6(a) |
| addendum 1.4, A1-G02 §7 (Control4); addendum 1.2, A1-G02 :96-100, A1-G03 disclosure 12 (T-ALL pools) | declare the sex coding and the T-ALL handling before any result | 2.1, 4.2, 4.3 |
| verify_G07 item 3 | "the vendor-versus-biobank question must be settled before the axis is interpreted" | cannot be met on disk; carried as a limitation (1.5 item 13; 7.7) |
| assessment G10(a) | non-AML comparators (systemic stress, CMV) | not run; carried as a limitation (1.5 item 12; 7.7) |
| A1-BATCH 12.5; A1-G01 11.5 | not T-specific; delay overlap; "beyond composition" is not "caused by AML" | 3.3; 7.6-7.7 |
| A1-G03 :466 | per-gene list status "stands", handed to G07 | used as is; claims are about sets (7.9) |
| A1-G02 :116-118 | the 25-gene always-selected core | 3.4, 3.6(b) |
| assessment 2.3 (:155-156) | generic labels may reopen the composition reading | 3.4; 7.4 |
| verify_G07 item 13 | GSE227903 diagnosis rows only | not used |
| verify_G08 items 4, 7 | overlap of the 370 with Lasry's lists; read the Author Correction first | 3.6(b), if on disk before the run (9.5); otherwise G08 |
| assessment G12 | "which genes are T-specific" for a bulk proxy | not built here; the direction is shared with B/NK/Mono (7.7) |

**Declared departures.**
1. **fry is not co-primary.** In D1's synthetic study, once a genome-wide donor component was present, fry "Mixed"
   flagged 71-89 % of sets that carried no effect. Only camera with estimated inter-gene correlation held its error rate
   in all three scenarios (1.4). fry directional is reported (4.4). D1's preset-correlation variant r* is not used.
2. **S5b: T only, with a reading pipeline and a rule.**
   - S5 4.2 registered S5b "per patient and compartment". It runs for T only. The NK and Mono parts are not run: S5's
     claim was narrowed to T (1.5 item 2), and both compartments carry known artefacts (the NK imprint vanished with the
     multiplexed covariate; Mono is clone-mixed). They are disclosed as an omission.
   - In D2's synthetic study, the two registered summaries could not separate a uniform shift from a sub-population, and
     AML's deeper cells moved them with no expression change (1.4).
   - They are reported as registered (pipeline S5b-raw). The question is read with a depth-equalised pipeline (S5b-1000),
     a detection floor and a shape rule (3.5).
   - The amendment is recorded here because S5 is tagged; S5's text is not edited.
3. **S5b runs after G01.** The verifier placed S5b before G01's result was interpreted. G01's verdict stands whatever S5b
   shows.
4. **No footprints and no TF activity** (if 9.3 as recommended). GATE_W1.1 and GATE_B0 failed in these data
   (GATES.md:24, 83). The paper outline lists the footprint control as a negative result (handoff :593).
5. **No dissociation list and no cryopreservation-stress list.** Neither is on disk. Freezing and cold-shock stress
   (e.g. CIRBP, RBM3) are not described; this is a limitation. The handling lists are Savage 2021 plus the project IEG
   list (2.3).
6. **Over-representation on the 370 is not a test.** Membership counts are reported with expected counts and no p value
   (they were already seen, 1.2).
7. **G06 leftovers, and two omissions.**
   - Only two description lines run (3.6(a)).
   - Disclosed omissions: de novo patients only, a 50-cell floor, the GATE_A1(4) re-read, a T-cell doublet
     drop-and-refit, and the capture-day arm of the T-ALL rule (3 healthy donors).
   - The median-UMI covariate already ran (A1-BATCH W2).
   - The per-cell GNA15/HMGB3 location is not run. The doublet-like stratum is concentrated in a few patients: 85 % of
     its 345 AML cells come from 3 patients, and 4 patients contribute none (review1/L2 t3). The doublet score is computed
     from the cells' own expression, so doublets and ambient-rich cells cannot be told apart. GNA15 already left the
     up-arm interpretation (A1-BATCH D4).
8. **Per-gene identity uses two references, by consensus.** The GSE185381 healthy donors are one arm of the contrast that
   selected the genes; Oetjen is out of sample. A class is stated only when both tag it (3.4).
9. **Handling removal is direction-matched** (3.2). The verifier asked for handling-artefact scoring; the set-level rule
   is new. The direction comes from A1-BATCH D8: the imprint overlaps the delay program in the same direction.

## 1. Already known before this draft (disclosed)

**1.1 Results G07 builds on** (GATES.md blocks and the outcome sections of the pre-registrations).

| step | result |
|---|---|
| S5 | T: 402 genes (171 up / 231 down), FDR < 0.10 and \|log2FC\| > 0.5. Fit: edgeR QL, `~ status + hb + gran + sex + CD8 share`, 7,379 genes tested. GATE_A1(2) failed (32 genes in the extended ambient ranking), giving P370 (S8). B 9 genes (not called), NK 91, Mono 275. GATE_A2 (coordination) failed. S5b not produced. |
| S8 | "GATE_S8 is read as T2 AND T1: PASS in GSE289435, consistent in GSE227903, not reproduced on Seq-Well GSE116256" (GATES.md:211-213) |
| S9 | f_age 0.102 (95 % CI −0.085 to 0.463); "bounded, not excluded"; older healthy donors do not score like AML |
| A1-G02 | cross-fitted gap 1.378 (1.98 out-of-fold healthy SD), p 0.008, AUC 0.93; 25 genes selected in all 25 folds (10 up, 15 down, all in P370) |
| A1-BATCH | PASS. Shared with B and NK in post-hoc relabellings and in G01 5.1; D10's declared reading is not usable. Overlaps the Savage delay program (D8: up 12 vs random 6; down 23 vs 10); delay score in MAIN gives R 0.82. IEG score +0.04 z (p 0.81). Post hoc: gap 1.318 without Savage genes; 1.3228 without Savage + IEG genes (`2026-11_a1_batch_posthoc/code_audit/a1_posthoc_d8_shared_genes.out`; unregistered) |
| A1-G01 | C 0.045 (−0.100 to 0.191); S_w 0.664 (0.497 to 0.772); out of sample, B 0.857, NK 0.893, Mono 0.732; 5.5 not fired; gamma-delta share 0.6 % (AML) vs 2.6 % (healthy); UMI per kept cell 1.15-2.07 times higher in AML in every label |
| A1-G04 | Chen2023 "not shown" (p 17/126; ρ_U 1.36); children "not readable"; the direction stays "shown in one study" |
| A1-G03 | held through R6; per-gene list "stands"; the count moves from 38 to 1,106 when one AML is left out; with Control4 coded female, 402 becomes 895 |

**1.2 Prior looks at the 370/402 against annotation** (ledger L1-L17 in A section 5; listed here because the
inventories are not in git).
- **L1, the S5 run (10-17).** 402 genes with logFC/FDR. The `ieg` flag marks IER2 only. Sensitivity arms:
  - IEG-removed 401;
  - multiplexed 405 (308 shared);
  - limma-voom 383 of 402;
  - subtype refits: CD4 Naive 85 hits, CD8 EM2 7 hits.
- **L2.** The identities and logFC of the 32 ambient-prone genes (15 up / 17 down; SELL, IL7R, FHIT, CD5, RUNX3, MAF,
  SRGN, CXCL8, ANXA2, EZR among them).
- **L3, S9.** Per-gene age effects of the 370 in Oetjen: 29 at raw p < 0.05, 0 at FDR 0.10; A2b alignment; the robust
  341 "record only".
- **L4, S10 review.**
  - The 370 self-identify better than about 94 % of expression-matched (about 82 % of composition-matched) random sets
    in healthy Griffin libraries.
  - chrX/Y and HLA-A/B membership (the 358-gene test-2 set).
- **L5, 10-28 conversation.** "No classic naive/effector marker in the 370; only 6 stress/IEG-list genes (IER2 down); 25
  of 214 down genes are lncRNA-like names". A 6-gene stress list was quoted and is in no table.
- **L6, 17 named-example lists (10-28).** The table is lost; the prose is in `inventory_prior_art.json` items[13].
  - Zero overlap: exhaustion, senescence, cytotoxic, Tpex, Knaus up, Lasry's named genes, Massoni-Badosa, Szabo
    tissue/proliferation/cytotoxic, Azimuth CD8 TEM / CD4 CTL / Treg.
  - Overlaps: Azimuth CD8 Naive (LINC02446, NOSIP, RGS10), CD4 Naive (NOSIP, LDHB), CD8 TCM (LDHB, LINC02446),
    gamma-delta/NK (TRDC), dnT (NUCB2), Knaus down (TMIGD2), Szabo IFN (IFIT3), Savage (DUSP4 up).
  - **Pinned HALLMARK membership counts:**

    | set | up | down |
    |---|---|---|
    | E2F targets | 7 | 3 |
    | UPR | 6 (ATF4, BAG3, EIF4E, MTHFD2, SHC1, ZBTB17) | 0 |
    | mTORC1 | 7 | 2 |
    | IFN-alpha | 0 | 4 (IFIT3, IFITM1, TRIM14, MVB12A) |
    | TNFA via NF-kB | 5 (incl. DUSP4, RELA) | 1 (IER2) |
    | OXPHOS | 3 | 6 |
  - Canonical markers absent from the 402: CCR7, TCF7, LEF1, GZMB, GZMK, PRF1, NKG7, GNLY, KLRG1, B3GAT1, CX3CR1, PDCD1,
    HAVCR2, TOX, FOXP3, CD69, HSPA1A, FOS, JUN, ISG15, MX1, STAT1.
- **L7, other prior-art items (10-28).**
  - DUSP4 +2.26 (FDR 0.040), CREM +1.53 (0.028), GMNN +1.51 (0.028), CXCR4 +1.30 (0.035), MTHFD2 +1.08 (0.046).
  - ZEB2 and DUSP4 shared with Vadakekolathu/Kazerani.
  - TSTR/heat-shock, IEG and tissue-residency genes (BAG3 up; none of HSPA1A/HSPA1B/DNAJB1/HSPH1).
  - Processing-delay partial resemblance; clock genes; FKBP5; HLA-A/B; 35 lncRNA-like genes.
  - Novelty judgements.
- **L8, objections inventory.**
  - Up-arm and down-arm classes, with logFC.
  - TRDC −1.86. CCR7, LEF1 and TCF7 "trend down but not significant" (they are K4's positive-control genes; K4 uses
    healthy donors only).
  - Clock and glucocorticoid genes; the lncRNA count (35 down vs 3 up); the ISG reading; mechanism hypotheses.
- **L9, evidence inventory.** Named genes with logFC; BAG3 and DNAJA3 (heat-shock-related); CXCR4, MTHFD2 and GMNN
  "broadly expressed outside T".
- **L10, the 10-28 synthesis.** "ISGs go down" (later corrected).
- **L11, verifiers.**
  - Flat: ISG15 +0.20, MX1 +0.09, IFI44L +0.01, IFI6 −0.16, STAT1 −0.56, IRF7 +0.07. Up: HLA-A +0.76, HLA-B +0.87.
  - Sub-threshold effector lean: GZMB +0.74, GZMH +0.79, FGFBP2 +0.99, all FDR > 0.10.
  - Flat: TOX −0.08, TIGIT −0.01, LAG3 +0.48, KLRG1 −0.54.
  - Clock genes CRY1, PER1 (only FKBP5 a hit). Per-channel ambient lists lack TRDC, DNTT, CD1A. Expected chance
    overlaps of 0.15-0.95 genes.
- **L12, addendum and assessment.**
  - Delay genes running against the delay direction: JUN −1.36, NFKBIA −0.87, IER2 −1.27. CD83 −0.32 and HIF1A −0.01
    are flat.
  - Running with it: ATF4, MTHFD2, BAG3, NDRG1, PLIN2 (up; "the reading is mixed").
  - GNA15 +2.27 and HMGB3 +1.59. Glucocorticoid program flat (TSC22D3, ZBTB16, KLF9). The clock line: DBP −1.18,
    FKBP5 −1.19, PER2 +1.05, ARNTL +0.89.
- **L13, A1-BATCH.**
  - Section-1 gene reads.
  - D8 Savage overlap counts (up 12, down 23, opposite 0).
  - **The fold-gene/Savage direction agreement, 0.994** (post hoc; `a1_posthoc_d8_shared_genes.out`).
  - D9 soup ratios (not interpreted); D4 GNA15/HMGB3.
  - Genome-wide T-versus-B/NK logFC correlations (r 0.60 / 0.70).
- **L14, A1-G01.** Named genes in its section 1; the 5.5 loading statistic; TRDC removed from the fold lists (gap
  1.24); AML4340's fold up genes broadly expressed (HLA-B, CREM, ARHGDIA, TANK, PIGT, MTRNR2L8).
- **L15, A1-G04.** Chen2023 gene-level sign concordance 0.832 (357 genes); halves; delay score.
- **L16, A1-G03.** P370 per rung (Sept-rule hits 73 / 1 / 1); the stored September table (326 genes, sign share
  1.000); AC007952.4; the S5 logFC of the 32 genes, seen by accident while drafting.
- **L17.** `signature_overlap_370.csv` is lost; its prose survives (L6).

**Not seen:**
- any camera, fry or fgsea run on the ranked statistic;
- any per-cell score;
- overlaps with full published lists other than Savage 2021;
- the healthy-subtype identity of the 370;
- set tests in B, NK or monocytes;
- the Gavish programs against the 370.

**1.3 Prior expectations** (stated now; they are not hypotheses).
- 2026-10-28 (verifier and objections inventory): "post-activation and biosynthetic-stress genes up; a quiescent,
  naive-like layer and lncRNAs down; a slight effector lean".
- 2026-11-09 (inventory C, from values already quoted): the same, plus "no exhaustion, no senescence, no ISG-up". The
  strongest competing reading is processing/handling.

A set-level result that matches these guesses is not a confirmation. The text quotes the result, not the guess.

**1.4 Seen while writing this draft** (no AML-versus-healthy expression value was computed).
- **Cell QC.**
  - Kept T cells: 14,609 in total. AML 6,911 in 15 patients (20-1,208 each); healthy 7,698 in 10 donors (123-1,917).
  - AML T cells are deeper: median UMI 2,952 vs 2,000.
  - Cells at ≤ 1,500 UMI: 27.6 % healthy vs 6.6 % AML (Control5 93.5 %).
  - Cells below 1,000 UMI: 7.1 % healthy vs 0.6 % AML.
  - W4 doublet score: median 0.095 vs 0.075; its Spearman correlation with log UMI is 0.62 in AML and 0.25 in healthy.
    The doublet-like stratum is concentrated (departure 7).
- **Membership counts (names only).**
  - Universes: 7,379 / 6,957; minus H_up 6,683; minus H_down 6,517; minus all of H 6,259.
  - 45 HALLMARK sets with ≥ 15 genes.
  - H has 698 genes in U_A (base rate 0.100).
  - **The delay list and the IEG list cover ≥ 25 % of the tested genes of 20 of the 22 Gavish T-cell programs. Only 9 of
    those 20 are naive, cytotoxic, regulatory or dysfunction programs.**
  - From the other side, only 178 of H's 698 genes are in any Gavish program.
  - Most naive-identity genes (LEF1, TCF7, CD27, PLAC8, TRAT1, SELL, IL7R) also fall with delay. So delay and naive loss
    cannot be separated by direction. Removal is still informative for NAIVE_2 and NAIVE_3.
  - After removing all of H, 7 Gavish sets drop below 15 genes: CD8 MEMORY_NAIVE_1 8, CD4 NAIVE_1 9, CD8 CELL_CYCLE 12,
    CD4 CELL_CYCLE 13, CD8 CYTOTOXIC 13, CD4 CYTOTOXIC 14, CD4 STRESS_HSP 14.
  - After direction-matched removal, every set keeps ≥ 10 genes (HALLMARK ≥ 16).
  - Ambient exclusion takes CCR7, SELL, IL7R, KLF2, CD69, LTB, JUNB, GZMB, GNLY and PRF1 out of the Gavish naive and
    cytotoxic programs.
- **Synthetic studies** (no project expression).
  - **D1, gene-set tests.** False-positive rate on sets without an effect, at p < 0.05 / FDR < 0.05:

    | method | no genome-wide component, ρ 0.05 | genome-wide, ρ 0.05 | genome-wide, ρ 0.20 |
    |---|---|---|---|
    | camera, estimated correlation | 0.041 / 0.004 | 0.024 / 0.001 | 0.011 / 0.000 |
    | camera, 0.01 | 0.173 / 0.084 | 0.093 / 0.020 | 0.291 / 0.189 |
    | fry, directional | 0.056 / 0.010 | 0.108 / 0.032 | 0.069 / 0.008 |
    | fry, Mixed | 0.116 / 0.044 | 0.917 / 0.892 | 0.771 / 0.714 |
    | fgsea, preranked | 0.235 / 0.175 | 0.167 / 0.091 | 0.377 / 0.310 |

    Power of camera with estimated correlation (FDR < 0.05) for a whole-set 0.30 log2 shift: 0.91, 0.32 and 0.05 across
    the three scenarios.
  - **L1, removal rule.**
    - A rule based on BH over the other sets called 58 of 103 uniform-effect sets "carried".
    - The own-p rule of 3.2 made 0 false "carried" calls and caught 75 of 78 handling-driven sets.
  - **L1, stacked T × X model.**
    - Full rank with the real covariates: B 46 × 30 (residual df 16), NK 42 × 28 (df 14).
    - Null sets: p < 0.05 in ≤ 0.056.
    - NK keeps 51 % of the information on `cT:AML`, B 72 %.
    - At 0.6 log2 the K5 code checks pass 20/20.
  - **D2 / L2, single-cell layer.** The scores are D2's synthetic up sets drawn from random tested genes.
    - S5b-raw shifted null AML patients by +0.45 healthy SD from depth alone.
    - Thinning to 1,500 UMI left a residual shift. In K6b replays it passed in 0.60 of runs in the base world; every
      failure was a false sub-population call in Control5.
    - Dropping cells below 1,000 UMI and thinning the rest to 1,000 passed 1.00 in the base world. It passed 0.05 when
      the set came from low-abundance genes.
    - Between-donor heterogeneity raised h_max (median 0.24 / 0.75 / 1.07) and cut power.
  - **Detection** (`drafting/detection_floor.out`). Random 156-gene tested sets detect a median of 8.9 genes per cell at
    1,000 UMI (P10 7.6) and 12.6 at 1,500. Low-abundance sets detect 2.7 and 4.1.

**1.5 Disclosures carried from earlier steps** (each enters the paper's limitations).
1. **The 32-gene removal was post hoc.** The 32 carry most recognisable lineage markers; they stay excluded from every
   statement (7.8).
2. **The narrowing to T cells chose where the axis is measured.** It was not a finding of T-specificity. B's
   "no imprint" (9 genes) is not evidence that B cells are unaffected.
3. **T3 is uninformative.** No CD4-versus-CD8 statement is made.
4. **The hashing covariate was mis-coded.** Gene membership belongs to one design's selection.
5. **Control4's sex.** Metadata male, cells female. Sample-level results barely move; the per-gene list moves a lot
   (402 → 895).
6. **Doublets.** The registered doublet arm never ran (replaced by A1-BATCH W4). decontX ran in A1-BATCH D4.
7. **T-ALL lanes.** Healthy lanes P02/P03 (Control1-5) were captured with T-ALL/ETP-ALL samples. Their soup was not
   removed and not bounded (D9).
8. **Age, and unrecorded factors.**
   - The healthy donors are 31.5 years younger in median. Age is bounded, not excluded (S9).
   - CMV, infection, transfusion, steroids, other medication and the time of aspiration are unrecorded in both arms.
9. **No public time stamp.** The pre-registrations are internal git tags.
10. **GATE_A1(4)** was registered as 3 vs 7 donors and run as 2 vs 8. S5's healthy ages were miscounted (true 5 / 5,
    median 32.5). The mixed channels hold 4 AML + 2 healthy distinct donors, not 7 + 3.
11. **S5's hit rule is lenient** (FDR 0.10, |log2FC| > 0.5): 252 of the 370 have |log2FC| ≤ 1. Claims are about sets
    (7.9). D10's declared reading is not usable, so "not T-specific" rests on post-hoc relabellings and G01 5.1.
12. **No comparator for acute illness.** The AML samples come from acutely ill patients at diagnosis. No non-AML patient
    comparator was run (G10(a)), so "AML-associated" includes being acutely ill.
13. **The source precondition is not met.** The verifier asked that the vendor-versus-biobank question be settled
    before the axis is interpreted. Data on disk cannot settle it (A1-BATCH section 0; A1-G04 :46-47).
14. **The September T_NK null** is reconciled by A1-G03. Its wording is bound by A1-G03 section 7.

## 2. Data and frozen definitions

**2.1 Samples, fits, codings.**
- **T**: S5's 25 samples (15 adult AML, 10 healthy).
  - Design: the T rows of `ws_s5_scores.csv`. Its `score` column (the in-sample score) is never read.
  - Counts: `PB$raw$T` of `/LARGE1/gr10634/gaozy/aml_niche_net/ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds`.
  - Fit (27's `fit_de`, as G03's K1 reproduced it, with `return_y = TRUE`):
    - S5's blacklist (307 genes) applied before filtering;
    - `filterByExpr(group = status)` plus CPM > 1 in at least half the samples;
    - TMM;
    - `estimateDisp(robust = TRUE)` and `glmQLFit(robust = TRUE)`.
  - The design columns are `(Intercept), statusAML, hb, gran, sexmale, subtype_share`. The coefficient is always passed
    by name. edgeR's `camera.DGEList` and `fry.DGEList` default to the last column.
- **B** (13 + 10), **NK** (12 + 9) and **Mono** (16 + 10; 21 of them in T's 25): each refitted with its own S5 design.
- **Sex**: metadata coding (as S5, G02 and G03) is primary. Control4 coded female is a sensitivity (4.2). Level names are
  `female` and `male`.
- **T-ALL pools**: Control1-5 are flagged in every table. Leaving them out is a sensitivity (4.3). The capture-day arm is
  not run (departure 7).
- **P370 and F402**: the text quotes P370. F402 values appear only as numbers in tables and in a 4.x parenthesis. No
  sentence names the 32 genes, reads them, or attributes a P370/F402 difference to them.

**2.2 Gene universes.**
- **U_T**: the 7,379 genes tested in S5's T fit.
- **U_A**: U_T minus the extended ambient ranking (`ws_s1_ambient_blacklist.csv`, `n_units ≥ 2`; 926 genes, 422 in
  U_T), 6,957 genes. U_A is primary for 3.1-3.3.
  - P370 lies inside U_A by construction.
  - U_A removes mostly immune and inflammatory HALLMARK members, and 5 IEG genes (EGR1, HSPA1A, JUNB, KLF2, SOCS3).
- B, NK and Mono use their own tested genes minus the extended ranking. The T × X interaction uses U_T ∩ X's tested
  genes, minus the extended ranking.
- **Implementation.** Fit on the full filtered DGEList, then subset rows with `keep.lib.sizes = TRUE` before camera. NB
  z-scores are computed per gene with fixed dispersions. Subsetting changes the background, and negligibly the variance
  prior (≤ 0.004; review1/L1).

**2.3 Gene-set families and gene lists (frozen).**
- **P-H, HALLMARK v2024.1 (pinned).** `/LARGE1/gr10634/gaozy/reference/h.all.v2024.1.Hs.symbols.gmt`, md5 `d3606c32…`,
  with its PROVENANCE file.
  - 45 sets have ≥ 15 genes in U_A.
  - Not testable: ANGIOGENESIS 5, PANCREAS_BETA_CELLS 9, HEDGEHOG 10, APICAL_SURFACE 13, NOTCH 13.
- **P-T, Gavish 2023 3CA T-cell metaprograms (decision 9.4).**
  - Source: the repo's MSigDB v2026.1 GMT (`LCC_proj/scripts/panels/msigdb.v2026.1.Hs.symbols.gmt`, md5
    `30893650e0c0735a39819f9f7691f0d6`; committed in aac665f; no provenance file).
  - If approved, the 22 sets are extracted to `/LARGE1/gr10634/gaozy/reference/gavish2023_3ca_T.v2026.1.gmt`, with a
    PROVENANCE note and an asserted md5.
  - **Interpreted: 19 sets.** The three UNASSIGNED programs (CD4_UNASSIGNED, CD8_UNASSIGNED_1, CD8_UNASSIGNED_2) are
    reported, not interpreted.
  - They are tumour-infiltrating T-cell programs, tested on whole-T pseudobulks.
  - For each set, the report lists its genes absent from U_A (untested or ambient-excluded) and its H share.
  - **Near-duplicate pairs** (Jaccard ≥ 0.3 in U_A) count as one unit in 7.3: CELL_CYCLE, CYTOTOXIC, GLYCOLYSIS_MYC,
    INTERFERON, CD4 NAIVE_1 / CD8 MEMORY_NAIVE_1, NAIVE_2, CD4 STRESS_HSP / CD8 HEAT_SHOCK.
- **Handling lists, split by delay direction.**
  - Savage 2021 T delay-18-h consistent set: 792 genes, 295 up / 497 down; md5 `5ee444a9…`.
  - Project IEG list: 26 genes, `config_ws.R` `WS_B0_IEG`. It has no delay direction.
  - **H_up** = IEG ∪ Savage-up; **H_down** = IEG ∪ Savage-down; **H** = both.
  - Each set's H fraction is reported as a description, next to the H fraction of size-matched random sets. No reading
    depends on it.
- **Per-compartment delay lists (decision 9.14).** These are derived without a download from the on-disk Table S4
  (`mmc5.xlsx`), using the T set's rule (PROVENANCE):
  - adjusted p < 0.05 at 18 h vs 2 h in at least one cell type of the compartment, in the healthy donor of both
    experiments;
  - the same direction in both;
  - no significant opposite call.

  | compartment | Savage cell types |
  |---|---|
  | B | B_Naive, B_Memory |
  | NK | NK_dim, NK_bright |
  | Mono | CD14+_Monocytes, CD16+_Monocytes |

  The files go next to the T set, with the PROVENANCE note extended. K9 checks the rule by reproducing the T set.
- **Opposite-sign class.** Genes in the Gavish CYTOTOXIC union that are also Savage-down, within U_A. There are 14:
  ANXA1, APOBEC3G, CCL5, CD27, CKLF, GZMA, GZMH, GZMK, HOPX, KLRD1, KLRG1, NKG7, TIGIT, TRGC2. An effector shift predicts
  them up; delay in the D8 orientation predicts them down.
  - Seen before (1.2): GZMH +0.79, KLRG1 −0.54 and TIGIT −0.01, all n.s.; GZMK is not among the 402.
- **lncRNA**: `gene_type "lncRNA"` and no other type in `refdata-gex-GRCh38-2020-A/genes/genes.gtf`, the annotation of
  the data's own gene axis. 305 genes in U_T, 274 in U_A.
- **Clock and glucocorticoid list (uncited; declared from canonical knowledge).**
  - Clock: ARNTL, CLOCK, NPAS2, PER1, PER2, PER3, CRY1, CRY2, NR1D1, NR1D2, DBP, TEF, HLF, CIART, BHLHE40, BHLHE41.
  - Glucocorticoid-induced: FKBP5, TSC22D3, ZBTB16, KLF9, DDIT4, PER1, SGK1, DUSP1, TXNIP, CXCR4.
  - In U_A: ARNTL, CLOCK, PER1, PER2, CRY1, CRY2, NR1D2, DBP, FKBP5, ZBTB16, KLF9, DUSP1, TXNIP, CXCR4.
  - Overlap with H: PER1, NR1D2 and CXCR4 are Savage-up; FKBP5 is Savage-down; DUSP1 is an IEG.
  - Already seen (1.2): the axis status of PER2, ARNTL, DBP, FKBP5, CXCR4, CRY1, PER1, TSC22D3, ZBTB16 and KLF9.

**2.4 Cells (S5b).**
- Kept T cells of the 25 samples: `PB$cells` with `kept` and `group4 == "T"`, 14,609 cells.
- Counts: the RNA `counts` layer of `02_seurat_objects/04_annotated/GSE185381/<sample>.rds` (36,601-gene axis).
- **UMI** = the column sum over the 36,601-gene axis (= `nCount_RNA`).
- Labels: the nine BoneMarrowMap T fine labels.

## 3. Analyses

### 3.1 Gene-set description (the interpreted families)

- **Test.**
  ```
  edgeR::camera.DGEList(y_A, index, design, contrast = "statusAML",
                        inter.gene.cor = NA, allow.neg.cor = FALSE, sort = FALSE)
  ```
  `y_A` is S5's T DGEList with its dispersions, restricted to U_A.
- **Families.**
  - P-H: 45 sets.
  - P-T: 19 sets (if 9.4).
  - BH within each family. **Flagged** = FDR < 0.05.
  - Direction = camera's `Direction`, relative to the average U_A gene.
- **Reported for every set:**
  - set size, genes in U_A, coverage;
  - NGenes, the estimated correlation, Direction, PValue, FDR;
  - **s**, the direction-signed set shift. It is the sign of Direction (+1 for Up) times (mean logFC of the set's genes
    minus mean logFC of the other U_A genes). If s ≤ 0, the set is marked "shift opposes the camera direction" and gets
    no further reading;
  - the H fraction with its random reference;
  - P370 members (up / down) with expected counts (members in U_A × arm size / 6,957);
  - F402 members, as numbers, with expected counts over U_T (/ 7,379).
- **Names of flagged Gavish sets.** A flagged Gavish set is named by its program name in quotes, with "as tested: n of
  50 genes, m of them in H".
  - "Exhausted", "dysfunctional" and "regulatory" are not used: PDCD1, HAVCR2, CTLA4 and CXCL13 are not in U_A, and the
    T_REG set lacks FOXP3, IL2RA and CTLA4 in U_A.
- **Leading genes of flagged sets.** fgsea's leading edge (4.4), as a description.

### 3.2 Handling (for every flagged set)

1. **H fraction.** The fraction of the set's U_A genes in H, against 1,000 random sets of the same size matched on
   status-free pooled logCPM (10 bins; seed SEED + 706). Description only.
2. **Direction-matched removal (primary).**
   - An Up-flagged set is re-tested in U_A minus H_up (6,683 genes); a Down-flagged set in U_A minus H_down (6,517).
   - Each run tests every set of the family that keeps ≥ 10 genes, with the settings of 3.1.
   - The assumption, taken from D8's orientation: AML samples resemble longer delay, so a delay artefact would sit in
     genes whose delay sign equals the flag.
3. **Reading** (each set on its own p; review1/L1 F2). Let p* be the largest primary p value among the family's flagged
   sets, which is the BH cut-off actually used.
   - **"consists mostly of the listed handling genes"**: fewer than 10 genes remain. K3 says none will.
   - **"carried by the listed handling genes"**: s after removal is below half of s before.
   - **"not carried by the listed handling genes"**: s after removal is at least half of s before, the direction is the
     same, and the p value is ≤ p*.
   - **"not resolved after removal"**: s holds but the p value exceeds p*.
4. **Full removal (sensitivity).** The same reading in U_A minus H (6,259 genes). CD8 MEMORY_NAIVE_1 (8 genes) and CD4
   NAIVE_1 (9) would read "consists mostly of the listed handling genes".
5. **Direction check** (descriptive; sets with ≥ 8 Savage members in U_A). The share of the set's Savage members whose
   AML-versus-healthy sign equals their delay direction:
   - ≥ 0.75: "runs with the delay direction";
   - ≤ 0.25: "runs against it";
   - otherwise: "mixed".

   A set that runs against the delay direction is never written "carried". The fold-gene-level agreement (0.994) was
   already seen (1.2, L13).
6. **The opposite-sign class** (2.3; descriptive). One camera test (estimated correlation; U_A background).
   - "Up": "runs against the delay direction, as an effector shift would".
   - "Down": "runs with the delay direction".
   - p reported.

### 3.3 Compartments

**(a) Each compartment alone (reported).**
- camera as in 3.1 on each compartment's refitted S5 fit (B, NK, Mono), in that compartment's tested genes minus the
  extended ranking.
- Sets: the P-H and P-T sets with ≥ 15 genes there, plus the compartment's own delay-up and delay-down lists and the T
  lists. BH within each compartment and family.
- Reading of the delay lists (description): tube-level delay predicts that every compartment aligns with its own delay
  list in the D8 orientation (delay-up "Up", delay-down "Down"). T aligning while B and NK do not would argue against
  tube-level delay; all aligning is consistent with it.
- Monocytes are reported, not read (often clonal in AML; A1-G01 5.1).

**(b) T versus X, for X = B and X = NK (read for flagged sets).**
- **Donors**: those with both T and X (B 23; NK 21).
- **Genes**: tested in both S5 fits, minus the extended ranking.
- **Model**: one edgeR model on the stacked raw counts, with TMM and `estimateDisp(robust = TRUE)` on the stacked design.
  - Columns:
    - one fixed effect per donor;
    - `cT`, the T-compartment indicator;
    - `cT:AML`;
    - `cT:hb`, `cT:gran`, `cT:sexmale`;
    - `cT:CD8share` and `cT:Xshare`.
  - `cT:AML` is T's adjusted AML effect minus X's, on the same donors.
- **Test**: camera, estimated correlation, contrast `cT:AML`. BH within each X and family (X-B: 41 HALLMARK + 18 Gavish;
  X-NK: 43 + 21).
- **ρ_X** = s in X / s in T, both computed on the paired donors and the common genes. It carries a 95 % donor-bootstrap
  interval (1,000 draws; seed SEED + 707).
- **Reading for each flagged set:**
  - "larger in T than in X": the interaction has FDR < 0.05 in the set's direction;
  - "larger in X than in T": FDR < 0.05 in the opposite direction;
  - "shared with X": otherwise, if ρ_X ≥ 0.5;
  - "not resolved": otherwise;
  - "not testable in X": fewer than 15 genes in common.
- NK keeps less information than B (51 % vs 72 %, review1/L1), so more NK calls will read "not resolved".
- K5 calibrates (b). If K5 fails, every reading of (b) is "not readable".

### 3.4 Gene identity in healthy marrow T subtypes (descriptive)

- **References.**
  - **G**: the 10 GSE185381 healthy donors.
  - **O**: Oetjen (GSE120221), 20 donors, libraries pooled per donor as in S9, genes mapped by symbol.
  - In each reference, a label is tagged if at least 5 donors have ≥ 20 kept cells in it. In G these labels are CD4
    Naive (8), CD4 Central Memory (7), CD4 Effector Memory (6), CD8 Naive (5), CD8 Central Memory (6), CD8 Effector
    Memory 1 (7) and CD8 Effector Memory 2 (9).
- **Statistic.** For gene g, label k and eligible donor d:

  L = log2((CPM of g in label k of donor d + 1) / (CPM of g in donor d's other kept T cells + 1)).

  M = the median of L over donors; c = the share of donors with L > 0. "k-high" if M ≥ 1, c ≥ 0.8 and the median CPM in
  k ≥ 10.
- **Classes**: naive (CD4 or CD8 Naive high); central memory; effector memory (CD4 EM, CD8 EM1 or EM2 high).
- **Consensus.** A class is stated for a gene only when both G and O tag it. A gene absent from O's axis is
  "not assessable". Each reference is also reported separately.
- **Readout.**
  - For the P370 up and down arms and G02's 25-gene core: the count per consensus class. The base rate comes from U_A,
    and the 1,000 matched random sets (10 logCPM bins; seed SEED + 704) are drawn from U_A.
  - F402: the same, with U_T.
- **Composition-only expectation.** For each gene:

  comp_g = log2((Σ_k s_AML,k · P_k,g + 1) / (Σ_k s_H,k · P_k,g + 1))

  - P_k,g is G's pooled healthy CPM profile of label k, over all nine labels.
  - s_AML,k and s_H,k are the mean kept-T label shares of the AML and healthy samples (as in A1-G01).

  This is the log-fold change expected from the label mix alone, with no AML expression used. Reported per class and
  arm: the median comp_g, and the share of genes whose comp_g has the arm's sign with |comp_g| ≥ 0.25.
- **Sensitivities.** Threshold 0.585; G without Control1-5 (5 donors; labels with ≥ 4 eligible donors; c ≥ 0.75).

### 3.5 S5b, the single-cell layer (T only; dated amendment of S5 4.2; descriptive)

**(a) Scores.**
- `UCell::ScoreSignatures_UCell` (2.10.1), maxRank 1,500, ties "average".
- Run on the full 36,601-gene count matrix of the kept T cells.
- The P370 up set and down set are scored separately. F402 is a sensitivity that reuses the P370 K6 states.

**(b) Pipelines.**
- **S5b-raw (as registered)**: raw counts.
- **S5b-1000 (reading)**:
  - cells below 1,000 UMI are dropped (healthy 7.1 %, AML 0.6 %);
  - every other cell is downsampled to exactly 1,000 UMI by sampling its UMIs without replacement, in base R.
  - Seeds: `set.seed(SEED + 7100 + i)`, where i is the sample's position in the T rows of `ws_s5_scores.csv`. Cells are
    in `PB$cells` order, one sample per call.
- **S5b-1500 (sensitivity)**: the same, at 1,500 UMI, with cells at or below 1,500 UMI kept unthinned (seed
  SEED + 7300 + i).
- The share of cells dropped or thinned is reported for every sample.

**(c) Reference.**
- For each fine label: the mean, SD, P5 and P95 of the healthy cells' scores under the same pipeline, and
  z = (score − mean_k) / SD_k.
- AML patients are compared with all 10 healthy donors. Each healthy donor e is compared with the other nine.
- For the K6 pseudo-patient d, the reference is the other nine, and h_max(−d) is the largest mean z of the other nine
  donors, each scored against the eight left after removing d and itself (nested leave-one-out).
- Under 4.3 (five donors), a label needs ≥ 20 reference cells from ≥ 2 donors. Below that, its cells use the pooled CD4
  or CD8 reference. The code asserts this.

**(d) Registered summaries (all pipelines).** Per patient:
- the mean UCell;
- the mean z;
- F95, the fraction of cells above their label's healthy P95 (strictly greater; quantile type 7);
- for the down set, also F5.

S5b-raw is reported as registered, with the depth warning of 1.4.

**(e) Reading rule (S5b-1000; per AML patient).**
1. **Eligibility.** ≥ 50 kept T cells after dropping. 13 of the 15 qualify; AML2451 (47) and AML0310 (20) are "too few
   cells".
2. **Detection floor** (K6e). The median number of detected up-set genes per healthy reference cell must be ≥ 7. That is
   the 10th percentile of random tested sets at 1,000 UMI (1.4). Below it, shape calls are "not readable" in advance;
   only the gate and the registered summaries are reported.
3. **Gate, up set.** The patient passes if

   mean z − 1.96 × sqrt(var(z) / (v_H × n)) > h_max.

   v_H is the variance of the pooled healthy reference z, and h_max the largest of the 10 healthy leave-one-out mean z.
   A patient who does not pass gets "no detectable shift".

   **Gate, down set.** The patient passes if mean z + 1.96 × sqrt(var(z) / (v_H × n)) < h_min, with h_min the smallest
   of the healthy leave-one-out means.
4. **Shape, up set only.**
   - f0 = the kernel density of the pooled healthy reference z. R `density()` defaults (Gaussian kernel, bw.nrd0),
     2,048 points, from min − 3 to max + 8, floored at 1e-8, `approx(rule = 2)`.
   - Uniform model: z ~ f0(z − δ), δ by maximum likelihood (a grid from −3 to 20 in steps of 1, then `optimize`
     within ± 1).
   - Mixture model: z ~ (1 − π) f0(z) + π f0(z − Δ), with π ∈ {0.1, …, 1.0} and Δ ∈ {0, 0.5, …, 20}.
   - LR = max(0, ℓ_mixture − ℓ_uniform). π̂ is the maximising π; ties go to the smaller π; π̂ = 1 if LR < 10^−6.
   - Calls:
     - "sub-population": π̂ ≤ 0.5 and LR ≥ 10;
     - "uniform": π̂ ≥ 0.9 or LR < 2;
     - "not separable": otherwise.
   - These functions are copied verbatim from D2's `13_sim.R` (`mk_ref`, `zf`, `f0`, `pstats`) and `17_final_rule.R`
     (md5s in 8.4).
   - A call is readable only as K6b allows (graded bars, section 6).
5. **Cohort reading** (first match wins; thresholds unvalidated):
   1. "not readable": fewer than 6 patients pass the gate, the detection floor fails, or K6b makes both calls
      unreadable.
   2. "mostly uniform": at least 2/3 of the patients with a detectable shift are "uniform" (written "uniform or not
      separable" if K6b leaves uniform calls unreadable).
   3. "mostly sub-population": at least 2/3 are "sub-population".
   4. "mixed": "uniform" and "sub-population" each occur in at least 2 patients.
   5. "mostly not separable": otherwise.

**(f) Per label (description).** CD4 Naive, CD4 Central Memory and CD8 Effector Memory 1 and 2 have ≥ 50 cells in at
least 7 AML patients. For each (patient, label) pair with ≥ 50 cells: the mean z, LR and π̂, with the pooled f0. No gate
and no call word.

**(g) Sensitivities.** S5b-1500; the healthy reference without Control1-5; F402.

**(h) Unit rule.** No per-cell p value is computed. Cells enter only through per-patient summaries (S5's principle P3).

### 3.6 Gene classes and the overlap table (descriptive)

**(a) Gene classes removed from the cross-fitted score.**
- Method: G02's stored fold lists, scored as in A1-BATCH's post-hoc script (log-CPM over the 25, out-of-fold z, up minus
  down). gap_oof = mean AML score − mean healthy score.
- Rows:

  | genes removed | gap_oof |
  |---|---|
  | none | known: 1.3782 |
  | Savage | known: 1.3184 |
  | Savage + IEG | known: 1.3228 |
  | lncRNA | |
  | clock + glucocorticoid | |
  | lncRNA + clock/glucocorticoid + H | |

- Also reported: the lncRNA share in the P370 arms against U_A (274 / 6,957); F402 arms against U_T, as numbers.
- Reading: "the score does not depend on [class] (gap g)", or "depends on it". The clock/glucocorticoid row adds: "this
  cannot exclude time of day or steroids, which act beyond these genes and are unrecorded".

**(b) Overlap table.**
- Lists on disk: Savage up / down; IEG; the 50 HALLMARK sets; the 22 Gavish sets; the cNMF normal MP10 (naive) and MP4
  (stress) programs, labelled "derived from this project's cells".
- Downloaded lists, only if 9.5 is approved and the files are on disk before the run (Lasry after reading the Author
  Correction): Lasry 2023 (246; iScore 38), Rutella 2022 (IED172, IED68, PI24), Mazziotta 2024 (top-100; the 48 genes
  if printed).
- Columns:
  - list and size;
  - in U_T and in U_A;
  - overlap with the P370 arms and G02's core, with expected counts over U_A;
  - F402 overlaps as numbers, over U_T;
  - the overlapping gene names (P370 genes only).
- No p values.
- Public repo: names, counts and citations only. Lasry's lists are malignant-cell signatures: positioning only.

## 4. Sensitivities (reported; none changes a reading)

4.1 All 7,379 genes (U_T) for 3.1-3.3.
4.2 Control4 coded female for 3.1-3.3. The recode asserts levels `female` and `male`.
4.3 Without Control1-5 (15 vs 5) for 3.1, 3.2, 3.4 and S5b's reference.
4.4 Other views of the same contrast:
    - camera 0.01;
    - `fgseaMultilevel` preranked on sign(logFC) × sqrt(F) from S5's QL refit (minSize 15, maxSize 500, eps 0, seed
      SEED + 705);
    - fry directional;
    - a joint BH over P-H ∪ P-T.

Reading: for a flagged set, the main text names in one parenthesis any of 4.1-4.3 under which the set loses FDR < 0.05
or flips direction. A set significant only in a 4.4 view is "not detected by the primary".

## 5. Multiplicity and what is read

| family | analysis | tests | FDR | status |
|---|---|---|---|---|
| P-H | 3.1, HALLMARK | 45 | BH within | interpreted |
| P-T | 3.1, Gavish (if 9.4) | 19 (+3 UNASSIGNED reported) | BH within | interpreted |
| removal | 3.2, direction-matched, two universes | ≤ 45 + ≤ 22 each | own p against p* (3.2) | read for flagged sets |
| full removal | 3.2(4) | 45 + 15 with ≥ 15 genes (all ≥ 10 except 2) | as 3.2 | sensitivity |
| X-B, X-NK | 3.3(b) | 41 + 18 / 43 + 21 | BH within | read for flagged sets |
| per compartment | 3.3(a) | each | BH within | reported |
| sensitivities and views | 4 | each | BH within | reported |

**Rule.** Nothing outside P-H and P-T is interpreted on its own. No secondary family can create a flagged set; it can
only qualify one. BH is never pooled across families, except in the joint-BH view of 4.4. Each family is read with its
own BH; the joint BH bounds the cost of having two.

## 6. Known answers (run first; no readout is printed until every check of a mode has a state)

There are three states, as in G04:
- **PASS**;
- **FAIL deterministic**: a code bug; fixed with a dated note, then rerun;
- **FAIL statistical**: the reading named in the row becomes "not readable" or "not interpretable"; the other readouts
  still run, labelled descriptive.

Rows marked "reported" carry no state.

| id | check | bar | type |
|---|---|---|---|
| K0 | md5s of the frozen inputs (2.x; Oetjen `ws_pb_GSE120221.rds` and its donor table; the cNMF tables; `mmc5.xlsx` and the derived delay lists; any downloaded file with its PROVENANCE sha256; the W4 scores; the GTF); the IEG list equals `WS_B0_IEG` | equal to recorded values | deterministic |
| K1 | S5 T refit: genes tested; hits up / down; the hit set equals the frozen list; max \|logFC − `ws_s5_de_T.csv`\|; `colnames(design)[2] == "statusAML"` | 7,379; 402 (171 / 231); TRUE; < 1e-9; TRUE | deterministic |
| K1b | B, NK, Mono refits: genes tested; max \|logFC diff\| | 4,292 / 4,151 / 9,252; < 1e-9 | deterministic |
| K2a | 500 random sets from U_A, each copying the logCPM-bin composition of a randomly drawn HALLMARK set (seed SEED + 701); share with p < 0.05; Up/Down split and median statistic reported | ≤ 0.10 | statistical (3.1 "not interpretable") |
| K2b | code check. 10 Up plants (y′ = y + Poisson(y(2^0.6 − 1)) in AML samples) and 10 Down plants (binomial thinning at 2^−0.6), each a random 50-gene U_A set; one set per refit on the original rows (no re-filtering); the plant joins the 45 HALLMARK sets for BH (seed SEED + 702) | flagged in its direction ≥ 8/10 per direction; never the wrong direction | statistical (3.1 "not interpretable") |
| K2c | real-fit power, as K2b: 20 whole-set plants at 0.3 log2, and 20 plants with 30 % of genes at 0.5 log2 | share flagged; quoted in the no-flag sentence (7.1) | reported |
| K3 | universes 6,683 / 6,517 / 6,259; every P-H and P-T set keeps ≥ 10 genes in both direction-matched universes; the 7 Gavish sets below 15 after full removal are those of 1.4; the Gavish pairs of 2.3 | as stated | deterministic |
| K4 | tags, positive control, in G and in O: of CCR7, LEF1, TCF7, SELL, and of GZMB, GNLY, NKG7, PRF1, GZMH (all in U_T) | ≥ 3 of 4 tagged naive and ≥ 3 of 5 tagged effector memory; none tagged in the opposite class | statistical (3.4 "not readable" for that reference) |
| K5 | code check: real stacked counts (seed SEED + 703); 20 random 50-gene sets per type, one per refit on the original rows; (i) T-only 0.6 log2 in AML T; (ii) shared 0.6 in AML T and AML X; (iii) half-shared (0.6 in T, 0.3 in X) and (iv) T-only at 0.3; the 3.3(b) reading applied to every planted set | (i) "larger in T" ≥ 14/20; (ii) "larger in T" ≤ 2/20 and "shared" ≥ 14/20; (iii) and (iv) reported | statistical (3.3(b) "not readable") |
| K6a | thinning and UCell | same seed twice gives identical scores; every kept cell sums to exactly 1,000; cells below 1,000 are dropped; a donor thinned alone equals the same donor in the full run | deterministic |
| K6b | planted shape on real healthy cells. Each healthy donor i in turn is a pseudo-patient (seed SEED + 7200 + i; (a) drawn before (b)). Plants on the real P370 up-set genes with synthetic per-gene β_g = 0.6 × min(2.5, 0.5 + Exp(mean 0.4)); S5's logFCs are never used. (a) uniform: y′ = y + Poisson(μ̂_g(2^β_g − 1)) in every cell. (b) exactly round(n/5) cells with multiplier m_g = (2^β_g − 1) Σ_c μ̂_cg / Σ_{c∈S} μ̂_cg. μ̂ = cell UMI × the pooled profile of the cell's label in the other nine donors. Then the exact 3.5 pipeline with the nested reference. Reported per donor (n; shares dropped and thinned). | graded: sub-population calls readable if (b) gives "sub-population" in ≥ 8/10 and (a) in ≤ 1/10; uniform calls readable if (a) gives "uniform" in ≥ 7/10, otherwise AML "uniform" calls read "uniform or not separable"; both fail → S5b cohort "not readable" | statistical |
| K6c | matching: per plant, log2[(Σ_up y′ / Σ_cells UMI′) / (Σ_up y / Σ_cells UMI)] on the donor's pre-thinning counts; the median over donors of \|(a) − (b)\| | ≤ 0.05 | statistical (as K6b) |
| K6d | reported. Each sample's dropped and thinned shares and median UMI; the unplanted held-out donors' gate results ("cannot fail by construction": the gate's specificity is ≤ 1/11 per exchangeable donor and does not cover AML-specific depth or handling); D2's synthetic null cited | none | reported |
| K6e | detection floor (3.5(e)2): median detected up-set genes per healthy reference cell under S5b-1000 | ≥ 7 | statistical (shape calls "not readable") |
| K7 | gene-class removal: rebuilt out-of-fold scores equal `ws_a1_g02_oof_scores.csv`; gap with no removal / without Savage / without Savage + IEG | max \|diff\| < 1e-9; 1.3782; 1.3184; 1.3228 | deterministic |
| K8 | per-cell inputs: kept-T sums equal `PB$raw$T`; label sums equal `PB$fine`; 14,609 cells; every set gene is in the matrix rownames; labels are exactly the nine T labels, with no NA; lncRNA count in U_T; median detected set genes per cell by status (reported) | exact; 305 | deterministic |
| K9 | delay-list rule: applied to Table S4, it reproduces `savage2021_T_delay18h_consistent.csv` | 792 genes (295 / 497), identical | deterministic |

Run-time blindness:
- K2, K5 and K6 read real counts after the tag. K6 uses healthy cells only.
- No AML-versus-healthy set or per-cell value is printed until every check of the mode has a state.

## 7. Reading rules and fixed sentences

**7.1 Sets.** Template:

> "In the S5 T-cell contrast (15 adult AML vs 10 healthy donors, GSE185381), set X had a higher [lower] average
> AML-versus-healthy statistic than the average tested gene (camera with estimated inter-gene correlation, FDR q; same
> contrast that defined the axis; descriptive)."

- A Gavish set adds: "(a program defined in tumour-infiltrating T cells; as tested here n of 50 genes, m of them in the
  handling lists; tested on whole-T pseudobulks, so the CD4/CD8 prefix does not localise the shift)".
- Then the handling reading (7.3) and the compartment reading (7.6).

**If no set is flagged:**

> "No HALLMARK [or T-cell-program, if 9.4] set differed from the average tested gene at FDR 0.05 in a competitive test
> that accounts for inter-gene correlation (n = 25). On this fit, planted whole-set shifts of 0.3 log2 were flagged in
> [K2c] of 20, so this is not evidence that the axis lacks set structure."

**7.2 Secondary views.** A set flagged only by a 4.4 view is "not detected by the primary". fry Mixed results never
appear.

**7.3 Handling.**
- **Every flagged set** carries its 3.2 reading, its H fraction with the random reference, and its direction check.
- "Not carried by the listed handling genes" is always followed by:

  > "(the lists represent 18 h ex-vivo delay of fresh blood from two donors and immediate-early stress; freezing,
  > thawing, shipping and sorting are not represented)".

  It is never shortened to "not a handling effect".
- **Consequence rule** (fixed now; decision 9.13).
  - Count each direction separately, over P-H and P-T flagged sets. A Gavish pair is one unit: it is bound if all its
    flagged members are bound.
  - "Bound" means "carried by" or "consists mostly of".
  - If at least one unit is flagged and at least half of the flagged units are bound, the paper states as a main
    limitation:

    > "The set-level description of the [up/down] arm rests on genes that also change with 18 h ex-vivo blood
    > processing delay (Savage 2021) or are immediate-early genes; without them it is [no longer carried / not
    > testable]. These lists cover at least a quarter of the tested genes of 20 of the 22 T-cell programs used here (9 of
    > them naive, cytotoxic, regulatory or dysfunction programs), and most naive-identity genes also fall with delay, so
    > a handling reading and a T-cell-state reading of this arm cannot be separated with these data."
  - Reported beside it: the share of the union of the flagged sets' U_A genes that is in H_up or H_down, against the
    base rate.

**7.4 Gene identity and composition (3.4).**
- Template:

  > "k of the [up/down] genes are [class]-high in healthy marrow T cells in both references (GSE185381 healthy donors
  > and Oetjen; paired within donor; base rate b; matched-random percentile p). Such genes are expected both from a
  > change within cells and from a lower [naive] share in AML (composition-only expectation: median log2 c); the tags
  > cannot tell these apart."
- Then G01's sentence, verbatim (A1-G01 11.5):

  > "In GSE185381 (15 adult AML, 10 healthy donors), the mix of BoneMarrowMap T subtypes accounts for about 5% of the
  > cross-fitted T-imprint gap (C 0.045, 95% CI −0.100 to 0.191; C_H 0.051, −0.107 to 0.193; C_A 0.039, −0.104 to
  > 0.196), and a shift is detected inside the four main T labels (S_w 0.664, 95% CI 0.497 to 0.772; a sign check, not a
  > share of the gap). The intervals are descriptive, with the fold lists held fixed. With gamma-delta and MAIT cells as
  > separate labels the composition share is 0.057 (the flags are lower bounds)."

  Followed by: "C is a share of the score gap, not of gene classes."
- **A naive-like description adds:**

  > "The healthy donors are 31.5 years younger in median; in healthy marrow the weak age association of the score runs
  > through the naive-T share (S9)."
- **Effector-like trigger.** Any effector-memory consensus tag in a reported arm, any flagged CYTOTOXIC set, any CD8 EM
  per-label sub-population pattern (3.5(f)), or any "shared with NK" adds:

  > "CMV serostatus is unrecorded in both arms; CMV expands TEMRA-like cells inside the CD8 effector-memory labels and
  > adaptive NK cells, and rises with age."

**7.5 S5b.**
- Template:

  > "With depth-equalised per-cell scores (cells thinned to 1,000 UMI), among the 13 AML patients with ≥ 50 eligible T
  > cells, u were uniform, s sub-population, n not separable and z showed no detectable shift (patient-level
  > summaries)."

  Then the cohort reading.
- A uniform reading adds: "a uniform shift is what a change in all cells would give, and also what handling, ambient RNA
  or a residual depth difference would give."
- A sub-population reading adds: "this can be composition below the label level, a clonal expansion, label error, cells
  rich in ambient RNA, or doublets."
- The wording "expansion of a T-cell state" (gaps_readable :456) is not used, because a minority call has several
  readings. No S5b result changes G01's category.

**7.6 Compartments (3.3(b)).**
- Allowed readings: "larger in T than in B [NK]", "larger in B [NK] than in T", "shared with B [NK] (ρ, 95 % CI)",
  "not resolved", "not testable in B [NK]".
- Always added: "(two compartments of the same donors; monocytes and other cells are not compared this way, so this is
  not a statement of T-specificity)".
- "Larger in T" does not mean absent in X. "Shared" means "not restricted to T", never "systemic biology": ambient RNA
  and handling are shared by every compartment of a sample.

**7.7 Context that travels with any G07 result** (verbatim).
- A1-BATCH (GATES.md:293):

  > "Allowed wording: none of the six within-study processing factors, tested one at a time, explained or weakened the
  > frozen direction; sample source, storage, sorting method and chemistry version are constant within arms or
  > unrecorded, and were not tested."
- A1-BATCH 12.5, must add:

  > "The imprint direction is shared with B and NK cells of the same donors (D10: R 0.69 and 0.93), so it is not
  > T-cell-specific. It overlaps a blood-processing-delay program in the same direction (D8), and adding the delay score
  > to MAIN lowers R to 0.82 (fhat 0.30; descriptive, no reference). Dropping the Control1-5 group lowers R to 0.81
  > (fhat 0.31; percentile 0.10), a share G1 could not have flagged."
- A1-G01 (GATES.md:319): "beyond T-subset composition (detected inside the main labels)". Monocytes: "monocytes in AML
  samples are often clonal, so a monocyte shift may be leukaemic biology rather than a shared environment or handling
  effect."
- A1-G04 (GATES.md:346):

  > "The direction was not shown in Chen2023 (exact p = 17/126); a common shift larger than 1.36 of GSE185381's
  > selection-corrected gap is excluded (one-sided 95.2%, common-shift model)." "5 vs 4 cannot exclude an effect of the
  > GSE185381 size." The direction stays "shown in one study".
- S9 (GATES.md:241, bound form):

  > "AML, with an age share bounded below one third in healthy marrow at the point estimate (one protocol; upper 95 %
  > bound 0.46)."
- Added by G07:

  > "Samples were taken at diagnosis; acute illness, infection, transfusion and medication are unrecorded, and no
  > non-AML patient comparator was run, so 'AML-associated' includes being acutely ill. The commercial-versus-biobank
  > source of the arms is not separated from disease."

**7.8 Never written, under any outcome.**
- **Blanket rule.** Every may-not item of A1-BATCH 12.5, A1-G01 11.5, A1-G03 section 7, A1-G04 11 and S9 applies
  unchanged.
- The 32 ambient-prone genes "stay excluded from every statement" (GATES.md:241).
- Mechanism and causation:
  - "pathway X is activated, induced or suppressed in AML T cells"; "drives"; "reflects signalling"; "caused by AML";
  - an RNA signal "equals" ligand binding or metabolite flux (handoff :466);
  - "AML-left imprint" (handoff :175).
- Evidence claims:
  - "independent evidence", "validates", "confirms the imprint";
  - "matches no published signature" as a novelty claim; "first to compare";
  - "the axis is not a pathway";
  - T3 as replication evidence; the in-sample 3.84 SD or AUC 1.0 as an effect size.
- T-specificity and handling:
  - "T-cell-specific program or state"; "T-cell-intrinsic"; "cell-intrinsic";
  - "not a handling effect"; "not an artefact"; "not a generic sample-level artefact"; "robust to processing";
  - "consistent within capture days" without the same-channel result.
- Composition: "not caused by subtype proportions" (from the in-sample subtype refits).
- Labels:
  - "exhausted", "dysfunctional", "regulatory" or "senescent" for a flagged set;
  - "a new cell state";
  - "opposite to published ISG-up reports".
- Other:
  - any per-cell p value;
  - any biological reading of GNA15;
  - any statement from PROGENy, CytoSig or ligand-receptor evidence;
  - the September null worded as underpowered, as the wrong model, or as tracing to pooled samples (A1-G03 section 7).

**7.9 Claims are about sets.** A single gene is named only as a member of a reported set, a leading edge, or G02's
25-gene core. No single gene gets a biological reading.

## 8. Code, outputs, compute (plan; the scripts are written after the tag)

**8.1 Scripts.**
- `scripts/ws/43_a1_g07.R`, modes:
  - `--known` (K0-K5, K7, K9);
  - `--sets` (3.1, 3.2, 4 for 3.1-3.2);
  - `--cross` (3.3 and its 4.x);
  - `--tags` (3.4);
  - `--classes` (3.6(a));
  - `--overlap` (3.6(b));
  - `--combine` (summary and GATES block).
- `scripts/ws/44_a1_g07_cells.R`, modes:
  - `--known` (K6, K8);
  - `--score` (3.5).
- `scripts/ws/submit_ws_a1_g07.sbatch`: dated read-only copies of the script, `OMP_PROC_BIND=false`, a smoke switch.
  `script_md5` in every table; a write guard on real outputs.

**8.2 Outputs.**
- `results/tables/ws/ws_a1_g07_{checks, sets, sets_removal, sets_sens, opposite, compartments, cross, tags, classes,
  overlap, s5b_patients, s5b_labels, s5b_checks, summary}.csv`.
- GATES.md block `A1-G07`, with a "Reading" line, not a "Verdict".
- `/LARGE1/gr10634/gaozy/aml_niche_net/ws_a1_g07/`: per-cell scores, code copies, smoke outputs.

**8.3 Seeds** (SEED = 491638; offsets ≥ +700 are free, B 8.2).

| offset | stream |
|---|---|
| +701 | K2a |
| +702 | K2b, K2c |
| +703 | K5 |
| +704 | 3.4 random sets |
| +705 | fgsea |
| +706 | 3.2(1) random sets |
| +707 | 3.3(b) bootstrap |
| +709 | smoke counts |
| +710 | smoke cells |
| +711-720 | reserved |
| +7100 + i | S5b-1000 thinning of sample i |
| +7300 + i | S5b-1500 thinning of sample i |
| +7200 + i | K6 plants of healthy donor i |

**8.4 Simulation provenance.** The rules rest on D1, D2 and review 1. The scripts and their md5s:
- D1, `inventory/scratch_D1/`: `02_sim_geneset.R` `94730713…`, `03_sim_crosscomp.R` `30079a4d…`, `04_summarise.R`
  `a2159db9…`, `05_crosscomp_readout.R` `76028d73…`.
- D2, `inventory/scratch_D2/`: `10_engine.R` `11573063…`, `13_sim.R` `bf631cdb…`, `14_summarise.R` `28c46716…`,
  `16_rule_table.R` `61e976a8…`, `17_final_rule.R` `e52077fd…`.
- Upstream inputs: D2's `01_real_structure.R`, `real_background_profile_nonF402.rds` and `real_T_percell_umi.csv.gz`
  (QC-derived).
- Review 1: `review1/L1/05_sim_removal.R`, `review1/L1/04_sim_k5_stacked.R`, `review1/L2/t2_k6b_sim.R`,
  `review1/L2/t6_k6b_sim_p3c.R`.

Decision 9.8: copy the scripts (not the data) into `scripts/ws/g07_sim/` with the code commit, with paths made relative.
The md5 list is frozen with the tag.

**8.5 Smoke runs.**
- Synthetic counts and cells on real IDs, with planted effects.
- No real covariate or status is written next to a smoke statistic.
- Outputs go under `/LARGE1/.../ws_a1_g07/smoke/` and are deleted after the real run reproduces.

**8.6 Compute.** About 60 T refits (K2b, K2c; 7-12 s each) and 160 stacked refits (K5; about 22 s each), plus the
UCell passes. Under 2 hours on one `c=16, m=64G` job.

**8.7 After the run.** Outcome section 11; GATES block; handoff S15 note.

## 9. Decisions for the user (recommended option first)

1. **Scope and timing.** Write and tag now, as specified. Alternative: wait for the late-November decision point.
2. **Primary test.** camera with estimated inter-gene correlation, FDR 0.05. Alternatives: camera 0.01 (2-19 % false set
   calls in simulation); fry as co-primary.
3. **Footprints and TF activity.** Not used. Alternative: record a descriptive use outside GATE_B0's verdict.
4. **T-cell programs.** Pin the 22 Gavish sets and interpret 19 (3 UNASSIGNED reported), with content-based naming.
   Alternatives: descriptive only; not used.
5. **Downloads, one batch.** Lasry 2023 tables (4.8 MB), Rutella 2022 tables (191 KB, CC BY 4.0), Mazziotta 2024
   supplement (2.3 MB; kept outside the repo).
   - Each gets a PROVENANCE note. They enter 3.6(b) only.
   - Options: approve now; at the decision point; not at all.
6. **S5b.** Run it for T as in 3.5: S5b-raw as registered, the S5b-1000 reading with the floor, gate and shape rule,
   and graded K6. Alternatives: S5b-raw only; disclose as never run.
7. **Codings.** Metadata sex primary; Control4 female and without Control1-5 as sensitivities.
8. **Simulation code.** Copied into the repo with the code commit.
9. **ADT (option c).** Deferred.
10. **G10(a).** A limitation now (1.5 item 12).
11. **The name "imprint".** Decided at writing.
12. **Time stamp.** An internal git tag.
13. **Consequence rule and generic-label wording** (7.3, 7.4): fixed now as written. Alternative: decide after the
    result. This binds the paper's wording.
14. **Per-compartment delay lists.** Derived from the on-disk Table S4 with the T rule (2.3), with no download.
    Alternative: T list only.

Real choices: 9.3-9.6 and 9.13. 9.1, 9.9, 9.11, 9.12 and 9.14 need confirmation. The rest are defaults.

## 10. Review history

- **2026-11-09, inventories.** A, B, C, D1, D2 and a completeness critic (workflow wf_309ce460-130). The critic checked
  44 citations and found no blindness breach.
- **2026-11-09, DRAFT v1** (c8b4f8f).
- **2026-11-09, review 1.** Workflow wf_90208e9a-ca1 was interrupted by a session end and finished as four agents.
  Outputs are in `review1/`.

  | lens | findings | blocking |
  |---|---|---|
  | L1 gene-set statistics | 11 | 1 |
  | L2 S5b | 13 | 1 |
  | L3 conformance | 39 | 4 |
  | L4 interpretation | 22 | 2 |

  Every blindness check held. Each blocking finding is fixed in v2:
  - removal judged on each set's own p and a direction-signed shift (L1 F2);
  - the per-cell doublet location dropped (L2-09);
  - G01's full sentence quoted (L3-01);
  - the 32-gene rule restated verbatim (L3-02);
  - the prior-looks ledger completed inline (L3-03);
  - 7.7 made verbatim (L3-04);
  - the handling sentence corrected (L4-01);
  - the 15-gene floor replaced by direction-matched removal with a 10-gene floor (L4-02).

  Main should-level changes:
  - S5b: thinning to exactly 1,000 UMI, a detection floor, graded K6 bars, an exact-match plant, nested
    leave-one-out, a down-set gate, a fifth cohort reading;
  - two-reference consensus tags with a composition-only expectation;
  - an opposite-sign gene class and per-compartment delay lists (L4-07);
  - ρ_X with a bootstrap interval on paired donors;
  - plants one per refit on fixed rows, plus real-fit power rows;
  - added limitations: acute illness, the source precondition, cryopreservation;
  - the never-written list completed;
  - decision 9.13 added.
