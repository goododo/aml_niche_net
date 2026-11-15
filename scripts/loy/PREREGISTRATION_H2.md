# Pre-registration v1: within-patient -Y vs Y+ leukaemic program (LOY line, H2) (2026-11-14)

Status: v1, user approved the design and the file split 2026-11-15 ("其他都可以"). Frozen copy; to be tagged prereg-loy-h2-v1 after Stage 0 counts and before 05 reads expression. No expression value from any library below has been read for this question.
What has been read: per-cell chrY UMI counts, fine labels, depth, the -Y fraction per label and timepoint
(REPORT_step2_external_2026-11-14.md section 3), and the diagnostic karyotypes of the Lambo patients
(external_clinical/lambo2023_target_aml1031_clinical.csv). On approval this file is copied to
`scripts/loy/PREREGISTRATION_H2.md` and tagged `prereg-loy-h2-v1` before any script reads expression for this contrast.

## 1. Question

Within one AML library, do leukaemic cells that have lost chromosome Y differ in transcriptional program from
Y-retaining leukaemic cells of the same patient, same fine label and same sequencing depth?

- "Differs": the program travels with a genomic branch of the leukaemia (the -Y subclone), i.e. part of the
  within-patient program variation is clonal.
- "No difference": this genomic branch does not carry a program; the within-patient variation seen in the field
  (state, stress, environment) is not explained by it.

The patient, library, capture day, ambient pool and handling are the block; nothing cross-sample enters the primary
contrast. The free second readout is the fate of the -Y clone at the patient's other timepoint, already observed
(section 2), which lets a program difference be set against a selection event without new data.

What -Y buys and what it costs: a per-cell genomic label at no genotyping cost; in exchange the label is defined by
absence of ~12 expressed chrY genes, so (a) chrY genes are removed from every feature, (b) the label is also a
selection on dropout and needs the technical-zero control (section 4), and (c) the contrast is "-Y branch versus the
rest", not "effect of losing Y": the branch may carry other lesions.

## 2. Cases and eligibility (counts only; chosen after the counts, disclosed as such)

Readable cases from the step-2 census (frozen three-state call, `scripts/loy/config_loy.R` section 5):

| case | source | age | library with -Y subclone | f_hat LEUK | suspected / retained | other timepoint |
|---|---|---|---|---|---|---|
| AML15 | Lambo 2023 (Zeng pack), 10x 3' | 9.6 y | AML15_DX | 0.169 | ~870 / ~3,900 | REL 0.000 (3,700 leukaemic), REM 0 |
| AML13 | Lambo 2023 (Zeng pack), 10x 3' | 3.9 y | AML13_REL | 0.132 | ~600 / ~4,100 | DX 0.000 (1,878), REM 0.019 (387) |
| PT43 | GSE269669 Daver 2025, 10x 5' | adult (59-73) | PT43A (pre-treatment) | 0.288 | 968 / 2,488 | PT43B on-treatment 0.000 (3,319), non-responder |

Fine labels with enough cells in both halves before matching (deep cells, from the step-2 tables): AML15_DX Early GMP
(2,240; 0.23), GMP-Neut (846; 0.28), GMP-Mono (285), LMPP (295), GMP-Cycle (161; 0.39), Cycling Prog (99; 0.40);
AML13_REL Early GMP (2,329; 0.17), EoBasoMast (762), LMPP (443), GMP-Neut (323), MEP (299), MPP-MkEry (261), Cycling
Prog (184); PT43A GMP (1,739; 0.35), HSC (1,259; 0.26), BaEoMa (560), CD14 Mono (334), plus 938 Azimuth-"Stromal"
cells that are leukaemic by depth and context (treated as their own label, reported, not pooled with the rest).

Selection-only cases (too few cells for a within-library contrast; used in section 9 only): AML5 (Lambo; DX 0.016 of
353, REL 0.858 of 141: the relapse is the -Y clone) and Mumme Sample20 (D 0.950 of 2,352; E 51 cells).

Genetics (TARGET AAML1031 open clinical data, diagnostic ISCN; USI map in external_clinical/):
- AML15 = PAWWXU, male, t(8;21) CBF-AML: `46,XY,t(8;21)[8] / idem,del(9)(q13q22)[4] / idem,dup(1)(p31p36.3),del(9)[4] /
  45,idem,-Y[3] / 45,idem,-Y,del(9)[5]`. -Y in 8 of 24 metaphases (0.33) as a sub-branch of the t(8;21) stemline.
  Not in CR after course 1 (MRD 29 %), relapsed. This is the one case with an orthogonal, cell-based -Y measurement.
- AML13 = PAWNPG, male, t(8;21)[12]/46,XY[3]; no -Y at diagnosis; CR, relapsed. The -Y clone is relapse-only.
- AML5 = PAWUAE, male, KMT2A-r (t(10;11)) with add(8)/+8 and del(9q) subclones; no -Y at diagnosis.
- PT43: TP53-mutant AML, non-responder (author metadata). Per-patient variant and karyotype are in Daver 2025
  Table S5, which could not be retrieved programmatically (PMC bot wall, AACR 403). Lookup is a Stage 0 item; it
  does not block the paediatric cases.

Eligibility rule (applied on counts only): per library, after depth matching (section 4), at least 50 suspected_LOY
and 50 Y_retained leukaemic cells inside at least one fine label; a technical-zero control compartment in the same
library with at least 50 Y-zero and 50 Y+ cells after matching (T/NK, or any normal-labelled compartment with
f_hat <= 0.03). Fewer than 2 eligible cases: STOP.

Order: AML15_DX first (paediatric, cytogenetic known answer), then AML13_REL, then PT43A. Adults carry the mosaic-LOY
confound shown by PT35B (step-2 report 3.2); PT43A enters the primary only if GATE_H2_CNV passes (section 6).

## 3. Per-cell labels and misclassification

- Labels are the frozen three-state call: suspected_LOY = chrY UMI 0 and lambda >= 4.605 (P(0 | intact) <= 0.01);
  Y_retained = (chrY UMI >= 2 on >= 2 genes) or >= 3; uncertain cells dropped. Lambda uses the step-2 per-study
  calibration (config 7b). No threshold is changed for this design.
- Y+ cells mislabelled -Y: <= 1 % by construction; measured per library by the Y-zero rate of the same library's deep
  T cells (PT43A 0.018; Lambo <= 0.02). -Y cells mislabelled Y+ (ambient chrY UMIs): bounded by the Y_retained
  fraction in clonal -Y leukaemias at similar depth (AML5_REL 0.14, Mumme 20D 0.05, GSM8791436 0.004); reported per
  case as c_amb. Both errors shrink a true difference toward zero; neither creates one.
- Fine labels: the pack's BoneMarrowMap `predicted_CellType` for Lambo; Azimuth `predicted.celltype.l2` for PT43.
  Labels are projections, not clusters, and do not use chrY (checked: no chrY gene in either reference's markers is
  assumed; verified in Stage 0 by re-running the label with chrY genes removed on a 2,000-cell subsample; agreement
  >= 0.95 required, else labels are coarsened to the LEUK bin).

## 4. Confounds and matching

1. Depth: within each fine label, each suspected_LOY cell is matched to up to 3 Y_retained cells by nearest log10
   nCount (caliper 0.1); unmatched cells dropped; the matched set is the analysis set. log(nUMI) stays in the model.
2. Fine label: stratum in the matching and a covariate in the model; per-label effects reported alongside.
3. Cell cycle: a primary feature (F4), not a covariate; a second model with the S/G2M score as covariate is a
   sensitivity arm, because AML15_DX showed a higher -Y fraction in cycling labels (0.39-0.40 vs 0.23).
4. chrY genes (the label genes plus every gene annotated to chrY in the gene order file) are removed from all
   feature sets and from the genome-wide list. XIST is not expressed in these male cells.
5. Dropout selection: a Y-zero cell is, among cells of equal depth, one that happened to capture fewer transcripts
   of a dozen low-expressed genes. The technical-zero control runs the identical pipeline on the same library's T/NK
   (or other normal) cells, where Y-zero is pure dropout. Any feature that "differs" there is a dropout artefact.
6. Doublets: cells flagged by the pack / author QC are already absent; a Y_retained cell with chrY UMIs from a
   second cell is a Y+ error of the kind bounded in section 3.

## 5. Features and statistic

Feature panel, frozen at the tag as `scripts/loy/gene_sets_h2_v1.tsv` (md5 in the tag message), assembled from files
already frozen in this project; a feature without a citable gene table at freeze time is dropped, not improvised:
- F1 LSC17 (Ng 2016), F2 van Galen malignant HSC/Prog-like, F3 HSPC core (`scripts/以前06_hierarchy/stemness_signatures.tsv`,
  as frozen in PREREGISTRATION_CLONAL_ENV.md).
- F4 cell cycle: Seurat S and G2M lists, reported as S, G2M and cycling fraction.
- F5 HALLMARK interferon-alpha, F6 interferon-gamma, F7 TNF-alpha via NF-kB (overlaps the handling program MP9;
  flagged), F8 oxidative phosphorylation, F9 MYC targets V1, F10 P53 pathway (relevant to PT43 only; reported for all).
  HALLMARK sets from the pinned local GMT `/LARGE1/gr10634/gaozy/reference/h.all.v2024.1.Hs.symbols.gmt` (as in config_ws.R), never fetched at run time.
- F11 this project's MP5 (HLA-II), MP6 (non-classical / IFN), MP10 (inflammatory secretory) and MP9 (dissociation
  stress; artefact control) from `results/tables/04_cnmf/malignant_metaprograms.tsv`.
- F12 chrY gene set: machinery check only (GATE_H2_K2), never a result.

Score: per library, log-normalised expression; genes detected in >= 5 % of the library's leukaemic cells (set before
the label is looked at); z-scored per gene over the matched set; mean over the set's genes.
Model, per library and feature, on the matched set: score ~ label + fine_label + log(nUMI). Effect d = label
coefficient / residual SD, then divided by (1 - c_amb) for the ambient mislabelling bound.
Per-feature p: 1,000 permutations of the label within fine label x depth bin (10 bins). Holm over F1-F11 per case.
No pooling across cases (three cases, three genetic backgrounds); the result is a cases x features table.

Genome-wide DE (secondary, S1): Wilcoxon on the matched set per gene, chrY excluded; the p-value distribution of the
technical-zero control in the same library is the empirical null; genes beyond it listed descriptively.

## 6. CNV co-segregation (is the -Y population a branch of the leukaemia?)

inferCNV with the library's own normal-labelled T/NK (+ B) cells as reference (>= 100 cells; otherwise the project's
external reference via `scripts/02_malignancy/00_infercnv_common.R`), run on the leukaemic cells of the library;
per-cell mean CNV signal in pre-specified regions:
- AML15_DX: 9q13-q22 (del in 8/16 Y+ abnormal and 5/8 -Y metaphases: expected in both halves, not discriminating) and
  1p31-p36.3 (dup in 4/16 Y+, 0/8 -Y metaphases: expected absent in -Y cells, present in a Y+ subset). The karyotype
  therefore predicts one region that separates and one that does not; both are reported.
- AML13_REL: no relapse karyotype available; descriptive scan of every arm for co-segregation with -Y (difference of
  per-cell CNV signal between halves beyond the T-cell null).
- PT43A: TP53-mutant AML; the usual lesions (5q, 7, 17p) plus the descriptive arm scan.
GATE_H2_CNV (PT43A only, adult): at least one leukaemic CNV present in -Y cells at a frequency within [0.7, 1.3] of
its frequency in Y+ cells. Pass: -Y is a branch of the leukaemia. Fail: PT43A is reclassified as "-Y clone of
uncertain lineage" (the PT35B pattern) and leaves the primary; its contrast is still reported as such.
The paediatric cases do not need the gate (no mosaic LOY at 4-10 y); their CNV readout is reported.

## 7. Gates (in order; any STOP is a stop for discussion)

- GATE_H2_0 (counts): >= 2 cases pass section 2 after matching. Else STOP.
- GATE_H2_K1 (technical-zero control): in every library, the full pipeline on the control compartment yields no
  feature with Holm p < 0.05 and |d| > 0.2. A feature that fails is removed from that case's primary readout and the
  failure recorded; >= 3 failing features in one library: STOP (the pipeline is reading dropout).
- GATE_H2_K2 (machinery): F12 chrY score is lower in suspected_LOY cells, permutation p < 0.001, in every case.
- GATE_H2_K3 (label independence): fine labels recomputed without chrY genes agree >= 0.95 (section 3).
- GATE_N1 (shuffle): 200 label shuffles within fine label x depth bin; per-feature rejection rate at alpha 0.05 in
  0.03-0.07.
- GATE_P (power, from the smoke run): d = 0.3 planted into a random 60-gene set in -Y cells of the real matched set
  is detected (Holm p < 0.05) in >= 80 % of 200 plants, per case. A case failing GATE_P can only yield "not
  readable" for a null feature.
- Consistency, not a gate (already observed): AML15_DX scRNA -Y fraction 0.17 (LEUK) to 0.31 (Mono_DC) versus 0.33 of
  diagnostic metaphases.

## 8. Primary result and its wording (fixed now)

Per case and feature, after Holm over F1-F11:
- "-Y branch differs": p < 0.05, |d| >= 0.2, and the sign of d is the same in every fine label with >= 50 / 50 cells.
- "No difference": 95 % CI of d inside (-0.2, 0.2) and GATE_P passed for that case.
- Otherwise "not resolved".
Across cases: the 3 x 11 table is the result. Two or more cases agreeing on a feature is reported as agreement, not
as a pooled estimate. The selection readout (AML15 lost its -Y clone at relapse; AML13 gained one; PT43 lost it under
aza-ven-magrolimab without responding; AML5's relapse is the -Y clone) was seen before this document and is
disclosed as such; programs are read without a hypothesis about direction, and the program-fate correspondence is
descriptive with n = 4 and no test.

## 9. Secondary (descriptive, no gate)

- S1 genome-wide DE per case (section 5), chrY excluded, control-calibrated.
- S2 -Y fraction by fine label and by cycling versus resting label (already seen for AML15; reported for all).
- S3 the depositors' own AUCell hierarchy scores (Zeng pack columns; Daver metadata) as a cross-check of F1-F3.
- S4 selection-only cases (AML5, Mumme Sample20): -Y fraction by timepoint and label; no contrast.
- S5 the denominator (h2_denominator_2026-11-14.csv): 50 unique male AML patients with >= 2 readable leukaemic
  timepoints (22 paediatric, 28 adult); 5 with a clonal or subclonal -Y class at one timepoint and a different class
  at another (PT43, AML13, AML15, AML5, Mumme 20). Reported as a table; no rate is claimed.
- S6 new cases from the second GEO batch (GSE315856, GSE223844, GSE343017, GSE198681; on disk 2026-11-14) and from
  ScPCA Gawad when downloaded: scanned with the frozen steps 1-2 first; a new subclonal case joins as a dated
  amendment with its results reported before and after.

## 10. Smoke before the real run

Two worlds on the real matched set of AML15_DX only: a null world (labels shuffled within fine label x depth bin)
that must traverse the "no difference" and "not resolved" branches and the control-compartment branch end to end,
and a planted world (d = 0.3) that must traverse "-Y branch differs". Outputs under `*_smoke` names, never next to
real covariates; deleted after the real run has reproduced, with parameters and results recorded in the report.

## 11. Known weaknesses (stated now)

- Three cases, two genetic backgrounds (t(8;21) x 2, TP53): nothing generalises; the result is per patient.
- The -Y branch may carry other lesions (AML15's -Y metaphases also carry del(9q) in 5 of 8); "differs" is a
  statement about the branch, not about chrY.
- Fine labels are projections onto references built from other marrows; label noise dilutes per-label effects.
- Lambo libraries are 3' v3 at 12-16k UMIs; PT43A is 5'; the two chemistries are never compared directly.
- The selection readout preceded this document.
- PT43's genotype is pending (Daver Table S5); the adult mosaic-LOY confound is handled by GATE_H2_CNV only.
- Ambient mislabelling (c_amb) is bounded from other samples, not measured in these libraries.

## 12. File split (for confirmation before any code is written)

| file | role | writes |
|---|---|---|
| `scripts/loy/PREREGISTRATION_H2.md` | frozen copy of this document; tag `prereg-loy-h2-v1` | - |
| `scripts/loy/gene_sets_h2_v1.tsv` | feature panel F1-F12, assembled from the files named in section 5; md5 in the tag | - |
| `scripts/loy/04_h2_stage0.R` | per-case counts, label-independence check (K3), depth matching, eligibility; the denominator table | `results/tables/loy/04_h2/stage0_counts.csv`, `stage0_matched_<case>.csv.gz`, `h2_denominator.csv` |
| `scripts/loy/05_h2_contrast.R` | features, model, permutations, K1 / K2 / N1 / P, `--world real|null|planted` | `results/tables/loy/04_h2/contrast_<case>.csv`, `control_<case>.csv`, `checks_<case>.csv` |
| `scripts/loy/06_h2_cnv.R` | within-library inferCNV via `00_infercnv_common.R::run_one`, region scores, GATE_H2_CNV | `results/tables/loy/04_h2/cnv_<case>.csv`; inferCNV dirs under LARGE1 |
| `scripts/loy/submit_loy_h2.sbatch` | array over cases, gr10634c, c=2 m=48G, OMP_PROC_BIND=false | logs `00_project/logs/loy_h2_%A_%a.{out,err}` |

Inputs: the step-2 per-cell tables (`LARGE1/loy/percell/<ds>__y_percell.csv.gz`) for labels and depth; counts from
the pack objects / GSE269669 objects via the step-2 adapter (`03_scan_external.R::counts_of`). No new dependency
(Seurat, infercnv, data.table already in use). Stage 0 runs before the tag only for counts; the tag is placed before
05 reads any expression value.

## 13. Amendments made before the freeze (2026-11-15, before any expression value is read)

1. **PT43 fine labels** (section 3): PT43A/B leukaemic and T cells are re-projected with this project's BoneMarrowMap
   pipeline (`04_h2_stage0.R` step [4]) so that all three cases share one label system; the authors' Azimuth l2 label is
   kept as `fine_depositor` (secondary, reported in the cross-table). Reason: the two label vocabularies are not
   comparable, and only the BMM projection can be rerun without chrY genes.
2. **GATE_H2_K3 procedure** (sections 3, 7): by construction when the scaffold's variable genes contain no chrY gene
   (checked and recorded); otherwise by re-projecting a 2,000-cell subsample of each Lambo contrast library with and
   without chrY genes. PT43 is always re-projected both ways.
3. **Technical-zero control** (section 4): the library's T cells, Y-zero at any depth versus Y_retained, depth-matched
   with the same procedure as the leukaemic contrast; stratified by T fine label and pooled.
4. **Shared code file** (section 12): `scripts/loy/00_loy_common.R` holds `counts_of()`, `read_gse269669()`,
   `loy_comp_of()`, `loy_cell_call()` and `loy_h2_is_leuk()` for 04-06. 02 and 03 are frozen and keep their own copies.
5. **Stage 0 outputs** (section 12): `stage0_k0_reproduction.csv` (K0: the re-derived call reproduces 02's LEUK counts
   for every H2 library; a mismatch stops the script), `stage0_cells_<case>.csv.gz` (the per-cell handoff to 05, with
   `match_id`), `stage0_pt43_bmm_percell.csv.gz`, `stage0_k3_label_independence.csv`, `stage0_gates.csv`.
6. **Daver Table S5** (section 2): the user downloads it by hand (PMID 40198272, PMC12165814; PMC lists 17 unnamed
   supplement PDFs); it is read from `reports/2026-11_loy/external_clinical/daver2025/` when present. Not blocking.
7. **CNV co-segregation operationalised** (section 6; `06_h2_cnv.R`, constants config_loy.R 8c): inferCNV on the
   library's leukaemic cells (suspected_LOY + Y_retained, one observation group so the clustering never sees the label)
   against the library's own normal T / NK / B cells (>= 100, >= 20 per lineage group) or, when the library has too few,
   the project's BoneMarrowMap external reference by lineage; project cutoff and denoising, HMM off and analysis mode
   "samples" (region means need the residual matrix, not HMM states). Per-cell region score = mean residual over the
   region's genes (>= 20 genes); regions = every chromosome arm from the UCSC hg38 cytoband table (fetched 2026-11-15
   into `/LARGE1/gr10634/gaozy/reference/hg38/`, provenance recorded) plus the pre-specified cytoband ranges per case.
   A cell carries a region's CNV when its score lies beyond 3 reference SDs in the direction of the Y+ leukaemic
   mean; a region is a leukaemic CNV when >= 30 % of Y_retained leukaemic cells carry it. GATE_H2_CNV (PT43A) passes
   when at least one leukaemic CNV has carrier frequency in -Y cells within [0.7, 1.3] of that in Y+ cells. chrY itself
   is scored as the machinery check of this arm (must separate the halves). Label permutation (1,000) within fine label
   gives the per-region p for the -Y vs Y+ difference of region scores.
8. **Stage 0 findings that changed the procedure, before the tag (2026-11-15, job 3689458, counts only)**:
   (a) the T cells of the contrast libraries hold too few Y-zero cells for the control floor (22 / 15 / 47 matched Y-zero
   cells in AML15_DX / AML13_REL / PT43A). The technical-zero control therefore uses every non-leukaemic cell of the
   compartments whose sample-level f_hat is <= 0.03 in 02's table (T, NK, B, erythroid, normal Mono/DC, megakaryocyte;
   as section 2 already allowed), and when the contrast library still gives < 50 matched Y-zero cells, the patient's
   other libraries join with library as a stratum (depth matching never crosses libraries; the library is a covariate
   in the control model). (b) Greedy 1:3 matching inside the 0.1 log10 caliper leaves median depth gaps of up to 9 %
   between the halves (PT43A HSC 0.91, AML15 Early GMP 1.07); the hard stop is set at half the caliper (0.05 log10) and
   log depth stays in the model as pre-registered. (c) GATE_H2_K3 is computed over cells that received a label in both
   projections; mapping-QC failures (1.2 % of PT43 cells) are reported, not counted as disagreement. (d) The scaffold's
   variable genes contain one chrY gene (PCDH11Y), so the Lambo subsample re-projection runs (section 13.2).
