# Figure legends — TP53-mut vs TP53-WT AML patch

Figures carry data only; all explanation lives here. Files are in `results/figures/` (PNG + PDF),
source tables in `results/tables/`.

## Design common to every figure

**Groups.** TP53-mut = 9 samples from 6 patients: 6 samples with a single-cell or literature-confirmed
*TP53* mutation (tier A) and 3 with Numbat allele-based 17p LOH (tier B). TP53-WT = 9 samples, each
matched 1:1 to one mutant sample **from the same dataset and the same timepoint**, nearest in cell
count, distinct patient within a timepoint. Study and timepoint are therefore eliminated by design
rather than modelled. Seed 491638; the matching is reproducible.

**Effect measure.** "In how many of the 9 matched pairs is TP53-mut higher." This is the quantity the
test operates on. The cell-pooled odds ratio is deliberately *not* used for any claim: it is
cell-weighted, so one 15,148-cell control sample can set it on its own, and it disagrees in sign with
the paired measure for 51 of 145 panel genes.

**Test.** One-sided exact paired signed-rank (direction fixed a priori: the hypothesis is *more*
fibrosis/inflammation in the mutant arm). Exact null by dynamic programming, no Monte Carlo error.
At 9 pairs the smallest attainable p is 1/512 = 0.002; at 8 pairs 1/256; at 6 pairs 1/64.

**Robustness ("share of control sets").** 500 alternative matchings were drawn at random under the
identical constraints, and the test re-run in each. The reported share is the fraction of those valid
control sets in which the feature reached p < 0.05. It answers "would another analyst, matching
equally correctly, have seen this?" — it is **not** a p-value and must not be quoted as one. 5% is
the chance level.

**Confounding checked.** Sequencing depth does not explain the results and runs the conservative way:
the mutant arm is deeper in only 2 of 9 pairs (median 1,020 vs 1,045 genes/cell; one-sided p = 0.951
for "mutant deeper").

---

## P1 — Samples used and the evidence that defined the groups
`P1_cohort_evidence.png` · `P1_cohort_evidence.csv`

All 18 samples, mutant arm above, wild-type below, carrying **only the evidence that actually defined
the grouping**: TP53 mutant cells read directly from the reads, Numbat 17p LOH fraction, and the cell
count the matching used. Each row names the evidence that made that sample mutant. Each column is
colour-scaled independently; printed values are raw.

Mutant arm = 6 samples with the mutation read directly (5 from single-cell reads, 809653 from the
study's targeted sequencing — hence its blank genotype-cell count) + 3 with Numbat 17p LOH
(1886_R 0.981, 6323_R 0.931, MLL_28830 0.986). Wild-type = neither.

An earlier version of this figure also carried two inferCNV columns. That was a presentation error:
inferCNV was tested and **rejected**, so showing it beside the evidence that was used implied it had
contributed. It now has its own figure.

## P1B — The inferCNV check, which failed
`P1B_infercnv_check.png` · `P1B_infercnv_check.csv`

The PI's first question was "use inferCNV to split TP53-mut from TP53-WT". This is the answer, and it
is negative, so it is reported rather than dropped.

Three samples have both a genotype and an inferCNV 17p call. **inferCNV got none of them right:**

| sample | genotype truth | inferCNV call | verdict |
|---|---|---|---|
| AML916-D0 (C238Y) | mutant | 17p negative | missed |
| 809653 (E286G) | mutant | 17p negative | missed |
| AML420B-D14 | wild-type | **17p positive** (0.964, 11 altered arms) | false alarm |

Sensitivity 0/2, specificity 0/1 (`07_proxy_vs_genotype.csv`). n = 3 is a very small validation set —
that limitation is real — but it is 0-for-3 in the direction that matters, and the mechanism is not
in doubt: all five confirmed variants (Q144P, P152R, R273L, C238Y, E286G) are DNA-binding-domain
**missense** mutations. They change one base and leave copy number untouched, so an expression-based
CNV caller cannot see them by construction.

Two further validation results point the same way (`02_validation_all_datasets.csv`): the
complex-karyotype track fails outright (healthy-donor false-positive rate 15/26 = 58%, AUC 0.473 =
chance), while the 17p track has clean specificity in healthy donors (0/26) but no testable positive
control and is underpowered.

## P2 — Myelofibrosis-related genes
`P2_myelofibrosis_genes.png` · `.csv`

Collagen/reticulin structural genes, ECM-remodelling enzymes, and the megakaryocyte axis. Dashed line
at 4.5/9 = no difference. Point size = robustness share; orange = p < 0.05 uncorrected.

**Result.** The signal is in the cross-linking machinery, not in collagen itself. Collagen structural
transcripts are flat or lower in the mutant arm (COL1A1 4/9, COL1A2 2/9, COL3A1 2/9, DCN 2/9,
POSTN 1/9), while the enzymes that build and cross-link collagen are higher — PLOD2 7/9 (p < 0.05),
LOX 6/9 (p < 0.05), MMP14 7/9, TIMP1 6/9 — as is the megakaryocyte marker GP9 8/9 (p < 0.05).

## P3 — Cytokines, growth factors, ECM components
`P3_cytokine_ecm_genes.png` · `.csv`

Same axes as P2, for the profibrotic cytokine/growth-factor axis, inflammatory cytokines, niche
factors, hypoxia and EMT gene sets. The profibrotic effectors (SPP1, OSM, PDGFB, PDGFRB, IL11) trend
upward; SLC2A1 (hypoxia) is among the most robust genes in the whole panel.

## P4 — Pathway activity
`P4_pathway_activity.png` · `P4_pathway_scores.csv`

Per-cell UCell scores averaged per sample, one line per matched pair, for TGF-β (3 sets),
ECM organisation / receptor interaction (3), EMT, hypoxia, and inflammatory response (3).

**Result: none of them differs, and this is a firm negative.** Across 500 alternative valid control
sets the best of the 14 reached p < 0.05 in 5.8% — exactly the chance rate — and 8 of 14 never
reached it in any draw. Re-picking controls cannot rescue it.

**Why gene-level hits coexist with pathway-level nulls: gene-set membership, not biology.** Six of the
13 significant genes (IL11, C1QB, C1QC, LYVE1, MAF, MRC1) belong to **none** of these 14 sets, and all
three TGF-β sets contain **none** of the significant genes. A 200-gene rank score cannot register a
5-gene shift. The conclusion is that these generic MSigDB sets are the wrong instrument for this
question, not that the biology is absent.

## P5 — Cell-cell communication
`P5_ccc_nodes_and_edges.png` · `P5_ccc_node_inventory.csv`

**P5a, node inventory.** The PI asked for leukaemia / HSC / macrophage / niche / MSC / adipocyte.
Present: leukaemia (malignant) 3,774 mut vs 7,416 WT; HSC/MPP 2,628 vs 2,216; Mono_DC 14,576 vs 5,446;
macrophage-like 72 vs 167. **Absent: MSC/fibroblast 5 vs 3, adipocyte 0 vs 0, endothelial 0 vs 0.**

This absence is a measurement, not an omission — see F5. Whole-marrow aspirates are not collagenase-
digested, so stroma is never released from bone; only physically enriched libraries contain it.
**No TP53-mut sample has a usable stromal population, so the leukaemia–stroma interaction question
has n = 0 on the mutant side and cannot be tested at any effect size.**

**P5b, edges among the compartments that do exist.** CellChat significant communication probability,
relative normalisation, 8 usable pairs (AML328-D0 lacks a CellChat object on the mutant side, so pair
1 is dropped; 3853_R was re-matched to 1216_R). **No edge survives FDR correction**; the best is
LMPP_GMP → Mono_DC at p = 0.031 (absolute) / 0.063 (relative), FDR ≈ 0.94 over 49 edges. A node must be
present in *both* members of a pair, which leaves most edges only 4–5 usable pairs, and
HSC_MPP → HSC_MPP only 4. Macrophages were not given their own node: 2,269 macrophage-like cells exist
across *all* AML samples.

## P6 — Nectin-4 (PVRL4) and the nectin/PVR axis
`P6_nectin4.png` · `P6_nectin_per_sample.csv` · `P6_nectin_axis_tests.csv`

**P6a**, per-sample detection, log scale, dashed line at 1% of cells. **P6b**, the same pair by pair,
with pairs-higher and the one-sided paired p printed per gene.

| gene | pairs higher in mut | p | median mut | median WT |
|---|---|---|---|---|
| NECTIN4 | 6/9 | 0.15 | 0.20% | 0.11% |
| NECTIN1 | 4/9 | 0.82 | 0.48% | 0.72% |
| NECTIN3 | 3/9 | 0.85 | 0.43% | 0.42% |
| NECTIN2 | 6/9 | 0.41 | 2.16% | 2.64% |
| PVR | 4/9 | 0.88 | 1.68% | 3.00% |
| CD226 | 5/9 | 0.33 | 3.50% | 1.14% |
| **TIGIT** | **6/9** | **0.10** | **3.35%** | **0.48%** |

**NECTIN4 sits at the detection floor in both arms.** The answer to the PI's question is that 10x 3′
marrow scRNA-seq *cannot measure* Nectin-4 — not that Nectin-4 is absent. Their IHC/protein finding is
not contradicted here.

**Receptor side, reported as requested even though not significant.** TIGIT is the strongest signal in
the axis: higher in 6 of 9 pairs with a ~7-fold median difference (3.35% vs 0.48%), p = 0.10. CD226 is
higher in 5 of 9 (3.50% vs 1.14%, p = 0.33). Both are receptors for PVR/CD155 and NECTIN2. Neither
clears significance and the robustness shares are low (TIGIT 4.8%, CD226 0.4% of control sets), so
these belong in the discussion as an observation to follow up, **not in a results claim**.

---

## Supporting figures

**F1 `F1_effect_robustness`** — the panel genes: x = median paired difference in log2 CPM,
y = robustness share. The cloud at the origin is the null result; genes leaving it are the candidates.
x is the *sample-level* effect on purpose — an earlier version plotted the cell-level odds ratio
against a sample-level y, which put genes such as VSIG4 on the wrong side of zero.

**F2 `F2_robustness_decay`** — the honest figure. Each line is a gene; left point counts samples as
the unit (so AML328's four timepoints count four times), right point counts patients. **Both panels
fall left to right**: part of the sample-level result was AML328 counted repeatedly. Adding tier C
(right facet, 11 patients) does not repair it — tier C is expression-CNV only, and expression-CNV scored 0/3 against genotype, so
most of what it adds is expected to be mislabelled, and the lines fall further. At 6 pairs the p-floor
is 1/64, so some of the fall is power rather than error; the two cannot be fully separated here.

**F3 `F3_pathway_null`** — all 14 pathway scores across all four designs against the 5% chance line.
No pathway is above chance in *every* panel: TNFA_SIGNALING runs 0.2% → 7.8% → 3.0% → 35.4% across
designs A/B/C/D, hypoxia 0.2% → 12.6% → 0.6% → 27.0%. A result that reverses with the analysis unit is
not a result.

**F4 `F4_macrophage_shift`** — the two effect measures side by side, because they disagree and the
disagreement is the lesson. VSIG4 and CD163 read as strongly *lower* when cells are pooled (OR 0.23,
0.28) yet are higher in 6 of 9 pairs. On the measure that matches the test, macrophage markers trend
**up broadly** (MAF 8/9, C1QB 7/9, MRC1 7/9) with STAB1 the one exception (2/9). **This is not a
polarisation switch** — an earlier reading of this figure claimed one and was wrong. Macrophage
frequency itself does not differ (0.36% vs 0.46%, p = 0.32).

**F5 `F5_stromal_reality`** — stromal content by library type, log scale. Stroma-enriched libraries
27.8%; whole-marrow aspirates 0.028%; CD34-sorted libraries 0.0085%. The threshold for calling a cell
stromal was calibrated so that CD34-sorted libraries — which contain no stroma by construction — yield
0.85 false positives per 10,000 cells (dashed line). Whole-marrow AML sits at 3× that floor: 56 stromal
cells in 200,096 cells across 36 samples. Calls require per-cell **co-expression** of a fibroblast/MSC
module with no haematopoietic module, so ambient RNA cannot produce them: 1886_R's 528 projection-
labelled "Stromal" cells yield 0 confirmed, as do GSE185991 M95 (132) and M93 (50).

**G1 `G1_ccc_map`** — 7×7 signalling map, median share of each sample's total significant
communication. Magnitude panels use a sequential ramp, the difference panel a diverging one centred at
zero; the hues differ so that blue never means "high" and "negative" on the same page.

**G2 `G2_focus_pairs`** — every pair plotted for HSC/Mono_DC signalling share. Crossing lines are what
a non-significant result looks like; plotted rather than summarised so that a single driving pair
would be visible.

**G3 `G3_ccc_pathways`** — pathway-level signalling into and out of HSC and Mono_DC. Labels mark where
to look next, not findings: nothing survives FDR.

**G4 `G4_hsc_mono_edges`** — the four HSC/Mono_DC edges, pair by pair.

---

## What may and may not be claimed

| Claim | Status |
|---|---|
| C1QB and macrophage markers higher in TP53-mut | Preliminary observation. Most robust single result (90% of control sets), but weakens to 13–32% at patient level. |
| SPP1 / OSM / SLC2A1 / LOX / PLOD2 / GP9 trend higher | Preliminary observation, same caveat. |
| Collagen transcripts themselves higher | **Not supported** — flat or lower. |
| TGF-β / ECM / EMT / hypoxia / inflammatory **pathway activity** higher | **Not supported**, and shown not to be rescuable by re-picking controls. |
| Nectin-4 differs between groups | **Cannot be measured** in 10x 3′ marrow scRNA-seq. |
| TIGIT / CD226 higher | Observation only (p = 0.10, 0.33). Discussion, not results. |
| Leukaemia ↔ stroma / MSC / adipocyte interaction | **Cannot be asked** — n = 0 on the mutant side. |
| CCC differences between haematopoietic compartments | **Underpowered** — most edges have 4–5 usable pairs; nothing survives FDR. |

**The binding constraint throughout is 6 TP53-aberrant patients**, one of whom (AML328) supplies 4 of
the 9 mutant samples. Materially advancing this requires more public datasets carrying TP53 genotype,
or in-house samples.

---

## Q1 / V1 / V2 / V3 / U1 — the expression views

Added after a review found the panel set had a single chart form for gene expression (the paired
"n of 9 pairs higher" lollipop) and no sample-quality figure at all. A paired count is the statistic
the test operates on; it is not a view of the expression.

**Q1 `Q1_sample_qc`** · `Q1_sample_qc.csv` — library quality of the 18 samples.

**Metrics are computed on the cells that entered the comparison.** A first version read
`00_ingest/00_MASTER_qc_summary.csv`, which is the *pre-filter* summary: it reported 21,731 cells for
809653 where 12,340 entered, and 84,603 across the 18 samples where 49,519 did. A QC figure
describing a different cell set from every other figure in the project is worse than none. The
ingest counts are retained, but as the "before" half of Q1c. Post-filter totals now agree exactly
with U1 (49,519).

**Q1a** pairs every metric so a systematic arm difference would show as parallel lines. **Q1b**
genes per cell per sample, median with min–max whiskers. **Q1c** cells at ingest vs cells retained,
per sample, with the retention rate. Q1a's y is square-root scaled — the haemoglobin facet contains
exact zeros, which log cannot take.

**The arms are comparable, measured rather than eyeballed.** Paired signed-rank on each metric,
mutant higher in *n* of 9 pairs:

| metric | mut higher | p |
|---|---|---|
| cells per sample | 4 / 9 | 0.203 |
| genes per cell | 2 / 9 | 0.098 |
| UMIs per cell | 4 / 9 | 0.426 |
| mitochondrial % | 3 / 9 | 0.164 |
| ribosomal % | 7 / 9 | 0.098 |
| haemoglobin % | 3 / 9 | 0.855 |

Nothing reaches p < 0.05. The two closest are worth stating rather than rounding to "comparable":
the mutant arm trends toward **fewer genes per cell** (2/9, p = 0.098) and **more ribosomal
reads** (7/9, p = 0.098). Neither is significant, and the first one matters for how V3 is read —
lower detection sensitivity in the mutant arm biases *against* the V3 finding that detection rates
rise there, so that finding is observed against the bias, not with it.

**Overall retention 59%** (49,519 / 84,603), range 40–81% per sample, with no arm pattern.

**V1 `V1_volcano`** — the 145 panel genes with a testable paired result (of 148 in the panel; 19 significant up, 12 down). x is the **median paired log2 fold change** across the 9
pairs, computed per pair from per-sample mean expression and then medianed. It is *not*
`median_delta_log2cpm`, which is a difference of mean log-normalised values and runs 0.001–0.05 for
the lowly-detected genes here — plotting that put almost every gene on the y axis. y is
−log10(p) made two-sided from the one-sided paired test (the test asks only "higher in mutant", so a
down gene scores p≈1, not a small p). Horizontal banding is real, not a rendering artefact: 9 pairs
admit only a handful of attainable p-values, the smallest being 1/512.

**The dashed line sits at two-sided 0.10, which is one-sided 0.05.** Points are coloured by the
project's one-sided p — the quantity every other figure and table reports — while y is the
two-sided mirror. A line drawn at two-sided 0.05 therefore sat above several coloured points and the
figure contradicted its own legend. Moving the line to the equivalent level makes colour and
geometry agree; the subtitle states which p the line is.

**V2 `V2_dotplot`** — genes × (compartment × arm); dot size = % of cells with a detected transcript,
colour = mean expression **scaled within each gene**. Within-gene scaling is necessary: absolute
expression spans three orders of magnitude across this panel and a shared scale renders every
low-expressed gene blank. Size is on a square-root scale for the same reason — most panel genes are
detected in 1–8% of cells.

**V3 `V3_violin`** — per-cell distributions. **Violins cover expressing cells only**, with the
detection rate over *all* cells printed above each. Drawn over all cells the figure was honest and
unreadable: 14 of the 16 genes are detected in 0.1–7% of cells, so every panel collapsed to a spike
at zero and only PGK1 and VIM had a visible shape.

**The two numbers together carry the finding.** The detection rate roughly doubles-to-quintuples in
the mutant arm — C1QB 0.1%→0.4%, SLC2A1 2.9%→7.3%, OSM 2.6%→5.5%, IL11 0.4%→2.0%, GP9 0.2%→0.7%,
PLOD2 0.1%→0.5% — while the level *among cells that already express* barely moves. **The mutant
effect is more cells switching these genes on, not expressing cells expressing more.** The two
well-detected controls behave as controls: PGK1 56.4%→55.8% and VIM 89.0%→81.5% do not move.

**U1 `U1_umap`** — 49,519 cells from the 18 samples, 20,438 genes shared across all of them.
CP10K + log1p → 2,000 HVGs → 30 PCs → **Harmony on dataset** → UMAP. Harmony is not optional: the
samples come from five studies with different chemistries, and an uncorrected embedding separates by
study rather than biology. Correction is on *dataset only* — correcting on arm would erase the
contrast being examined. **U1a** cell type (malignant blasts pulled out of their haematopoietic bin,
8-slot categorical cap with the remainder folded into "other"), **U1b** TP53 group, **U1c** four
genes painted on, expressing cells drawn last so they are never buried.

**Caveat that must travel with U1b.** Harmony corrects *study*, not *patient*. The large right-hand
island is close to arm-pure, and that reflects one 15,000-cell sample rather than a TP53 effect.
U1b shows that the arms otherwise co-occupy the same territory — it is not evidence of a global
transcriptional difference, and no claim in this project rests on it.

---

## C1 / C2 / C3 — cell–cell communication, macrophage and HSC axes

Produced by `20_ccc_macrophage.R`. These extend P5 and G1–G4 in two directions the earlier figures
did not cover: the macrophage node the brief asked for, and the ligand–receptor resolution under the
compartment-level edges.

**C1 `C1_macrophage_node`** · `20_macrophage_node_availability.csv`
Per matched pair, cells available for a macrophage node, log scale, dashed line at CellChat's
10-cell floor. A pair counts only when **both** dots clear it, because the test is paired.

- **Macrophage-like (strict call): 1 of 8 pairs qualifies.** Counts run 0–158, median 1.
- **Mono_DC compartment: 5 of 8 pairs qualify.**

This is the answer to "run CCC on macrophages", and it is a negative one: a strict macrophage node
is not estimable in this cohort. What *is* estimable is the Mono_DC compartment that contains the
macrophages, and every macrophage statement below is at that resolution. The figure exists so the
limitation is visible rather than asserted in a footnote — "we skipped it" and "we measured that it
cannot be done" are different claims and only the second one is true here.

**C3 `C3_ccc_axis_totals`** · `20_ccc_axis_totals.csv`
Total significant communication on each axis as a share of the sample's total, bar = median across
qualifying pairs, one line per pair. Four tests, one-sided exact paired signed-rank:

| axis | direction | pairs | higher in mut | median Δ | p |
|---|---|---|---|---|---|
| Macrophage | **receives** | 5 | **5 / 5** | **+7.3 pp** | **0.031** |
| Macrophage | sends | 5 | 1 / 5 | −3.4 pp | 0.906 |
| HSC | receives | 4 | 1 / 4 | −2.5 pp | 0.938 |
| HSC | sends | 4 | 2 / 4 | −1.3 pp | 0.812 |

**What may be claimed.** In every one of the 5 qualifying pairs, the monocyte/macrophage compartment
of the TP53-mutant marrow receives a larger share of the marrow's signalling. p = 0.031 is the exact
floor at 5 pairs — the result is as strong as this design can produce, and no stronger. Across the
four tests two-sided BH gives FDR = 0.25, so **this is a consistent direction, not a corrected
significant finding.**

**Two structural checks, both of which the result survives.**

1. *Structural zeros.* A sample whose compartment holds fewer than 10 cells gets ~zero communication
   by construction, not because signalling fell. A first run left those rows in and returned 7/7
   pairs, p = 0.008 — an artefact, because two of the seven "increases" were WT samples carrying 1
   and 9 Mono_DC cells. Those rows are excluded; the table above is post-exclusion.
2. *Compartment size.* CellChat's permutation p depends on group size, so a smaller node yields less
   measured signal. **The mutant node is the smaller one in 4 of the 5 macrophage pairs** (14 vs
   1,225; 16 vs 554; 1,427 vs 12,058; 213 vs 833). The "receives" result is therefore observed
   *against* that bias. The mutant Mono_DC fraction is also *lower* in 3 of 5 pairs while receives
   is higher in 5 of 5.

**What may not be claimed.** "Mutant macrophages send less" is exactly what a smaller node would
manufacture on its own, and it is not separable from that here. Relative normalisation also forces
the shares to compete, so "receives up" and "sends down" are not independent observations.

**C2 `C2_ccc_ligand_receptor`** · `20_ccc_lr_tests.csv`
Every ligand–receptor pair on the two axes, paired mut vs WT, top 12 per panel by effect size; text
gives pairs-higher-in-mut over pairs measurable. Ligand names contain hyphens (HLA-DRA), so only the
spaced separator is rendered as an arrow.

**No single LR pair carries the C3 result.** The strongest on the macrophage-receives axis are
LGALS9→CD44 (4/4 pairs, p = 0.063) and the MHC-II→CD4 family (3/4, p = 0.125); nothing reaches
FDR < 1. The axis-level shift is distributed across many small contributions rather than driven by
one signal, which is what the 12 similarly-sized bars in that panel show. Read C2 as *where to look
next*, not as a list of findings.

**MSC / stromal axis.** Not computed, at any resolution. P5a and F5 give the reason: no TP53-mut
sample has a usable stromal population, so the node has n = 0 on the mutant side.
