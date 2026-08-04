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

## P1 — Samples used and the evidence behind the grouping
`P1_cohort_evidence.png` · `P1_cohort_evidence.csv`

All 18 samples, mutant arm above, wild-type below, with every number that entered the grouping:
genotyped mutant cells, Numbat 17p LOH fraction, inferCNV 17p fraction, inferCNV altered-arm count.
Each column is colour-scaled independently; printed values are raw.

**The point of this figure is the contrast between columns.** Not one of the nine TP53-mut samples is
called 17p-positive by inferCNV, while the *wild-type* sample AML420B-D14 is (17p = 0.964, 11 altered
arms). Numbat cleanly recovers 1886_R (0.981), 6323_R (0.931) and MLL_28830 (0.986). This is why
expression-CNV was never allowed to define the mutant group: against single-cell genotype truth its
sensitivity is 1/3, and all five confirmed variants (Q144P, P152R, R273L, C238Y, E286G) are
DNA-binding-domain **missense** mutations, which an expression-based caller cannot see by construction.

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

**F1 `F1_effect_robustness`** — all 160 panel genes: x = median paired difference in log2 CPM,
y = robustness share. The cloud at the origin is the null result; genes leaving it are the candidates.
x is the *sample-level* effect on purpose — an earlier version plotted the cell-level odds ratio
against a sample-level y, which put genes such as VSIG4 on the wrong side of zero.

**F2 `F2_robustness_decay`** — the honest figure. Each line is a gene; left point counts samples as
the unit (so AML328's four timepoints count four times), right point counts patients. **Both panels
fall left to right**: part of the sample-level result was AML328 counted repeatedly. Adding tier C
(right facet, 11 patients) does not repair it — tier C is expression-CNV only at 1/3 sensitivity, so
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
