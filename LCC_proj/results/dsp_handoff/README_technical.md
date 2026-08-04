# DSP handoff — technical README

Files produced by `LCC_proj/scripts/16_dsp_handoff.R`.
Plain-language summary for clinicians: `README_for_clinicians_CN.md` (same folder).

| file | what it is |
|---|---|
| `16_profile_matrix_AMLmarrow.tsv` | genes × cell types, CP10K units. The SpatialDecon `X` matrix. |
| `16_profile_celltype_meta.csv` | per column: contributing cells/samples/datasets, RNA-content scaler, honesty flags |
| `16_profile_qc_markers.csv` | canonical marker → highest-expressing column, and that column's share |
| `16_profile_collinearity.csv` | pairwise correlation of log profiles, descending |
| `16_genesets_LCC.gmt` | the six pre-specified gene sets |
| `16_genesets_LCC_evidence.csv` | every set member with the statistic it was selected on |

---

## 1. The profile matrix

### How it was built

- 220 AML / normal marrow samples across 13 studies (`whole_MNC` 153, `sorted_CD34` 46,
  `stroma_enriched` 21) that have both a per-cell label table and a QC'd counts object.
- Haematopoietic labels from the BoneMarrowMap projection (`hierarchy_bin` / `bmm_broad`);
  malignant calls from the main-line pipeline; **stromal labels called de novo per cell** by module
  co-expression (`stromal_modules.R`), not by projection label — the projected `Stromal` bin is
  unreliable in both directions and 1886_R's 528 projection-"Stromal" cells yield 0 confirmations.
- Each cell scaled to CP10K, averaged within (sample, cell type), then **median across samples**.
  The median rather than a pooled mean is deliberate: cell-weighted pooling lets one 15,148-cell
  sample set the answer alone, and in the main analysis it disagreed in sign with the sample-level
  measure for 51 of 145 genes.
- Gene universe: symbols present in the reference of ≥80% of contributing samples. A gene missing
  from a sample's reference is treated as *not measured there*, not as zero — the zeros are restored
  exactly at aggregation, so the median is exact, not approximated.
- Columns are rescaled to equal sums, i.e. the matrix is in **RNA-composition units**.

### Using it

```r
library(SpatialDecon)
X <- as.matrix(read.delim("16_profile_matrix_AMLmarrow.tsv", row.names = 1))
res <- spatialdecon(norm = your_Q3_normalised_AOI_matrix,
                    bg   = your_background_matrix,
                    X    = X[intersect(rownames(X), rownames(your_norm)), ])
```

**Abundance units.** Equal column sums means `res$beta` is RNA abundance, not cell abundance.
To convert, divide each cell type's estimate by `rna_content_rel` from the meta table —
an erythroid progenitor and a megakaryocyte differ several-fold in RNA content and both are in
these AOIs.

### What this reference cannot see — read this first

**Droplet scRNA-seq loses the two most mature marrow lineages, and both matter here.** Same root
cause: mature erythroblasts and mature granulocytes have very low RNA content and are removed by the
droplet prep and the QC gene-count filter. Numbers are max CP10K across all 16 columns:

| lineage | present | absent |
|---|---|---|
| granulocytic | promyelocyte/GMP: `MPO` 13.1, `ELANE` 4.4, `AZU1` 8.5, `PRTN3` 3.0 (all in `LMPP_GMP`) | myelocyte onward: `DEFA1B` **0**, `CEACAM8` **0**, `FCGR3B` 0.005, `CXCR2` 0.016, `DEFA4` 0.02, `CAMP` 0.03, `DEFA3` 0.09, `LTF` 0.11, `LCN2` 0.16 |
| erythroid | MEP / early progenitors | erythroblasts: `GYPA` 0.001, `ALAS2` 0.006, `AHSP` 0.08, `SLC4A1` 0.03; `HBB` peaks at 1.05 in `AML_blast`, i.e. ambient |

The genes in the "absent" column top essentially random rare-stromal columns — that is noise, not
signal, and should not be read as attribution.

**Direct consequences for the DSP analysis:**

1. Marrow trephine AOIs contain erythroid islands and mature granulocytes. This reference cannot
   represent them, so that signal will be pushed onto whichever column is nearest. Either accept the
   bias, or splice in erythroblast/neutrophil profiles from an external reference.
2. `DEFA1B` is a top-4 "higher in nonfibrotic" gene in two of the four contrasts, and it is
   **exactly zero in every column here** — deconvolution cannot attribute it. It *can* attribute
   `AZU1` / `PRTN3` / `ELANE` / `CTSG` / `MPO`, which is the bulk of that signal, all cleanly to
   `LMPP_GMP` (shares 0.55–0.74). So the composition hypothesis remains testable; just note that the
   mature-neutrophil part of it is out of reach on this reference.

`16_profile_qc_markers.csv` records both gaps explicitly (`granulocyte_mature` and `erythroid`
marker sets) so the limitation travels with the data rather than living only in this file.

### Read the meta table before trusting a column

`flags` is not decoration. Specifically:

- **`no_erythroblasts_see_README` on `Erythroid_MEP`** — the column holds MEP and early progenitors
  only; see the table above.
- **`single_dataset`** — marrow stroma survives only in physically enriched libraries, so several
  stromal columns necessarily come from one study. Their profile is that study's and cannot be
  cross-checked.
- **`few_cells`** (<1000) — `Endothelial`, `Pericyte`, and `Adipo_MSC` are built on a few hundred
  cells each.
- **`adipogenic_MSC_not_mature_adipocyte`** — `ADIPOQ` reaches only ~1.8 CP10K. Droplet scRNA-seq
  destroys mature marrow adipocytes (large, buoyant, over the droplet size limit). This column is
  the adipogenic / adipo-CAR MSC, which is what droplet data *can* see. A low estimate here is
  **not** evidence that the marrow lacked fat.

### Collinearity — this decides how many columns are usable

`16_profile_collinearity.csv`. Log-profile correlations reach ~0.99. High baseline correlation is
normal for scRNA reference matrices (shared housekeeping/ribosomal expression) and SpatialDecon is
built for it, but the top pairs should be **read jointly rather than as independent estimates**:

- `HSC_MPP` / `LMPP_GMP` / `Erythroid_MEP` — three progenitor columns, r ≈ 0.99
- `DC` / `Macrophage` — r ≈ 0.98; treat as one mononuclear-phagocyte readout unless the split is
  clean in your data
- `AML_blast` / `Monocyte` — r ≈ 0.98; expected, much AML is monocytic
- `Adipo_MSC` / `MSC_fibroblast` — r ≈ 0.98; not independently identifiable

With three AOI masks, do not expect all columns to be separately identifiable. The columns that
carry the questions in the deck and that pass marker QC cleanly are:
**`LMPP_GMP`, `Macrophage`, `Plasma`, `Megakaryocyte`, `MSC_fibroblast`, `Osteolineage`,
`AML_blast`, `T_NK`.**

### Marker QC — what it does and does not confirm

`16_profile_qc_markers.csv` gives each canonical marker's highest-expressing column and that
column's share of the total. Passing cleanly:

| set | result |
|---|---|
| granule / GMP (AZU1, PRTN3, ELANE, CTSG, MPO) | → `LMPP_GMP`, shares 0.55–0.74 |
| macrophage (C1QA/B/C) | → `Macrophage`, shares 0.93–0.99 |
| plasma (IGHG1, IGKC, MZB1, JCHAIN) | → `Plasma`, shares 0.80–0.90 |
| megakaryocyte (PF4, PPBP, GP9, ITGA2B) | → `Megakaryocyte`, shares 0.92–1.00 |
| T/NK (CD3D, CD3E, IL7R) | → `T_NK`, shares 0.93–0.98 |
| osteolineage (BGLAP, IBSP, SP7) | → `Osteolineage`, shares 0.97–1.00 |
| endothelial (CDH5, EMCN, VWF) | → `Endothelial`, shares 0.84–0.95 |
| fibroblast/MSC (DCN, LUM) | → `MSC_fibroblast`, shares 0.38–0.39 |
| ribosomal control (RPL13, RPS23, …) | uniform ~0.10 across columns, as it should be |

Apparent failures that are correct biology, not errors: `COL1A1`/`COL1A2` peak in `Osteolineage`
(osteoblasts are the dominant type-I collagen source), `CD34`/`CRHBP` peak in `Endothelial` (both
are genuinely endothelial as well as HSPC genes), and `CXCL12`/`COL3A1` peak in `Adipo_MSC`
(adipo-CAR cells are the CXCL12-abundant reticular population).

One caveat that matters for the gene sets below: **`LYVE1` peaks in `Pericyte` (share 0.78), not in
`Macrophage`.** LYVE1 is genuinely expressed by perivascular cells as well as by resident
macrophages. It is a member of `LCC_TP53_MACROPHAGE_C1Q`, so if that set moves in an AOI-level test,
check `C1QA/B/C` (shares 0.93–0.99, unambiguous) before attributing the change to macrophages.

**The first three rows are the point of the whole exercise.** They are the columns that test the
three things driving the DSP DE lists: whether `AZU1`/`PRTN3`/`ELANE` reflect GMP fraction, whether
the CD68+ signal is macrophage abundance, and whether `IGHG1`/`IGKC` is plasma-cell content.

---

## 2. The pre-specified gene sets

`16_genesets_LCC.gmt`, six sets, 6–8 genes each. Curated from `11_gene_results.csv` (design A,
9 matched pairs, scRNA) — a **different platform and cohort, fixed before seeing the DSP per-gene
tables**. `16_genesets_LCC_evidence.csv` carries each member's pairs-higher, p, robustness share and
effect size, so nothing has to be taken on trust.

| set | strongest members |
|---|---|
| `LCC_TP53_MACROPHAGE_C1Q` | C1QB 7/9 p=0.008 (robustness 0.90), MAF 8/9, LYVE1 6/9, C1QC 6/9 |
| `LCC_COLLAGEN_CROSSLINK` | PLOD2 7/9 p=0.039, LOX 6/9 p=0.049 (0.62) |
| `LCC_PROFIBROTIC_EFFECTOR` | OSM 8/9, SPP1 7/9 (0.65), PDGFB 7/9, IL11 7/9 |
| `LCC_HYPOXIA_GLYCOLYSIS` | SLC2A1 8/9 p=0.004 (0.68), PGK1 7/9 |
| `LCC_MEGAKARYOCYTE_AXIS` | GP9 8/9 p=0.004 (0.51) |
| `LCC_COLLAGEN_STRUCTURAL_NEGCTRL` | **negative control** — COL1A1 4/9, COL1A2 2/9, COL3A1 2/9, POSTN 1/9 |

The last set is a control, not a hypothesis. We find collagen structural transcripts flat or
*lower*, and the DSP top-20 lists contain no collagen gene at all. **If an analysis reports this set
moving while the other five do not, suspect the pipeline before the biology.**

Note that several members are individually weak (LDHA, NDRG1, ADM, BNIP3 in the hypoxia set; PPBP,
VWF, NFE2 in the megakaryocyte set). They are kept because membership was fixed on biology, not
thresholded on the data — a set whose membership moves with the p-value cutoff is not pre-specified.

### Why small directional sets rather than MSigDB collections

Their slide 19 plans Reactome ECM organisation / collagen formation, NABA matrisome, Hallmark
TGF-beta. Our P4 is the empirical case against those here: across 500 alternative valid control sets
the best of 14 such collections reached p<0.05 in **5.8%** — the chance rate — 8 of 14 never reached
it in any draw, and **all three TGF-beta collections contain none of our significant genes**.
Six of our 13 significant genes (IL11, C1QB, C1QC, LYVE1, MAF, MRC1) belong to none of the 14.
A 200-gene rank score cannot register a 5-gene shift.

### Testing them at n = 3 cases

Use rotation-based self-contained tests — `limma::fry` or `roast`, which are valid at very small n —
or `camera`, which adjusts for inter-gene correlation. **Do not run standard preranked GSEA on the
H3 interaction ranked list.** Their interaction volcanoes are strongly directional (CD68+ almost
entirely negative, CD45+ and Leftover almost entirely positive) and the top genes are marrow-
irrelevant (`OR1J4` an olfactory receptor, `OLIG1`, `SLC6A2`, `ADAM2`, `RORB`, `GNAO1`, `FRMPD4`).
A global directional offset of that kind is invisible to gene-wise FDR but makes a competitive GSEA
return a wall of significant sets that are pure artefact. Diagnose the offset first — a GeoMx
limit-of-quantification filter from the negative probes, then check whether the shift tracks AOI
area / nuclei count.

---

## 3. Reproducing

```
Rscript LCC_proj/scripts/16_dsp_handoff.R \
    [--max_cells_per_type 3000] [--min_cells 25] [--min_samples 3] [--gene_presence 0.8]
```

Seed 491638. Peak memory is the binding constraint: per-sample profiles are held as nonzero values
against an integer gene index (a dense character-keyed version was OOM-killed at sample 180/220).

`stromal_modules.R` is shared with `10_stromal_denovo.R` so the per-cell stromal calls here are
identical to the ones behind the project's stromal numbers; the operating point
(`fibro_msc >= 3 & haematopoietic <= 2`, CD34-sorted FPR 8.5e-5) is read back from
`10_stromal_calibration.csv` rather than hard-coded, so a re-calibration propagates automatically.
