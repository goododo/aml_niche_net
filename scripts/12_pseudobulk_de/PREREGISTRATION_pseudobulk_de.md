# Pre-registration: pseudobulk AML-vs-healthy differential expression, per hierarchy bin

Written and committed **2026-09-16, before any pseudobulk matrix or p-value
exists**. That is the only thing that makes it a pre-registration rather than a
description. Nothing below was chosen after seeing a result.

Supersedes the design in `*_project_context/PSEUDOBULK_DE_HANDOFF_v1.md`. Where
this document and the handoff disagree, this document is what was run; §12
records every disagreement and why.

---

## 1. The question

Per hierarchy bin, at the gene level: does AML differ from healthy bone marrow?

The project's node vector is 150 features, of which ~141 are bin-level means of
signature, pathway, cNMF-program, pseudotime or CNV scores, and 9 are bin-level
mean expression of three genes (`BCL2`, `MCL1`, `BCL2L1`) across three strata
(named `mt_expr_*` in the node table, `expr_*` in `config_ccc.R:65`). Expression
therefore enters the analysis only through gene sets chosen upstream. No
cohort-wide, full-transcriptome, model-based differential expression has ever
been run in this repo. This closes that gap.

It is **not** part of the graph/FGW line. It is the simpler question underneath
it, and it is the analysis a reviewer asks for first.

## 2. Stated in advance: a null is the expected outcome and is the result

Every route to the AML-vs-healthy contrast through **graph topology** has
returned null: 11 cost x mass configurations at alpha=1 (p 0.305-0.966), the
transport-free Frobenius baseline (p=0.679), and the per-edge regression (0 of 49
edges within dataset).

Routes through **node features** are *not* null. They are significant and
depth-confounded, which is why `DECISIONS_pending.md` O6 and the covariate gate
in `PREREGISTRATION_panel_screen.md` exist. Do not read the topology nulls as
evidence that nothing is there; read them as evidence that the graph summary
discards it.

A clean null here is reportable and will be reported as one. It will not be
reframed later as "underpowered" or "needs a different model". §8.3's positive
control exists precisely so that "we found nothing" and "we could not have found
anything" are separable, and §11 states which one applies.

---

## 3. Sample set, fixed now

Ground truth: `results/tables/07_fgw/fgw_input_index.csv` (138 rows, 10 datasets,
115 AML / 23 healthy; the arm column is the boolean `healthy`, TRUE = control).

### 3.1 Only both-arm datasets enter

| dataset | AML | healthy | arms | used |
| --- | --- | --- | --- | --- |
| Chen2023 | 5 | 4 | both | **yes** |
| GSE116256 | 13 | 3 | both | **yes** |
| GSE185381 | 37 | 10 | both | **yes** |
| E-MTAB-11536 | 0 | 6 | healthy only | no |
| GSE201966 | 5 | 0 | AML only | no |
| GSE207356 | 3 | 0 | AML only | no |
| GSE227903 | 26 | 0 | AML only | no |
| GSE239721 | 10 | 0 | AML only | no |
| GSE289435 | 11 | 0 | AML only | no |
| Petti2019 | 5 | 0 | AML only | no |

For the excluded 66 samples, dataset and disease status are **perfectly
collinear**. No covariate and no batch-correction method separates perfectly
collinear factors. Including them adds batch, not information about the contrast.
72 both-arm samples remain.

### 3.2 Discovery / Validation

`results/tables/01_preprocess/02_study_split.csv` is a 13-dataset, **five-level**
assignment (`split` in {Discovery, Validation, Exploratory, Healthy, Reference}),
not a two-level 70/30 split. Rule-based assignment
(`scripts/01_preprocess/02_study_split.R:40-52`) fixes Healthy, Reference,
Exploratory and every longitudinal dataset first; only the residual 5-dataset
pool goes to the stratified 70/30 draw (`:57-78`, `disc_frac <- 0.70`, with the
non-empty-Validation guard at `:80-87`). **GSE116256 is Validation because it is
longitudinal (forced), not by the draw.** Stated because "dataset-level 70/30"
overstates what the file is.

Restricted to the three both-arm datasets, the split gives:

- **Discovery**: Chen2023 (9) + GSE185381 (47) = **56 samples**, 42 AML / 14 healthy.
- **Validation**: GSE116256 = **16 samples**, 13 AML / 3 healthy, reduced to 11 by §3.4.

No analysis in this project has ever respected this split
(`FINDINGS_topology_null.md` §16 limitation 5). This one does.

**Run Discovery first. Freeze the hit list. Only then open Validation.**

### 3.3 Discovery is clean at the sample level; §5 is where it stops being clean

All 42 Discovery AML samples are `timepoint == Diagnosis`; all 14 controls are
`timepoint == Healthy`. **Every Discovery sample is a distinct `Patient_ID`** —
zero repeats.

A scan of every column of `00_curated_manifest.csv` for a variable that is
constant within each arm but differs between them, inside each both-arm dataset,
returns only `disease`, `control_type`, `timepoint`, `timepoint_class`,
`cfg_timepoint` and `days_from_diagnosis` — all restatements of the arm itself.
**No hidden technical confounder (chemistry, platform, collection site, year)
separates the arms within a dataset.** Sorting, cell prep and sample state are
constant within each Discovery dataset — in particular Chen2023 applies its
FACS-recombination protocol to both arms alike. That negative result is recorded
because it is the one thing in this design that came out clean, and §5 and §7
are about everything that did not.

### 3.4 Validation is reduced, deliberately, and is not an independent replication

Two defects in GSE116256 were found before any run and are handled here rather
than discovered at write-up.

**(a) `BM5-34p` is dropped.** Per
`results/tables/01_preprocess/00_curated_manifest.csv`, it is
`sorting_r=CD34pos`, `cell_prep_r=sorted`, `sample_state=fresh`. All 13 AML
samples and the other two controls (`BM3`, `BM4`) are `viability_only`,
`mononuclear_ficoll`, `cryopreserved`. A CD34 enrichment plus a fresh/cryo
difference dominates any AML-vs-healthy logFC. Its distortion is already visible
in `ccc_node_features.csv`: `BM5-34p` has T_NK `n_cells=0` and Mono_DC 43,
against `BM4`'s 812 and 927. **The Validation control arm is n = 2.**

**(b) Validation is restricted to `timepoint == Diagnosis`.** Its 16 samples are
12 patients. `AML328`, `AML329`, `AML556` and `AML707B` each contribute a
Diagnosis sample **and** a post-treatment sample (D29 On_treatment, D20 and D31
Post_induction, D113 Post_consolidation), and all four post-treatment marrows
would otherwise be pooled into the "AML" arm. A regenerating or aplastic marrow
is a different biological state, not a second observation of the same one.

Validation is therefore **9 AML (9 distinct patients, all Diagnosis) vs 2
healthy** = 11 samples, each in its own `library_id` (verified: 11 distinct
values, max block size 1).

**(c) Validation is a directional-concordance check, fully specified here.**

- Model: **`~ arm`** on GSE116256 alone. No `dataset` term — Validation is one
  dataset and `~ dataset + arm` is rank-deficient there. No
  `duplicateCorrelation` — every Validation block has size 1.
- Gene universe: the frozen Discovery hit list, intersected with the genes that
  clear §6.3's filter in Validation. Genes lost to that intersection are counted
  and reported, not silently dropped.
- Statistic: per bin, the number of hit genes whose Validation `logFC` has the
  same sign as Discovery's, against a one-sided binomial test at p = 0.5,
  **alpha = 0.05, no multiplicity correction across bins** (three bins at most,
  fixed in advance, reported together with all three p-values).
- **If a bin's frozen hit list is empty, no test is run for that bin** and the
  result is recorded as "not evaluable — no Discovery hits". It is not pooled
  with other bins to manufacture an n. If **every** bin's list is empty, §9
  step 6 does not execute and the study's result is the Discovery null.
- Sign agreement is the claim; a Validation q-value is not produced. Decided
  before the arm is opened so it cannot be upgraded to "replication" if the signs
  agree, nor downgraded to "exploratory" if they do not.

Controls are **disjoint** between arms (Discovery 14, Validation 2; no
`Patient_ID` appears in two datasets). That is a genuine improvement over
`PREREGISTRATION_panel_screen.md`, where only 23 healthy samples exist
cohort-wide and the controls had to be shared. It is also why Validation is this
small: the 6 E-MTAB-11536 controls sit in a healthy-only dataset and cannot enter
a within-dataset contrast.

---

## 4. Cell set and aggregation (Stage A)

### 4.1 The reconciled labels are used, and the LCC fast path is not

The 220-file, 408 MB pseudobulk at `LCC_PB_DIR`
(`/LARGE1/.../LCC_proj/pseudobulk/`, built by
`LCC_proj/scripts/03_percell_pass.R:163-176`) is **not used**, not even as a
feasibility read. Reasons, all verified before the decision:

1. Its bins come from the raw projection (`LCC_BMM_DIR/<ds>/<sample>__bmm_percell.csv`),
   not `ANNO_RECONCILED_DIR`. `scripts/05_ccc/03_node_features.R` records 6.8%
   of GSE116256 cells changing `hierarchy_bin` after reconciliation.
2. It applies no `in_ccc_graph` / `high_error` filter, so it includes cells the
   production cell set drops.
3. `03_percell_pass.R:171` writes a literal `0L` into the per-bin `n_cells` slot,
   so any cell-count gate on it must read `n_cells_bin` instead — a trap, not a
   feature.
4. `NA`/empty bins are recoded to `"unassigned"` at line 91 **before** `by_bin`
   is keyed, so `by_bin` carries up to **nine** keys (the 8 raw hierarchy bins
   plus `"unassigned"`), not the seven `CCC_NODES`.
5. It spans 220 samples over 13 datasets against this analysis's 72 over 3.

A null from that object would be uninterpretable — it could only ever justify
spending the compute, never justify not spending it — so it is skipped.

### 4.2 The cell set

Reuses the load/join/filter block of `scripts/05_ccc/02_run_cellchat.R`, with two
changes: (i) bins join from **`ANNO_RECONCILED_DIR`** (`config_hierarchy.R:51`),
not the raw projection directory; (ii) the `CCC_EXCLUDE_FINE` exclusion of
`scripts/05_ccc/03_node_features.R:125-127` is applied, so the cell set is the
one behind `ccc_node_features.csv` — which is what §5 counts — and **not** the
one CellChat scored. `ANNO_RECONCILED_DIR` is **not** in `config_paths.R`;
`config_hierarchy.R` must be sourced too.

- Objects: `CCC_QC_OBJ_DIR` (`config_ccc.R:14`).
- Cell filter: `in_ccc_graph & !high_error & hierarchy_bin %in% CCC_NODES &
  !(bmm_broad %in% CCC_EXCLUDE_FINE)`, with `CCC_EXCLUDE_FINE = "Early Lymphoid"`
  (`config_ccc.R:130`). It filters on `bmm_broad`, the projection's own label, so
  `bmm_broad` must be added to the `fread(select=)` list at
  `02_run_cellchat.R:96`, which does not read that column.
- Why the fourth condition is declared rather than inherited: CellChat kept those
  cells (`02_run_cellchat.R:104` applies only the first three), the §5 table does
  not, and §5 is binding. Dropping the condition changes only `T_NK` — Discovery
  AML 39/24 -> 41/25, Validation AML 4/4 -> 5/5, the crossing units being
  `GSE185381/AML0160` 22->31, `GSE185381/AML2910` 11->139,
  `GSE116256/AML916-D0` 5->49. `Early Lymphoid` is 0.95% of the healthy
  non-malignant T_NK pool against 10.64% of the AML one (`config_ccc.R:119-123`),
  i.e. exactly the arm-asymmetric composition §7.2 is about, so it is excluded
  rather than left to the model.
- Aggregation: `Matrix::rowSums` over bin-split column indices on the **counts**
  layer. Never on normalized or log values.

### 4.3 Per-(sample, bin) depth is computed at aggregation time

`ccc_node_features.csv` carries per-(sample, bin) **cell** counts — `n_cells`,
`n_evaluable`, `n_malignant`, `n_vg_eval`, `n_vg_mal`, `n_vg_scored`,
`n_vg_pos` — but **no UMI, gene-count or depth column**;
`03_qc_report__ALL.csv`'s `med_ncount_final` is per **sample**. Per-bin depth is
free at aggregation time and unrecoverable later without rerunning, so Stage A
also writes, per (sample, bin): `n_cells`, `median_umi`, `mean_umi`,
`total_counts`, `median_genes`. §7.1's depth covariate is the **per-bin** value.

---

## 5. Bin inclusion, fixed before any model runs

### 5.1 The 30-cell unit gate is a new use of an existing constant

A (sample, bin) enters only if it holds **>= 30 cells** after §4.2 filtering.

Stated plainly: **30 is not the repo's existing per-(sample, bin) gate.** The
repo's per-node presence gate is `CCC_MIN_CELLS_PER_NODE = 10L`
(`config_ccc.R:132`, used at `02_run_cellchat.R:117`), and
`CCC_MIN_CELLS_PER_OCCUPIED_BIN = 30L` (`config_ccc.R:146`) is a **sample-level**
occupancy threshold: a sample needs `>= CCC_MIN_OCCUPIED_BINS` (=3) bins holding
30+ cells to get a graph at all. Using 30 as a unit-inclusion gate is a decision
made here. It is stricter than CellChat's: a bin with 10-29 cells is a present
node under CellChat and is excluded here. The value is reused rather than
invented so that no third threshold enters the project.

### 5.2 Bin-level rules

A **bin** is analysed on Discovery only if, among (sample, bin) units clearing
5.1:

1. each arm draws on **>= 2 of the both-arm Discovery datasets**;
2. the healthy arm spans **>= 4 distinct `library_id`** values;
3. the **cell-mass composition gap** — |Chen2023 share of healthy cells minus
   Chen2023 share of AML cells|, in cells — is **<= 25 percentage points**.

Rules 1 and 2 count **libraries, not samples**, because §6.2 shows samples are
not the independent unit. Rule 3 counts **cells**, because rules 1 and 2 can pass
while one dataset still supplies almost all of one arm's cells.

### 5.3 The table

From `ccc_node_features.csv` joined to `fgw_input_index.csv`, at the 30-cell gate,
as `samples/libraries`:

| bin | Chen AML | Chen hty | 185381 AML | 185381 hty | **Disc AML** | **Disc hty** | Val AML | Val hty | Chen cell-share hty / AML | gap | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Mono_DC | 5/5 | 4/4 | 35/19 | 10/7 | **40/24** | **14/11** | 7/7 | 2/2 | 31.4% / 24.9% | 6.5 | **primary** |
| B_Plasma | 5/5 | 4/4 | 35/20 | 10/7 | **40/25** | **14/11** | 3/3 | 1/1 | 46.3% / 29.9% | 16.4 | **primary** |
| T_NK | 5/5 | 4/4 | 34/19 | 10/7 | **39/24** | **14/11** | 4/4 | 2/2 | 43.4% / 24.4% | 19.0 | **primary** |
| LMPP_GMP | 5/5 | 4/4 | 34/19 | 2/2 | **39/24** | **6/6** | 9/9 | 2/2 | 97.4% / 13.3% | **84.1** | secondary |
| HSC_MPP | 4/4 | 4/4 | 24/17 | 4/2 | **28/21** | **8/6** | 6/6 | 1/1 | 95.0% / 9.1% | **85.8** | secondary |
| Erythroid | 3/3 | 3/3 | 18/13 | 1/1 | **21/16** | **4/4** | 6/6 | 2/2 | 97.5% / 4.7% | **92.8** | secondary |
| Megakaryocyte | 4/4 | 3/3 | 3/3 | 0/0 | **7/7** | **3/3** | 0/0 | 0/0 | 100% / 4.6% | 95.4 | **dropped** |

### 5.4 What the three tiers mean

**Primary — Mono_DC, T_NK, B_Plasma (3 bins).** Pooled two-dataset analysis per
§6. These are the only bins where both arms draw comparable cell mass from both
datasets, so `dataset` as a blocking term does what it is there to do.

**Secondary — LMPP_GMP, HSC_MPP, Erythroid (3 bins), Chen2023 only.** Rule 3
fails catastrophically: the healthy arm is 95-97.5% Chen2023 cells while the AML
arm is 4.7-13.3%. Pooled, the `arm` coefficient in these bins would largely be
*Chen2023 healthy vs GSE185381 AML* — a between-dataset contrast that the
`dataset` term cannot absorb, because there is almost no within-GSE185381 healthy
cell mass to anchor it. Chen2023 is a five-fraction FACS-recombined protocol
(`sorting_detail`: `CD45-CD235a-CD71-CD31-`, `CD45-CD235a-CD71-CD31+`,
`CD45+CD34+`, `CD45+CD34-CD117+/CD33+`, `CD45+CD34-CD117-CD33-`), pooled into two
libraries and merged per donor — its CD34+ pool is why it dominates the
progenitor bins' healthy mass, and its uniform exclusion of CD235a+/CD71+ cells
is why its Erythroid bin cannot be read as erythroid precursor biology at all.

They are therefore **not dropped and not pooled**. They are analysed
**within Chen2023 only** (`~ arm`, no `dataset` term, no
`duplicateCorrelation` — Chen2023 `library_id` is 1:1 with sample), where the
protocol is identical across arms by §3.3 and the contrast is internally valid.
n per arm is 3-5. They get their **own BH correction, within bin, and are never
pooled with the primary bins in any denominator.** They are reported as
underpowered single-dataset contrasts and any hit from them requires the §8.3
positive control to have passed in that bin.

**Dropped — Megakaryocyte.** Rule 1: zero GSE185381 controls clear 30 cells (its
best is `Control3` at 10), so its healthy arm is Chen2023-only *and* it has 0
Validation samples in either arm. Unlike the secondary bins it has no viable
within-dataset arm either: Chen2023 gives 4 AML vs 3 healthy with 154 healthy
cells total.

### 5.5 The threshold does not drive the outcome

The observed gaps are 6.5, 16.4, 19.0, 84.1, 85.8, 92.8, 95.4. **Nothing lies
between 19.0 and 84.1.** Any rule-3 threshold in [20, 84] produces this exact
partition. The 25-point value is declared now for the record, not because the
answer is sensitive to it.

### 5.6 Where marker-vs-projection agreement is low, and why it is not a gate

`results/tables/03_hierarchy/annotation_agreement_by_bin_dataset.csv` (Chen2023 /
GSE185381 / GSE116256): Mono_DC 0.853/0.722/0.814, T_NK 0.986/0.773/0.748,
B_Plasma 0.892/0.791/0.705 — versus LMPP_GMP 0.487/0.329/0.451, HSC_MPP
0.849/0.266/0.205, Erythroid 0.403/0.169/0.054.

**This is not used to drop bins**, and the reason is recorded so the number is
not misread later: `06_reconcile_annotation.R:20-24` measured both methods against
van Galen 2019's independent typing and found the projection better (0.875 vs
0.680; where they disagree, BMM right 68.4%, marker right 9.9%). Low agreement
means the **weaker** method disagreed. It is a caution flag, not evidence the bin
is wrong. It is reported because the three primary bins are also the three with
high agreement, and that coincidence should be visible rather than claimed as
independent support.

---

## 6. The model, fixed now

### 6.1 Primary bins

```r
DGEList(counts)
  -> filterByExpr(design)                    # see 6.4
  -> calcNormFactors(method = "TMM")
  -> voom(design)
  -> duplicateCorrelation(block = library_id)
  -> voom(design, block = library_id, correlation = cor1$consensus)
  -> duplicateCorrelation(block = library_id)     # re-estimated on the reweighted fit
  -> lmFit(design, block = library_id, correlation = cor2$consensus)
  -> eBayes

design <- model.matrix(~ dataset + arm)      # arm: healthy = reference level
```

Tested coefficient: `arm`. `dataset` is a blocking fixed effect. The two
voom/`duplicateCorrelation` passes are limma's documented idiom for a blocked
voom fit, not an extra precaution.

### 6.2 Secondary bins and Validation

Both are single-dataset, singleton-library designs:

```r
design <- model.matrix(~ arm)                # no dataset term, no block
```

Secondary: Chen2023 only (§5.4). Validation: GSE116256 only (§3.4c). Neither
uses `duplicateCorrelation`; in both, `library_id` is 1:1 with sample.

### 6.3 `dataset` fixed, `library_id` random; O3 is closed here and only here

`DECISIONS_pending.md` **O3 is an open item**, not a settled decision: it records
that the blueprint asked for random effects and the repo uses fixed effects
everywhere. This analysis follows current practice for `dataset` and departs for
`library_id`. O3 is closed **for this analysis only**, in writing, before the
run. It is not claimed to be settled repo-wide.

`library_id` is not optional in the primary bins. GSE185381 — 84% of Discovery
at sample level — is HTO-multiplexed: **37 of its 47 samples share a 10x library
with other samples in this analysis** (`00_curated_manifest.csv`,
`is_multiplexed`, `demux_method=hashing_HTO`, `n_donors_in_library` 2-5), and
**all four of `Control1`, `Control2`, `Control3`, `Control4` sit in one library,
`GSE185381__lib__P02`**. The §5.3 table shows the damage: B_Plasma's 14 Discovery
controls are 11 libraries, T_NK's and Mono_DC's likewise.

`library_id`'s provenance, stated because it is not what it looks like: for
GSE185381 multiplexed samples it is `GSE185381__lib__P<nn>` and genuinely shared;
for unmultiplexed GSE185381 and all GSE116256 samples it is a
`<dataset>__<sample>` string, unique per sample; for Chen2023 it is a
**semicolon-joined pair of deposit library names** (e.g.
`filtered_feature_bc_matrix.7z_AML103_CD34;filtered_feature_bc_matrix.7z_AML103_Niche_Immune`),
also unique per sample. So every non-GSE185381-multiplexed sample is a singleton
block, which is what `duplicateCorrelation` needs.

That Chen2023 string encodes a second nesting the column cannot express and the
model therefore does **not** absorb: each donor is two separately sequenced pools
(HSPC+myeloid | stromal+lymphoid) merged per donor, so **within Chen2023, library
is nested in bin**. It is recorded as a limitation (§7.3), not fixed. It is also
a direct argument against pooling bins into one correction (§6.5).

### 6.4 The gene universe

`filterByExpr(design)` derives its detection threshold from the design's
**leverage** (`edgeR` 4.4.2: `h <- hat(design); MinSampleSize <- 1/max(h)`). That
is design-dependent, not data-dependent — it depends only on the sample metadata
fixed in §3 — so it is reproducible and is stated here in advance. Per-bin
`MinSampleSize` on the real primary designs: **Mono_DC / T_NK / B_Plasma 7.04**.
(The secondary bins under `~ arm` take `MinSampleSize` from their own designs and
report it.)

`filterByExpr(group = arm)` was considered and **rejected**: it sets
`MinSampleSize = min(arm size)`, shrunk above `large.n = 10` by `min.prop = 0.7`,
which is still different per bin and preserves the same ordering — it does not
fix what it was proposed to fix.

Because §6.5 makes BH **within bin** primary, a per-bin gene universe is correct,
not a defect: each bin's correction is over exactly the genes tested in that bin.
`n_genes_tested` and the realised `MinSampleSize` are written into every output
table. **The across-bin sensitivity arm (§6.5) uses the intersection of the three
primary bins' universes**, so its denominator is a single number, computable
before the run and reported with it.

Two gene classes are removed from every universe before filtering:

- **Sex-chromosome genes** — all chrY genes plus `XIST`. **No sex or age column
  exists anywhere in the repo** (`00_curated_manifest.csv`,
  `ccc_sample_manifest.csv`, `00_MASTER_qc_summary.csv` all lack one), so donor
  sex cannot be balanced, covaried, or even measured from metadata. With healthy
  arms of n = 4-14 a chance sex imbalance is near-certain, and `XIST`, `RPS4Y1`,
  `DDX3Y`, `UTY`, `KDM5D`, `EIF1AY`, `TXLNGY` would top any such list while
  passing every other gate in this document. Excluding them is the only
  defensible option and it is declared now. Inferred sex per sample, read off
  those same genes in the pseudobulk, is reported as a **diagnostic** cross-tab
  by arm — never as a covariate, because it is derived from the data being
  tested.
- Nothing else. In particular mitochondrial and ribosomal genes stay in; see
  §7.1 on why ribosomal fraction is not a second gate.

### 6.5 Multiple testing: BH **within bin** is primary

This **reverses** the handoff's recommendation (`HANDOFF_v1.md` G3), for reasons
fixed before any p-value exists:

1. Within Chen2023, bins come from **two physically separate sequencing
   libraries** (§6.3). Pooling p-values across bins treats two separate
   experiments as one.
2. The bins have different `n_genes_tested`, different cell depths and different
   effective n. They are separate questions.
3. The primary and secondary tiers use **different models on different cohorts**
   (§6.1 vs §6.2). Pooling across tiers would be indefensible under any
   denominator.

This matches `PREREGISTRATION_panel_screen.md`, which corrects within family for
the same reason ("pooling them would spend the correction on unrelated
hypotheses").

- **Primary**: BH within bin, `q < 0.05`, over the three primary bins separately.
- **Secondary tier**: BH within bin, separately again, never pooled with primary.
- **Sensitivity, reported alongside and never instead of**: BH across the three
  primary bins pooled, over the intersected universe of §6.4.

A gene is a **hit** if `q_within_bin < 0.05` **and** `|logFC| >= 1`. Both
conditions are declared now and neither is tuned later.

**Disagreement between primary and sensitivity is defined**: any gene that is a
hit under one and not the other. The count of such genes, and their identities,
are reported for every bin. The primary is the claim.

---

## 7. What makes a hit non-reportable

### 7.1 Depth dependence — mandatory, same rule as the panel screen

AML libraries in this cohort are systematically deeper than healthy ones:
`med_ncount_final` **4874 vs 3006.5** (AUC 0.695, p = 0.0033;
`DECISIONS_pending.md` O6 `:302-326` rounds to 3007,
`PREREGISTRATION_panel_screen.md` truncates to 3006).

The partial correlation quoted in the panel screen — depth still tracking the
label at **r = +0.284, p = 0.012** after dataset fixed effects — is reproducible,
but **on a different sample set from the one this document calls Discovery**. It
is the panel screen's Discovery arm: `split_sample %in% {Discovery, Healthy}`,
77 samples over 5 datasets (57 AML + 20 healthy, including E-MTAB-11536,
GSE239721 and Petti2019, controls shared between arms). §3.2's Discovery is 56
samples from 2 datasets. The number is carried here as **prior evidence that the
confound exists in this cohort**, not as a measurement on this analysis's
Discovery arm; it is also sample-level, not per-bin, which is why §4.3 computes
per-bin depth rather than reusing it. (O6 itself carries 4874/3007, AUC 0.695,
p=0.0033 and a rho table; the partial correlation is from the panel screen's
2026-08-29 amendment, not from O6.)

TMM plus voom's precision weights handle unequal library size by construction — a
real advantage over the signature-score route, where depth entered through the
scoring function itself. They do not handle it fully.

**Fixed in advance: a hit that fails EITHER condition of §6.5's hit definition —
`q < 0.05` or `|logFC| >= 1` — when per-bin `median_umi` is added as a covariate
is recorded as depth-dependent and is NOT REPORTABLE.** The gate tests the same
two conditions the hit definition uses; testing only `q` would leave the looser
half unpoliced, and `|logFC|` is the half depth moves. No hit is re-tested with a
different covariate set to recover significance.

`median_pct_ribo` exists at `results/tables/00_ingest/00_MASTER_qc_summary.csv`
(populated 138/138) but is an **ingest-stage, pre-QC-filter** median while
`med_ncount_final` is post-filter. They are not on the same footing, and O6's
warning about over-controlling stands. **Ribosomal fraction is not added to the
primary model and is not a second gate.**

### 7.2 Arm-asymmetric mapping error — recorded, not gated

Bin labels come from projection onto a healthy bone-marrow reference, and
`results/tables/03_hierarchy/bmm_projection_summary.csv` shows `frac_high_error`
is higher in AML in all three both-arm datasets (Chen2023 0.0394 vs 0.0344;
GSE116256 0.0646 vs 0.0557; GSE185381 0.0401 vs 0.0282). The `!high_error` filter
therefore removes a larger share of AML cells, so each bin's AML pseudobulk is a
differently-selected subset of its healthy counterpart — and by construction the
surviving AML cells are the ones that most resemble the healthy reference, which
biases `logFC` **toward zero**. A null here is therefore conservative, and a hit
is not inflated by this mechanism.

This is not something the model can fix and it is not used to discard hits. It is
reported as a standing table — per bin, per arm: `frac_high_error`, pre-filter
and post-filter cell counts — and named in the limitations. §8.2's per-arm marker
check is the detector for the severe version of it.

### 7.3 Limitations recorded now, so they cannot be presented as findings later

- **Chen2023's two-pool nesting is not in the model** (§6.3). Within Chen2023,
  library is nested in bin, so the secondary tier's three bins each sit in a
  single pool and their contrasts are not independent of it.
- **The secondary tier is a single dataset with n = 3-5 per arm.** Its results
  are reported with the §8.3 power number attached, or not reported.
- **No ambient-RNA correction exists anywhere in `scripts/`** (no SoupX,
  CellBender or decontX). Haemoglobin genes outside Erythroid are expected.
- **Pseudobulk averages over sub-states within a bin.** A composition shift
  inside a bin reads as differential expression. `ccc_node_features.csv` carries
  per-(sample, bin) cNMF program means (`mp_*`) and pseudotime (`pt_*`) that
  bear on this; the per-bin arm difference in those is reported as a
  **diagnostic** alongside the DE result, and any bin where it is large has its
  hits labelled composition-confounded. It is not a gate — the two cannot be
  separated with this design, and saying so is the honest position.
- **Validation carries a treatment axis and Discovery does not**, even after
  §3.4(b), because §3.4(b) removes the post-treatment samples but cannot remove
  the fact that GSE116256 was selected into Validation *for* being longitudinal.
  "Failed to replicate" has "different cohort design" as a competing explanation.
- **Donor sex is unknown and unknowable from this metadata** (§6.4).

---

## 8. Sanity checks, in the code, with numeric pass criteria

Correctness here is established by checking outputs, not by reading code. Each
check has a criterion that can **fail**; a check that cannot fail is not a check.

### 8.1 Conservation, not order of magnitude

The handoff proposed "within 10x of cells x median UMI". That tolerance is wider
than the variation it would need to notice, and it can only catch summing the
wrong layer. Replaced with exact equalities:

- `sum(over bins) pseudobulk_counts[sample]` **==** total counts over that
  sample's §4.2-filtered cells. Exact.
- `sum(over bins) n_cells[sample]` **==** §4.2-filtered cell count for that
  sample. Exact.

**Fails on**: dropped cells, double-counted barcodes, bad barcode joins,
normalized-layer aggregation. The handoff's version fails on none of these
except the last.

### 8.2 Lineage markers, **per arm**, with a numeric threshold

Run **before any DE**. Pooled across arms the check cannot detect the failure that
would actually manufacture DE (§7.2), so it is run per (bin, arm).

Criterion: the bin's markers must rank in the **top 50 genes by mean CPM,
computed within each arm separately, in both arms** — `CD3E`/`CD3D` in T_NK,
`LYZ`/`S100A8` in Mono_DC, `MS4A1`/`CD79A` in B_Plasma; for the secondary tier
`HBB`/`HBA1` in Erythroid and `CD34` above the bin's own median in both arms for
HSC_MPP and LMPP_GMP. At least one marker of each named pair must clear it, in
each arm.

**If this fails for a primary bin, that bin's labels are wrong, the bin is
dropped, and nothing from it is reported.** A failure is not itself a finding.

`HBB`/`HBA1` appearing in a **non**-Erythroid bin is expected from ambient RBC
mRNA (§7.3) and is explicitly **not** a failure of that bin. Declared now so it
cannot be invoked either way after the fact.

### 8.3 Permutation null and positive control, fully specified

**`N_PERM = 1000`** for every permutation reported in this document.

**The exchangeable unit.** Arm is a *sample* attribute, not a library attribute:
5 of GSE185381's 22 libraries contain both an AML and a healthy sample
(`lib__P03` 2+1, `lib__P05` 1+1, `lib__P06` 2+1, `lib__P12` 2+1, `lib__P18` 1+1),
holding 5 of its 10 controls. A pure "permute whole libraries" null is therefore
**undefined** for those five, and a pure sample-level shuffle inside `lib__P02`
is not exchangeable because all four of its members are controls. Both nulls are
run and both are reported:

- **Primary null — permute `arm` among samples, within dataset.** This matches
  the exchangeability the blocked model of §6.1 actually assumes (samples
  exchangeable *conditional on* block), which is why it is primary.
- **Conservative null — permute whole libraries, within dataset, restricted to
  arm-pure libraries**, with mixed libraries held fixed. The number of permutable
  units and the minimum attainable p are computed and reported **per bin**, not
  assumed.

Reported as the full **hit-count distribution** over the 1000 permutations, not
its mean. `observed = 0` against `mean(null) = 0` carries no information, which
is why the positive control below is mandatory and not optional.

**Positive control, pre-registered numerically.** Per bin, on the real count
matrix, a planted effect is spiked into **200 genes** drawn from the middle
expression tertile of that bin's universe (seeded with `.perm_seed("planted|"
+ bin)`), at **logFC = 1.0** applied to the AML arm, and the full §6 pipeline is
rerun.

- **Realised power** for a bin = the fraction of those 200 planted genes
  recovered under §6.5's full hit definition (`q_within_bin < 0.05` **and**
  `|logFC| >= 1`).
- **A bin with realised power < 0.50 is declared power-limited**, and its null is
  reported as *uninformative*, not as a null.
- **A bin with realised power >= 0.80 supports the strong statement**: no effect
  of |logFC| >= 1 at this prevalence is present.
- Between 0.50 and 0.80, the null is reported with the power number attached and
  no strong statement.

These thresholds are fixed now, before the effect size is known, precisely so
that a post-hoc choice of planted logFC cannot make the power read high or low at
will.

### 8.4 Seeding

Per-test streams via the existing `.perm_seed(label)` helper
(defined at `scripts/05_ccc/04_stemness_purity_sweep.R:111-116`, with the
rationale at `:95-110`), **not** a bare `set.seed(SEED)` at the top.
`SEED = 491638L` (`config_paths.R:28`). A single shared stream makes a test's
p-value depend on how many tests ran before it — that bug is written up mid-script
at `:95-100` and is not repeated here.

### 8.5 Leave-one-out on the marginal control

The handoff's "bin with 3 cells" edge case tests a state the 30-cell gate makes
unreachable. The reachable edge case is a unit sitting just above 30.

Applies to **LMPP_GMP and Erythroid in the secondary tier**, whose GSE185381
controls are `Control0005` (35), `Control1` (41) and, for Erythroid, `Control1`
alone at exactly 30. Because the secondary tier is analysed **within Chen2023
only** (§5.4), those marginal GSE185381 controls do not enter the secondary model
at all — so the §5.2 rule-2 circularity that a Discovery-pooled leave-one-out
would have created (dropping `Control1` would drop Erythroid below its own
inclusion rule) does not arise. Stated explicitly because it was a real hole in
the pooled design and its absence here is a consequence of §5.4, not an oversight.

What is still run: within the secondary tier, **leave-one-out over each healthy
sample in turn** (3-4 per bin). Criterion, fixed now: a secondary hit is reported
as **robust** only if it remains a hit (§6.5, both conditions) in **every**
leave-one-out fit; otherwise it is reported as **contingent**, with the number of
folds in which it survived. "Survive" is thereby a count, not a judgement.

**If a bin's hit list is empty, no leave-one-out is run for it** and the result
is recorded as "not evaluable — no hits".

### 8.6 The new-arm-reproduces-the-old-result tell

`FINDINGS_topology_null.md` §14 records the recurring failure: a new analysis arm
silently scoring the old model and exiting 0. The baseline arm is run first as a
gate; its `logFC` must reproduce bit-identically across reruns. **The depth arm
(§7.1), the permutation arms (§8.3) and the planted-effect arm (§8.3) must each
differ from the baseline somewhere.** If any of them reproduces the baseline to
the last digit, that is a bug, not a result: stop and find the branch that did not
take.

---

## 9. Order of operations. Nothing skips ahead.

1. `PREREGISTRATION_pseudobulk_de.md` committed. **No p-value exists yet.**
2. Stage A aggregation (`01_aggregate.R`), 72 both-arm samples, §4.2 cell set,
   with per-bin depth (§4.3).
3. Checks 8.1 and 8.2. **Both must pass or the run is void**; a primary bin
   failing 8.2 is dropped.
4. Primary-tier DE (§6.1), three bins, within-bin BH (§6.5). Sensitivity arm on
   the intersected universe. Hit lists **frozen**.
5. Secondary-tier DE (§6.2), three bins, Chen2023 only, own BH, never pooled.
   Hit lists **frozen**.
6. §8.3: primary null, conservative null, and the planted-effect positive
   control with realised power **per bin**. §7.1 depth arm. §8.5 leave-one-out on
   the secondary tier. Non-reportable hits struck out.
7. Only then: Validation sign concordance (§3.4c) on the frozen **primary** hit
   lists. The secondary tier has no Validation arm — §5.3 shows GSE116256 gives
   it 1-2 controls — and this is fixed now, not decided after seeing step 5.

Nothing that fails step 3 proceeds. Nothing struck at step 6 is re-tested with a
different covariate set. Nothing that fails Discovery is looked at again on
Validation. If a step's input list is empty, that step is recorded as not
evaluable and the run continues to the next.

---

## 10. What this cannot settle, and must not be used to claim

- **Anything about Megakaryocyte.** Dropped (§5.4); the cohort cannot support it.
- **Anything about stroma.** `Stromal` is outside `CCC_NODES` and outside this
  analysis.
- **Erythroid biology.** Its healthy arm is 97.5% Chen2023 cells, and every
  Chen2023 FACS fraction gates out CD235a+/CD71+ cells (§5.4). Whatever the
  secondary tier finds in that bin is not a statement about erythroid precursors.
- **Independent replication.** §3.4(c) — Validation is 9 vs 2 and gives a sign
  check on three primary bins, nothing more. A concordant sign is not a
  replicated effect.
- **That a null means no biological difference.** This is settled per bin by
  §8.3's realised power, not by assertion. The design's floor is the **secondary
  tier**: Erythroid at **4 healthy samples in 4 libraries**, Chen2023-only with
  3 controls once restricted (§5.4) — not the 14-control figure that applies only
  to the three primary bins. A null is a null **at the reported power**, and the
  power is reported next to it.
- **Treatment response.** Removed by §3.4(b), not studied.
- **Malignant vs non-malignant.** Not a split here. inferCNV under-calls ~9x
  (`PREREGISTRATION_panel_screen.md`); every §4.2-surviving cell in a bin enters
  its pseudobulk regardless of CNV call.
- **Anything sex-linked.** chrY and `XIST` are excluded (§6.4) and donor sex is
  not recorded anywhere in this project.

---

## 11. What "done" looks like

1. This document committed before any p-value exists.
2. §8.2 passing, printed, recorded per (bin, arm).
3. Per-bin gene tables: three primary bins, three secondary bins, each with
   `n_genes_tested`, `MinSampleSize`, and both correction arms.
4. Both permutation nulls (§8.3) as distributions, with observed counts beside
   them, and the per-bin minimum attainable p for the conservative null.
5. Realised power per bin from the planted-effect control, and each bin's null
   classified informative / power-limited by §8.3's thresholds.
6. The depth arm (§7.1) with an explicit list of which hits it removes.
7. The secondary tier's leave-one-out, each hit labelled robust or contingent.
8. Only then: Validation sign concordance on the frozen primary lists, with the
   three binomial p-values and the count of hit genes lost to the universe
   intersection.

A null with realised power >= 0.80 is a publishable result and is reported as
one. A null with realised power < 0.50 is reported as an uninformative test, not
as a null. §2 pre-commits to this distinction so it cannot be blurred afterwards.

---

## 12. Every place this departs from `PSEUDOBULK_DE_HANDOFF_v1.md`

| # | handoff | here | why |
| --- | --- | --- | --- |
| 1 | G4: check whether limma/edgeR are installed; may need approval | already installed | limma 3.62.2, edgeR 4.4.2, statmod, locfit, DESeq2 in `ENV_PREFIX` (`/FAST/gr10634/gaozy/general_env`). The handoff checked `~/R/.../4.4/`. Scripts must run through the conda env; login-node `Rscript` is R 4.5 and cannot see limma. |
| 2 | G2: run the LCC fast path first | skipped | §4.1 — wrong bins, wrong cohort, `n_cells` slot is a literal 0, up to 9 `by_bin` keys. Could only justify spending compute, never justify not. |
| 3 | G3: across-bin BH primary | within-bin BH primary, per tier | §6.5 |
| 4 | `~ dataset + arm` | `+ duplicateCorrelation(block = library_id)` in the primary tier; `~ arm` in the secondary tier and Validation | §6.1-§6.3 — `Control1/2/3/4` are one lane; single-dataset arms are rank-deficient with a `dataset` term. |
| 5 | O3 cited as settled precedent | O3 is open; closed here for this analysis only | §6.3 |
| 6 | seven bins, one analysis | **three primary, three secondary (Chen2023-only), one dropped** | §5.4 — the decisive finding. In LMPP_GMP, HSC_MPP and Erythroid the healthy arm is 95-97.5% Chen2023 cells against an AML arm that is 4.7-13.3%, so pooled they would be a between-dataset contrast wearing an `arm` label. |
| 7 | ">= 5 samples per arm" | >= 2 datasets per arm, >= 4 healthy libraries, **and cell-mass gap <= 25 pp** | §5.2 — the sample-level rule is unachievable within-dataset and counts the wrong unit; libraries and cells are the units that bind. |
| 8 | 30 cells described as an existing gate | declared as a **new use** of an existing constant | §5.1 — the repo's per-(sample,bin) gate is `CCC_MIN_CELLS_PER_NODE = 10`. |
| 9 | Validation = GSE116256, 16 samples | 9 AML vs 2 healthy, `~ arm`, sign check only, empty-list behaviour defined | §3.4 — `BM5-34p` is CD34+ sorted; 4 samples are post-treatment; 16 samples are 12 patients. |
| 10 | `filterByExpr(design)` | `filterByExpr(design)` **kept**, plus chrY/`XIST` exclusion and an intersected universe for the sensitivity arm | §6.4 — `group=arm` was considered and rejected; it does not fix the ordering it was proposed to fix, and per-bin universes are correct under within-bin BH. |
| 11 | no mention of donor sex | chrY + `XIST` excluded; inferred sex reported as a diagnostic | §6.4 — no sex column exists anywhere in the repo, so it cannot be balanced or covaried. |
| 12 | depth covariate = sample-level median UMI | per-(sample, bin) median UMI, computed in Stage A | §4.3 |
| 13 | depth gate on significance | gate on **both** `q < 0.05` and `\|logFC\| >= 1` | §7.1 — a gate looser than the hit definition does not police it. |
| 14 | sanity 1: within 10x | exact conservation equalities | §8.1 |
| 15 | sanity 2: markers "top-ranked" | per arm, top 50 by mean CPM within arm, both arms, named pairs | §8.2 |
| 16 | sanity 3: a 3-cell bin | leave-one-out over every secondary-tier control, hits labelled robust/contingent by fold count | §8.5 |
| 17 | sanity 4: "hit count collapses to the BH expectation" | two nulls (sample-level primary, arm-pure-library conservative), full distributions, `N_PERM = 1000`, plus a numerically specified planted-effect control with power thresholds | §8.3 — 0-vs-0 is not a passing check, and "permute libraries" is undefined for the 5 mixed lanes. |
| 18 | "every other route returned null" | every route through **graph topology** returned null; node-feature routes are significant and depth-confounded | §2 |
| 19 | split described as "dataset-level 70/30" | 13 datasets, 5 levels; 70/30 applied to a 5-dataset pool at `:57-78` | §3.2 |
| 20 | median healthy depth 3006 | 3006.5 | §7.1 |
| 21 | r = +0.284, p = 0.012 "on the Discovery arm" | reproducible, but on the **panel screen's** 77-sample Discovery+Healthy arm, not this document's 56-sample Discovery | §7.1 — the word was redefined in §3.2 and the number did not move with it. |
| 22 | "no DE has ever been run" | no **cohort-wide, full-transcriptome, model-based** DE has been run | A TP53-mut vs WT detection-rate comparison exists in `LCC_proj` over a 148-gene curated panel, using no DE package. |
| 23 | node vector "every one a signature mean", 3 gene features | ~141 of 150 are score means; **9** are gene features (3 genes x 3 strata), named `mt_expr_*` | The 27 `mp_*` features are this project's own cNMF programs, not externally chosen gene sets. |

---

## 13. Files

```
scripts/12_pseudobulk_de/
  PREREGISTRATION_pseudobulk_de.md   this file
  01_aggregate.R                     Stage A: (sample, bin) -> raw summed counts + per-bin depth
  02_de_limma.R                      Stage B: primary tier (blocked voom) + secondary tier (~ arm)
  03_permute.R                       both nulls + planted-effect positive control
```

**Amendment 2026-09-16, implementation only, no design content.** This section
originally listed `01_aggregate.sbatch`, a SLURM array modeled on
`05_ccc/02_run_cellchat.sbatch`. Stage A turned out to run the whole 67-sample
roster serially in ~7 minutes (6-9 s per sample; it is a `rowSums` over a counts
matrix, not the CellChat fit the array was sized for), so the array was not
written. Nothing about the cell set, the model, or any decision rule changes.
Recorded here rather than silently dropped.

Paths come from `scripts/config/config_paths.R` (+ the `.sh` mirror, edited
together), `config_ccc.R`, and `config_hierarchy.R`. Nothing is hardcoded.
