# Handoff: pseudobulk differential expression, AML vs healthy, per hierarchy bin

**Written** 2026-09-16. **For** a fresh Claude Code session in `aml_niche_net`.
**Status** not started. No script exists yet.

---

## 0. What this task is

Run the standard baseline analysis this project has never run: per cell type,
compare AML against healthy at the gene level, using pseudobulk and limma-voom,
with dataset as a blocking covariate.

Seven hierarchy bins, so seven independent DE analyses. The output is a gene
table per bin, plus a permutation null.

**This is not part of the graph/FGW line.** It is a separate, simpler question
that sits underneath it, and it is the analysis a reviewer will ask for first.

---

## 1. Why it exists (read this, it determines the design)

A collaborator (Eric Verbeke, Yachie lab) pointed out two things in September 2026:

1. He is not sure how much of the single-cell data actually enters the current
   analysis, because everything is reduced to 7 nodes.
2. Node vectors should probably be built from real expression, not derived scores.

He is right. The current node vector is **150 features, every one a bin-level mean
of a signature or pathway score**. Only three are genes (`expr_BCL2`, `expr_MCL1`,
`expr_BCL2L1`). The expression data enters only through gene sets someone else chose.

A repo-wide search confirms: **no `edgeR`, `limma`, `DESeq2`, `muscat`,
`FindMarkers` call exists anywhere in `scripts/`.** No AML-vs-healthy differential
expression has ever been run, pseudobulk or otherwise.

So this task closes a real gap. Whatever it finds, including nothing, is
reportable and strengthens the write-up.

---

## 2. USER-CONFIRM gates

**Do not write code before the user answers these.** Project rule: propose the
file layout and wait (`/FAST/gr10634/gaozy/CLAUDE.md`, 原则 2).

**G1 — Where do the scripts live?** Module boundaries follow pipeline stages and
are the user's call, not the AI's. Stage A (aggregation) is naturally a sibling of
`05_ccc/02_run_cellchat.R` because it reuses that exact cell set. Stage B (the DE
model) is a statistics step, but `08_scoring/` is all Python and this is R.
Options to put to the user: (a) `05_ccc/04_pseudobulk.R` + a new
`scripts/12_pseudobulk_de/`; (b) both in a new `scripts/12_pseudobulk_de/`;
(c) something else. Do **not** create `utils/`, `helpers/`, `common/`.

**G2 — Fast path or correct path for the pseudobulk?** See §4. Recommendation:
run the fast path first as a feasibility read that is declared non-reportable in
advance, then decide whether to spend the SLURM round on the correct path.

**G3 — Multiple-testing denominator.** BH within bin, or BH across all bins and
genes together? This project has already been burned once by a looser denominator
(`mp` family, panel screen). Recommendation: **pre-declare the across-bin
correction as primary** and report the within-bin version as the sensitivity arm.
Whatever is chosen must be written down before any p-value is computed.

**G4 — Dependencies.** This needs `limma` + `edgeR` (voom needs edgeR's `DGEList`
and TMM). Neither is currently called anywhere in the repo. Check whether they are
installed in `/home/b/b39170/R/x86_64-conda-linux-gnu-library/4.4/` and report
before installing anything. New dependencies require the user's approval
(原则: 不要擅自引入新依赖).

---

## 3. Step 1 is a pre-registration document, not code

House style, and it is what makes the result reportable. Precedents to match:

- `scripts/08_scoring/PREREGISTRATION_panel_screen.md`
- `scripts/08_scoring/PREREGISTRATION_paired_gate.md`

Write `PREREGISTRATION_pseudobulk_de.md` and commit it **before running anything
that produces a p-value.** It must fix, in advance:

- the sample set (§5) and the Discovery → Validation order
- the bin inclusion rule (§6)
- the model formula and what counts as a hit (§7)
- the multiple-testing denominator (G3)
- the rule for what makes a hit non-reportable (depth sensitivity arm, §7)

The project's existing disclosure problem was that the `mp` feature family was
added after a result existed. Do not repeat it.

---

## 4. Stage A — pseudobulk per (sample, hierarchy_bin)

### Path 1 (fast, same day, NOT reportable)

A full-gene pseudobulk **already exists on disk** from the LCC side project:

```
/LARGE1/gr10634/gaozy/aml_niche_net/LCC_proj/pseudobulk/<dataset>/<sample>__pseudobulk.rds
```

- 220 files, 408 MB. Built by `LCC_proj/scripts/03_percell_pass.R:163-175`.
- Structure: `list(dataset, sample, genes, n_cells, n_cells_bin, all, malignant,
  nonmalignant, by_bin)`. `by_bin` is keyed by `hierarchy_bin`.
- **Raw summed counts, full gene set** — correct input for voom, which needs
  library sizes for its precision weights.
- Path constant: `LCC_PB_DIR` in `LCC_proj/scripts/config_lcc.R:34`.

**Three caveats that make it non-reportable:**

1. Bins come from the raw projection (`LCC_BMM_DIR/<ds>/<sample>__bmm_percell.csv`),
   not the reconciled labels the CCC/FGW stack uses
   (`ANNO_RECONCILED_DIR/<ds>/<sample>__anno_percell.csv`).
   `scripts/05_ccc/03_node_features.R` records the disagreement: on GSE116256
   alone, **6.8% of cells carry a different `hierarchy_bin` after reconciliation.**
2. No `in_ccc_graph` / `high_error` filter, so stromal and high-mapping-error
   cells are included where the production cell set drops them.
3. `NA` bins are folded into an `"unassigned"` stratum that does not exist downstream.

Use it to answer one question only: **is there any DE signal at all, and is it
worth the compute to do this properly?** Declare that in the pre-registration.

### Path 2 (correct, one SLURM round, reportable)

Rebuild with reconciled bins and the production filters.

Cheapest insertion point is the load/join/filter block already written in
`scripts/05_ccc/02_run_cellchat.R:93-127` — it reads the QC `.rds`, joins bins by
exact barcode, applies `in_ccc_graph & !high_error & hierarchy_bin %in% CCC_NODES`,
and normalizes. **Note that script joins raw bins too**; for this task the join
must come from `ANNO_RECONCILED_DIR` instead.

- Aggregate with `Matrix::rowSums` over bin-split column indices on the **counts**
  layer (match the LCC idiom), not on normalized values.
- 138 CCC-eligible samples (`ccc_eligible == TRUE` in
  `results/tables/05_ccc/ccc_sample_manifest.csv`).
- Reusable SLURM array: `scripts/05_ccc/02_run_cellchat.sbatch`
  (`--array=1-${N}%8`, `p=1:t=2:c=2:m=36G`, `-t 02:00:00`, per-task TMPDIR).
  A pure aggregation pass is much cheaper than the CellChat run it was sized for.

---

## 5. Sample set — this is the part that decides everything

Ground truth from `results/tables/07_fgw/fgw_input_index.csv` (138 rows, 10 datasets,
115 AML / 23 healthy):

| dataset | AML | healthy | arms |
|---|---|---|---|
| Chen2023 | 5 | 4 | **both** |
| GSE116256 | 13 | 3 | **both** |
| GSE185381 | 37 | 10 | **both** |
| E-MTAB-11536 | 0 | 6 | healthy only |
| GSE201966 | 5 | 0 | AML only |
| GSE207356 | 3 | 0 | AML only |
| GSE227903 | 26 | 0 | AML only |
| GSE239721 | 10 | 0 | AML only |
| GSE289435 | 11 | 0 | AML only |
| Petti2019 | 5 | 0 | AML only |

**Only 3 of 10 datasets contain both arms. For the other 66 samples, dataset and
disease status are perfectly collinear** — no covariate, and no batch-correction
method, can separate perfectly collinear factors.

**Therefore: the analysis runs on the 72 both-arm samples only.** Including the
other 66 does not add information about the contrast; it adds batch.

### Use the Discovery / Validation split

It exists (`results/tables/01_preprocess/02_study_split.csv`, dataset-level 70/30
from 2026-08-04) and **no analysis in this project has ever respected it**
(`FINDINGS_topology_null.md` §16 limitation 5). This task is a clean chance to.

- **Discovery, both-arm:** Chen2023 (9) + GSE185381 (47) = **56 samples.**
  GSE185381 supplies 47/56 = 84% of it. Note this concentration in the write-up.
- **Validation, both-arm:** GSE116256 = **16 samples.**

Unlike the panel screen — where only 23 healthy samples exist cohort-wide so the
controls had to be shared between arms — here the controls are **disjoint**
(Discovery 14, Validation 3). That is a genuine improvement over the previous
design and is worth stating.

**Run Discovery first. Freeze the hit list. Only then touch Validation.**

---

## 6. Bin inclusion rule (fix before running)

Per-bin cell counts are very uneven. Medians per sample:
Mono_DC 476, LMPP_GMP 259, T_NK 245, HSC_MPP 90, B_Plasma 77, Erythroid 34,
**Megakaryocyte 2** (103 of 138 samples fall below 10 cells).

Existing gates in `scripts/config/config_ccc.R`: `CCC_MIN_CELLS_PER_NODE = 10`,
`CCC_MIN_CELLS_PER_OCCUPIED_BIN = 30`.

Proposed rule, to be fixed in the pre-registration: a (sample, bin) enters only if
it has **≥ 30 cells**, and a bin is analysed only if **≥ 5 samples per arm** clear
that in the arm being analysed. Megakaryocyte will almost certainly drop out.
Report which bins dropped and why, rather than silently analysing a bin built from
2-cell pseudobulks.

---

## 7. The model

Per bin, independently:

```
DGEList(counts) -> filterByExpr(design) -> calcNormFactors(method="TMM")
  -> voom(design) -> lmFit -> eBayes
design <- model.matrix(~ dataset + arm)      # arm: healthy = reference level
```

Test the `arm` coefficient. `dataset` is a blocking factor, not a random effect —
consistent with every other platform control in this repo (see `DECISIONS_pending.md`
O3, which records that the blueprint asked for random effects and the repo uses
fixed effects everywhere).

### The depth sensitivity arm is mandatory

AML libraries in this cohort are systematically deeper than healthy ones:
**median 4874 vs 3006 counts per cell, AUC 0.695, p = 0.0033**
(`DECISIONS_pending.md` O6, `PREREGISTRATION_panel_screen.md` amendment 2026-08-29).
On the Discovery arm, depth still tracks the label at r = +0.284 (p = 0.012)
**after** dataset fixed effects.

Pseudobulk partly handles this by construction — TMM plus voom's precision weights
are exactly how bulk RNA-seq handles unequal library sizes, which is a real
advantage over the signature-score approach where depth entered through the scoring
function itself. But it does not handle it fully.

**Pre-declare:** a hit that loses significance when median UMI per cell is added
as a covariate is recorded as depth-dependent and is **not reportable.** Same rule
as the panel screen. Decide this before seeing numbers.

---

## 8. Sanity checks — required, and they go in the code

Correctness in this project is established by checking outputs, not by reading
code (`/FAST/gr10634/gaozy/CLAUDE.md`, 原則 1). Four kinds, all four needed:

1. **Order of magnitude.** A bin's pseudobulk library size should be roughly
   (cells in that bin) × (median UMI per cell). If it is off by 10x, the
   aggregation is wrong.
2. **Known answer.** Lineage markers must be top-ranked in their own bin:
   `CD3E`/`CD3D` in T_NK, `LYZ`/`S100A8` in Mono_DC, `HBB`/`HBA1` in Erythroid,
   `MS4A1`/`CD79A` in B_Plasma, `PF4`/`PPBP` in Megakaryocyte. **If this fails,
   the bin labels are wrong and nothing downstream means anything.** Run this check
   before any DE.
3. **Edge case.** A bin with 3 cells in one sample: does it fail loudly, or does it
   silently produce a pseudobulk of near-zero counts that voom then weights
   confidently? Test it explicitly.
4. **Permutation.** Shuffle the AML label **within dataset**, rerun the whole DE,
   and confirm the hit count collapses to the BH expectation. This is the check that
   distinguishes signal from batch, and this project's history says it is the one
   that matters. Use `SEED` from `scripts/config/config_paths.R`.

---

## 9. Failure modes this project has actually hit

Read `FINDINGS_topology_null.md` §14 ("One bug found and fixed during this step")
before starting. The recurring failure here is **a new analysis arm silently
scoring the old model and exiting 0**. Concretely, watch for:

- A new code branch gated on a string prefix that the new arm does not match, so it
  falls back to the default path and reproduces the old numbers exactly.
  **If a new arm reproduces an old result to the last digit, that is a tell, not a
  result.**
- A rank or fraction written into an integer column and truncated.
- Scripts sharing one random stream, so rerunning with a different argument changes
  every p-value.
- Significance read from the uncorrected p while the corrected value is computed
  and discarded.

Every one of these produced correctly shaped output and exit code 0.

---

## 10. Operational

- **Git / the 23:30 cron.** `scripts/99_admin/daily_commit.sh` runs `git add -A`
  on the **current branch** at 23:30 and pushes. If a working branch is checked out
  at that moment, the whole day's output gets swept into it. Either stay on `main`,
  or switch back to `main` before 23:30. This is why the daily four-step loop commits
  directly to `main`.
- **Never `Read` or `cat` anything under `/LARGE1/`.** Loading a `.rds` inside R is
  fine; reading it into the conversation is not. See `aml_niche_net/CLAUDE.md` for
  the full off-limits list, including the 350 plain-text CSVs under
  `02_seurat_objects/03_bmm_projected/`.
- **Paths** come from `scripts/config/config_paths.R` (+ the `.sh` mirror, edit both
  together). Do not hardcode.
- **Script headers** are the only documentation in this repo. Every new script gets
  the house block: purpose / INPUT / OUTPUT / how invoked / WHY.
- **Do not** propose restructuring, renaming, or migrating anything. One script at
  a time, verified, committed.

---

## 11. What "done" looks like

1. `PREREGISTRATION_pseudobulk_de.md` committed, before any p-value exists.
2. Sanity check 2 (lineage markers) passing, printed, and recorded.
3. A per-bin gene table for the Discovery arm, with the pre-declared correction.
4. The within-dataset permutation null, with its hit count next to the observed one.
5. The depth sensitivity arm, and an explicit list of which hits it removes.
6. Only then: the Validation arm on the frozen hit list.

A clean null here is a publishable result and should be reported as one. It is
also the most likely outcome — every other route to the AML-vs-healthy contrast in
this cohort has returned null. Say so in advance, in the pre-registration, so that
a null does not get quietly reframed later.
