# AML Niche Network — project guide

Language rule: **conversation in Chinese, everything else in English** (code, comments, figures,
tables, docs — this file included). Per CODING_STANDARDS.md §10.

Integrated analysis pipeline over 13 public AML / healthy bone-marrow scRNA-seq datasets:
**ingest → QC → malignant-cell labelling → hematopoietic hierarchy projection (7 bins) → CellChat cell-cell communication graphs → FGW graph alignment → statistical testing.**

- Project root: `/FAST/gr10634/gaozy/aml_niche_net` (scripts + tables + figures, purgeable tier)
- Large objects: `/LARGE1/gr10634/gaozy/aml_niche_net` (Seurat objects, CNV, cNMF — **785G**)
- All paths and the global `SEED=491638` are defined in one place, `scripts/config/config_paths.R`,
  mirrored for the cluster in `config_paths.sh`. **Edit both together.**
- No formal docs. The header comment block of each script is the only spec — read it before changing anything.
  (Exceptions: `.claude/settings.local.json` is this repo's agent permission config;
  `00_project/tmp/mp_bin_check.txt` and the root-level `recon_*` files are one-off check notes, not specs.)

---

## Established premises (do not re-litigate)

### 1. Malignancy labelling runs on inferCNV as the main arm; Numbat covers only a few datasets

**Decision**: Numbat needs FASTQ/BAM to run, and **only two or three datasets in the cohort qualify**.
It therefore covers only a small subset of samples by nature and **cannot be a cohort-wide evidence
source** — this is a data-availability limit, not a reliability problem with Numbat itself.
The main labelling arm is inferCNV (the only arm that covers the whole cohort):

```
20_refnorm_identify.R → 44_infercnv_run_one.R (+45 array) → 41_infercnv_to_percell.R → 50_consensus_malignancy.R
```

`50_consensus_malignancy.R` is the **sole malignancy labeller** (multi-arm union vote + evidence tier).

**But the actual state on disk and in code does not fully match this decision — know this before taking over:**

- **125 of 130 samples are `arms=infercnv` / `tier=C_single` (single evidence).
  5 GSE227903 samples already carry Numbat evidence in their labels**:
  `1886_Dg`, `1886_R`, `3853_Dg`, `6323_Dg` are `B_multi_partial`, `6323_R` is `A_concordant`.
  Their `conflict_frac` runs 0.14–0.76.
- **The Numbat arm is not turned off in code**: `run_consensus_all.sh:50` appends `--numbat`
  whenever a numbat percell file exists; `50_consensus_malignancy.R:137` feeds it into the union vote as `allele_call`.
- **bamchain / Numbat is still running** (as of 2026-07-24 a SLURM array `bamchain` is R/PD;
  `logs/bamchain.*` spans 07-16 ~ 07-24). GSE289435 has produced 5 `__numbat_percell.csv` (07-16 ~ 07-23),
  **but they were never merged into consensus** — the dataset's 12 sample labels are all still the 07-14 `arms=infercnv`.

> **Decided (2026-07-24)**: do **not** kill the bamchain array — let it finish; a few more samples
> with allelic evidence is net positive. But Numbat is not the cohort-wide strategy.
>
> ⚠️ **Still open**: when to merge GSE289435's 5 finished Numbat results.
> Re-running `50_consensus` changes that dataset's labels and invalidates everything downstream (03–08),
> so it **must be folded into one run together with** the post-metadata-fix rerun, never triggered alone.

**Next step to improve label quality (decided, run after the lab meeting)**: run `42_exprcnv_run.R`
(CopyKAT/SCEVAN). It reads the per-sample QC `.rds` directly, **needs no FASTQ/BAM**, and is therefore
the only second independent evidence arm that can cover the whole cohort — it can lift a large batch of
samples out of `tier=C_single` into multi-evidence, directly attacking HSC_MPP's high false-positive rate.
This is the most realistic approximation of the blueprint's Phase 1 ">= 2/3 three-way consensus" under
the current data conditions.

Status of the other arms:
- `42_exprcnv_run.R` (CopyKAT/SCEVAN): never run.
- `43_author_to_percell.R`: **has run** (15 `__author_percell.csv` under `/LARGE1/.../03_cnv_snv/author/`,
  Chen2023 x1 + GSE239721 x14, 06-19/06-25), but was never passed to `50_consensus` via `--author`
  (logs show `[author] absent` for all 130 samples).

Most labels are single-method conclusions — treat them as "single evidence" in downstream interpretation.
Two negative controls exist: `96_malignancy_fpr_healthy.R` (healthy-donor false-positive rate),
`94_diagnose_label_quality.R` (sequencing-depth confounder).
`51_malignancy_plan.py` scans each sample's currently reachable evidence tier.

### 2. Niche redefined as the hematopoietic-immune microenvironment (stroma excluded)

**Two datasets carry stromal content, but they cannot form a contrast**:

| Dataset | Role (`config_qc.R` ROLE_TABLE verbatim) | Stromal cells |
|---|---|---|
| Chen2023 | `Discovery-AML+AuxStroma` | ~6,284 (19 samples) |
| GSE253355 | `Reference-scaffold+AuxStroma` | ~25,370 (12 samples) |

**Chen2023 is the only AML dataset with stroma; GSE253355 has more stroma but is a pure healthy-donor
reference** (Bandyopadhyay 2024, `Role : healthy reference + auxiliary stroma`, all baseline, no timepoints).
So a cross-cohort "AML vs healthy" stroma-hematopoietic communication contrast has only one arm and cannot
be done. The project therefore **redefines "niche" as the hematopoietic-immune microenvironment**, with
stromal cells excluded from the analysis.

This premise is encoded in code, not just a verbal agreement:

- `scripts/03_hierarchy/bmm_bin_map.tsv`: `Stromal → in_ccc_graph = FALSE`
  (comment verbatim: `non-hematopoietic; NOT a CCC-graph node (B-layer add-on, Chen2023 later)`)
- `scripts/config/config_ccc.R:20`: `CCC_NODES` is fixed to 7 hematopoietic bins —
  `HSC_MPP, LMPP_GMP, Mono_DC, Erythroid, Megakaryocyte, T_NK, B_Plasma`

The FGW 7x7 graph structure, the barycenters, and the 49-edge tests in 08_scoring all rest on these 7 nodes.
**Adding or removing a node invalidates every output of stages 06/07/08.**

`scripts/判断样本MSC含量.r` is the one-off screening script (non-pipeline) that originally drove this decision.

---

## Directory guide

### `scripts/` — all code (205 files)

~103 active files; the other 102 are archived under `以前06_hierarchy/` (85 of them old SLURM logs).

| Directory | Contents |
|---|---|
| `config/` | 10 config libraries, only `source()`d, never run standalone. `config_paths.R`/`.sh` is the single source of truth for paths |
| `00_ingest/` | 13 `ingest_<dataset>.R` + shared reader `00_common_readers.R` (**where files are actually written**) + MASTER rollup |
| `01_preprocess/` | role manifest → 70/30 study-level split → `03_per_sample_qc.R` (core, SLURM array of 13) → landing check |
| `02_malignancy/` | **24 files**, the most complex stage. Audit/plan tools: `51_malignancy_plan.py` (scans each sample's reachable tier and emits `_generated/run_consensus.sh`), `91_check_infercnv.R`, `92_check_numbat.R`, `93_check_bams.sh` |
| `03_hierarchy/` | BoneMarrowMap projection → per-bin malignant fraction → distribution-shift test → stemness scoring. Strictly linear 01→02→03, with 04 hanging off 01 |
| `04_cnmf/` | GeneNMF meta-programs (branch A). `mp_labels.tsv` is a hand-curated table |
| `05_ccc/` | CellChat communication graphs (branch B, main line). `02_run_cellchat.R` runs as a SLURM array |
| `06_distance/` | single script: LR tensor → directed 7x7 edge weights → within-sample rank distance `C = 1 - rank` |
| `07_fgw/` | FGW graph alignment (R builds inputs → Python computes barycenter/HDS/ATS → H3 relapse test) |
| `08_scoring/` | **11 files** (9 `.py` + 2 `.sbatch`), the statistical-inference layer, **the current active frontier** |
| `10_figures/` | 1 shared config (`f00_fig_config.R`) + 5 figure scripts (incl. `f03_*_test.R` trial version), **written but never run**. NB the assessment figures live in `g_assessment_figures.R` → `results/figures/11_assessment/` |
| `以前06_hierarchy/` | **archived legacy** ("以前" = previous), superseded by `03_hierarchy/` + `04_cnmf/`. Do not modify, do not place in the execution order |

Two top-level oddments: `00.0_Job_submit.sh` (SLURM/tmux command cheat-sheet, **not a driver script**),
`判断样本MSC含量.r` (one-off MSC screen, see premise 2 above).

### Scattered `recon_*` probes (one-off diagnostics, non-pipeline)

`00_ingest/recon_GSE185381_structure.R`, `05_ccc/recon_ccc_probe.R`, `05_ccc/recon_node_feasibility.R`,
`07_fgw/recon_fgw_probe.py`, `07_fgw/recon_fgw_directed_probe.py`.

**Their output is written to the project root, not `results/tables/`**: `recon_node_counts.csv` (100K),
`recon_ccc_probe.txt`, `recon_fgw_probe.txt`, `recon_fgw_directed_probe.txt`,
`recon_GSE185381_{donor_qc,donor_span,library_mix}.csv`.
Do not interpret them as stage outputs, and do not place them in the execution order.

### `results/tables/` — per-stage outputs

`00_ingest/` `01_preprocess/` `02_malignancy/` `03_hierarchy/` `04_cnmf/`
`05_ccc/` `06_distance/` `07_fgw/` `08_scoring/` `09_robustness/` (empty) `side_tp53_pvrl4/`

The pipeline head is currently `results/tables/08_scoring/`.
Post-hoc negative-control outputs: `02_malignancy/{sample_quality_vs_malignancy,malignancy_fpr_healthy,malignancy_fpr_by_bin}.csv`,
`00_ingest/timepoint_rule0_dryrun.csv`.

**Note: `03_hierarchy` per-cell outputs are NOT here** — `__bmm_percell.csv` /
`__stemness_percell.csv` live under `/LARGE1/.../02_seurat_objects/03_bmm_projected/` (see below).

### `results/figures/`

Only `01_preprocess/` (13 pngs), `04_cnmf/` (12 files = 5 pdf+png pairs + 2 pngs), and `11_assessment/`
(the 5 lab-meeting figures, g01–g05) have content.
`00_ingest/` and `02_malignancy/` are **empty directories** (placeholders, not evidence a figure script ran).

### Logs

`logs/` (344 files, 11M) + `00_project/logs/` (26 files, 3.7M) = 370 files / 14M.
A further 85 old logs live under `scripts/以前06_hierarchy/logs/` (archived, ignore).

---

## Large-data directories NOT to read

**Never open any path below with Read, and never `cat` files inside them.**
To learn what is there, use `ls` / `du` / `find -printf` for names, sizes, mtimes, or read the script header schema.

### Absolutely off-limits (a single file can blow up the context window)

| Path | Size | Note |
|---|---|---|
| `/LARGE1/gr10634/gaozy/aml_niche_net/` | **785G** | whole tree. Enumerated below |
| `└─ 05_cnv_snv/` | 374G | inferCNV run dirs + burden |
| `└─ 03_cnv_snv/` | 365G | Numbat / cellsnp / STARsolo workspace (incl. BAMs); also `author/` with 15 percell csv |
| `└─ 00_raw/` | 23G | raw public downloads (GEO/Zenodo/ArrayExpress) |
| `└─ 01_processed_counts/rds/` | 8.8G | 13 merged Seurat objects |
| `└─ 06_cnmf/` | 6.7G | legacy cNMF input/runs |
| `└─ 02_seurat_objects/` | 4.7G | 220 per-sample QC objects + **350 plain-text CSVs** (see warning below) |
| `└─ 05_ccc_graphs/` | 2.6G | 148 CellChat objects |
| `└─ reference/` | 200M | SingleR/inferCNV external refs. **`gencode_GRCh38_gene_order.txt` is 1.2M plain text, the extension list won't catch it — also do not read** |
| `└─ 04_cnmf/` | 309M | NMF results `nmf_res.rds` |
| `└─ objects/`, `└─ 99_archive/` | empty | listed for completeness |
| `/FAST/.../SRR32323369/` | **30G** | `MLL_16703.bam` — leftover download, not a pipeline output |
| `/FAST/.../SRR32323368/` | 2.4G | `MLL_17746.bam.tmp` — **incomplete download** (has a .lock) |

### Skip by extension (in any directory)

`.rds` `.RDS` `.h5ad` `.h5` `.mtx` `.bam` `.bam.tmp` `.npz`

> ⚠️ **The extension list only blocks binaries. The following are plain text but must also not be read:**
> - **350 `.csv`** under `/LARGE1/.../02_seurat_objects/03_bmm_projected/`
>   (`__bmm_percell.csv` / `__stemness_percell.csv`, largest single file **2.3M**)
> - `/LARGE1/.../reference/gencode_GRCh38_gene_order.txt` (1.2M, ~60k lines)

### FAST-side, handle with care (not binary, but large)

| Path | Size | How to handle |
|---|---|---|
| `results/tables/00_ingest/*_qc_percell.csv.gz` | 90M total | largest single file **49M** (GSE185381). Read `_qc_summary.csv` only, do not decompress percell |
| `results/tables/02_malignancy/` 9 subdirs | 49M / 265 files | largest single 1.4M. Do not iterate per-sample; read the 4 root summary tables |
| `results/tables/05_ccc/ccc_edge_distance.csv` | **728K** | duplicate copy, known-issue 5 |
| `results/tables/05_ccc/ccc_node_features.csv` | 104K | sample with `head`/`awk` when needed |
| `results/tables/05_ccc/tensors/` | 148 LR tensors | read a single sample only, do not glob all |
| `results/tables/06_distance/edge_distance.csv` | **728K** | sample only |
| `results/tables/07_fgw/fgw_edges_long.csv` | **447K** | sample only |
| `results/tables/07_fgw/fgw_nodes_long.csv` | 135K | sample only |
| `logs/` | 344 files / 11M | group-count by prefix; read at most the first 30 lines of a log |
| `00_project/logs/` | 26 files / 3.7M | **largest single 521K** (`01_qc_3622821_*.err`). Never read whole; use `tail -30` / `grep` |

### Safe read entry points (all small tables, read directly)

```
results/tables/{01_preprocess,03_hierarchy,04_cnmf,08_scoring}/*.csv
results/tables/00_ingest/*_qc_summary.csv          # excludes *_qc_percell.csv.gz
results/tables/02_malignancy/*.csv                 # 4 root tables only, all <11K, non-recursive
results/tables/06_distance/edge_qc.csv             # 9K
results/tables/07_fgw/{fgw_input_index,patient_scores,paired_rls_scs,rls_grouped}.csv
results/tables/05_ccc/ccc_edge_qc.csv              # 9K
```

---

## Execution order

**There are no SLURM `--dependency`/`afterok` anywhere. Stage-to-stage sequencing is entirely manual.**
Only two drivers actually chain multiple scripts: `run_consensus_all.sh` (41→50) and
`10_prealign_one_sample.sh` (download→STARsolo→numbat pileup_and_phase→30_numbat_run.R,
submitted by `11_submit_array.slurm:74`). Additionally `51_malignancy_plan.py` emits (but does not run) a 41/42/43→50 plan.

```
[1] 00_ingest/ingest_<13 datasets>.R → ingest_MASTER_summary.R ; (side audit) 95_test_timepoint_rule0.R
[2] 01_preprocess: 01 → 02 → 04 ; (preflight 91/90) → 03_per_sample_qc.R [sbatch array] → 05
        ↓ per-sample QC .rds  ← the project's central hub, read by 7 downstream stages
[3] 02_malignancy: 01_gene_order ; 20_refnorm → 90_preflight → 44+45[array] → 91_check ↺
                   → run_consensus_all.sh (41 → 50) → 60_rollup
                   (side branch) 11+10[array] → 30_numbat  ← still running, see premise 1
                   (post-hoc negative controls, rerunnable) 94_diagnose_label_quality / 96_malignancy_fpr_healthy
[4] 03_hierarchy:  01_bmm_project → 02_per_bin → 03_shift ; 01 → 04_stemness
[5A] 04_cnmf:  01(normal) → 02(normal) → 01(malignant) → 02(malignant) → 03 → 04
[5B] 05_ccc:   01_define → 02_run_cellchat[array] → 03_node_features
[6] 06_distance/01  →  [7] 07_fgw: 01(R) → 02(py) → 03 → 04
[8] 08_scoring: 01 → 02 → 03 → 04 → 05 → 06 → 07/08/09   ← current frontier
```

### Easy traps

- **The numbering in `02_malignancy` is not the execution order**: `20` must run before `40/44`;
  `30_numbat_run.R` is never called directly (invoked inside `10_prealign_one_sample.sh:272`);
  `02/03/04/10/11/30` are a parallel FASTQ side branch, not a step on the main line.
- **`04_cnmf/02` has a self-dependency**: the normal run's output is the malignant run's input, so order matters.
- **`08_scoring` numbering is narrative order, not a data dependency**: besides `01→05`,
  `02/03/04/06/07` read only the three `07_fgw` CSVs; **`08/09` additionally read `06_distance/edge_distance.csv`**
  (using only nodes/index, not edges). These scripts can be rerun independently and in parallel.
- **`07_fgw/02_fgw_align.py` does not source `config_fgw.R`** — it re-declares the same constants inline
  (comment says "keep in sync"). Changing an FGW parameter means editing two places.
- `config_paths.R`'s `PROJ_OBJ_DIR` points at `04_bmm_projected` and is a **dead constant**;
  the real path is `03_bmm_projected` in `config_hierarchy.R`.

---

## Known issues (as of 2026-07-24, not yet addressed)

1. **`results/tables/02_malignancy/ref_norm_summary.csv` does not exist on disk**, yet `40/44/91`
   `stopifnot` it, and `90_preflight` checks then `quit(status=1)` with "run 20_refnorm_identify.R first".
   The routing table itself can be rebuilt with `20_refnorm_identify.R` (inputs complete: 220 QC rds +
   106 `ref_norm_cells.txt`), but the QC objects were refreshed on 07-15, so **the rebuilt routing table
   no longer equals the one the 07-14 labels were based on — what is lost is reproducibility of the original run.**
2. **Timestamp vs stage-numbering contradiction**: `02_malignancy` per-sample consensus outputs are 07-14
   (the rollup `ALL_consensus_summary.csv` is 07-17, the negative-control CSVs are 07-22), while
   `00_ingest`/`01_preprocess` are 07-15 — ingest and QC were rerun **after** malignancy labelling, so the
   on-disk malignant labels correspond to an older batch of QC cells. `90_preflight_infercnv.R`'s comment confirms this is known.
3. **GSE289435 has 5 new Numbat results that were never merged** (07-16 ~ 07-23), and `50_consensus`
   was never rerun: the dataset's 12 `__consensus_percell.csv` are all still 07-14 18:12.
   So **downstream (03/04/05/06/07/08) is consistent with the on-disk labels**; what is inconsistent is the
   "produced but not merged" Numbat evidence. See the open item under premise 1.
4. **Incomplete cohort**: all 13 datasets have ingest + QC (220 QC objects), but only 9 datasets / 130 samples
   have malignancy consensus (missing E-MTAB-11536, GSE147989, GSE185991, GSE253355).
   Coverage comparison: 220 `__bmm_percell` / 130 `__stemness_percell` / 148 CCC tensors.
   **stemness = 130 is not a gap**: `04_stemness_score.R` takes the intersection QC ∩ consensus ∩ projection,
   which is exactly the 130 samples with consensus. The 148 CCC tensors are not gated on consensus — a separate filter.
5. **`results/tables/05_ccc/ccc_edge_distance.csv` and `ccc_edge_qc.csv` are duplicate copies**,
   byte-identical (matching md5) to `06_distance/edge_distance.csv` / `edge_qc.csv`, with an mtime 19 minutes
   earlier; no script anywhere references these two filenames. Safe to delete.
6. **`scripts/10_figures/` (f01–f04) never ran**: they do not source the project config and hard-code legacy
   paths (`results/tables/01_qc/`, `03_malignancy/`), so they currently fail at `stopifnot`.
   Also `03_qc_report__ALL.csv` does not exist (only the 13 per-dataset versions).
   (The English assessment figures `g_assessment_figures.R` → `results/figures/11_assessment/` do work.)
7. **`results/tables/09_robustness/` is an empty directory with no corresponding script directory.**
8. **`results/tables/side_tp53_pvrl4/` are orphan outputs** — a tree-wide grep finds no script that produces them.

---

## Working conventions

- Before changing a config constant, check whether it is redeclared elsewhere (`config_paths.R`/`.sh` are a pair;
  `config_fgw.R` and `02_fgw_align.py` are a pair; `config_cnmf.R` copies `HIER_PROJ_DIR`).
- When adding/editing a script, keep the existing header comment style: purpose, INPUT, OUTPUT, how invoked, WHY.
  This is the project's only documentation.
- Before rerunning any stage, check whether the script has a `--force` flag and "skip if output exists" logic —
  most stages are resume-safe.
- **Rerunning `50_consensus_malignancy.R` invalidates everything downstream** (03/04/05/06/07/08) and changes
  the existing labels. Read the open item under premise 1 before touching it.
