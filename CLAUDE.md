# CLAUDE.md

Working context for Claude Code sessions in this repository.

## What this project is

Cross-cohort single-cell analysis of AML bone marrow: malignant-cell
annotation, hierarchy projection onto a healthy reference, and comparison of
cell–cell communication network topology across public cohorts (~339 samples,
spanning 10x 3'/5', BD Rhapsody, and Seq-Well).

Design document lives under `*_project_context/`.

## Environment

- Primary language: **R**. The project uses the `here` package (`.here` at root).
- Python where tooling requires it.
- Runs on a remote server over SSH.

---

## How to work in this repo

These are the rules that matter most. Follow them even when a different
approach would be technically cleaner.

### Never work directly on `main`

`main` holds work that already runs. All changes go on a branch.

### Incremental, never wholesale

Do **not** propose restructuring the repository, renaming directories, or
migrating to a new framework. Change the one thing currently being worked on
and leave everything else alone. A working messy pipeline beats a clean broken
one.

### Module boundaries follow pipeline stages, not software convention

Split along: `qc` / `align` / `cnv` / `projection` / `ccc`.

Do **not** create `utils.R`, `helpers.py`, `common/`, `core/`, or any
catch-all module. Those boundaries don't map onto anything in this project and
make the code harder to navigate, not easier.

### Plan before writing code

For anything beyond ~30 lines, describe the approach and the file layout
first, then wait for confirmation. Don't produce several interdependent files
in one go.

### No new dependencies without asking

If a new package is needed, say what it is, why, and what the alternative is.
Wait for an answer.

### Make outputs inspectable

Print dimensions and the first few rows of intermediate results after each
step. Prefer a runnable minimal example over a large block that can only be
verified by running the whole pipeline.

### Explain, don't just deliver

When using a language feature or idiom that isn't obvious, add a one-line
comment saying what it does. The goal here is not only working code — it is
code I understand well enough to defend in a technical interview.

---

## Daily working loop

When I open a script I do four things. Help with all four, in this order:

1. Write a three-line **English** header: `input` / `output` / `what it does`
2. Identify the 3–5 lines I don't understand — explain only those, not the
   whole script
3. Add **one** sanity check on the output
4. Commit with an English message

Do not expand this into a refactor. One script, four steps, done.

---

## Verification over reading

Correctness in this project is established by **checking outputs**, not by
reading code. When asked "is this right", reach for these first:

| Check | Example here |
| --- | --- |
| Order of magnitude | 8000 cells in, 200 after filtering? investigate |
| Known answer | Petti samples carry mutation labels — malignant calls should agree |
| Edge case | A sample with 3 malignant cells: does it fail loudly or silently? |
| Permutation | Shuffle labels, rerun — signal should disappear |

---

## Method decisions — do not silently change these

- **inferCNV is the primary CNV caller.** Numbat is *not* usable: several of
  the public datasets do not distribute usable raw data. Do not suggest
  switching back to it.
- **VarTrix supplements** on the subset of samples with known somatic SNVs
  (e.g. Petti, Riether), giving two independent lines of evidence there.
- **No cross-sample expression integration** before per-sample analysis, so
  batch structure is never baked into per-cell calls.
- **Malignancy and hierarchy are orthogonal axes.** Malignancy comes from
  genomic evidence; hierarchy comes from reference projection. Keep them
  separate — a cell can be `GMP-like + malignant` or `GMP-like + normal`.
- **"Microenvironment" here means the hematopoietic-immune and paracrine
  compartments**, not structural stroma. Structural stroma (MSC, endothelial)
  is cohort-level auxiliary evidence only. Don't phrase things as structural
  niche claims.
- Which evidence types were available per sample is **an output**, not a
  footnote — it feeds the Methods section directly.

---

## Known repo debt

Fix incrementally, only when already touching the relevant part. Do not
schedule a cleanup sprint.

- ~15 loose `.tsv` / `.csv` / `.txt` files at repo root belong under `config/`
- No README yet
- `infercnv_tasks*.tsv` (`_retry`, `_sorted`, `_sens_eval`,
  `_matchedref_eval`) are a **hand-maintained job queue** — this is the thing
  Snakemake will replace

---

## Snakemake migration (in progress)

Five steps. One at a time. Each must actually run before moving on, and each
gets its own commit.

1. One rule, one sample, hardcoded paths
2. Add `{sample}` wildcard; read the sample list from `infercnv_tasks.tsv`
3. Add the upstream rule so two stages chain
4. Move parameters into `config.yaml`
5. Per-rule environment declarations

**Do not scaffold all five at once.** A skeleton where every rule is a stub
teaches nothing and cannot be verified.

---

## Do not put in this repository

This repo is **public**. Never add personal, career, visa, salary, contact, or
collaborator-assessment information to any file here, including commit
messages. Keep that material outside the repository entirely.
