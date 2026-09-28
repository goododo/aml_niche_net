# Handoff: split-half noise floor for the `node_feat` ladder rung

**Written** 2026-09-24. **For** a fresh Claude Code session in `aml_niche_net`.
**Status** not started. The rung itself is built and scored; only its noise floor is missing.

---

## 0. The one-sentence job

`node_feat` won GATE 1 and cannot yet be believed, because it is the only winning arm with no
split-half retest floor. Build that floor, the same way every other arm has one.

---

## 1. What already exists (do not rebuild any of this)

| thing | where |
|---|---|
| the rung's distance matrices | `results/tables/06_distance/scaccordion_distance__node_feat{,150}.csv` |
| the script that built them | `scripts/06_distance/04_node_feature_distance.py` |
| its registration | `scripts/08_scoring/PREREGISTRATION_scaccordion.md` **AMENDMENT 1** |
| GATE 1 scores, all 8 arms | `results/tables/08_scoring/scaccordion_gate1.csv` |
| the floor machinery for the other arms | `scripts/08_scoring/16_scaccordion_noise_floor.py` |
| split-half **cell lists**, 37 samples x 2 | `results/tables/05_ccc/split_half/<ds>__<sample>__<A\|B>.csv` (74 files + `tasks.csv`) |

### The result this exists to qualify

| arm | contrast | median pct | p_raw | p_BH | floor |
|---|---|---|---|---|---|
| `node_feat` | Dx_to_Relapse | 0.043 | 0.000488 | 0.0039 | **missing** |
| `node_feat150` | Dx_to_Relapse | 0.045 | 0.000488 | 0.0039 | **missing** |
| `prop7` | Dx_to_Relapse | 0.174 | 0.0049 | 0.026 | n/a by design |
| every communication arm | Dx_to_Relapse | 0.125–0.130 | 0.022–0.035 | 0.068–0.070 | 2.18–2.35 |

`p_raw = 0.000488` is the minimum attainable at n = 11. The registered reading in
`17_scaccordion_gate1.py` is explicit: **`d_true_over_floor < 1` means the cell cannot be believed
whatever its p.** Without a floor, `node_feat` has no such number.

`prop7` is excluded from `16` for a stated reason that does **not** transfer: `05_ccc/06` splits
cells **stratified by hierarchy_bin**, so both halves have identical composition by construction and
prop7's retest distance would be ~0 by design. Node features are bin-level *means over different
cells within* each bin, so they do vary between halves. That is why this floor is computable and
prop7's is not.

---

## 2. Why this is not a one-liner (read before estimating)

The obvious route — add `--cell_subset` to `05_ccc/03_node_features.R` — does not work:

1. **`03_node_features.R` is cohort-level.** It loops the whole manifest in one process
   (`rbindlist(lapply(seq_len(nrow(man)), ...))`) and writes one file. There is no per-sample entry
   point to subset.
2. **The scaling has to match production, and it is not persisted.** The distances in §1 were
   computed on `results/tables/07_fgw/fgw_nodes_long.csv`, which is **globally z-scored** over the
   966-row cohort by `07_fgw/01_build_fgw_inputs.R` (`FGW_FEATURE_SCALE <- "global_z"`, NAs
   mean-imputed first). Half-sample features must be scaled with **those same mu/sigma**, not with
   parameters re-derived from the halves. Otherwise numerator and denominator of
   `d_true / floor` live on different scales and the ratio is meaningless.

Both problems are avoidable. The route below does not touch `03_node_features.R` at all.

---

## 3. The route, with its two proof gates

Write **one** new script. Proposed: `scripts/06_distance/05_node_feature_noise_floor.py` — but
module placement is the user's call, ask first (project rule: propose the file layout and wait).

### Step 1 — a focused aggregator for the 47 all-cells features

Only the base stratum is needed. No `_normal` / `_malignant`, no van Galen, no CNV burden — that is
most of what makes `03_node_features.R` long.

Replicate, from `scripts/05_ccc/03_node_features.R` lines 175-205 and `config/config_ccc.R`:

- `CCC_PANELS` gives, per panel key, the directory constant, the filename suffix, and the columns:
  `st` → `__stemness_percell.csv` (4), `pg` → `__progeny_percell.csv` (14),
  `cs` / `mt` → `__cellstate_percell.csv` (13 / 6), `pt` → `__bmm_percell.csv` (1),
  `mp` → `__mp_usage.csv` (9). Total 47.
- Column naming is `paste0(pk, "_", gsub("[^A-Za-z0-9_]", "_", col))` — e.g. `JAK-STAT` → `pg_JAK_STAT`.
- **The `pt` special case is load-bearing**: cells whose `bmm_broad` is in `BMM_PSEUDOTIME_OFFTRAJ`
  are set to NA *before* averaging. The reference pins six terminally differentiated classes at
  exactly 0, and 30.7% of it sits there. Averaging over them inverts the column's meaning.
- Cell set and bins come from the **reconciled** per-cell table, with `in_ccc_graph & !high_error &
  hierarchy_bin %in% CCC_NODES`, and `CCC_EXCLUDE_FINE = "Early Lymphoid"` removed.

**PROOF GATE 1 — do not proceed past this.** Run the aggregator on the **full** cell set of the 37
paired samples and reproduce `results/tables/05_ccc/ccc_node_features.csv`'s 47 all-cells columns.
Require max |difference| < 1e-9 on every (sample, bin, feature) cell. If it does not reproduce, the
aggregator is wrong and everything downstream is a complete-looking wrong table. This is the same
discipline as `12_pseudobulk_de`'s "449 / 449 by an independent path" check.

### Step 2 — apply the production scaling

Recompute mu / sigma per feature exactly as `07_fgw/01_build_fgw_inputs.R` lines ~130-146 does:
mean-impute NAs, then mean and sd over the full 966-row cohort table. Apply to the half tables.

**PROOF GATE 2.** Re-scale the **production** raw table with those same mu / sigma and reproduce
`fgw_nodes_long.csv`'s 47 columns to < 1e-9. Only then are the half values on the right scale.

### Step 3 — retest distances

Per sample: build half-A and half-B vectors (7 bins x 47, `CCC_NODES` order, flattened), Euclidean
distance between them. 37 numbers.

### Step 4 — extend `16_scaccordion_noise_floor.py`

Add `node_feat` to its `ARMS`, matching how it reports the existing five. Keep its existing
structure; it already carries the "both terms computed at half depth" convention in its header,
which applies here unchanged.

### Step 5 — rerun `17_scaccordion_gate1.py`

It reads the floor from `results/tables/08_scoring/scaccordion_noise_floor.csv` and will fill
`retest_floor` and `d_true_over_floor` for the rung. Nothing else in that script changes.

---

## 4. What the answer licenses

Fix this before seeing the number, in the amendment block, the way AMENDMENT 1 fixed its prediction:

- **`d_true_over_floor` comfortably above 1** (the communication arms sit at 2.18–2.35, so that is
  the natural comparison) — the GATE 1 win stands, and the ladder's ordering becomes
  *node state > composition > communication*, all without transport.
- **At or below 1** — the win is the 329-dimensional representation remembering measurement noise,
  not patient biology. The rung is reported as uninterpretable and `prop7` remains the best
  believable arm.
- Either way, the floor is computed **at half depth on both terms**, so it is conservative in the
  same direction as every other arm's, and that must be said rather than assumed.

---

## 5. Two follow-ups this unblocks (do not do them in the same session)

1. **`FINDINGS_project_status.md` §4 and `FINDINGS_topology_null.md` currently assert "in the limit,
   counting cells, are sufficient". That is refuted by §1 above** and must be narrowed to the
   wording AMENDMENT 1 fixed in advance. Both files had uncommitted edits from another session on
   2026-09-24 — check with the user before editing, to avoid a collision.
2. `node_feat150` tracks `node_feat` almost exactly (0.045 vs 0.043, identical p). The extra 101
   features add nothing, which is worth one line in the findings and no further work.

---

## 6. Operational

- **`main` only.** `scripts/99_admin/daily_commit.sh` runs `git add -A` on the *current* branch at
  23:30 and pushes. A working branch checked out at that moment swallows the day's output.
- **Never `Read` or `cat` anything under `/LARGE1/`.** The per-cell score CSVs live there; loading
  them inside R or Python is fine, reading them into the conversation is not.
- Run R through `ENV_PREFIX` (`/FAST/gr10634/gaozy/general_env`); the login-node `Rscript` is R 4.5
  and cannot see the project's libraries. Python: `/FAST/gr10634/gaozy/general_env/bin/python`.
- `SEED = 491638`, from `scripts/config/config_paths.R`. Paths come from config, never hardcoded.
- Script headers are this repo's only documentation: purpose / INPUT / OUTPUT / how invoked / WHY.

## 7. The failure mode to watch for

From `FINDINGS_topology_null.md` §14: **a new arm silently scoring the old model and exiting 0.**
If the half-sample features come out equal to the full-sample features, or the floor comes out
suspiciously close to another arm's, that is a tell and not a result. Both proof gates in §3 exist
to catch exactly this before it reaches a p-value.
