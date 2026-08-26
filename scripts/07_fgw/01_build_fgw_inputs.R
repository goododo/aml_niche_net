# 01_build_fgw_inputs.R ----
# Phase 6 (stage 07_fgw), R side. Assemble per-sample FGW inputs from 06 (edge distance) + 03 (node
# features) into TIDY LONG CSVs that 02_fgw_align.py reads and pivots to arrays (C 7x7 directed, F 7x3
# z-scored, p 7 masses). CSV handoff = no R<->Python bridge dependency.
#
# LOCKED assembly rules (config_fgw.R):
#   C  : edge_distance C column (already 49 directed edges incl. self-loops). DIAGONAL KEPT = self-loop
#        rank-distance (autocrine, e.g. HSC_MPP->HSC_MPP CD99 is real). [decision ii] Python builds the
#        7x7 by reindexing on FGW_NODES; no diagonal override.
#   F  : frac_malignant / mean_stemness / n_cells. Healthy frac_malignant forced 0 (FGW_ZERO_HEALTHY_MAL).
#        z-scored GLOBALLY across all node x sample rows (FGW_SCALE_FEATURES) so the 3 disparate scales are
#        comparable (n_cells won't dominate); MUST be global (per-sample scaling erases between-sample
#        signal). Missing-node feature NA -> imputed to global mean BEFORE z-score (neutral -> 0 post-scale).
#   p  : from RAW n_cells (NOT the z-scored feature). FGW_MASS_MODE "ncells" ~ n_cells / "uniform" = present
#        equal. Absent node -> FGW_EPS_MASS. Renormalized to sum 1 within sample.
#   Only has_graph==TRUE samples assembled. sparse_flag carried (02 excludes from barycenter but still
#   computes pairwise HDS/ATS).
#
# INPUT  : DIR_DISTANCE/edge_distance.csv , DIR_DISTANCE/edge_qc.csv , DIR_CCC/ccc_node_features.csv ,
#          DIR_CCC/ccc_sample_manifest.csv (Disease_state/Timepoint for Healthy zeroing)
# OUTPUT : DIR_FGW/fgw_nodes_long.csv  (per sample x node: scaled feats + mass + present + flags)
#          DIR_FGW/fgw_edges_long.csv  (per sample x directed edge: C)
#          DIR_FGW/fgw_input_index.csv (per sample: timepoint, sparse_flag, healthy)
# Usage  : Rscript scripts/07_fgw/01_build_fgw_inputs.R [--force]
suppressPackageStartupMessages({ library(optparse); library(data.table); library(here); library(jsonlite) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_fgw.R"))   # FGW_* ; pulls CCC_NODES, DIR_DISTANCE, DIR_CCC
source(here::here("scripts", "config", "utils.R"))

# --out_dir and --features exist so an ALTERNATIVE FEATURE SET can be built and scored without
# touching the production results. The Python side was already parameterised (08_scoring/* take
# --root, 07_fgw/02 takes --input_dir); this stage was the only one that could not be redirected,
# which is why every previous feature question had to be answered from 08/07's decomposition alone
# -- and that only reports the alpha=0 term, never the alpha sweep or the barycenter.
opt <- parse_args(OptionParser(option_list = list(
  make_option("--force",    action = "store_true", default = FALSE),
  make_option("--out_dir",  type = "character", default = DIR_FGW,
              help = "where to write; default the production DIR_FGW"),
  make_option("--features", type = "character", default = "",
              help = "comma-separated feature set to put IN the FGW distance; default FGW_FEATURES")
)))

# A feature named here must exist in ccc_node_features.csv and must not be silently dropped later.
FEATURES <- if (nzchar(opt$features)) trimws(strsplit(opt$features, ",")[[1]]) else FGW_FEATURES
if (!identical(FEATURES, FGW_FEATURES))
  message("[cfg] NON-DEFAULT feature set: ", paste(FEATURES, collapse = " + "),
          "\n      (production set is ", paste(FGW_FEATURES, collapse = " + "), ")")
if (!identical(normalizePath(opt$out_dir, mustWork = FALSE), normalizePath(DIR_FGW, mustWork = FALSE)))
  message("[cfg] writing to ", opt$out_dir, " -- production DIR_FGW untouched")
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

out_index <- file.path(opt$out_dir, "fgw_input_index.csv")
out_nodes <- file.path(opt$out_dir, "fgw_nodes_long.csv")
out_edges <- file.path(opt$out_dir, "fgw_edges_long.csv")
# FRESHNESS, not existence. This guard printed "[skip]" and exited 0 on a superseded cohort:
# results/tables/07_fgw/patient_scores.csv holds 148 rows of which 55 name samples that have
# left the cohort, and 47 current samples have never entered CCC at all. Re-running the chain
# hit five of these guards in a row and reported success.
# the exact four files read at lines below -- a guard that names a file the stage does not read
# certifies nothing, and a non-existent path is silently ignored by is_stale()
.ins <- c(file.path(DIR_DISTANCE, "edge_distance.csv"),
          file.path(DIR_DISTANCE, "edge_qc.csv"),
          file.path(DIR_CCC, "ccc_node_features.csv"),
          file.path(DIR_CCC, "ccc_sample_manifest.csv"))
if (!is_stale(c(out_index, out_nodes, out_edges), .ins, force = opt$force)) {
  message("[skip] FGW inputs are current"); quit(status = 0)
}
if (file.exists(out_index)) message("[recompute] ", stale_reason(c(out_index, out_nodes, out_edges), .ins, force = opt$force))

## -- Step 1. load inputs ----
message("[1] loading edge distance / node features / qc / manifest")
ed  <- fread(file.path(DIR_DISTANCE, "edge_distance.csv"))
qc  <- fread(file.path(DIR_DISTANCE, "edge_qc.csv"))
nf  <- fread(file.path(DIR_CCC, "ccc_node_features.csv"))
man <- fread(file.path(DIR_CCC, "ccc_sample_manifest.csv"))
if (!"sample" %in% names(man) && "Sample" %in% names(man)) setnames(man, "Sample", "sample")

ds_col <- if ("Disease_state" %in% names(man)) "Disease_state" else "Timepoint"
man_ds <- unique(man[, .(dataset, sample, healthy = get(ds_col) == "Healthy")])
qc_key <- unique(qc[, .(dataset, sample, sparse_flag)])

samples <- unique(nf[has_graph == TRUE, .(dataset, sample)])
setorder(samples, dataset, sample)
message("[1] samples with graph: ", nrow(samples))

## -- Step 2. node MASS from RAW n_cells (before any scaling) ----
mass <- nf[has_graph == TRUE, .(dataset, sample, hierarchy_bin, n_cells_raw = n_cells)]
mass[, present := n_cells_raw >= 1L]
if (FGW_MASS_MODE == "ncells")  mass[, m := as.numeric(n_cells_raw)] else
if (FGW_MASS_MODE == "uniform") mass[, m := as.numeric(present)] else stop("bad FGW_MASS_MODE")
mass[present == FALSE | m == 0, m := FGW_EPS_MASS]
mass[, mass := m / sum(m), by = .(dataset, sample)]     # normalize to 1 within sample

## -- Step 3. GLOBAL feature z-score (Healthy frac_malignant -> 0 first) ----
feat <- merge(nf[has_graph == TRUE], man_ds, by = c("dataset","sample"), all.x = TRUE)

# ZEROING APPLIES TO frac_malignant ONLY. frac_malignant_vg is deliberately left untouched: it is the
# non-circular counterpart, and zeroing it for healthy would rebuild the very circularity that
# 08_scoring/07 exists to test. Keep the two apart.
if (FGW_ZERO_HEALTHY_MAL) feat[healthy == TRUE, frac_malignant := 0]

# Candidates ride the same scaling path but never reach FGW (see config_fgw.R).
missing_cand <- setdiff(FGW_CANDIDATE_FEATURES, names(feat))
if (length(missing_cand))
  stop("ccc_node_features.csv lacks candidate feature column(s): ", paste(missing_cand, collapse = ", "),
       "\n  Rebuild scripts/05_ccc/03_node_features.R first.")
ALL_FEATURES <- unique(c(FEATURES, FGW_CANDIDATE_FEATURES))

zpar <- list()
for (fcol in ALL_FEATURES) {
  mu <- mean(feat[[fcol]], na.rm = TRUE); sdv <- sd(feat[[fcol]], na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) sdv <- 1
  n_na <- sum(is.na(feat[[fcol]]))                      # report it: a candidate that is 90% imputed
  feat[is.na(get(fcol)), (fcol) := mu]                  # neutral impute (pre-scale)
  if (FGW_SCALE_FEATURES) feat[, (fcol) := (get(fcol) - mu) / sdv]
  zpar[[fcol]] <- c(mean = mu, sd = sdv, n_imputed = n_na)
}
message("[3] global z-score params:")
for (fcol in ALL_FEATURES)
  message("    ", if (fcol %in% FEATURES) "[use] " else "[cand] ", fcol,
          ": mean=", round(zpar[[fcol]]['mean'], 4), " sd=", round(zpar[[fcol]]['sd'], 4),
          " imputed=", zpar[[fcol]]['n_imputed'], "/", nrow(feat))

## -- Step 4. assemble long tables ----
nodes_long <- merge(
  feat[, c("dataset","sample","timepoint","hierarchy_bin","healthy", ALL_FEATURES), with = FALSE],
  mass[, .(dataset, sample, hierarchy_bin, n_cells_raw, present, mass)],
  by = c("dataset","sample","hierarchy_bin"))
nodes_long <- merge(nodes_long, qc_key, by = c("dataset","sample"), all.x = TRUE)
# enforce node order via a factor so downstream row order is deterministic
nodes_long[, hierarchy_bin := factor(hierarchy_bin, levels = FGW_NODES)]
setorder(nodes_long, dataset, sample, hierarchy_bin)

edges_long <- ed[, .(dataset, sample, sender_bin, receiver_bin, C)]
edges_long[, `:=`(sender_bin   = factor(sender_bin,   levels = FGW_NODES),
                  receiver_bin = factor(receiver_bin, levels = FGW_NODES))]
setorder(edges_long, dataset, sample, sender_bin, receiver_bin)

index <- unique(nodes_long[, .(dataset, sample, timepoint, sparse_flag,
                               healthy, mass_mode = FGW_MASS_MODE)])

## -- VALIDATE BEFORE WRITING ANYTHING ---------------------------------------------------------
# Every timepoint present must be a label CANONICAL_TIMEPOINTS covers. If curation introduces a
# label nobody wired up -- or, as happened here, the UPSTREAM input still carries labels that were
# RETIRED on 2026-08-04 -- that is an error, not 34 silently deleted rows six stages downstream.
# This runs before the writes on purpose: the first version of this check sat after them and left
# a complete, well-formed, stale-vocabulary output set on disk next to a non-zero exit code, which
# is exactly the state a resume-guard would then treat as done.
seen_tp <- setdiff(unique(as.character(index$timepoint)), NA_character_)
unknown_tp <- setdiff(seen_tp, CANONICAL_TIMEPOINTS)
if (length(unknown_tp))
  stop(sprintf(paste0("timepoint(s) not in CANONICAL_TIMEPOINTS: %s\n",
                      "  These are retired or unrecognised labels carried by the INPUT, so the node\n",
                      "  features are older than the 2026-08-04 vocabulary migration. Rebuild\n",
                      "  05_ccc/03_node_features.R before rebuilding FGW inputs; writing them now\n",
                      "  would bake a dead vocabulary into every downstream stage."),
               paste(unknown_tp, collapse = ", ")))

## -- Step 5. write ----
fwrite_safe(nodes_long, out_nodes)
fwrite_safe(edges_long, out_edges)
fwrite_safe(index,      out_index)

## -- EMIT THE VOCABULARY THE PYTHON STAGES MUST USE ------------------------------------------
# The R side derives FGW_BARY_GROUPS$aml from CANONICAL_TIMEPOINTS precisely so a vocabulary change
# cannot silently shrink B_AML, and config_fgw.R says so in a comment that names the failure:
# "Spelling out c(\"Diagnosis\",\"MRD\",\"Post_treatment\",...) is exactly how that would have
# happened". Nine Python files then did exactly that. CANONICAL_TIMEPOINTS was migrated on
# 2026-08-04 -- MRD and Post_treatment retired, six labels added -- and the Python mirror was not,
# so 34 of 214 samples (Post_induction 17, Post_treatment_unspecified 8, On_treatment 7,
# Post_consolidation 1, Refractory 1) fall to grp=="other" and are deleted by the AML/healthy
# filters in the H2 regression, both permutation pools, the alpha sweep, the edge regression, the
# feature decomposition and the CLR test. That is 64% of the treated arm -- the arm H3 is about --
# and nothing prints an "other" count.
# So the vocabulary is EMITTED here, next to the index it describes, and the Python loader aborts
# if it is missing rather than falling back to a literal. One source, no mirror to drift.
vocab <- list(
  aml_timepoints     = sort(setdiff(CANONICAL_TIMEPOINTS, c("Healthy", "Unknown"))),
  healthy_timepoints = "Healthy",
  excluded_timepoints= "Unknown",
  all_timepoints     = CANONICAL_TIMEPOINTS,
  nodes              = FGW_NODES,
  features           = FEATURES,
  # A candidate that has been promoted INTO the distance for this run must not also be listed as a
  # candidate, or 08/07 would test it as if it were still outside the model it is now inside.
  candidate_features = setdiff(FGW_CANDIDATE_FEATURES, FEATURES),
  generated_by       = "scripts/07_fgw/01_build_fgw_inputs.R"
)
out_vocab <- file.path(opt$out_dir, "fgw_vocab.json")
writeLines(jsonlite::toJSON(vocab, auto_unbox = TRUE, pretty = TRUE), out_vocab)

message(sprintf("[vocab] %s : %d AML timepoints, %d present in this index",
                out_vocab, length(vocab$aml_timepoints),
                length(intersect(seen_tp, vocab$aml_timepoints))))
message("[5] wrote ", out_nodes, "  (", nrow(nodes_long), " rows = ", nrow(samples), " x 7 nodes)")
message("[5] wrote ", out_edges, "  (", nrow(edges_long), " rows = ", nrow(samples), " x 49 edges)")
message("[5] wrote ", out_index, "  (", nrow(index), " samples)")
message("[5] barycenter-eligible (non-sparse): ", index[sparse_flag == FALSE | is.na(sparse_flag), .N],
        " | sparse (pairwise-only): ", index[sparse_flag == TRUE, .N])
message("[done]")
