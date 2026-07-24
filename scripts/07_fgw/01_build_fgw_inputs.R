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
suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_fgw.R"))   # FGW_* ; pulls CCC_NODES, DIR_DISTANCE, DIR_CCC
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--force", action = "store_true", default = FALSE)
)))

out_index <- file.path(DIR_FGW, "fgw_input_index.csv")
out_nodes <- file.path(DIR_FGW, "fgw_nodes_long.csv")
out_edges <- file.path(DIR_FGW, "fgw_edges_long.csv")
if (all(file.exists(out_index, out_nodes, out_edges)) && !opt$force) {
  message("[skip] outputs exist; --force to rebuild"); quit(status = 0)
}

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
if (FGW_ZERO_HEALTHY_MAL) feat[healthy == TRUE, frac_malignant := 0]
zpar <- list()
for (fcol in FGW_FEATURES) {
  mu <- mean(feat[[fcol]], na.rm = TRUE); sdv <- sd(feat[[fcol]], na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) sdv <- 1
  feat[is.na(get(fcol)), (fcol) := mu]                  # neutral impute (pre-scale)
  if (FGW_SCALE_FEATURES) feat[, (fcol) := (get(fcol) - mu) / sdv]
  zpar[[fcol]] <- c(mean = mu, sd = sdv)
}
message("[3] global z-score params:")
for (fcol in FGW_FEATURES) message("    ", fcol, ": mean=", round(zpar[[fcol]]['mean'],4), " sd=", round(zpar[[fcol]]['sd'],4))

## -- Step 4. assemble long tables ----
nodes_long <- merge(
  feat[, c("dataset","sample","timepoint","hierarchy_bin","healthy", FGW_FEATURES), with = FALSE],
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

## -- Step 5. write ----
fwrite_safe(nodes_long, out_nodes)
fwrite_safe(edges_long, out_edges)
fwrite_safe(index,      out_index)
message("[5] wrote ", out_nodes, "  (", nrow(nodes_long), " rows = ", nrow(samples), " x 7 nodes)")
message("[5] wrote ", out_edges, "  (", nrow(edges_long), " rows = ", nrow(samples), " x 49 edges)")
message("[5] wrote ", out_index, "  (", nrow(index), " samples)")
message("[5] barycenter-eligible (non-sparse): ", index[sparse_flag == FALSE | is.na(sparse_flag), .N],
        " | sparse (pairwise-only): ", index[sparse_flag == TRUE, .N])
message("[done]")
