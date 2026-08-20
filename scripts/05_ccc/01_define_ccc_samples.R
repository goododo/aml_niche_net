# 01_define_ccc_samples.R ----
# Define the CCC-eligible sample set and per-node presence, BEFORE running any CellChat.
# Gate is METADATA-DRIVEN, not composition-based [R9]:
#   eligible = projected(has bmm) & l2_capable(study) & NOT sorted-sublibrary(Chen CD34).
# A lymphoid-cell floor is deliberately NOT a gate: blast-dominated whole-MNC AML has ~0 lymphoid
# cells yet is valid; its missing nodes are biological (kept, unbalanced FGW). Sorted CD34 missingness
# is technical (excluded). Composition cannot distinguish them -> use protocol metadata only.
# INPUT  : CCC_ROLE_MANIFEST (01_preprocess) ; CCC_BMM_DIR/<ds>/<sample>__bmm_percell.csv
# OUTPUT : DIR_CCC/ccc_sample_manifest.csv   (one row/sample: eligibility + reason + node summary)
#          DIR_CCC/ccc_node_presence.csv     (one row/sample x node: n_cells + present flag)
# Usage  : Rscript scripts/05_ccc/01_define_ccc_samples.R [--manifest=<path>] [--force]
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(here)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--manifest", type = "character", default = CCC_ROLE_MANIFEST),
  make_option("--force",    action = "store_true", default = FALSE)
)))

out_manifest <- file.path(DIR_CCC, "ccc_sample_manifest.csv")
out_presence <- file.path(DIR_CCC, "ccc_node_presence.csv")
# FRESHNESS, not existence -- see is_stale() in utils.R. Existence-only guards pinned this whole
# chain to the July cohort: re-running it printed "[skip]" five times and exited 0 while
# patient_scores.csv still named 55 samples that had left the cohort.
.ins <- c(opt$manifest,
          list.files(CCC_BMM_DIR, pattern = "__bmm_percell\\.csv$", recursive = TRUE, full.names = TRUE))
if (!is_stale(c(out_manifest, out_presence), .ins, force = opt$force)) {
  message("[skip] CCC sample manifest is current"); quit(status = 0)
}
if (file.exists(out_manifest)) message("[recompute] ", stale_reason(c(out_manifest, out_presence), .ins, force = opt$force))

## -- Step 1. per-(dataset,sample,bin) cell counts from BMM projection ----
# Only cells that will actually seed CCC: in_ccc_graph & !high_error.
message("[1] scanning BMM per-cell files under ", CCC_BMM_DIR)
bmm_files <- list.files(CCC_BMM_DIR, pattern = "__bmm_percell\\.csv$",
                        recursive = TRUE, full.names = TRUE)
message("[1] found ", length(bmm_files), " projected samples")
stopifnot(length(bmm_files) > 0L)

node_counts <- rbindlist(lapply(bmm_files, function(f) {
  d <- fread(f, select = c("cell","dataset","sample","hierarchy_bin","in_ccc_graph","high_error"))
  d[, `:=`(in_ccc_graph = as.logical(in_ccc_graph), high_error = as.logical(high_error))]
  d[in_ccc_graph == TRUE & high_error == FALSE & hierarchy_bin %in% CCC_NODES,
    .(n_cells = .N), by = .(dataset, sample, hierarchy_bin)]
}))

## -- Step 2. per-node presence table (fixed node vocabulary; absent bins -> n_cells=0) ----
smp_keys <- unique(node_counts[, .(dataset, sample)])
grid <- CJ(idx = seq_len(nrow(smp_keys)), hierarchy_bin = CCC_NODES)
grid <- cbind(smp_keys[grid$idx], hierarchy_bin = grid$hierarchy_bin)
presence <- merge(grid, node_counts, by = c("dataset","sample","hierarchy_bin"), all.x = TRUE)
presence[is.na(n_cells), n_cells := 0L]
presence[, present := n_cells >= CCC_MIN_CELLS_PER_NODE]
# OCCUPIED is a stricter bar than PRESENT: present (>=10 cells) means the node appears in the
# graph, occupied (>=30) means its mean expression is estimated from enough cells to be worth
# comparing across samples. The eligibility gate counts occupied nodes.
presence[, occupied := n_cells >= CCC_MIN_CELLS_PER_OCCUPIED_BIN]
setorder(presence, dataset, sample, hierarchy_bin)

## -- Step 3. per-sample rollups (total cells, nodes present, lymphoid count for eyeball) ----
wide <- dcast(presence, dataset + sample ~ hierarchy_bin, value.var = "n_cells", fill = 0)
roll <- presence[, .(total_cells     = sum(n_cells),
                     n_nodes_present = sum(present),
                     n_bins_occupied = sum(occupied)),
                 by = .(dataset, sample)]
roll <- merge(roll, wide[, .(dataset, sample,
                             lymphoid_cells = T_NK + B_Plasma,
                             immune_recv    = T_NK + B_Plasma + Mono_DC)],
              by = c("dataset","sample"))

## -- Step 4. eligibility gate (metadata) ----
message("[4] reading role manifest: ", opt$manifest)
if (!file.exists(opt$manifest)) stop("role manifest not found: ", opt$manifest)
man <- fread(opt$manifest)
# harmonize the sample key column name (manifest uses 'Sample')
if (!"sample" %in% names(man) && "Sample" %in% names(man)) setnames(man, "Sample", "sample")
man_cols <- intersect(c("dataset","sample","Study","Patient_ID","uid_patient",
                        "Timepoint","study_role","sample_role","is_longitudinal",
                        "subtype_stratum","l2_capable"), names(man))
man <- man[, ..man_cols]

m <- merge(roll, man, by = c("dataset","sample"), all.x = TRUE)   # keep all projected samples

# sorted-sublibrary flag (dataset -> regex, case-insensitive). Chen CD34 is the only case.
m[, sorted_sublib := FALSE]
for (ds in names(CCC_SORTED_SUBLIB_PATTERN)) {
  pat <- CCC_SORTED_SUBLIB_PATTERN[[ds]]
  m[dataset == ds & grepl(pat, sample, ignore.case = TRUE), sorted_sublib := TRUE]
}

# robust coercion: manifest l2_capable may be logical or "TRUE"/"FALSE"/"True"; NA -> FALSE for gating.
m[, l2_val := as.logical(as.character(l2_capable))]
m[is.na(l2_val), l2_val := FALSE]
m[, l2_flag := (!CCC_ELIGIBLE_REQUIRES_L2) | l2_val]
m[, occupancy_ok := n_bins_occupied >= CCC_MIN_OCCUPIED_BINS]
m[, ccc_eligible := l2_flag & !sorted_sublib &
                    (total_cells >= CCC_MIN_TOTAL_CELLS) & occupancy_ok]
# human-readable exclusion reason (first failing rule)
m[, exclude_reason := fifelse(ccc_eligible, "",
                       fifelse(!l2_flag,                       "not_l2_capable",
                       fifelse(sorted_sublib,                  "sorted_sublibrary",
                       fifelse(total_cells < CCC_MIN_TOTAL_CELLS, "below_min_cells",
                       fifelse(!occupancy_ok,                  "below_min_occupied_bins", "other")))))]

setcolorder(m, c("dataset","sample","ccc_eligible","exclude_reason",
                 "Timepoint","study_role","sample_role","l2_capable","sorted_sublib",
                 "total_cells","n_nodes_present","n_bins_occupied","lymphoid_cells"))
setorder(m, -ccc_eligible, dataset, sample)

## -- Step 5. write + report ----
fwrite_safe(m,        out_manifest)
fwrite_safe(presence, out_presence)
message("[5] wrote ", out_manifest)
message("[5] wrote ", out_presence)

n_elig <- m[ccc_eligible == TRUE, .N]
message("[done] eligible samples: ", n_elig, " / ", nrow(m), " projected")
message("       excluded breakdown:")
print(m[ccc_eligible == FALSE, .N, by = exclude_reason])
message("       eligible by dataset x timepoint:")
print(dcast(m[ccc_eligible == TRUE], dataset ~ Timepoint, value.var = "sample", fun.aggregate = length))
