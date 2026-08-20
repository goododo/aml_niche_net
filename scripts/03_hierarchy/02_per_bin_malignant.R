#!/usr/bin/env Rscript
# 02_per_bin_malignant.R ----
# Phase 2.2 -- the quantitative crux of the central hypothesis. Per sample, join the BMM projection
# per-cell table (cell -> hierarchy_bin, from 01) with the 02_malignancy consensus per-cell table
# (cell -> malignant), then aggregate the MALIGNANT FRACTION within each of the 8 hierarchy bins.
# This yields, for every sample, where malignant cells sit along the hematopoietic hierarchy -- the
# substrate for "LSC-like states converge on conserved bins" and "deviation grows Dx->MRD->Relapse".
#
# JOIN: exact on `cell` (both tables use the same QC-object barcodes: <sample>_<16bp>-1). The script
#   REPORTS the hit-rate per sample; a low rate means a format mismatch (investigate, don't trust).
#
# INPUT  : HIER_PROJ_DIR/<ds>/<sample>__bmm_percell.csv   (01 output)
#          DIR_MALIGNANCY/<ds>/<sample>__consensus_percell.csv  (02_malignancy output)
# OUTPUT : HIER_TAB_DIR/per_bin_malignant.csv     (long: one row per sample x bin)
#            dataset, sample, timepoint, hierarchy_bin, in_ccc_graph, n_cells, n_evaluable,
#            n_malignant, malignant_frac, frac_high_error
#          HIER_TAB_DIR/per_bin_malignant_joinQC.csv   (one row per sample: join hit-rate)
#
# NOTE: samples with no consensus percell (e.g. healthy references with no malignancy call) are
#   projection-only -> skipped here (they have no malignant layer, by design).

suppressPackageStartupMessages({ library(data.table); library(optparse) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--datasets", type = "character", default = "", help = "comma-separated dataset filter (default: all)"),
  make_option("--min_join", type = "double",    default = 0.90, help = "warn if join hit-rate below this [0.90]")
)))

dir.create(HIER_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# CURATED timepoint, not a name regex. The private tp_from_name() this replaces resolved 27 of 214
# samples -- all GSE227903 -- and returned NA for 34 post-treatment and 19 relapse samples in five
# other datasets, which 03 then dropped without a message. See TP_AXIS_LEVELS in config_hierarchy.R.
TPMAP <- sample_timepoints()

## ---- samples that have BOTH a projection and a consensus per-cell table ----
proj <- data.table(path = list.files(HIER_PROJ_DIR, pattern = "__bmm_percell\\.csv$", recursive = TRUE, full.names = TRUE))
proj[, sample := sub("__bmm_percell\\.csv$", "", basename(path))][, dataset := basename(dirname(path))]
cons <- data.table(path = list.files(DIR_MALIGNANCY, pattern = "__consensus_percell\\.csv$", recursive = TRUE, full.names = TRUE))
cons[, sample := sub("__consensus_percell\\.csv$", "", basename(path))][, dataset := basename(dirname(path))]
man <- merge(proj[, .(dataset, sample, proj_path = path)],
            cons[, .(dataset, sample, cons_path = path)], by = c("dataset", "sample"))
if (nzchar(opt$datasets)) man <- man[dataset %in% trimws(strsplit(opt$datasets, ",")[[1]])]
stopifnot(nrow(man) > 0)
message(sprintf("[0] %d sample(s) with BOTH projection + consensus (of %d projected, %d consensus)",
                nrow(man), nrow(proj), nrow(cons)))

bins_all <- HIERARCHY_BINS
per_bin <- vector("list", nrow(man)); join_qc <- vector("list", nrow(man))

for (i in seq_len(nrow(man))) {
  ds <- man$dataset[i]; sid <- man$sample[i]
  p <- fread(man$proj_path[i]); c <- fread(man$cons_path[i], select = c("cell", "malignant"))

  m <- merge(p[, .(cell, hierarchy_bin, in_ccc_graph, high_error)], c, by = "cell")   # exact join
  # BMM masks the label of cells that fail its mapping QC -> empty bmm_broad -> empty/NA bin here.
  # Make that explicit as 'Unassigned' (NOT a CCC node); such cells also can't be validly called
  # malignant (they don't fit the hematopoietic scaffold), same caveat as Stromal.
  m[is.na(hierarchy_bin) | hierarchy_bin == "", `:=`(hierarchy_bin = "Unassigned", in_ccc_graph = FALSE)]
  hit <- nrow(m) / max(nrow(p), 1L)
  join_qc[[i]] <- data.table(dataset = ds, sample = sid, n_proj = nrow(p), n_cons = nrow(c),
                             n_joined = nrow(m), join_rate = round(hit, 4))
  if (hit < opt$min_join) message(sprintf("  [WARN] %s::%s low join rate %.3f (proj=%d cons=%d joined=%d)",
                                           ds, sid, hit, nrow(p), nrow(c), nrow(m)))
  if (!nrow(m)) next

  # ensure every bin present (0-fill) so downstream sees a complete vector per sample
  agg <- m[, .(n_cells = .N,
               n_evaluable = sum(!is.na(malignant)),
               n_malignant = sum(malignant == 1, na.rm = TRUE),
               frac_high_error = round(mean(high_error), 4)), by = hierarchy_bin]
  allbins <- c(bins_all, "Unassigned")
  full <- agg[data.table(hierarchy_bin = allbins), on = "hierarchy_bin"]   # right join -> all bins present
  full[is.na(n_cells), `:=`(n_cells = 0L, n_evaluable = 0L, n_malignant = 0L, frac_high_error = NA_real_)]
  full[, in_ccc_graph := !(hierarchy_bin %in% c("Stromal", "Unassigned"))]  # authoritative: only these two are non-nodes
  full[, malignant_frac := fifelse(n_evaluable > 0, round(n_malignant / n_evaluable, 4), NA_real_)]
  .tp <- TPMAP[dataset == ds & sample == sid]
  full[, `:=`(dataset = ds, sample = sid,
              timepoint = if (nrow(.tp)) .tp$timepoint[1] else NA_character_,
              tp_axis   = if (nrow(.tp)) .tp$tp_axis[1]   else NA_character_,
              patient   = if (nrow(.tp)) .tp$uid_patient[1] else NA_character_)]
  per_bin[[i]] <- full[, .(dataset, sample, timepoint, tp_axis, patient, hierarchy_bin, in_ccc_graph,
                           n_cells, n_evaluable, n_malignant, malignant_frac, frac_high_error)]
}

pb <- rbindlist(per_bin, fill = TRUE)
jq <- rbindlist(join_qc, fill = TRUE)
# A silent NA here is what shrank the longitudinal cohort from 28 patients to 6.
.cov <- mean(!is.na(pb$timepoint))
if (.cov < 0.95) stop(sprintf("curated timepoint covers only %.1f%% of per-bin rows", 100 * .cov))
message(sprintf("[timepoint] curated for %.1f%% of rows | axis: %s", 100 * .cov,
                paste(sprintf("%s=%d", names(table(pb$tp_axis)), table(pb$tp_axis)), collapse = " ")))
fwrite_safe(pb, file.path(HIER_TAB_DIR, "per_bin_malignant.csv"))
fwrite_safe(jq, file.path(HIER_TAB_DIR, "per_bin_malignant_joinQC.csv"))

message("[done] per-bin malignant -> ", file.path(HIER_TAB_DIR, "per_bin_malignant.csv"))
cat("\n== join hit-rate (should be ~1.0; low = barcode format mismatch) ==\n")
print(jq[order(join_rate)][1:min(12, .N)])
cat("\n== malignant fraction by bin, pooled over samples (evaluable-weighted) ==\n")
cat("   (Stromal & Unassigned are NON-hematopoietic / QC-fail: their 'malignant' is an inferCNV\n")
cat("    artifact -- reference is hematopoietic -- so ccc=FALSE rows are EXCLUDED from interpretation)\n")
pooled <- pb[, .(n_evaluable = sum(n_evaluable), n_malignant = sum(n_malignant)), by = hierarchy_bin]
pooled[, `:=`(malignant_frac = round(n_malignant / pmax(n_evaluable, 1), 4),
              ccc = !(hierarchy_bin %in% c("Stromal", "Unassigned")))]
print(pooled[order(-ccc, -malignant_frac)])
cat("\n-- CCC bins only (this is the interpretable LSC-distribution signal) --\n")
print(pooled[ccc == TRUE][order(-malignant_frac), .(hierarchy_bin, n_evaluable, n_malignant, malignant_frac)])
cat("\n== sanity: Stromal-bin malignant fraction per sample (HIGH => AML cells misprojecting to Stromal, e.g. 1886_R) ==\n")
print(pb[hierarchy_bin == "Stromal" & n_evaluable >= 20][order(-malignant_frac)][1:min(10, .N),
        .(dataset, sample, n_cells, n_malignant, malignant_frac)])