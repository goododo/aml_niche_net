#!/usr/bin/env Rscript
# d40_rollup.R ----
# Phase 2 PILOT close-out: consolidate per-sample x per-bin tables across all pilot samples into ONE
# cross-sample table -- the L1/L2 entry point. Pure aggregation: reads d30 per-bin + d25 stemness;
# adds full-bin stemness means and honest caution flags. Computes nothing new about malignancy.
#
# INPUTS per sample:
#   PERBIN_DIR/<ds>/<sample>_perbin.csv          (d30, run against the d35 inferCNV labels)
#   PERCELL_STEM_DIR/<ds>/<sample>_stem.csv       (d25; per-cell, z-scored WITHIN sample)
#   PERCELL_BINNED_DIR/<ds>/<sample>_binned.csv   (d20; for cell->bin to aggregate stemness)
#
# STEMNESS AGGREGATION = full-bin mean (DECISION 1): mean over ALL cells in a bin (malignant +
#   normal pooled). Chosen because the pilot malignant label is single-evidence inferCNV and not yet
#   reliable enough to split malignant/normal stemness; revisit (mean_stem_malignant/_normal) once
#   the three-way consensus is restored. NOTE: d25 scores are z-scored within sample, so bin means
#   are SAMPLE-RELATIVE (compare bins within a sample, not absolute values across samples).
#
# CAUTION flags (DECISION c; ';'-joined, empty = none):
#   caution_cnv_single             : sample where only ONE CNV method was valid (evidence_tier C);
#                                     the other method was blind (no_CNV_detected) or not run
#   caution_tierB_saturated        : GSE239721 (expr-only consensus bimodally saturated -- known limit)
#   caution_low_malignant_primitive: Dg/R sample whose HSC_MPP+LMPP_GMP pooled malignant_frac < 0.10
#                                     (contradicts expected blast enrichment in primitive bins at
#                                      diagnosis/relapse -> malignant label suspect for this sample)
#
# Usage:
#   Rscript d40_rollup.R --datasets=GSE227903,GSE239721
#   Rscript d40_rollup.R --datasets=GSE227903,GSE239721 --out=/path/rollup.csv

suppressPackageStartupMessages({ library(optparse); library(data.table) })
source("d00_config.R")

PERCELL_STEM_DIR   <- file.path(FAST_DIR, "results/tables/04_hierarchy/percell_stem")
PERCELL_BINNED_DIR <- file.path(FAST_DIR, "results/tables/04_hierarchy/percell_binned")
ROLLUP_DEFAULT     <- file.path(FAST_DIR, "results/tables/04_hierarchy/phase2_pilot_rollup.csv")

# tunables for the biological-plausibility caution
PRIMITIVE_BINS         <- c("HSC_MPP", "LMPP_GMP")
LOW_PRIMITIVE_MAL_FRAC <- 0.10
TIERB_SATURATED_DATASETS <- c("GSE239721")   # expr-only saturated (blueprint known limitation)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--datasets", type = "character", help = "comma-separated dataset ids"),
  make_option("--out", type = "character", default = ROLLUP_DEFAULT, help = "output csv path")
)))
if (is.null(opt$datasets)) stop("--datasets is required (e.g. GSE227903,GSE239721)")
datasets <- strsplit(opt$datasets, ",")[[1]]

## -- timepoint from sample name: trailing _Dg / _MRD / _R / _R2 ... ----
parse_timepoint <- function(samp) {
  m <- regmatches(samp, regexpr("_(Dg|MRD|R[0-9]*)$", samp))
  if (length(m)) sub("^_", "", m) else NA_character_
}

stem_cols <- NULL  # discovered from first stem file

## -- per-sample: attach full-bin stemness means to the d30 per-bin table ----
load_sample <- function(ds, samp) {
  pb_csv  <- file.path(PERBIN_DIR, ds, paste0(samp, "_perbin.csv"))
  st_csv  <- file.path(PERCELL_STEM_DIR, ds, paste0(samp, "_stem.csv"))
  bn_csv  <- file.path(PERCELL_BINNED_DIR, ds, paste0(samp, "_binned.csv"))
  if (!file.exists(pb_csv)) { message("    [miss perbin] ", ds, "/", samp); return(NULL) }
  pb <- fread(pb_csv)

  # full-bin stemness mean (pool all cells)
  if (file.exists(st_csv) && file.exists(bn_csv)) {
    st <- fread(st_csv); bn <- fread(bn_csv, select = c("cell", "hierarchy_bin"))
    sc <- grep("_score$", names(st), value = TRUE)
    if (is.null(stem_cols)) stem_cols <<- sc
    m  <- merge(bn, st[, c("cell", sc), with = FALSE], by = "cell")
    sm <- m[, lapply(.SD, mean, na.rm = TRUE), by = hierarchy_bin, .SDcols = sc]
    setnames(sm, sc, paste0("mean_", sc))
    pb <- merge(pb, sm, by = "hierarchy_bin", all.x = TRUE)
  } else {
    message("    [no stem] ", ds, "/", samp)
  }
  pb[, timepoint := parse_timepoint(samp)]
  pb[]
}

all_pb <- list()
for (ds in datasets) {
  bdir <- file.path(PERBIN_DIR, ds)
  pbf  <- list.files(bdir, pattern = "_perbin\\.csv$", full.names = TRUE)
  if (!length(pbf)) { message("[!] no perbin csv for ", ds, " (run d30 first)"); next }
  samps <- sub("_perbin\\.csv$", "", basename(pbf))
  message("[", ds, "] ", length(samps), " sample(s)")
  for (s in samps) all_pb[[paste(ds, s)]] <- load_sample(ds, s)
}
roll <- rbindlist(all_pb, fill = TRUE)
if (!nrow(roll)) stop("[d40] nothing to roll up")

## -- per-sample primitive malignant fraction (pooled HSC_MPP + LMPP_GMP) for the caution ----
prim <- roll[hierarchy_bin %in% PRIMITIVE_BINS,
             .(prim_called = sum(n_called, na.rm = TRUE),
               prim_mal    = sum(n_malignant, na.rm = TRUE)),
             by = .(Dataset, Sample)]
prim[, prim_mal_frac := fifelse(prim_called > 0, prim_mal / prim_called, NA_real_)]
roll <- merge(roll, prim[, .(Dataset, Sample, prim_mal_frac)], by = c("Dataset", "Sample"), all.x = TRUE)

## -- caution flags (';'-joined; per ROW, driven by sample-level facts) ----
roll[, caution := {
  out <- vector("character", .N)
  for (i in seq_len(.N)) {
    fl <- character(0)
    # single-method CNV evidence (evidence_tier starting with "C")
    if (!is.na(evidence_tier[i]) && grepl("^C", evidence_tier[i]))
      fl <- c(fl, "caution_cnv_single")
    if (Dataset[i] %in% TIERB_SATURATED_DATASETS) fl <- c(fl, "caution_tierB_saturated")
    tp <- timepoint[i]
    if (!is.na(tp) && grepl("^(Dg|R)", tp) &&
        !is.na(prim_mal_frac[i]) && prim_mal_frac[i] < LOW_PRIMITIVE_MAL_FRAC) {
      fl <- c(fl, "caution_low_malignant_primitive")
    }
    out[i] <- paste(fl, collapse = ";")
  }
  out
}]

## -- order & write ----
lead <- c("Dataset", "Sample", "timepoint", "hierarchy_bin", "in_ccc_graph",
          "n_cells", "n_called", "n_malignant", "malignant_frac",
          "malignant_frac_lo", "malignant_frac_hi", "n_malignant_high_error",
          "frac_high_error", "frac_robust", "evidence_tier", "bin_confidence",
          "prim_mal_frac", "caution")
stem_means <- grep("^mean_.*_score$", names(roll), value = TRUE)
setcolorder(roll, c(intersect(lead, names(roll)), sort(stem_means)))
setorder(roll, Dataset, Sample, hierarchy_bin)
fwrite_safe(roll, opt$out)

## -- console summary ----
message("\n[done] d40 rollup -> ", opt$out)
message("  rows (sample x bin): ", nrow(roll),
        " | samples: ", uniqueN(roll[, .(Dataset, Sample)]),
        " | datasets: ", uniqueN(roll$Dataset))
message("  bin_confidence: ",
        paste(sprintf("%s=%d", names(table(roll$bin_confidence, useNA = "ifany")),
                      as.integer(table(roll$bin_confidence, useNA = "ifany"))), collapse = "  "))
cf <- sort(table(unlist(strsplit(roll$caution, ";"))), decreasing = TRUE)
message("  caution tallies: ", paste(sprintf("%s=%d", names(cf), as.integer(cf)), collapse = "  "))
n_low_prim <- uniqueN(roll[grepl("caution_low_malignant_primitive", caution), .(Dataset, Sample)])
message("  samples flagged low_malignant_primitive: ", n_low_prim,
        "  (if a known-CNV Dg sample lands here, its inferCNV label is still under-calling)")
