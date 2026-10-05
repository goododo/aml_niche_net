#!/usr/bin/env Rscript
# 02_doublet_scores.R ----
# INPUT  : 02_seurat_objects/04_annotated/<ds>/<sample>.rds ; 00_curated_manifest.csv
#          results/tables/01_preprocess/03_qc_report__<ds>.csv  (n_raw, for the expected rate)
#          scripts/01_preprocess/03_per_sample_qc.R  (call_sc, call_doubletfinder -- REUSED, not rewritten,
#          so these scores stay comparable with the 28-sample calibration 07_doublet_calibration.R)
# OUTPUT : results/tables/ws/ws_doublet_scores_percell.csv.gz   cell, sc_score, sc_class, df_pANN, df_class
#          results/tables/ws/ws_doublet_libraries.csv           one row per library: rate used and its basis
# WHAT IT DOES : task T3a of HANDOFF_within_sample_v1.5. Both callers are run PER LIBRARY, never per
#          patient object, because a doublet is a property of the droplet emulsion. GSE185381 is
#          multiplexed both ways (32 of 42 eligible patients span 2-3 libraries; 27 of 47 libraries
#          hold 2-5 patients), so its cells are regrouped across patient objects into whole libraries
#          before either caller sees them.
# Usage  : Rscript scripts/ws/02_doublet_scores.R
#
# TWO THINGS THE FROZEN SPEC LEFT OPEN, settled here before any score exists:
#   1. WHICH N in "D = N * 0.008 * N/1000". The existing pipeline used the PRE-QC count: I verified
#      dbl_rate_exp == n_raw * 8e-6 exactly (sd/mean = 0.0000) across 03_qc_report__GSE185381.csv.
#      So n_raw it is. A regrouped library has no n_raw of its own, so it is estimated as
#      assembled_cells / median(n_final/n_raw) for that dataset, and rate_basis records that.
#   2. THE POOLED CORRECTION. demux already removed cross-donor doublets, so only same-donor ones
#      survive. The exact surviving fraction is sum(p_i^2) over the donors' shares of the library,
#      not 1/k. Both are recorded so the choice is visible.
# Nothing is written back into 04_annotated: the output is a side table keyed by cell barcode.

suppressPackageStartupMessages({ library(data.table); library(Seurat) })
args <- commandArgs(TRUE)
cfg <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "..", "config", "config_paths.R")
source(cfg); set.seed(SEED)
# call_doubletfinder lives in 03_per_sample_qc.R and is reused verbatim so the pANN/class match the
# 28-sample calibration. call_sc does NOT live there (it is defined in 07_doublet_calibration.R, a
# top-level analysis script that must not be sourced), so scDblFinder is called directly below --
# once, returning score and class together, instead of the two separate runs the first version did.
suppressMessages(source(file.path(SCRIPTS_DIR, "01_preprocess", "03_per_sample_qc.R")))
stopifnot(exists("call_doubletfinder"))        # fail loudly, not 95 times in a row

AN_DIR  <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated")
OUT_DIR <- file.path(TAB_DIR, "ws"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
STUDIES <- c("GSE185381", "GSE239721", "GSE116256", "GSE227903", "GSE289435")
MULTIPLEXED <- "GSE185381"              # the only study in this cohort with a Library column that matters
MIN_LIB_CELL <- 100L                    # below this DoubletFinder's pK sweep is not meaningful
RATE_PER_CELL <- 8e-6                   # 0.008 per 1000 recovered cells, the pipeline's own constant

man <- fread(file.path(TAB_DIR, "01_preprocess", "00_curated_manifest.csv"))
elig <- man[disease == "AML" & timepoint == "Diagnosis" & tissue_r == "BM" &
            sorting_r %in% c("viability_only", "none") & dataset %in% STUDIES, .(dataset, sample)]
message(sprintf("[1] eligible diagnostic samples: %d across %d studies", nrow(elig), uniqueN(elig$dataset)))

## -- per-dataset QC retention, for estimating a regrouped library's raw count -------------------
ret <- rbindlist(lapply(STUDIES, function(ds) {
  f <- file.path(TAB_DIR, "01_preprocess", paste0("03_qc_report__", ds, ".csv"))
  if (!file.exists(f)) return(NULL)
  q <- fread(f)
  data.table(dataset = ds, retention = median(q$n_final / q$n_raw, na.rm = TRUE))
}))
qcr <- rbindlist(lapply(STUDIES, function(ds) {
  f <- file.path(TAB_DIR, "01_preprocess", paste0("03_qc_report__", ds, ".csv"))
  if (!file.exists(f)) return(NULL)
  fread(f)[, .(dataset = ds, sample = Sample, n_raw)]
}))
print(ret)

## -- assemble libraries -------------------------------------------------------------------------
# A library is the unit. For a non-multiplexed dataset that is the sample itself, and the sample's
# own n_raw is the rate basis. For GSE185381 the library's cells are collected across patient objects.
LIBS <- list()
for (ds in unique(elig$dataset)) {
  rows <- elig[dataset == ds]
  if (!ds %in% MULTIPLEXED) {
    for (sm in rows$sample) LIBS[[paste(ds, sm, sep = "|")]] <-
      list(dataset = ds, library = sm, members = sm)
    next
  }
  # multiplexed: read metadata only, build library -> list of (patient object, barcodes)
  map <- rbindlist(lapply(rows$sample, function(sm) {
    f <- file.path(AN_DIR, ds, paste0(sm, ".rds")); if (!file.exists(f)) return(NULL)
    s <- readRDS(f); md <- s@meta.data
    r <- data.table(sample = sm, cell = rownames(md),
                    library = if ("Library" %in% names(md)) as.character(md$Library) else sm)
    rm(s); gc(FALSE); r
  }))
  for (lb in unique(map$library)) LIBS[[paste(ds, lb, sep = "|")]] <-
    list(dataset = ds, library = lb, members = unique(map[library == lb, sample]),
         cells = map[library == lb, cell])
  message(sprintf("    %s: %d libraries assembled from %d patient objects", ds, uniqueN(map$library), uniqueN(map$sample)))
}
message(sprintf("[2] libraries to score: %d", length(LIBS)))

## -- score ---------------------------------------------------------------------------------------
PC <- list(); LB <- list(); skipped <- character()
for (k in names(LIBS)) {
  L <- LIBS[[k]]; ds <- L$dataset
  obj <- tryCatch({
    parts <- lapply(L$members, function(sm) {
      f <- file.path(AN_DIR, ds, paste0(sm, ".rds")); if (!file.exists(f)) return(NULL)
      s <- readRDS(f)
      if (!is.null(L$cells)) { keep <- intersect(colnames(s), L$cells); s <- s[, keep] }
      DefaultAssay(s) <- "RNA"; s
    })
    parts <- parts[!vapply(parts, is.null, TRUE)]
    if (!length(parts)) NULL else if (length(parts) == 1L) parts[[1]] else merge(parts[[1]], parts[-1])
  }, error = function(e) { message("  assemble fail ", k, ": ", conditionMessage(e)); NULL })
  if (is.null(obj) || ncol(obj) < MIN_LIB_CELL) {
    skipped <- c(skipped, sprintf("%s(n=%s)", k, if (is.null(obj)) "NA" else ncol(obj))); next }
  # GSE239721 stores a Seurat v3 Assay; the other four studies store Assay5. JoinLayers only has a
  # method for Assay5, and a v3 assay has a single counts slot with nothing to join.
  if (inherits(obj[["RNA"]], "Assay5") && length(SeuratObject::Layers(obj[["RNA"]])) > 1L)
    obj <- suppressWarnings(JoinLayers(obj))
  n <- ncol(obj)

  # donor shares inside this library -> the exact post-demux surviving-doublet fraction
  shares <- if (!is.null(L$cells)) {
    tb <- table(obj@meta.data$Sample); as.numeric(tb) / sum(tb)
  } else 1
  same_donor <- sum(shares^2)
  if (is.null(L$cells)) {
    nraw <- qcr[dataset == ds & sample == L$library, n_raw][1]
    basis <- "sample n_raw from 03_qc_report"
    if (is.na(nraw)) { nraw <- n / ret[dataset == ds, retention]; basis <- "estimated from retention (n_raw missing)" }
  } else {
    nraw <- n / ret[dataset == ds, retention]
    basis <- "assembled cells / dataset median retention"
  }
  rate_raw <- nraw * RATE_PER_CELL
  rate <- min(0.4, rate_raw * same_donor)

  scr <- tryCatch({ sce <- Seurat::as.SingleCellExperiment(obj)
                    sce <- scDblFinder::scDblFinder(sce, dbr = rate, dbr.sd = 0)   # dbr.sd=0: dbr is the rate
                    list(score = setNames(sce$scDblFinder.score, colnames(sce)),
                         class = setNames(as.character(sce$scDblFinder.class) == "doublet", colnames(sce))) },
                  error = function(e) { message("  sc fail ", k, ": ", conditionMessage(e)); NULL })
  sc   <- if (is.null(scr)) NULL else scr$class
  sc_s <- if (is.null(scr)) NULL else scr$score
  df <- tryCatch(call_doubletfinder(obj, rate), error = function(e) { message("  df fail ", k, ": ", conditionMessage(e)); NULL })
  if (is.null(sc) && is.null(df)) { skipped <- c(skipped, paste0(k, "(both callers failed)")); rm(obj); gc(FALSE); next }

  cn <- colnames(obj)
  PC[[k]] <- data.table(dataset = ds, library = L$library, cell = cn,
                        sample = as.character(obj@meta.data$Sample),
                        sc_score = if (is.null(sc_s)) NA_real_ else as.numeric(sc_s[cn]),
                        sc_class = if (is.null(sc))   NA        else as.logical(sc[cn]),
                        df_class = if (is.null(df))   NA        else as.logical(df[cn]))
  LB[[k]] <- data.table(dataset = ds, library = L$library, assay_class = class(obj[["RNA"]])[1],
                        n_donors = length(shares), n_cells = n,
                        n_raw_used = round(nraw), same_donor_frac = round(same_donor, 4),
                        rate_before_demux = round(rate_raw, 5), rate_used = round(rate, 5),
                        rate_basis = basis,
                        rate_sc = if (is.null(sc)) NA_real_ else mean(sc),
                        rate_df = if (is.null(df)) NA_real_ else mean(df))
  message(sprintf("  %-34s n=%5d donors=%d rate=%.4f | sc=%.3f df=%.3f", k, n, length(shares), rate,
                  if (is.null(sc)) NA_real_ else mean(sc), if (is.null(df)) NA_real_ else mean(df)))
  rm(obj); gc(FALSE)
}
P <- rbindlist(PC, fill = TRUE); B <- rbindlist(LB, fill = TRUE)
fwrite(P, file.path(OUT_DIR, "ws_doublet_scores_percell.csv.gz"))
fwrite(B, file.path(OUT_DIR, "ws_doublet_libraries.csv"))

message(sprintf("\n[3] scored %d cells in %d libraries | skipped %d: %s",
                nrow(P), nrow(B), length(skipped), paste(skipped, collapse = ", ")))
message("\n[4] observed / expected by caller, per dataset")
B[, `:=`(ratio_sc = rate_sc / rate_used, ratio_df = rate_df / rate_used)]
print(B[, .(libs = .N, med_rate_exp = round(median(rate_used), 4),
            med_ratio_sc = round(median(ratio_sc, na.rm = TRUE), 2),
            med_ratio_df = round(median(ratio_df, na.rm = TRUE), 2)), by = dataset][order(dataset)])

## SELF-CHECKS ----------------------------------------------------------------------------------
message("\n[SELF-CHECK]")
message(sprintf("  (a) every scored cell appears exactly once -> %s",
                if (!anyDuplicated(P$cell)) "PASS" else "** FAIL **"))
message(sprintf("  (b) pooled libraries have same_donor_frac < 1 and single-donor libraries == 1 -> %s",
                if (all(B[n_donors == 1L, same_donor_frac] == 1) &&
                    all(B[n_donors > 1L, same_donor_frac] < 1)) "PASS" else "** FAIL **"))
message(sprintf("  (c) Jaccard(sc, df) per library, median %.2f (the 28-sample calibration got 0.18)",
                median(P[, .(j = if (sum(sc_class | df_class, na.rm = TRUE) > 0)
                              sum(sc_class & df_class, na.rm = TRUE) / sum(sc_class | df_class, na.rm = TRUE)
                            else NA_real_), by = library]$j, na.rm = TRUE)))
message(sprintf("  (d) GSE185381 cells recovered: %d scored vs %d in the eligible objects",
                nrow(P[dataset == "GSE185381"]), sum(B[dataset == "GSE185381", n_cells])))
message("[done] ", OUT_DIR)
