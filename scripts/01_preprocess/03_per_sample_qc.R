# =============================================================================
# 03_per_sample_qc.R   (Task 3)
# Per-sample QC on the merged RDS objects:
#   1. Split each merged object by Sample.
#   2. Per-sample data-driven MAD outlier filtering (scater::isOutlier) on
#      log10(nCount_RNA), log10(nFeature_RNA) [both tails] and percent.mt [upper].
#   3. Doublet removal via dual-R consensus (scDblFinder + DoubletFinder),
#      intersection by default.
#   4. Sample-level gate: keep sample only if >= MIN_CELLS_SAMPLE cells
#      remain after QC + doublet removal; otherwise -> P2-exclude.
#   5. Save per-sample QC'd objects + a per-sample QC report.
#
# DUAL RUN MODE (chosen workflow):
#   (a) Check locally / serially, ONE dataset at a time:
#         Rscript 03_per_sample_qc.R --dataset Chen2023
#       or all datasets serially:
#         Rscript 03_per_sample_qc.R --all
#   (b) SLURM array (one array task per dataset) -- see submit_03_per_sample_qc.sbatch:
#         the script reads SLURM_ARRAY_TASK_ID and picks the i-th dataset.
#
# Within a dataset, samples are processed serially by default; set
# cfg_parallel_workers > 1 to parallelise samples with future (single node).
#
# Outputs:
#   <obj_out_dir>/<dataset>/<Sample>.rds            (PASS samples only)
#   <tab_dir>/03_qc_report__<dataset>.csv           (one row per sample)
#   <tab_dir>/03_P2_exclude__<dataset>.csv          (failed samples)
#   <fig_dir>/03_qc_violin__<dataset>.png           (pre/post QC overview)
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(scater)         # isOutlier
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
})

cfg_parallel_workers <- 1L   # set >1 to parallelise SAMPLES within a dataset
options(future.globals.maxSize = 32000 * 1024^2)

# ----------------------------------------------------------------------------
# Argument / dataset selection (supports --dataset / --all / SLURM array)
# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
all_datasets <- list_datasets()

pick_datasets <- function() {
  # SLURM array task id takes priority if present.
  aid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  if (nzchar(aid)) {
    i <- as.integer(aid)
    stopifnot(i >= 1, i <= length(all_datasets))
    return(all_datasets[i])
  }
  if ("--all" %in% args) return(all_datasets)
  k <- which(args == "--dataset")
  if (length(k) == 1 && length(args) >= k + 1) return(args[k + 1])
  stop("Specify --dataset <name>, --all, or run as a SLURM array job.")
}
# NOTE: dataset selection is deferred to the run block below so that sourcing
# this file (e.g. from the dry-run script) only DEFINES functions.

# ----------------------------------------------------------------------------
# Helper: ensure a counts layer is joined (Seurat v5 may store split layers).
# ----------------------------------------------------------------------------
ensure_joined <- function(obj) {
  if (inherits(obj[["RNA"]], "Assay5")) {
    obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  }
  obj
}

# ----------------------------------------------------------------------------
# Helper: per-sample MAD outlier flags using existing per-cell metrics.
# Returns a logical vector "keep" (TRUE = passes QC).
# ----------------------------------------------------------------------------
mad_keep <- function(md) {
  nc <- md[[COL_NCOUNT]]; nf <- md[[COL_NFEAT]]; mt <- md[[COL_MT]]
  out_nc <- isOutlier(nc, nmads = MAD_NMADS, log = MAD_LOG_COUNTS, type = "both")
  out_nf <- isOutlier(nf, nmads = MAD_NMADS, log = MAD_LOG_FEATS,  type = "both")
  out_mt <- isOutlier(mt, nmads = MAD_NMADS, log = FALSE, type = MT_TAIL)
  out_hb <- rep(FALSE, length(nc)); hb_hi <- NA_real_
  if (isTRUE(MAD_HB_ENABLE) && COL_HB %in% names(md)) {
    hb <- md[[COL_HB]]
    # NOTE: hb MAD is a RELATIVE upper-tail flag, so genuine erythroid-program
    # cells are not blanket-removed the way a hard hb cap would do.
    out_hb <- isOutlier(hb, nmads = MAD_NMADS, log = FALSE, type = MAD_HB_TAIL)
    hb_hi  <- attr(out_hb, "thresholds")[["higher"]]
  }
  # ABSOLUTE floors/caps, applied IN ADDITION to the per-sample MAD. isOutlier alone is purely
  # RELATIVE, so the worse a sample is the looser its own thresholds become (observed nfeat_lo
  # down to 0.08 genes, percent.mt cutoffs computed above 100%) -- filter strength ended up
  # inversely proportional to sample quality. A cell must clear BOTH gates.
  # NA-safe: a missing metric must not silently propagate into the keep vector.
  na0 <- function(x, fill) { x[is.na(x)] <- fill; x }
  abs_fail <- (na0(nf, Inf) < ABS_MIN_NFEAT) |
              (na0(nc, Inf) < ABS_MIN_NCOUNT) |
              (na0(mt, -Inf) > ABS_MAX_MT)

  out_mad <- out_nc | out_nf | out_mt | out_hb
  keep <- !(out_mad | abs_fail)
  attr(keep, "thresh") <- list(
    # EFFECTIVE thresholds (MAD combined with the absolute bound), so the report shows the
    # bound that actually applied rather than the relative one that was overridden.
    ncount_lo = max(attr(out_nc, "thresholds")[["lower"]],  ABS_MIN_NCOUNT),
    ncount_hi = attr(out_nc, "thresholds")[["higher"]],
    nfeat_lo  = max(attr(out_nf, "thresholds")[["lower"]],  ABS_MIN_NFEAT),
    nfeat_hi  = attr(out_nf, "thresholds")[["higher"]],
    mt_hi     = min(attr(out_mt, "thresholds")[["higher"]], ABS_MAX_MT),
    hb_hi     = hb_hi,
    # provenance: how much of the filtering each gate is responsible for (QC-covariate table)
    ncount_lo_mad = attr(out_nc, "thresholds")[["lower"]],
    nfeat_lo_mad  = attr(out_nf, "thresholds")[["lower"]],
    mt_hi_mad     = attr(out_mt, "thresholds")[["higher"]],
    n_fail_mad    = sum(out_mad, na.rm = TRUE),
    n_fail_abs    = sum(abs_fail & !out_mad, na.rm = TRUE))
  keep
}

# ----------------------------------------------------------------------------
# Helper: scDblFinder doublet calls on a (QC-passed) Seurat object.
# ----------------------------------------------------------------------------
call_scdblfinder <- function(obj, rate) {
  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(sce, dbr = rate)
  cls <- as.character(sce$scDblFinder.class)
  setNames(cls == "doublet", colnames(sce))
}

# ----------------------------------------------------------------------------
# Helper: DoubletFinder doublet calls. Standard pipeline + pK sweep + homotypic
# adjustment. Handles both old (_v3 suffix) and new DoubletFinder APIs.
# ----------------------------------------------------------------------------
call_doubletfinder <- function(obj, rate) {
  suppressPackageStartupMessages(library(DoubletFinder))
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)

  has_v3 <- exists("paramSweep_v3")
  ps  <- if (has_v3) paramSweep_v3(obj, PCs = 1:30, sct = FALSE)
         else        paramSweep(obj, PCs = 1:30, sct = FALSE)
  sw  <- summarizeSweep(ps, GT = FALSE)
  bcmvn <- find.pK(sw)
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  homo    <- modelHomotypic(Idents(obj))
  nExp    <- round(rate * ncol(obj))
  nExp.adj<- max(0, round(nExp * (1 - homo)))

  obj <- if (has_v3)
           doubletFinder_v3(obj, PCs = 1:30, pN = 0.25, pK = pK,
                            nExp = nExp.adj, sct = FALSE)
         else
           doubletFinder(obj, PCs = 1:30, pN = 0.25, pK = pK,
                         nExp = nExp.adj, sct = FALSE)
  cls_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)[1]
  cls <- obj@meta.data[[cls_col]]
  setNames(cls == "Doublet", colnames(obj))
}

# ----------------------------------------------------------------------------
# Per-sample QC worker. Returns a one-row report data.table and (on PASS)
# writes the cleaned object to disk.
# ----------------------------------------------------------------------------
qc_one_sample <- function(obj, ds, samp, integrity_fail) {
  obj <- ensure_joined(obj)
  md  <- obj@meta.data
  n0  <- ncol(obj)

  rep <- data.table(dataset = ds, Sample = samp, n_raw = n0,
                    n_after_mad = NA_integer_, n_doublet = NA_integer_,
                    n_final = NA_integer_, status = NA_character_,
                    reason = NA_character_,
                    ncount_lo = NA_real_, ncount_hi = NA_real_,
                    nfeat_lo = NA_real_, nfeat_hi = NA_real_,
                    mt_hi = NA_real_, hb_hi = NA_real_,
                    med_ncount_final = NA_real_, med_nfeat_final = NA_real_,
                    lowcomplexity_flag = FALSE, dbl_method = NA_character_,
                    # -- sample-level QC-STRENGTH covariates (Phase 8 confounding check) --
                    # Per-sample filter strength varies 10-100x across this cohort, and the whole
                    # project compares networks BETWEEN samples: if "what counts as a cell" differs
                    # by 10x, that difference propagates straight into the edge weights. These are
                    # carried forward so HDS/ATS/RLS/SCS can be regressed on them; a score that
                    # tracks mad_loss or nfeat_lo is a QC artefact, not biology.
                    ncount_lo_mad = NA_real_, nfeat_lo_mad = NA_real_, mt_hi_mad = NA_real_,
                    n_fail_mad = NA_integer_, n_fail_abs = NA_integer_,
                    mad_loss = NA_real_, dbl_rate_exp = NA_real_, dbl_rate_obs = NA_real_,
                    platform = NA_character_, upstream_filtered = NA)

  # --- integrity short-circuit (e.g. GSE227903 counts not loaded) ---
  if (isTRUE(integrity_fail)) {
    rep[, `:=`(status = "P2-exclude", reason = "INTEGRITY_FAIL: counts not loaded")]
    return(rep)
  }

  # --- 1. MAD cell filtering ---
  keep <- mad_keep(md); th <- attr(keep, "thresh")
  obj  <- obj[, keep]
  rep[, `:=`(n_after_mad = ncol(obj),
             ncount_lo = th$ncount_lo, ncount_hi = th$ncount_hi,
             nfeat_lo = th$nfeat_lo, nfeat_hi = th$nfeat_hi,
             mt_hi = th$mt_hi, hb_hi = th$hb_hi,
             ncount_lo_mad = th$ncount_lo_mad, nfeat_lo_mad = th$nfeat_lo_mad,
             mt_hi_mad = th$mt_hi_mad,
             n_fail_mad = th$n_fail_mad, n_fail_abs = th$n_fail_abs,
             mad_loss = if (n0 > 0) 1 - ncol(obj) / n0 else NA_real_)]

  # --- 2. doublet consensus (skip if too few cells) ---
  # Expected rate is driven by how many barcodes were LOADED, so it is computed on n_raw, not on
  # the MAD-survivors: MAD removes ~24% of cells on average, which systematically deflated the
  # expected rate (and therefore the callers' dbr prior) in the previous run.
  # Platform matters here: the 0.8%/1k law is a 10x droplet model and is wrong for Seq-Well.
  plat     <- platform_of(ds)$platform
  rate_exp <- expected_dbl_rate(n0, plat)
  rep[, `:=`(dbl_rate_exp = rate_exp, platform = plat,
             upstream_filtered = is_upstream_filtered(ds))]

  if (ncol(obj) < DOUBLET_MIN_CELLS) {
    rep[, `:=`(n_doublet = 0L, dbl_method = "skipped_lowN")]
  } else {
    sc <- tryCatch(call_scdblfinder(obj, rate_exp),   error = function(e) { message_ts("scDblFinder fail ", samp, ": ", conditionMessage(e)); NULL })
    df <- tryCatch(call_doubletfinder(obj, rate_exp), error = function(e) { message_ts("DoubletFinder fail ", samp, ": ", conditionMessage(e)); NULL })
    if (is.null(sc) && is.null(df)) {
      is_dbl <- rep(FALSE, ncol(obj)); rep[, dbl_method := "both_failed"]
    } else if (is.null(df)) {
      is_dbl <- sc[colnames(obj)]; rep[, dbl_method := "scDblFinder_only"]
    } else if (is.null(sc)) {
      is_dbl <- df[colnames(obj)]; rep[, dbl_method := "DoubletFinder_only"]
    } else {
      sc <- sc[colnames(obj)]; df <- df[colnames(obj)]
      is_dbl <- if (DOUBLET_CONSENSUS == "union") (sc | df) else (sc & df)
      rep[, dbl_method := paste0("consensus_", DOUBLET_CONSENSUS)]
    }
    rep[, n_doublet := sum(is_dbl, na.rm = TRUE)]
    rep[, dbl_rate_obs := n_doublet / n_after_mad]
    obj <- obj[, !is_dbl]
  }

  # --- 3. sample-level gate ---
  rep[, n_final := ncol(obj)]

  # low-complexity soft flag: catches uniformly poor libraries that per-sample
  # MAD cannot (e.g. GSE201966 AML01M4_CR/_P). FLAG by default; drop only if
  # FLAG_LOWCOMPLEXITY_DROP is TRUE.
  if (ncol(obj) > 0) {
    mednc <- median(obj@meta.data[[COL_NCOUNT]], na.rm = TRUE)
    mednf <- median(obj@meta.data[[COL_NFEAT]],  na.rm = TRUE)
    rep[, `:=`(med_ncount_final = mednc, med_nfeat_final = mednf,
               lowcomplexity_flag = (mednc < FLAG_MIN_MED_NCOUNT |
                                     mednf < FLAG_MIN_MED_NFEAT))]
  }
  drop_lowc <- isTRUE(FLAG_LOWCOMPLEXITY_DROP) && isTRUE(rep$lowcomplexity_flag)

  if (ncol(obj) >= MIN_CELLS_SAMPLE && !drop_lowc) {
    rep[, `:=`(status = "PASS",
               reason = fifelse(lowcomplexity_flag, "PASS_lowcomplexity_flag", "OK"))]
    out_d <- file.path(QC_RDS_DIR, ds)
    if (!dir.exists(out_d)) dir.create(out_d, recursive = TRUE, showWarnings = FALSE)
    saveRDS(obj, file.path(out_d, paste0(samp, ".rds")))
  } else {
    rsn <- if (drop_lowc) "LOW_COMPLEXITY (median below floor)"
           else sprintf("n_final=%d < %d", ncol(obj), MIN_CELLS_SAMPLE)
    rep[, `:=`(status = "P2-exclude", reason = rsn)]
  }
  rep
}

# ----------------------------------------------------------------------------
# Dataset driver.
# ----------------------------------------------------------------------------
process_dataset <- function(ds) {
  message_ts("==== dataset: ", ds, " ====")
  integrity <- ROLE_TABLE[dataset == ds, integrity_flag]
  integrity_fail <- length(integrity) == 1 && integrity == "INTEGRITY_FAIL"

  seu <- readRDS(file.path(RDS_INGEST_DIR, paste0(ds, ".rds")))
  seu <- ensure_joined(seu)

  # Demux pre-filter (Lasry GSE185381): keep only successfully demultiplexed
  # cells. Drops the un-demuxed composite/fallback libraries before splitting,
  # so they never become per-sample objects and fail the cell-count gate naturally.
  #
  # GUARD: fire ONLY where demultiplexing was actually ATTEMPTED, i.e. at least one
  # cell carries a donor call. "Demuxed == FALSE everywhere" does not mean every cell
  # failed demux -- it means the deposit ships no per-cell donor call at all, which is
  # true of GSE185991 (3 two-patient libraries, curated demux_method = "unknown"). The
  # unguarded version read that as "0/159183 cells kept" and killed the whole dataset.
  # The distinction matters: an un-demuxed cell in a demuxed library belongs to an
  # unknown one of 5 donors and is junk; a cell in a never-demuxed library still belongs
  # to a known library whose donor composition is recorded in N_donors_in_library.
  if (isTRUE(DEMUX_PREFILTER) && COL_DEMUXED %in% colnames(seu@meta.data)) {
    dmx      <- seu@meta.data[[COL_DEMUXED]]
    keep_dmx <- dmx %in% c(TRUE, "TRUE", "True", "true", 1, "1")
    if (!any(keep_dmx)) {
      message_ts(ds, ": Demuxed pre-filter SKIPPED -- no cell carries a donor call, so ",
                 "demultiplexing was never attempted for this deposit. ",
                 "Multi-donor libraries stay flagged via Patient_resolved instead.")
    } else {
      message_ts(ds, ": Demuxed pre-filter ", sum(keep_dmx), "/", length(keep_dmx), " cells kept")
      seu <- seu[, keep_dmx]
    }
  }

  # Dataset-level integrity floor (auto-detect broken counts even if not flagged).
  med_nc <- median(seu@meta.data[[COL_NCOUNT]], na.rm = TRUE)
  if (med_nc < INTEGRITY_MIN_MED_NCOUNT) {
    message_ts("INTEGRITY: ", ds, " median nCount = ", round(med_nc, 1),
         " -> flagging all samples as P2-exclude (re-import required)")
    integrity_fail <- TRUE
  }

  samp_levels <- unique(as.character(seu@meta.data[[COL_SAMPLE]]))
  obj_list <- SplitObject(seu, split.by = COL_SAMPLE)
  rm(seu); gc(verbose = FALSE)

  run_one <- function(samp) {
    set.seed(SEED)   # reproducible doublet sampling per sample
    qc_one_sample(obj_list[[samp]], ds, samp, integrity_fail)
  }

  if (cfg_parallel_workers > 1L) {
    suppressPackageStartupMessages(library(future.apply))
    future::plan(future::multicore, workers = cfg_parallel_workers)
    reports <- future.apply::future_lapply(samp_levels, run_one,
                                            future.seed = TRUE)
    future::plan(future::sequential)
  } else {
    reports <- lapply(samp_levels, run_one)
  }
  report <- rbindlist(reports, use.names = TRUE, fill = TRUE)

  # write tables
  fwrite_safe(report, file.path(DIR_PREPROCESS, paste0("03_qc_report__", ds, ".csv")))
  excl <- report[status == "P2-exclude"]
  fwrite_safe(excl, file.path(DIR_PREPROCESS, paste0("03_P2_exclude__", ds, ".csv")))

  # quick overview figure (raw counts distribution per sample)
  p <- ggplot(report, aes(x = Sample)) +
    geom_col(aes(y = n_final), fill = "#3b7dd8") +
    geom_hline(yintercept = MIN_CELLS_SAMPLE, linetype = 2, colour = "red") +
    coord_flip() +
    labs(title = paste0("Post-QC cell counts: ", ds),
         y = "cells after QC + doublet removal", x = NULL) +
    theme_bw(base_size = 9)
  ggsave(file.path(FIG_PREPROCESS, paste0("03_qc_cellcounts__", ds, ".png")),
         p, width = 6, height = max(3, 0.18 * nrow(report)), dpi = 150)

  message_ts(ds, ": PASS=", sum(report$status == "PASS"),
       " / P2-exclude=", sum(report$status == "P2-exclude"))
  report
}

# ----------------------------------------------------------------------------
# Run  (guarded: skipped when sourced with options(aml_qc.source_only = TRUE),
#       e.g. by the dry-run script, which reuses the functions above)
# ----------------------------------------------------------------------------
if (!isTRUE(getOption("aml_qc.source_only", FALSE))) {
  datasets <- pick_datasets()
  all_reports <- rbindlist(lapply(datasets, process_dataset),
                           use.names = TRUE, fill = TRUE)

  # When run with --all (or after all array tasks finish, re-run this
  # aggregation), write a combined report.
  if ("--all" %in% args || length(datasets) > 1) {
    fwrite_safe(all_reports, file.path(DIR_PREPROCESS, "03_qc_report__ALL.csv"))
  }
  message_ts("Task 3 done for: ", paste(datasets, collapse = ", "))
}
