#!/usr/bin/env Rscript
# d30_per_bin_malignant.R ----
# Phase 2 step 3: overlay the (OPTIONAL) Phase-1b per-cell malignant consensus onto the binned
# backbone, cross high_error x malignant per cell, and emit the per-sample x per-bin summary with
# DECOMPOSED confidence (three components + a derived categorical -- no single opaque scalar).
#
# Malignant overlay is OPTIONAL: if the consensus csv is absent, malignant columns are NA and the
# backbone per-bin summary (n_cells, frac_high_error, bin membership) is still produced.
#
# Per-sample, serial, resume-skip.
# Usage:
#   Rscript d30_per_bin_malignant.R --dataset=GSE227903
#   Rscript d30_per_bin_malignant.R --dataset=GSE227903 --sample=6323_R

suppressPackageStartupMessages({ library(optparse); library(data.table) })
source("d00_config.R")

PERCELL_BINNED_DIR <- file.path(FAST_DIR, "results/tables/04_hierarchy/percell_binned")
PERCELL_LABELED_DIR <- file.path(FAST_DIR, "results/tables/04_hierarchy/percell_labeled")  # backbone + malignant overlay

# Knob (promote to d00 if you like): a bin whose malignant cells are mostly single-type evidence
# is capped at 'low' confidence even if the sample's evidence_tier is A.
MIN_FRAC_ROBUST <- 0.50   # frac of malignant cells with per-cell confidence == 'high' to avoid the cap

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", help = "dataset id"),
  make_option("--sample",  type = "character", default = NULL, help = "optional single sample"),
  make_option("--malignancy_dir", type = "character", default = NULL,
              help = "override dir holding <sample>__consensus_percell.csv (default: MALIGNANCY_DIR/<ds>)")
)))
if (is.null(opt$dataset)) stop("--dataset is required")
mal_dir <- if (!is.null(opt$malignancy_dir)) opt$malignancy_dir else file.path(MALIGNANCY_DIR, opt$dataset)

## -- Helpers ----
# 16bp barcode core: strip sample prefix (everything up to & incl. last '_') and trailing -N.
core16 <- function(x) sub("-\\d+$", "", sub("^.*_", "", x))

# Wilson score interval for a binomial proportion (no extra packages).
wilson_ci <- function(x, n, z = 1.96) {
  if (is.na(n) || n == 0) return(list(lo = NA_real_, hi = NA_real_))
  p <- x / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(lo = max(0, ctr - hw), hi = min(1, ctr + hw))
}

## -- Samples ----
bdir <- file.path(PERCELL_BINNED_DIR, opt$dataset)
f <- list.files(bdir, pattern = "_binned\\.csv$", full.names = TRUE)
if (!is.null(opt$sample)) f <- f[basename(f) == paste0(opt$sample, "_binned.csv")]
if (!length(f)) stop("[d30] no binned per-cell csv for ", opt$dataset, " (run d20 first)")
message("[1] ", length(f), " sample(s); malignancy dir = ", mal_dir)

process_one <- function(binned_csv) {
  samp     <- sub("_binned\\.csv$", "", basename(binned_csv))
  out_pb   <- file.path(PERBIN_DIR, opt$dataset, paste0(samp, "_perbin.csv"))
  out_pc   <- file.path(PERCELL_LABELED_DIR, opt$dataset, paste0(samp, "_labeled.csv"))
  if (file.exists(out_pb)) { message("    [skip] ", samp); return(invisible()) }

  dt <- fread(binned_csv)
  dt[, core := core16(cell)]

  ## -- optional malignant overlay ----
  cons_pc <- file.path(mal_dir, paste0(samp, "__consensus_percell.csv"))
  cons_sm <- file.path(mal_dir, paste0(samp, "__consensus_summary.csv"))
  have_mal <- file.exists(cons_pc)
  evidence_tier <- NA_character_

  if (have_mal) {
    cp <- fread(cons_pc, select = intersect(c("cell", "malignant", "confidence"),
                                            names(fread(cons_pc, nrows = 0))))
    cp[, core := core16(cell)]
    cov <- mean(dt$core %in% cp$core)
    message("    [overlay] ", samp, " consensus join coverage = ", round(100 * cov, 1), "%")
    dt <- merge(dt, cp[, .(core, malignant, mal_conf = confidence)], by = "core", all.x = TRUE)
    if (file.exists(cons_sm)) {
      sm <- fread(cons_sm)
      if ("evidence_tier" %in% names(sm)) evidence_tier <- as.character(sm$evidence_tier[1])
    }
  } else {
    message("    [overlay] ", samp, " NO consensus -> backbone only (malignant = NA)")
    dt[, `:=`(malignant = NA_integer_, mal_conf = NA_character_)]
  }

  ## -- per-cell aberrant class (cross high_error x malignant; never drops cells) ----
  dt[, aberrant_class := fcase(
    is.na(malignant),                                              NA_character_,
    malignant == 1 & mapping_qc_flag == "high_error",             "malignant_high_error",
    malignant == 1,                                               "malignant",
    malignant == 0 & mapping_qc_flag == "high_error",             "normal_high_error",
    default = "normal"
  )]
  fwrite_safe(dt, out_pc)

  ## -- per-bin summary with decomposed confidence ----
  tier_letter <- if (!is.na(evidence_tier)) substr(evidence_tier, 1, 1) else NA_character_
  pb <- dt[, {
    n   <- .N
    fhe <- mean(mapping_qc_flag == "high_error")
    if (have_mal) {
      called <- sum(!is.na(malignant))
      nmal   <- sum(malignant == 1, na.rm = TRUE)
      frac   <- if (called > 0) nmal / called else NA_real_
      ci     <- wilson_ci(nmal, called)
      nmal_he <- sum(malignant == 1 & mapping_qc_flag == "high_error", na.rm = TRUE)
      # robustness of the malignant calls in this bin: frac with per-cell confidence == 'high'
      mc     <- mal_conf[which(malignant == 1)]
      frob   <- if (length(mc)) mean(mc == "high", na.rm = TRUE) else NA_real_
    } else {
      called <- NA_integer_; nmal <- NA_integer_; frac <- NA_real_
      ci <- list(lo = NA_real_, hi = NA_real_); nmal_he <- NA_integer_; frob <- NA_real_
    }
    .(n_cells = n, n_called = called, n_malignant = nmal, malignant_frac = frac,
      malignant_frac_lo = ci$lo, malignant_frac_hi = ci$hi,
      n_malignant_high_error = nmal_he, frac_high_error = fhe,
      frac_robust = frob, in_ccc_graph = in_ccc_graph[1])
  }, by = hierarchy_bin]

  pb[, `:=`(Dataset = opt$dataset, Sample = samp, evidence_tier = evidence_tier)]

  ## -- derived categorical bin_confidence (conservative; components stay separate) ----
  tier_conf <- if (!is.na(tier_letter) && tier_letter %in% names(TIER2CONF)) TIER2CONF[[tier_letter]] else NA_character_
  pb[, bin_confidence := {
    out <- rep(NA_character_, .N)
    for (i in seq_len(.N)) {
      if (n_cells[i] < N_MIN_BIN || frac_high_error[i] > FRAC_ERR_MAX) { out[i] <- "low"; next }
      if (!have_mal) { out[i] <- NA_character_; next }           # no malignant labels -> no fraction confidence
      base <- if (!is.na(tier_conf)) tier_conf else "medium"
      # frac_robust penalty applies ONLY to single-evidence samples (tier C). For tier A/B the
      # malignant set legitimately mixes high(both-methods-agree) + medium(single-method) cells --
      # that's the normal shape of a union call, so frac_robust must NOT veto the sample-level tier.
      # frac_robust stays in the table as an independent column (it flags how many "hard" malignant
      # cells a bin holds; useful for picking high-purity cNMF input in Phase 3).
      is_single_evidence <- is.na(tier_letter) || identical(tier_letter, "C")
      if (is_single_evidence && !is.na(frac_robust[i]) && frac_robust[i] < MIN_FRAC_ROBUST)
        base <- "low"
      out[i] <- base
    }
    out
  }]

  setcolorder(pb, c("Dataset", "Sample", "hierarchy_bin", "in_ccc_graph", "n_cells", "n_called",
                    "n_malignant", "malignant_frac", "malignant_frac_lo", "malignant_frac_hi",
                    "n_malignant_high_error", "frac_high_error", "frac_robust",
                    "evidence_tier", "bin_confidence"))
  setorder(pb, hierarchy_bin)
  fwrite_safe(pb, out_pb)
  message("    [write] ", out_pb, "  (", nrow(pb), " bins)")
  invisible()
}

for (x in f) process_one(x)
message("[done] d30 for ", opt$dataset)