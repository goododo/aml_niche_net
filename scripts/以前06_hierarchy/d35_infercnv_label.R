#!/usr/bin/env Rscript
# d35_infercnv_label.R ----
# Phase 2: derive a per-cell malignant label as the UNION of two CNV methods (inferCNV + Numbat),
# because each method has an orthogonal blind spot:
#   - inferCNV misses CN-normal / weak-CNV AML (e.g. 3853_Dg, 6323_Dg -> ~0% while Numbat ~50-66%)
#   - Numbat misses samples where clone reconstruction fails (e.g. 1216_Dg -> "no_CNV_detected"
#     while inferCNV finds 56.8%)
# A cell is malignant if EITHER valid method calls it. Where both methods are valid AND agree, the
# call is high-confidence (evidence_tier A); where only one method is valid (the other was blind or
# not run), it is single-evidence (tier C). Thresholds of each method are UNCHANGED (no loosening:
# relaxing inferCNV's reference quantile cannot conjure CNVs in CN-normal samples and only adds
# false positives in T/NK -- the fix is method complementarity, not a looser cutoff).
#
# Numbat per-cell has TWO schemas (schema-adaptive read):
#   full     : ...,compartment_opt,sample   -> malignant = compartment_opt == "tumor"
#   degraded : cell,sample,p_cnv,note        -> Numbat found nothing -> this method NA for the sample
#
# BARCODE MATCH: inferCNV cell = "<sample>_<16bp>-1"; Numbat cell = "<16bp>". Join on 16bp core.
#
# OUTPUT (in d30's consumed schema so d30 runs UNCHANGED via --malignancy_dir):
#   <MALIGNANCY_INFERCNV_DIR>/<ds>/<sample>__consensus_percell.csv   cell, malignant, confidence
#   <MALIGNANCY_INFERCNV_DIR>/<ds>/<sample>__consensus_summary.csv   sample, evidence_tier, ...
# Per-sample, serial, resume-skip.
# Usage:
#   Rscript d35_infercnv_label.R --dataset=GSE227903
#   Rscript d35_infercnv_label.R --dataset=GSE239721 --sample=PT10A
#   Rscript d35_infercnv_label.R --dataset=GSE227903 --no_numbat   # inferCNV-only (e.g. GSE239721)

suppressPackageStartupMessages({ library(optparse); library(data.table) })
source("d00_config.R")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", help = "dataset id"),
  make_option("--sample",  type = "character", default = NULL, help = "optional single sample"),
  make_option("--infercnv_dir", type = "character", default = NULL,
              help = "override inferCNV dir (default: INFERCNV_PERCELL_DIR/<ds>)"),
  make_option("--numbat_dir", type = "character", default = NULL,
              help = "override Numbat dir (default: NUMBAT_PERCELL_DIR/<ds>)"),
  make_option("--no_numbat", action = "store_true", default = FALSE,
              help = "ignore Numbat entirely (inferCNV-only; e.g. datasets without Numbat runs)")
)))
if (is.null(opt$dataset)) stop("--dataset is required")

in_dir   <- if (!is.null(opt$infercnv_dir)) opt$infercnv_dir else file.path(INFERCNV_PERCELL_DIR, opt$dataset)
nb_dir   <- if (!is.null(opt$numbat_dir))   opt$numbat_dir   else file.path(NUMBAT_PERCELL_DIR, opt$dataset)
out_dir  <- file.path(MALIGNANCY_INFERCNV_DIR, opt$dataset)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -- helpers ----
# 16bp barcode core: strip sample prefix (up to & incl last '_') and trailing -N.
core16 <- function(x) sub("-\\d+$", "", sub("^.*_", "", x))

# Read Numbat per-cell, schema-adaptive. Returns data.table(core, nb_malignant in {0,1}) or NULL
# if Numbat is unusable for this sample (degraded schema / missing file / no compartment_opt).
read_numbat <- function(samp) {
  nf <- file.path(nb_dir, samp, "numbat", paste0(samp, "__numbat_percell.csv"))
  if (!file.exists(nf)) return(NULL)
  d <- fread(nf)
  if (!"compartment_opt" %in% names(d)) return(NULL)   # degraded schema -> Numbat blind here
  if (!"cell" %in% names(d)) return(NULL)
  d[, core := core16(cell)]
  d[, nb_malignant := as.integer(compartment_opt == "tumor")]
  unique(d[, .(core, nb_malignant)], by = "core")
}

## -- discover inferCNV per-cell files ----
pc_files <- list.files(in_dir, pattern = "__infercnv_percell\\.csv$",
                       recursive = TRUE, full.names = TRUE)
if (!is.null(opt$sample)) {
  pc_files <- pc_files[basename(pc_files) == paste0(opt$sample, "__infercnv_percell.csv")]
}
if (!length(pc_files)) stop("[d35] no __infercnv_percell.csv under ", in_dir, " (run c41 first)")
message("[1] ", length(pc_files), " inferCNV per-cell file(s) in ", opt$dataset)

process_one <- function(pc_csv) {
  samp    <- sub("__infercnv_percell\\.csv$", "", basename(pc_csv))
  out_pc  <- file.path(out_dir, paste0(samp, "__consensus_percell.csv"))
  out_sm  <- file.path(out_dir, paste0(samp, "__consensus_summary.csv"))
  if (file.exists(out_pc) && file.exists(out_sm)) { message("    [skip] ", samp); return(invisible()) }

  d <- fread(pc_csv)  # cell, malignant, score, method, sample
  if (!all(c("cell", "malignant") %in% names(d))) {
    warning("[d35] ", samp, " missing cell/malignant cols; skipping"); return(invisible())
  }
  d[, ic_malignant := as.integer(malignant)]
  d[, core := core16(cell)]

  ## -- attach Numbat (if usable for this sample) ----
  nb <- if (opt$no_numbat) NULL else read_numbat(samp)
  numbat_usable <- !is.null(nb)
  if (numbat_usable) {
    d <- merge(d, nb, by = "core", all.x = TRUE)   # cells absent in Numbat -> nb_malignant NA
  } else {
    d[, nb_malignant := NA_integer_]
  }

  ## -- UNION call: malignant if EITHER valid method calls it ----
  # treat each method's NA as "no opinion" (not a normal vote), so union is blind-spot robust.
  d[, malignant := as.integer((ic_malignant %in% 1L) | (nb_malignant %in% 1L))]

  ## -- per-cell confidence from method agreement ----
  # both methods present a call for the cell:
  #   agree malignant (1,1) -> high ; disagree (1,0)/(0,1) -> medium ; agree normal (0,0) -> NA(normal)
  # only one method present -> low (single evidence)
  d[, confidence := fcase(
    malignant == 0L,                                          NA_character_,
    !is.na(ic_malignant) & !is.na(nb_malignant) &
      ic_malignant == 1L & nb_malignant == 1L,               "high",
    !is.na(ic_malignant) & !is.na(nb_malignant),             "medium",  # one says 1, other 0
    default = "low"                                                       # only one method available
  )]

  pc <- d[, .(cell, malignant, confidence)]
  fwrite_safe(pc, out_pc)

  ## -- sample-level evidence tier (drives d30 TIER2CONF) ----
  n_called <- nrow(pc)
  n_mal    <- sum(pc$malignant == 1L, na.rm = TRUE)
  ic_n     <- sum(d$ic_malignant == 1L, na.rm = TRUE)
  nb_n     <- if (numbat_usable) sum(d$nb_malignant == 1L, na.rm = TRUE) else NA_integer_
  # concordance among malignant calls: fraction of union-malignant cells where both methods agree
  both_n   <- sum(d$ic_malignant == 1L & d$nb_malignant == 1L, na.rm = TRUE)
  if (!numbat_usable) {
    tier <- TIER_CNV_SINGLE; src <- "infercnv_only"
  } else {
    conc <- if (n_mal > 0) both_n / n_mal else 0
    tier <- if (conc >= 0.5) TIER_CNV_CONCORDANT else TIER_CNV_UNION
    src  <- "infercnv_union_numbat"
  }
  sm <- data.table(
    sample          = samp,
    Dataset         = opt$dataset,
    evidence_tier   = tier,
    label_source    = src,
    numbat_usable   = numbat_usable,
    n_called        = n_called,
    n_malignant     = n_mal,
    n_infercnv_mal  = ic_n,
    n_numbat_mal    = nb_n,
    n_both_mal      = both_n,
    malignant_frac  = if (n_called > 0) n_mal / n_called else NA_real_
  )
  fwrite_safe(sm, out_sm)
  message(sprintf("    [write] %-10s union %d/%d (%.1f%%)  [ic=%d nb=%s both=%d] tier=%s",
                  samp, n_mal, n_called, 100 * n_mal / max(1, n_called),
                  ic_n, ifelse(is.na(nb_n), "NA", nb_n), both_n, tier))
  invisible()
}

for (x in pc_files) process_one(x)
message("[done] d35 (inferCNV union Numbat) for ", opt$dataset,
        "  -> ", out_dir, "  (feed d30 with --malignancy_dir=", out_dir, ")")
