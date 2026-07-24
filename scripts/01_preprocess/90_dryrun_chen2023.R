# =============================================================================
suppressPackageStartupMessages({ library(here) })
# 90_dryrun_chen2023.R
# 5-minute INTEGRATION test of the real QC chain on ONE Chen2023 sample.
# Reuses the actual functions from 03_per_sample_qc.R (no logic duplication):
#   Demux pre-filter (n/a for Chen) -> per-sample MAD -> scDblFinder +
#   DoubletFinder consensus -> >=500 gate -> save .rds.
#
# Writes to a *_dryrun output dir so real outputs are never touched.
#
# Run:
#   Rscript 90_dryrun_chen2023.R                  # auto-picks largest sample
#   Rscript 90_dryrun_chen2023.R --sample AML103_Niche_Immune
#   Rscript 90_dryrun_chen2023.R --maxcells 1500  # subsample for speed
# =============================================================================

t0 <- Sys.time()

# Source 03 in "functions only" mode (guard prevents the full run block).
options(aml_qc.source_only = TRUE)
source(here::here("scripts", "01_preprocess", "03_per_sample_qc.R"))   # this also sources the config layer and loads libs

suppressPackageStartupMessages({ library(Seurat); library(data.table) })

# ---- args ----
args     <- commandArgs(trailingOnly = TRUE)
get_arg  <- function(flag, default = NA) {
  k <- which(args == flag); if (length(k) == 1 && length(args) >= k + 1) args[k + 1] else default
}
ds        <- "Chen2023"
want_samp <- get_arg("--sample", NA)
maxcells  <- suppressWarnings(as.integer(get_arg("--maxcells", NA)))

stopifnot(file.exists(file.path(RDS_INGEST_DIR, paste0(ds, ".rds"))))

# ---- redirect outputs to a dry-run sandbox ----
orig_obj_dir   <- QC_RDS_DIR
QC_RDS_DIR <- file.path(orig_obj_dir, "_dryrun")
dir.create(QC_RDS_DIR, recursive = TRUE, showWarnings = FALSE)

message_ts("DRY-RUN on dataset: ", ds)
seu <- readRDS(file.path(RDS_INGEST_DIR, paste0(ds, ".rds")))
seu <- ensure_joined(seu)

# Demux pre-filter (no-op for Chen, but exercises the code path if present).
if (isTRUE(DEMUX_PREFILTER) && COL_DEMUXED %in% colnames(seu@meta.data)) {
  dmx <- seu@meta.data[[COL_DEMUXED]]
  keep_dmx <- dmx %in% c(TRUE, "TRUE", "True", "true", 1, "1")
  message_ts("Demuxed pre-filter ", sum(keep_dmx), "/", length(keep_dmx), " kept")
  seu <- seu[, keep_dmx]
}

# ---- pick the sample ----
samp_counts <- sort(table(as.character(seu@meta.data[[COL_SAMPLE]])), decreasing = TRUE)
if (is.na(want_samp)) {
  samp <- names(samp_counts)[1]                 # largest -> exercises save path
  message_ts("auto-picked largest sample: ", samp, " (", samp_counts[1], " cells)")
} else {
  stopifnot(want_samp %in% names(samp_counts))
  samp <- want_samp
}

obj_list <- SplitObject(seu, split.by = COL_SAMPLE)
obj <- obj_list[[samp]]
rm(seu, obj_list); gc(verbose = FALSE)

# optional subsample for a guaranteed-fast check
if (!is.na(maxcells) && ncol(obj) > maxcells) {
  set.seed(SEED)
  obj <- obj[, sample(colnames(obj), maxcells)]
  message_ts("subsampled to ", ncol(obj), " cells for speed")
}

# ---- run the REAL per-sample QC function ----
set.seed(SEED)
rep <- qc_one_sample(obj, ds, samp, integrity_fail = FALSE)

# ---- assertions / integration checks ----
ok <- TRUE
check <- function(cond, msg) { if (!isTRUE(cond)) { ok <<- FALSE; message_ts("FAIL: ", msg) }
                               else message_ts("pass: ", msg) }
check(rep$n_raw >= rep$n_after_mad,                 "n_after_mad <= n_raw")
check(rep$n_after_mad >= rep$n_final,               "n_final <= n_after_mad")
check(!is.na(rep$status),                           "status assigned")
check(grepl("consensus|skipped|only|failed",
            ifelse(is.na(rep$dbl_method), "", rep$dbl_method)),
      "doublet method recorded")
if (rep$status == "PASS") {
  f <- file.path(QC_RDS_DIR, ds, paste0(samp, ".rds"))
  check(file.exists(f), paste0("saved object exists: ", f))
  if (file.exists(f)) {
    o2 <- readRDS(f)
    check(ncol(o2) == rep$n_final, "saved object cell count == n_final")
  }
}

cat("\n================ DRY-RUN REPORT ================\n")
print(as.data.frame(rep))
cat("\nElapsed: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min\n")
cat(ifelse(ok, ">>> INTEGRATION TEST PASSED\n", ">>> INTEGRATION TEST FAILED (see FAIL lines)\n"))

# restore (in case of interactive reuse)
QC_RDS_DIR <- orig_obj_dir
