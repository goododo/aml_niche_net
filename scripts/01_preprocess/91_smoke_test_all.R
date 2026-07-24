# =============================================================================
suppressPackageStartupMessages({ library(here) })
# 91_smoke_test_all.R
# Pre-flight smoke test BEFORE the full sbatch array.
# For each of the 13 datasets: load -> (demux prefilter) -> split by Sample ->
# pick ONE sample (subsampled for speed) -> run the REAL stage functions from
# 03 (mad_keep, call_scdblfinder, call_doubletfinder) inside per-stage tryCatch.
#
# Purpose: confirm every dataset can be split, MAD runs, and BOTH doublet tools
# return without error. It does NOT save per-sample objects (smoke test only).
#
# Run:
#   Rscript 91_smoke_test_all.R                 # maxcells = 1200 (fast)
#   Rscript 91_smoke_test_all.R --maxcells 800
# Output:
#   <tab_dir>/91_smoke_test.csv
# =============================================================================

t0 <- Sys.time()
options(aml_qc.source_only = TRUE)
source(here::here("scripts", "01_preprocess", "03_per_sample_qc.R"))     # reuse real functions; guard prevents full run

suppressPackageStartupMessages({ library(Seurat); library(data.table) })

args     <- commandArgs(trailingOnly = TRUE)
k        <- which(args == "--maxcells")
maxcells <- if (length(k) == 1 && length(args) >= k + 1) as.integer(args[k + 1]) else 1200L

datasets <- list_datasets()
message_ts("smoke test over ", length(datasets), " datasets; maxcells = ", maxcells)

# Run a stage and capture (ok, message). Returns list(ok=, msg=).
try_stage <- function(expr) {
  res <- tryCatch(list(ok = TRUE, val = force(expr), msg = ""),
                  error = function(e) list(ok = FALSE, val = NULL,
                                           msg = conditionMessage(e)))
  res
}

rows <- list()

for (ds in datasets) {
  td <- Sys.time()
  r  <- data.table(dataset = ds, sample = NA_character_,
                   n_raw = NA_integer_, n_after_mad = NA_integer_,
                   split_ok = FALSE, mad_ok = FALSE,
                   scdbl_ok = NA, df_ok = NA,
                   n_sc_dbl = NA_integer_, n_df_dbl = NA_integer_,
                   note = "", secs = NA_real_)

  ld <- try_stage({
    seu <- readRDS(file.path(RDS_INGEST_DIR, paste0(ds, ".rds")))
    ensure_joined(seu)
  })
  if (!ld$ok) { r[, note := paste0("load FAIL: ", ld$msg)]; rows[[ds]] <- r; next }
  seu <- ld$val

  # demux prefilter (no-op unless a Demuxed column exists, e.g. Lasry)
  if (isTRUE(DEMUX_PREFILTER) && COL_DEMUXED %in% colnames(seu@meta.data)) {
    dmx <- seu@meta.data[[COL_DEMUXED]]
    keep_dmx <- dmx %in% c(TRUE, "TRUE", "True", "true", 1, "1")
    seu <- seu[, keep_dmx]
    r[, note := paste0("demux ", sum(keep_dmx), "/", length(keep_dmx), "; ")]
  }

  # integrity heads-up (e.g. counts-loading failure)
  med_nc <- suppressWarnings(median(seu@meta.data[[COL_NCOUNT]], na.rm = TRUE))
  if (is.finite(med_nc) && med_nc < INTEGRITY_MIN_MED_NCOUNT) {
    r[, note := paste0(note, "INTEGRITY: median nCount=", round(med_nc, 1), "; ")]
  }

  # split + pick the largest sample, then subsample for speed
  sp <- try_stage({
    cnts <- sort(table(as.character(seu@meta.data[[COL_SAMPLE]])), decreasing = TRUE)
    samp <- names(cnts)[1]
    obj  <- SplitObject(seu, split.by = COL_SAMPLE)[[samp]]
    list(samp = samp, obj = obj)
  })
  rm(seu); gc(verbose = FALSE)
  if (!sp$ok) { r[, `:=`(note = paste0(note, "split FAIL: ", sp$msg))]
                rows[[ds]] <- r; next }
  r[, `:=`(split_ok = TRUE, sample = sp$val$samp, n_raw = ncol(sp$val$obj))]
  obj <- sp$val$obj
  if (!is.na(maxcells) && ncol(obj) > maxcells) {
    set.seed(SEED); obj <- obj[, sample(colnames(obj), maxcells)]
  }

  # MAD
  md <- try_stage({
    keep <- mad_keep(obj@meta.data)
    obj2 <- obj[, keep]; obj2
  })
  if (!md$ok) { r[, note := paste0(note, "MAD FAIL: ", md$msg)]; rows[[ds]] <- r; next }
  obj <- md$val
  r[, `:=`(mad_ok = TRUE, n_after_mad = ncol(obj))]

  # doublet (skip if too few cells; that is not a failure)
  if (ncol(obj) < DOUBLET_MIN_CELLS) {
    r[, note := paste0(note, "doublet skipped (n=", ncol(obj), " < ",
                       DOUBLET_MIN_CELLS, ")")]
  } else {
    set.seed(SEED)
    sc <- try_stage(call_scdblfinder(obj))
    df <- try_stage(call_doubletfinder(obj))
    r[, `:=`(scdbl_ok = sc$ok, df_ok = df$ok,
             n_sc_dbl = if (sc$ok) sum(sc$val, na.rm = TRUE) else NA_integer_,
             n_df_dbl = if (df$ok) sum(df$val, na.rm = TRUE) else NA_integer_)]
    if (!sc$ok) r[, note := paste0(note, "scDblFinder ERR: ", sc$msg, "; ")]
    if (!df$ok) r[, note := paste0(note, "DoubletFinder ERR: ", df$msg, "; ")]
  }

  r[, secs := round(as.numeric(difftime(Sys.time(), td, units = "secs")), 1)]
  message_ts(sprintf("%-14s split=%s mad=%s sc=%s df=%s (%ss)",
               ds, r$split_ok, r$mad_ok, r$scdbl_ok, r$df_ok, r$secs))
  rows[[ds]] <- r
}

smoke <- rbindlist(rows, use.names = TRUE, fill = TRUE)
fwrite_safe(smoke, file.path(DIR_PREPROCESS, "91_smoke_test.csv"))

# overall verdict: every dataset must split + MAD; doublet must not ERROR where
# it was attempted (NA = skipped-lowN is acceptable).
dbl_failed <- smoke[(scdbl_ok == FALSE) | (df_ok == FALSE)]
all_ok <- all(smoke$split_ok) && all(smoke$mad_ok) && nrow(dbl_failed) == 0

cat("\n================ SMOKE TEST SUMMARY ================\n")
print(as.data.frame(smoke[, .(dataset, sample, n_raw, n_after_mad,
                              split_ok, mad_ok, scdbl_ok, df_ok, secs, note)]))
cat("\nTotal elapsed: ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min\n")
cat(ifelse(all_ok,
           ">>> SMOKE TEST PASSED -- safe to launch the full array\n",
           ">>> SMOKE TEST FAILED -- inspect 'note' / FAIL rows before sbatch\n"))
