# utils.R ----
# Shared helpers used across stages. Source AFTER config_paths.R.
# Единственная реализация core16 (was duplicated 5x with 2 implementations) lives here.

suppressPackageStartupMessages({ library(data.table) })

## -- default assay name (10x RNA) ----
if (!exists("ASSAY")) ASSAY <- "RNA"

## -- logging: unified message("[N] ...") with optional timestamp ----
# message_ts() is the timestamped logger (replaces the legacy .log() from 02_preprocess).
# Plain message("[N] ...") per CODING_STANDARDS is also fine for step markers.
message_ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                            paste0(..., collapse = "")))
msg <- message_ts  # alias

## -- counts extraction, Seurat v5 layer / v4 slot tolerant ----
get_counts <- function(obj, assay = ASSAY) {
  ct <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
                 error = function(e) NULL)
  if (is.null(ct) || nrow(ct) == 0L) {
    ct <- SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
  }
  ct
}

## -- safe fwrite: always create the parent dir first (FAST is purgeable) ----
fwrite_safe <- function(x, file, ...) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file, ...)
}

## -- core16: extract the 16bp 10x cell-barcode core, tool-agnostic ----
# Different tools decorate barcodes differently:
#   inferCNV : "<sample>_<BC>-1"     Numbat : "<BC>"     STARsolo : "<BC>-1"
# We strip a leading "<sample>_" and a trailing "-<lane>", then VALIDATE the result is a
# 16bp ACGT string. A length/content check guards against silently mis-joining when a
# barcode format is unexpected (previously two different regex implementations disagreed).
core16 <- function(x, strict = TRUE) {
  core <- sub("-\\d+$", "", sub("^.*_", "", x))
  if (strict) {
    ok <- nchar(core) == 16L & grepl("^[ACGT]{16}$", core)
    if (!all(ok)) {
      n_bad <- sum(!ok)
      warning(sprintf("core16: %d/%d barcodes are not clean 16bp ACGT after stripping (e.g. '%s'). ",
                      n_bad, length(core), core[which(!ok)[1]]),
              "Check the barcode format for this tool/dataset before joining.")
    }
  }
  core
}

## -- list per-sample QC RDS for a dataset -> data.table(dataset, sample, rds) ----
list_qc_samples <- function(dataset) {
  d <- file.path(QC_RDS_DIR, dataset)
  if (!dir.exists(d)) stop("[utils] QC RDS dir not found: ", d)
  f <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  data.table(dataset = dataset, sample = sub("\\.rds$", "", basename(f)), rds = f)
}

## -- THE COHORT ROSTER IS THE QC REPORT, NOT THE DIRECTORY LISTING ----
# Eight scripts built their sample list with list.files(QC_RDS_DIR), which quietly
# equates "an object is on disk" with "this sample is in the study". The two
# diverged the moment the v2 ingest stopped emitting samples the v1 ingest had
# written: 77 objects from 2026-07-15 stayed behind, and J3 processed 74 of them,
# so ref_norm_summary described a cohort that no longer exists. Everything keyed
# off that summary -- including the malignancy calibration measurement, which is
# supposed to be the gate -- would have been computed over the wrong cohort and
# would have looked entirely normal doing it.
#
# This returns the PASS roster and REFUSES to guess when disk and roster disagree,
# because the failure is silent in both directions: extra objects inflate the
# cohort, missing ones shrink it, and neither raises an error on its own.
qc_rds_roster <- function(datasets = NULL, on_extra = c("error", "warn", "ignore")) {
  on_extra <- match.arg(on_extra)
  rep_csv <- file.path(DIR_PREPROCESS, "03_qc_report__ALL.csv")
  if (!file.exists(rep_csv))
    stop("[utils] QC report not found: ", rep_csv,
         "\n  the roster comes from the report, so run the J2c combiner first.")
  Q <- data.table::fread(rep_csv)
  scol <- intersect(c("Sample", "sample", "sample_id"), names(Q))[1]
  if (is.na(scol) || !all(c("dataset", "status") %in% names(Q)))
    stop("[utils] QC report lacks dataset/status/sample columns: ", rep_csv)
  R <- Q[status == "PASS", .(dataset, sample = as.character(get(scol)))]
  if (!is.null(datasets)) R <- R[dataset %in% datasets]
  R[, rds := file.path(QC_RDS_DIR, dataset, paste0(sample, ".rds"))]

  miss <- R[!file.exists(rds)]
  if (nrow(miss))
    stop(sprintf("[utils] %d PASS sample(s) have no QC object -- the QC stage is incomplete: %s",
                 nrow(miss), paste(head(miss$sample, 10), collapse = ", ")))

  # Only meaningful for a whole-cohort call; a per-dataset call legitimately sees
  # the other datasets' objects.
  if (is.null(datasets)) {
    on_disk <- list.files(QC_RDS_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
    on_disk <- on_disk[!grepl("/_", on_disk)]                 # scratch dirs (e.g. /_dryrun/)
    extra   <- setdiff(on_disk, R$rds)
    if (length(extra)) {
      msg <- sprintf(paste0(
        "[utils] %d QC object(s) under QC_RDS_DIR are NOT in the PASS roster -- stale output ",
        "from an earlier ingest. Quarantine them (99_admin/14_quarantine_stale_qc_objects.R) ",
        "before running any cohort-level stage.\n  e.g. %s"),
        length(extra), paste(head(basename(extra), 3), collapse = ", "))
      if (on_extra == "error") stop(msg) else if (on_extra == "warn") warning(msg, call. = FALSE)
    }
  }
  R[]
}

## -- list merged-RDS dataset basenames present in the ingest output ----
list_datasets <- function() {
  f <- list.files(RDS_INGEST_DIR, pattern = "\\.rds$", full.names = FALSE)
  sub("\\.rds$", "", f)
}
