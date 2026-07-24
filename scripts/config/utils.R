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

## -- list merged-RDS dataset basenames present in the ingest output ----
list_datasets <- function() {
  f <- list.files(RDS_INGEST_DIR, pattern = "\\.rds$", full.names = FALSE)
  sub("\\.rds$", "", f)
}
