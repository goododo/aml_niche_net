#!/usr/bin/env Rscript
# 05_qc_landing_check.R ----
# Cross-check per-sample QC Seurat objects (directories on /LARGE1) against the
# per-dataset QC report CSVs (on /FAST). Decide which datasets landed cleanly,
# whether outputs are mutually complete (dir <-> csv), per-sample count consistency,
# and flag datasets that need a re-run.
#
# WHERE TO RUN : an HPC login/compute node where /LARGE1 and /FAST are mounted.
# DEPENDS ON   : R + data.table only (no Seurat load required; we just list .rds files).
# DESIGN NOTE  : column names of the QC report CSVs are unknown a priori, so this
#                script is introspective: it prints the detected schema and uses
#                configurable candidate column names with graceful degradation.

suppressPackageStartupMessages({ library(here); library(data.table) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))

## ---- 0. Paths (from config layer) ----
QC_OBJ_DIR  <- QC_RDS_DIR        # per-sample QC objects on LARGE1 (01_preprocess output)
QC_TAB_DIR  <- DIR_PREPROCESS    # per-sample QC report CSVs
OUT_TAB_DIR <- DIR_PREPROCESS    # verification summary lands beside QC tables

# Report file naming: 03_qc_report__<DATASET>.csv
REPORT_PREFIX <- "03_qc_report__"
REPORT_REGEX  <- paste0("^", REPORT_PREFIX, "(.+)\\.csv$")

# Directory names to ignore inside QC_OBJ_DIR (scaffolding / dry runs / hidden).
IGNORE_DIRS   <- c("_dryrun")

# Per-sample object file extension (case-insensitive).
RDS_REGEX     <- "\\.rds$"

# Sample-level cell gate from the blueprint (>= 500 total cells per sample).
MIN_CELLS     <- 500L

# Expected active dataset roster (13). Used only to surface datasets that are
# absent from BOTH the object dir and the report tables (i.e. never produced).
EXPECTED_DATASETS <- c(
  "Chen2023", "E-MTAB-11536", "GSE116256", "GSE147989", "GSE185381",
  "GSE185991", "GSE201966", "GSE207356", "GSE227903", "GSE239721",
  "GSE253355", "GSE289435", "Petti2019"
)

# Candidate column names for introspection. The first match (case-insensitive)
# in each group is used; if none match, that metric is reported as NA and we
# fall back to row counts only.
SAMPLE_COL_CANDIDATES <- c("Sample", "sample", "Sample_ID", "sample_id",
                           "Patient_ID", "orig.ident", "GSM", "library", "lib")
NCELL_COL_CANDIDATES  <- c("n_cells_pass", "n_cells_post", "nCells_post",
                           "n_cells", "nCells", "cells_pass", "n_pass",
                           "n_cells_after", "ncell_post")
FLAG_COL_CANDIDATES   <- c("integrity_flag", "integrity", "qc_flag",
                           "status", "flag")

## ---- 1. Helper functions ----

# Pick the first candidate present in `cols` (case-insensitive); return the
# actual column name as it appears in the data, or NA_character_.
pick_col <- function(cols, candidates) {
  lc <- tolower(cols)
  for (cand in candidates) {
    hit <- which(lc == tolower(cand))
    if (length(hit)) return(cols[hit[1]])
  }
  NA_character_
}

# Robust CSV read; returns NULL on failure.
safe_fread <- function(path) {
  out <- tryCatch(fread(path, showProgress = FALSE),
                  error = function(e) NULL)
  out
}

# Count .rds files (non-recursive, files only) in a directory.
count_rds <- function(dir) {
  if (!dir.exists(dir)) return(character(0))
  ff <- list.files(dir, pattern = RDS_REGEX, ignore.case = TRUE,
                   full.names = FALSE, recursive = FALSE)
  ff
}

## ---- 2. Discover dataset dirs and report CSVs ----
message("[1] Scanning object dirs and report CSVs ...")

# 2a. Dataset directories on /LARGE1 (exclude ignore list + hidden).
obj_dirs_all <- list.dirs(QC_OBJ_DIR, full.names = FALSE, recursive = FALSE)
obj_dirs     <- setdiff(obj_dirs_all, IGNORE_DIRS)
obj_dirs     <- obj_dirs[!grepl("^\\.", obj_dirs)]

# 2b. Report CSVs on /FAST -> extract dataset name from filename.
csv_files <- list.files(QC_TAB_DIR, pattern = REPORT_REGEX, full.names = FALSE)
csv_dsets <- sub(REPORT_REGEX, "\\1", csv_files)

# 2c. Union of all dataset names seen anywhere (plus the expected roster, so a
#     dataset missing from BOTH sources still shows up as a row).
all_dsets <- sort(unique(c(obj_dirs, csv_dsets, EXPECTED_DATASETS)))

message(sprintf("    object dirs found : %d  (ignored: %s)",
                length(obj_dirs),
                paste(intersect(obj_dirs_all, IGNORE_DIRS), collapse = ", ")))
message(sprintf("    report CSVs found : %d", length(csv_files)))
message(sprintf("    datasets to check : %d", length(all_dsets)))

## ---- 3. Per-dataset cross-check ----
message("[2] Cross-checking each dataset (dir <-> csv, counts, flags) ...")

rows        <- vector("list", length(all_dsets))
schema_log  <- list()  # store detected colnames per dataset for the user to verify

for (i in seq_along(all_dsets)) {
  ds       <- all_dsets[i]
  this_dir <- file.path(QC_OBJ_DIR, ds)
  this_csv <- file.path(QC_TAB_DIR, paste0(REPORT_PREFIX, ds, ".csv"))

  dir_exists <- dir.exists(this_dir)
  csv_exists <- file.exists(this_csv)

  # --- object side ---
  rds_files <- if (dir_exists) count_rds(this_dir) else character(0)
  n_rds     <- length(rds_files)

  # --- report side ---
  n_csv_rows    <- NA_integer_
  n_csv_samples <- NA_integer_
  sample_col    <- NA_character_
  ncell_col     <- NA_character_
  flag_col      <- NA_character_
  n_pass_gate   <- NA_integer_
  n_flag_notok  <- NA_integer_
  csv_cols      <- character(0)

  if (csv_exists) {
    dt <- safe_fread(this_csv)
    if (!is.null(dt) && nrow(dt) > 0L) {
      csv_cols   <- names(dt)
      n_csv_rows <- nrow(dt)

      sample_col <- pick_col(csv_cols, SAMPLE_COL_CANDIDATES)
      ncell_col  <- pick_col(csv_cols, NCELL_COL_CANDIDATES)
      flag_col   <- pick_col(csv_cols, FLAG_COL_CANDIDATES)

      if (!is.na(sample_col)) {
        n_csv_samples <- length(unique(dt[[sample_col]]))
      } else {
        n_csv_samples <- n_csv_rows  # assume 1 row per sample if no id column
      }

      if (!is.na(ncell_col)) {
        v <- suppressWarnings(as.numeric(dt[[ncell_col]]))
        n_pass_gate <- sum(!is.na(v) & v >= MIN_CELLS)
      }

      if (!is.na(flag_col)) {
        fv <- toupper(trimws(as.character(dt[[flag_col]])))
        n_flag_notok <- sum(!(fv %in% c("OK", "PASS", "TRUE", "1")), na.rm = TRUE)
      }
    } else {
      csv_cols <- "<empty-or-unreadable>"
    }
  }

  schema_log[[ds]] <- csv_cols

  # --- count consistency (only meaningful when both sides present) ---
  count_consistent <- NA
  if (dir_exists && csv_exists && !is.na(n_csv_samples)) {
    count_consistent <- (n_rds == n_csv_samples)
  }

  # --- status decision ---
  status <- if (!dir_exists && !csv_exists) {
    "ABSENT_BOTH"          # expected but never produced
  } else if (dir_exists && !csv_exists) {
    "MISSING_CSV"          # objects landed, summary report missing
  } else if (!dir_exists && csv_exists) {
    "MISSING_DIR"          # report exists, objects missing
  } else if (n_rds == 0L) {
    "EMPTY_DIR"            # dir exists but no .rds inside
  } else if (isFALSE(count_consistent)) {
    "COUNT_MISMATCH"       # #rds != #samples in report
  } else {
    "OK"
  }

  # --- secondary review flag (integrity not OK) ---
  review <- if (!is.na(n_flag_notok) && n_flag_notok > 0L) "REVIEW_FLAG" else ""

  rows[[i]] <- data.table(
    dataset          = ds,
    dir_exists       = dir_exists,
    csv_exists       = csv_exists,
    n_rds_files      = n_rds,
    n_csv_rows       = n_csv_rows,
    n_csv_samples    = n_csv_samples,
    count_consistent = count_consistent,
    sample_col       = sample_col,
    ncell_col        = ncell_col,
    flag_col         = flag_col,
    n_pass_min_cells = n_pass_gate,
    n_flag_not_ok    = n_flag_notok,
    status           = status,
    review           = review
  )
}

summary_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
setorder(summary_dt, status, dataset)

## ---- 4. Console report ----
message("[3] Summary (status per dataset):")
print(summary_dt[, .(dataset, dir_exists, csv_exists,
                     n_rds_files, n_csv_samples, count_consistent,
                     status, review)])

# 4a. Detected CSV schemas (so the user can correct the candidate column lists).
message("\n[4] Detected CSV column names per dataset (verify SAMPLE/NCELL/FLAG picks):")
for (ds in names(schema_log)) {
  cols <- schema_log[[ds]]
  if (length(cols) && !identical(cols, character(0))) {
    message(sprintf("    %-14s : %s", ds, paste(cols, collapse = ", ")))
  }
}

# 4b. Cross-direction mismatches.
dir_no_csv <- summary_dt[status == "MISSING_CSV", dataset]
csv_no_dir <- summary_dt[status == "MISSING_DIR", dataset]
absent     <- summary_dt[status == "ABSENT_BOTH", dataset]
empty      <- summary_dt[status == "EMPTY_DIR",  dataset]
mismatch   <- summary_dt[status == "COUNT_MISMATCH", dataset]
review_ds  <- summary_dt[review == "REVIEW_FLAG", dataset]

message("\n[5] Action items:")
message(sprintf("    OK (clean)            : %s",
                paste(summary_dt[status == "OK", dataset], collapse = ", ")))
message(sprintf("    NEEDS RERUN (objects) : %s", paste(c(csv_no_dir, empty), collapse = ", ")))
message(sprintf("    NEEDS RERUN (report)  : %s", paste(dir_no_csv, collapse = ", ")))
message(sprintf("    COUNT MISMATCH        : %s", paste(mismatch, collapse = ", ")))
message(sprintf("    ABSENT FROM BOTH      : %s", paste(absent, collapse = ", ")))
message(sprintf("    INTEGRITY REVIEW      : %s", paste(review_ds, collapse = ", ")))

# 4c. For mismatches, print the two lists side by side so the user can eyeball
#     which sample is missing on which side.
if (length(mismatch)) {
  message("\n[6] Count-mismatch detail (rds files vs report rows):")
  for (ds in mismatch) {
    this_dir <- file.path(QC_OBJ_DIR, ds)
    rds_files <- count_rds(this_dir)
    message(sprintf("    --- %s ---", ds))
    message(sprintf("      rds files (%d): %s", length(rds_files),
                    paste(rds_files, collapse = ", ")))
  }
}

## ---- 5. Persist summary ----
if (!dir.exists(OUT_TAB_DIR)) dir.create(OUT_TAB_DIR, recursive = TRUE)
out_path <- file.path(OUT_TAB_DIR, "04_qc_landing_check__summary.csv")
fwrite_safe(summary_dt, out_path)
message(sprintf("\n[7] Written: %s", out_path))

message("[done] QC landing check complete.")
