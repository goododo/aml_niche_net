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


## ---------------------------------------------------------------------------------------------
## curated per-sample timepoint / patient -- the ONE source
## ---------------------------------------------------------------------------------------------
# Replaces three private tp_from_name() copies (02_per_bin_malignant.R, 04_stemness_score.R,
# 04_cnmf/04_cnmf_figures.R) that parsed the timepoint out of the SAMPLE NAME. Between them they
# resolved 27 of 214 samples, all from one dataset, and used a vocabulary ("Dx"/"MRD"/"Relapse")
# that matched neither CANONICAL_TIMEPOINTS nor the Python stages. The patient id was derived the
# same way -- by stripping a suffix -- when the manifest carries uid_patient directly.
#
# Returns dataset, sample, timepoint (canonical), tp_axis (3-level), patient, uid_patient.
sample_timepoints <- function(manifest_csv = NULL) {
  if (is.null(manifest_csv))
    manifest_csv <- file.path(FAST_DIR, "results/tables/01_preprocess/01_sample_role_manifest.csv")
  if (!file.exists(manifest_csv)) stop("no curated manifest at ", manifest_csv)
  M <- data.table::fread(manifest_csv)
  need <- c("dataset", "Sample", "Timepoint")
  miss <- setdiff(need, names(M))
  if (length(miss)) stop("manifest is missing column(s): ", paste(miss, collapse = ", "))
  out <- M[, .(dataset, sample = Sample, timepoint = as.character(Timepoint),
               patient    = if ("Patient_ID"  %in% names(M)) as.character(Patient_ID)  else NA_character_,
               uid_patient= if ("uid_patient" %in% names(M)) as.character(uid_patient) else NA_character_)]
  out[, tp_axis := tp_axis(timepoint)]
  unique(out, by = c("dataset", "sample"))
}

# Join the curated timepoint onto a table and REFUSE to proceed on poor coverage. A silent NA here
# is what turned a 28-patient longitudinal cohort into a 6-patient one.
add_timepoint <- function(DT, min_cov = 0.95, manifest_csv = NULL) {
  stopifnot(all(c("dataset", "sample") %in% names(DT)))
  TP <- sample_timepoints(manifest_csv)
  for (cc in c("timepoint", "tp_axis", "patient", "uid_patient"))
    if (cc %in% names(DT)) DT[, (cc) := NULL]
  out <- merge(DT, TP, by = c("dataset", "sample"), all.x = TRUE)
  cov <- mean(!is.na(out$timepoint))
  if (cov < min_cov)
    stop(sprintf(paste0("curated timepoint covers only %.1f%% of rows (need >= %.0f%%). ",
                        "Unmatched: %s. A name-derived fallback is exactly the bug this replaces."),
                 100 * cov, 100 * min_cov,
                 paste(utils::head(unique(out[is.na(timepoint)]$sample), 5), collapse = ", ")))
  message(sprintf("[timepoint] curated labels joined for %d of %d rows (%.1f%%); axis: %s",
                  sum(!is.na(out$timepoint)), nrow(out), 100 * cov,
                  paste(sprintf("%s=%d", names(table(out$tp_axis)), table(out$tp_axis)), collapse = " ")))
  out[]
}


## ---------------------------------------------------------------------------------------------
## resume guards must test FRESHNESS, not existence
## ---------------------------------------------------------------------------------------------
# Generalises the idiom already used at 02_malignancy/44_infercnv_run_one.R:62-73 and
# 03_hierarchy/01_bmm_project.R:114, which every stage after 03 was missing. Existence-only guards
# pinned published tables to superseded inputs: 68 stemness files dated 2026-07-17 whose three
# inputs all postdate them (QC objects 08-05, projection 08-07, consensus 08-13/14/17) and which
# reference 48,087 cells that no longer exist; 33 cellstate files written before their consensus
# existed and therefore carrying malignant=NA permanently; and the whole 05_ccc -> 06 -> 07 -> 08
# chain, where results/tables/07_fgw/patient_scores.csv holds 148 rows of which 55 name samples
# that left the cohort. Re-running the chain printed "[skip]" five times and exited 0.
#
# Returns TRUE when `out` must be rebuilt. Missing inputs are IGNORED rather than treated as fresh:
# an input that is not there yet cannot certify an output as current.
is_stale <- function(out, inputs, force = FALSE) {
  if (isTRUE(force)) return(TRUE)
  out <- out[nzchar(out)]
  if (!length(out) || !all(file.exists(out))) return(TRUE)
  inputs <- unique(inputs[nzchar(inputs)])
  inputs <- inputs[file.exists(inputs)]
  if (!length(inputs)) return(FALSE)
  min(file.mtime(out)) < max(file.mtime(inputs))
}

# Print WHICH input made it stale. "[recompute]" with no reason is how a stale-output problem
# becomes invisible again the next time someone reads the log.
stale_reason <- function(out, inputs, force = FALSE) {
  # `force` must be reported as force. Without it the log reads "RECOMPUTE -- up to date", which is
  # a contradiction, and a reader cannot tell a forced rebuild from a genuine staleness detection.
  if (isTRUE(force)) return("forced (--force)")
  out <- out[nzchar(out) & file.exists(out)]
  inputs <- unique(inputs[nzchar(inputs)]); inputs <- inputs[file.exists(inputs)]
  if (!length(out) || !length(inputs)) return("output missing")
  newer <- inputs[file.mtime(inputs) > min(file.mtime(out))]
  if (!length(newer)) return("up to date")
  sprintf("%d newer input(s), e.g. %s (%s) > output (%s)", length(newer), basename(newer[1]),
          format(file.mtime(newer[1]), "%m-%d %H:%M"), format(min(file.mtime(out)), "%m-%d %H:%M"))
}
