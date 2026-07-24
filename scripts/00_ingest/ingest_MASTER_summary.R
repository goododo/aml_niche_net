# ============================================================================
# ingest_MASTER_summary.R   --   combine all per-dataset QC summaries
# Reads every "<DATASET>_qc_summary.csv" in the 01_qc table dir, stacks them,
# and writes a master CSV + XLSX plus a per-dataset overview to the console.
# Run this LAST, after all 01_qc_<DATASET>.R scripts have finished.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

suppressPackageStartupMessages(library(writexl))

# ----------------------------------------------------------------------------
# [1] Collect per-dataset summary tables ----
# ----------------------------------------------------------------------------
message("[1] collecting per-dataset summaries")
summary_files <- list.files(DIR_INGEST, pattern = "_qc_summary\\.csv$", full.names = TRUE)
summary_files <- summary_files[basename(summary_files) != "00_MASTER_qc_summary.csv"]
message("    ", length(summary_files), " summary files found")
stopifnot(length(summary_files) > 0)

master <- dplyr::bind_rows(lapply(summary_files, function(f) {
  readr::read_csv(f, col_types = readr::cols(
    Dataset    = readr::col_character(),
    Sample     = readr::col_character(),
    Patient_ID = readr::col_character(),   # "1216" stays a string, not a double
    Timepoint  = readr::col_character(),
    .default   = readr::col_double()
  ))
}))

# ----------------------------------------------------------------------------
# [2] Per-dataset overview ----
# ----------------------------------------------------------------------------
message("[2] per-dataset overview")
overview <- master %>%
  dplyr::group_by(Dataset) %>%
  dplyr::summarise(
    n_samples   = dplyr::n(),
    n_patients  = dplyr::n_distinct(Patient_ID),
    total_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Dataset)

print(as.data.frame(overview))
message("    grand total: ", sum(overview$total_cells), " cells across ",
        sum(overview$n_samples), " samples / ", nrow(overview), " datasets")

# ----------------------------------------------------------------------------
# [3] Write master CSV + XLSX ----
# ----------------------------------------------------------------------------
message("[3] writing master tables")
fwrite_safe(master, file.path(DIR_INGEST, "00_MASTER_qc_summary.csv"))
writexl::write_xlsx(
  list(per_sample = as.data.frame(master), overview = as.data.frame(overview)),
  path = file.path(DIR_INGEST, "00_MASTER_qc_summary.xlsx")
)

message("[ok] master summary written to ", DIR_INGEST)
