# =============================================================================
# 04_export_role_split_xlsx.R
# Reads the REAL outputs of Task 1 / Task 2 and writes a colour-coded workbook:
#   Sheet 1 "Study_roles_split"  : study-level role + 70/30 split (one row/study)
#   Sheet 2 "Sample_manifest"    : sample-level role + split (one row/sample)
#   Sheet 3 "Legend"             : colour key + column definitions
#
# Inputs (produced by 01_dataset_roles.R and 02_study_split.R):
#   <tab_dir>/01_study_role_summary.csv
#   <tab_dir>/02_study_split.csv
#   <tab_dir>/02_sample_split.csv   (manifest + split, sample level)
# Output:
#   <tab_dir>/AML_Phase1_roles_and_split.xlsx
# Run:
#   Rscript 04_export_role_split_xlsx.R     (after 01 and 02 have run)
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

# ----------------------------------------------------------------------------
# Lane palette (shared by both sheets). Lane = operational train/val category.
# ----------------------------------------------------------------------------
pal <- c(
  Discovery          = "#C6EFCE",  # green
  Validation         = "#BDD7EE",  # blue
  `Treatment-axis`   = "#C2E0D6",  # teal
  Healthy            = "#D9D9D9",  # grey
  Reference          = "#E1D5E7",  # light purple
  Exploratory        = "#FFF2CC",  # yellow
  Exclude            = "#F4CCCC",  # red
  Other              = "#FFFFFF"
)
flag_fill   <- "#FCE4D6"           # peach: integrity / attention flag
header_fill <- "#1F4E78"           # dark blue header

# Map any role/split string to a lane key for colouring.
to_lane <- function(x) {
  x <- as.character(x)
  fcase(
    grepl("Exclude", x, ignore.case = TRUE),       "Exclude",
    grepl("Healthy", x, ignore.case = TRUE),       "Healthy",
    grepl("Reference", x, ignore.case = TRUE),     "Reference",
    grepl("Treatment-axis", x, ignore.case = TRUE),"Treatment-axis",
    grepl("Exploratory", x, ignore.case = TRUE),   "Exploratory",
    grepl("Discovery", x, ignore.case = TRUE),     "Discovery",
    grepl("Validation", x, ignore.case = TRUE),    "Validation",
    default = "Other")
}

read_if <- function(p) if (file.exists(p)) fread(p) else {
  warning("missing: ", p); NULL
}

study_role  <- read_if(file.path(DIR_PREPROCESS, "01_study_role_summary.csv"))
study_split <- read_if(file.path(DIR_PREPROCESS, "02_study_split.csv"))
sample_tbl  <- read_if(file.path(DIR_PREPROCESS, "02_sample_split.csv"))

wb <- createWorkbook()
hdr_style <- createStyle(fgFill = header_fill, fontColour = "white",
                         textDecoration = "bold", halign = "center",
                         border = "TopBottom")

# Generic: write a data.frame to a sheet, colour rows by `lane_col`, freeze header.
write_coloured <- function(wb, sheet, df, lane_col, flag_col = NULL, widths = NULL) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df, headerStyle = hdr_style)
  n <- nrow(df)
  if (n > 0) {
    lanes <- to_lane(df[[lane_col]])
    for (ln in unique(lanes)) {
      rws <- which(lanes == ln) + 1L                 # +1 for header
      st  <- createStyle(fgFill = pal[[ln]])
      addStyle(wb, sheet, st, rows = rws, cols = seq_len(ncol(df)),
               gridExpand = TRUE, stack = TRUE)
    }
    if (!is.null(flag_col) && flag_col %in% names(df)) {
      fl <- which(!is.na(df[[flag_col]]) & df[[flag_col]] != "OK") + 1L
      if (length(fl)) {
        st <- createStyle(fgFill = flag_fill, textDecoration = "bold")
        addStyle(wb, sheet, st, rows = fl, cols = seq_len(ncol(df)),
                 gridExpand = TRUE, stack = TRUE)
      }
    }
  }
  freezePane(wb, sheet, firstActiveRow = 2, firstActiveCol = 2)
  if (!is.null(widths)) setColWidths(wb, sheet, cols = seq_along(widths), widths = widths)
  else setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = "auto")
}

# ----- Sheet 1: study level -----
if (!is.null(study_split)) {
  s1 <- copy(study_split)
  if (!is.null(study_role)) {
    add <- setdiff(c("n_samples","n_discovery_pool","n_healthy","n_excluded","integrity"),
                   names(s1))
    s1 <- merge(s1, study_role[, c("dataset", intersect(add, names(study_role))), with = FALSE],
                by = "dataset", all.x = TRUE)
  }
  # colour by a finer lane: treatment-axis distinguished from plain Validation
  s1[, lane_key := fifelse(grepl("Treatment-axis", study_role), "Treatment-axis", split)]
  setcolorder(s1, c("dataset","study_role","split","lane_key"))
  write_coloured(wb, "Study_roles_split", as.data.frame(s1),
                 lane_col = "lane_key", flag_col = "integrity")
}

# ----- Sheet 2: sample level -----
if (!is.null(sample_tbl)) {
  s2 <- copy(sample_tbl)
  lane_src <- if ("split_sample" %in% names(s2)) "split_sample" else "sample_role"
  s2[, lane_key := to_lane(get(lane_src))]
  write_coloured(wb, "Sample_manifest", as.data.frame(s2),
                 lane_col = "lane_key", flag_col = "integrity_flag")
}

# ----- Sheet 3: legend -----
legend <- data.table(
  Lane = names(pal),
  Meaning = c("Discovery (random 70% pool)","Validation (held-out 30%)",
              "Treatment-axis (longitudinal, forced to val side)",
              "Healthy control (B_healthy)","Reference scaffold / stroma aux",
              "Exploratory (hypothesis only)","Excluded (cell line / PDX / fail)",
              "Unclassified"))
addWorksheet(wb, "Legend")
writeData(wb, "Legend", legend, headerStyle = hdr_style)
for (i in seq_len(nrow(legend))) {
  addStyle(wb, "Legend", createStyle(fgFill = pal[[legend$Lane[i]]]),
           rows = i + 1L, cols = 1:2, gridExpand = TRUE, stack = TRUE)
}
addStyle(wb, "Legend", createStyle(fgFill = flag_fill, textDecoration = "bold"),
         rows = nrow(legend) + 3L, cols = 1:2, gridExpand = TRUE)
writeData(wb, "Legend", data.frame(Note = "Peach + bold row = integrity/attention flag (e.g. re-import, sorted-only)"),
          startRow = nrow(legend) + 3L, colNames = FALSE)
setColWidths(wb, "Legend", cols = 1:2, widths = c(18, 55))

out <- file.path(DIR_PREPROCESS, "AML_Phase1_roles_and_split.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)
message_ts("wrote workbook: ", out)
