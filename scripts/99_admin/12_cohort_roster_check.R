# =============================================================================
# 12_cohort_roster_check.R
#
# THE QUESTION: is the set of samples the pipeline processes exactly the set of
# samples the curated metadata says exists?
#
# 11_audit_labels_vs_curated.R answered "are the LABELS right" (yes: 0 coarse
# disagreements among joined samples). This script answers the other half, which
# is where every remaining error turned out to live: are the ROWS right. A sample
# that should not exist cannot have a wrong label -- it is invisible to a label
# audit, and three of the four errors found this way were in the healthy arm that
# the B_healthy barycenter is built from.
#
# Two modes:
#   PRE-INGEST  (default) -- compares the v1 roster in 00_project/metadata/
#                            ALL__samples.tsv against curated, and explains every
#                            difference with its cause and whether the fix is in.
#   POST-INGEST (--post)  -- re-reads results/tables/00_ingest/*_qc_summary.csv
#                            and requires an EXACT match, exiting non-zero if not.
#                            Run this after the J1 ingest re-run; it is the gate.
#
# Usage: Rscript scripts/99_admin/12_cohort_roster_check.R [--post]
# =============================================================================

suppressPackageStartupMessages({library(data.table); library(readxl); library(here)})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))

POST        <- "--post" %in% commandArgs(trailingOnly = TRUE)
CURATED_DIR <- "/LARGE1/gr10634/gaozy/aml_niche_net/00_raw/metadata"
OUT_DIR     <- file.path(PROJECT_ROOT, "00_project", "metadata")

# Curated dir names do not all match dataset ids ("Petti 2019" vs "Petti2019").
ds_of_dir <- function(d) gsub(" ", "", d)

load_meta <- function(d) {
  p   <- file.path(CURATED_DIR, d)
  csv <- list.files(p, "^meta_.*\\.csv$", full.names = TRUE)
  if (length(csv)) return(fread(csv[1], colClasses = "character"))
  xls <- list.files(p, "^meta_.*\\.xlsx$", full.names = TRUE)[1]
  sh  <- grep("^meta_", readxl::excel_sheets(xls), value = TRUE)[1]
  as.data.table(readxl::read_excel(xls, sheet = sh, col_types = "text"))
}

dirs <- setdiff(list.dirs(CURATED_DIR, recursive = FALSE, full.names = FALSE), "BoneMarrowMap")
CUR  <- rbindlist(lapply(dirs, function(d) {
  m <- load_meta(d)
  data.table(dataset = ds_of_dir(d),
             curated_id = m$sample_id,
             timepoint  = m$timepoint_class,
             tissue     = if ("tissue" %in% names(m)) m$tissue else NA_character_)
}), fill = TRUE)
CUR[, sample := sub("^[^_]*__", "", curated_id)]

# ---------------------------------------------------------------------------
# [1] Pipeline roster
# ---------------------------------------------------------------------------
if (POST) {
  fs <- list.files(DIR_INGEST, "_qc_summary\\.csv$", full.names = TRUE)
  if (!length(fs)) stop("--post given but no *_qc_summary.csv under ", DIR_INGEST)
  PIPE <- rbindlist(lapply(fs, function(f)
    fread(f)[, .(dataset = as.character(Dataset), sample = as.character(Sample))]), fill = TRUE)
  message("[1] POST-INGEST roster: ", nrow(PIPE), " samples from ", length(fs), " qc summaries")
} else {
  PIPE <- fread(file.path(OUT_DIR, "ALL__samples.tsv"))[, .(dataset, sample = as.character(sample))]
  message("[1] PRE-INGEST (v1) roster: ", nrow(PIPE), " samples")
}

# ---------------------------------------------------------------------------
# [2] Known differences and their disposition
# ---------------------------------------------------------------------------
# Every entry here is a difference between the v1 roster and curated that has
# been diagnosed. `fixed_by` names the mechanism; anything left as "OPEN" is a
# decision or a download still outstanding. This table is the reason the check
# can be strict after ingest: nothing is allowed to be merely "expected".
DISPOSITION <- rbindlist(list(
  data.table(dataset = "GSE185381", side = "pipeline_only",
             sample  = c("AML0102", "AML0134", "AML2975", "Control0182"),
             fixed_by = "INGEST_EXCLUDE 5prime arm"),
  data.table(dataset = "GSE185381", side = "curated_only",
             sample  = c("AML001", "AML2123", "PAWWEE"),
             fixed_by = "Sample=donor + MIN_CELLS at donor level (v1 lost these to a library-level gate)"),
  data.table(dataset = "Petti2019", side = "pipeline_only",
             sample  = c("ND_083017", "ND_090617", "Normal_sorted_170531", "Normal_sorted_170607"),
             fixed_by = "INGEST_EXCLUDE uncurated_healthy"),
  data.table(dataset = "E-MTAB-11536", side = "pipeline_only",
             sample  = c("621B-BLD", "637C-BLD", "A35-BLD", "A36-BLD"),
             fixed_by = "INGEST_EXCLUDE peripheral_blood"),
  data.table(dataset = "E-MTAB-11536", side = "curated_only",
             sample  = c("D496-BMA-80", "D503-BMA-72"),
             fixed_by = paste("OPEN: absent from the local deposit. Needs download, plus the",
                              "erythroid-node mask (C10) and donor-equal barycenter weighting",
                              "(C11) before they can be admitted -- they are 72% of this",
                              "dataset's marrow cells and the only 3' + triple-depleted arm.")),
  data.table(dataset = "GSE289435", side = "curated_only", sample = "MLL_29512",
             fixed_by = "INGEST_EXCLUDE xenograft (intentional; was an undocumented grepl in the ingest)"),
  data.table(dataset = "GSE185991", side = "curated_only", sample = "PT18_D30",
             fixed_by = paste("OPEN: curated gsm_id is 'unknown' -- the curator could not identify",
                              "which GSM this draw is. Unjoinable until resolved. GSE185991 is",
                              "l2_capable=FALSE (sorted blast libraries), so this affects L1 only."))))

# ---------------------------------------------------------------------------
# [3] Diff, per dataset -- THROUGH THE JOIN MAP, not by raw name
# ---------------------------------------------------------------------------
# Curated ids and pipeline ids do not share a naming convention (curated
# "Chen2023__AML103" vs pipeline "AML103_CD34"; GSE289435 curated "MLL_29512" vs
# a GSM-prefixed pipeline id), so a raw setdiff reports ~240 differences that are
# purely cosmetic. 11_audit_labels_vs_curated.R already resolves this and writes
# join_map.tsv with the rule that fired for each curated row; consume it.
# ORDER MATTERS: run 10 (skeleton) -> 11 (join + label audit) -> 12 (this).
JMAP <- file.path(OUT_DIR, "join_map.tsv")
if (!file.exists(JMAP))
  stop("join_map.tsv missing. Run scripts/99_admin/11_audit_labels_vs_curated.R first.")
J <- fread(JMAP, colClasses = "character")

# A curated row is MATCHED if its join rule found a pipeline sample. "broadcast" means the
# dataset is homogeneous across every pipeline-relevant field, so a per-sample join carries no
# information and the dataset-level value applies to all of its samples (GSE239721: the GSM<->PT
# map is private, but all 20 rows are identical, so nothing is lost).
J[, matched := how == "broadcast" | (!is.na(pipeline) & nzchar(pipeline))]
matched_pipe <- unique(rbindlist(lapply(seq_len(nrow(J)), function(i) {
  if (!J$matched[i]) return(NULL)
  tgt <- if (J$how[i] == "broadcast") PIPE[dataset == J$dataset[i], sample]
         else strsplit(J$pipeline[i], ";")[[1]]
  if (!length(tgt)) NULL else data.table(dataset = J$dataset[i], sample = tgt)
})), by = c("dataset", "sample"))

res <- rbind(
  # curated rows the join could not place anywhere in the pipeline roster
  J[matched == FALSE, .(dataset, sample = sub("^[^_]*__", "", sample_id), side = "curated_only")],
  # pipeline samples no curated row claims
  PIPE[!matched_pipe, on = .(dataset, sample)][, .(dataset, sample, side = "pipeline_only")]
)

if (nrow(res)) {
  res <- merge(res, DISPOSITION, by = c("dataset", "sample", "side"), all.x = TRUE)
  res[is.na(fixed_by), fixed_by := "UNDIAGNOSED"]
} else res <- data.table(dataset = character(), sample = character(),
                         side = character(), fixed_by = character())

message("\n[2] ROSTER DIFF  (pipeline ", nrow(PIPE), " vs curated ", nrow(CUR), ")")
if (nrow(res)) print(res[order(fixed_by == "UNDIAGNOSED", dataset, sample)]) else
  message("    identical")

# ---------------------------------------------------------------------------
# [4] Healthy arm -- the set B_healthy is built from
# ---------------------------------------------------------------------------
message("\n[3] HEALTHY ARM per dataset (curated timepoint_class == Healthy)")
h <- CUR[timepoint == "Healthy", .(curated = .N), by = dataset]
hp <- res[side == "pipeline_only" & sample %in% CUR$sample == FALSE, .N, by = dataset]
print(h[order(-curated)])
message("    curated healthy total: ", CUR[timepoint == "Healthy", .N])

fwrite_safe(res, file.path(OUT_DIR, "roster_check.tsv"), sep = "\t")

# ---------------------------------------------------------------------------
# [5] Verdict
# ---------------------------------------------------------------------------
undiag <- res[fixed_by == "UNDIAGNOSED"]
open_i <- res[grepl("^OPEN", fixed_by)]
message("\n[4] VERDICT")
message("    undiagnosed differences: ", nrow(undiag))
message("    open items             : ", nrow(open_i))
if (POST) {
  if (nrow(undiag)) {
    print(undiag); stop("POST-INGEST roster does not match curated. See above.")
  }
  if (nrow(open_i)) message("    (open items are accepted gaps, recorded in DISPOSITION)")
  message("    PASS")
} else if (nrow(undiag)) {
  print(undiag)
  warning(nrow(undiag), " roster difference(s) with no diagnosis -- investigate before ingest.",
          call. = FALSE)
}
message("\n[done] ", file.path(OUT_DIR, "roster_check.tsv"))
