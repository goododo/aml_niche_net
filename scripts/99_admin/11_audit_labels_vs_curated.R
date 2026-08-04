# =============================================================================
# 11_audit_labels_vs_curated.R
#
# Compares every label the pipeline currently assigns against the curated
# per-sample metadata, and reports EVERY disagreement.
#
# WHY: the pipeline derives Timepoint and disease state from dataset-level rules
# plus sample-name regexes. That inference is what put 8 Chen2023 normal marrows
# into the AML arm. This script answers the only question that matters before a
# re-run: is that the ONLY case, or are there others?
#
# Join strategy, in order of decreasing confidence. Each row records which rule
# fired so a reviewer can audit the join itself, not just its result:
#   exact          curated sample_id (prefix stripped) == pipeline Sample
#   gsm_prefix     pipeline Sample starts with the curated gsm_id
#   gsm_multi      curated row carries several ';'-separated GSMs -> several
#                  pipeline samples (one curated row merges >1 library)
#   suffix_strip   curated id has a trailing -<digits> the pipeline id lacks
#   prefix_N       curated id prefixes N pipeline ids (compartment split)
#   broadcast      every curated row in the dataset is IDENTICAL across all
#                  pipeline-relevant fields, so a per-sample join carries no
#                  information -- the dataset-level value is applied directly.
#                  Used for GSE239721, whose GSM<->PT map is private and
#                  therefore unrecoverable, but whose 20 rows are homogeneous.
#
# OUTPUT: 00_project/metadata/join_map.tsv        one row per curated row
#         00_project/metadata/label_audit.tsv     one row per pipeline sample
#         console: the disagreements, grouped by severity
#
# Usage: Rscript scripts/99_admin/11_audit_labels_vs_curated.R
# =============================================================================

suppressPackageStartupMessages({library(data.table); library(readxl); library(here)})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))

CURATED_DIR <- "/LARGE1/gr10634/gaozy/aml_niche_net/00_raw/metadata"
OUT_DIR     <- file.path(PROJECT_ROOT, "00_project", "metadata")

# Fields the pipeline actually consumes. `broadcast` is only legitimate when a
# dataset is homogeneous across ALL of them.
PIPE_FIELDS <- c("disease", "control_type", "timepoint_class", "tissue",
                 "sorting", "cell_prep", "material_source", "platform")

load_meta <- function(d) {
  p   <- file.path(CURATED_DIR, d)
  csv <- list.files(p, "^meta_.*\\.csv$", full.names = TRUE)
  if (length(csv)) return(fread(csv[1], colClasses = "character"))
  xls <- list.files(p, "^meta_.*\\.xlsx$", full.names = TRUE)[1]
  sh  <- grep("^meta_", readxl::excel_sheets(xls), value = TRUE)[1]
  as.data.table(readxl::read_excel(xls, sheet = sh, col_types = "text"))
}

dirs <- setdiff(list.dirs(CURATED_DIR, recursive = FALSE, full.names = FALSE), "BoneMarrowMap")
META <- lapply(setNames(dirs, gsub(" ", "", dirs)), load_meta)

skel <- fread(file.path(OUT_DIR, "ALL__samples.tsv"))

# ---------------------------------------------------------------------------
# [1] Join map
# ---------------------------------------------------------------------------
build_join <- function(d) {
  cu   <- META[[d]]
  real <- skel[dataset == d, sample]
  sid  <- sub("^[^_]*__", "", cu$sample_id)
  gsm  <- if ("gsm_id" %in% names(cu)) cu$gsm_id else rep(NA_character_, nrow(cu))

  # homogeneity test -> broadcast
  fld <- intersect(PIPE_FIELDS, names(cu))
  homogeneous <- length(fld) > 0 &&
                 all(vapply(fld, function(c) length(unique(cu[[c]])) == 1L, logical(1)))

  res <- data.table(dataset = d, sample_id = cu$sample_id, curated = sid,
                    gsm = gsm, pipeline = NA_character_, how = NA_character_)
  for (i in seq_len(nrow(res))) {
    h <- real[real == sid[i]]
    if (length(h) == 1) { res[i, `:=`(pipeline = h, how = "exact")]; next }

    gs <- if (is.na(gsm[i])) character(0) else trimws(strsplit(gsm[i], ";")[[1]])
    gs <- gs[grepl("^GSM", gs)]
    if (length(gs)) {
      h <- unlist(lapply(gs, function(g) real[grepl(paste0("^", g, "[_-]"), real)]))
      if (length(h) == 1) { res[i, `:=`(pipeline = h, how = "gsm_prefix")]; next }
      if (length(h) > 1)  { res[i, `:=`(pipeline = paste(h, collapse = ";"),
                                        how = paste0("gsm_multi_", length(h)))]; next }
    }
    # curated id longer than the pipeline id by a trailing -<digits>
    st <- sub("-[0-9]+$", "", sid[i])
    h  <- real[real == st]
    if (length(h) == 1) { res[i, `:=`(pipeline = h, how = "suffix_strip")]; next }

    h <- real[grepl(paste0("^", sid[i], "([_-]|$)"), real, ignore.case = TRUE)]
    if (length(h) >= 1) { res[i, `:=`(pipeline = paste(h, collapse = ";"),
                                      how = paste0("prefix_", length(h)))]; next }

    if (homogeneous) { res[i, how := "broadcast"]; next }
    res[i, how := "UNRESOLVED"]
  }
  attr(res, "homogeneous") <- homogeneous
  res
}

J <- rbindlist(lapply(names(META), build_join))
fwrite_safe(J, file.path(OUT_DIR, "join_map.tsv"), sep = "\t")

message("[1] join map:")
print(J[, .N, by = .(dataset, how)][order(dataset, -N)])

# ---------------------------------------------------------------------------
# [2] Curated label per PIPELINE sample (explode the 1:N joins)
# ---------------------------------------------------------------------------
get <- function(d, sid, col) {
  x <- META[[d]]; if (!col %in% names(x)) return(NA_character_)
  v <- x[[col]][x$sample_id == sid]; if (length(v) == 1) v else NA_character_
}

rows <- list()
for (i in seq_len(nrow(J))) {
  d <- J$dataset[i]
  tgt <- if (J$how[i] == "broadcast") skel[dataset == d, sample]
         else if (is.na(J$pipeline[i])) character(0)
         else strsplit(J$pipeline[i], ";")[[1]]
  if (!length(tgt)) next
  rows[[length(rows) + 1]] <- data.table(
    dataset = d, sample = tgt, how = J$how[i], sample_id = J$sample_id[i],
    cur_timepoint = get(d, J$sample_id[i], "timepoint_class"),
    cur_disease   = get(d, J$sample_id[i], "disease"),
    cur_control   = get(d, J$sample_id[i], "control_type"),
    cur_sorting   = get(d, J$sample_id[i], "sorting"),
    cur_material  = get(d, J$sample_id[i], "material_source"))
}
C <- unique(rbindlist(rows), by = c("dataset", "sample"))

# ---------------------------------------------------------------------------
# [3] Compare against what the pipeline assigns
# ---------------------------------------------------------------------------
# The curated vocabulary is RICHER than CANONICAL_TIMEPOINTS. The audit therefore
# compares at the coarse level where a mislabel actually changes an analysis arm:
# healthy vs diseased, and diagnosis vs post-treatment vs relapse. Finer values
# (Post_induction vs Post_consolidation) are a vocabulary upgrade, not an error.
coarse <- function(x) fifelse(is.na(x), NA_character_,
  fifelse(x == "Healthy", "Healthy",
  fifelse(x == "Diagnosis", "Diagnosis",
  fifelse(x %in% c("Relapse", "Relapse2"), "Relapse",   # our canonical splits 1st/2nd relapse
  fifelse(x %in% c("Refractory"), "Refractory", "Post_treatment")))))

A <- merge(skel[, .(dataset, sample, pipe_timepoint = timepoint, pipe_v1 = timepoint_v1)],
           C, by = c("dataset", "sample"), all.x = TRUE)
A[, cur_coarse  := coarse(cur_timepoint)]
A[, pipe_coarse := coarse(pipe_timepoint)]
A[, mismatch := !is.na(cur_coarse) & !is.na(pipe_coarse) & cur_coarse != pipe_coarse]
A[, unjoined := is.na(cur_timepoint)]
setorder(A, -mismatch, dataset, sample)
fwrite_safe(A, file.path(OUT_DIR, "label_audit.tsv"), sep = "\t")

message("\n[3] AUDIT RESULT")
message("    pipeline samples: ", nrow(A),
        " | joined: ", A[unjoined == FALSE, .N], " | unjoined: ", A[unjoined == TRUE, .N])

mm <- A[mismatch == TRUE]
if (nrow(mm)) {
  message("\n    *** ", nrow(mm), " LABEL DISAGREEMENTS ***")
  print(mm[, .(dataset, sample, pipeline = pipe_coarse, curated = cur_coarse,
               curated_detail = cur_timepoint, disease = cur_disease)])
} else message("\n    no coarse-level disagreements")

message("\n    healthy-arm membership:")
print(A[, .(pipeline = sum(pipe_coarse == "Healthy", na.rm = TRUE),
            curated  = sum(cur_coarse  == "Healthy", na.rm = TRUE)), by = dataset][
        pipeline > 0 | curated > 0])

if (A[unjoined == TRUE, .N])
  print(A[unjoined == TRUE, .N, by = dataset])
message("\n[done] ", file.path(OUT_DIR, "label_audit.tsv"))
