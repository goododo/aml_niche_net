#!/usr/bin/env Rscript
# =============================================================================
# 00_curated_manifest.R   --   the curated metadata becomes the single source
#
# WHY THIS EXISTS. Until now every sample-level fact lived in TWO places: the
# curated per-sample tables under /LARGE1/.../00_raw/metadata/, and a rule in
# scripts/config/ that re-derived the same fact from a sample name. Every defect
# found during the v2 re-run was that duplication failing:
#
#   - the timepoint vocabulary was migrated in TIMEPOINT_MAP but not in the
#     hardcoded GSE116256 day-code branch, so 19 samples carried a retired value
#   - TISSUE_OVERRIDE matched "^PT06_DX$", the CURATED id, while the pipeline id
#     is GSM5628167_M07, so a peripheral-blood draw stayed labelled marrow
#   - the metadata skeleton read a frozen 2026-07-22 side-car and silently
#     reinstated the v1 GSE185381 roster on every rebuild
#
# None of those were wrong RULES. They were rules that had drifted from the
# curation they were meant to encode. The fix is not better rules; it is to stop
# deriving what has already been curated.
#
# PRECEDENCE: curated value > config rule > NA. Every field carries its source,
# so "where did this come from" is answerable per sample per field.
#
# The config rules are NOT deleted. They stay as the fallback for samples the
# curation does not cover, and this script DIFFS them against curated so a drift
# becomes a reported disagreement instead of a silent divergence.
#
# OUTPUTS
#   00_curated_manifest.csv   one row per PIPELINE sample; resolved values +
#                             <field>_src provenance + conflict flags
#   00_curated_join_map.tsv   one row per CURATED row, with the derived
#                             pipeline_sample_id. The curation ships that column
#                             empty; this is the paste-back so it can be filled
#                             once and be authoritative from then on.
#   00_curated_diff.csv       every field where curated and the config rule
#                             disagree, with both values
#
# Usage: Rscript scripts/01_preprocess/00_curated_manifest.R
# =============================================================================

suppressPackageStartupMessages({library(data.table); library(readxl); library(here)})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))

CURATED_DIR <- "/LARGE1/gr10634/gaozy/aml_niche_net/00_raw/metadata"
OUT_MANIFEST <- file.path(DIR_PREPROCESS, "00_curated_manifest.csv")
OUT_JOIN     <- file.path(DIR_PREPROCESS, "00_curated_join_map.tsv")
OUT_DIFF     <- file.path(DIR_PREPROCESS, "00_curated_diff.csv")

# Fields carried forward. Deliberately a named list rather than "everything":
# a column nothing consumes is a column nobody validates.
CURATED_FIELDS <- c(
  # identity / provenance of the biological unit
  "patient_id", "gsm_id", "library_id", "n_libraries_merged",
  # clinical -- drives the analysis arms
  "disease", "control_type", "aml_ontogeny", "tissue", "timepoint_class",
  "days_from_diagnosis", "cycle_day_label", "blast_pct_clinical",
  # genetic -- drives subtype strata
  "karyotype", "fusion", "driver_mutations", "fab_class", "eln_risk",
  # technical -- drives node_status (which nodes CAN exist in this sample) and
  # the platform covariate. sorting / cell_prep / material_source are the three
  # that decide whether a missing node is biological or technical.
  "material_source", "platform", "sorting", "sorting_detail", "cell_prep",
  "sample_state", "is_multiplexed", "n_donors_in_library", "demux_method",
  # curation state -- an unverified row must be visible downstream
  "curation_status", "open_questions", "curator_note")

# Curated placeholders that mean "no value", normalised so downstream logic does
# not have to know the vocabulary of absence. The verbatim string is preserved in
# the join map; only the resolved manifest is normalised.
NULLISH <- c("", "NA", "na", "unknown", "not_applicable", "not_reported", "none_reported")
nz <- function(x) { x <- as.character(x); x[x %in% NULLISH] <- NA_character_; x }

# ---------------------------------------------------------------------------
# [1] Load every curated table
# ---------------------------------------------------------------------------
load_meta <- function(d) {
  p   <- file.path(CURATED_DIR, d)
  csv <- list.files(p, "^meta_.*\\.csv$", full.names = TRUE)
  if (length(csv)) return(fread(csv[1], colClasses = "character"))
  xls <- list.files(p, "^meta_.*\\.xlsx$", full.names = TRUE)[1]
  sh  <- grep("^meta_", readxl::excel_sheets(xls), value = TRUE)[1]
  as.data.table(readxl::read_excel(xls, sheet = sh, col_types = "text"))
}

dirs <- setdiff(list.dirs(CURATED_DIR, recursive = FALSE, full.names = FALSE), "BoneMarrowMap")
META <- lapply(setNames(dirs, gsub(" ", "", dirs)), load_meta)   # "Petti 2019" -> "Petti2019"
message("[1] loaded ", length(META), " curated tables, ",
        sum(vapply(META, nrow, integer(1))), " rows")

# ---------------------------------------------------------------------------
# [2] Pipeline roster
# ---------------------------------------------------------------------------
roster <- rbindlist(lapply(names(META), function(d) {
  f <- file.path(DIR_INGEST, paste0(d, "_qc_summary.csv"))
  if (!file.exists(f)) { message("    [warn] no ingest summary for ", d); return(NULL) }
  r <- fread(f)
  # UNDEMUX__* are cells the author metadata could not assign to a donor. The
  # demux prefilter removes them before the per-sample split, so they never
  # become samples and must not be offered to the join.
  data.table(dataset = d, sample = as.character(r$Sample))[!grepl("^UNDEMUX__", sample)]
}), fill = TRUE)
message("[2] pipeline roster: ", nrow(roster), " samples across ", uniqueN(roster$dataset), " datasets")

# ---------------------------------------------------------------------------
# [3] Derive the join, curated row -> pipeline sample(s)
# ---------------------------------------------------------------------------
# Rules in decreasing confidence. Each row records WHICH rule fired, so the join
# itself is auditable rather than just its result.
#   exact         curated sample_id (dataset prefix stripped) == pipeline sample
#   gsm_prefix    pipeline sample starts with the curated gsm_id
#   gsm_multi_N   curated row lists several ';'-separated GSMs -> N samples
#   suffix_strip  curated id carries a trailing -<digits> the pipeline id lacks
#   prefix_N      curated id prefixes N pipeline ids
#   broadcast     every curated row in the dataset is IDENTICAL across all
#                 pipeline-relevant fields, so a per-sample join carries no
#                 information (GSE239721: the GSM<->PT map is private, but all 20
#                 rows agree on every field the pipeline reads)
PIPE_FIELDS <- c("disease", "control_type", "timepoint_class", "tissue",
                 "sorting", "cell_prep", "material_source", "platform")

build_join <- function(d) {
  cu   <- META[[d]]
  real <- roster[dataset == d, sample]
  sid  <- sub("^[^_]*__", "", cu$sample_id)
  gsm  <- if ("gsm_id" %in% names(cu)) cu$gsm_id else rep(NA_character_, nrow(cu))

  fld <- intersect(PIPE_FIELDS, names(cu))
  homogeneous <- length(fld) > 0 &&
                 all(vapply(fld, function(c) uniqueN(cu[[c]]) == 1L, logical(1)))

  res <- data.table(dataset = d, sample_id = cu$sample_id, curated = sid, gsm = gsm,
                    pipeline = NA_character_, how = NA_character_)
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
    st <- sub("-[0-9]+$", "", sid[i]); h <- real[real == st]
    if (length(h) == 1) { res[i, `:=`(pipeline = h, how = "suffix_strip")]; next }

    h <- real[grepl(paste0("^", sid[i], "([_-]|$)"), real, ignore.case = TRUE)]
    if (length(h) >= 1) { res[i, `:=`(pipeline = paste(h, collapse = ";"),
                                      how = paste0("prefix_", length(h)))]; next }

    if (homogeneous) { res[i, `:=`(pipeline = paste(real, collapse = ";"), how = "broadcast")]; next }
    res[i, how := "UNRESOLVED"]
  }
  res
}

J <- rbindlist(lapply(names(META), build_join))
message("[3] join rules fired:")
print(J[, .N, by = .(how)][order(-N)])
unres <- J[how == "UNRESOLVED"]
if (nrow(unres)) {
  message("    ", nrow(unres), " curated row(s) UNRESOLVED (no data in the deposit, or a join gap):")
  print(unres[, .(dataset, sample_id, gsm)])
}

# ---------------------------------------------------------------------------
# [4] Explode to one row per PIPELINE sample, resolving N:1 collisions
# ---------------------------------------------------------------------------
# N curated rows CAN legitimately map to one pipeline sample: GSE185991 pools two
# patients into three libraries, so e.g. GSM5628191_M44 is claimed by both
# PT12_D14 and PT13_D14. Where those rows agree on a field, the value is
# unambiguous. Where they disagree, there is no defensible single value -- the
# field is set NA and flagged, rather than silently taking the first row.
long <- rbindlist(lapply(seq_len(nrow(J)), function(i) {
  if (J$how[i] == "UNRESOLVED" || is.na(J$pipeline[i])) return(NULL)
  tgt <- strsplit(J$pipeline[i], ";")[[1]]
  cu  <- META[[J$dataset[i]]]
  row <- cu[cu$sample_id == J$sample_id[i]]
  if (!nrow(row)) return(NULL)
  vals <- lapply(CURATED_FIELDS, function(f) if (f %in% names(row)) as.character(row[[f]][1]) else NA_character_)
  names(vals) <- CURATED_FIELDS
  cbind(data.table(dataset = J$dataset[i], sample = tgt,
                   curated_sample_id = J$sample_id[i], join_rule = J$how[i]),
        as.data.table(vals))
}), fill = TRUE)

collapse_field <- function(v) {
  u <- unique(nz(v)); u <- u[!is.na(u)]
  if (length(u) == 0) return(NA_character_)
  if (length(u) == 1) return(u)
  paste0("CONFLICT:", paste(sort(u), collapse = "|"))
}

M <- long[, c(list(n_curated_rows = .N,
                   curated_sample_id = paste(unique(curated_sample_id), collapse = ";"),
                   join_rule = paste(unique(join_rule), collapse = ";")),
              lapply(.SD, collapse_field)),
          by = .(dataset, sample), .SDcols = CURATED_FIELDS]

# Free-text fields conflict on every N:1 join by construction (two curators' notes
# are never byte-identical) and printing them buries the structured conflicts that
# actually need a decision. Collapsed and flagged like everything else; just not
# printed.
FREETEXT <- c("curator_note", "open_questions", "sorting_detail", "driver_mutations", "karyotype")
conflict_cols <- setdiff(
  CURATED_FIELDS[vapply(CURATED_FIELDS, function(f) any(grepl("^CONFLICT:", M[[f]])), logical(1))],
  FREETEXT)
if (length(conflict_cols)) {
  message("[4] structured fields with N:1 conflicts: ", paste(conflict_cols, collapse = ", "))
  for (f in conflict_cols)
    print(M[grepl("^CONFLICT:", get(f)),
            .(dataset, sample, curated_sample_id, field = f, value = substr(get(f), 1, 60))])
} else message("[4] no structured N:1 conflicts")

# Samples the curation does not cover at all.
M <- merge(roster, M, by = c("dataset", "sample"), all.x = TRUE)
M[is.na(n_curated_rows), `:=`(n_curated_rows = 0L, join_rule = "uncurated")]
message("    ", M[n_curated_rows == 0, .N], " pipeline sample(s) with no curated row")

# ---------------------------------------------------------------------------
# [5] Diff curated against the config rules
# ---------------------------------------------------------------------------
# The rules stay in place as the fallback. Diffing them here is what turns a
# silent divergence into a reported one.
M[, cfg_tissue := vapply(seq_len(.N), function(i) tissue_of(dataset[i], sample[i]), character(1))]
M[, cfg_platform := vapply(seq_len(.N), function(i) platform_of(dataset[i], sample[i])$platform, character(1))]

ing <- rbindlist(lapply(unique(M$dataset), function(d) {
  f <- file.path(DIR_INGEST, paste0(d, "_qc_summary.csv"))
  if (!file.exists(f)) return(NULL)
  fread(f)[, .(dataset = d, sample = as.character(Sample),
               cfg_timepoint = as.character(Timepoint))]
}), fill = TRUE)
M <- merge(M, ing, by = c("dataset", "sample"), all.x = TRUE)

# timepoint_class is the CURATED vocabulary and CANONICAL_TIMEPOINTS is now the
# same vocabulary, so these are directly comparable -- that alignment was the
# point of the v2 vocabulary migration.
diffs <- rbindlist(list(
  M[!is.na(timepoint_class) & !is.na(cfg_timepoint) & timepoint_class != cfg_timepoint,
    .(dataset, sample, field = "timepoint", curated = timepoint_class, config = cfg_timepoint)],
  M[!is.na(tissue) & tissue != cfg_tissue,
    .(dataset, sample, field = "tissue", curated = tissue, config = cfg_tissue)]
), fill = TRUE)

message("\n[5] curated vs config-rule disagreements: ", nrow(diffs))
if (nrow(diffs)) print(diffs)

# ---------------------------------------------------------------------------
# [6] Resolve with precedence, and record the source of every resolved field
# ---------------------------------------------------------------------------
resolve <- function(cur, cfg, name) {
  out <- nz(cur); src <- rep("curated", length(out))
  miss <- is.na(out); out[miss] <- cfg[miss]; src[miss] <- "config_rule"
  src[is.na(out)] <- "none"
  list(value = out, src = src)
}
r_tp <- resolve(M$timepoint_class, M$cfg_timepoint, "timepoint")
r_ti <- resolve(M$tissue,          M$cfg_tissue,    "tissue")
M[, `:=`(timepoint = r_tp$value, timepoint_src = r_tp$src,
         tissue_r  = r_ti$value, tissue_src    = r_ti$src)]

# PLATFORM IS DELIBERATELY NOT MERGED. The two sources answer different questions
# and are not the same vocabulary:
#   platform_coarse  config's platform_of() -> "10x" / "SeqWell". This is the
#                    DROPLET-vs-not distinction the doublet model keys on
#                    (DBL_RATE_FLAT has a SeqWell entry; the 0.8%/1k law is a 10x
#                    loading model and is simply wrong for Seq-Well).
#   platform_fine    curated -> "10x_3p_v3", "10x_5p_v1.1", ... This is the
#                    chemistry covariate, and it is what distinguishes the two
#                    arms inside E-MTAB-11536 or the 3'/5' split in GSE185381.
# Collapsing them into one column would have silently mixed "10x" and "10x_3p_v3"
# in the same field -- which is how a covariate stops stratifying anything.
M[, `:=`(platform_coarse = cfg_platform, platform_fine = nz(platform))]
M[, platform_fine_src := fifelse(is.na(platform_fine), "missing", "curated")]

# sorting / cell_prep get the same precedence chain, because node_status reads them
# and a blank there costs more than a wrong-but-checkable value (see SORTING_FALLBACK).
M[, `:=`(sorting_r = nz(sorting), cell_prep_r = nz(cell_prep))]
M[, `:=`(sorting_src = fifelse(is.na(sorting_r), "none", "curated"),
         cell_prep_src = fifelse(is.na(cell_prep_r), "none", "curated"))]
for (i in seq_len(nrow(SORTING_FALLBACK))) {
  d <- SORTING_FALLBACK$dataset[i]
  M[dataset == d & is.na(sorting_r),
    `:=`(sorting_r = SORTING_FALLBACK$sorting[i], sorting_src = "config_rule")]
  M[dataset == d & is.na(cell_prep_r),
    `:=`(cell_prep_r = SORTING_FALLBACK$cell_prep[i], cell_prep_src = "config_rule")]
}
message("    sorting resolved: ", paste(sprintf("%s=%d", names(table(M$sorting_src)),
                                                table(M$sorting_src)), collapse = "  "))
if (M[is.na(sorting_r), .N])
  warning(M[is.na(sorting_r), .N], " sample(s) still have no sorting value -- node_status ",
          "must treat these as NA_technical for every node, which is very costly. Datasets: ",
          paste(unique(M[is.na(sorting_r), dataset]), collapse = ", "), call. = FALSE)

message("\n[6] provenance of resolved fields:")
for (f in c("timepoint", "tissue"))
  print(M[, .N, by = c(paste0(f, "_src"))][order(-N)])
message("    platform_fine missing for ", M[is.na(platform_fine), .N], " sample(s): ",
        paste(unique(M[is.na(platform_fine), dataset]), collapse = ", "))

# ---------------------------------------------------------------------------
# [7] Write
# ---------------------------------------------------------------------------
fwrite_safe(M, OUT_MANIFEST)
# Paste-back: the curation ships pipeline_sample_id empty. This is the derived
# value, so it can be filled once and stop being derived.
fwrite_safe(J[, .(dataset, sample_id, gsm, pipeline_sample_id = pipeline, join_rule = how)],
            OUT_JOIN, sep = "\t")
fwrite_safe(diffs, OUT_DIFF)

message("\n[done]")
message("  manifest : ", OUT_MANIFEST, "  (", nrow(M), " samples x ", ncol(M), " fields)")
message("  join map : ", OUT_JOIN, "  <- paste pipeline_sample_id back into the curation")
message("  diff     : ", OUT_DIFF, "  (", nrow(diffs), " disagreements)")
