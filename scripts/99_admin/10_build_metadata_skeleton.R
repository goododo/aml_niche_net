# =============================================================================
# 10_build_metadata_skeleton.R
#
# Emits one per-sample metadata TSV per dataset, pre-filled with the sample ids
# the pipeline ACTUALLY produces plus every value the code can currently infer,
# for manual curation.
#
# WHY: the sample-level facts the pipeline needs (timepoint, disease state,
# subtype, platform, sorted-vs-whole-marrow compartment) are currently spread
# across hardcoded tables in config_qc.R -- TIMEPOINT_MAP, TP_HEALTHY_SAMPLE_
# PATTERNS, PLATFORM_TABLE, ROLE_TABLE$subtype_stratum -- and inferred from
# sample-name regexes. That is how 8 Chen2023 normal marrows ended up labelled
# Diagnosis, and how a 50-donor deposit ended up keyed by library. Curated
# per-sample metadata replaces inference with fact.
#
# THE JOIN KEY IS `sample`. It must match the id ingest_<dataset>.R writes into
# the Sample column, or the curated row silently fails to join and the pipeline
# falls back to the inferred value. So the skeleton is GENERATED from the real
# sample lists rather than typed by hand.
#
# NOTE on GSE185381: this dataset's Sample key changed from library to DONOR
# (a library pools up to 5 donors and mixes AML with Healthy). The skeleton
# therefore lists the 53 donors, not the 50 libraries -- the v1 ingest QC
# summary for this dataset is keyed the old way and must NOT be used.
#
# review_flag marks every field the code could not establish, so curation
# effort goes where it is actually needed.
#
# Usage: Rscript scripts/99_admin/10_build_metadata_skeleton.R
# =============================================================================

suppressPackageStartupMessages({ library(here); library(data.table) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))

OUT_DIR <- file.path(PROJECT_ROOT, "00_project", "metadata")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Compartment inference from the sample id. This is the field that SHOULD drive
# l2_capable: a CD34/CD117-sorted library has no stroma/immune compartment, so
# no whole-marrow CCC graph can be built from it. Today that decision is a
# hand-written l2_capable column plus a Chen2023-only regex
# (CCC_SORTED_SUBLIB_PATTERN) -- two mechanisms for one fact.
# ---------------------------------------------------------------------------
COMPARTMENT_VALUES <- c("whole_marrow", "niche_immune", "CD34_sorted",
                        "CD117_sorted", "lineage_sorted", "unknown")

infer_compartment <- function(dataset, sample) {
  s <- tolower(as.character(sample))
  d <- as.character(dataset)
  if (grepl("cd34", s))                          return("CD34_sorted")
  if (grepl("cd117", s))                         return("CD117_sorted")
  if (grepl("niche|immune", s))                  return("niche_immune")
  if (grepl("34p|38n", s))                       return("CD34_sorted")   # GSE116256 BM5-34p / -34p38n
  if (grepl("sorted", s))                        return("lineage_sorted")
  # dataset-level protocol facts recorded in ROLE_TABLE$l2_reason
  if (d %in% c("GSE185991", "GSE147989"))        return("CD34_sorted")
  "unknown"
}

DISEASE_VALUES  <- c("Healthy", "AML", "MDS", "Other", "Unknown")
SUBTYPE_VALUES  <- c("NPM1", "FLT3-ITD", "TP53", "KMT2A-r", "monocytic",
                     "complex_karyotype", "other", "NA")

# ---------------------------------------------------------------------------
# Sample lists: from the ingest QC summaries, except GSE185381 (key changed).
# ---------------------------------------------------------------------------
build_rows <- function(ds) {
  if (ds == "GSE185381") {
    f <- file.path(PROJECT_ROOT, "recon_GSE185381_donor_qc.csv")
    if (!file.exists(f)) {
      message("[skip] GSE185381: donor recon table missing (", f, ")"); return(NULL)
    }
    r <- fread(f)
    return(data.table(
      dataset          = ds,
      sample           = r$Patient_ID,        # donor IS the sample now
      patient_id       = r$Patient_ID,
      n_cells_ingest   = r$n_cells,
      timepoint        = r$timepoints,        # from per-cell Disease_state
      timepoint_detail = r$disease_states,
      disease_state    = r$disease_states,
      n_libraries      = r$n_libraries,
      med_nFeature     = r$med_nFeature))
  }
  f <- file.path(DIR_INGEST, paste0(ds, "_qc_summary.csv"))
  if (!file.exists(f)) { message("[skip] ", ds, ": no ingest QC summary"); return(NULL) }
  r <- fread(f)
  data.table(dataset          = ds,
             sample           = as.character(r$Sample),
             patient_id       = as.character(r$Patient_ID),
             n_cells_ingest   = r$n_cells,
             timepoint        = as.character(r$Timepoint),
             timepoint_detail = NA_character_,
             disease_state    = NA_character_,
             n_libraries      = NA_integer_,
             med_nFeature     = r$median_nFeature)
}

datasets <- ROLE_TABLE$dataset
all_rows <- rbindlist(lapply(datasets, build_rows), use.names = TRUE, fill = TRUE)

# ---------------------------------------------------------------------------
# Fill everything the code can establish; flag the rest for curation.
# ---------------------------------------------------------------------------
all_rows[, compartment := mapply(infer_compartment, dataset, sample)]
all_rows[, timepoint_days := timepoint_days(timepoint_detail)]

# The timepoint carried in the ingest QC summaries is the v1 (2026-07-15) value, produced BEFORE
# rule 0 existed -- so the 8 Chen2023 NBM samples still read "Diagnosis" there. Keep it as
# timepoint_v1 and overlay rule 0, the one rule that decides from sample IDENTITY and is a proven
# bug. The other canonicalisation rules are deliberately NOT re-applied: they need the RAW label,
# which the summaries no longer carry (re-running GSE116256's day-code rule against an
# already-canonical "Post_treatment" would map it to Unknown). Everything else stays at the v1
# value for curation to confirm or correct.
setnames(all_rows, "timepoint", "timepoint_v1")
all_rows[, rule0_healthy := mapply(is_healthy_sample, sample, patient_id)]
all_rows[, timepoint := fifelse(rule0_healthy, "Healthy", timepoint_v1)]
all_rows[, changed_vs_v1 := timepoint != timepoint_v1]
if (all_rows[changed_vs_v1 == TRUE, .N] > 0) {
  message("[rule0] sample-level healthy donors reclassified out of the AML arm:")
  print(all_rows[changed_vs_v1 == TRUE, .(dataset, sample, timepoint_v1, timepoint)])
}

plat <- rbindlist(lapply(seq_len(nrow(all_rows)), function(i) {
  p <- platform_of(all_rows$dataset[i])
  data.table(platform = p$platform, chemistry = p$chemistry)
}))
all_rows[, `:=`(platform = plat$platform, chemistry = plat$chemistry)]

all_rows <- merge(all_rows, ROLE_TABLE[, .(dataset, subtype_stratum, study_role)],
                  by = "dataset", all.x = TRUE)
# subtype_stratum is a DATASET-level guess; per-sample subtype is what H1's
# conditional-SRRS strata actually need, so it is carried as a starting value only.
setnames(all_rows, "subtype_stratum", "subtype_dataset_level")
all_rows[, subtype := NA_character_]
all_rows[, karyotype := NA_character_]
all_rows[, mutations := NA_character_]
all_rows[, upstream_doublet_removed := is_upstream_filtered(dataset)]
all_rows[, notes := NA_character_]

# review_flag: the fields curation must resolve, most-consequential first.
all_rows[, review_flag := ""]
all_rows[timepoint %in% c("Unknown", "", NA), review_flag := paste0(review_flag, "timepoint;")]
all_rows[is.na(disease_state) | disease_state == "Unknown",
         review_flag := paste0(review_flag, "disease_state;")]
all_rows[compartment == "unknown",  review_flag := paste0(review_flag, "compartment;")]
all_rows[chemistry  == "unknown",   review_flag := paste0(review_flag, "chemistry;")]
all_rows[is.na(subtype),            review_flag := paste0(review_flag, "subtype;")]

col_order <- c("dataset","sample","patient_id","timepoint","timepoint_v1","changed_vs_v1",
               "timepoint_detail","timepoint_days",
               "disease_state","subtype","subtype_dataset_level","karyotype","mutations",
               "compartment","platform","chemistry","upstream_doublet_removed",
               "n_cells_ingest","med_nFeature","n_libraries","study_role","notes","review_flag")
setcolorder(all_rows, intersect(col_order, names(all_rows)))
setorder(all_rows, dataset, sample)

for (ds in unique(all_rows$dataset)) {
  f <- file.path(OUT_DIR, paste0(ds, "__samples.tsv"))
  fwrite_safe(all_rows[dataset == ds], f, sep = "\t")
  message(sprintf("  %-14s %3d samples -> %s", ds, all_rows[dataset == ds, .N], basename(f)))
}
fwrite_safe(all_rows, file.path(OUT_DIR, "ALL__samples.tsv"), sep = "\t")

# Controlled vocabularies, written beside the tables so curation has one place to look.
vocab <- rbindlist(list(
  data.table(field = "timepoint",     allowed = CANONICAL_TIMEPOINTS),
  data.table(field = "disease_state", allowed = DISEASE_VALUES),
  data.table(field = "subtype",       allowed = SUBTYPE_VALUES),
  data.table(field = "compartment",   allowed = COMPARTMENT_VALUES),
  data.table(field = "platform",      allowed = unique(PLATFORM_TABLE$platform)),
  data.table(field = "chemistry",     allowed = c("3prime","5prime","NA","unknown"))))
fwrite_safe(vocab, file.path(OUT_DIR, "VOCABULARY.tsv"), sep = "\t")

message("\n[summary] ", nrow(all_rows), " samples across ", uniqueN(all_rows$dataset), " datasets")
message("[summary] fields needing curation:")
print(all_rows[, .(n = .N), by = review_flag][order(-n)])
message("\n[done] ", OUT_DIR)
