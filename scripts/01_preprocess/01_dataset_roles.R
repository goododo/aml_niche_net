# =============================================================================
# 01_dataset_roles.R   (Task 1)
# Build a sample-level role manifest by joining:
#   (a) the curated STUDY-LEVEL role table (ROLE_TABLE), and
#   (b) per-sample metadata extracted from each merged RDS object,
# then applying sample-level refinement (healthy / cell-line / PDX / longitudinal).
#
# Output: results/tables/01_qc/01_sample_role_manifest.csv
# Run:    Rscript 01_dataset_roles.R           (reads RDS metadata, then caches)
#         Rscript 01_dataset_roles.R --cache   (reuse cached manifest, skip RDS I/O)
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

args        <- commandArgs(trailingOnly = TRUE)
use_cache   <- "--cache" %in% args
manifest_raw_path <- file.path(DIR_PREPROCESS, "01_sample_metadata_raw.csv")

# ----------------------------------------------------------------------------
# Step 1. Per-sample metadata (one row per Sample) from each RDS.
#         meta.data is small relative to counts but we still load each object
#         once; result is cached so this only runs when --cache is absent.
# ----------------------------------------------------------------------------
extract_sample_meta <- function() {
  datasets <- list_datasets()
  out <- vector("list", length(datasets))
  need <- c(COL_SAMPLE, COL_DATASET, COL_STUDY,
            COL_PATIENT, COL_TP, COL_DISEASE)
  for (i in seq_along(datasets)) {
    ds   <- datasets[i]
    message_ts("reading metadata: ", ds)
    seu  <- readRDS(file.path(RDS_INGEST_DIR, paste0(ds, ".rds")))
    md   <- as.data.table(seu@meta.data)
    # Guard: ensure expected columns exist; fill missing with NA.
    for (cn in need) if (!cn %in% names(md)) md[[cn]] <- NA_character_
    md   <- md[, ..need]
    # Collapse to one row per Sample; n_cells from object (raw, pre-QC).
    smry <- md[, .(.n_cells_obj = .N,
                   Disease_state = data.table::first(get(COL_DISEASE)),
                   Patient_ID    = data.table::first(get(COL_PATIENT)),
                   Timepoint     = data.table::first(get(COL_TP)),
                   Study         = data.table::first(get(COL_STUDY))),
                by = c(Sample = COL_SAMPLE)]
    smry[, dataset := ds]
    out[[i]] <- smry
    rm(seu); gc(verbose = FALSE)
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

if (use_cache && file.exists(manifest_raw_path)) {
  message_ts("using cached sample metadata: ", manifest_raw_path)
  smeta <- fread(manifest_raw_path)
} else {
  smeta <- extract_sample_meta()
  fwrite_safe(smeta, manifest_raw_path)
  message_ts("wrote raw sample metadata: ", manifest_raw_path)
}

# ----------------------------------------------------------------------------
# Step 2. Join the curated study-level role table.
# ----------------------------------------------------------------------------
roles <- copy(ROLE_TABLE)
m <- merge(smeta, roles, by = "dataset", all.x = TRUE)
if (anyNA(m$study_role)) {
  warning("Datasets without a curated role: ",
          paste(unique(m[is.na(study_role), dataset]), collapse = ", "))
}

# ----------------------------------------------------------------------------
# Step 3. Sample-level role refinement.
#   - healthy/normal/control/donor  -> Healthy-control (pooled for B_healthy)
#   - cell line / PDX               -> Exclude
#   - everything else inherits the study role
# Note: healthy donors inside mixed studies (Petti/Lasry/Chen/vanGalen) are
# different individuals, so pooling them as controls causes no train/val leak.
# ----------------------------------------------------------------------------
ds_text <- tolower(paste(m$Disease_state, m$Sample))
m[, sample_role := study_role]
m[grepl(HEALTHY_REGEX,  ds_text), sample_role := "Healthy-control"]
m[grepl(CELLLINE_REGEX, ds_text), sample_role := "Exclude-cellline"]
m[grepl(PDX_REGEX,      ds_text), sample_role := "Exclude-PDX"]

# Discovery-candidacy is a property of NON-healthy, NON-excluded AML samples
# within discovery_candidate studies (longitudinal studies are never candidates).
m[, in_discovery_pool :=
    discovery_candidate == TRUE &
    is_longitudinal == FALSE &
    !sample_role %in% c("Healthy-control", "Exclude-cellline", "Exclude-PDX")]

# Surface integrity warnings at the sample level too.
m[, integrity_flag := fifelse(is.na(integrity_flag), "OK", integrity_flag)]

# Namespaced patient key: prevents cross-dataset ID collisions
# (e.g. GSE239721 "PT10A" is NOT GSE185991 "PT10"). Use uid_patient for any
# cross-dataset / leakage-control logic; never join on raw Patient_ID alone.
m[, uid_patient := paste(dataset, Patient_ID, sep = ":")]

# ----------------------------------------------------------------------------
# Step 4. Write manifest + a compact study x role summary.
# ----------------------------------------------------------------------------
keep <- c("dataset","Study","Sample","Patient_ID","uid_patient","Timepoint",
          "Disease_state",".n_cells_obj","study_role","sample_role",
          "is_longitudinal","discovery_candidate","in_discovery_pool",
          "subtype_stratum","l2_capable","integrity_flag")
manifest <- m[, ..keep]
setorder(manifest, dataset, Sample)
fwrite_safe(manifest, file.path(DIR_PREPROCESS, "01_sample_role_manifest.csv"))

study_summary <- manifest[, .(n_samples = .N,
                              n_discovery_pool = sum(in_discovery_pool, na.rm = TRUE),
                              n_healthy = sum(sample_role == "Healthy-control"),
                              n_excluded = sum(grepl("^Exclude", sample_role)),
                              integrity = data.table::first(integrity_flag)),
                          by = .(dataset, study_role, is_longitudinal,
                                 discovery_candidate, subtype_stratum)]
setorder(study_summary, -discovery_candidate, dataset)
fwrite_safe(study_summary, file.path(DIR_PREPROCESS, "01_study_role_summary.csv"))

message_ts("Task 1 done.")
print(study_summary)
