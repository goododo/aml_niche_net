# ============================================================================
# ingest_GSE185991.R   --   Naldini 2023 Nat Commun (10x 3'; sorted blast/prog)
# Format : prefixed 10x mtx  "<GSM>_<M##>_{barcodes,features,matrix}"
# Patient: hardcoded GSM -> label map (provided by PI). Label parsed as
#          patient + (DX|D14|D30|REL). Pooled labels (contain "_", e.g.
#          "PT01_PT10D30") have no per-cell donor file in this dataset, so
#          per the agreed fallback they keep orig.ident = composition string
#          and are flagged Pooled = TRUE.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE185991"
STUDY   <- "Naldini 2023 Nat Commun"
RAW     <- file.path(GEO_RAW_DIR, "GSE185991_RAW")

# ----------------------------------------------------------------------------
# [1] Hardcoded GSM -> patient/timepoint label map ----
# ----------------------------------------------------------------------------
PT_MAP <- c(
  GSM5628163 = "PT01DX",  GSM5628164 = "PT01DX",  GSM5628165 = "PT01D14", GSM5628166 = "PT01D14",
  GSM5628167 = "PT06DX",  GSM5628168 = "PT06DX",  GSM5628169 = "PT06D14", GSM5628170 = "PT06D14",
  GSM5628171 = "PT07DX",  GSM5628172 = "PT07DX",  GSM5628173 = "PT07D14", GSM5628174 = "PT07D14",
  GSM5628175 = "PT02DX",  GSM5628176 = "PT02D30", GSM5628177 = "PT01D30",
  GSM5628178 = "PT08DX",  GSM5628179 = "PT08D14", GSM5628180 = "PT08REL", GSM5628181 = "PT08REL",
  GSM5628182 = "PT09DX",  GSM5628183 = "PT10DX",  GSM5628184 = "PT09D14", GSM5628185 = "PT10D14",
  GSM5628186 = "PT08D30", GSM5628187 = "PT12DX",  GSM5628188 = "PT13DX",  GSM5628189 = "PT10REL",
  GSM5628190 = "PT01_PT10D30", GSM5628191 = "PT12_PT13D14", GSM5628192 = "PT12D30",
  GSM5628195 = "PT15DX",  GSM5628196 = "PT15REL", GSM5628197 = "PT09_PT13D30",
  GSM6412447 = "PT11DX",  GSM6412448 = "PT11D30", GSM6412449 = "PT17DX",  GSM6412450 = "PT17D30",
  GSM6412451 = "PT18DX",  GSM6412452 = "PT18D14", GSM6412453 = "PT17D14",
  GSM6412454 = "PT19REL", GSM6412455 = "PT20REL"
)

# Split a label into patient prefix + trailing timepoint token.
parse_label <- function(lab) {
  tp  <- stringr::str_extract(lab, "(DX|D\\d+|REL)$")
  pid <- sub("(DX|D\\d+|REL)$", "", lab)
  list(patient = pid, timepoint = tp)
}

# ----------------------------------------------------------------------------
# [2] Discover library prefixes ----
# ----------------------------------------------------------------------------
message("[2] listing library prefixes")
bc_files <- list.files(RAW, pattern = "_barcodes\\.tsv\\.gz$", full.names = FALSE)
prefixes <- sub("_barcodes\\.tsv\\.gz$", "", bc_files)   # e.g. "GSM5628163_M03"
message("    ", length(prefixes), " libraries found")

# ----------------------------------------------------------------------------
# [3] Build one Seurat object per library (no filtering) ----
# ----------------------------------------------------------------------------
message("[3] building per-library Seurat objects")
obj_list <- list()

for (prefix in prefixes) {
  gsm <- stringr::str_extract(prefix, "^GSM\\d+")
  lab <- PT_MAP[[gsm]]
  if (is.null(lab) || is.na(lab)) { message("    ! no map for ", gsm, " (", prefix, ") - skipped"); next }

  pt <- parse_label(lab)
  message("    - ", prefix, "  (patient=", pt$patient, ", timepoint=", pt$timepoint, ")")

  counts <- read_mtx_prefixed(RAW, prefix, sep = "_", features_name = "features")
  obj_list[[prefix]] <- make_seurat(
    counts    = counts,
    sample    = prefix,
    dataset   = DATASET,
    patient   = pt$patient,
    timepoint = pt$timepoint,
    study     = STUDY
  )
}

# ----------------------------------------------------------------------------
# [4] Merge + flag pooled libraries ----
# ----------------------------------------------------------------------------
message("[4] merging ", length(obj_list), " libraries")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

seu$Pooled <- grepl("_", seu$Patient_ID)   # multiplexed (multi-donor) libraries

# ----------------------------------------------------------------------------
# [5] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[5] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] GSE185991 finished")
