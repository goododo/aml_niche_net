# config_qc.R ----
# Stage config for 00_ingest + 01_preprocess (Phase 1A). Sources paths, adds QC constants.
# Style: ALL_CAPS constants (was cfg$list in legacy 02_preprocess/00_config.R -- now unified).

suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
suppressPackageStartupMessages({ library(data.table) })

## -- metadata schema (uniform across all merged RDS) ----
COL_SAMPLE  <- "Sample";     COL_DATASET <- "Dataset";  COL_STUDY   <- "Study"
COL_PATIENT <- "Patient_ID"; COL_TP      <- "Timepoint"; COL_DISEASE <- "Disease_state"
COL_NCOUNT  <- "nCount_RNA"; COL_NFEAT   <- "nFeature_RNA"
COL_MT      <- "percent.mt"; COL_HB      <- "percent.hb"; COL_RIBO <- "percent.ribo"
COL_DEMUXED <- "Demuxed"

## -- gene patterns (also the ingest QC metrics) ----
PAT_MITO <- "^MT[-.]"       # + MITO_SHORT fallback for prefix-stripped datasets
PAT_RIBO <- "^RP[SL]"
PAT_HB   <- "^HB[^P]"
MITO_SHORT <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6","CYTB","COX1","COX2","COX3","ATP6","ATP8")

## -- CANONICAL TIMEPOINT VOCABULARY (C1: unify Timepoint across datasets for H3) ----
# Every ingest_<dataset>.R sets its RAW timepoint as before, then calls canonicalize_timepoint()
# to map (dataset, raw_timepoint, disease_state) -> ONE canonical value. Raw value is preserved
# in a new Timepoint_detail column, so nothing is lost.
CANONICAL_TIMEPOINTS <- c("Healthy", "Diagnosis", "MRD", "Post_treatment", "Relapse", "Relapse2", "Unknown")

# WHY a table (not a generic regex): the string "Baseline" means HEALTHY in E-MTAB/GSE253355/
# GSE116256 (healthy donors) but AML DIAGNOSIS in Chen2023/Petti2019. Mapping must be dataset-aware.
# Resolution order in canonicalize_timepoint():
#   0. sample-level healthy donor    -> Healthy (NBM*/ND_*/Normal_*/BM<n> inside an AML dataset)
#   1. whole-healthy datasets       -> Healthy (ignore raw string)
#   2. disease_state-driven datasets -> from Disease_state (Petti, Lasry)
#   3. GSE116256 day-codes           -> D0=Diagnosis, D<n>=Post_treatment, Baseline=Healthy
#   4. explicit per-dataset raw->canonical map (TIMEPOINT_MAP)
#   5. fallback                      -> Unknown (with a warning)

# datasets whose every cell is a healthy donor (raw timepoint irrelevant)
TP_HEALTHY_DATASETS <- c("E-MTAB-11536", "GSE253355")
# datasets where the canonical timepoint comes from Disease_state, not Timepoint
TP_DISEASE_DRIVEN   <- c("Petti2019", "GSE185381")

# SAMPLE-level healthy donors inside an otherwise-AML dataset. The dataset-level rules above cannot
# see these, so they silently inherited "Diagnosis" (and the dataset's subtype_stratum): Chen2023
# NBM1-4, Petti2019 ND_*/Normal_sorted_*, GSE116256 BM4/BM5-*. That put healthy marrow into the AML
# arm of every downstream comparison. Rule 0 in canonicalize_timepoint() applies this FIRST.
# Patterns are anchored and deliberately narrow; every hit is logged so the list stays auditable.
# Cross-check with a name-independent signal before trusting a new dataset here.
TP_HEALTHY_SAMPLE_PATTERNS <- c(
  "^NBM[0-9_-]",        # Chen2023   normal bone marrow
  "^ND[_-][0-9]",       # Petti2019  normal donor (ND_083017)
  "^Normal[_-]",        # Petti2019  Normal_sorted_*
  "^BM[0-9]+([_-]|$)",  # GSE116256  van Galen healthy BM (BM4, BM5-34p, BM5-34p38n)
  "^HD[_-]?[0-9]"       # generic    healthy donor
)
TP_HEALTHY_SAMPLE_REGEX <- paste(TP_HEALTHY_SAMPLE_PATTERNS, collapse = "|")

is_healthy_sample <- function(sample_id, patient_id = NA_character_) {
  x <- c(as.character(sample_id), as.character(patient_id))
  x <- x[!is.na(x) & nzchar(x)]
  length(x) > 0 && any(grepl(TP_HEALTHY_SAMPLE_REGEX, x, ignore.case = TRUE))
}

# explicit raw -> canonical, keyed per dataset (only datasets needing a lookup)
TIMEPOINT_MAP <- data.table::fread(text = '
dataset|raw|canonical
Chen2023|Baseline|Diagnosis
Chen2023|Diagnosis|Diagnosis
GSE239721|Diagnosis|Diagnosis
GSE289435|Diagnosis|Diagnosis
GSE227903|Diagnosis|Diagnosis
GSE227903|MRD|MRD
GSE227903|Relapse|Relapse
GSE227903|Relapse2|Relapse2
GSE185991|DX|Diagnosis
GSE185991|D14|Post_treatment
GSE185991|D30|Post_treatment
GSE185991|REL|Relapse
GSE147989|SCR|Diagnosis
GSE147989|C1D1|Post_treatment
GSE201966|Primary|Diagnosis
GSE201966|Complete_remission|Post_treatment
GSE201966|Relapse|Relapse
GSE207356|Sample_D|Diagnosis
GSE207356|Sample_E|Post_treatment
GSE207356|Sample_F|Post_treatment
', sep = "|")

# Disease_state -> canonical (for TP_DISEASE_DRIVEN datasets). "Mixed" = un-demuxed composite
# library (removed by DEMUX_PREFILTER); map to Unknown but preserve "Mixed" in Timepoint_detail.
DISEASE_TO_TP <- c(Healthy = "Healthy", AML = "Diagnosis", Mixed = "Unknown")

# Map a single sample's raw timepoint to canonical. Returns canonical string; caller keeps raw.
# sample_id/patient_id are optional so older call sites keep working, but SHOULD be passed: without
# them rule 0 cannot fire and healthy donors inside AML datasets are mislabelled as Diagnosis.
canonicalize_timepoint <- function(ds_id, raw_tp = NA_character_, dis = NA_character_,
                                   sample_id = NA_character_, patient_id = NA_character_) {
  ds_id <- as.character(ds_id); raw_tp <- as.character(raw_tp); dis <- as.character(dis)
  # 0. sample-level healthy donor inside an AML dataset (highest priority: a normal marrow is
  #    normal regardless of what the dataset-level rule would say). Logged, because overriding a
  #    dataset rule from a NAME is exactly the kind of decision that must stay visible.
  if (is_healthy_sample(sample_id, patient_id)) {
    would_be <- if (ds_id %in% TP_HEALTHY_DATASETS) "Healthy" else NA_character_
    if (!identical(would_be, "Healthy"))
      message(sprintf("[tp] rule0 healthy sample: %s::%s (raw='%s', disease='%s') -> Healthy",
                      ds_id, sample_id, raw_tp, dis))
    return("Healthy")
  }
  # 1. whole-healthy datasets
  if (ds_id %in% TP_HEALTHY_DATASETS) return("Healthy")
  # 2. disease-state-driven
  if (ds_id %in% TP_DISEASE_DRIVEN) {
    d <- DISEASE_TO_TP[dis]
    return(if (!is.na(d)) unname(d) else "Unknown")
  }
  # 3. GSE116256 day-codes
  if (ds_id == "GSE116256") {
    if (grepl("^Baseline$", raw_tp, ignore.case = TRUE)) return("Healthy")   # van Galen healthy BM
    if (grepl("^D0$", raw_tp))                           return("Diagnosis")
    if (grepl("^D[0-9]+$", raw_tp))                      return("Post_treatment")
    return("Unknown")
  }
  # 4. explicit table (plain vector match; no data.table subset to avoid name/scoping issues)
  hit <- TIMEPOINT_MAP$canonical[TIMEPOINT_MAP$dataset == ds_id & TIMEPOINT_MAP$raw == raw_tp]
  if (length(hit) == 1L) return(hit)
  # 5. fallback
  warning(sprintf("canonicalize_timepoint: no mapping for dataset=%s raw=%s ds=%s -> Unknown",
                  ds_id, raw_tp, dis))
  "Unknown"
}

## -- study-level role table (curated; drives split + L2 capability) ----
ROLE_TABLE <- data.table::fread(text = '
dataset|study_role|is_longitudinal|discovery_candidate|subtype_stratum|l2_capable|integrity_flag
Chen2023|Discovery-AML+AuxStroma|FALSE|TRUE|NPM1|TRUE|OK
GSE239721|Discovery-AML|FALSE|TRUE|other|TRUE|OK
GSE289435|Discovery-AML|FALSE|TRUE|other|TRUE|OK
Petti2019|Discovery-AML|FALSE|TRUE|other|TRUE|OK
GSE185381|Discovery-AML-immuneReceiver|FALSE|TRUE|other|TRUE|OK
E-MTAB-11536|Healthy-control|FALSE|FALSE|NA|TRUE|OK
GSE253355|Reference-scaffold+AuxStroma|FALSE|FALSE|NA|FALSE|OK
GSE116256|Treatment-axis|TRUE|FALSE|mixed|TRUE|OK
GSE185991|Treatment-axis-L1only|TRUE|FALSE|NPM1|FALSE|OK
GSE147989|Treatment-axis-L1only|TRUE|FALSE|mixed|FALSE|OK
GSE227903|Validation-relapse+treatment|TRUE|FALSE|mixed|TRUE|OK
GSE201966|Validation-relapse|TRUE|FALSE|monocytic|TRUE|OK
GSE207356|Exploratory|TRUE|FALSE|NA|TRUE|OK
', sep = "|")

## -- sample-level role refinement regex ----
HEALTHY_REGEX  <- "healthy|normal|control|donor"
CELLLINE_REGEX <- "cell line|cell-line|cellline"
PDX_REGEX      <- "pdx|xenograft"

## -- QC thresholds: per-sample data-driven MAD (scater) ----
MAD_NMADS        <- 3L        # 3 * MAD
MAD_LOG_COUNTS   <- TRUE      # MAD on log10(nCount_RNA), both tails
MAD_LOG_FEATS    <- TRUE      # MAD on log10(nFeature_RNA), both tails
MT_TAIL          <- "higher"  # percent.mt outliers: upper tail only
MAD_HB_ENABLE    <- TRUE
MAD_HB_TAIL      <- "higher"

## -- ABSOLUTE bounds, applied IN ADDITION to the per-sample MAD ----
# MAD alone is purely relative, so the worse a sample is the looser its own thresholds become:
# observed nfeat_lo as low as 0.08 genes and percent.mt cutoffs computed above 100%. Filter
# strength ended up inversely proportional to sample quality. These floors/caps are the fixed
# biological sanity limits MAD cannot provide; a cell must clear BOTH to be kept.
ABS_MIN_NFEAT <- 200L    # below this a barcode is debris/empty, whatever the sample median says
ABS_MIN_NCOUNT <- 500L
ABS_MAX_MT    <- 20      # percent.mt hard cap (MAD-derived caps exceeded 100% on 2 samples)
MIN_CELLS_SAMPLE <- 300L      # sample gate AFTER cell QC + doublet removal
# Lowered from 500. Cell count alone was the ONLY exclusion criterion, which kept large low-quality
# samples (GSE201966 AML01M4_P: 6000 cells at 182 genes/cell) and dropped small high-quality ones
# (GSE227903 3822_R: 462 cells at 1650 genes/cell -- 38 cells short of the line). Quality is now
# enforced by FLAG_LOWCOMPLEXITY_DROP above, so this gate only needs to guarantee enough cells to
# populate the hierarchy bins. NOTE: the binding constraint for L2 is per-bin occupancy (>=3 bins
# with >=30 cells), which can only be evaluated after 03_hierarchy -- apply it there, not here.
DEMUX_PREFILTER  <- TRUE      # keep Demuxed==TRUE cells before per-sample QC (Lasry)

## -- low-complexity gate (was a flag-only soft check; now an actual filter) ----
# 250 genes/cell is far too low for CCC: a ligand and its receptor must be co-detected in the same
# cell, so libraries at ~200 genes/cell yield noise, not communication. Raised to 500 and turned ON.
# Consequence: samples like GSE201966 AML01M4_P (6000 cells at 182 genes/cell) are now excluded,
# while good small samples are retained by the complexity-aware gate below.
FLAG_MIN_MED_NCOUNT     <- 500
FLAG_MIN_MED_NFEAT      <- 500
FLAG_LOWCOMPLEXITY_DROP <- TRUE
INTEGRITY_MIN_MED_NCOUNT <- 100   # hard integrity floor -> flag + skip doublet

## -- doublet detection: dual-R consensus (scDblFinder + DoubletFinder) ----
# UNION, not intersection. Requiring BOTH callers recovered only ~37% of the expected doublets
# (15,488 removed of 41,573 expected; 54 samples with >=100 cells had none removed at all). A
# surviving blast+T-cell doublet expresses a ligand and its receptor in one barcode, which CellChat
# reads as a genuine communication edge -- the exact artefact this project cannot afford.
# Union is less conservative on real cells; that trade is deliberate.
DOUBLET_CONSENSUS <- "union"
DOUBLET_MIN_CELLS <- 100L
DBL_RATE_PER_1K   <- 0.008           # ~0.8% per 1000 cells (10x)
DBL_RATE_CAP      <- 0.10

expected_dbl_rate <- function(n_cells) min(DBL_RATE_PER_1K * (n_cells / 1000), DBL_RATE_CAP)