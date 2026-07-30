# config_hierarchy.R ----
# Stage config for 03_hierarchy (Phase 2). Sources paths, adds projection/bin/stemness constants.
# Extracted from legacy d00_config.R (Phase 2 portion).
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
# PLATFORM_TABLE / platform_of() live in config_qc.R and are the single source of truth for the
# platform stratum. Sourcing is safe: config_qc.R only pulls config_paths.R.
source(here::here("scripts", "config", "config_qc.R"))
suppressPackageStartupMessages({ library(data.table) })

## -- reference tables (live beside the hierarchy scripts) ----
HIER_SCRIPT_DIR  <- file.path(SCRIPTS_DIR, "03_hierarchy")
BIN_MAP_TSV      <- file.path(HIER_SCRIPT_DIR, "bmm_bin_map.tsv")       # broad(24) -> hierarchy_bin(8) + in_ccc_graph
STEMNESS_TSV     <- file.path(HIER_SCRIPT_DIR, "stemness_signatures.tsv")
# REMOVED: DATASET_PLAT_TSV pointed at scripts/03_hierarchy/dataset_platform.tsv, which DOES NOT
# EXIST -- a dangling constant left behind by the refactor out of scripts/以前06_hierarchy/ (the
# only consumers, d15/d20, still live there and carry their own copy of the path). Consequence:
# the per-platform stratification that STRATA_KEY/MAPPING_QC_QUANTILE describe below was silently
# dropped in the current pipeline. Platform now comes from PLATFORM_TABLE (config_qc.R), resolved
# at LIBRARY level because GSE185381 mixes 3'/5' inside one dataset.

## -- BoneMarrowMap projection scaffold (added: needed by 01_bmm_project.R) ----
# Symphony reference + its uwot UMAP model, both from the BMM Zenodo deposit. The scaffold path is
# config_paths.R::BMM_REF_RDS (BoneMarrowMap_SymphonyReference.rds); the uwot model sits beside it.
BMM_SCAFFOLD_RDS <- BMM_REF_RDS
BMM_UWOT         <- file.path(ZENODO_DIR, "BoneMarrowMap_uwot_model.uwot")
BMM_REF_LABEL    <- "CellType_Annotation"   # fine label predict_CellTypes returns (include_broad=TRUE adds broad)
BMM_KNN_K        <- 30L                      # predict_CellTypes k

## -- projection I/O (added) ----
HIER_PROJ_DIR    <- file.path(LARGE1_DIR, "02_seurat_objects", "03_bmm_projected")  # per-sample per-cell projection tables
HIER_TAB_DIR     <- file.path(FAST_DIR, "results", "tables", "03_hierarchy")        # summaries / derived tables
HIER_PROJ_SUMMARY_CSV <- file.path(HIER_TAB_DIR, "bmm_projection_summary.csv")

## -- 8 hierarchy bins (nodes for the future CCC graph) ----
HIERARCHY_BINS <- c("HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid",
                    "Megakaryocyte", "T_NK", "B_Plasma", "Stromal")
# Stromal is annotated but NOT a CCC-graph node (in_ccc_graph = FALSE in bmm_bin_map.tsv).

## -- mapping-QC flag (PRIMARY, scale-free): prob-based ----
# high_error if bmm_celltype_fine_prob < MIN_PROJ_PROB. prob is a 0-1 KNN label-agreement
# proportion, comparable across depth/platform (unlike absolute mapping_error, whose per-platform
# healthy P95 proved unreliable). [DECISION: prob-flag replaced abs-error flag]
MIN_PROJ_PROB <- 0.50

## -- absolute-error thresholds (RETAINED as optional reference column only, NOT the flag) ----
STRATA_KEY          <- "platform"
MAPPING_QC_QUANTILE <- 0.95     # P95 of same-platform healthy cells (reference only)
MAD_THRESHOLD_BMM   <- 2.5      # BMM built-in MAD QC fallback for platforms w/o healthy control

## -- per-bin confidence (drives the derived categorical label) ----
N_MIN_BIN    <- 20L    # bins with fewer cells -> 'low' (binomial fraction too imprecise)
FRAC_ERR_MAX <- 0.50   # bins with > this fraction high_error -> 'low'
TIER2CONF    <- c(A = "high", B = "medium", C = "low")  # evidence_tier -> confidence

## -- healthy datasets anchoring per-platform thresholds ----
HEALTHY_DATASETS <- c("E-MTAB-11536", "GSE253355")

## -- stemness scoring: 4 signatures (from van Galen 2019 Table S3 + LSC17) ----
# vanGalen_HSC_Prog (PRIMARY normal stemness), HSPC_core (cross-val), vanGalen_HSCprog_like
# (tumor LSC readout), LSC17 (weighted bulk control). Loaded from STEMNESS_TSV by the stemness step.

load_bin_map  <- function() {
  if (!file.exists(BIN_MAP_TSV)) stop("[config] bin map not found: ", BIN_MAP_TSV)
  fread(BIN_MAP_TSV, sep = "\t")
}
load_stemness <- function() {
  if (!file.exists(STEMNESS_TSV)) stop("[config] stemness tsv not found: ", STEMNESS_TSV)
  fread(STEMNESS_TSV, sep = "\t")
}
