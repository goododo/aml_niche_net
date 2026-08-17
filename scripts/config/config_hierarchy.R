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

## =====================================================================================
## MARKER-BASED ANNOTATION (05) -- the fallback for datasets the BMM projection cannot type
## =====================================================================================
# Used where the projection is untrustworthy. GSE253355 is the case that forced this: its median
# mapping_error runs ~50% above every other dataset in EVERY bin, so its projected bins are not
# usable -- and it supplies 80% of the cohort's stromal cells, so it cannot simply be dropped.
MARKER_TABLE_CSV <- file.path(REF_DIR, "cell_annotation_markers_for_code_260206.csv")
MARKER_ANNO_DIR  <- file.path(LARGE1_DIR, "02_seurat_objects", "03_marker_annotated")
# 06 joins the projection and the marker call per cell and emits ONE corrected bin while keeping
# both inputs as their own columns (bin_bmm / bin_marker). Downstream should read THIS, not either
# source table, so that a result which depends on the typing method is visible rather than latent.
ANNO_RECONCILED_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "03_annotation_reconciled")
# 07 writes Seurat objects carrying both annotations HERE, not back into 01_per_sample_qc: the
# inferCNV freshness guard keys on the QC object's mtime, so re-saving those would mark all 179
# finished runs stale and trigger a full recompute for what is only a metadata edit.
ANNO_OBJ_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "04_annotated")
# OBSOLETE as a bin selector [2026-08-14]. It used to name the datasets whose hierarchy_bin came
# from 05 instead of 01. No dataset does that any more: 08_validate_annotation.R measured the
# marker bins against van Galen 2019 and they lost badly (0.664 vs 0.828; marker overrides wrong
# 5.5x more often than right), so 06 now takes the bin from the projection everywhere and uses the
# markers only for cell_subtype WITHIN that bin. Kept as a non-empty declaration of which datasets
# most depend on marker subtyping -- GSE253355 supplies 80.1% of the cohort's stromal cells, so its
# 17 resolved niche subtypes are the main thing the marker route still buys.
MARKER_ANNO_DATASETS <- c("GSE253355")

MARKER_ANNO_RESOLUTION  <- 1.0    # Louvain resolution; stroma needs enough clusters to resolve
MARKER_ANNO_MIN_CELLS   <- 200L   # below this a per-sample clustering is not meaningful
MARKER_ANNO_MIN_MARKERS <- 2L     # a cell type needs >=2 positive markers PRESENT to be scorable
MARKER_ANNO_NEG_WEIGHT  <- 1.0    # weight on the negative-marker penalty
MARKER_ANNO_MIN_MARGIN  <- 0.25   # top-vs-second score gap below which the call is low_confidence
# "Proliferation state" / "Blast state" describe a STATE, not a lineage: a cycling or blast-like
# cell still belongs to its lineage bin, and letting a state win the argmax would put cells in a
# bin that does not exist in the hierarchy.
MARKER_ANNO_DROP_CATEGORIES <- c("Proliferation state")

# marker-table category -> hierarchy bin. Aligned to bmm_bin_map.tsv so the two typing methods
# produce the SAME vocabulary: MEP -> Erythroid, EoBasoMast Precursor -> LMPP_GMP, Early Lymphoid
# -> T_NK there, so Ery/Meg progenitor, Granulocyte and the T/NK categories follow suit here.
MARKER_CATEGORY_TO_BIN <- c(
  "Stem/Progenitor"         = "HSC_MPP",
  "Blast state"             = "HSC_MPP",   # blast-like states sit in the progenitor compartment
  "Lymphoid progenitor"     = "LMPP_GMP",
  "Granulocyte"             = "LMPP_GMP",  # cf. EoBasoMast Precursor -> LMPP_GMP in the BMM map
  "Monocyte/Macrophage"     = "Mono_DC",
  "Dendritic cell"          = "Mono_DC",
  "Myeloid (general)"       = "Mono_DC",
  "Erythroid"               = "Erythroid",
  "Ery/Meg progenitor"      = "Erythroid", # cf. MEP -> Erythroid in the BMM map
  "Megakaryocyte/Platelet"  = "Megakaryocyte",
  "T cell (general)"        = "T_NK",
  "T cell (CD4)"            = "T_NK",
  "T cell (CD8)"            = "T_NK",
  "T helper"                = "T_NK",
  "Innate-like T"           = "T_NK",
  "NK"                      = "T_NK",
  "B cell"                  = "B_Plasma",
  "Plasma cell"             = "B_Plasma",
  "Bone marrow stroma"      = "Stromal",
  "Skeletal lineage"        = "Stromal",
  "Adipo lineage"           = "Stromal",
  "Bone marrow vasculature" = "Stromal"    # BMM carries a single "Stromal" broad type; endothelium
)                                          # folds into it so the two vocabularies stay comparable

# CELL-TYPE overrides, applied ON TOP of the category map. Several categories are heterogeneous in
# lineage, so mapping them wholesale puts cells in the wrong bin and manufactures a disagreement
# with the projection that is entirely our own doing. Measured: with the category map alone,
# LMPP_GMP agreement was 0.00-0.17 in EVERY dataset -- a cohort-wide constant, i.e. a vocabulary
# bug, not a data property. Each line below is anchored to bmm_bin_map.tsv so the two methods
# partition the same way:
#   * "Blast state" names its own lineage per cell type ("Blast-like: GMP" is a GMP), so it cannot
#     collapse to one progenitor bin.
#   * BMM maps LMPP -> LMPP_GMP and Early Lymphoid -> T_NK, so the HSPC types must follow.
#   * The marker table files a "T/NK" type under category "Stem/Progenitor"; taking the category
#     there would bin T cells as stem cells.
#   * The 8-bin hierarchy has no granulocyte node. BMM sends EoBasoMast Precursor -> LMPP_GMP, so
#     baso/eo/mast follow it; mature neutrophils are not progenitors and BMM assigns them Mono_DC.
MARKER_CELLTYPE_TO_BIN <- c(
  "Blast-like: HSC/MPP"             = "HSC_MPP",
  "Blast-like: GMP"                 = "LMPP_GMP",
  "Blast-like: Pro-monocyte"        = "Mono_DC",
  "Blast-like: Monocyte"            = "Mono_DC",
  "Blast-like: cDC"                 = "Mono_DC",
  "HSPC: HSC/MPP"                   = "HSC_MPP",
  "HSPC: LMPP"                      = "LMPP_GMP",
  "HSPC: CLP/Early lymphoid"        = "T_NK",
  "T/NK"                            = "T_NK",
  "Granulocyte: Neutrophil"         = "Mono_DC",
  "Granulocyte: Neutrophil (CST7+)" = "Mono_DC",
  "Granulocyte: Basophil"           = "LMPP_GMP",
  "Granulocyte: Eosinophil"         = "LMPP_GMP",
  "Granulocyte: Mast"               = "LMPP_GMP",
  # Genuinely bi-lineage by name (MEP is erythroid-megakaryocytic). The category map already sent
  # it to Erythroid, matching BMM's MEP -> Erythroid, but leaving that implicit means the choice
  # is invisible; stated here so it is a decision on the record rather than a side effect.
  "Ery/Meg prog: MEP/MkP"           = "Erythroid"
)

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

## =====================================================================================
## CELL-STATE FEATURE PANEL (09) -- functional axes the node features currently lack
## =====================================================================================
# The FGW feature term today is frac_malignant + mean_stemness + n_cells, and the measured
# decomposition (08_scoring/feature_decomposition.csv) says mean_stemness carries NO global signal
# (beta -0.044, p 0.396) while adding it DILUTES the combined estimate (all3 beta 0.135 vs
# frac_malignant alone 0.334). So the panel below is built to be tested one axis at a time, and
# nothing enters FGW_FEATURES until it is individually significant. Most of it is for mechanism
# and for validating CCC edges, NOT for graph alignment.
CELLSTATE_TABLE_TSV <- file.path(HIER_SCRIPT_DIR, "cellstate_signatures.tsv")
CELLSTATE_DIR       <- HIER_PROJ_DIR   # per-cell scores sit beside the projection they are keyed to

# UCELL, NOT AddModuleScore. UCell ranks genes WITHIN each cell, so a score depends on nothing but
# that cell -- which is what makes it comparable across 214 samples. AddModuleScore subtracts the
# mean of control-gene bins drawn from whatever object it is handed, and 04_stemness_score.R calls
# it PER SAMPLE, so part of the between-sample signal is removed by construction: measured
# between/within sd ratio 0.300 for HSPC_core vs 0.689 for LSC17, which is a raw weighted sum and
# is not recentred. (The node feature in use is LSC17, so this is not the explanation for the null
# stemness result -- but it is a reason not to build NEW cross-sample features that way.)
CELLSTATE_MAXRANK   <- 1500L   # UCell rank horizon; genes below this rank contribute nothing
CELLSTATE_NCORES    <- 8L      # measured 646.6s (1 core) vs 38.5s (8) on an 8,469-cell sample

# Single genes reported as mean log-normalised expression rather than a UCell score. A rank-based
# score over one gene is meaningless, but the MCL1-vs-BCL2 balance is exactly the quantity the
# venetoclax-resistance literature is about, so it is emitted directly.
CELLSTATE_SINGLE_GENES <- c("BCL2", "MCL1", "BCL2L1", "CXCR4", "KIT", "CD36", "HLA-DRA", "CD74")

# POSITIVE CONTROLS, AS ORDERED CONTRASTS. A signature with a stale gene symbol or cross-axis
# bleed still produces a perfectly ordinary-looking number for every cell, so each panel has to be
# checked against biology that is not in question.
#
# These were argmax tests ("this signature must PEAK in bin X") and that was the wrong instrument.
# mhc_ii failed it by peaking in HSC_MPP -- which is CORRECT biology, not a bug: human CD34+
# progenitors are largely HLA-DR positive (CD34+CD38-HLA-DR- is how you deplete them to enrich
# HSC). Argmax asks which compartment is highest across the whole marrow, which is often a genuine
# empirical question; a pairwise contrast asks something that is actually settled.
#
# Each row is high_bin > low_bin, and the pair is chosen so that a violation can only mean the
# panel is broken. T cells do not express MHC-II; stroma is not cytotoxic; T cells are not stromal.
CELLSTATE_CONTRASTS <- data.table::data.table(
  signature = c("cytotoxicity", "cytotoxicity", "exhaustion", "naive_memory",
                "mhc_ii", "mhc_ii", "retention_ligand", "retention_ligand", "sasp"),
  high_bin  = c("T_NK", "T_NK", "T_NK", "T_NK",
                "B_Plasma", "Mono_DC", "Stromal", "Stromal", "Stromal"),
  low_bin   = c("Erythroid", "Stromal", "B_Plasma", "Mono_DC",
                "T_NK", "T_NK", "T_NK", "B_Plasma", "T_NK"),
  why       = c("erythroid cells carry no cytotoxic granules",
                "stroma is not a cytotoxic lineage",
                "PDCD1/LAG3/TOX are T-cell exhaustion markers, not B-cell",
                "CCR7/TCF7/LEF1 are T-naive markers",
                "T cells do not express MHC-II; B cells do",
                "T cells do not express MHC-II; monocytes/DC do",
                "CXCL12/KITLG/VCAM1 are stromal niche ligands",
                "CXCL12/KITLG/VCAM1 are stromal niche ligands",
                "SASP is a stromal/myeloid inflammatory programme")
)

## -- PROGENy pathway footprints (10) ----
# The model ships INSIDE the progeny package, so nothing here needs the network at run time. top=500
# is the PROGENy default and the setting its benchmarks were run at.
PROGENY_TOP          <- 500L
PROGENY_MINSIZE      <- 5L      # a pathway needs >=5 of its genes present to be scored at all
PROGENY_MIN_OVERLAP  <- 2000L   # gene-symbol sanity: a near-empty overlap still returns numbers
# Interferon signals through JAK-STAT, so PROGENy's JAK-STAT footprint and the interferon_I/II
# UCell panels in 09 must agree despite sharing no genes and no method. Set from the observed
# separation against the unrelated-pathway baseline, not from a rule of thumb.
PROGENY_MIN_IFN_COR  <- 0.15

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
