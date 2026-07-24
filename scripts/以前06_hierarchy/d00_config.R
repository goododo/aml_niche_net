# d00_config.R ----
# Phase 2 (BoneMarrowMap projection + hierarchy bins) shared configuration.
# Source at the top of every d-script:  source("d00_config.R")
# Conventions: ALL_CAPS config constants; make_/load_ helpers; message("[N] ...") logging.

suppressPackageStartupMessages({
  library(data.table)
})

## -- Roots ----
LARGE1_DIR <- "/LARGE1/gr10634/gaozy/aml_niche_net"   # large / persistent objects (BAM, CNV, big RDS)
FAST_DIR   <- "/FAST/gr10634/gaozy/aml_niche_net"     # small files, scripts, tables (PURGEABLE!)

## -- Inputs ----
QC_RDS_DIR <- file.path(LARGE1_DIR, "02_seurat_objects/01_per_sample_qc")  # <dataset>/<sample>.rds
ASSAY      <- "RNA"

# BoneMarrowMap canonical Symphony projection files (downloaded from S3, NOT the Zenodo
# "Annotated_Dataset" Seurat object). map_Query() needs SymphonyReference + uwot model.
BMM_REF_RDS    <- file.path(LARGE1_DIR, "00_raw/public/Zenodo/BoneMarrowMap_SymphonyReference.rds")
BMM_UWOT_MODEL <- file.path(LARGE1_DIR, "00_raw/public/Zenodo/BoneMarrowMap_uwot_model.uwot")
# Fine reference label transferred by predict_CellTypes (55 states). Broad (24) is DERIVED in d20
# from the reference's own meta_data fine->broad pairs (author-defined; do not invent).
BMM_LABEL_COL  <- "CellType_Annotation"

# Optional Phase-1b consensus per-cell labels. Malignant overlay is OPTIONAL: Phase 2 runs on the
# malignancy-INDEPENDENT backbone whether or not these exist yet.
MALIGNANCY_DIR <- file.path(FAST_DIR, "results/tables/03_malignancy")  # <dataset>/...consensus percell csv

# inferCNV per-cell calls from c41 (<ds>/<sample>/<sample>__infercnv_percell.csv: cell,malignant,score,method,sample)
INFERCNV_PERCELL_DIR <- file.path(LARGE1_DIR, "05_cnv_snv", "infercnv")
# Numbat per-cell from c30 (<ds>/<sample>/numbat/<sample>__numbat_percell.csv). TWO schemas:
#   full     : has compartment_opt ("tumor"/"normal") -> usable
#   degraded : cell,sample,p_cnv,note  (note == "no_CNV_detected") -> Numbat found nothing -> NA
NUMBAT_PERCELL_DIR <- file.path(LARGE1_DIR, "03_cnv_snv", "numbat")
# d35 writes a UNIFIED malignant label here, in the SAME schema d30 consumes
# (<sample>__consensus_percell.csv + __consensus_summary.csv), so d30 can run unchanged
# via --malignancy_dir. Label = inferCNV UNION Numbat (each method's negative may be a blind
# spot, so a cell is malignant if EITHER method calls it; agreement -> higher evidence tier).
MALIGNANCY_INFERCNV_DIR <- file.path(FAST_DIR, "results/tables/03_malignancy_infercnv")
# evidence_tier letters drive TIER2CONF (A=high,B=medium,C=low):
#   A_cnv_concordant : inferCNV & Numbat both valid AND agree on this sample's malignant set
#   B_cnv_union      : both methods valid but only partial agreement (union > intersection)
#   C_cnv_single     : only one method valid for this sample (other was blind / not run)
TIER_CNV_CONCORDANT <- "A_cnv_concordant"
TIER_CNV_UNION      <- "B_cnv_union"
TIER_CNV_SINGLE     <- "C_cnv_single"
PILOT_EVIDENCE_TIER <- TIER_CNV_SINGLE   # legacy fallback when only inferCNV exists

## -- Outputs ----
PROJ_OBJ_DIR        <- file.path(LARGE1_DIR, "02_seurat_objects/04_bmm_projected")       # <ds>/<sample>.rds (big)
PERCELL_DIR         <- file.path(FAST_DIR, "results/tables/04_hierarchy/percell")        # <ds>/<sample>_percell.csv
PERBIN_DIR          <- file.path(FAST_DIR, "results/tables/04_hierarchy/perbin")         # <ds>/<sample>_perbin.csv
THRESH_TSV          <- file.path(FAST_DIR, "results/tables/04_hierarchy/mapping_error_thresholds.tsv")
THRESH_OVERRIDE_TSV <- file.path(FAST_DIR, "results/tables/04_hierarchy/mapping_error_thresholds_override.tsv")

## -- Version-controlled lookups (live WITH the scripts) ----
SCRIPT_DIR       <- file.path(FAST_DIR, "scripts/06_hierarchy")
BIN_MAP_TSV      <- file.path(SCRIPT_DIR, "bmm_bin_map.tsv")       # broad(24) -> hierarchy_bin(8) + in_ccc_graph
DATASET_PLAT_TSV <- file.path(SCRIPT_DIR, "dataset_platform.tsv")  # dataset -> platform (strata key for thresholds)

## -- Phase 3 (cNMF meta-programs) ----
# e05 (R) exports per-sample mtx of selected cells (malignant or bin-fallback) with technical genes
# removed; e10 (Python) reads the mtx, lets cNMF pick HVGs, factorizes; e20 (R) clusters programs.
CNMF_INPUT_DIR  <- file.path(LARGE1_DIR, "06_cnmf", "input")    # <ds>/<sample>/{counts.mtx,barcodes.tsv,genes.tsv,meta.json}
CNMF_RUN_DIR    <- file.path(LARGE1_DIR, "06_cnmf", "runs")     # <ds>/<sample>/  (cnmf working dir, big)
CNMF_PROG_DIR   <- file.path(FAST_DIR, "results/tables/05_cnmf/programs")  # <ds>/<sample>__programs.csv (small)
# Input-selection thresholds (DECISION Q2a, three-stage gate):
MIN_MALIGNANT_CNMF <- 200L   # >= -> use union-malignant cells (input_mode="malignant")
MIN_CELLS_CNMF     <- 200L   # else if bin-fallback cells >= this -> input_mode="bin_fallback"; else SKIP
CNMF_FALLBACK_BINS <- c("HSC_MPP", "LMPP_GMP", "Mono_DC")   # myeloid axis bins for the fallback
# cNMF factorization (DECISION Q2b / Q4):
CNMF_K_RANGE  <- 4:9    # keep ALL K programs into the aggregation pool (no single-K selection)
CNMF_N_ITER   <- 50L    # random restarts per K
CNMF_NUM_HVG  <- 2000L  # cNMF picks this many HVGs itself (--numhvgenes)
CNMF_TOP_N    <- 50L    # top-N genes per program for Jaccard (e20)
# Technical gene patterns to REMOVE before cNMF. Conservative set: mito / ribo / hemoglobin, plus
# UNAMBIGUOUS non-coding RNA host transcripts (miRNA, sno/scaRNA, RNU/RN7 small RNAs) -- these are
# essentially always technical noise. We deliberately KEEP lncRNAs (AC*/AL*/LINC*/*-AS) and AP-1/IEG
# (FOS/JUN/EGR...) genes IN the input: some are real leukemia programs, so we let cross-sample SRRS
# reproducibility filter the noise ones rather than hand-removing them (avoids cherry-picking).
# Cell-cycle genes are also KEPT (proliferation may be a real malignant program; flagged at annotation).
TECH_GENE_PATTERNS <- c("^MT[-.]", "^RP[SL]", "^HB[^P]",
                        "^MIR[0-9]", "^SNOR", "^SCARNA", "^RNU[0-9]", "^RN7S")

## -- Tunable knobs (edit here, re-run; no code changes needed) ----
SEED <- 20260605L

# Mapping-QC flag (PRIMARY, scale-free): a cell is high_error if its BMM KNN projection
# confidence bmm_celltype_fine_prob < MIN_PROJ_PROB. prob is a 0-1 proportion (KNN label
# agreement), so it is comparable across depth/platform/dataset -- unlike absolute mapping_error,
# whose per-platform healthy P95 thresholds proved unreliable (healthy anchors too sparse/divergent).
MIN_PROJ_PROB <- 0.50   # prob below this -> high_error (healthy P5 ~= 0.5; AML BMMC tail much lower)

# Absolute-error thresholds (d15) are RETAINED only as an optional reference column, NOT the flag.
STRATA_KEY          <- "platform"  # threshold stratum (shared across healthy & AML datasets)
MAPPING_QC_QUANTILE <- 0.95        # P95 of same-platform healthy cells (reference only)
MAD_THRESHOLD_BMM   <- 2.5         # BMM built-in MAD QC; fallback for platforms w/o healthy control

# Per-bin confidence (DECOMPOSED; these only drive the derived categorical label)
N_MIN_BIN    <- 20L    # bins with fewer cells -> 'low' (binomial fraction too imprecise)
FRAC_ERR_MAX <- 0.50   # bins with > this fraction high_error -> 'low' (bin membership unreliable)
TIER2CONF    <- c(A = "high", B = "medium", C = "low")  # Phase-1b evidence_tier -> confidence

## -- Healthy datasets that ANCHOR the per-platform thresholds ----
# (Authoritative roles live in Phase-1b dataset_roles.tsv; here we only need the healthy set.)
HEALTHY_DATASETS <- c("E-MTAB-11536", "GSE253355")  # extend if more healthy cohorts are added

## -- Helpers ----
# Cross-version counts extraction (Seurat v5 layer / v4 slot).
get_counts <- function(obj, assay = ASSAY) {
  ct <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
                 error = function(e) NULL)
  if (is.null(ct) || nrow(ct) == 0L) {
    ct <- SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
  }
  ct
}

# Safe fwrite: always create the parent dir first (FAST scratch is purgeable).
fwrite_safe <- function(x, file, ...) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file, ...)
}

# Per-sample QC RDS for a dataset -> data.table(dataset, sample, rds).
list_qc_samples <- function(dataset) {
  d <- file.path(QC_RDS_DIR, dataset)
  if (!dir.exists(d)) stop("[config] QC RDS dir not found: ", d)
  f <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  data.table(dataset = dataset, sample = sub("\\.rds$", "", basename(f)), rds = f)
}

load_bin_map <- function() {
  if (!file.exists(BIN_MAP_TSV)) stop("[config] bin map not found: ", BIN_MAP_TSV)
  fread(BIN_MAP_TSV, sep = "\t")
}
