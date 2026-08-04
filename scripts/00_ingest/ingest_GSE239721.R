# ============================================================================
# ingest_GSE239721.R   --   IFNgamma-AML single-cell (already a Seurat object)
# Format : a single Seurat .RDS  "GSE239721_aml_seurat.RDS"
# Action : load, update legacy object, derive Patient_ID/Timepoint from the
#          existing sample column, recompute QC with the shared rules,
#          preserve original meta as "anno_<col>". No merge (single object).
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE239721"
STUDY   <- "GSE239721 (IFN-gamma AML)"
RDS_IN  <- file.path(GEO_RAW_DIR, "GSE239721", "GSE239721_aml_seurat.RDS")

# ----------------------------------------------------------------------------
# [1] Load and update the Seurat object ----
# ----------------------------------------------------------------------------
message("[1] loading ", basename(RDS_IN))
seu <- readRDS(RDS_IN)
seu <- UpdateSeuratObject(seu)   # harmonise to the installed Seurat version
message("    ", ncol(seu), " cells x ", nrow(seu), " genes")
message("    existing meta columns: ", paste(colnames(seu@meta.data), collapse = ", "))

# ----------------------------------------------------------------------------
# [2] Preserve original metadata as anno_<col> ----
# ----------------------------------------------------------------------------
orig_meta <- seu@meta.data
keep_as_anno <- setdiff(colnames(orig_meta),
                        c("nCount_RNA", "nFeature_RNA"))   # QC recomputed below
for (cn in keep_as_anno) seu[[paste0("anno_", cn)]] <- orig_meta[[cn]]

# ----------------------------------------------------------------------------
# [3] Derive Patient_ID / Timepoint from an existing sample column ----
# ----------------------------------------------------------------------------
# Probe likely sample columns; first hit wins. Cross-sectional -> Diagnosis.
SAMPLE_COL_CANDIDATES <- c("sample","Sample","sample_id","sampleID","patient",
                           "patient_id","donor","orig.ident")
hit <- intersect(SAMPLE_COL_CANDIDATES, colnames(orig_meta))
sample_col <- if (length(hit)) hit[1] else "orig.ident"
message("[3] using '", sample_col, "' as Patient_ID source")

pid <- as.character(orig_meta[[sample_col]])

seu$Dataset    <- DATASET
seu$Study      <- STUDY
seu$Sample     <- pid
seu$Patient_ID <- pid
seu$Timepoint  <- "Diagnosis"
seu$orig.ident <- pid

# This ingest builds its metadata by hand (the input is an author Seurat object, not a count
# matrix), so it never calls make_seurat() and therefore never got make_seurat's standard
# columns. That is invisible until something downstream filters on one of them -- Tissue came
# through blank for all 20 samples. Set the same defaults here, explicitly.
seu$uid_patient         <- paste0(DATASET, ":", pid)
seu$Disease_state       <- "AML"
seu$Timepoint_detail    <- "Diagnosis"
seu$Timepoint_days      <- NA_integer_
.plat                   <- platform_of(DATASET)
seu$Platform            <- .plat$platform
seu$Chemistry           <- .plat$chemistry
seu$Tissue              <- TISSUE_DEFAULT     # all 20 curated rows are BM
seu$N_donors_in_library <- 1L
seu$Patient_resolved    <- TRUE

# ----------------------------------------------------------------------------
# [4] Recompute QC metrics with the shared rules ----
# ----------------------------------------------------------------------------
message("[4] recomputing QC metrics")
DefaultAssay(seu) <- "RNA"
seu <- add_qc_metrics(seu)

# ----------------------------------------------------------------------------
# [5] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[5] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] GSE239721 finished")
