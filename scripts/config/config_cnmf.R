# config_cnmf.R ----
# Stage config for 04_cnmf: per-sample NMF (GeneNMF, Gavish-style) -> cross-sample meta-programs.
# Run SEPARATELY per compartment: malignant vs normal (two parallel program-discovery arms), then
# compare (malignant MPs vs normal MPs) to keep tumor-specific programs and drop normal-like ones.
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
suppressPackageStartupMessages({ library(data.table) })

## -- upstream inputs consumed by 04_cnmf (declared here so this stage doesn't source config_hierarchy) --
# projection per-cell tables (from 03_hierarchy/01_bmm_project.R) -- for hierarchy_bin per cell
HIER_PROJ_DIR   <- file.path(LARGE1_DIR, "02_seurat_objects", "03_bmm_projected")
# consensus per-cell malignant labels come from DIR_MALIGNANCY (defined in config_paths.R)

CNMF_SCRIPT_DIR <- file.path(SCRIPTS_DIR, "04_cnmf")

## -- compartment cell selection (nodes = the 7 hematopoietic CCC bins; Stromal/Unassigned are
##    inferCNV false-positive-prone / non-hematopoietic, so excluded from BOTH arms) ----
CNMF_CCC_BINS   <- c("HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma")
CNMF_MIN_CELLS  <- 100L        # min compartment cells per sample to enter NMF discovery (per-sample NMF needs enough cells)

## -- per-sample NMF (GeneNMF::multiNMF) ----
CNMF_K_RANGE    <- 4:9         # k values per sample (Gavish-style range of program numbers)
CNMF_NFEATURES  <- 2000L       # HVG per sample fed to NMF
CNMF_L1         <- c(0, 0)     # NMF L1 penalty (GeneNMF default; no sparsity constraint)
# keep NMF programs OFF technical / lineage-nuisance axes (pipeline gene patterns)
CNMF_BLOCK_PATTERNS <- c("^MT[-.]", "^RP[SL]", "^HB[^P]")   # mitochondrial / ribosomal / hemoglobin

## -- meta-program aggregation (GeneNMF::getMetaPrograms) ----
CNMF_NMP        <- 10L         # target number of meta-programs per compartment (tune from the tree)
CNMF_METRIC     <- "jaccard"   # program-similarity metric (gene-set overlap; Gavish default)
CNMF_MAX_GENES  <- 200L        # max genes per meta-program signature
CNMF_SPEC_WT    <- 5           # getMetaPrograms specificity.weight
CNMF_MIN_CONF   <- 0.5         # getMetaPrograms min.confidence (drop unstable clusters)

## -- malignant-vs-normal comparison (tumor-specificity) ----
CNMF_MP_SIM_MAX <- 0.20        # a malignant MP with Jaccard > this vs any normal MP is "normal-like" (co-opted/contaminant)

## -- I/O (two arms live under <compartment>/) ----
CNMF_RES_DIR    <- file.path(LARGE1_DIR, "04_cnmf")                        # heavy: nmf_res.rds, per-sample programs
CNMF_TAB_DIR    <- file.path(FAST_DIR, "results", "tables", "04_cnmf")     # meta-program signatures / usage / comparison
CNMF_FIG_DIR    <- file.path(FAST_DIR, "results", "figures", "04_cnmf")

for (d in c(CNMF_RES_DIR, CNMF_TAB_DIR, CNMF_FIG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)