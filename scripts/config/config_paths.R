# config_paths.R ----
# SINGLE SOURCE OF TRUTH for all filesystem paths (R side).
# ┌──────────────────────────────────────────────────────────────────────────┐
# │ MIRROR FILE: config_paths.sh  --  EDIT BOTH TOGETHER. If you change a root │
# │ here, change the same line in config_paths.sh, or bash/R will disagree.    │
# └──────────────────────────────────────────────────────────────────────────┘
# Sourced (directly or transitively) by every R script. Path RESOLUTION uses the
# `here` package anchored on a `.here` file at the project root, so this file needs
# NO hard-coded path to itself and scripts work from any working directory.
# Output dirs are numbered by the NEW phase scheme (00_ingest .. 09_robustness).

suppressPackageStartupMessages({ library(here) })

## -- project root + config dir, resolved via the `.here` anchor (zero hard-coding) ----
# A file named `.here` must exist at the project root (scripts/ sits under it). here() walks
# up from the running script to find it. PROJECT_ROOT is the code tree root, which may differ
# from the DATA roots (FAST_DIR/LARGE1_DIR) below -- data lives on separate mounts.
PROJECT_ROOT <- here::here()
CONFIG_DIR   <- here::here("scripts", "config")

## -- Roots (the only hard-coded absolute paths in the project) ----
FAST_DIR   <- "/FAST/gr10634/gaozy/aml_niche_net"      # small files, tables, figures (PURGEABLE scratch)
LARGE1_DIR <- "/LARGE1/gr10634/gaozy/aml_niche_net"    # large / persistent objects (BAM, CNV, big RDS)
REF_DIR    <- "/LARGE1/gr10634/gaozy/reference"        # shared references (STAR index, panels, whitelists)
ENV_PREFIX <- "/FAST/gr10634/gaozy/general_env"        # conda env prefix

## -- Global reproducibility ----
SEED <- 491638L

## -- Scripts / config (code tree; may live on a different mount than data) ----
SCRIPTS_DIR <- here::here("scripts")

## -- Big-object roots (LARGE1) ----
RAW_DIR         <- file.path(LARGE1_DIR, "00_raw")                              # raw public deposits
RAW_PUBLIC_DIR  <- file.path(RAW_DIR, "public")
GEO_RAW_DIR     <- file.path(RAW_PUBLIC_DIR, "GEO")                             # 00_ingest: per-dataset GEO deposits
ZENODO_DIR      <- file.path(RAW_PUBLIC_DIR, "Zenodo")                          # 00_ingest: Zenodo + BoneMarrowMap
# Canonical GRCh38-2020-A symbol reference (36601 genes) for make.names() de-sanitising.
REF_FEATURES_2020A <- file.path(GEO_RAW_DIR, "GSE289435_RAW", "GSM8791432_MLL_14666.features.tsv.gz")
RDS_INGEST_DIR  <- file.path(LARGE1_DIR, "01_processed_counts/rds")             # 00_ingest out: <dataset>.rds (merged)
QC_RDS_DIR      <- file.path(LARGE1_DIR, "02_seurat_objects/01_per_sample_qc")  # 01_preprocess out: <ds>/<sample>.rds
PROJ_OBJ_DIR    <- file.path(LARGE1_DIR, "02_seurat_objects/04_bmm_projected")  # 03_hierarchy: projected <ds>/<sample>.rds
CNV_ROOT        <- file.path(LARGE1_DIR, "03_cnv_snv")                          # numbat/cellsnp/vartrix/author/copykat
INFERCNV_ROOT   <- file.path(LARGE1_DIR, "05_cnv_snv/infercnv")                 # inferCNV lives on a separate root (legacy)
INFERCNV_BURDEN_ROOT <- file.path(LARGE1_DIR, "05_cnv_snv/infercnv_burden")
NUMBAT_ROOT     <- file.path(CNV_ROOT, "numbat")                                # <ds>/<sample>/numbat/<sample>__numbat_percell.csv
CNMF_INPUT_DIR  <- file.path(LARGE1_DIR, "06_cnmf/input")                       # 04_cnmf: <ds>/<sample>/{counts.mtx,...}
CNMF_RUN_DIR    <- file.path(LARGE1_DIR, "06_cnmf/runs")                        # cnmf working dir (big scratch)

## -- BoneMarrowMap reference (Zeng 2025) ----
BMM_REF_RDS    <- file.path(ZENODO_DIR, "BoneMarrowMap_SymphonyReference.rds")
BMM_UWOT_MODEL <- file.path(ZENODO_DIR, "BoneMarrowMap_uwot_model.uwot")

## -- Results roots (FAST) : tables + figures, numbered by NEW phase scheme ----
RES_DIR    <- file.path(FAST_DIR, "results")
TAB_DIR    <- file.path(RES_DIR, "tables")
FIG_DIR    <- file.path(RES_DIR, "figures")

DIR_INGEST      <- file.path(TAB_DIR, "00_ingest")       # (migrated from legacy 01_qc)
DIR_PREPROCESS  <- file.path(TAB_DIR, "01_preprocess")   # (migrated from legacy 01_qc mixed outputs)
DIR_MALIGNANCY  <- file.path(TAB_DIR, "02_malignancy")   # (migrated from legacy 03_malignancy) -- c50 consensus, SOLE labeler
DIR_HIERARCHY   <- file.path(TAB_DIR, "03_hierarchy")    # (migrated from legacy 04_hierarchy)
DIR_CNMF        <- file.path(TAB_DIR, "04_cnmf")         # (migrated from legacy 05_cnmf)
DIR_ROBUSTNESS  <- file.path(TAB_DIR, "09_robustness")

## -- stage figure dirs + key cross-stage files ----
FIG_INGEST      <- file.path(FIG_DIR, "00_ingest")
FIG_PREPROCESS  <- file.path(FIG_DIR, "01_preprocess")
QC_MASTER_CSV   <- file.path(DIR_INGEST, "00_MASTER_qc_summary.csv")   # 00_ingest master summary (feeds 01_preprocess)

## -- Frequently used derived table dirs ----
CNMF_PROG_DIR    <- file.path(DIR_CNMF, "programs")                    # <ds>/<sample>__programs.csv
PERCELL_DIR      <- file.path(DIR_HIERARCHY, "percell")                # <ds>/<sample>_percell.csv
PERBIN_DIR       <- file.path(DIR_HIERARCHY, "perbin")                 # <ds>/<sample>_perbin.csv
PERCELL_BINNED_DIR <- file.path(DIR_HIERARCHY, "percell_binned")       # <ds>/<sample>_binned.csv

## -- ensure output table + figure dirs exist (idempotent; FAST is purgeable scratch) ----
# Central creation so any script writing with a bare fwrite() won't fail on a missing dir.
for (.d in c(DIR_INGEST, DIR_PREPROCESS, DIR_MALIGNANCY, DIR_HIERARCHY, DIR_CNMF, DIR_ROBUSTNESS,
             FIG_INGEST, FIG_PREPROCESS)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE, showWarnings = FALSE)
}
rm(.d)
