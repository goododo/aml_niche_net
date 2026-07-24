#!/usr/bin/env bash
# config_paths.sh ----
# SINGLE SOURCE OF TRUTH for all filesystem paths (bash side).
# ┌──────────────────────────────────────────────────────────────────────────┐
# │ MIRROR FILE: config_paths.R  --  EDIT BOTH TOGETHER. If you change a root  │
# │ here, change the same line in config_paths.R, or bash/R will disagree.     │
# └──────────────────────────────────────────────────────────────────────────┘
# Source from every cluster (.sh/.slurm) step:  source "<...>/config/config_paths.sh"

## -- Roots (the only hard-coded absolute paths in the project) ----
export FAST_DIR="/FAST/gr10634/gaozy/aml_niche_net"      # small files, tables, figures (PURGEABLE scratch)
export LARGE1_DIR="/LARGE1/gr10634/gaozy/aml_niche_net"  # large / persistent objects (BAM, CNV, big RDS)
export REF_DIR="/LARGE1/gr10634/gaozy/reference"         # shared references
export ENV_PREFIX="/FAST/gr10634/gaozy/general_env"      # conda env prefix

## -- Global reproducibility ----
export SEED=491638

## -- Big-object roots (LARGE1) ----
export QC_RDS_DIR="${LARGE1_DIR}/02_seurat_objects/01_per_sample_qc"
export CNV_ROOT="${LARGE1_DIR}/03_cnv_snv"              # cellsnp/numbat/vartrix/author/copykat
export INFERCNV_ROOT="${LARGE1_DIR}/05_cnv_snv/infercnv"
export WORK_ROOT="${CNV_ROOT}/_work"                    # per-sample scratch (transient)

## -- Results (FAST), NEW phase numbering ----
export MAL_TAB="${FAST_DIR}/results/tables/02_malignancy"   # c50 consensus output
export MAL_FIG="${FAST_DIR}/results/figures/02_malignancy"

## ===================== USER-CONFIRM (3 items, per new dataset) =====================
# See docs/NEW_DATASET_RUNBOOK.md for the decision rationale on each.
#
# (1) STAR index. hg38/star_index and refdata-gex-GRCh38-2024-A/star both exist; neither is
#     the GRCh38-2020-A used by QC. Numbat CNV is position-based (annotation-tolerant).
export STAR_INDEX="${REF_DIR}/hg38/star_index"
#
# (2) SNP VCF for cellsnp -- contig naming MUST match the STARsolo BAM (chr vs non-chr).
export NB_SNP_VCF="${REF_DIR}/hg38/vcf/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
# export NB_SNP_VCF="${REF_DIR}/hg38/vcf/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.CHR.sorted.vcf.gz"
#
# (3) STAR FIFO tmp -- Lustre (/LARGE1) may reject FIFOs; default node-local /tmp.
export STAR_TMP_BASE="/tmp"
## ==================================================================================

## -- 10x whitelists ----
export WHITELIST_V3="${REF_DIR}/10x/whitelist/3M-february-2018.txt"
export WHITELIST_V2="${REF_DIR}/10x/whitelist/737K-august-2016.txt"

## -- Numbat references ----
export NB_PANEL_DIR="${REF_DIR}/1000G_hg38/1000G_hg38"
export NB_EAGLE_BIN="${REF_DIR}/Eagle_v2.4.1/eagle"
export NB_GMAP="${REF_DIR}/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
export NB_GENOME_VERSION="hg38"
export NB_REF_RDS="${REF_DIR}/numbat/numbat_ref_hg38.rds"

## -- Compute / behavior ----
export THREADS=16
export KEEP_BAM="false"
export CLEAN_FASTQ="true"

mkdir -p "${WORK_ROOT}" "${MAL_TAB}" "${MAL_FIG}" "${CNV_ROOT}/numbat" "${CNV_ROOT}/cellsnp"
echo "[cfg] ENV=${ENV_PREFIX}  STAR_INDEX=${STAR_INDEX}  THREADS=${THREADS}  SEED=${SEED}"
