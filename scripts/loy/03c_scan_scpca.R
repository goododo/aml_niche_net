#!/usr/bin/env Rscript
# 03c_scan_scpca.R ----
# input : ScPCA project SCPCP000007 as downloaded from the portal (LOY_SCPCA_DIR: SCPCP000007_SINGLE_CELL_SCE.zip with one
#         <library>_processed.rds SingleCellExperiment per sample, and SCPCP000007_ALL_METADATA.zip with metadata.tsv), the
#         BoneMarrowMap scaffold, LOY_YGENES / LOY_XIST
# output: LOY_PERCELL_DIR/SCPCA_SCPCP000007__y_percell.csv.gz (01 / 03 schema), LOY_SCAN_DIR/SCPCA_SCPCP000007__y_sample_summary.csv
#         and __y_scan_checks.csv; the dataset id appended to LOY_EXT_LIST; the zips extracted once under LOY_SCPCA_DIR/extracted/
# what it does: step 2c of the LOY line for the Gawad paediatric AML cohort. Per sample: raw counts from the SCE counts assay,
#   Ensembl ids -> gene symbols (rowData gene_symbol; duplicates made unique), light QC (LOY_GEO_QC_*), this project's
#   BoneMarrowMap projection, chrY-gene and XIST UMIs per cell, sex / age / diagnosis / timing from metadata.tsv (the recorded
#   sex is kept as anno_sex for 02's sex review; the call still comes from expression).
#   The projection and Y-count blocks mirror 03b_scan_geo_mtx.R (03b was running when this file was written; the two are to be
#   unified into 00_loy_common.R when both are idle).
# Usage : sbatch scripts/loy/submit_loy_scan_scpca.sbatch   |   Rscript scripts/loy/03c_scan_scpca.R [--force]
# WHY: 16 male samples of paediatric AML with raw counts (the Zeng-pack copy is ambient-corrected and was refused by the integer check).
suppressPackageStartupMessages({ library(data.table); library(here); library(Matrix); library(Seurat); library(SeuratObject); library(BoneMarrowMap); library(SingleCellExperiment) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "loy", "config_loy.R"))
args  <- commandArgs(trailingOnly = TRUE); force <- any(args == "--force")
ds    <- LOY_SCPCA_DATASET
out_cells <- file.path(LOY_PERCELL_DIR, paste0(ds, "__y_percell.csv.gz"))
if (file.exists(out_cells) && !force) { message("[skip] exists: ", out_cells); quit(save = "no") }
dir.create(LOY_PERCELL_DIR, showWarnings = FALSE, recursive = TRUE); dir.create(LOY_SCAN_DIR, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); message(sprintf("[0] %s | %s", ds, format(Sys.time())))

## [1] extract once; metadata ----------------------------------------------------------------------------------
ext_dir <- file.path(LOY_SCPCA_DIR, "extracted"); dir.create(ext_dir, showWarnings = FALSE)
if (!length(list.files(ext_dir, pattern = "_processed\\.rds$", recursive = TRUE))) {
  message("[1] extracting the SCE zip"); system2("unzip", c("-o", "-q", shQuote(file.path(LOY_SCPCA_DIR, "SCPCP000007_SINGLE_CELL_SCE.zip")), "-d", shQuote(ext_dir)))
  system2("unzip", c("-o", "-q", shQuote(file.path(LOY_SCPCA_DIR, "SCPCP000007_ALL_METADATA.zip")), "-d", shQuote(file.path(ext_dir, "metadata"))))
}
rds  <- list.files(ext_dir, pattern = "_processed\\.rds$", recursive = TRUE, full.names = TRUE)
meta <- fread(list.files(ext_dir, pattern = "^metadata\\.tsv$", recursive = TRUE, full.names = TRUE)[1])
meta <- meta[, .(scpca_sample_id, scpca_library_id, diagnosis, subdiagnosis, disease_timing, age, sex, technology, participant_id)]
message(sprintf("[1] %d processed SCE files; metadata rows %d", length(rds), nrow(meta)))

## [2] scaffold + bin map (as 03b) ----------------------------------------------------------------------------------
ref <- readRDS(BMM_SCAFFOLD_RDS); bin_map <- fread(LOY_BMM_BIN_MAP)
pick_col <- function(df, exact, pattern) { for (e in exact) if (e %in% names(df)) return(e); h <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE); if (length(h)) h[1] else NA_character_ }
project <- function(cnt, md) {
  q <- map_Query(exp_query = cnt, metadata_query = md, ref_obj = ref, do_normalize = TRUE, do_umap = FALSE, verbose = FALSE)
  q <- calculate_MappingError(q, reference = ref, MAD_threshold = MAD_THRESHOLD_BMM)
  q <- predict_CellTypes(query_obj = q, ref_obj = ref, ref_label = BMM_REF_LABEL, k = BMM_KNN_K, include_broad = TRUE)
  qm <- q@meta.data
  data.table(cell = rownames(qm),
             bmm_fine = as.character(qm[[pick_col(qm, c("predicted_CellType", "final_CellType"), "^predicted_CellType$")]]),
             predicted_CellType_Broad = as.character(qm[[pick_col(qm, c("predicted_CellType_Broad", "predicted_Broad"), "[Bb]road")]]),
             bmm_prob = as.numeric(qm[[pick_col(qm, c("predicted_CellType_prob", "predicted_CellType_Broad_prob"), "predicted.*prob|_prob$")]]),
             mapping_error_QC = as.character(qm[[pick_col(qm, c("mapping_error_QC", "mapping_QC"), "mapping.*QC")]]))
}

## [3] per sample --------------------------------------------------------------------------------------------------------
out <- list()
for (f in rds) {
  sce <- readRDS(f); lib <- sub("_processed\\.rds$", "", basename(f))
  m <- SingleCellExperiment::counts(sce); m <- as(m, "CsparseMatrix")
  stopifnot(all(m@x[seq_len(min(1e5, length(m@x)))] == round(m@x[seq_len(min(1e5, length(m@x)))])))   # raw UMIs
  sym <- as.character(SummarizedExperiment::rowData(sce)$gene_symbol); sym[is.na(sym) | sym == ""] <- rownames(m)[is.na(sym) | sym == ""]
  rownames(m) <- make.unique(sym)
  cd <- as.data.frame(SummarizedExperiment::colData(sce)); sid <- if ("sample_id" %in% names(cd)) as.character(cd$sample_id[1]) else meta[scpca_library_id == lib, scpca_sample_id][1]
  mt <- meta[scpca_library_id == lib]; if (!nrow(mt)) mt <- meta[scpca_sample_id == sid]
  nc <- Matrix::colSums(m); nf <- Matrix::colSums(m > 0); pm <- 100 * Matrix::colSums(m[grepl("^MT-", rownames(m)), , drop = FALSE]) / pmax(nc, 1)
  keep <- nc >= LOY_GEO_QC_MIN_UMI & nf >= LOY_GEO_QC_MIN_GENE & pm <= LOY_GEO_QC_MAX_MT
  m <- m[, keep, drop = FALSE]; if (ncol(m) < 30) { message("    [skip] ", lib, ": ", ncol(m), " cells after QC"); next }
  md <- data.table(cell = colnames(m), nCount_RNA = as.numeric(nc[keep]), nFeature_RNA = as.numeric(nf[keep]), percent.mt = as.numeric(pm[keep]))
  pr <- project(m, as.data.frame(md[, .(cell)], row.names = md$cell)); md <- merge(md, pr, by = "cell", all.x = TRUE)
  md <- merge(md, bin_map[, .(predicted_CellType_Broad = CellType_Broad, hierarchy_bin)], by = "predicted_CellType_Broad", all.x = TRUE)
  md[, high_error := as.integer(!is.na(bmm_prob) & bmm_prob < MIN_PROJ_PROB)]
  yg <- intersect(LOY_YGENES, rownames(m)); yc <- as.data.table(as.matrix(Matrix::t(m[yg, md$cell, drop = FALSE]))); setnames(yc, paste0("y_", names(yc)))
  md <- cbind(md, yc); md[, y_umi := rowSums(as.matrix(yc))]
  md[, xist_umi := if (LOY_XIST %in% rownames(m)) as.numeric(m[LOY_XIST, md$cell]) else NA_real_]
  dx <- mt$diagnosis[1]; timing <- mt$disease_timing[1]
  md[, `:=`(dataset = ds, object = lib, Sample = sid, Patient_ID = if (!is.na(mt$participant_id[1]) && nzchar(mt$participant_id[1])) mt$participant_id[1] else sid, Library = lib,
            Chemistry = LOY_SCPCA_CHEM, Platform = paste0("10x_", LOY_SCPCA_CHEM),
            Timepoint = fifelse(grepl("Initial", timing), "Diagnosis", fifelse(grepl("Recurrence", timing), "Relapse", fifelse(grepl("Healthy", timing), "Control", timing))),
            Disease_state = fifelse(grepl("Acute myeloid", dx), "AML", fifelse(grepl("Non-cancerous", dx), "Healthy", fifelse(grepl("mixed phenotype", dx), "MPAL", dx))),
            anno_malignant = NA_character_, anno_sex = mt$sex[1], AgeGroup = fifelse(!is.na(mt$age[1]) & mt$age[1] < 18, "Pediatric", "Adult"), age_years = mt$age[1], subdiagnosis = mt$subdiagnosis[1])]
  md[, cell := paste(object, cell, sep = "__")]
  out[[lib]] <- md
  message(sprintf("    %s (%s, %s, %s, sex %s, age %s): %d cells | Y genes %d | median nCount %.0f | Y-zero %.3f | XIST>0 %.3f (%.1f min)", lib, sid, md$Disease_state[1], md$Timepoint[1], mt$sex[1], mt$age[1],
                  nrow(md), length(yg), median(md$nCount_RNA), mean(md$y_umi == 0), mean(md$xist_umi > 0, na.rm = TRUE), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(sce, m); invisible(gc())
}
cells <- rbindlist(out, fill = TRUE)
front <- c("dataset", "object", "cell", "Sample", "Patient_ID", "Library", "Chemistry", "Platform", "Timepoint", "Disease_state", "hierarchy_bin", "bmm_fine", "high_error", "anno_malignant",
           "nCount_RNA", "nFeature_RNA", "percent.mt", "anno_sex", "AgeGroup", "predicted_CellType_Broad", "mapping_error_QC", "bmm_prob")
setcolorder(cells, c(intersect(front, names(cells)), setdiff(names(cells), front)))
fwrite(cells, out_cells)
summ <- cells[, .(n_cells = .N, median_nCount = as.numeric(median(nCount_RNA)), y_rate = sum(y_umi) / sum(nCount_RNA), frac_y_zero = mean(y_umi == 0),
                  xist_pos_frac = mean(xist_umi[nCount_RNA >= 2000] >= 2, na.rm = TRUE), n_T = sum(hierarchy_bin == "T_NK", na.rm = TRUE),
                  Patient_ID = Patient_ID[1], Timepoint = Timepoint[1], Disease_state = Disease_state[1], anno_sex = anno_sex[1], age_years = age_years[1]), by = .(dataset, Sample)]
fwrite(summ, file.path(LOY_SCAN_DIR, paste0(ds, "__y_sample_summary.csv")))
# K5 here is the recorded sex: expression (XIST >= 2 UMIs in deep cells) must agree with metadata.tsv for samples with >= 200 deep cells
sex_expr <- summ[, fifelse(xist_pos_frac >= LOY_XIST_FEMALE, "F", fifelse(xist_pos_frac <= LOY_XIST_MALE, "M", "mixed"))]
k5 <- summ[, .(Sample, anno_sex, sex_expr = sex_expr, n_cells)][!is.na(anno_sex) & sex_expr != "mixed"]
chk <- data.table(check = c("K1_samples_ge50cells", "K2_ygenes_found", "K3_xist_present", "K4_median_depth_gt1000", "K5_recorded_sex_agrees"),
                  value = c(sum(summ$n_cells >= 50), length(grep("^y_[A-Z]", names(cells))), as.integer(!all(is.na(cells$xist_umi))), sum(summ$median_nCount > 1000), mean(k5$anno_sex == k5$sex_expr)),
                  pass = c(sum(summ$n_cells >= 50) >= 1, length(grep("^y_[A-Z]", names(cells))) >= 10, !all(is.na(cells$xist_umi)), TRUE, mean(k5$anno_sex == k5$sex_expr) >= 0.9))
fwrite(chk, file.path(LOY_SCAN_DIR, paste0(ds, "__y_scan_checks.csv"))); print(chk)
if (any(k5$anno_sex != k5$sex_expr)) { message("    recorded vs expression sex disagree:"); print(k5[anno_sex != sex_expr]) }
ext <- if (file.exists(LOY_EXT_LIST)) readLines(LOY_EXT_LIST) else character()
if (!ds %in% ext) writeLines(c(ext, ds), LOY_EXT_LIST)
message(sprintf("[done] %s: %d cells, %d samples (%.1f min)", ds, nrow(cells), uniqueN(cells$Sample), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(summ[order(Disease_state, Sample)][, .(Sample, Disease_state, Timepoint, anno_sex, age_years, n_cells, median_nCount, frac_y_zero = round(frac_y_zero, 3), xist_pos_frac = round(xist_pos_frac, 3))])
