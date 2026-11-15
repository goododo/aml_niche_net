#!/usr/bin/env Rscript
# 03b_scan_geo_mtx.R ----
# input : one GEO series of config section 7e (LOY_GEO_MTX): its RAW tar under RAW_PUBLIC_DIR/GEO/<GSE>_RAW/, the series matrix,
#         for GSE198681 the depositors' CMO sheet; the BoneMarrowMap scaffold (BMM_SCAFFOLD_RDS); LOY_YGENES / LOY_XIST
# output: LOY_PERCELL_DIR/<GSE>__y_percell.csv.gz (same schema as 01 / 03), LOY_SCAN_DIR/<GSE>__y_sample_summary.csv and
#         <GSE>__y_scan_checks.csv, the series id appended to LOY_EXT_LIST (so 02 picks it up); LOY_GEO_EXTRACT_DIR/<GSE>/ (tar extracted once)
# what it does: step 2b of the LOY line -- the per-cell chrY / XIST count table for series deposited as raw cellranger matrices
#   without labels. Per library: read the 10x triplet, keep Gene Expression features, light QC (LOY_GEO_QC_*), this project's
#   BoneMarrowMap projection (fine label, broad label -> hierarchy bin, mapping QC), chrY-gene and XIST UMIs per cell, patient /
#   timepoint / disease from the series matrix (one_donor layout) or from the CMO sheet after CMO demultiplexing (cmo layout).
#   Nothing here calls -Y: that is 02's frozen rule.
# Usage : sbatch scripts/loy/submit_loy_scan_geo.sbatch            # array over the LOY_GEO_MTX rows
#         Rscript scripts/loy/03b_scan_geo_mtx.R --dataset=GSE315856 [--force]
# WHY: the 2026-11-14 batch holds adults with two or more timepoints (GSE315856, GSE223844) and 13 untreated AML (GSE198681):
#   more chances of a subclonal -Y case with a second timepoint for the H2 design.
suppressPackageStartupMessages({ library(data.table); library(here); library(Matrix); library(Seurat); library(SeuratObject); library(BoneMarrowMap) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "loy", "config_loy.R"))
args  <- commandArgs(trailingOnly = TRUE)
ds    <- sub("^--dataset=", "", grep("^--dataset=", args, value = TRUE)); force <- any(args == "--force")
stopifnot(length(ds) == 1, ds %in% LOY_GEO_MTX$dataset)
row   <- LOY_GEO_MTX[dataset == ds]
raw_dir <- file.path(RAW_PUBLIC_DIR, "GEO", paste0(ds, "_RAW")); ext_dir <- file.path(LOY_GEO_EXTRACT_DIR, ds)
out_cells <- file.path(LOY_PERCELL_DIR, paste0(ds, "__y_percell.csv.gz"))
if (file.exists(out_cells) && !force) { message("[skip] exists: ", out_cells); quit(save = "no") }
dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE); dir.create(LOY_PERCELL_DIR, showWarnings = FALSE, recursive = TRUE); dir.create(LOY_SCAN_DIR, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); message(sprintf("[0] %s | layout %s | chem %s | %s", ds, row$layout, row$chem_lit, format(Sys.time())))

## [1] extract the tar once; list the matrices ---------------------------------------------------------------
tar <- file.path(raw_dir, paste0(ds, "_RAW.tar")); stopifnot(file.exists(tar))
if (!length(list.files(ext_dir, pattern = "matrix\\.mtx\\.gz$"))) { message("[1] extracting ", basename(tar)); system2("tar", c("xf", shQuote(tar), "-C", shQuote(ext_dir))) }
mtx <- list.files(ext_dir, pattern = "matrix\\.mtx\\.gz$", full.names = TRUE)
stem_of <- function(f) sub("_?matrix\\.mtx\\.gz$", "", basename(f))
libs <- data.table(mtx = mtx, stem = stem_of(mtx)); libs[, gsm := sub("_.*$", "", stem)]
message(sprintf("[1] %d matrices", nrow(libs)))

## [2] series matrix: GSM -> title / characteristics -------------------------------------------------------------
sm <- readLines(file.path(raw_dir, "series_matrix", paste0(ds, "_series_matrix.txt.gz")))
field <- function(key) { l <- grep(paste0("^!Sample_", key, "\\b"), sm, value = TRUE); lapply(l, function(x) gsub("\"", "", strsplit(x, "\t")[[1]][-1])) }
gsm_ids <- field("geo_accession")[[1]]; titles <- field("title")[[1]]
chars <- do.call(cbind, field("characteristics_ch1")); rownames(chars) <- gsm_ids
meta <- data.table(gsm = gsm_ids, title = titles)
char_val <- function(key) apply(chars, 1, function(v) { h <- grep(paste0("^", key, ": "), v, value = TRUE); if (length(h)) sub(paste0("^", key, ": "), "", h[1]) else NA_character_ })
libs <- merge(libs, meta, by = "gsm", all.x = TRUE)

## [3] per-dataset sample metadata (patient, timepoint, disease) ------------------------------------------------
if (ds == "GSE315856") {                       # title "..., AML1, TP1" ; library AML_Vtx<k>
  m <- regmatches(libs$title, regexec(", (AML[0-9]+), (TP[0-9])$", libs$title))
  libs[, `:=`(Patient_ID = vapply(m, `[`, "", 2), Timepoint = vapply(m, `[`, "", 3), Disease_state = "AML", Sample = sub("^GSM[0-9]+_", "", sub("_S[0-9]+$", "", stem)))]
  libs[, Timepoint := fifelse(Timepoint == "TP1", "Diagnosis", fifelse(Timepoint == "TP2", "On_treatment", Timepoint))]
} else if (ds == "GSE223844") {                # title "Screening [AML1007.1.GEX]"
  m <- regmatches(libs$title, regexec("^([A-Za-z0-9 ]+) \\[(AML[0-9]+)\\.([0-9]+)\\.GEX\\]$", libs$title))
  libs[, `:=`(Timepoint = vapply(m, `[`, "", 2), Patient_ID = vapply(m, `[`, "", 3), Disease_state = "AML", Sample = sub("^GSM[0-9]+_", "", sub("\\.GEX$", "", stem)))]
  libs[, treatment := char_val("treatment")[gsm]]
} else if (ds == "GSE343017") {                # per-donor RNA matrices: <GSM>_<donor>_<WBM|NK>_RNA ; HD = healthy donor
  m <- regmatches(libs$stem, regexec("^GSM[0-9]+_([A-Za-z0-9]+)_(WBM|NK)_(RNA|ADT_HTO)$", libs$stem))
  libs[, `:=`(donor = vapply(m, `[`, "", 2), fraction = vapply(m, `[`, "", 3), assay = vapply(m, `[`, "", 4))]
  libs <- libs[assay == "RNA" & fraction == "WBM"]                                   # NK-sorted lanes and antibody matrices are not scanned
  libs[, `:=`(Patient_ID = donor, Sample = paste0(donor, "_WBM"), Timepoint = "Diagnosis", Disease_state = fifelse(grepl("^HD", donor), "Healthy", "AML"))]
} else if (ds == "GSE198681") {                # one matrix per lane; donors by CMO (sheet)
  libs[, lane := sub("^GSM[0-9]+_", "", stem)]
  libs[, lib_name := fifelse(grepl("^scAML", lane), sub("^sc", "", lane), lane)]         # scAML4 -> AML4 ; scRNAseq7 stays
  cmo_sheet <- data.table::rbindlist(list(
    data.table(patient = c("HEMA-090", "HEMA-096", "HEMA-100", "HEMA-102", "HEMA-103", "HEMA-105", "HEMA-121", "HEMA-131", "HEMA-085", "HEMA-109", "HEMA-117", "HEMA-129", "HEMA-136"),
               group = c(rep("CR", 8), rep("non-CR", 5)),
               lib_name = c("AML4", "AML4", "AML4", "AML1", "AML1", "AML1", "AML2", "AML4", "AML1", "AML1", "AML0", "AML4", "scRNAseq7"),
               cmo = c("CMO301", "CMO302", "CMO309", "CMO303,CMO304", "CMO305,CMO306", "CMO307,CMO308", "CMO309,CMO310", "CMO312", "CMO301,CMO302", "CMO311,CMO312", "CMO305", "CMO310", "CMO309"))))
  # transcribed from GSE198681_CMO_conditions.xlsx (sheet 1, rows 3-15) on 2026-11-15; the xlsx stays the source of truth
  cmo_map <- cmo_sheet[, .(cmo = strsplit(cmo, ",")[[1]]), by = .(patient, group, lib_name)]
}
message(sprintf("[3] libraries to scan: %d", nrow(libs))); print(libs[, intersect(c("gsm", "Sample", "Patient_ID", "Timepoint", "Disease_state", "lane"), names(libs)), with = FALSE][1:min(.N, 8)])

## [4] scaffold + bin map ------------------------------------------------------------------------------------------
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
read_triplet <- function(f) {               # cellranger triplet -> list(gex = dgCMatrix genes x cells, cmo = CMO counts or NULL)
  # sibling files are found by glob on the GSM stem: one GSE198681 lane is named "<GSM>_scAML2matrix.mtx.gz" (no underscore)
  # while its features/barcodes files keep the underscore, so paste0(base, "features.tsv.gz") does not exist for it
  sibling <- function(what) { hit <- list.files(dirname(f), paste0("^", stem_of(f), "_?", what, "\\.tsv\\.gz$"), full.names = TRUE)
    if (length(hit) != 1L) stop("expected one ", what, " file for ", basename(f), ", found ", length(hit)); hit }
  feats <- fread(sibling("features"), header = FALSE); cells <- fread(sibling("barcodes"), header = FALSE)[[1]]
  m <- Matrix::readMM(f); m <- as(m, "CsparseMatrix"); colnames(m) <- cells
  sym <- if (ncol(feats) >= 2) feats[[2]] else feats[[1]]; type <- if (ncol(feats) >= 3) feats[[3]] else rep("Gene Expression", nrow(feats))
  rownames(m) <- make.unique(sym)
  list(gex = m[type == "Gene Expression", , drop = FALSE], cmo = if (any(type == "Multiplexing Capture")) m[type == "Multiplexing Capture", , drop = FALSE] else NULL)
}

## [5] per library ---------------------------------------------------------------------------------------------------
out <- list(); summ <- list()
for (i in seq_len(nrow(libs))) {
  L <- libs[i]; tr <- read_triplet(L$mtx); m <- tr$gex
  stopifnot(all(m@x[seq_len(min(1e5, length(m@x)))] == round(m@x[seq_len(min(1e5, length(m@x)))])))
  nc <- Matrix::colSums(m); nf <- Matrix::colSums(m > 0); mt <- 100 * Matrix::colSums(m[grepl("^MT-", rownames(m)), , drop = FALSE]) / pmax(nc, 1)
  keep <- nc >= LOY_GEO_QC_MIN_UMI & nf >= LOY_GEO_QC_MIN_GENE & mt <= LOY_GEO_QC_MAX_MT
  m <- m[, keep, drop = FALSE]; nc <- nc[keep]; nf <- nf[keep]; mt <- mt[keep]
  if (ncol(m) < 30) { message(sprintf("    [skip] %s: %d cells after QC", L$stem, ncol(m))); next }
  md <- data.table(cell = colnames(m), nCount_RNA = as.numeric(nc), nFeature_RNA = as.numeric(nf), percent.mt = as.numeric(mt))
  if (row$layout == "cmo") {                # demultiplex: dominant CMO with >= LOY_GEO_CMO_MIN_UMI UMIs and >= LOY_GEO_CMO_RATIO x the runner-up
    cm <- as.matrix(tr$cmo[, colnames(m), drop = FALSE]); o <- apply(cm, 2, function(v) { s <- sort(v, decreasing = TRUE); c(names(s)[1], s[1], s[2]) })
    md[, `:=`(cmo = o[1, ], cmo_top = as.numeric(o[2, ]), cmo_second = as.numeric(o[3, ]))]
    md[, cmo_ok := cmo_top >= LOY_GEO_CMO_MIN_UMI & cmo_top >= LOY_GEO_CMO_RATIO * pmax(cmo_second, 1)]
    md <- merge(md, cmo_map[lib_name == L$lib_name, .(cmo, Patient_ID = patient, response = group)], by = "cmo", all.x = TRUE)
    message(sprintf("    %s: %d cells after QC; CMO assigned %d (%.0f%%), to a listed patient %d", L$stem, nrow(md), sum(md$cmo_ok), 100 * mean(md$cmo_ok), sum(md$cmo_ok & !is.na(md$Patient_ID))))
    md <- md[cmo_ok == TRUE & !is.na(Patient_ID)]; m <- m[, md$cell, drop = FALSE]
    md[, `:=`(Sample = Patient_ID, Timepoint = "Diagnosis", Disease_state = "AML", Library = L$lib_name)]
  } else {
    md[, `:=`(Sample = L$Sample, Patient_ID = L$Patient_ID, Timepoint = L$Timepoint, Disease_state = L$Disease_state, Library = L$Sample)]
    if ("treatment" %in% names(L)) md[, treatment := L$treatment]
  }
  pr <- project(m, as.data.frame(md[, .(cell)], row.names = md$cell))
  md <- merge(md, pr, by = "cell", all.x = TRUE)
  md <- merge(md, bin_map[, .(predicted_CellType_Broad = CellType_Broad, hierarchy_bin)], by = "predicted_CellType_Broad", all.x = TRUE)
  md[, high_error := as.integer(!is.na(bmm_prob) & bmm_prob < MIN_PROJ_PROB)]
  yg <- intersect(LOY_YGENES, rownames(m)); yc <- as.data.table(as.matrix(Matrix::t(m[yg, md$cell, drop = FALSE]))); setnames(yc, paste0("y_", names(yc)))
  md <- cbind(md, yc); md[, y_umi := rowSums(as.matrix(yc))]
  md[, xist_umi := if (LOY_XIST %in% rownames(m)) as.numeric(m[LOY_XIST, md$cell]) else NA_real_]
  md[, `:=`(dataset = ds, object = L$stem, Chemistry = row$chem_lit, Platform = paste0("10x_", row$chem_lit), anno_malignant = NA_character_, anno_sex = NA_character_, AgeGroup = "Adult")]
  md[, cell := paste(object, cell, sep = "__")]
  out[[L$stem]] <- md
  message(sprintf("    %s: %d cells | samples %s | Y genes %d | median nCount %.0f | Y-zero %.3f | XIST>0 %.3f (%.1f min)", L$stem, nrow(md), paste(unique(md$Sample), collapse = ","),
                  length(yg), median(md$nCount_RNA), mean(md$y_umi == 0), mean(md$xist_umi > 0, na.rm = TRUE), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cells <- rbindlist(out, fill = TRUE)
front <- c("dataset", "object", "cell", "Sample", "Patient_ID", "Library", "Chemistry", "Platform", "Timepoint", "Disease_state", "hierarchy_bin", "bmm_fine", "high_error", "anno_malignant",
           "nCount_RNA", "nFeature_RNA", "percent.mt", "anno_sex", "AgeGroup", "predicted_CellType_Broad", "mapping_error_QC", "bmm_prob")
setcolorder(cells, c(intersect(front, names(cells)), setdiff(names(cells), front)))
fwrite(cells, out_cells)
summ <- cells[, .(n_cells = .N, median_nCount = as.numeric(median(nCount_RNA)), y_rate = sum(y_umi) / sum(nCount_RNA), frac_y_zero = mean(y_umi == 0),
                  xist_pos_frac = mean(xist_umi[nCount_RNA >= 2000] >= 2, na.rm = TRUE), n_T = sum(hierarchy_bin == "T_NK", na.rm = TRUE),
                  Patient_ID = Patient_ID[1], Timepoint = Timepoint[1], Disease_state = Disease_state[1]), by = .(dataset, Sample)]
fwrite(summ, file.path(LOY_SCAN_DIR, paste0(ds, "__y_sample_summary.csv")))
chk <- data.table(check = c("K1_samples_ge50cells", "K2_ygenes_found", "K3_xist_present", "K4_median_depth_gt1000", "K5_bins_assigned"),
                  value = c(sum(summ$n_cells >= 50), length(intersect(LOY_YGENES, sub("^y_", "", grep("^y_[A-Z]", names(cells), value = TRUE)))), as.integer(!all(is.na(cells$xist_umi))), sum(summ$median_nCount > 1000), mean(!is.na(cells$hierarchy_bin))),
                  pass = c(sum(summ$n_cells >= 50) >= 1, length(intersect(LOY_YGENES, sub("^y_", "", grep("^y_[A-Z]", names(cells), value = TRUE)))) >= 10, !all(is.na(cells$xist_umi)), TRUE, mean(!is.na(cells$hierarchy_bin)) >= 0.9))
fwrite(chk, file.path(LOY_SCAN_DIR, paste0(ds, "__y_scan_checks.csv"))); print(chk)
ext <- if (file.exists(LOY_EXT_LIST)) readLines(LOY_EXT_LIST) else character()
if (!ds %in% ext) writeLines(c(ext, ds), LOY_EXT_LIST)
message(sprintf("[done] %s: %d cells, %d samples, %d patients (%.1f min)", ds, nrow(cells), uniqueN(cells$Sample), uniqueN(cells$Patient_ID), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(summ[order(Patient_ID, Timepoint)][, .(Sample, Patient_ID, Timepoint, Disease_state, n_cells, median_nCount, frac_y_zero = round(frac_y_zero, 3), xist_pos_frac = round(xist_pos_frac, 3))])
