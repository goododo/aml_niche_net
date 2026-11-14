#!/usr/bin/env Rscript
# 03_scan_external.R ----
# INPUT  : one Zeng-pack object  00_raw/public/Zenodo/zenodo.18528442/<file_stem>_projected_cleaned.rds, or the two
#          GSE269669 author objects (config section 7c)
#          (LOY_EXT_STUDIES in config_loy.R section 7; --study=<study> or --task=<n> over external_datasets.txt)
# OUTPUT : LARGE1/loy/percell/ZENG_<study>__y_percell.csv.gz   same columns as 01_scan_y_counts.R output, plus
#            AgeGroup, Differentiation_State, mapping_error_QC, predicted_CellType_Broad
#          results/tables/loy/01_scan/ZENG_<study>__y_sample_summary.csv, ZENG_<study>__y_scan_checks.csv
# what it does: the step-0 count scan for the external pack. Maps the pack's metadata onto the project's column
#   names so that 02_call_loy_status.R reads it unchanged: Sample <- SampleID_Unique (the pack's batch variable),
#   hierarchy_bin <- bmm_bin_map.tsv[predicted_CellType_Broad], bmm_fine <- predicted_CellType, high_error <-
#   mapping_error_QC == "Fail", anno_malignant <- NormalvsMalignant, Disease_state <- Diagnosis. Chemistry is
#   left "unknown" on purpose (the pack has none); 02 takes chem_lit from config section 7.
# Usage  : Rscript scripts/loy/03_scan_external.R --study=vanGalen_2019     |  --dataset=GSE269669 (section 7c branch)
#          Rscript scripts/loy/03_scan_external.R --task=3      (3rd line of results/tables/loy/01_scan/external_datasets.txt)
# Checks : K1 every SampleID_Unique produced rows; K2 >= 8 Y genes on the gene axis; K3 rows == ncol(object);
#          K4 every predicted_CellType_Broad value (except NA) is in bmm_bin_map.tsv (unmapped labels are listed).
# WHY a separate script and not a flag in 01: 01 walks a directory of per-sample objects with the project's
#   metadata; the pack is one object per study with different column names. Keeping the mapping in one place
#   makes the vocabulary translation reviewable (K4) instead of buried in a branch.

suppressPackageStartupMessages({ library(here); library(Seurat); library(Matrix); library(data.table) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "loy", "config_loy.R"))
set.seed(SEED)

args <- commandArgs(TRUE)
getarg <- function(key) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), args, value = TRUE)); if (length(v)) v[1] else NA }
study_id <- getarg("study"); task <- getarg("task"); ds_arg <- getarg("dataset")
if (!is.na(task)) ds_arg <- readLines(LOY_EXT_LIST)[as.integer(task)]
if (!is.na(ds_arg) && grepl("^ZENG_", ds_arg)) study_id <- sub("^ZENG_", "", ds_arg)
is_geo <- !is.na(ds_arg) && ds_arg %in% LOY_EXT_GEO$dataset          # section 7c objects (GSE269669)
if (!is_geo) {
  stopifnot(!is.na(study_id), study_id %in% LOY_EXT_STUDIES$study)
  row <- LOY_EXT_STUDIES[study == study_id] ; ds <- row$dataset   # data.table i-expression: column `study` vs. the R variable
} else { ds <- ds_arg; row <- LOY_EXT_GEO[dataset == ds] }
dir.create(LOY_PERCELL_DIR, showWarnings = FALSE, recursive = TRUE); dir.create(LOY_SCAN_DIR, showWarnings = FALSE, recursive = TRUE)

out_cells <- file.path(LOY_PERCELL_DIR, paste0(ds, "__y_percell.csv.gz"))
if (file.exists(out_cells) && is.na(getarg("force"))) { message("[skip] exists: ", out_cells); quit(save = "no") }

# y_block(): the counting step shared by both branches -- Y-gene and XIST UMIs per cell, appended to a metadata table
y_block <- function(m, dt) {
  yg <- intersect(LOY_YGENES, rownames(m))
  ycounts <- as.data.table(as.matrix(Matrix::t(m[yg, , drop = FALSE]))); setnames(ycounts, paste0("y_", colnames(ycounts)))
  dt <- cbind(dt, ycounts)
  dt[, y_umi := rowSums(as.matrix(ycounts))]
  dt[, xist_umi := if (LOY_XIST %in% rownames(m)) as.numeric(m[LOY_XIST, ]) else NA_real_]
  attr(dt, "ygenes") <- yg
  dt
}
counts_of <- function(s) {
  m <- tryCatch(LayerData(s, assay = "RNA", layer = "counts"), error = function(e) GetAssayData(s, assay = "RNA", slot = "counts"))
  stopifnot(all(m@x[seq_len(min(1e5, length(m@x)))] == round(m@x[seq_len(min(1e5, length(m@x)))])))   # raw UMI counts
  m
}

if (!is_geo) {
## ================= Zeng pack branch =================================================================
f <- file.path(LOY_ZENG_RDS_DIR, paste0(row$file_stem, "_projected_cleaned.rds"))
stopifnot(file.exists(f))
message(sprintf("[1] %s: reading %s (%.2f GB)", ds, basename(f), file.size(f) / 1e9))
s <- readRDS(f)
m <- counts_of(s)
md <- as.data.table(s@meta.data, keep.rownames = "cell")
stopifnot(identical(md$cell, colnames(m)))
rm(s); invisible(gc())
message(sprintf("    cells %d  genes %d  samples %d", ncol(m), nrow(m), uniqueN(md$SampleID_Unique)))

## -- vocabulary translation ----------------------------------------------------------------------
binmap <- fread(LOY_BMM_BIN_MAP)
md[, hierarchy_bin := binmap$hierarchy_bin[match(predicted_CellType_Broad, binmap$CellType_Broad)]]
unmapped <- setdiff(unique(na.omit(md$predicted_CellType_Broad)), binmap$CellType_Broad)
dt <- data.table(
  dataset = ds, object = md$SampleID_Unique, cell = md$cell,
  Sample = md$SampleID_Unique, Patient_ID = md$PatientID, Library = md$SampleID_Unique,
  Chemistry = "unknown", Platform = row$chem_lit, Timepoint = md$Timepoint,
  Disease_state = md$Diagnosis, hierarchy_bin = md$hierarchy_bin, bmm_fine = md$predicted_CellType,
  high_error = as.integer(md$mapping_error_QC == "Fail"), anno_malignant = md$NormalvsMalignant,
  nCount_RNA = as.numeric(md$nCount_RNA), nFeature_RNA = md$nFeature_RNA, percent.mt = md$pct_Mito,
  anno_sex = NA_character_, AgeGroup = md$AgeGroup, Differentiation_State = md$Differentiation_State,
  mapping_error_QC = md$mapping_error_QC, predicted_CellType_Broad = md$predicted_CellType_Broad,
  sample_label = md$Sample)
if (all(is.na(dt$nCount_RNA))) dt[, nCount_RNA := Matrix::colSums(m)]
dt <- y_block(m, dt); yg <- attr(dt, "ygenes"); n_expected <- nrow(md); n_samples_expected <- uniqueN(md$SampleID_Unique)
rm(m); invisible(gc())
} else {
## ================= GSE269669 branch (section 7c) =====================================================
# two author objects; the tme object lacks Patient/time, so both come from the library id (PTnnX) and the aml object
parts <- list(); yg <- NULL; n_expected <- 0L; unmapped <- character()
tp_map <- NULL
for (obj in c("avm_aml", "avm_tme")) {
  f <- file.path(LOY_GSE269669_DIR, paste0("GSE269669_", obj, ".rds.gz")); stopifnot(file.exists(f))
  message(sprintf("[1] %s: reading %s (%.2f GB)", ds, basename(f), file.size(f) / 1e9))
  con <- gzcon(gzfile(f, "rb")); s <- readRDS(con); close(con)          # gzip of a gzip-compressed RDS
  m <- counts_of(s); md <- as.data.table(s@meta.data, keep.rownames = "cell"); rm(s); invisible(gc())
  stopifnot(identical(md$cell, colnames(m)))
  if (obj == "avm_aml") tp_map <- unique(md[, .(Sample = as.character(orig.ident), Timepoint = as.character(time), BestResponse = as.character(BestResponse2))])
  md[, hierarchy_bin := unname(LOY_AZIMUTH_L2_TO_BIN[as.character(predicted.celltype.l2)])]
  unmapped <- union(unmapped, setdiff(unique(as.character(md$predicted.celltype.l2)), names(LOY_AZIMUTH_L2_TO_BIN)))
  dt <- data.table(
    dataset = ds, object = obj, cell = md$cell,
    Sample = as.character(md$orig.ident), Patient_ID = substr(as.character(md$orig.ident), 1, 4), Library = as.character(md$orig.ident),
    Chemistry = "unknown", Platform = row$chem_lit, Timepoint = NA_character_,
    Disease_state = "AML", hierarchy_bin = md$hierarchy_bin, bmm_fine = as.character(md$predicted.celltype.l2),
    high_error = 0L, anno_malignant = if (obj == "avm_aml") "malignant" else "normal",
    nCount_RNA = as.numeric(md$nCount_RNA), nFeature_RNA = md$nFeature_RNA, percent.mt = 100 * md$mitoRatio,
    anno_sex = NA_character_, AgeGroup = "Adult", author_class = as.character(md$class1),
    symphony_granja = as.character(md$Symphony_granja), l2_score = md$predicted.celltype.l2.score)
  dt <- y_block(m, dt); yg <- union(yg, attr(dt, "ygenes")); n_expected <- n_expected + nrow(md)
  parts[[obj]] <- dt; rm(m, md, dt); invisible(gc())
}
dt <- rbindlist(parts, fill = TRUE)
dt[, Timepoint := tp_map$Timepoint[match(Sample, tp_map$Sample)]]          # Pre / Post from the aml object; NA for tme-only libraries
dt[, BestResponse := tp_map$BestResponse[match(Sample, tp_map$Sample)]]
n_samples_expected <- uniqueN(dt$Sample)
message(sprintf("    libraries %d (with Pre/Post: %d) | aml cells %d | tme cells %d", n_samples_expected, uniqueN(tp_map$Sample), nrow(parts$avm_aml), nrow(parts$avm_tme)))
}
fwrite(dt, out_cells)
message("[2] wrote ", out_cells, " rows ", nrow(dt), " | Y genes on axis: ", paste(yg, collapse = ","))
print(dt[, .(n = .N, med_ncount = median(nCount_RNA), frac_y0 = round(mean(y_umi == 0), 3), xist_pos = round(mean(xist_umi > 0), 3)),
         by = .(Sample, Disease_state, Timepoint)][order(Sample)])

## -- per sample x bin x depth summary (same shape as 01) ----------------------------------------------
dt[, depth_bin := cut(nCount_RNA, LOY_DEPTH_BREAKS, right = FALSE)]
dt[, bin := fifelse(is.na(hierarchy_bin), "NA", hierarchy_bin)]
summ <- dt[, .(n = .N, n_y0 = sum(y_umi == 0), n_y1 = sum(y_umi == 1), n_y2p = sum(y_umi >= 2),
               med_ncount = median(nCount_RNA), sum_y = sum(y_umi), sum_ncount = sum(nCount_RNA),
               xist_pos_frac = mean(xist_umi > 0, na.rm = TRUE), n_high_error = sum(high_error == 1L)),
           by = .(dataset, object, Sample, Library, Chemistry, Platform, Timepoint, Disease_state, bin, depth_bin)]
summ[, y_rate := sum_y / sum_ncount]
fwrite(summ, file.path(LOY_SCAN_DIR, paste0(ds, "__y_sample_summary.csv")))

## -- known answers --------------------------------------------------------------------------------
chk <- data.table(check = c("K1_samples_with_rows", "K2_ygenes_found", "K3_rows_equal_ncol", "K4_broad_labels_mapped"),
                  expected = c(n_samples_expected, ">=8", n_expected, "0 unmapped"),
                  observed = c(uniqueN(dt$Sample), length(yg), nrow(dt), paste0(length(unmapped), " unmapped", if (length(unmapped)) paste0(": ", paste(unmapped, collapse = ";")) else "")))
chk[, pass := c(observed[1] == expected[1], length(yg) >= 8, observed[3] == expected[3], length(unmapped) == 0)]
fwrite(chk, file.path(LOY_SCAN_DIR, paste0(ds, "__y_scan_checks.csv")))
print(chk)
message("[done] ", ds)
