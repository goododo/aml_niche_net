#!/usr/bin/env Rscript
# 03_prepare_integration_input.R ----
# INPUT  : /LARGE1/gr10634/gaozy/leukemia_ml/03.2_anno/inte_rough_anno.RDS   (42.5 GB, 161 samples)
#          /FAST/gr10634/gaozy/leukemia_ml/05_CCC_per_sample/sample_disease_corrected.csv
#          /LARGE1/gr10634/gaozy/aml_niche_net/02_seurat_objects/04_annotated/GSE185381/*.rds (52)
#          results/tables/01_preprocess/00_curated_manifest.csv
# OUTPUT : /LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/01_input/<Dataset_ID>.rds
#          /FAST/gr10634/gaozy/aml_niche_net/results/tables/ws/ws_reintegration_manifest.csv
# WHAT IT DOES : step 1 of the re-integration the user approved on 2026-10-05. The leukemia_ml
#          integration inputs were deleted from LARGE0, so the 42.5 GB annotated object is the only
#          surviving copy of the QC'd cells for 161 samples; it does still carry integer RNA counts.
#          This script takes counts out of it per dataset, applies two approved corrections, adds
#          GSE185381, and writes per-dataset objects so no later step has to touch 42.5 GB again.
# Usage  : Rscript scripts/ws/03_prepare_integration_input.R
#
# THE TWO APPROVED CORRECTIONS
#   1. DROP four datasets that contribute nothing to a within-sample diagnostic design:
#      GSE185991 (100% sorted, 0 diagnostic unsorted, worst QC retention 0.47 / MAD loss 0.51),
#      E-MTAB-11536 (most dissociation stress, worst doublet over-removal), GSE207356 (1 patient,
#      0 diagnostic), GSE147989 (4 sorted PB samples, 0 diagnostic). GSE116256 is not in the atlas
#      and is not added (Seq-Well: 505 cells/sample, 1,917 UMI, 819 genes/cell).
#   2. FIX the Chen2023 (PMID37470778) label bug: 8 normal-bone-marrow fraction samples
#      (NBM1-4 x CD34/Niche_Immune, 51,613 cells) are recorded as SampleType "AML" in the atlas.
#      sample_disease_corrected.csv is the authority; the script asserts that exactly 8 rows change.
# BoneMarrowMap is KEPT (user reversed the earlier decision on 2026-10-05). Note it is kept for
# integration and annotation only -- the standing decision to exclude it from ANALYSIS is unchanged.

suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(data.table) })
cfg <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "..", "config", "config_paths.R")
source(cfg); set.seed(SEED)

ATLAS   <- "/LARGE1/gr10634/gaozy/leukemia_ml/03.2_anno/inte_rough_anno.RDS"
CORR    <- "/FAST/gr10634/gaozy/leukemia_ml/05_CCC_per_sample/sample_disease_corrected.csv"
OUT_RDS <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/01_input"
OUT_TAB <- file.path(TAB_DIR, "ws")
AN_185  <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated/GSE185381")
DROP    <- c("GSE185991", "E-MTAB-11536", "GSE207356", "GSE147989")
KEEP_MD <- c("Sample_Name", "Dataset_ID", "SampleType", "TimePoint",
             "Manual_anno_2", "Sub_anno", "DataType", "Original_CellType", "percent.mt")
dir.create(OUT_RDS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
CHK <- list(); chk <- function(name, expected, observed) {
  pass <- identical(as.character(expected), as.character(observed))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected),
                                        observed = as.character(observed), pass = pass)
  message(sprintf("  %-46s expect %-10s got %-10s %s", name, expected, observed,
                  if (pass) "PASS" else "** FAIL **"))
}

## -- 1. the atlas ---------------------------------------------------------------------------------
message("[1] reading the atlas object (42.5 GB, several minutes) ...")
A <- readRDS(ATLAS)
md <- as.data.table(A@meta.data, keep.rownames = "cell")
chk("atlas cells", 808539, ncol(A))
chk("atlas samples", 161, uniqueN(md$Sample_Name))
chk("atlas datasets", 12, uniqueN(md$Dataset_ID))
if (!"percent.mt" %in% colnames(A@meta.data)) {
  message("    percent.mt absent -- computing from MT- genes (needed by ScaleData vars.to.regress)")
  A[["percent.mt"]] <- PercentageFeatureSet(A, pattern = "^MT-", assay = "RNA")
  md <- as.data.table(A@meta.data, keep.rownames = "cell")
}

## -- 2. correction 2: the Chen2023 label bug ------------------------------------------------------
message("\n[2] applying the Chen2023 SampleType correction")
co <- fread(CORR)[, .(Sample_Name, Dataset_ID, disease, disease_class)]
md <- merge(md, co, by = c("Sample_Name", "Dataset_ID"), all.x = TRUE, sort = FALSE)
md[, SampleType_new := SampleType]
md[!is.na(disease_class) & disease_class == "healthy", SampleType_new := "Control"]
md[!is.na(disease_class) & disease_class == "malignant", SampleType_new := "AML"]
changed <- md[SampleType_new != SampleType, .(cells = .N), by = .(Dataset_ID, Sample_Name, SampleType, SampleType_new)]
print(changed)
chk("samples whose SampleType changes", 8, nrow(changed))
chk("cells relabelled", 51613, sum(changed$cells))
chk("all relabelled are PMID37470778", "TRUE", as.character(all(changed$Dataset_ID == "PMID37470778")))
chk("all relabelled go AML -> Control", "TRUE",
    as.character(all(changed$SampleType == "AML" & changed$SampleType_new == "Control")))
setkey(md, cell); A$SampleType <- md[J(colnames(A)), SampleType_new]

## -- 3. correction 1: drop four datasets ----------------------------------------------------------
message("\n[3] dropping ", paste(DROP, collapse = ", "))
keep_cells <- colnames(A)[!A$Dataset_ID %in% DROP]
chk("cells kept from the atlas", 728092, length(keep_cells))
A <- A[, keep_cells]
chk("datasets kept", 8, uniqueN(A$Dataset_ID))

## -- 4. write one object per dataset, counts only -------------------------------------------------
message("\n[4] writing per-dataset objects (RNA counts only; integrated/progeny/dorothea dropped)")
man <- list()
for (ds in sort(unique(A$Dataset_ID))) {
  j <- which(A$Dataset_ID == ds)
  cnt <- SeuratObject::LayerData(A, assay = "RNA", layer = "counts")[, j, drop = FALSE]
  m <- A@meta.data[j, intersect(KEEP_MD, colnames(A@meta.data)), drop = FALSE]
  o <- CreateSeuratObject(counts = cnt, meta.data = m, project = ds)
  saveRDS(o, file.path(OUT_RDS, paste0(ds, ".rds")))
  man[[length(man) + 1]] <- data.table(Dataset_ID = ds, source = "atlas", samples = uniqueN(m$Sample_Name),
                                       cells = length(j), genes = nrow(cnt),
                                       types = paste(sort(unique(m$SampleType)), collapse = "/"),
                                       annotated = sum(!is.na(m$Manual_anno_2)))
  message(sprintf("    %-16s %6d cells | %3d samples | %d genes",
                  ds, length(j), uniqueN(m$Sample_Name), nrow(cnt)))
  rm(cnt, m, o); gc(FALSE)
}
rm(A); gc(FALSE)

## -- 5. add GSE185381 -----------------------------------------------------------------------------
message("\n[5] ingesting GSE185381 (52 post-QC objects, zero QC exclusions)")
cm <- fread(file.path(TAB_DIR, "01_preprocess", "00_curated_manifest.csv"))[dataset == "GSE185381"]
fs <- list.files(AN_185, pattern = "\\.rds$", full.names = TRUE)
chk("GSE185381 objects on disk", 52, length(fs))
parts <- vector("list", length(fs))
for (i in seq_along(fs)) {
  sm <- sub("\\.rds$", "", basename(fs[i]))
  s <- readRDS(fs[i]); r <- cm[sample == sm]
  m <- data.frame(
    Sample_Name = sm, Dataset_ID = "GSE185381",
    SampleType = if (nrow(r) && r$disease[1] == "AML") "AML" else "Control",
    TimePoint = if (nrow(r)) r$timepoint[1] else NA_character_,
    Manual_anno_2 = NA_character_, Sub_anno = NA_character_,   # to be annotated in 03.2
    DataType = if (nrow(r)) r$platform_fine[1] else NA_character_,
    Original_CellType = if ("bmm_fine" %in% colnames(s@meta.data)) as.character(s$bmm_fine) else NA_character_,
    percent.mt = if ("percent.mt" %in% colnames(s@meta.data)) s$percent.mt else NA_real_,
    row.names = colnames(s))
  parts[[i]] <- CreateSeuratObject(counts = SeuratObject::LayerData(s, assay = "RNA", layer = "counts"),
                                   meta.data = m, project = sm)
  rm(s); gc(FALSE)
  if (i %% 10 == 0) message(sprintf("    ... %d/%d", i, length(fs)))
}
G <- if (length(parts) == 1L) parts[[1]] else merge(parts[[1]], parts[-1])
G <- suppressWarnings(JoinLayers(G))
rm(parts); gc(FALSE)
chk("GSE185381 cells", 167136, ncol(G))
chk("GSE185381 samples", 52, uniqueN(G$Sample_Name))
saveRDS(G, file.path(OUT_RDS, "GSE185381.rds"))
man[[length(man) + 1]] <- data.table(Dataset_ID = "GSE185381", source = "aml_niche_net",
  samples = uniqueN(G$Sample_Name), cells = ncol(G), genes = nrow(G),
  types = paste(sort(unique(G$SampleType)), collapse = "/"), annotated = 0L)
rm(G); gc(FALSE)

## -- 6. manifest and the gene space the integration will have -------------------------------------
M <- rbindlist(man)
message("\n[6] re-integration input")
print(M)
message(sprintf("    TOTAL %d datasets | %d samples | %d cells",
                nrow(M), sum(M$samples), sum(M$cells)))
gl <- lapply(file.path(OUT_RDS, paste0(M$Dataset_ID, ".rds")),
             function(f) { o <- readRDS(f); r <- rownames(o); rm(o); gc(FALSE); r })
names(gl) <- M$Dataset_ID
M[, genes_in_file := vapply(gl, length, 1L)]
inter <- Reduce(intersect, gl); uni <- Reduce(union, gl)
message(sprintf("    gene intersection across the 9 datasets: %d | union: %d", length(inter), length(uni)))
M[, genes_missing_vs_union := vapply(gl, function(g) length(setdiff(uni, g)), 1L)]
fwrite(M, file.path(OUT_TAB, "ws_reintegration_manifest.csv"))
fwrite(data.table(gene = inter), file.path(OUT_TAB, "ws_reintegration_common_genes.csv"))
chk("datasets prepared", 9, nrow(M))
chk("samples prepared", 182, sum(M$samples))
chk("cells prepared", 895228, sum(M$cells))

R <- rbindlist(CHK); fwrite(R, file.path(OUT_TAB, "ws_reintegration_input_checks.csv"))
message(sprintf("\n[7] %d checks, %d PASS, %d FAIL", nrow(R), sum(R$pass), sum(!R$pass)))
if (any(!R$pass)) print(R[pass == FALSE])
message("[done] ", OUT_RDS)
