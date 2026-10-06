#!/usr/bin/env Rscript
# INPUT  : /LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/01_input/<Dataset_ID>.rds  (9 datasets)
#          results/tables/ws/ws_reintegration_common_genes.csv                           (32,042 genes)
# OUTPUT : /LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated/ws_rpca_integrated.rds
#          results/tables/ws/ws_rpca_checks.csv          every assertion, expected vs got
#          results/tables/ws/ws_rpca_mixing.csv          neighbour dataset-mixing, before vs after
#          results/figures/ws/ws_rpca_umap.pdf           UMAP by Dataset_ID / SampleType / phase
# what it does: re-runs 02.1_Compare_inte_methods.r's method2_seurat_rpca on the new 9-dataset
#   manifest (BoneMarrowMap kept, GSE185991/E-MTAB-11536/GSE207356/GSE147989 dropped, GSE185381
#   added, Chen2023 SampleType corrected). Split variable is Dataset_ID, matching that script's
#   one-object-per-dataset list. Anchor parameters are left at the defaults it used.
#
# Three departures from 02.1, each deliberate:
#   1. the 13 MT- genes are dropped from the gene axis before NormalizeData. They carry UMIs in
#      only 3 of the 9 datasets -- the other 6 were downloaded already QC-filtered, with the MT
#      rows zeroed by the original authors -- so keeping them makes the library-size denominator
#      dataset-dependent and lets a perfectly dataset-discriminating feature into the anchor set.
#      (approved 2026-10-06)
#   2. RunTSNE is skipped: hours at 895k cells, and UMAP answers the same batch question.
#   3. the duplicate-embedding guard is rewritten. 02.1 does colnames(obj)[-duplicate_cells] with
#      duplicate_cells a CHARACTER vector, which is an R error ("invalid argument to unary
#      operator"); it never fired only because no duplicate was ever found. setdiff() here.
#
# This script deletes nothing. The step-3 checkpoint is kept on purpose so that an OOM in
# anchoring or integration does not cost the normalize/scale/PCA work again.

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(Seurat); library(SeuratObject)
})
.a   <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
RESUME <- "--resume" %in% .a           # reuse the step-3 checkpoint if it is on disk
SMOKE  <- as.integer(.get("smoke", "0"))       # >0: subsample each dataset to this many cells
MIX_N  <- as.integer(.get("mix_n", if (SMOKE > 0) "3000" else "30000"))
SFX    <- if (SMOKE > 0) "_smoke" else ""      # a smoke run never touches a production path
set.seed(1234)
if (SMOKE > 0) message("** SMOKE RUN: ", SMOKE, " cells per dataset; outputs suffixed ", SFX, " **")

ROOT <- "/FAST/gr10634/gaozy/aml_niche_net"
BASE <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration"
IN   <- file.path(BASE, "01_input")
OUT  <- file.path(BASE, "02_integrated");  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CKPT <- file.path(BASE, paste0("02_integrated/ws_prepared_list", SFX, ".rds"))
TAB  <- file.path(ROOT, "results/tables/ws");  dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
FIG  <- file.path(ROOT, "results/figures/ws"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

CHK <- list()
chk <- function(name, expected, got) {
  # A check that ERRORS must not kill the run -- it reports FAIL and we carry on. Job 3641242
  # died here because vapply(sl, ncol, 1L) throws on Seurat v5 (ncol returns a double), which
  # threw while all.equal was lazily evaluating `got`, aborting after step 1 had completed.
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected),
                                        got = as.character(got), pass = ok)
  message(sprintf("  %-46s expect %-12s got %-12s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
  invisible(ok)
}
chk_n <- function(...) if (SMOKE == 0) chk(...) else invisible(TRUE)   # count checks: full runs only
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)

# both ship with Seurat; used only by the mixing check. Asserted so a missing one fails here
# rather than two hours into the run.
stopifnot(requireNamespace("RANN",  quietly = TRUE),
          requireNamespace("irlba", quietly = TRUE))

## ---- step 0: gene axis -------------------------------------------------------------------------
ts("[0] gene axis")
RAW <- fread(file.path(TAB, "ws_reintegration_common_genes.csv"))$gene
chk("intersection genes read", 32042L, length(RAW))
MT  <- grep("^MT-", RAW, value = TRUE)
GEN <- setdiff(RAW, MT)
chk("MT- genes dropped", 13L, length(MT))
chk("gene axis used", 32029L, length(GEN))
chk("drop accounts for every gene", TRUE, length(GEN) + length(MT) == length(RAW))
message("    dropped: ", paste(MT, collapse = ", "))

FILES <- sort(list.files(IN, "\\.rds$", full.names = TRUE))
chk("input objects found", 9L, length(FILES))

if (RESUME && file.exists(CKPT)) {
  ts("[1-3] resuming from the step-3 checkpoint (nothing is overwritten or removed)")
  z <- readRDS(CKPT); sl <- z$sl; features <- z$features; rm(z); gc(FALSE)
  chk("checkpoint datasets", 9L, length(sl))
  chk("checkpoint anchor features", 3000L, length(features))
} else {
  ## ---- step 1: per dataset normalize / HVF / cell cycle -----------------------------------------
  ts("[1] per dataset: NormalizeData -> FindVariableFeatures(5000) -> CellCycleScoring")
  s.genes <- cc.genes$s.genes; g2m.genes <- cc.genes$g2m.genes
  chk("cc s.genes on the axis", 42L, sum(s.genes %in% GEN))
  chk("cc g2m.genes on the axis", 52L, sum(g2m.genes %in% GEN))

  sl <- vector("list", length(FILES))
  names(sl) <- sub("\\.rds$", "", basename(FILES))
  for (i in seq_along(FILES)) {
    ds <- names(sl)[i]
    o  <- readRDS(FILES[i])
    if (SMOKE > 0 && ncol(o) > SMOKE) o <- subset(o, cells = sample(colnames(o), SMOKE))
    stopifnot(all(GEN %in% rownames(o)))              # the axis must be a subset of every dataset
    o  <- subset(o, features = GEN)
    DefaultAssay(o) <- "RNA"
    o  <- NormalizeData(o, verbose = FALSE)
    o  <- FindVariableFeatures(o, nfeatures = 5000, verbose = FALSE)
    o  <- CellCycleScoring(o, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
    sl[[ds]] <- o
    ts(sprintf("    %-15s %7d cells | %5d genes | pmt median %.2f", ds, ncol(o), nrow(o),
               median(o$percent.mt, na.rm = TRUE)))
    rm(o); gc(FALSE)
  }
  chk_n("cells after the gene subset", 895228, sum(vapply(sl, function(x) as.numeric(ncol(x)), numeric(1))))
  chk("every dataset on the same axis", TRUE, all(vapply(sl, function(x) as.numeric(nrow(x)), numeric(1)) == length(GEN)))

  ## ---- step 2: anchor features ------------------------------------------------------------------
  ts("[2] SelectIntegrationFeatures(nfeatures = 3000)")
  features <- SelectIntegrationFeatures(object.list = sl, nfeatures = 3000)
  chk("anchor features", 3000L, length(features))
  chk("no MT- gene became an anchor feature", 0L, sum(grepl("^MT-", features)))

  ## ---- step 3: per dataset scale + PCA ----------------------------------------------------------
  ts("[3] per dataset: ScaleData(regress percent.mt + S.Score + G2M.Score) -> RunPCA(50)")
  for (ds in names(sl)) {
    x <- sl[[ds]]
    DefaultAssay(x) <- "RNA"
    LayerData(x[["RNA"]], layer = "scale.data") <- NULL
    x <- ScaleData(x, features = features, assay = "RNA", verbose = FALSE,
                   vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
    x <- RunPCA(x, npcs = 50, features = features, verbose = FALSE)
    sl[[ds]] <- x; ts(sprintf("    %-15s scaled + PCA", ds)); rm(x); gc(FALSE)
  }
  nan_ds <- names(sl)[vapply(sl, function(x) anyNA(Embeddings(x, "pca")), TRUE)]
  chk("datasets with NA in their PCA", 0L, length(nan_ds))
  if (length(nan_ds)) message("    ** ", paste(nan_ds, collapse = ", "))

  ts("[3c] writing the step-3 checkpoint (kept, not cleaned up): ", CKPT)
  saveRDS(list(sl = sl, features = features, genes = GEN), CKPT)
}

## ---- mixing, BEFORE correction ------------------------------------------------------------------
# fraction of a cell's 30 nearest neighbours that come from a different dataset, divided by the
# fraction expected if neighbours were drawn at random (1 - that dataset's share). 1.0 = fully
# mixed, 0 = fully segregated. Measured on an uncorrected PCA of the same anchor features, so the
# before/after pair is comparable.
mixing <- function(emb, ds, k = 30L) {
  nn   <- RANN::nn2(emb, k = k + 1L)$nn.idx[, -1, drop = FALSE]
  frac <- 1 - rowMeans(matrix(ds[nn], nrow = nrow(nn)) == ds)
  share <- table(ds) / length(ds)
  frac / (1 - as.numeric(share[ds]))
}

ts("[4] mixing BEFORE correction (uncorrected PCA on the anchor features, ", MIX_N, " cells)")
sub_cells <- unlist(lapply(sl, function(x) {
  n <- max(1L, round(MIX_N * ncol(x) / sum(vapply(sl, function(y) as.numeric(ncol(y)), numeric(1)))))
  sample(colnames(x), min(n, ncol(x)))
}), use.names = FALSE)
pre <- do.call(cbind, lapply(sl, function(x)
  as.matrix(LayerData(x, assay = "RNA", layer = "data")[features, intersect(colnames(x), sub_cells), drop = FALSE])))
pre_ds <- rep(names(sl), vapply(sl, function(x) length(intersect(colnames(x), sub_cells)), numeric(1)))
vr  <- apply(pre, 1, var)
if (any(vr <= 0)) ts(sprintf("    %d of %d anchor features have zero variance in the subsample -> excluded",
                             sum(vr <= 0), length(vr)))
pre_pca <- irlba::prcomp_irlba(t(pre[vr > 0, , drop = FALSE]), n = 30, center = TRUE, scale. = TRUE)$x
mix_pre <- mixing(pre_pca, pre_ds)
ts(sprintf("    before: median %.3f  mean %.3f", median(mix_pre), mean(mix_pre)))
rm(pre, pre_pca); gc(FALSE)

## ---- step 5: anchors + integration --------------------------------------------------------------
ts("[5] FindIntegrationAnchors(reduction = 'rpca', dims = 1:30)  -- 9 datasets, 36 pairs")
anchors <- FindIntegrationAnchors(object.list = sl, anchor.features = features,
                                  reduction = "rpca", dims = 1:30, verbose = TRUE)
n_in <- sum(vapply(sl, function(x) as.numeric(ncol(x)), numeric(1))); rm(sl); gc(FALSE)

ts("[6] IntegrateData(dims = 1:30)  -- this is the memory peak")
obj <- IntegrateData(anchorset = anchors, dims = 1:30, verbose = TRUE)
rm(anchors); gc(FALSE)
DefaultAssay(obj) <- "integrated"
chk("cells survive integration", as.numeric(n_in), ncol(obj))
chk("integrated assay features", 3000L, nrow(obj[["integrated"]]))

ts("[7] ScaleData -> RunPCA(50) -> RunUMAP(dims = 1:30)")
obj <- ScaleData(obj, vars.to.regress = NULL, verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

# 02.1's duplicate guard, rewritten: the original subscripts with a character vector, which errors.
emb <- Embeddings(obj, reduction = "pca")
dup <- rownames(emb)[duplicated(as.data.frame(emb))]
chk("cells with a duplicate PCA embedding", 0L, length(dup))
if (length(dup)) {
  ts("    dropping ", length(dup), " duplicate-embedding cells")
  obj <- subset(obj, cells = setdiff(colnames(obj), dup))
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
}
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)

## ---- mixing, AFTER correction -------------------------------------------------------------------
ts("[8] mixing AFTER correction")
keep <- intersect(sub_cells, colnames(obj))
mix_post <- mixing(Embeddings(obj, "pca")[keep, 1:30], obj$Dataset_ID[keep])
ts(sprintf("    after : median %.3f  mean %.3f", median(mix_post), mean(mix_post)))
chk("mixing improved (median after > median before)", TRUE, median(mix_post) > median(mix_pre))

MIX <- rbind(
  data.table(stage = "before", Dataset_ID = pre_ds, mixing = as.numeric(mix_pre)),
  data.table(stage = "after",  Dataset_ID = as.character(obj$Dataset_ID[keep]), mixing = as.numeric(mix_post))
)[, .(cells = .N, median_mixing = round(median(mixing), 4), mean_mixing = round(mean(mixing), 4)),
  by = .(stage, Dataset_ID)]
fwrite(MIX[order(Dataset_ID, stage)], file.path(TAB, paste0("ws_rpca_mixing",  SFX, ".csv")))
print(MIX[order(Dataset_ID, stage)])

## ---- output -------------------------------------------------------------------------------------
ts("[9] checks and output")
cnt <- LayerData(obj, assay = "RNA", layer = "counts")
chk("RNA counts are still integers", TRUE, all(cnt@x == floor(cnt@x)))
chk("datasets in the output", 9L, uniqueN(obj$Dataset_ID))
chk_n("samples in the output", 182L, uniqueN(obj$Sample_Name))
rm(cnt); gc(FALSE)

pdf(file.path(FIG, paste0("ws_rpca_umap",    SFX, ".pdf")), width = 22, height = 6)
print(DimPlot(obj, reduction = "umap", group.by = "Dataset_ID", raster = TRUE, shuffle = TRUE) |
      DimPlot(obj, reduction = "umap", group.by = "SampleType", raster = TRUE, shuffle = TRUE) |
      DimPlot(obj, reduction = "umap", group.by = "Phase", raster = TRUE, shuffle = TRUE))
dev.off()

saveRDS(obj, file.path(OUT, paste0("ws_rpca_integrated", SFX, ".rds")))
C <- rbindlist(CHK)
fwrite(C, file.path(TAB, paste0("ws_rpca_checks",  SFX, ".csv")))
ts(sprintf("[10] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", file.path(OUT, paste0("ws_rpca_integrated", SFX, ".rds")))
