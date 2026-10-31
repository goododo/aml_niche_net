#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../01_processed_counts/rds/GSE185381.rds (ingest object: every deposited barcode, counts layers);
#          results/tables/00_ingest/GSE185381_qc_percell.csv.gz (Sample, Demuxed and QC metrics per barcode);
#          results/tables/01_preprocess/03_qc_report__GSE185381.csv (per-sample QC thresholds, n_doublet);
#          /LARGE1/.../02_seurat_objects/04_annotated/GSE185381/<sample>.rds of the 25 S5 T samples (cell names only);
#          ws_pb_GSE185381.rds ($cells: kept T flags); results/tables/ws/ws_s5_scores.csv (the 25 T samples);
#          GEO RAW matrix/features/barcodes of the six P02/P03 lanes (GSM5613747-52)
# OUTPUT : /LARGE1/.../ws_a1_batch/w4_doublet_scores.csv.gz  cell, sample, library, kept_T, score (the 25 samples' 04_annotated cells)
#          /LARGE1/.../ws_a1_batch/w4_libraries.csv          per library: labelled barcodes, warnings, K5 (QC-removed doublets vs kept cells)
#          /LARGE1/.../ws_a1_batch/w4_k5_samples.csv         per sample: reconstructed doublets vs the QC report's n_doublet
#          /LARGE1/.../ws_a1_batch/d9_soup_profile.csv.gz     gene, counts summed over T-ALL/ETP barcodes (top share >= 0.8), six lanes
#          /LARGE1/.../ws_a1_batch/d9_lanes.csv               per lane: called barcodes, T-ALL/ETP top-hashtag barcodes (K6), used, share
# what it does: PREREGISTRATION_A1_batch.md (tag prereg-a1-batch-v1) W4 (doublet scores) and D9 (soup profile), nothing
#   else. Reads no expression of any paediatric, post-treatment or relapse sample beyond using labelled barcodes of the same
#   libraries as scoring context; writes scores only for the 25 S5 T samples' cells.
#   usage: Rscript scripts/ws/38_a1_batch_cells.R [--smoke]   (--smoke: 2 libraries and 1 lane; outputs end in _smoke)
# This script deletes nothing and does not modify its inputs.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here) })
SMOKE <- "--smoke" %in% commandArgs(TRUE)
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R"))
TABWS <- file.path(TAB_DIR, "ws"); BDIR <- file.path(LARGE1_DIR, "ws_a1_batch"); dir.create(BDIR, showWarnings = FALSE)
SFX <- if (SMOKE) "_smoke" else ""
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
wfile <- function(name) {   # the only way this script names an output: a new file under ws_a1_batch/, never overwriting a real one
  f <- file.path(BDIR, sub("(\\.csv(\\.gz)?)$", paste0(SFX, "\\1"), name))
  if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f)
  f }
CHK <- list(); chk <- function(name, expected, got) { ok <- isTRUE(all.equal(expected, got))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-72s expect %-14s got %-14s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }

## -- 1. the 25 samples, their 04_annotated cells, kept T cells ---------------------------------------------------------
S25 <- fread(file.path(TABWS, "ws_s5_scores.csv"))[compartment == "T", sample]
AN_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "04_annotated", "GSE185381")
ANN <- rbindlist(lapply(S25, function(s) data.table(cell = colnames(readRDS(file.path(AN_DIR, paste0(s, ".rds")))), sample = s)))
ANN[, library := sub("_.*$", "", cell)]                       # cell names start with the GSM accession
PB <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"))
KT <- PB$cells[kept == TRUE & group4 == "T" & sample %in% S25, cell]
ANN[, kept_T := cell %in% KT]
chk("25 samples", 25L, length(S25))
chk("kept T cells of the 25 (S4)", 14609L, length(KT))
chk("every kept T cell is an 04_annotated cell", TRUE, all(KT %in% ANN$cell))
LIBS <- sort(unique(ANN$library))
chk("libraries holding the 25 samples' cells", 35L, length(LIBS))
ts(sprintf("[1] %d 04_annotated cells of the 25 samples in %d libraries; %d kept T cells", nrow(ANN), length(LIBS), sum(ANN$kept_T)))

## -- 2. W4: scDblFinder per library on the author-labelled barcodes --------------------------------------------------------
QCC <- fread(file.path(TAB_DIR, "00_ingest", "GSE185381_qc_percell.csv.gz"),
             select = c("cell_barcode", "Sample", "Demuxed", "GSM", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hb"))
setnames(QCC, c("cell", "sample", "demuxed", "library", "nc", "nf", "mt", "hb"))
QCR <- fread(file.path(TAB_DIR, "01_preprocess", "03_qc_report__GSE185381.csv"))
LIB_RUN <- if (SMOKE) LIBS[1:2] else LIBS
ts("[2] loading the ingest object (all deposited barcodes)")
obj <- readRDS(file.path(LARGE1_DIR, "01_processed_counts", "rds", "GSE185381.rds"))
assay <- obj[["RNA"]]; LAY <- SeuratObject::Layers(assay, search = "counts")
ts(sprintf("    %d cells, %d counts layers", ncol(obj), length(LAY)))
cell_layer <- rbindlist(lapply(LAY, function(l) data.table(cell = SeuratObject::Cells(assay, layer = l), layer = l)))
get_counts <- function(cells) {   # counts of these cells from whichever layer(s) hold them, in the given order
  m <- cell_layer[cell %in% cells]; stopifnot(nrow(m) == length(cells))
  M <- do.call(cbind, lapply(unique(m$layer), function(l) SeuratObject::LayerData(assay, layer = l)[, m[layer == l, cell], drop = FALSE]))
  M[, cells, drop = FALSE] }
pass_qc <- function(d, s) {   # 03_per_sample_qc thresholds of sample s (effective MAD + absolute bounds, from the QC report)
  q <- QCR[Sample == s]; stopifnot(nrow(q) == 1)
  d$nc >= q$ncount_lo & d$nc <= q$ncount_hi & d$nf >= q$nfeat_lo & d$nf <= q$nfeat_hi & d$mt <= q$mt_hi & (is.na(q$hb_hi) | d$hb <= q$hb_hi) }
SC <- list(); LIBTAB <- list(); K5S <- list()
for (L in LIB_RUN) {
  lab <- QCC[library == L & demuxed == TRUE, cell]
  M <- get_counts(lab)
  warns <- character(0); msgs <- character(0)
  set.seed(SEED)
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = M))
  sce <- withCallingHandlers(scDblFinder::scDblFinder(sce),     # default arguments; in this environment the score is the pre-training heuristic
           warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
           message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  sc <- data.table(cell = colnames(sce), score = sce$scDblFinder.score)
  # K5: labelled barcodes of the 25 samples that pass their sample's QC thresholds but are absent from 04_annotated = union doublets removed at QC
  q <- QCC[library == L & demuxed == TRUE & sample %in% S25]
  q[, passes := pass_qc(.SD, sample[1]), by = sample]
  q[is.na(passes), passes := FALSE]                             # a missing QC metric counts as failing QC, as at QC
  q[, removed_dbl := passes & !(cell %in% ANN$cell)]
  k <- merge(q[, .(cell, sample, removed_dbl, in_ann = cell %in% ANN$cell)], sc, by = "cell")
  a <- k[removed_dbl == TRUE, score]; b <- k[in_ann == TRUE, score]
  auc <- if (length(a) && length(b)) mean(outer(a, b, ">") + 0.5 * outer(a, b, "==")) else NA_real_
  K5S[[L]] <- k[, .(n_removed_dbl = sum(removed_dbl)), by = sample][, library := L]
  LIBTAB[[L]] <- data.table(library = L, n_labelled = length(lab), n_scored = nrow(sc), n_cells_25 = ANN[library == L, .N],
                            n_removed_dbl = length(a), k5_auc = auc, n_warnings = length(warns),
                            xgboost_mentioned = any(grepl("xgb", c(warns, msgs), ignore.case = TRUE)),
                            first_warning = if (length(warns)) substr(warns[1], 1, 160) else "")
  SC[[L]] <- merge(ANN[library == L], sc, by = "cell", all.x = TRUE)
  ts(sprintf("    %s: %d labelled barcodes, %d cells of the 25, K5 AUC %.3f (%d removed doublets), %d warnings", L, length(lab), ANN[library == L, .N], auc, length(a), length(warns)))
  rm(M, sce); gc(FALSE)
}
SCORES <- rbindlist(SC)
chk("every 04_annotated cell of the scored libraries has a score", 0L, SCORES[is.na(score), .N])
LIBTAB <- rbindlist(LIBTAB)
K5S <- rbindlist(K5S)[, .(n_removed_dbl = sum(n_removed_dbl)), by = sample]
K5S <- merge(K5S, QCR[, .(sample = Sample, n_doublet)], by = "sample", all.x = TRUE)
if (!SMOKE) chk("K5: median per-library AUC (QC-removed doublets vs kept cells) > 0.5", TRUE, isTRUE(median(LIBTAB$k5_auc, na.rm = TRUE) > 0.5))
fwrite(SCORES[, .(cell, sample, library, kept_T, score)], wfile("w4_doublet_scores.csv.gz"))
fwrite(LIBTAB, wfile("w4_libraries.csv")); fwrite(K5S, wfile("w4_k5_samples.csv"))
print(LIBTAB[, .(library, n_labelled, n_cells_25, n_removed_dbl, k5_auc, n_warnings, xgboost_mentioned)])
print(K5S)
rm(obj, assay); gc(FALSE)

## -- 3. D9: T-ALL/ETP soup profile from the six P02/P03 lanes (GEO matrices; HTO rows for the hashtags) -------------------
RAW <- file.path(RAW_DIR, "public", "GEO", "GSE185381_RAW")
LANES <- c("GSM5613747_2019-08-01-count-1", "GSM5613748_2019-08-01-count-2", "GSM5613749_2019-08-01-count-3",
           "GSM5613750_2019-08-01-count-4", "GSM5613751_2019-08-29-count-1", "GSM5613752_2019-08-29-count-2")
K6_TOP <- c(403L, 2757L, 2408L, 2558L, 1911L, 1924L)          # 2026-10-28 HTO tally (VERIFICATION_G03-G08 section 1.2)
TALL <- c("TALL3783", "TALL3862", "ETP9153", "TALL4954")
read_lane <- function(pre) {   # 10x matrix.mtx (genes + antibody rows) -> sparse matrix; gene names made unique as Read10X does
  ft <- fread(paste0(pre, "_features.tsv.gz"), header = FALSE); bc <- fread(paste0(pre, "_barcodes.tsv.gz"), header = FALSE)$V1
  m <- fread(cmd = paste("zcat", shQuote(paste0(pre, "_matrix.mtx.gz")), "| grep -v '^%'"), header = FALSE)
  M <- sparseMatrix(i = m$V1[-1], j = m$V2[-1], x = m$V3[-1], dims = c(m$V1[1], m$V2[1]))
  list(M = M, ft = ft, bc = bc) }
LANE_RUN <- if (SMOKE) LANES[1] else LANES
PROF <- NULL; LT <- list()
for (k in seq_along(LANE_RUN)) {
  x <- read_lane(file.path(RAW, LANE_RUN[k])); ge <- x$ft$V3 == "Gene Expression"
  hto <- which(x$ft$V3 == "Antibody Capture" & grepl("^HTO", x$ft$V1))
  H <- as.matrix(x$M[hto, , drop = FALSE]); rownames(H) <- sub("^HTO-", "", x$ft$V1[hto])
  top <- rownames(H)[max.col(t(H), ties.method = "first")]; tot <- colSums(H); share <- apply(H, 2, max) / pmax(tot, 1)
  is_tall <- top %in% TALL; use <- is_tall & share >= 0.8
  genes <- make.unique(x$ft$V2[ge])
  p <- Matrix::rowSums(x$M[ge, use, drop = FALSE]); names(p) <- genes
  PROF <- if (is.null(PROF)) p else { stopifnot(identical(names(PROF), genes)); PROF + p }
  LT[[k]] <- data.table(lane = LANE_RUN[k], n_called = length(x$bc), n_top_tall = sum(is_tall), n_used = sum(use), tall_share = sum(is_tall) / length(x$bc))
  ts(sprintf("    %s: %d called barcodes, %d with a T-ALL/ETP top hashtag (%.1f %%), %d used (top share >= 0.8)", LANE_RUN[k], length(x$bc), sum(is_tall), 100 * sum(is_tall) / length(x$bc), sum(use)))
}
LT <- rbindlist(LT)
chk("D9: soup-profile gene axis equals the pseudobulk gene axis", TRUE, identical(names(PROF), PB$genes))
if (!SMOKE) chk("K6: T-ALL/ETP top-hashtag barcodes per lane equal the 2026-10-28 tally", K6_TOP, LT$n_top_tall)
# K6b (prereg section 11, replacing K6's top-1 % marker check): the soup profile is closest to healthy T among healthy T, B, monocytes
H10 <- fread(file.path(TABWS, "ws_s5_scores.csv"))[compartment == "T" & status == "healthy", sample]
cpm1 <- function(v) v / sum(v) * 1e6
REFP <- sapply(c(T = "T", B = "B", Mono = "Mono"), function(g) Matrix::rowSums(PB$raw[[g]][, intersect(H10, colnames(PB$raw[[g]])), drop = FALSE]))
ALLP <- cbind(soup = cpm1(PROF[PB$genes]), apply(REFP[PB$genes, ], 2, cpm1))
keep6 <- !grepl("^(RP[LS][0-9]|MT-)", rownames(ALLP)) & apply(ALLP >= 10, 1, any)        # housekeeping genes removed; expressed somewhere
RHO <- cor(log1p(ALLP[keep6, ]), method = "spearman")["soup", c("T", "B", "Mono")]
ts(sprintf("    K6b: Spearman of the soup profile with healthy T %.3f, B %.3f, Mono %.3f (%d genes)", RHO["T"], RHO["B"], RHO["Mono"], sum(keep6)))
chk("K6b: the soup profile correlates best with healthy T (vs B, monocytes)", "T", names(which.max(RHO)))
fwrite(data.table(gene = names(PROF), counts = as.numeric(PROF)), wfile("d9_soup_profile.csv.gz"))
fwrite(cbind(LT, rho_T = RHO["T"], rho_B = RHO["B"], rho_Mono = RHO["Mono"]), wfile("d9_lanes.csv"))
print(LT)

C <- rbindlist(CHK); fwrite(C, wfile("cells_checks.csv"))
ts(sprintf("[done] %d checks, %d FAIL", nrow(C), sum(!C$pass)))
if (any(!C[!grepl("^K6b", check), pass])) quit(status = 1)                       # K6b concerns D9 only (prereg section 11)
