#!/usr/bin/env Rscript
# 01_extract_barcodes.R ----
# Doublet-pairing SMOKE, stage 1: pull the smoke lanes out of the GSE185381 ingest object and put every deposited barcode
# into class A / B_union / B_mad / C (PROPOSAL_doublet_pairing_external_v1.3 section 4.1), then apply the line QC to B and C.
# INPUT  : RDS_INGEST_DIR/GSE185381.rds                         every deposited barcode, one counts layer per lane
#          TAB_DIR/00_ingest/GSE185381_qc_percell.csv.gz         Sample, Demuxed, nCount/nFeature/percent.mt/.hb per barcode
#          TAB_DIR/01_preprocess/03_qc_report__GSE185381.csv     per-donor effective QC thresholds and n_doublet
#          LARGE1/02_seurat_objects/04_annotated/GSE185381/<donor>.rds   class A = these cells; their projection labels
#          TAB_DIR/ws/ws_percell_labels.csv.gz                   ws atlas cluster (cl) of class-A cells
#          GEO_RAW_DIR/GSE185381_RAW/<lane>_{matrix.mtx,features.tsv,barcodes.tsv}.gz   hashtag rows (hashed lanes: pooled, P01)
# OUTPUT : DBL_LANE_DIR/<lane>.rds  list(lane, counts = genes x analysis barcodes (A + line-QC-passing B and C),
#                                        bc = one row per deposited barcode, ann = class-A labels, hto = hashtags x barcodes
#                                        or NULL, thr = line-QC thresholds used)
#          DBL_TAB_DIR/dbl_barcode_classes_smoke.csv   lane x donor x class: deposited, passing line QC
#          DBL_TAB_DIR/dbl_checks_01_smoke.csv          known answers and invariants (PASS/FAIL)
#          --lanes=p01_2,p01_3 (stage 0a, P01 count-2/3): the same two tables with suffix _p01 instead of _smoke
#          --share --dataset=<ds> (share design Stage 1 S1b): every analysis lane of <ds> in scripts/dbl/dbl_lanes_share.tsv;
#            tables suffixed _share_<ds>; a lane object that exists is not re-extracted, its stored classes enter the table.
#            GSE289435 / Petti2019: per-barcode metadata from the ingest object (one sample per lane, all author-deposited,
#            so no class C), no hashtag rows.
# Usage  : Rscript scripts/dbl/01_extract_barcodes.R [--lanes=single,pooled | --lanes=p01_2,p01_3 | --share --dataset=<ds>] [--force]
# Nothing here types a doublet, so nothing needs masking. Deletes nothing; refuses to overwrite a lane object without --force.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
ARGS  <- commandArgs(TRUE)
FORCE <- "--force" %in% ARGS
la    <- sub("^--lanes=", "", grep("^--lanes=", ARGS, value = TRUE))
SHARE <- "--share" %in% ARGS
DS    <- if (SHARE) sub("^--dataset=", "", grep("^--dataset=", ARGS, value = TRUE)) else DBL_DATASET
if (SHARE) { stopifnot(length(DS) == 1, !length(la))
  LT <- fread(here::here("scripts", "dbl", "dbl_lanes_share.tsv"))[dataset == DS & role == "analysis"]
  LANES <- setNames(sort(unique(LT$lane)), sort(unique(LT$lane))); SFX <- paste0("_share_", DS)
} else {
KEYS  <- c(DBL_SMOKE_LANES, DBL_P01_LANES[c("p01_2", "p01_3")])   # p01_1 IS "pooled": never a key, never re-extracted
LANES <- if (length(la)) KEYS[strsplit(la, ",")[[1]]] else DBL_SMOKE_LANES
P01   <- !all(LANES %in% DBL_SMOKE_LANES)
if (P01) stopifnot(setequal(names(LANES), c("p01_2", "p01_3")), length(LANES) == 2)   # P01 lanes only both, never mixed
SFX   <- if (P01) "_p01" else "_smoke" }                           # a P01 run never touches the _smoke tables
stopifnot(length(LANES) > 0, !anyNA(LANES))                        # an unknown key (e.g. p01_1) indexes to NA and stops here
G185  <- DS == "GSE185381"                                         # the only dataset with Library / Demuxed / hashtag rows
for (d in c(DBL_LANE_DIR, DBL_TAB_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {   # record and print one check; a FAIL never stops the run, it is read at the end
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-78s expect %-12s got %-12s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }

## -- [1] per-barcode metadata and QC thresholds ----
ts("[1] per-barcode metadata and QC thresholds")
GSMS <- sub("_.*$", "", LANES)                                 # lane names start with the GSM accession (Petti: the sample id)
if (G185) {
QCC <- fread(file.path(TAB_DIR, "00_ingest", "GSE185381_qc_percell.csv.gz"),
             select = c("cell_barcode", "Sample", "Demuxed", "GSM", "Library", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hb"))
setnames(QCC, c("cell", "donor", "demuxed", "gsm", "lane", "nc", "nf", "mt", "hb"))
QCC <- QCC[gsm %in% GSMS] }
QCR <- fread(file.path(TAB_DIR, "01_preprocess", paste0("03_qc_report__", DS, ".csv")))[, Sample := as.character(Sample)]   # Petti ids are numeric
pass_full <- function(d, s) {   # 03_per_sample_qc thresholds of donor s (effective MAD + absolute bounds) -- copied from ws/38
  q <- QCR[Sample == s]; stopifnot(nrow(q) == 1)
  d$nc >= q$ncount_lo & d$nc <= q$ncount_hi & d$nf >= q$nfeat_lo & d$nf <= q$nfeat_hi & d$mt <= q$mt_hi & (is.na(q$hb_hi) | d$hb <= q$hb_hi) }
if (G185) for (L in LANES) chk(sprintf("%s: lane name agrees with the Library column", L), TRUE, all(QCC[gsm == sub("_.*$", "", L), unique(lane)] == L))

## -- [2] ingest object: the lanes' counts layers ----
ts("[2] loading the ingest object (all deposited barcodes)")
obj   <- readRDS(file.path(RDS_INGEST_DIR, paste0(DS, ".rds")))
assay <- obj[["RNA"]]
LAY   <- SeuratObject::Layers(assay, search = "counts")
if (G185) { META <- as.data.table(obj@meta.data, keep.rownames = "cell")[, .(cell, Sample, Demuxed, Library)]
} else {   # one sample per lane, all barcodes author-deposited: per-barcode metadata straight from the ingest object
  META <- as.data.table(obj@meta.data, keep.rownames = "cell")[, .(cell, Sample = as.character(Sample), Demuxed = TRUE, Library = as.character(Sample))]
  QCC  <- as.data.table(obj@meta.data, keep.rownames = "cell")[, .(cell, donor = as.character(Sample), demuxed = TRUE, gsm = sub("_.*$", "", Sample),
            lane = as.character(Sample), nc = nCount_RNA, nf = nFeature_RNA, mt = percent.mt, hb = percent.hb)][lane %in% LANES] }
cell_layer <- rbindlist(lapply(LAY, function(l) data.table(cell = SeuratObject::Cells(assay, layer = l), layer = l)))
GENES <- rownames(assay)
ts(sprintf("    %d barcodes, %d genes, %d counts layers", ncol(obj), length(GENES), length(LAY)))
get_counts <- function(cells) {   # counts of these cells from whichever layer(s) hold them, in the given order (as ws/38)
  m <- cell_layer[cell %in% cells]; stopifnot(nrow(m) == length(cells))
  M <- do.call(cbind, lapply(unique(m$layer), function(l) SeuratObject::LayerData(assay, layer = l)[, m[layer == l, cell], drop = FALSE]))
  M <- M[, cells, drop = FALSE]; rownames(M) <- GENES; M }

## -- [3] class-A labels: 04_annotated meta + ws atlas cluster ----
ts("[3] class-A cells and their labels (04_annotated) + ws atlas cluster")
AN_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "04_annotated", DS)
ATL <- fread(file.path(TAB_DIR, "ws", "ws_percell_labels.csv.gz"), select = c("cell", "cl"))
ATL <- ATL[sub("_.*$", "", cell) %in% GSMS]
CLASS_TAB <- list()
for (L in LANES) {
  out <- file.path(DBL_LANE_DIR, paste0(L, ".rds"))
  if (file.exists(out) && !FORCE) { ts("    [skip] ", out, " exists (use --force)")
    if (SHARE) { b0 <- readRDS(out)$bc                              # share tables cover every analysis lane, extracted earlier or now
      CLASS_TAB[[L]] <- b0[, .(deposited = .N, pass_line_qc = sum(pass_line_qc)), by = .(lane, donor, class)][order(donor, class)]; rm(b0) }
    next }
  ts("[4] lane ", L)
  bc <- merge(META[Library == L, .(cell, donor_meta = Sample, demuxed_meta = Demuxed)], QCC[lane == L], by = "cell", all = TRUE)
  chk(sprintf("%s: ingest barcodes == qc_percell barcodes", L), TRUE, bc[, !anyNA(donor_meta) && !anyNA(donor)])
  chk(sprintf("%s: Sample and Demuxed agree between ingest meta and qc_percell", L), TRUE,
      bc[, all(donor_meta == donor) && all(demuxed_meta == demuxed)])
  bc[, c("donor_meta", "demuxed_meta") := NULL]
  donors <- sort(unique(bc[demuxed == TRUE, donor]))
  chk(sprintf("%s: undemultiplexed barcodes carry UNDEMUX__<lane>", L), TRUE, bc[demuxed == FALSE, all(donor == paste0("UNDEMUX__", L))])
  ## class-A labels from each donor's 04_annotated object (only this lane's cells)
  ANN <- rbindlist(lapply(donors, function(s) {
    o <- readRDS(file.path(AN_DIR, paste0(s, ".rds")))
    md <- as.data.table(o@meta.data, keep.rownames = "cell")
    md[startsWith(cell, paste0(L, "_")), .(cell, donor = s, bmm_fine, bmm_prob, high_error, anno_source, hierarchy_bin)] }))
  ANN <- merge(ANN, ATL, by = "cell", all.x = TRUE)
  setnames(ANN, "cl", "atlas_cl")
  ## classes
  bc[, class := fifelse(cell %in% ANN$cell, "A", fifelse(demuxed == FALSE, "C", NA_character_))]
  for (s in donors) { i <- bc[, which(is.na(class) & donor == s)]
    pf <- pass_full(bc[i], s); pf[is.na(pf)] <- FALSE          # a missing QC metric counts as failing QC, as at QC
    bc[i, class := fifelse(pf, "B_union", "B_mad")] }
  chk(sprintf("%s: every barcode has exactly one class", L), 0L, bc[is.na(class), .N])
  chk(sprintf("%s: every class-A cell is author-assigned", L), TRUE, bc[class == "A", all(demuxed)])
  chk(sprintf("%s: every 04_annotated cell of the lane is a deposited barcode", L), TRUE, all(ANN$cell %in% bc$cell))
  chk(sprintf("%s: class-A donors agree with 04_annotated donors", L), TRUE, all(bc[class == "A"][ANN, on = "cell", donor == i.donor]))
  ## known answer: B_union per donor never exceeds the QC report's n_doublet (equality for a donor captured in this lane only)
  bu <- bc[class == "B_union", .N, by = donor][QCR[, .(donor = Sample, n_doublet, n_raw)], on = "donor", nomatch = NULL]
  n_dep <- bc[demuxed == TRUE, .N, by = donor]
  for (s in donors) {
    nb <- bu[donor == s, if (.N) N else 0L]; q <- QCR[Sample == s]
    if (n_dep[donor == s, N] == q$n_raw) chk(sprintf("%s: %s B_union == QC report n_doublet (donor only in this lane)", L, s), q$n_doublet, nb)
    else chk(sprintf("%s: %s B_union <= QC report n_doublet (donor spans lanes)", L, s), TRUE, nb <= q$n_doublet) }
  ## line QC for B and C: lower tails of nCount/nFeature, upper tail of percent.mt  [H 4.4.3]
  thr <- QCR[Sample %in% donors, .(donor = Sample, ncount_lo, nfeat_lo, mt_hi)]
  thr_c <- thr[, .(donor = paste0("UNDEMUX__", L), ncount_lo = median(ncount_lo), nfeat_lo = median(nfeat_lo), mt_hi = median(mt_hi))]
  thr <- rbind(thr, thr_c)
  bc <- merge(bc, thr, by = "donor", all.x = TRUE)
  bc[, pass_line_qc := class == "A" | (!is.na(nc) & nc >= ncount_lo & nf >= nfeat_lo & !is.na(mt) & mt <= mt_hi)]
  bc[, c("ncount_lo", "nfeat_lo", "mt_hi") := NULL]
  ## counts: class A + line-QC-passing B and C
  keep <- bc[pass_line_qc == TRUE, cell]
  M <- get_counts(keep)
  chk(sprintf("%s: qc_percell nCount == column sums of the ingest counts (analysis barcodes)", L), TRUE,
      isTRUE(all.equal(as.numeric(Matrix::colSums(M)), bc[match(keep, cell), as.numeric(nc)])))
  ## hashtags: Antibody Capture rows of the raw GEO matrix (pooled lane); ingest dropped them
  HTO <- NULL
  pre <- file.path(GEO_RAW_DIR, "GSE185381_RAW", L)
  hrows <- integer(0)
  if (G185) { ft <- fread(paste0(pre, "_features.tsv.gz"), header = FALSE)
    hrows <- which(ft$V3 == "Antibody Capture" & grepl(DBL_HTO_PREFIX, ft$V1)) }
  if (length(hrows)) {
    rbc <- fread(paste0(pre, "_barcodes.tsv.gz"), header = FALSE)$V1
    m <- fread(cmd = paste("zcat", shQuote(paste0(pre, "_matrix.mtx.gz")), "| grep -v '^%'"), header = FALSE)
    RAW <- sparseMatrix(i = m$V1[-1], j = m$V2[-1], x = m$V3[-1], dims = c(m$V1[1], m$V2[1])); rm(m)
    colnames(RAW) <- paste0(L, "_", rbc)
    chk(sprintf("%s: raw GEO barcodes == ingest barcodes", L), TRUE, setequal(colnames(RAW), bc$cell))
    ge <- which(ft$V3 == "Gene Expression")
    rn <- make.unique(ft$V2[ge])                               # gene symbols made unique as Read10X does
    common <- intersect(rn, GENES)
    chk(sprintf("%s: every ingest gene is a raw GEO gene (by name)", L), length(GENES), length(common))
    probe <- head(keep, 50)                                    # same deposit: raw gene counts equal ingest counts, gene by gene
    chk(sprintf("%s: raw GEO gene counts == ingest counts (50 barcodes, matched by name)", L), TRUE,
        isTRUE(all.equal(unname(as.matrix(RAW[ge[match(common, rn)], probe])), unname(as.matrix(M[common, probe])))))
    HTO <- RAW[hrows, , drop = FALSE]; rownames(HTO) <- ft$V1[hrows]; rm(RAW); gc(FALSE)
    ts(sprintf("    %d hashtag rows: %s", nrow(HTO), paste(rownames(HTO), collapse = ", ")))
  } else ts("    no hashtag rows in this lane")
  ## save + class table
  bc[, in_counts := cell %in% keep]
  saveRDS(list(lane = L, counts = M, bc = bc, ann = ANN, hto = HTO, thr = thr), out)
  ts("    [write] ", out)
  CLASS_TAB[[L]] <- bc[, .(deposited = .N, pass_line_qc = sum(pass_line_qc)), by = .(lane, donor, class)][order(donor, class)]
  print(CLASS_TAB[[L]])
  print(bc[, .(deposited = .N, pass_line_qc = sum(pass_line_qc)), by = class][order(class)])
  rm(M, bc, ANN, HTO); gc(FALSE)
}
CHK <- rbindlist(CHK)
if (length(CLASS_TAB) == length(LANES)) {   # write only when every requested lane was extracted: a skip never shrinks a table
  fwrite(rbindlist(CLASS_TAB), file.path(DBL_TAB_DIR, paste0("dbl_barcode_classes", SFX, ".csv")))
  fwrite(CHK, file.path(DBL_TAB_DIR, paste0("dbl_checks_01", SFX, ".csv")))
} else ts("    [no table written] ", length(LANES) - length(CLASS_TAB), " requested lane(s) skipped (use --force to redo all)")
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
