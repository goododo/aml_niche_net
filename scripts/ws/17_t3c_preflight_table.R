#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../02_seurat_objects/04_annotated/GSE185381/<sample>.rds, 00_curated_manifest.csv
# OUTPUT : results/tables/ws/ws_t3c_preflight_ligands.csv    A[L,p,bin]: per patient x sender bin
#          results/tables/ws/ws_t3c_preflight_receptors.csv  receptor pos/neg counts per patient x bin
#          results/tables/ws/ws_t3c_preflight_grid.csv       qualifying patients per (L, sender, receiver, R)
#          results/tables/ws/ws_t3c_preflight_checks.csv
# what it does: T3c-0 needs two cells of the GATE_W1.1 table that the original pre-flight
#   (scratchpad t3c_preflight.R, 2026-10-02) computed but never persisted -- only the §0A summary
#   survived, and it does not contain them: the receptor for IFNG -> LMPP_GMP (rule 2: the one with
#   most qualifying patients) and the sender for TNF -> Mono_DC after Mono_DC itself is excluded
#   (rule 1: best-detected compartment, same standard as the pre-flight). This re-runs the pre-flight
#   with its definitions copied verbatim, writes the FULL grid this time, and first reproduces the six
#   numbers already in §0A (33, 32, 32, 19, 28, 15) as a check that the re-run matches the original
#   before the two new cells are trusted.
#   Definitions (from the original script, unchanged):
#     eligible samples : GSE185381, AML, Diagnosis, BM, sorting viability_only|none
#     state unit       : hierarchy_bin, cells with high_error == 1 dropped
#     sender usable    : >= 30 cells of that bin in the patient              (5.4 step 0)
#     receptor usable  : >= 50 receptor-positive AND >= 50 negative cells in the receiver bin
#     "best sender"    : largest cross-patient IQR of mean ligand expression, >= 10 patients
#     qualifying patient for a combo: sender usable AND receptor usable in the same patient
# This script deletes nothing and does not modify the 04_annotated objects.
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix); library(data.table) })
ROOT <- "/FAST/gr10634/gaozy/aml_niche_net"; TAB <- file.path(ROOT, "results/tables/ws")
AN  <- "/LARGE1/gr10634/gaozy/aml_niche_net/02_seurat_objects/04_annotated/GSE185381"
LIG <- c("IFNG","TNF","TGFB1")
REC <- list(IFNG = c("IFNGR1","IFNGR2"), TNF = c("TNFRSF1A","TNFRSF1B"), TGFB1 = c("TGFBR1","TGFBR2"))
RECV_BINS <- c("LMPP_GMP","Mono_DC")
MIN_SEND <- 30L; MIN_POS <- 50L; MIN_NEG <- 50L; MIN_PAT <- 10L
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-58s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

ts("[1] eligible samples")
man  <- fread(file.path(ROOT, "results/tables/01_preprocess/00_curated_manifest.csv"))
elig <- man[dataset == "GSE185381" & disease == "AML" & timepoint == "Diagnosis" & tissue_r == "BM" &
            sorting_r %in% c("viability_only","none"), sample]
ts("    ", length(elig), " eligible samples (§15.2 item 12 asks whether this is 41 or 42)")

ts("[2] per patient x bin: ligand availability and receptor counts")
A <- list(); R <- list()
for (sm in elig) {
  f <- file.path(AN, paste0(sm, ".rds")); if (!file.exists(f)) next
  s <- readRDS(f); md <- s@meta.data
  b <- md$hierarchy_bin; b[!is.na(md$high_error) & md$high_error == 1] <- NA
  cts <- LayerData(s, assay = "RNA", layer = "counts"); tot <- Matrix::colSums(cts)
  gl <- intersect(LIG, rownames(cts)); gr <- intersect(unlist(REC), rownames(cts))
  for (bb in unique(b[!is.na(b)])) {
    j <- which(b == bb); n <- length(j)
    if (n >= MIN_SEND && length(gl)) {
      m <- cts[gl, j, drop = FALSE]
      A[[length(A)+1]] <- data.table(sample = sm, bin = bb, n = n, gene = gl,
        mean_expr = as.numeric(Matrix::rowSums(m) / sum(tot[j]) * 1e4),       # CPM-like, per bin
        det_frac  = as.numeric(Matrix::rowSums(m > 0) / n))
    }
    if (bb %in% RECV_BINS && n >= MIN_POS + MIN_NEG && length(gr)) {
      pos <- as.numeric(Matrix::rowSums(cts[gr, j, drop = FALSE] > 0))
      R[[length(R)+1]] <- data.table(sample = sm, bin = bb, n = n, gene = gr, pos = pos, neg = n - pos)
    }
  }
  rm(s, cts); gc(FALSE)
}
A <- rbindlist(A); R <- rbindlist(R)
R[, ok := pos >= MIN_POS & neg >= MIN_NEG]
fwrite(A, file.path(TAB, "ws_t3c_preflight_ligands.csv"))
fwrite(R, file.path(TAB, "ws_t3c_preflight_receptors.csv"))

ts("[3] sender ranking per ligand: cross-patient IQR of mean expression, >= 10 patients")
SEND <- A[, .(patients = .N, IQR = round(IQR(mean_expr), 4), med_det = round(median(det_frac), 3),
              pct_zero = round(100*mean(mean_expr == 0), 0)), by = .(gene, bin)][patients >= MIN_PAT][order(gene, -IQR)]
print(SEND)
best <- SEND[, .SD[1], by = gene]
chk("IFNG best sender is T_NK (as §0A)",     "T_NK",    best[gene == "IFNG", bin])
chk("TNF best sender is Mono_DC (as §0A)",   "Mono_DC", best[gene == "TNF",  bin])
chk("TGFB1 best sender is T_NK (as §0A)",    "T_NK",    best[gene == "TGFB1", bin])
chk("IFNG T_NK patients = 36 (as §0A)",      36L, best[gene == "IFNG", patients])
chk("TNF Mono_DC patients = 40 (as §0A)",    40L, best[gene == "TNF",  patients])

ts("[4] the full grid: qualifying patients per (ligand, sender, receiver, receptor)")
GRID <- rbindlist(lapply(LIG, function(L) {
  snd <- SEND[gene == L, bin]
  rbindlist(lapply(snd, function(sb) {
    sp <- A[gene == L & bin == sb, sample]                      # sender usable in these patients
    rbindlist(lapply(RECV_BINS, function(rb) {
      if (rb == sb) return(NULL)                                 # rule 1: sender != receiver
      rbindlist(lapply(REC[[L]], function(rg) {
        rp <- R[gene == rg & bin == rb & ok == TRUE, sample]     # receptor usable in these patients
        data.table(ligand = L, sender = sb, receiver = rb, receptor = rg,
                   sender_IQR = SEND[gene == L & bin == sb, IQR],
                   qualifying_patients = length(intersect(sp, rp)))
      }))
    }))
  }))
}))
setorder(GRID, ligand, receiver, -qualifying_patients)
fwrite(GRID, file.path(TAB, "ws_t3c_preflight_grid.csv"))
print(GRID)

ts("[5] reproduce the six §0A numbers before trusting the two new cells")
g <- function(L, s, r, rg) GRID[ligand == L & sender == s & receiver == r & receptor == rg, qualifying_patients]
chk("IFNG  T_NK    -> Mono_DC  IFNGR2  = 33", 33L, g("IFNG","T_NK","Mono_DC","IFNGR2"))
chk("IFNG  T_NK    -> Mono_DC  IFNGR1  = 32", 32L, g("IFNG","T_NK","Mono_DC","IFNGR1"))
chk("TNF   Mono_DC -> LMPP_GMP TNFRSF1A = 19", 19L, g("TNF","Mono_DC","LMPP_GMP","TNFRSF1A"))
chk("TGFB1 T_NK    -> Mono_DC  TGFBR2  = 28", 28L, g("TGFB1","T_NK","Mono_DC","TGFBR2"))
chk("TGFB1 T_NK    -> LMPP_GMP TGFBR2  = 15", 15L, g("TGFB1","T_NK","LMPP_GMP","TGFBR2"))

ts("[6] the two cells T3c-0 must fill")
r2 <- GRID[ligand == "IFNG" & sender == "T_NK" & receiver == "LMPP_GMP"][order(-qualifying_patients)][1]
cat(sprintf("\n  ROW 2  IFNG | T_NK -> LMPP_GMP | receptor by rule 2 = %s | qualifying patients = %d\n",
            r2$receptor, r2$qualifying_patients))
snd4 <- SEND[gene == "TNF" & bin != "Mono_DC"][order(-IQR)][1]
r4 <- GRID[ligand == "TNF" & sender == snd4$bin & receiver == "Mono_DC" & receptor == "TNFRSF1B"]
cat(sprintf("  ROW 4  TNF  | best sender excluding Mono_DC = %s (IQR %.4f, %d patients) -> Mono_DC | TNFRSF1B | qualifying patients = %d\n",
            snd4$bin, snd4$IQR, snd4$patients, if (nrow(r4)) r4$qualifying_patients else NA_integer_))
cat(sprintf("  rule: combos with < %d qualifying patients are not run\n\n", MIN_PAT))
FILL <- data.table(row = c(2L, 4L),
  ligand = c("IFNG","TNF"), sender = c("T_NK", snd4$bin), receiver = c("LMPP_GMP","Mono_DC"),
  receptor = c(r2$receptor, "TNFRSF1B"),
  qualifying_patients = c(r2$qualifying_patients, if (nrow(r4)) r4$qualifying_patients else NA_integer_))
FILL[, runnable := !is.na(qualifying_patients) & qualifying_patients >= MIN_PAT]
fwrite(FILL, file.path(TAB, "ws_t3c_preflight_fill.csv")); print(FILL)

C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_t3c_preflight_checks.csv"))
ts(sprintf("[7] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done]")
