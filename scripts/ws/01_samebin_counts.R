#!/usr/bin/env Rscript
# 01_samebin_counts.R ----
# INPUT  : results/tables/01_preprocess/00_curated_manifest.csv
#          LARGE1/02_seurat_objects/04_annotated/<ds>/<sample>.rds  (hierarchy_bin, high_error)
#          CellChatDB.human (installed CellChat; no new dependency)
# OUTPUT : results/tables/ws/ws_w1_samebin_counts.csv        per (study, bin, ligand, receptor) patient counts
#          results/tables/ws/ws_w1_samebin_summary.csv       GATE_W1.2b verdict per study
#          results/tables/ws/ws_w1_samebin_ligands.csv       measurable-ligand census per patient x bin
# WHAT IT DOES : task T3d of HANDOFF_within_sample_v1.5, the GATE_W1.2b count. Counts the SAME-BIN
#          (ligand, receptor, bin) combinations that have enough cells, which the already-passing
#          GATE_W1.2 does not cover because that count pooled all sender compartments.
# Usage  : Rscript scripts/ws/01_samebin_counts.R
#
# THE RULE, frozen by the user on 2026-10-02 BEFORE this file existed (v1.5 section 5.6 + 15.2 #5):
#   ligand measurable in bin b of patient p : detected in 1% to 50% of that bin's cells
#   receivers                               : cells of bin b that do NOT detect the ligand
#   per patient                             : >= 50 receptor-positive AND >= 50 receptor-negative
#                                             receivers  (50 is primary; 30 is the pre-registered
#                                             fallback written into section 5.6 before any count)
#   combination counts for a study           : >= 10 patients satisfy all of the above
#   gate                                     : >= 30 combinations in >= 2 studies
# Nothing here chooses a threshold after seeing a number. The 30-cell arm is reported because the
# handoff already fixed it as the fallback, not because 50 failed.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(CellChat) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "..", "config", "config_paths.R"))
set.seed(SEED)

AN_DIR   <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated")   # config's PROJ_OBJ_DIR is stale
OUT_DIR  <- file.path(TAB_DIR, "ws"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LIG_LO   <- 0.01; LIG_HI <- 0.50          # "measurable" window, frozen
MIN_CELL <- 100L                          # a bin must allow 50 + 50 at all
MIN_PAT  <- 10L                           # patients per combination
ARMS     <- c(primary = 50L, fallback = 30L)
STUDIES  <- c("GSE185381", "GSE239721", "GSE116256", "GSE227903", "GSE289435")
MAIN_BINS <- c("LMPP_GMP", "Mono_DC")     # section 5.4 step 0; others reported but not the gate

## -- ligand-receptor resource -------------------------------------------------------------------
data("CellChatDB.human", package = "CellChat")
inter <- as.data.table(CellChatDB.human$interaction)[, .(ligand, receptor)]
cplx  <- CellChatDB.human$complex
expand <- function(r) if (r %in% rownames(cplx))
  unlist(strsplit(gsub(" ", "", paste(cplx[r, ], collapse = ",")), ",")) else r
pairs <- unique(inter[, .(gene = unlist(lapply(receptor, expand))), by = .(ligand, receptor)][gene != ""])
LIG <- sort(unique(pairs$ligand)); REC <- sort(unique(pairs$gene))
message(sprintf("[1] CellChatDB: %d interactions | %d ligands | %d receptor genes | %d (ligand,receptor-gene) pairs",
                nrow(inter), length(LIG), length(REC), uniqueN(pairs[, .(ligand, gene)])))

## -- cohort -------------------------------------------------------------------------------------
man <- fread(file.path(TAB_DIR, "01_preprocess", "00_curated_manifest.csv"))
elig <- man[disease == "AML" & timepoint == "Diagnosis" & tissue_r == "BM" &
            sorting_r %in% c("viability_only", "none") & dataset %in% STUDIES, .(dataset, sample)]
message(sprintf("[2] eligible diagnostic samples: %d across %d studies", nrow(elig), uniqueN(elig$dataset)))

## -- accumulate ---------------------------------------------------------------------------------
# The measurable-ligand set differs per patient, so nothing dense is accumulated across patients.
# Instead the receptor x ligand grid is evaluated per patient x bin and read ONLY at the positions
# of real CellChatDB pairs, then stacked long. Dimensions can therefore never disagree.
HIT <- list()                              # long: one row per (patient, bin, pair) that clears >=30
LIGCEN <- list()                           # measurable-ligand census
GL <- GR <- NULL; PIDX <- NULL             # gene axes + pair index, fixed PER DATASET
CUR_DS <- NA_character_                    # GSE116256 is Seq-Well: its gene universe differs
MISSING <- character()                     # manifest rows with no object on disk

elig <- elig[order(dataset, sample)]       # dataset-major, so the gene axes are fixed once per study
for (i in seq_len(nrow(elig))) {
  ds <- elig$dataset[i]; sm <- elig$sample[i]
  if (!identical(ds, CUR_DS)) { CUR_DS <- ds; GL <- GR <- NULL; PIDX <- NULL }
  f <- file.path(AN_DIR, ds, paste0(sm, ".rds"))
  if (!file.exists(f)) { MISSING <<- c(MISSING, paste0(ds, "/", sm)); next }
  s <- readRDS(f); md <- s@meta.data
  b <- md$hierarchy_bin
  if (!is.null(md$high_error)) b[!is.na(md$high_error) & md$high_error == 1] <- NA
  cts <- SeuratObject::LayerData(s, assay = "RNA", layer = "counts")
  if (is.null(GL)) {                                         # fix the gene axes once
    GL <- intersect(LIG, rownames(cts)); GR <- intersect(REC, rownames(cts))
    pk <- unique(pairs[, .(ligand, gene)])[ligand %in% GL & gene %in% GR]
    PIDX <- data.table(ligand = pk$ligand, gene = pk$gene,
                       li = match(pk$ligand, GL), ri = match(pk$gene, GR))
    message(sprintf("[2b] %-10s gene axes: %d ligands x %d receptor genes | %d real pairs evaluable",
                    ds, length(GL), length(GR), nrow(PIDX)))
  }
  stopifnot(all(GL %in% rownames(cts)), all(GR %in% rownames(cts)))   # same gene universe within a dataset
  for (bb in setdiff(unique(b[!is.na(b)]), NA)) {
    j <- which(b == bb); n <- length(j); if (n < MIN_CELL) next
    Dl <- as(cts[GL, j, drop = FALSE] > 0, "dMatrix")        # ligands x cells (FULL fixed axis)
    Dr <- as(cts[GR, j, drop = FALSE] > 0, "dMatrix")        # receptors x cells
    nl <- Matrix::rowSums(Dl); fpos <- nl / n
    meas <- fpos >= LIG_LO & fpos <= LIG_HI
    LIGCEN[[length(LIGCEN) + 1]] <- data.table(dataset = ds, sample = sm, bin = bb, n_cells = n,
                                               n_ligands_present = length(GL),
                                               n_ligands_measurable = sum(meas),
                                               n_receptors_present = length(GR))
    if (!any(meas)) next
    SIZE <- n - nl                                           # receiver-set size per ligand
    POS  <- matrix(Matrix::rowSums(Dr), length(GR), length(GL)) - as.matrix(Matrix::tcrossprod(Dr, Dl))
    NEG  <- matrix(SIZE, length(GR), length(GL), byrow = TRUE) - POS
    k <- cbind(PIDX$ri, PIDX$li)
    p <- POS[k]; q <- NEG[k]; mm <- meas[PIDX$li]
    keep <- mm & p >= min(ARMS) & q >= min(ARMS)             # >=30 implies >=50 is a subset
    if (!any(keep)) next
    HIT[[length(HIT) + 1]] <- data.table(dataset = ds, sample = sm, bin = bb,
                                         ligand = PIDX$ligand[keep], gene = PIDX$gene[keep],
                                         pos = p[keep], neg = q[keep])
  }
  rm(s, cts); gc(FALSE)
  if (i %% 10 == 0) message(sprintf("    ... %d/%d samples", i, nrow(elig)))
}

## -- melt, keep only real ligand-receptor pairs, count ------------------------------------------
dbp <- unique(pairs[, .(ligand, gene)])
H <- rbindlist(HIT)
long <- rbindlist(lapply(names(ARMS), function(a)
  H[pos >= ARMS[[a]] & neg >= ARMS[[a]],
    .(arm = a, n_patients = uniqueN(sample)), by = .(dataset, bin, ligand, gene)]))
fwrite(long[order(dataset, bin, arm, -n_patients)], file.path(OUT_DIR, "ws_w1_samebin_counts.csv"))
cen <- rbindlist(LIGCEN); fwrite(cen, file.path(OUT_DIR, "ws_w1_samebin_ligands.csv"))

sa <- long[n_patients >= MIN_PAT, .(combos = uniqueN(paste(ligand, gene, bin))), by = .(dataset, arm)]
sa[, scope := "all_bins"]
sm2 <- long[n_patients >= MIN_PAT & bin %in% MAIN_BINS,
            .(combos = uniqueN(paste(ligand, gene, bin))), by = .(dataset, arm)]
sm2[, scope := "main_bins"]
summ <- rbind(sa, sm2)
fwrite(summ, file.path(OUT_DIR, "ws_w1_samebin_summary.csv"))

message("\n[3] measurable-ligand census (patient x bin units with >= ", MIN_CELL, " cells)")
print(cen[, .(units = .N, median_cells = median(n_cells),
              median_measurable_ligands = median(n_ligands_measurable)), by = .(dataset, bin)][order(dataset, bin)])
message("\n[4] GATE_W1.2b: same-bin combinations with >= ", MIN_PAT, " patients")
print(dcast(summ, dataset + scope ~ arm, value.var = "combos", fill = 0L)[order(scope, dataset)])
for (sc in c("main_bins", "all_bins")) for (a in names(ARMS)) {
  k <- summ[scope == sc & arm == a & combos >= 30L]
  message(sprintf("    %-10s arm %-8s (>=%d cells): %d studies with >=30 combos -> %s",
                  sc, a, ARMS[[a]], nrow(k), if (nrow(k) >= 2L) "PASS" else "FAIL"))
}
message("\n[5] top same-bin combinations in the main bins, primary arm")
print(long[arm == "primary" & bin %in% MAIN_BINS][order(-n_patients)][1:15,
      .(dataset, bin, ligand, receptor_gene = gene, n_patients)])

## SELF-CHECKS ----------------------------------------------------------------------------------
message("\n[SELF-CHECK]")
s50 <- summ[arm == "primary" & scope == "all_bins"]; s30 <- summ[arm == "fallback" & scope == "all_bins"]
m <- merge(s50, s30, by = c("dataset", "scope"), suffixes = c("_50", "_30"))
message(sprintf("  (a) monotone: combos(30) >= combos(50) in every study -> %s",
                if (all(m$combos_30 >= m$combos_50)) "PASS" else "** FAIL **"))
cmp <- merge(summ[scope == "main_bins", .(dataset, arm, main = combos)],
             summ[scope == "all_bins",  .(dataset, arm, all = combos)], by = c("dataset", "arm"))
message(sprintf("  (b) main_bins <= all_bins in every (study, arm): %d rows compared -> %s",
                nrow(cmp), if (nrow(cmp) > 0L && all(cmp$main <= cmp$all)) "PASS" else "** FAIL **"))
message(sprintf("  (c) every reported pair is in CellChatDB -> %s",
                if (nrow(long[!paste(ligand, gene) %in% dbp[, paste(ligand, gene)]]) == 0L) "PASS" else "** FAIL **"))
message(sprintf("  (d) measurable ligands all inside [%.2f, %.2f] by construction; census max = %d of %d present",
                LIG_LO, LIG_HI, max(cen$n_ligands_measurable), max(cen$n_ligands_present)))
message(sprintf("  (e) every retained unit satisfies pos+neg == receiver-set size -> %s",
                if (nrow(H[pos + neg < MIN_CELL - 1L]) == 0L) "PASS" else "** FAIL **"))
message(sprintf("  (f) patients contributing a >=%d-cell bin: %d of %d eligible (%d had no object: %s)",
                MIN_CELL, uniqueN(cen$sample), nrow(elig), length(MISSING),
                if (length(MISSING)) paste(MISSING, collapse = ", ") else "none"))
message("[done] ", OUT_DIR)
