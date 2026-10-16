#!/usr/bin/env Rscript
# INPUT  : results/tables/ws/decontx/<dataset>/{ws_decontx_index.csv.gz, ws_decontx_units.csv, <unit>__percell.csv.gz}
#          /LARGE1/.../aml_niche_net/ws_decontx/<dataset>/<unit>.rds            (22's decontXcounts, per unit)
#          /LARGE1/.../02_seurat_objects/04_annotated/<dataset>/<sample>.rds     (raw counts, read once per unit)
#          results/tables/ws/ws_doublet_scores_percell.csv.gz                   (scDblFinder score, union flag)
#          /LARGE1/.../00_raw/public/GEO/GSE185381_RAW/<GSM>_*_RNA_soupx_processed.csv.gz  (depositors' SoupX matrices;
#          GSE185381 only, used as an external reference, never ingested into the pipeline)
# OUTPUT : results/tables/ws/decontx/<dataset>/ws_s1_contamination_per_cell.csv.gz   cell, sample, unit, group, contamination
#          results/tables/ws/decontx/<dataset>/ws_s1_contamination_per_sample.csv    median/IQR, HB and granule shares, doublets
#          results/tables/ws/decontx/<dataset>/ws_s1_ambient_genes.csv              per unit: top 200 genes by removed fraction
#          results/tables/ws/decontx/<dataset>/ws_s1_ambient_blacklist.csv          genes in a unit's top 100 (>= 50 UMIs), with frequency
#          results/tables/ws/decontx/<dataset>/ws_s1_soupx_concordance.csv          per unit: DecontX vs SoupX removed fractions
#          results/tables/ws/decontx/<dataset>/ws_s1_validation.pdf                 the two S1 check metrics, doublets, SoupX
#          results/tables/ws/decontx/<dataset>/ws_s1_checks.csv
# what it does: v2.0 S1, step 2 of 2. Collects 22's per-channel DecontX output into the four S1 deliverables and
#   checks them three ways: (a) the two metrics v2.0 §5.3 S1 names (haemoglobin UMI share in non-erythroid cells;
#   MPO/ELANE/AZU1/PRTN3 share in T/NK and B cells) must rise with the DecontX contamination; (b) cells decontX
#   calls heavily contaminated should be enriched for doublets (a mixed profile is what both look like);
#   (c) the per-gene fraction of UMIs removed should agree with the depositors' own SoupX run (GSE185381).
#   Ambient genes are ranked by the FRACTION of a gene's UMIs removed, not by removed UMIs: the latter is a list of
#   the most abundant genes in the channel (MALAT1, MT-, RP-), which is what ambient RNA is made of but not what
#   singles out contamination-driven genes for the S5 blacklist.
#   usage: --dataset=GSE185381 [--no-soupx]
# This script deletes nothing and does not modify 22's outputs or the 04_annotated objects.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
.a <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
DATASET <- .get("dataset", "GSE185381"); DO_SOUPX <- !("--no-soupx" %in% .a) && DATASET == "GSE185381"

suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
AN   <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated", DATASET)
OBJD <- file.path(LARGE1_DIR, "ws_decontx", DATASET)
TABD <- file.path(TAB_DIR, "ws", "decontx", DATASET)
RAW  <- file.path(LARGE1_DIR, "00_raw/public/GEO/GSE185381_RAW")
EARLY_LYMPHOID_FINE <- c("CLP", "MLP", "MLP-II", "Pre-ProB")
HB <- c("HBB", "HBA1", "HBA2"); GRAN <- c("MPO", "ELANE", "AZU1", "PRTN3")
MIN_UMI_GENE <- 50L      # a gene enters the removed-fraction ranking only with this many raw UMIs in the unit
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-60s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- 1. per-cell table -------------------------------------------------------------------------------
ts("[1] per-cell contamination, ", DATASET)
IX <- fread(file.path(TABD, "ws_decontx_index.csv.gz")); U <- fread(file.path(TABD, "ws_decontx_units.csv"))
pcf <- file.path(TABD, paste0(U$unit, "__percell.csv.gz"))
chk("every unit has a per-cell file", nrow(U), sum(file.exists(pcf)))
PC <- rbindlist(lapply(pcf[file.exists(pcf)], fread))
chk("per-cell rows equal the index", nrow(IX), nrow(PC))
chk("one zmode across units", 1L, uniqueN(PC$zmode))
DBL <- fread(file.path(TAB_DIR, "ws", WS_DBL_SCORES_TABLE))[dataset == DATASET, .(cell, sc_score, dbl_union = (sc_class %in% c(TRUE, "TRUE")) | (df_class %in% c(TRUE, "TRUE")))]
PC <- merge(PC, DBL, by = "cell", all.x = TRUE)
chk("doublet scores cover >= 95% of cells", TRUE, mean(!is.na(PC$sc_score)) >= 0.95)
fwrite(PC[, .(cell, sample, unit, group, decontx_z, contamination, umi_before, umi_after, sc_score, dbl_union)],
       file.path(TABD, "ws_s1_contamination_per_cell.csv.gz"), compress = "gzip")

## -- 2. per unit: raw counts back, removed fraction per gene, HB / granule shares, SoupX ---------------
ts("[2] per unit: removed fraction per gene, check metrics", if (DO_SOUPX) ", SoupX" else "")
soupx_files <- if (DO_SOUPX) list.files(RAW, pattern = "_RNA_soupx_processed\\.csv\\.gz$", full.names = TRUE) else character(0)
GENES <- list(); CELLM <- list(); SX <- list()
for (u in U$unit) {
  R <- readRDS(file.path(OBJD, paste0(u, ".rds"))); cells <- R$cells
  mats <- lapply(R$samples, function(sm) { o <- readRDS(file.path(AN, paste0(sm, ".rds")))
    m <- LayerData(o, assay = "RNA", layer = "counts"); m[, intersect(colnames(m), cells), drop = FALSE] })
  X <- do.call(cbind, mats)[, cells, drop = FALSE]; DX <- R$decontXcounts; rm(mats)
  stopifnot(identical(dim(X), dim(DX)), identical(rownames(X), rownames(DX)))
  raw_g <- Matrix::rowSums(X); dx_g <- Matrix::rowSums(DX); rem <- raw_g - dx_g
  G <- data.table(unit = u, gene = rownames(X), raw_umi = raw_g, removed_umi = rem,
                  removed_fraction = ifelse(raw_g > 0, rem / raw_g, NA_real_), removed_share = rem / sum(rem))
  GENES[[u]] <- G[raw_umi >= MIN_UMI_GENE][order(-removed_fraction)][seq_len(min(.N, 200L))][, rank := .I]
  # per-cell check metrics on RAW counts (what S5 will use), plus the decontaminated version
  grp <- R$group; lymph <- grp %in% c("T_NK", "B_Plasma")        # EarlyLymphoid is its own group, so it is excluded here
  hbg <- intersect(HB, rownames(X)); grg <- intersect(GRAN, rownames(X))
  CELLM[[u]] <- data.table(cell = cells, unit = u, nonery = grp != "Erythroid", lymph = lymph,
    hb_raw = as.numeric(Matrix::colSums(X[hbg, , drop = FALSE])), hb_dx = as.numeric(Matrix::colSums(DX[hbg, , drop = FALSE])),
    gran_raw = as.numeric(Matrix::colSums(X[grg, , drop = FALSE])), gran_dx = as.numeric(Matrix::colSums(DX[grg, , drop = FALSE])),
    tot_raw = as.numeric(Matrix::colSums(X)), tot_dx = as.numeric(Matrix::colSums(DX)))
  if (DO_SOUPX) {
    f <- soupx_files[startsWith(basename(soupx_files), paste0(u, "_"))]
    if (length(f) == 1L) {
      S <- fread(f, header = TRUE); setnames(S, 1L, "gene")
      bc <- substring(colnames(S)[-1], nchar(colnames(S)[-1]) - 15L)          # the 16-bp barcode core
      ours <- sub("-[0-9]+$", "", sub("^.*_", "", cells)); j <- match(ours, bc)
      ok <- !is.na(j); gi <- match(rownames(X), S$gene)
      if (mean(ok) >= 0.95 && mean(!is.na(gi)) >= 0.95) {
        SM <- as.matrix(S[gi[!is.na(gi)], (j[ok] + 1L), with = FALSE]); sx_g <- rowSums(SM); rm(SM)
        raw_sub <- Matrix::rowSums(X[!is.na(gi), ok, drop = FALSE]); dx_sub <- Matrix::rowSums(DX[!is.na(gi), ok, drop = FALSE])
        keep <- raw_sub >= MIN_UMI_GENE
        rf_dx <- 1 - dx_sub[keep] / raw_sub[keep]; rf_sx <- pmax(1 - sx_g[keep] / raw_sub[keep], 0)
        top <- function(v) names(sort(v, decreasing = TRUE))[seq_len(min(50L, length(v)))]
        SX[[u]] <- data.table(unit = u, cells_matched = sum(ok), cells = length(cells), genes_compared = sum(keep),
          spearman_removed_fraction = cor(rf_dx, rf_sx, method = "spearman"),
          top50_overlap = length(intersect(top(rf_dx), top(rf_sx))),
          total_removed_decontx = 1 - sum(dx_sub) / sum(raw_sub), total_removed_soupx = 1 - sum(sx_g) / sum(raw_sub))
        if (u == U$unit[1]) EX <- data.table(gene = names(rf_dx), decontx = rf_dx, soupx = rf_sx, raw = raw_sub[keep])
      } else SX[[u]] <- data.table(unit = u, cells_matched = sum(ok), cells = length(cells), genes_compared = NA_integer_)
      rm(S)
    }
  }
  ts(sprintf("    %-11s %5d cells | removed overall %.3f | top by fraction: %s", u, length(cells), 1 - sum(dx_g) / sum(raw_g),
             paste(head(GENES[[u]]$gene, 6), collapse = " ")))
  rm(X, DX, R); gc(FALSE)
}
GENES <- rbindlist(GENES); fwrite(GENES, file.path(TABD, "ws_s1_ambient_genes.csv"))
BL <- GENES[rank <= 100L, .(n_units = .N, median_removed_fraction = round(median(removed_fraction), 3),
                            median_rank = median(rank)), by = gene][order(-n_units, -median_removed_fraction)]
BL[, frac_units := round(n_units / nrow(U), 3)]
fwrite(BL, file.path(TABD, "ws_s1_ambient_blacklist.csv"))
ts(sprintf("    blacklist: %d genes in some unit's top 100; %d in at least half the units: %s", nrow(BL), BL[frac_units >= 0.5, .N],
           paste(head(BL[frac_units >= 0.5, gene], 15), collapse = " ")))
CELLM <- rbindlist(CELLM); PC <- merge(PC, CELLM[, -"unit"], by = "cell")

## -- 3. per-sample table and checks ------------------------------------------------------------------
ts("[3] per sample")
SS <- PC[, .(cells = .N, units = uniqueN(unit), median_contamination = median(contamination),
             q25 = quantile(contamination, .25), q75 = quantile(contamination, .75), frac_gt_0.5 = mean(contamination > 0.5),
             hb_share_nonery_raw = sum(hb_raw[nonery]) / sum(tot_raw[nonery]), hb_share_nonery_dx = sum(hb_dx[nonery]) / sum(tot_dx[nonery]),
             gran_share_lymph_raw = if (any(lymph)) sum(gran_raw[lymph]) / sum(tot_raw[lymph]) else NA_real_,
             gran_share_lymph_dx = if (any(lymph)) sum(gran_dx[lymph]) / sum(tot_dx[lymph]) else NA_real_,
             n_lymph = sum(lymph), doublet_score_mean = mean(sc_score, na.rm = TRUE),
             spearman_contam_doublet = suppressWarnings(cor(contamination, sc_score, method = "spearman", use = "complete.obs"))),
         by = sample][order(-median_contamination)]
for (k in grep("share|contamination|q25|q75|frac|doublet", names(SS), value = TRUE)) SS[, (k) := signif(get(k), 4)]
fwrite(SS, file.path(TABD, "ws_s1_contamination_per_sample.csv"))
chk("every sample in the index has a row", uniqueN(IX$sample), nrow(SS))
r_hb <- cor(SS$median_contamination, SS$hb_share_nonery_raw, method = "spearman")
r_gr <- cor(SS$median_contamination, SS$gran_share_lymph_raw, method = "spearman", use = "complete.obs")
ts(sprintf("    Spearman(median contamination, HB share in non-erythroid cells) = %.2f; (.., granule share in T/NK+B) = %.2f (n = %d)", r_hb, r_gr, nrow(SS)))
chk("HB check metric rises with contamination (rho > 0)", TRUE, r_hb > 0)
chk("granule check metric rises with contamination (rho > 0)", TRUE, r_gr > 0)
d_med <- PC[!is.na(dbl_union), median(contamination), by = dbl_union]
ts(sprintf("    median contamination: doublet-flagged %.3f vs other %.3f; cells > 0.5 that are flagged doublets: %.1f%% (flag rate overall %.1f%%)",
           d_med[dbl_union == TRUE, V1], d_med[dbl_union == FALSE, V1], 100 * PC[contamination > 0.5, mean(dbl_union, na.rm = TRUE)], 100 * PC[, mean(dbl_union, na.rm = TRUE)]))
chk("doublet-flagged cells carry more contamination", TRUE, d_med[dbl_union == TRUE, V1] > d_med[dbl_union == FALSE, V1])
hb_removed <- PC[nonery == TRUE, 1 - sum(hb_dx) / sum(hb_raw)]; all_removed <- PC[nonery == TRUE, 1 - sum(tot_dx) / sum(tot_raw)]
ts(sprintf("    non-erythroid cells: fraction of HB UMIs removed %.3f vs all UMIs removed %.3f", hb_removed, all_removed))
chk("HB removed preferentially in non-erythroid cells", TRUE, hb_removed > all_removed)
if ("AML3133" %in% PC$sample) {                                   # known answer from the 2026-10-16 audit: 790 HBB UMIs in 355 of 405 T cells
  a <- PC[sample == "AML3133" & group == "T_NK", .(hb_raw = sum(hb_raw), hb_dx = sum(hb_dx), cells = .N)]
  ts(sprintf("    AML3133 T_NK: HB UMIs %d in %d cells -> %.0f after DecontX (removed %.2f)", as.integer(a$hb_raw), a$cells, a$hb_dx, 1 - a$hb_dx / a$hb_raw))
  chk("AML3133 T_NK haemoglobin mostly removed (> 0.5)", TRUE, (1 - a$hb_dx / a$hb_raw) > 0.5)
}
if (length(SX)) {
  SX <- rbindlist(SX, fill = TRUE); fwrite(SX, file.path(TABD, "ws_s1_soupx_concordance.csv"))
  ok <- SX[!is.na(genes_compared)]
  ts(sprintf("    SoupX: %d of %d units compared; Spearman of per-gene removed fraction median %.2f (min %.2f); top-50 overlap median %d/50; total removed DecontX %.3f vs SoupX %.3f (medians)",
             nrow(ok), nrow(SX), median(ok$spearman_removed_fraction), min(ok$spearman_removed_fraction), as.integer(median(ok$top50_overlap)),
             median(ok$total_removed_decontx), median(ok$total_removed_soupx)))
  chk("SoupX comparison done for >= 90% of units", TRUE, nrow(ok) >= 0.9 * nrow(SX))
  chk("DecontX and SoupX agree on which genes are ambient (median rho > 0.3)", TRUE, median(ok$spearman_removed_fraction) > 0.3)
}

## -- 4. figure -----------------------------------------------------------------------------------------
ts("[4] figure")
pdf(file.path(TABD, "ws_s1_validation.pdf"), width = 10, height = 7.5); par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 1))
plot(SS$median_contamination, SS$hb_share_nonery_raw, log = "y", pch = 19, col = "#2a78d6", xlab = "DecontX contamination, sample median",
     ylab = "HB UMI share, non-erythroid cells (raw)", main = sprintf("S1 check (a): rho = %.2f", r_hb))
plot(SS$median_contamination, SS$gran_share_lymph_raw, log = "y", pch = 19, col = "#eb6834", xlab = "DecontX contamination, sample median",
     ylab = "MPO/ELANE/AZU1/PRTN3 share, T/NK + B (raw)", main = sprintf("S1 check (b): rho = %.2f", r_gr))
boxplot(contamination ~ dbl_union, data = PC[!is.na(dbl_union)], outline = FALSE, col = c("grey85", "#1baf7a"),
        names = c("singlet", "doublet (union)"), ylab = "DecontX contamination", main = "per cell vs doublet call")
hist(SS$median_contamination, breaks = 20, col = "grey80", main = "sample medians", xlab = "DecontX contamination")
if (exists("EX")) { plot(EX$soupx, EX$decontx, pch = ".", col = "#2a78d6", xlab = "SoupX removed fraction", ylab = "DecontX removed fraction",
     main = sprintf("%s: per gene (%d genes)", U$unit[1], nrow(EX))); abline(0, 1, lty = 2) }
if (length(SX) && nrow(ok)) hist(ok$spearman_removed_fraction, breaks = 15, col = "grey80", main = "DecontX vs SoupX, per unit", xlab = "Spearman, per-gene removed fraction")
dev.off()
C <- rbindlist(CHK); fwrite(C, file.path(TABD, "ws_s1_checks.csv"))
ts(sprintf("[5] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", TABD)
