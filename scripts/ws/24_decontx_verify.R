#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../aml_niche_net/ws_decontx/<dataset>/_summary_step2.rds   (23's per-cell HB / granule / total UMIs)
#          results/tables/ws/decontx/<dataset>/ws_s1_contamination_per_cell.csv.gz (22's per-cell DecontX estimate)
#          /LARGE1/.../02_seurat_objects/04_annotated/<dataset>/<sample>.rds      (raw counts of ONE channel, for the spike-in)
# OUTPUT : results/tables/ws/decontx/<dataset>/ws_s1_verify_marker_estimate.csv  per sample: DecontX vs a model-free estimate
#          results/tables/ws/decontx/<dataset>/ws_s1_verify_spikein.csv          planted vs recovered contamination, one channel
#          results/tables/ws/decontx/<dataset>/ws_s1_verify.pdf, ws_s1_verify_checks.csv
# what it does: answers "how do we know DecontX is right?" two ways, neither of which uses DecontX's own model.
#   A. Model-free estimate. T, NK and B cells do not transcribe MPO, ELANE, AZU1 or PRTN3, so in those cells every
#      granule UMI is ambient. If a cell's UMIs are a fraction r ambient and the ambient pool has granule share s
#      (approximated by the channel's average profile), the cell's granule share is r * s. Per sample:
#      r_marker = (granule share in its T/NK + B cells) / (channel granule share), averaged over channels by cell
#      number. This is compared with DecontX's sample median (all cells, and lymphoid cells only) across samples.
#      The same estimate is NOT made with haemoglobin: red-cell lysate is not a captured cluster, so the channel
#      average under-represents it and the ratio would overshoot.
#   B. Planted known answer. One channel's raw counts; a random half of its cells receive extra UMIs drawn from the
#      channel-average profile so that the added UMIs are a known fraction f of each spiked cell's final total
#      (f = 0.05, 0.10, 0.20, 0.30). DecontX (auto grouping, as in 22) runs on the original and the spiked matrix;
#      recovered = estimate(spiked) - estimate(original) for spiked cells, which should track f; unspiked cells
#      should move little. This tests calibration for ambient RNA whose source IS captured; it cannot test the
#      red-cell case.
#   usage: --dataset=GSE185381 [--unit=GSM5613744] [--marker-only | --spike-only]
# This script deletes nothing and does not modify 22's or 23's outputs.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
.a <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
DATASET <- .get("dataset", "GSE185381"); UNIT <- .get("unit", "GSM5613744")
DO_A <- !("--spike-only" %in% .a); DO_B <- !("--marker-only" %in% .a)

suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
set.seed(SEED)
AN   <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated", DATASET)
OBJD <- file.path(LARGE1_DIR, "ws_decontx", DATASET)
TABD <- file.path(TAB_DIR, "ws", "decontx", DATASET)
F_LEVELS <- c(0.05, 0.10, 0.20, 0.30)
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-60s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}
pdf(file.path(TABD, "ws_s1_verify.pdf"), width = 10, height = 5); par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))

## -- A. model-free estimate from granule genes in lymphoid cells ----------------------------------------
if (DO_A) {
  ts("[A] model-free estimate, ", DATASET)
  S2 <- readRDS(file.path(OBJD, "_summary_step2.rds")); CM <- rbindlist(S2$CELLM)
  PC <- fread(file.path(TABD, "ws_s1_contamination_per_cell.csv.gz"), select = c("cell", "sample", "unit", "group", "contamination"))
  D  <- merge(CM[, .(cell, gran_raw, tot_raw)], PC, by = "cell")
  chk("per-cell tables align", nrow(PC), nrow(D))
  soup <- D[, .(soup_gran_share = sum(gran_raw) / sum(tot_raw)), by = unit]           # channel-average profile
  L <- merge(D[group %in% c("T_NK", "B_Plasma")], soup, by = "unit")
  su <- L[, .(cells = .N, gran_umi = sum(gran_raw), r_marker = (sum(gran_raw) / sum(tot_raw)) / soup_gran_share[1],
              decontx_lymph = median(contamination)), by = .(sample, unit)]
  SS <- su[, .(cells_lymph = sum(cells), gran_umi = sum(gran_umi), r_marker = pmin(weighted.mean(r_marker, cells), 1),
               decontx_lymph_median = weighted.mean(decontx_lymph, cells)), by = sample]
  SS <- merge(SS, D[, .(decontx_all_median = median(contamination)), by = sample], by = "sample")
  ok <- SS[gran_umi >= 30]                                                               # below this the ratio is noise
  r_all <- cor(ok$r_marker, ok$decontx_all_median, method = "spearman"); r_ly <- cor(ok$r_marker, ok$decontx_lymph_median, method = "spearman")
  ts(sprintf("    %d of %d samples with >= 30 granule UMIs in lymphoid cells | Spearman(marker estimate, DecontX sample median) = %.2f; vs DecontX lymphoid median = %.2f",
             nrow(ok), nrow(SS), r_all, r_ly))
  ts(sprintf("    marker estimate median %.3f (IQR %.3f-%.3f) vs DecontX lymphoid median %.3f (IQR %.3f-%.3f)",
             median(ok$r_marker), quantile(ok$r_marker, .25), quantile(ok$r_marker, .75), median(ok$decontx_lymph_median), quantile(ok$decontx_lymph_median, .25), quantile(ok$decontx_lymph_median, .75)))
  chk("model-free estimate and DecontX rank samples alike (rho > 0.5)", TRUE, r_all > 0.5)
  for (k in c("r_marker", "decontx_lymph_median", "decontx_all_median")) SS[, (k) := signif(get(k), 4)]
  fwrite(SS[order(-r_marker)], file.path(TABD, "ws_s1_verify_marker_estimate.csv"))
  plot(ok$r_marker, ok$decontx_all_median, pch = 19, col = "#2a78d6", xlab = "model-free estimate (granule genes in T/NK + B)",
       ylab = "DecontX, sample median", main = sprintf("A: %d samples, rho = %.2f", nrow(ok), r_all)); abline(0, 1, lty = 2)
}

## -- B. planted known answer on one channel ---------------------------------------------------------------
if (DO_B) {
  ts("[B] spike-in on ", UNIT)
  IX <- fread(file.path(TABD, "ws_decontx_index.csv.gz"))[unit == UNIT]
  mats <- lapply(unique(IX$sample), function(sm) { o <- readRDS(file.path(AN, paste0(sm, ".rds")))
    m <- LayerData(o, assay = "RNA", layer = "counts"); m[, intersect(colnames(m), IX$cell), drop = FALSE] })
  X <- do.call(cbind, mats)[, IX$cell]; rm(mats); n <- ncol(X); tot <- Matrix::colSums(X)
  p <- Matrix::rowSums(X) / sum(X)                                                         # channel-average profile
  spiked <- sample(n, floor(n / 2)); f <- sample(F_LEVELS, length(spiked), replace = TRUE)
  add <- round(tot[spiked] * f / (1 - f))                                                  # so that added / (tot + added) = f
  E <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = dim(X))
  cols <- integer(0); rows <- integer(0); vals <- numeric(0)
  for (k in seq_along(spiked)) { v <- rmultinom(1, add[k], p)[, 1]; nz <- which(v > 0)
    rows <- c(rows, nz); cols <- c(cols, rep(spiked[k], length(nz))); vals <- c(vals, v[nz]) }
  E <- Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = dim(X)); Xs <- X + E; dimnames(Xs) <- dimnames(X)
  planted <- numeric(n); planted[spiked] <- add / (tot[spiked] + add)
  chk("planted fraction equals the target within 1%", TRUE, max(abs(planted[spiked] - f)) < 0.01)
  t0 <- Sys.time(); r0 <- celda::decontX(X, seed = SEED, verbose = FALSE); r1 <- celda::decontX(Xs, seed = SEED, verbose = FALSE)
  ts(sprintf("    two decontX runs in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  c0 <- as.numeric(r0$contamination); c1 <- as.numeric(r1$contamination); rec <- c1 - c0
  R <- data.table(f = planted[spiked], recovered = rec[spiked], base = c0[spiked])
  out <- R[, .(cells = .N, recovered_mean = mean(recovered), recovered_median = median(recovered), recovered_sd = sd(recovered),
               baseline_mean = mean(base)), by = .(planted = round(f, 2))][order(planted)]
  out <- rbind(out, data.table(planted = 0, cells = n - length(spiked), recovered_mean = mean(rec[-spiked]), recovered_median = median(rec[-spiked]),
                               recovered_sd = sd(rec[-spiked]), baseline_mean = mean(c0[-spiked])))[order(planted)]
  slope <- coef(lm(recovered ~ f, data = R))[2]
  ts(sprintf("    recovered / planted slope %.2f; by level: %s | unspiked cells move %.3f (sd %.3f)", slope,
             paste(out[planted > 0, sprintf("%.2f -> %.3f", planted, recovered_mean)], collapse = ", "), out[planted == 0, recovered_mean], out[planted == 0, recovered_sd]))
  chk("planted ambient recovered monotonically", TRUE, all(diff(out[planted > 0, recovered_mean]) > 0))
  chk("recovery slope within 0.5-1.5", TRUE, slope > 0.5 && slope < 1.5)
  chk("unspiked cells move less than 0.02 on average", TRUE, abs(out[planted == 0, recovered_mean]) < 0.02)
  out[, unit := UNIT]; for (k in 3:6) out[[k]] <- signif(out[[k]], 4)
  fwrite(out, file.path(TABD, "ws_s1_verify_spikein.csv")); print(out)
  boxplot(recovered ~ round(f, 2), data = R, col = "grey85", xlab = "planted ambient fraction", ylab = "DecontX estimate, spiked - original",
          main = sprintf("B: %s, slope %.2f", UNIT, slope)); abline(h = 0, lty = 3)
  points(seq_along(F_LEVELS), F_LEVELS, pch = 4, col = "#eb6834", lwd = 2)
}
dev.off()
C <- rbindlist(CHK); fwrite(C, file.path(TABD, "ws_s1_verify_checks.csv"))
ts(sprintf("[done] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
