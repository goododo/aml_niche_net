# 05b_h2_pseudotime_sensitivity.R
# input : results/tables/loy/04_h2/stage0_cells_<case>.csv.gz (depth-matched leukaemic set, from Stage 0)
#         results/tables/loy/04_h2/stage0_counts.csv           (eligible strata per case)
#         the case's Seurat object (Lambo pack RDS / GSE269669), carrying predicted_Pseudotime in meta.data
#         scripts/loy/gene_sets_h2_v2.tsv                       (v2 feature sets, PAR excluded from scores)
# output: results/tables/loy/04_h2/sens_pseudotime_<case>.csv  per feature: baseline vs +pseudotime (d, ci, p_perm)
#         results/tables/loy/04_h2/sens_pseudotime_summary.csv  the immune panel across all cases + the -Y pseudotime gap
# what it does: POST-HOC sensitivity, NOT the pre-registered model (PREREGISTRATION_H2.md freezes score ~ label + fine +
#   log_depth). It adds BoneMarrowMap predicted_Pseudotime as a continuous covariate to that design and re-fits the immune
#   panel, to test whether the -Y HLA-II / IFN signal is a within-fine-label differentiation-state effect (the mundane
#   explanation that -Y cells are slightly more mature). It also reports whether -Y (A) cells sit later in pseudotime than
#   Y+ (B) cells. Baseline (no pseudotime, on the pseudotime-available cells) should reproduce contrast_<case>.csv.

suppressPackageStartupMessages({ library(data.table); library(here); library(Matrix); library(Seurat); library(SeuratObject) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "loy", "config_loy.R"))
source(here::here("scripts", "loy", "00_loy_common.R"))

set.seed(LOY_H2_SEED)
N_PERM <- 999L
CASES  <- c("AML15", "AML13", "PT43")
t0 <- Sys.time()

## --- scoring (identical to 05): z-score each gene over the given cells, mean over the set's genes present (>=5) ---------
score_sets <- function(Xc, sets) {
  Xd <- as.matrix(Matrix::t(Xc)); mu <- colMeans(Xd); sdv <- apply(Xd, 2, sd); sdv[sdv == 0] <- 1
  Z <- sweep(sweep(Xd, 2, mu), 2, sdv, "/")
  sc <- lapply(sets, function(g) { g <- intersect(g, colnames(Z)); if (length(g) < 5) return(NULL); rowMeans(Z[, g, drop = FALSE]) })
  sc <- sc[!vapply(sc, is.null, TRUE)]; do.call(cbind, sc)
}

## --- fit: same Frisch-Waugh residualisation as 05's fit_arm, but the design W optionally carries pseudotime ------------
# Residualising scores S and label g on W, then regressing, gives the coefficient of g in `score ~ g + W` (Frisch-Waugh).
# Adding pseudotime to W removes its linear contribution before the -Y coefficient is read.
fit_arm <- function(d, S, n_perm, use_pt) {
  g <- as.numeric(d$group == "A")
  base <- model.matrix(if (uniqueN(d$fine) > 1) ~ fine + log_depth else ~ log_depth, data = d[, .(fine = factor(fine), log_depth)])
  W    <- if (use_pt) cbind(base, pseudotime = d$pseudotime) else base
  Q <- qr.Q(qr(W)); resid <- function(v) v - Q %*% crossprod(Q, v)
  Sr <- resid(S); gr <- as.numeric(resid(g)); gg <- sum(gr^2)
  coef <- as.numeric(crossprod(gr, Sr)) / gg
  res  <- Sr - outer(gr, coef); dfr <- nrow(S) - ncol(W) - 1
  sigma <- sqrt(colSums(res^2) / dfr); se <- sigma / sqrt(gg); d_hat <- coef / sigma
  strata <- d[, paste(fine, cut(log_depth, unique(quantile(log_depth, 0:LOY_H2_N_DEPTH_BINS / LOY_H2_N_DEPTH_BINS)), include.lowest = TRUE))]
  idx <- split(seq_along(g), strata); perm_d <- matrix(NA_real_, n_perm, ncol(S))
  for (p in seq_len(n_perm)) {
    gp <- g; for (ix in idx) if (length(ix) > 1) gp[ix] <- g[ix][sample.int(length(ix))]
    gpr <- as.numeric(resid(gp)); cp <- as.numeric(crossprod(gpr, Sr)) / sum(gpr^2)
    perm_d[p, ] <- cp / sqrt(colSums((Sr - outer(gpr, cp))^2) / dfr)
  }
  p_perm <- (1 + colSums(abs(perm_d) >= matrix(abs(d_hat), n_perm, ncol(S), byrow = TRUE))) / (n_perm + 1)
  data.table(feature = colnames(S), n = nrow(S), d = d_hat, ci_lo = (coef - 1.96 * se) / sigma, ci_hi = (coef + 1.96 * se) / sigma, p_perm = p_perm)
}

## --- gene sets ---------------------------------------------------------------------------------------------------------
gs    <- fread(LOY_H2_GENESETS)
chrY  <- gs[set == "chrY_exclude", gene]; par <- gs[set == "PAR_exclude", gene]
PANEL <- c("MP5", "F5_IFNA", "F6_IFNG", "MP6", "MP9", "F12_chrY")            # immune panel + artefact control + machinery
sets  <- split(gs[set %in% PANEL, gene], gs[set %in% PANEL, set])
pt_of <- function(md) {                                                      # pull a pseudotime vector (named by cell) from a meta.data
  if (is.null(md)) return(NULL)
  col <- grep("pseudotime", names(md), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(col)) return(NULL)
  setNames(as.numeric(md[[col]]), rownames(md))
}

## --- one case ----------------------------------------------------------------------------------------------------------
run_case <- function(cs, s_full) {
  lib   <- LOY_H2_CASES[case == cs & role == "contrast"]; stopifnot(nrow(lib) == 1)
  cl    <- fread(file.path(LOY_H2_DIR, sprintf("stage0_cells_%s.csv.gz", cs)))
  elig  <- fread(file.path(LOY_H2_DIR, "stage0_counts.csv"))[case == cs & eligible == TRUE]
  leukS <- cl[arm == "leukaemic" & fine != "ALL" & matched == TRUE & fine %in% elig[arm == "leukaemic", fine]]
  all_leuk_cells <- cl[arm == "leukaemic" & fine == "ALL", cell]
  need  <- unique(c(leukS$cell, all_leuk_cells))

  if (lib$dataset == "GSE269669") {
    ptv <- numeric(0); mats <- list()
    for (obj in c("avm_aml", "avm_tme")) {                                    # an object with none of the needed leukaemic cells is skipped
      s <- read_gse269669(obj); v <- pt_of(s@meta.data); if (!is.null(v)) ptv <- c(ptv, v)
      hit <- intersect(colnames(s), need); if (length(hit)) mats[[obj]] <- counts_of(subset(s, cells = hit)); rm(s); gc()
    }
    gmn <- Reduce(intersect, lapply(mats, rownames)); cnt <- do.call(cbind, lapply(mats, function(m) m[gmn, , drop = FALSE])); rm(mats); gc()
    if (!length(ptv)) ptv <- NULL
  } else {
    s <- subset(s_full, cells = intersect(colnames(s_full), need)); cnt <- counts_of(s); ptv <- pt_of(s@meta.data); rm(s); gc()
  }

  det  <- Matrix::rowMeans(cnt[, all_leuk_cells, drop = FALSE] > 0)
  keep <- setdiff(names(det)[det >= LOY_H2_MIN_DETECT & !names(det) %in% chrY], par)
  chrY_present <- intersect(rownames(cnt), gs[set == "F12_chrY", gene])
  lognorm <- function(m) { sf <- 1e4 / Matrix::colSums(m); m@x <- log1p(m@x * rep.int(sf, diff(m@p))); m }
  X <- lognorm(cnt[, need, drop = FALSE]); rm(cnt); gc()

  d <- copy(leukS)[order(cell)]
  d[, pseudotime := if (is.null(ptv)) NA_real_ else ptv[cell]]
  pt_exists <- !is.null(ptv) && any(!is.na(d$pseudotime)) && sd(d$pseudotime, na.rm = TRUE) > 0
  d_common  <- if (pt_exists) d[!is.na(pseudotime)] else d
  pt_ok     <- pt_exists && nrow(d_common) >= LOY_H2_MIN_CELLS && sd(d_common$pseudotime) > 0
  message(sprintf("[%s] matched leukaemic cells: %d (A %d, B %d); pseudotime present: %d; usable: %s",
                  cs, nrow(d), sum(d$group=="A"), sum(d$group=="B"), sum(!is.na(d$pseudotime)), pt_ok))

  score_on <- function(dd) {
    S <- score_sets(X[keep, dd$cell, drop = FALSE], sets[setdiff(names(sets), "F12_chrY")])
    if (length(chrY_present) >= 1) { Y <- as.matrix(Matrix::t(X[chrY_present, dd$cell, drop = FALSE])); S <- cbind(S, F12_chrY = rowMeans(Y)) }
    S
  }
  Sc <- score_on(d_common)
  base_fit <- fit_arm(d_common, Sc, N_PERM, use_pt = FALSE)
  res <- base_fit[, .(feature, n, d_base = d, p_base = p_perm, ci_lo_base = ci_lo, ci_hi_base = ci_hi)]
  pt_gap <- NA_real_; pt_raw <- NA_real_
  if (pt_ok) {
    pt_fit <- fit_arm(d_common, Sc, N_PERM, use_pt = TRUE)
    res <- merge(res, pt_fit[, .(feature, d_pt = d, p_pt = p_perm, ci_lo_pt = ci_lo, ci_hi_pt = ci_hi)], by = "feature")
    res[, delta_d := d_pt - d_base]
    # is -Y later in pseudotime? coefficient of group on pseudotime, residualised on fine + log_depth (pseudotime as outcome)
    Wd <- model.matrix(if (uniqueN(d_common$fine) > 1) ~ fine + log_depth else ~ log_depth, data = d_common[, .(fine = factor(fine), log_depth)])
    Qd <- qr.Q(qr(Wd)); rr <- function(v) v - Qd %*% crossprod(Qd, v)
    gres <- as.numeric(rr(as.numeric(d_common$group == "A"))); ptres <- as.numeric(rr(d_common$pseudotime))
    pt_gap <- as.numeric(crossprod(gres, ptres)) / sum(gres^2)
    pt_raw <- d_common[group == "A", mean(pseudotime)] - d_common[group == "B", mean(pseudotime)]
  } else {
    res[, c("d_pt","p_pt","ci_lo_pt","ci_hi_pt","delta_d") := NA_real_]
  }
  res[, `:=`(case = cs, n_cells = nrow(d_common), pt_gap_adj = pt_gap, pt_gap_raw = pt_raw, pt_usable = pt_ok)]
  fwrite(res, file.path(LOY_H2_DIR, sprintf("sens_pseudotime_%s.csv", cs)))
  message(sprintf("[%s] -Y minus Y+ pseudotime: raw %.3f, fine+depth-adjusted %.3f", cs, pt_raw, pt_gap))
  res
}

## --- run, loading the shared Lambo object once -------------------------------------------------------------------------
all <- list()
lambo_cases <- CASES[LOY_H2_CASES[match(CASES, case), dataset] == "ZENG_Lambo_2023"]
if (length(lambo_cases)) {
  stem <- LOY_EXT_STUDIES[dataset == "ZENG_Lambo_2023", file_stem]
  f <- file.path(LOY_ZENG_RDS_DIR, paste0(stem, "_projected_cleaned.rds")); message("reading shared Lambo object ", basename(f))
  s_full <- readRDS(f); message(sprintf("  Lambo object: %d cells; pseudotime col: %s", ncol(s_full), paste(grep("seudotime", names(s_full@meta.data), value = TRUE), collapse = ",")))
  for (cs in lambo_cases) all[[cs]] <- run_case(cs, s_full)
  rm(s_full); gc()
}
for (cs in setdiff(CASES, lambo_cases)) all[[cs]] <- run_case(cs, NULL)

summ <- rbindlist(all, fill = TRUE)[feature %in% PANEL]
setcolorder(summ, c("case","feature","n_cells","d_base","d_pt","delta_d","p_base","p_pt","pt_gap_raw","pt_gap_adj","pt_usable"))
fwrite(summ[order(case, feature)], file.path(LOY_H2_DIR, "sens_pseudotime_summary.csv"))
message("\n=== sens_pseudotime_summary ===")
print(summ[order(case, feature), .(case, feature, n_cells, d_base = round(d_base,3), d_pt = round(d_pt,3), delta_d = round(delta_d,3), p_pt = round(p_pt,4), pt_gap_adj = round(pt_gap_adj,3), pt_usable)])
message(sprintf("done (%.1f min)", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
