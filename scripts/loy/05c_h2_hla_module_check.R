# 05c_h2_hla_module_check.R
# input : results/tables/loy/04_h2/stage0_cells_<case>.csv.gz, stage0_counts.csv, the case's Seurat object,
#         scripts/loy/gene_sets_h2_v2.tsv (for MP5, F5_IFNA, F6_IFNG)
# output: results/tables/loy/04_h2/hla_module_<case>.csv        per-gene d for the HLA-II structural + CIITA-regulatory modules
#         results/tables/loy/04_h2/hla_mp5_designs_<case>.csv    MP5 d under baseline / +pseudotime / +IFN / +both
#         results/tables/loy/04_h2/hla_module_summary.csv        all cases, the key rows
# what it does: GATE 01 computational check (post-hoc, not pre-registered). For the leukaemic matched set of each case:
#   (1) re-fit MP5 (HLA-II) adding the IFN score (F5_IFNA+F6_IFNG) and/or pseudotime as covariates -> is the HLA-II up
#   downstream of an interferon program? (2) per-gene d for the HLA-II structural genes and the CIITA regulatory module
#   -> is it a coordinated module (CIITA/RFX/STAT1 move with HLA-DR) or a few highly expressed genes (ambient/lineage)?
#   (3) report -Y vs Y+ fine-label composition and the MP5<->IFN gene overlap. Decisive case = AML15.

suppressPackageStartupMessages({ library(data.table); library(here); library(Matrix); library(Seurat); library(SeuratObject) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "loy", "config_loy.R")); source(here::here("scripts", "loy", "00_loy_common.R"))
set.seed(LOY_H2_SEED); N_PERM <- 999L; CASES <- c("AML15", "AML13", "PT43"); t0 <- Sys.time()

HLA2_STRUCT <- c("HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DMA","HLA-DMB","HLA-DOA","CD74")
CIITA_MOD   <- c("CIITA","RFX5","RFXAP","RFXANK","STAT1","IRF1","NLRC5")

score_sets <- function(Xc, sets) {                                           # z-score each gene over these cells, mean over the set (>=5 present)
  Xd <- as.matrix(Matrix::t(Xc)); mu <- colMeans(Xd); sdv <- apply(Xd, 2, sd); sdv[sdv == 0] <- 1
  Z <- sweep(sweep(Xd, 2, mu), 2, sdv, "/")
  sc <- lapply(sets, function(g) { g <- intersect(g, colnames(Z)); if (length(g) < 5) return(NULL); rowMeans(Z[, g, drop = FALSE]) })
  do.call(cbind, sc[!vapply(sc, is.null, TRUE)])
}
zcols <- function(Xc) { Xd <- as.matrix(Matrix::t(Xc)); mu <- colMeans(Xd); sdv <- apply(Xd, 2, sd); sdv[sdv == 0] <- 1; sweep(sweep(Xd, 2, mu), 2, sdv, "/") }
fit_arm <- function(d, S, n_perm, Zextra = NULL) {                           # Frisch-Waugh as in 05; Zextra = extra covariate columns appended to W
  g <- as.numeric(d$group == "A")
  base <- model.matrix(if (uniqueN(d$fine) > 1) ~ fine + log_depth else ~ log_depth, data = d[, .(fine = factor(fine), log_depth)])
  W <- if (!is.null(Zextra)) cbind(base, Zextra) else base
  Q <- qr.Q(qr(W)); resid <- function(v) v - Q %*% crossprod(Q, v)
  Sr <- resid(S); gr <- as.numeric(resid(g)); gg <- sum(gr^2); coef <- as.numeric(crossprod(gr, Sr)) / gg
  res <- Sr - outer(gr, coef); dfr <- nrow(S) - ncol(W) - 1; sigma <- sqrt(colSums(res^2) / dfr); d_hat <- coef / sigma
  strata <- d[, paste(fine, cut(log_depth, unique(quantile(log_depth, 0:LOY_H2_N_DEPTH_BINS / LOY_H2_N_DEPTH_BINS)), include.lowest = TRUE))]
  idx <- split(seq_along(g), strata); perm_d <- matrix(NA_real_, n_perm, ncol(S))
  for (p in seq_len(n_perm)) { gp <- g; for (ix in idx) if (length(ix) > 1) gp[ix] <- g[ix][sample.int(length(ix))]
    gpr <- as.numeric(resid(gp)); cp <- as.numeric(crossprod(gpr, Sr)) / sum(gpr^2); perm_d[p, ] <- cp / sqrt(colSums((Sr - outer(gpr, cp))^2) / dfr) }
  p_perm <- (1 + colSums(abs(perm_d) >= matrix(abs(d_hat), n_perm, ncol(S), byrow = TRUE))) / (n_perm + 1)
  data.table(feature = colnames(S), d = d_hat, p_perm = p_perm)
}

gs <- fread(LOY_H2_GENESETS); chrY <- gs[set == "chrY_exclude", gene]; par <- gs[set == "PAR_exclude", gene]
score_gsets <- split(gs[set %in% c("MP5","F5_IFNA","F6_IFNG"), gene], gs[set %in% c("MP5","F5_IFNA","F6_IFNG"), set])
mp5_overlap_ifn <- length(intersect(score_gsets$MP5, c(score_gsets$F5_IFNA, score_gsets$F6_IFNG)))
pt_of <- function(md) { if (is.null(md)) return(NULL); col <- grep("pseudotime", names(md), ignore.case = TRUE, value = TRUE)[1]; if (is.na(col)) NULL else setNames(as.numeric(md[[col]]), rownames(md)) }

run_case <- function(cs, s_full) {
  lib <- LOY_H2_CASES[case == cs & role == "contrast"]
  cl <- fread(file.path(LOY_H2_DIR, sprintf("stage0_cells_%s.csv.gz", cs)))
  elig <- fread(file.path(LOY_H2_DIR, "stage0_counts.csv"))[case == cs & eligible == TRUE]
  leukS <- cl[arm == "leukaemic" & fine != "ALL" & matched == TRUE & fine %in% elig[arm == "leukaemic", fine]]
  all_leuk_cells <- cl[arm == "leukaemic" & fine == "ALL", cell]; need <- unique(c(leukS$cell, all_leuk_cells))
  if (lib$dataset == "GSE269669") {
    ptv <- numeric(0); mats <- list()
    for (obj in c("avm_aml","avm_tme")) { s <- read_gse269669(obj); v <- pt_of(s@meta.data); if (!is.null(v)) ptv <- c(ptv, v)
      hit <- intersect(colnames(s), need); if (length(hit)) mats[[obj]] <- counts_of(subset(s, cells = hit)); rm(s); gc() }
    gmn <- Reduce(intersect, lapply(mats, rownames)); cnt <- do.call(cbind, lapply(mats, function(m) m[gmn, , drop = FALSE])); rm(mats); gc(); if (!length(ptv)) ptv <- NULL
  } else { s <- subset(s_full, cells = intersect(colnames(s_full), need)); cnt <- counts_of(s); ptv <- pt_of(s@meta.data); rm(s); gc() }
  lognorm <- function(m) { sf <- 1e4 / Matrix::colSums(m); m@x <- log1p(m@x * rep.int(sf, diff(m@p))); m }
  X <- lognorm(cnt[, need, drop = FALSE]); det <- Matrix::rowMeans(cnt[, all_leuk_cells, drop = FALSE] > 0); rm(cnt); gc()
  keep <- setdiff(names(det)[det >= LOY_H2_MIN_DETECT & !names(det) %in% chrY], par)

  d <- copy(leukS)[order(cell)]; d[, pseudotime := if (is.null(ptv)) NA_real_ else ptv[cell]]
  has_pt <- !is.null(ptv) && all(!is.na(d$pseudotime)) && sd(d$pseudotime, na.rm = TRUE) > 0
  cells <- d$cell
  Ssc <- score_sets(X[intersect(keep, rownames(X)), cells, drop = FALSE], score_gsets)   # MP5 / IFN-a / IFN-g scores on the matched set
  mp5 <- Ssc[, "MP5", drop = FALSE]
  ifn <- Ssc[, intersect(c("F5_IFNA","F6_IFNG"), colnames(Ssc)), drop = FALSE]
  designs <- list(baseline = NULL, plus_IFN = ifn)
  if (has_pt) { designs$plus_pseudotime <- cbind(pseudotime = d$pseudotime); designs$plus_both <- cbind(pseudotime = d$pseudotime, ifn) }
  mp5res <- rbindlist(lapply(names(designs), function(nm) { r <- fit_arm(d, mp5, N_PERM, designs[[nm]]); r[, design := nm]; r }), fill = TRUE)
  mp5res[, `:=`(case = cs, mp5_ifn_gene_overlap = mp5_overlap_ifn)]
  fwrite(mp5res, file.path(LOY_H2_DIR, sprintf("hla_mp5_designs_%s.csv", cs)))

  mod <- c(HLA2_STRUCT, CIITA_MOD); present <- intersect(mod, rownames(X))
  Sg <- zcols(X[present, cells, drop = FALSE]); colnames(Sg) <- present
  gfit <- fit_arm(d, Sg, N_PERM)
  gfit[, `:=`(module = fifelse(feature %in% HLA2_STRUCT, "HLA2_struct", "CIITA_reg"),
              det_frac = round(det[feature], 3), case = cs)]
  fwrite(gfit[order(module, -abs(d))], file.path(LOY_H2_DIR, sprintf("hla_module_%s.csv", cs)))

  comp <- d[, .(nA = sum(group == "A"), nB = sum(group == "B")), by = fine][, `:=`(fracA = round(nA/sum(nA),3), fracB = round(nB/sum(nB),3))]
  message(sprintf("[%s] cells %d (A %d, B %d); pseudotime usable %s; MP5<->IFN gene overlap %d; module genes present %d/%d",
                  cs, nrow(d), sum(d$group=="A"), sum(d$group=="B"), has_pt, mp5_overlap_ifn, length(present), length(mod)))
  message(sprintf("[%s] MP5 d by design: %s", cs, paste(sprintf("%s=%.3f(p%.3f)", mp5res$design, mp5res$d, mp5res$p_perm), collapse = "  ")))
  message(sprintf("[%s] CIITA-module d: %s", cs, paste(sprintf("%s=%.2f", gfit[module=="CIITA_reg", feature], gfit[module=="CIITA_reg", d]), collapse = "  ")))
  message(sprintf("[%s] fine-label composition (-Y vs Y+): %s", cs, paste(sprintf("%s A%.2f/B%.2f", comp$fine, comp$fracA, comp$fracB), collapse = "  ")))
  list(mp5 = mp5res, genes = gfit)
}

all <- list()
lambo <- CASES[LOY_H2_CASES[match(CASES, case), dataset] == "ZENG_Lambo_2023"]
if (length(lambo)) { stem <- LOY_EXT_STUDIES[dataset == "ZENG_Lambo_2023", file_stem]
  s_full <- readRDS(file.path(LOY_ZENG_RDS_DIR, paste0(stem, "_projected_cleaned.rds"))); message("Lambo loaded: ", ncol(s_full), " cells")
  for (cs in lambo) all[[cs]] <- run_case(cs, s_full); rm(s_full); gc() }
for (cs in setdiff(CASES, lambo)) all[[cs]] <- run_case(cs, NULL)

mp5all <- rbindlist(lapply(all, `[[`, "mp5"), fill = TRUE)
gall   <- rbindlist(lapply(all, `[[`, "genes"), fill = TRUE)
fwrite(dcast(mp5all, case + mp5_ifn_gene_overlap ~ design, value.var = "d"), file.path(LOY_H2_DIR, "hla_mp5_designs_summary.csv"))
fwrite(gall[order(case, module, -abs(d))], file.path(LOY_H2_DIR, "hla_module_summary.csv"))
message("\n=== MP5 by design (d) ==="); print(dcast(mp5all, case ~ design, value.var = "d"))
message("\n=== CIITA regulatory module (d, p) ==="); print(gall[module == "CIITA_reg", .(case, feature, d = round(d,3), p = round(p_perm,3), det = det_frac)][order(case, -abs(d))])
message("\n=== HLA-II structural genes (d) ==="); print(gall[module == "HLA2_struct", .(case, feature, d = round(d,3), p = round(p_perm,3), det = det_frac)][order(case, -abs(d))])
message(sprintf("done (%.1f min)", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
