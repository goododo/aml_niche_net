#!/usr/bin/env Rscript
# 05_h2_contrast.R ----
# input : results/tables/loy/04_h2/stage0_cells_<case>.csv.gz and stage0_counts.csv (from 04), scripts/loy/gene_sets_h2_v1.tsv,
#         raw counts of the case's library (Zeng pack Lambo object / GSE269669 author objects), 02's loy_compartment_calls.csv (c_amb)
# output: results/tables/loy/04_h2/
#           contrast_<case><sfx>.csv     per feature: n_A, n_B, d, d_corr, ci, permutation p, Holm p, sign consistency, verdict
#           control_<case><sfx>.csv      the same pipeline on the library's T cells (GATE_H2_K1)
#           perstratum_<case><sfx>.csv   d per fine label, leukaemic arm
#           checks_<case><sfx>.csv       GATE_H2_K1, K2, N1, P (P only in the planted world)
#           de_<case><sfx>.csv.gz        S1: per-gene Wilcoxon on the matched set, control-calibrated empirical FDR (real world only)
#         <sfx> = "" for --world real, "_smoke_null" / "_smoke_planted" otherwise (smoke outputs never sit next to real ones)
# what it does: the pre-registered contrast of PREREGISTRATION_H2.md sections 4-8 for ONE case. Within the depth-matched set
#   written by 04: log-normalise, keep genes detected in >= LOY_H2_MIN_DETECT of the library's leukaemic cells, drop chrY genes,
#   z-score per gene, average per feature set; per feature the OLS coefficient of the -Y label in score ~ label + fine + log depth
#   (by Frisch-Waugh residualisation, so permutations are cheap), d = coefficient / residual SD, corrected by 1/(1 - c_amb);
#   permutation p from LOY_H2_N_PERM label shuffles within fine label x depth bin; Holm over the primary panel; the T-cell control
#   arm runs the identical pipeline; wording per section 8.
# Usage : sbatch scripts/loy/submit_loy_h2.sbatch 05 --case=AML15 --world=real|null|planted
#   --world=null     labels shuffled within strata once, then the full pipeline (must traverse the no-result branches)
#   --world=planted  real labels; LOY_H2_N_PLANTS random 60-gene sets each shifted by d = LOY_H2_PLANT_D in -Y cells -> power (GATE_P)
# WHY: the tag prereg-loy-h2-v1 precedes the first run of this script on any world; the smoke worlds exist because the G07 and
#   dbl lines each hid a bug in the branch the planted signal never visits.
suppressPackageStartupMessages({ library(data.table); library(here); library(Matrix); library(Seurat); library(SeuratObject) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "loy", "config_loy.R"))
source(here::here("scripts", "loy", "00_loy_common.R"))
args  <- commandArgs(trailingOnly = TRUE)
argv  <- function(k, default) { v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE)); if (length(v)) v[1] else default }
case  <- argv("case", NA_character_); world <- argv("world", "real")
stopifnot(case %in% LOY_H2_CASES$case, world %in% c("real", "null", "planted"))
sfx   <- if (world == "real") "" else paste0("_smoke_", world)
set.seed(LOY_H2_SEED + nchar(world))
t0 <- Sys.time(); out <- function(stem) file.path(LOY_H2_DIR, paste0(stem, "_", case, sfx, ".csv"))
message(sprintf("[0] case %s | world %s | %s", case, world, format(Sys.time())))

## [1] matched cells from Stage 0 ------------------------------------------------------------------------
lib   <- LOY_H2_CASES[case == get("case", envir = globalenv()) & role == "contrast"]
stopifnot(nrow(lib) == 1)
cl    <- fread(file.path(LOY_H2_DIR, sprintf("stage0_cells_%s.csv.gz", case)))
elig  <- fread(file.path(LOY_H2_DIR, "stage0_counts.csv"))[case == get("case", envir = globalenv()) & eligible == TRUE]
leukS <- cl[arm == "leukaemic" & fine != "ALL" & matched == TRUE & fine %in% elig[arm == "leukaemic", fine]]   # per-stratum matched set, eligible strata only
ctrlS <- cl[arm == "control" & fine == "ALL" & matched == TRUE]                                                # pooled control, one "ALL" stratum per library
ctrlS[, fine := Sample]                                                                                       # library = the control's stratum (depth is library-specific)
stopifnot(nrow(leukS) > 0)
if (world == "null") {   # the null world: the real -Y label destroyed inside fine x depth strata, before anything else is computed
  shuffle_within <- function(d) { d <- copy(d); d[, group := group[sample.int(.N)], by = .(fine, dbin = cut(log_depth, unique(quantile(log_depth, 0:LOY_H2_N_DEPTH_BINS / LOY_H2_N_DEPTH_BINS)), include.lowest = TRUE))]; d[, dbin := NULL]; d }
  leukS <- shuffle_within(leukS); ctrlS <- shuffle_within(ctrlS)
}
message(sprintf("[1] leukaemic matched set: %d cells (A %d, B %d) in %d strata | control: %d cells (A %d, B %d)",
                nrow(leukS), sum(leukS$group == "A"), sum(leukS$group == "B"), uniqueN(leukS$fine), nrow(ctrlS), sum(ctrlS$group == "A"), sum(ctrlS$group == "B")))

## [2] counts of the library ----------------------------------------------------------------------------------
all_leuk_cells <- cl[arm == "leukaemic" & fine == "ALL", cell]        # every leukaemic cell with a call (for the detection filter)
need <- unique(c(leukS$cell, ctrlS$cell, all_leuk_cells))
if (lib$dataset == "GSE269669") {
  mats <- lapply(c("avm_aml", "avm_tme"), function(obj) { s <- read_gse269669(obj); s <- subset(s, cells = intersect(colnames(s), need)); m <- counts_of(s); rm(s); gc(); m })
  g <- Reduce(intersect, lapply(mats, rownames)); cnt <- do.call(cbind, lapply(mats, function(m) m[g, , drop = FALSE])); rm(mats)
} else {
  stem <- LOY_EXT_STUDIES[dataset == lib$dataset, file_stem]
  f <- file.path(LOY_ZENG_RDS_DIR, paste0(stem, "_projected_cleaned.rds")); message("    reading ", basename(f))
  s <- readRDS(f); s <- subset(s, cells = intersect(colnames(s), need)); cnt <- counts_of(s); rm(s); invisible(gc())
}
stopifnot(all(need %in% colnames(cnt)))
message(sprintf("[2] counts: %d genes x %d cells (%.1f min)", nrow(cnt), ncol(cnt), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
gs    <- fread(LOY_H2_GENESETS)
chrY  <- gs[role == "exclude", gene]
det   <- Matrix::rowMeans(cnt[, all_leuk_cells, drop = FALSE] > 0)                         # detection among ALL leukaemic cells, label not looked at
keep  <- names(det)[det >= LOY_H2_MIN_DETECT & !names(det) %in% chrY]
chrY_present <- intersect(rownames(cnt), gs[set == "F12_chrY", gene])
message(sprintf("    genes kept: %d of %d (detected >= %.0f%% of %d leukaemic cells, chrY removed); chrY genes present for F12: %d",
                length(keep), nrow(cnt), 100 * LOY_H2_MIN_DETECT, length(all_leuk_cells), length(chrY_present)))
lognorm <- function(m) { sf <- 1e4 / Matrix::colSums(m); m@x <- log1p(m@x * rep.int(sf, diff(m@p))); m }   # log1p(counts / depth * 1e4), per column
X <- lognorm(cnt[, need, drop = FALSE]); rm(cnt); invisible(gc())

## [3] feature scores and the OLS contrast by residualisation --------------------------------------------------
# score_sets(): z-score each gene over the given cells, mean over the set's genes present; sets with < 5 genes present are dropped
score_sets <- function(Xc, sets) {
  Xd <- as.matrix(Matrix::t(Xc)); mu <- colMeans(Xd); sdv <- apply(Xd, 2, sd); sdv[sdv == 0] <- 1
  Z <- sweep(sweep(Xd, 2, mu), 2, sdv, "/")
  sc <- lapply(sets, function(g) { g <- intersect(g, colnames(Z)); if (length(g) < 5) return(NULL); rowMeans(Z[, g, drop = FALSE]) })
  sc <- sc[!vapply(sc, is.null, TRUE)]; do.call(cbind, sc)
}
# fit_arm(): for one matched set, every feature at once. Design W = [fine dummies, log_depth, 1]; residualise the scores S and the
# label g on W (Frisch-Waugh: the coefficient of g in score ~ g + W equals the OLS of residualised score on residualised g).
design_of <- function(d) {        # fine label enters only when the set spans > 1 label (a one-level factor breaks model.matrix)
  if (uniqueN(d$fine) > 1) model.matrix(~ fine + log_depth, data = d[, .(fine = factor(fine), log_depth)]) else model.matrix(~ log_depth, data = d[, .(log_depth)])
}
fit_arm <- function(d, S, n_perm) {
  g <- as.numeric(d$group == "A")
  W <- design_of(d)
  Q <- qr.Q(qr(W)); resid <- function(v) v - Q %*% crossprod(Q, v)
  Sr <- resid(S); gr <- as.numeric(resid(g)); gg <- sum(gr^2)
  coef <- as.numeric(crossprod(gr, Sr)) / gg
  res  <- Sr - outer(gr, coef); dfr <- nrow(S) - ncol(W) - 1
  sigma <- sqrt(colSums(res^2) / dfr); se <- sigma / sqrt(gg)
  d_hat <- coef / sigma
  # permutations: shuffle the label within fine x depth bin (10 quantile bins of log depth), recompute the residualised coefficient
  strata <- d[, paste(fine, cut(log_depth, unique(quantile(log_depth, 0:LOY_H2_N_DEPTH_BINS / LOY_H2_N_DEPTH_BINS)), include.lowest = TRUE))]
  idx <- split(seq_along(g), strata)
  perm_d <- matrix(NA_real_, n_perm, ncol(S))
  for (p in seq_len(n_perm)) {
    gp <- g; for (ix in idx) if (length(ix) > 1) gp[ix] <- g[ix][sample.int(length(ix))]
    gpr <- as.numeric(resid(gp)); cp <- as.numeric(crossprod(gpr, Sr)) / sum(gpr^2)
    perm_d[p, ] <- cp / sqrt(colSums((Sr - outer(gpr, cp))^2) / dfr)
  }
  p_perm <- (1 + colSums(abs(perm_d) >= matrix(abs(d_hat), n_perm, ncol(S), byrow = TRUE))) / (n_perm + 1)
  data.table(feature = colnames(S), n_A = sum(g == 1), n_B = sum(g == 0), coef = coef, sigma = sigma, d = d_hat,
             ci_lo = (coef - 1.96 * se) / sigma, ci_hi = (coef + 1.96 * se) / sigma, p_perm = p_perm)
}
sets  <- split(gs[role %in% c("primary", "artefact_control", "machinery"), gene], gs[role %in% c("primary", "artefact_control", "machinery"), set])
roles <- gs[, .(role = role[1]), by = set]
# c_amb from 02's table: Y_retained fraction among LEUK cells of the clonal -Y reference of this chemistry
cc02  <- fread(file.path(LOY_CALL_DIR, "loy_compartment_calls.csv"))[comp == "LEUK"]
camb_ref <- LOY_H2_C_AMB_REF[chem == LOY_EXT_ALL$chem_lit[match(lib$dataset, LOY_EXT_ALL$dataset)]]
camb_row <- cc02[dataset == camb_ref$dataset & Sample == camb_ref$Sample]
c_amb <- if (nrow(camb_row)) camb_row$n_ret / (camb_row$n_susp + camb_row$n_ret) else NA_real_
message(sprintf("[3] c_amb = %.3f (from %s)", c_amb, camb_ref$Sample))

run_arm <- function(d, label) {
  d <- copy(d)[order(cell)]
  S <- score_sets(X[keep, d$cell, drop = FALSE], sets[setdiff(names(sets), "F12_chrY")])
  if (length(chrY_present) >= 1) {                       # F12 from the chrY genes themselves (not in `keep`), machinery check only
    Y <- as.matrix(Matrix::t(X[chrY_present, d$cell, drop = FALSE])); S <- cbind(S, F12_chrY = rowMeans(Y))
  }
  r <- fit_arm(d, S, LOY_H2_N_PERM)
  r[, role := roles$role[match(feature, roles$set)]]
  r[, d_corr := d / (1 - c_amb)]
  fam <- r$role %in% c("primary", "artefact_control")
  r[, p_holm := NA_real_]; r[fam, p_holm := p.adjust(p_perm, "holm")]
  # sign consistency across eligible strata with >= LOY_H2_MIN_CELLS matched cells on both sides (leukaemic arm only)
  r[, arm := label]
  list(res = r, S = S, d = d)
}
message("[3] leukaemic arm ..."); L <- run_arm(leukS, "leukaemic")
message("[3] control arm ...");   C <- run_arm(ctrlS, "control")
# per-stratum d (leukaemic arm): the same residualised contrast inside each eligible fine label (log depth only)
per <- rbindlist(lapply(split(seq_len(nrow(L$d)), L$d$fine), function(ix) {
  d <- L$d[ix]; if (sum(d$group == "A") < LOY_H2_MIN_CELLS || sum(d$group == "B") < LOY_H2_MIN_CELLS) return(NULL)
  d[, fine := "one"]; r <- fit_arm(d, L$S[ix, , drop = FALSE], 0L); r[, fine := L$d$fine[ix[1]]]; r[, .(fine, feature, n_A, n_B, d)] }))
if (nrow(per)) { per[, case := case]; fwrite(per, out("perstratum")) }
sign_ok <- per[, .(n_strata = .N, sign_consistent = all(sign(d) == sign(d[1]))), by = feature]
L$res <- merge(L$res, sign_ok, by = "feature", all.x = TRUE)

## [4] gates ----------------------------------------------------------------------------------------------------
checks <- list()
k1_fail <- C$res[role %in% c("primary", "artefact_control") & p_holm < LOY_H2_ALPHA & abs(d) > LOY_H2_D_SMALL, feature]
checks[["K1"]] <- data.table(gate = "GATE_H2_K1", value = length(k1_fail), detail = paste(k1_fail, collapse = ","), pass = length(k1_fail) < 3)
k2 <- L$res[feature == "F12_chrY"]
checks[["K2"]] <- data.table(gate = "GATE_H2_K2", value = if (nrow(k2)) k2$d else NA_real_, detail = if (nrow(k2)) sprintf("chrY score d %.2f, p %.4f", k2$d, k2$p_perm) else "no chrY genes",
                             pass = nrow(k2) > 0 && k2$d < 0 && k2$p_perm < 0.001)
# GATE_N1: LOY_H2_N_SHUFFLE whole-arm shuffles with a short permutation run each; the per-feature rejection rate at alpha must sit in 0.03-0.07
n1 <- sapply(seq_len(LOY_H2_N_SHUFFLE), function(i) {
  d <- copy(L$d); d[, group := group[sample.int(.N)], by = .(fine, dbin = cut(log_depth, unique(quantile(log_depth, 0:LOY_H2_N_DEPTH_BINS / LOY_H2_N_DEPTH_BINS)), include.lowest = TRUE))]
  fit_arm(d, L$S, 199L)$p_perm < LOY_H2_ALPHA })
n1_rate <- rowMeans(n1); names(n1_rate) <- colnames(L$S)
checks[["N1"]] <- data.table(gate = "GATE_N1", value = mean(n1_rate), detail = paste(sprintf("%s=%.3f", names(n1_rate), n1_rate), collapse = " "),
                             pass = all(n1_rate >= 0.03 & n1_rate <= 0.07))
if (world == "planted") {   # GATE_P: plant d = LOY_H2_PLANT_D into a random gene set's score for the A cells; Holm within panel + planted set
  Sz <- L$S[, setdiff(colnames(L$S), "F12_chrY"), drop = FALSE]; gpool <- setdiff(keep, unlist(sets)); hits <- logical(LOY_H2_N_PLANTS)
  for (i in seq_len(LOY_H2_N_PLANTS)) {
    g60 <- sample(gpool, LOY_H2_PLANT_NGENES)
    s60 <- score_sets(X[g60, L$d$cell, drop = FALSE], list(planted = g60))
    Qw <- qr.Q(qr(design_of(L$d))); rs <- as.numeric(s60 - Qw %*% crossprod(Qw, s60))
    s60[L$d$group == "A"] <- s60[L$d$group == "A"] + LOY_H2_PLANT_D * sd(rs)          # shift of d residual SDs in the -Y cells
    r <- fit_arm(L$d, cbind(Sz, planted = as.numeric(s60)), 199L)
    hits[i] <- p.adjust(r$p_perm, "holm")[r$feature == "planted"] < LOY_H2_ALPHA
    if (i %% 50 == 0) message(sprintf("    plant %d/%d: power so far %.2f", i, LOY_H2_N_PLANTS, mean(hits[1:i])))
  }
  checks[["P"]] <- data.table(gate = "GATE_P", value = mean(hits), detail = sprintf("%d plants, d = %.1f, %d genes", LOY_H2_N_PLANTS, LOY_H2_PLANT_D, LOY_H2_PLANT_NGENES), pass = mean(hits) >= 0.80)
}
checks <- rbindlist(checks); checks[, `:=`(case = case, world = world)]
fwrite(checks, out("checks")); message("[4] checks:"); print(checks[, .(gate, value = round(value, 3), pass, detail = substr(detail, 1, 80))])

## [5] wording (section 8) and outputs ----------------------------------------------------------------------------
gate_p_pass <- if (world == "planted") checks[gate == "GATE_P", pass] else {   # the real world reads GATE_P from the planted smoke file
  f <- file.path(LOY_H2_DIR, sprintf("checks_%s_smoke_planted.csv", case)); gp <- if (file.exists(f)) fread(f)[gate == "GATE_P", pass] else NA; if (length(gp) != 1) NA else gp }
L$res[, k1_flag := feature %in% k1_fail]
L$res[, verdict := fifelse(!(role %in% c("primary", "artefact_control")), "machinery",
                   fifelse(k1_flag, "removed_by_K1",
                   fifelse(p_holm < LOY_H2_ALPHA & abs(d) >= LOY_H2_D_SMALL & (is.na(sign_consistent) | sign_consistent), "differs",
                   fifelse(ci_lo > -LOY_H2_D_SMALL & ci_hi < LOY_H2_D_SMALL & isTRUE(gate_p_pass), "no_difference", "not_resolved"))))]
L$res[, `:=`(case = case, world = world, c_amb = c_amb, gate_p_pass = gate_p_pass)]
C$res[, `:=`(case = case, world = world)]
fwrite(L$res[order(role, p_holm)], out("contrast")); fwrite(C$res[order(role, p_holm)], out("control"))
message("[5] leukaemic arm:"); print(L$res[, .(feature, role, n_A, n_B, d = round(d, 3), ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi), p_perm, p_holm = round(p_holm, 4), sign_consistent, verdict)])
message("    control arm (K1):"); print(C$res[, .(feature, d = round(d, 3), p_perm, p_holm = round(p_holm, 4))])

## [6] S1 genome-wide Wilcoxon on the matched set, calibrated by the control arm (real world only) -------------------
if (world == "real") {
  wil <- function(d) { Xd <- X[keep, d$cell, drop = FALSE]; a <- d$group == "A"
    p <- apply(as.matrix(Xd), 1, function(v) suppressWarnings(wilcox.test(v[a], v[!a])$p.value))
    data.table(gene = rownames(Xd), mean_A = Matrix::rowMeans(Xd[, a, drop = FALSE]), mean_B = Matrix::rowMeans(Xd[, !a, drop = FALSE]), p = p) }
  de_l <- wil(L$d); de_c <- wil(C$d)
  de_l[, diff := mean_A - mean_B]; de_l[, p_ctrl_le := vapply(p, function(x) mean(de_c$p <= x), 0)]       # fraction of control genes at least as extreme
  de_l[, efdr := pmin(1, p_ctrl_le / pmax(1 / .N, vapply(p, function(x) mean(de_l$p <= x), 0)))]        # empirical FDR = control rate / observed rate
  fwrite(de_l[order(p)], paste0(out("de"), ".gz"))
  message(sprintf("[6] S1: %d genes; %d with empirical FDR < 0.05", nrow(de_l), sum(de_l$efdr < 0.05)))
}
message(sprintf("[done] %s %s in %.1f min", case, world, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
