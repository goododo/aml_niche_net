#!/usr/bin/env Rscript
# 07_h2_lox_contrast.R ----
# input : LARGE1/loy/percell/GSE185381__y_percell.csv.gz (XIST / chrY UMIs, depth, labels per cell), the 04_annotated objects of
#         the case (PAUMTZ = donor PAXMIJ), the two control samples and AML005, scripts/loy/gene_sets_h2_v2.tsv, the gene order file
# output: results/tables/loy/05_lox/
#           stage0_cells_lox.csv.gz, stage0_counts_lox.csv   labels, matching, eligibility, misclassification bound (counts only)
#           contrast_lox<sfx>.csv     every feature: n, d, CI, permutation p, role; primary rows carry Holm p and the verdict
#           control_lox<sfx>.csv      the same pipeline on X-intact samples, where XIST-zero is dropout (GATE_LX_K1)
#           perstratum_lox<sfx>.csv   d per eligible fine label
#           checks_lox<sfx>.csv       GATE_LX_0, K1, K2, N1, P (P only in the planted world)
#           sens_lox.csv, de_lox.csv.gz, dose_lox.csv         real world only
#         <sfx> = "" for --world=real, "_smoke_null" / "_smoke_planted" otherwise
# what it does: the pre-registered contrast of PREREGISTRATION_LOX.md for the one female t(8;21) case with a subclonal -X
#   branch. -X = 0 XIST UMIs, X-retained = >= 2, inside a depth window, malignant-called cells without chrY UMIs; each
#   X-retained cell matched to up to 3 -X cells of the same fine label by depth; score ~ label + fine label + lane + log depth;
#   permutation p within fine label x depth bin; Holm over the two primary features (HLA-II expected higher, F7 expected lower).
# Usage : sbatch scripts/loy/submit_loy_h2.sbatch 07 --world=planted|null|real     (in that order; the tag prereg-lox-1-v1 comes first)
suppressPackageStartupMessages({ library(data.table); library(here); library(Matrix); library(Seurat); library(SeuratObject) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "loy", "config_loy.R")); source(here::here("scripts", "loy", "00_loy_common.R"))
args  <- commandArgs(trailingOnly = TRUE)
argv  <- function(k, default) { v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE)); if (length(v)) v[1] else default }
world <- argv("world", "real"); stopifnot(world %in% c("real", "null", "planted"))
sfx   <- if (world == "real") "" else paste0("_smoke_", world)
set.seed(LOY_H2_SEED); t0 <- Sys.time()                 # the matching below is identical in every world; the world seed is set after it
dir.create(LOY_LOX_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(stem, ext = ".csv") file.path(LOY_LOX_DIR, paste0(stem, "_lox", sfx, ext))
message(sprintf("[0] world %s | %s", world, format(Sys.time())))

## [1] labels and matching: XIST UMIs, chrY UMIs, depth and labels only ---------------------------------------------------
pc <- fread(file.path(LOY_PERCELL_DIR, paste0(LOY_LOX_DATASET, "__y_percell.csv.gz")),
            select = c("cell", "Sample", "Library", "bmm_fine", "anno_malignant", "nCount_RNA", "y_umi", "xist_umi"))
pc <- pc[Sample %in% c(LOY_LOX_CASE_SAMPLE, LOY_LOX_CTRL_SAMPLES, LOY_LOX_AMB_SAMPLE)]
pc[, leuk := !is.na(anno_malignant) & tolower(anno_malignant) == "malignant" & y_umi == 0]   # depositors' call; a chrY UMI = another donor of the pool
pc[, inwin := nCount_RNA >= LOY_LOX_DEPTH[1] & nCount_RNA < LOY_LOX_DEPTH[2]]
pc[, xs := fcase(xist_umi == 0, "zero", xist_umi >= LOY_LOX_XIST_RET, "ret", default = "one")]
pc[, `:=`(log_depth = log10(nCount_RNA), fine = bmm_fine)]
use <- pc[leuk & inwin & !is.na(fine) & fine != ""]
message(sprintf("[1] cells under the filter (malignant-called, no chrY UMI, %d <= depth < %d): %s", LOY_LOX_DEPTH[1], LOY_LOX_DEPTH[2],
                paste(sprintf("%s %d", use[, .N, by = Sample]$Sample, use[, .N, by = Sample]$N), collapse = " | ")))

# match_depth(): as 04 -- greedy nearest neighbour on log10 depth, without replacement, up to `ratio` partners within `caliper`
match_depth <- function(ld_s, ld_r, ratio, caliper) {
  ord <- order(ld_r); rs <- ld_r[ord]; used <- logical(length(rs)); res <- vector("list", length(ld_s))
  for (i in seq_along(ld_s)) {
    pos <- findInterval(ld_s[i], rs); lo <- pos; hi <- pos + 1L; picked <- integer()
    while (length(picked) < ratio) {
      while (lo >= 1L && used[lo]) lo <- lo - 1L
      while (hi <= length(rs) && used[hi]) hi <- hi + 1L
      dlo <- if (lo >= 1L) ld_s[i] - rs[lo] else Inf
      dhi <- if (hi <= length(rs)) rs[hi] - ld_s[i] else Inf
      if (min(dlo, dhi) > caliper) break
      if (dlo <= dhi) { picked <- c(picked, lo); used[lo] <- TRUE; lo <- lo - 1L } else { picked <- c(picked, hi); used[hi] <- TRUE; hi <- hi + 1L }
    }
    res[[i]] <- ord[picked]
  }
  res
}
match_stratum <- function(d, anchor) {     # one stratum; `anchor` = the smaller side, each of its cells gets up to 3 partners from the other side
  d <- copy(d); d[, `:=`(matched = FALSE, match_id = NA_integer_)]
  s <- which(d$group == anchor); r <- which(d$group != anchor)
  if (!length(s) || !length(r)) return(d)
  s <- s[sample.int(length(s))]                          # random processing order: no cell monopolises the nearest partners
  res <- match_depth(d$log_depth[s], d$log_depth[r], LOY_H2_MATCH_RATIO, LOY_H2_CALIPER)
  for (i in seq_along(s)) if (length(res[[i]])) set(d, i = c(s[i], r[res[[i]]]), j = c("matched", "match_id"), value = list(TRUE, i))
  d
}
# case: group A = -X (XIST 0), B = X-retained (XIST >= 2); anchors = B
cs <- use[Sample == LOY_LOX_CASE_SAMPLE & xs != "one"]; cs[, `:=`(group = fifelse(xs == "zero", "A", "B"), stratum = fine, arm = "case")]
cs <- rbindlist(lapply(split(cs, cs$stratum), match_stratum, anchor = "B"))
cnt_case <- cs[, .(n_A = sum(group == "A"), n_B = sum(group == "B"), matched_A = sum(matched & group == "A"), matched_B = sum(matched & group == "B"),
                   med_depth_A_m = as.numeric(median(nCount_RNA[matched & group == "A"])), med_depth_B_m = as.numeric(median(nCount_RNA[matched & group == "B"]))), by = .(arm, stratum)]
cnt_case[, eligible := matched_A >= LOY_H2_MIN_CELLS & matched_B >= LOY_H2_MIN_CELLS]
elig <- cnt_case[eligible == TRUE, stratum]
# control: X-intact samples, group A = XIST 0 (dropout), B = XIST >= 2; anchors = A; stratum = sample x fine label
ct <- use[Sample %in% LOY_LOX_CTRL_SAMPLES & xs != "one" & fine %in% elig]; ct[, `:=`(group = fifelse(xs == "zero", "A", "B"), stratum = paste(Sample, fine, sep = "::"), arm = "control")]
nz <- ct[, .(n_zero = sum(group == "A")), by = Sample]; message("    control XIST-zero cells in the eligible labels: ", paste(sprintf("%s %d", nz$Sample, nz$n_zero), collapse = " | "))
stopifnot(all(LOY_LOX_CTRL_SAMPLES %in% nz$Sample), all(nz$n_zero >= LOY_LOX_CTRL_MIN_ZERO))
ct <- rbindlist(lapply(split(ct, ct$stratum), match_stratum, anchor = "A"))
cnt_ctrl <- ct[, .(n_A = sum(group == "A"), n_B = sum(group == "B"), matched_A = sum(matched & group == "A"), matched_B = sum(matched & group == "B"),
                   med_depth_A_m = as.numeric(median(nCount_RNA[matched & group == "A"])), med_depth_B_m = as.numeric(median(nCount_RNA[matched & group == "B"]))), by = .(arm, stratum)]
cnt_ctrl[, eligible := matched_A >= LOY_LOX_CTRL_STRATUM_MIN]
# misclassification bound (section 10 item 12), counts only
k  <- use[fine %in% elig]
e  <- k[Sample %in% LOY_LOX_CTRL_SAMPLES, mean(xs == "zero")]; r2 <- k[Sample %in% LOY_LOX_CTRL_SAMPLES, mean(xs == "ret")]
a  <- k[Sample == LOY_LOX_AMB_SAMPLE, mean(xs == "ret")]
n0 <- k[Sample == LOY_LOX_CASE_SAMPLE, sum(xs == "zero")]; n2 <- k[Sample == LOY_LOX_CASE_SAMPLE, sum(xs == "ret")]
c_A <- (n2 / r2) * e / n0; c_B <- n0 * a / n2; c_mis <- c_A + c_B
stage_counts <- rbind(cnt_case, cnt_ctrl); stage_counts[, `:=`(e_dropout = e, a_ambient = a, c_A = c_A, c_B = c_B, c_mis = c_mis)]
fwrite(stage_counts, file.path(LOY_LOX_DIR, "stage0_counts_lox.csv"))
fwrite(rbind(cs, ct)[, .(cell, arm, Sample, Library, fine, stratum, group, matched, match_id, nCount_RNA, log_depth, xist_umi)], file.path(LOY_LOX_DIR, "stage0_cells_lox.csv.gz"))
message("[1] strata:"); print(stage_counts[, .(arm, stratum, n_A, n_B, matched_A, matched_B, med_depth_A_m, med_depth_B_m, eligible)])
message(sprintf("    misclassification: dropout rate e %.3f, ambient rate a %.4f -> c_A %.3f, c_B %.3f, c %.3f", e, a, c_A, c_B, c_mis))
caseS <- cs[matched == TRUE & stratum %in% elig][order(cell)]
ctrlS <- ct[matched == TRUE & stratum %in% cnt_ctrl[eligible == TRUE, stratum]][order(cell)]
# sanity check: matched halves sit at the same depth (hard stop at half the caliper, as 04)
dr <- caseS[, .(ratio = median(nCount_RNA[group == "A"]) / median(nCount_RNA[group == "B"])), by = stratum]; print(dr)
stopifnot(all(abs(log10(dr$ratio)) <= LOY_H2_CALIPER / 2))
gate0 <- length(elig) >= 2 && sum(caseS$group == "B") >= LOY_LOX_MIN_RET
message(sprintf("    GATE_LX_0: %d eligible labels, matched -X %d / X-retained %d -> %s | control matched zero %d / retained %d in %d strata",
                length(elig), sum(caseS$group == "A"), sum(caseS$group == "B"), if (gate0) "pass" else "FAIL", sum(ctrlS$group == "A"), sum(ctrlS$group == "B"), uniqueN(ctrlS$stratum)))
stopifnot(gate0)

dbin_of <- function(d) cut(d$log_depth, unique(quantile(d$log_depth, 0:LOY_H2_N_DEPTH_BINS / LOY_H2_N_DEPTH_BINS)), include.lowest = TRUE)   # 10 depth bins over the set
shuffle_within <- function(d) { d <- copy(d); d[, dbin := dbin_of(d)]; d[, group := group[sample.int(.N)], by = .(stratum, dbin)]; d[, dbin := NULL]; d }
set.seed(LOY_H2_SEED + nchar(world))
if (world != "real") {   # smoke worlds: the XIST label of BOTH arms is destroyed here, before any expression value is loaded
  caseS <- shuffle_within(caseS); ctrlS <- shuffle_within(ctrlS); message("    labels shuffled within stratum x depth bin (smoke world)")
}

## [2] counts ------------------------------------------------------------------------------------------------------------
read_counts <- function(smp) { s <- readRDS(file.path(LOY_ANNOT_DIR, LOY_LOX_DATASET, paste0(smp, ".rds"))); m <- counts_of(s); rm(s); invisible(gc()); m }
lognorm <- function(m) { sf <- 1e4 / Matrix::colSums(m); m@x <- log1p(m@x * rep.int(sf, diff(m@p))); m }   # log1p(counts / depth * 1e4), per column
m <- read_counts(LOY_LOX_CASE_SAMPLE); all_leuk <- pc[Sample == LOY_LOX_CASE_SAMPLE & leuk == TRUE, cell]
stopifnot(all(caseS$cell %in% colnames(m)), all(all_leuk %in% colnames(m)))
det <- Matrix::rowMeans(m[, all_leuk, drop = FALSE] > 0)                              # detection among the donor's leukaemic cells, label not looked at
raw_case <- m[, caseS$cell, drop = FALSE]; X <- lognorm(raw_case); rm(m); invisible(gc())
mc <- do.call(cbind, lapply(LOY_LOX_CTRL_SAMPLES, function(smp) { mm <- read_counts(smp); mm[, intersect(colnames(mm), ctrlS$cell), drop = FALSE] }))
stopifnot(all(ctrlS$cell %in% colnames(mc)), identical(rownames(mc), rownames(raw_case)))
raw_ctrl <- mc[, ctrlS$cell, drop = FALSE]; Xc <- lognorm(raw_ctrl); rm(mc); invisible(gc())
gs   <- fread(LOY_H2_GENESETS); chrY <- gs[set == "chrY_exclude", gene]; par <- gs[set == "PAR_exclude", gene]
gord <- fread(LOY_H2_GENE_ORDER, header = FALSE, col.names = c("gene", "chr", "start", "end")); chrX <- union(gord[chr == "chrX", gene], par)
keep_s1 <- names(det)[det >= LOY_H2_MIN_DETECT & !names(det) %in% chrY & names(det) != LOY_XIST]   # S1: chrX stays in, flagged
keep    <- setdiff(keep_s1, chrX)                                                                   # feature scores: no chrX, no chrY gene
par_keep <- intersect(par, keep_s1); xnp_keep <- setdiff(intersect(chrX, keep_s1), par)
message(sprintf("[2] genes: %d detected in >= %.0f%% of %d leukaemic cells; %d autosomal for the scores; %d PAR and %d other chrX genes for the label check",
                length(keep_s1), 100 * LOY_H2_MIN_DETECT, length(all_leuk), length(keep), length(par_keep), length(xnp_keep)))

## [3] scores and the model -------------------------------------------------------------------------------------------------
score_sets <- function(Xm, sets) {      # as 05: z-score each gene over these cells, mean over the set's genes present; sets with < 5 genes dropped
  Xd <- as.matrix(Matrix::t(Xm)); mu <- colMeans(Xd); sdv <- apply(Xd, 2, sd); sdv[sdv == 0] <- 1
  Z <- sweep(sweep(Xd, 2, mu), 2, sdv, "/")
  sc <- lapply(sets, function(g) { g <- intersect(g, colnames(Z)); if (length(g) < 5) return(NULL); rowMeans(Z[, g, drop = FALSE]) })
  sc <- sc[!vapply(sc, is.null, TRUE)]; do.call(cbind, sc)
}
design_of <- function(d) {              # stratum dummies + lane dummies + log depth; a one-level factor is left out (it breaks model.matrix)
  dd <- d[, .(stratum = factor(stratum), lane = factor(Library), log_depth)]
  model.matrix(reformulate(c(if (nlevels(dd$stratum) > 1) "stratum", if (nlevels(dd$lane) > 1) "lane", "log_depth")), data = dd)
}
fit_arm <- function(d, S, n_perm, Zextra = NULL) {   # as 05 (Frisch-Waugh residualisation), with the design reduced to its rank
  g <- as.numeric(d$group == "A"); W <- design_of(d); if (!is.null(Zextra)) W <- cbind(W, Zextra)
  q <- qr(W); Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]; resid <- function(v) v - Q %*% crossprod(Q, v)   # first `rank` columns of Q span the design
  Sr <- resid(S); gr <- as.numeric(resid(g)); gg <- sum(gr^2); coef <- as.numeric(crossprod(gr, Sr)) / gg
  res <- Sr - outer(gr, coef); dfr <- nrow(S) - q$rank - 1; sigma <- sqrt(colSums(res^2) / dfr); se <- sigma / sqrt(gg); d_hat <- coef / sigma
  idx <- split(seq_along(g), paste(d$stratum, dbin_of(d))); perm_d <- matrix(NA_real_, n_perm, ncol(S))
  for (p in seq_len(n_perm)) { gp <- g; for (ix in idx) if (length(ix) > 1) gp[ix] <- g[ix][sample.int(length(ix))]
    gpr <- as.numeric(resid(gp)); cp <- as.numeric(crossprod(gpr, Sr)) / sum(gpr^2); perm_d[p, ] <- cp / sqrt(colSums((Sr - outer(gpr, cp))^2) / dfr) }
  p_perm <- if (n_perm > 0) (1 + colSums(abs(perm_d) >= matrix(abs(d_hat), n_perm, ncol(S), byrow = TRUE))) / (n_perm + 1) else rep(NA_real_, ncol(S))
  data.table(feature = colnames(S), n_A = sum(g == 1), n_B = sum(g == 0), coef = coef, sigma = sigma, d = d_hat,
             ci_lo = (coef - 1.96 * se) / sigma, ci_hi = (coef + 1.96 * se) / sigma, p_perm = p_perm)
}
f7 <- gs[set == "F7_TNFA_NFKB", gene]; mp9 <- gs[set == "MP9", gene]
panel <- gs[!set %in% c("chrY_exclude", "PAR_exclude", "F12_chrY", "F7_TNFA_NFKB")]
sets  <- c(list(P1_HLA2 = LOY_LOX_HLA2, P2_F7 = f7, F7_noMP9 = setdiff(f7, mp9), NFKB_core = LOY_LOX_NFKB_CORE, MHC1 = LOY_LOX_MHC1), split(panel$gene, panel$set))
role_of <- function(f) fifelse(f %in% c("P1_HLA2", "P2_F7"), "primary", fifelse(f %in% c("PAR", "X_nonPAR"), "machinery", "descriptive"))
run_arm <- function(d, Xm, label) {
  S <- score_sets(Xm[intersect(keep, unlist(sets)), d$cell, drop = FALSE], sets)
  M <- score_sets(Xm[c(par_keep, xnp_keep), d$cell, drop = FALSE], list(PAR = par_keep, X_nonPAR = xnp_keep))   # label check: chrX dose, never a result
  S <- cbind(S, M); r <- fit_arm(d, S, LOY_H2_N_PERM); r[, `:=`(role = role_of(feature), arm = label)]
  list(res = r, S = S, d = d)
}
message("[3] case arm ...");    L <- run_arm(caseS, X,  "case")
message("[3] control arm ..."); C <- run_arm(ctrlS, Xc, "control")
stopifnot(all(c("P1_HLA2", "P2_F7", "F7_noMP9", "MP9", "PAR") %in% L$res$feature))
per <- rbindlist(lapply(split(seq_len(nrow(L$d)), L$d$stratum), function(ix) {   # per fine label, for the sign rule
  d <- L$d[ix]; if (sum(d$group == "A") < LOY_H2_MIN_CELLS || sum(d$group == "B") < LOY_H2_MIN_CELLS) return(NULL)
  r <- fit_arm(d, L$S[ix, , drop = FALSE], 0L); r[, stratum := d$stratum[1]]; r[, .(stratum, feature, n_A, n_B, d)] }))
fwrite(per, out("perstratum"))
sign_ok <- per[, .(n_strata = .N, sign_consistent = all(sign(d) == sign(d[1]))), by = feature]

## [4] gates -------------------------------------------------------------------------------------------------------------------
checks <- list()
checks[["LX0"]] <- data.table(gate = "GATE_LX_0", value = sum(caseS$group == "B"), detail = sprintf("%d eligible labels: %s", length(elig), paste(elig, collapse = ", ")), pass = gate0)
k2 <- L$res[feature == "PAR"]
checks[["K2"]] <- data.table(gate = "GATE_LX_K2", value = k2$d, detail = sprintf("PAR score (%d genes) d %.3f, p %.4f; other chrX (%d genes) d %.3f", length(par_keep), k2$d, k2$p_perm, length(xnp_keep), L$res[feature == "X_nonPAR", d]),
                             pass = k2$d < 0 && k2$p_perm < 0.001)
mp9_moves <- abs(L$res[feature == "MP9", d]) >= LOY_H2_D_SMALL
prim_feats <- c("P1_HLA2", if (mp9_moves) "F7_noMP9" else "P2_F7")      # section 10 item 10
checks[["MP9"]] <- data.table(gate = "MP9_clause", value = L$res[feature == "MP9", d], detail = sprintf("second primary feature = %s", prim_feats[2]), pass = !mp9_moves)
Sp <- L$S[, c("P1_HLA2", "P2_F7"), drop = FALSE]
n1 <- sapply(seq_len(LOY_H2_N_SHUFFLE), function(i) fit_arm(shuffle_within(L$d), Sp, 199L)$p_perm < LOY_H2_ALPHA)   # 200 shuffles x 199 permutations
n1_rate <- rowMeans(n1); names(n1_rate) <- colnames(Sp)
checks[["N1"]] <- data.table(gate = "GATE_N1", value = mean(n1_rate), detail = paste(sprintf("%s=%.3f", names(n1_rate), n1_rate), collapse = " "),
                             pass = all(n1_rate >= 0.020 & n1_rate <= 0.080) && mean(n1_rate) >= 0.03 && mean(n1_rate) <= 0.07)
if (world == "planted") {   # GATE_P: d = 0.3 planted into a random 60-gene set in the (shuffled) -X cells; Holm over the two primary features + the planted set
  gpool <- setdiff(keep, unlist(sets)); hits <- logical(LOY_H2_N_PLANTS); Qw <- { q <- qr(design_of(L$d)); qr.Q(q)[, seq_len(q$rank), drop = FALSE] }
  for (i in seq_len(LOY_H2_N_PLANTS)) {
    g60 <- sample(gpool, LOY_H2_PLANT_NGENES); s60 <- as.numeric(score_sets(X[g60, L$d$cell, drop = FALSE], list(planted = g60)))
    rs <- as.numeric(s60 - Qw %*% crossprod(Qw, s60)); s60[L$d$group == "A"] <- s60[L$d$group == "A"] + LOY_H2_PLANT_D * sd(rs)
    r <- fit_arm(L$d, cbind(Sp, planted = s60), LOY_H2_N_PERM); hits[i] <- p.adjust(r$p_perm, "holm")[r$feature == "planted"] < LOY_H2_ALPHA
    if (i %% 50 == 0) message(sprintf("    plant %d/%d: power so far %.2f", i, LOY_H2_N_PLANTS, mean(hits[1:i])))
  }
  checks[["P"]] <- data.table(gate = "GATE_P", value = mean(hits), detail = sprintf("%d plants, d = %.1f, %d genes", LOY_H2_N_PLANTS, LOY_H2_PLANT_D, LOY_H2_PLANT_NGENES), pass = mean(hits) >= 0.80)
}
gate_p_pass <- if (world == "planted") checks[["P"]]$pass else {   # the other worlds read GATE_P from the planted smoke file
  f <- file.path(LOY_LOX_DIR, "checks_lox_smoke_planted.csv"); gp <- if (file.exists(f)) fread(f)[gate == "GATE_P", pass] else NA; if (length(gp) != 1) NA else gp }

## [5] wording (section 7) -------------------------------------------------------------------------------------------------------
res <- merge(L$res, sign_ok, by = "feature", all.x = TRUE)
res <- merge(res, C$res[, .(feature, d_ctrl = d, p_ctrl = p_perm, n_A_ctrl = n_A, n_B_ctrl = n_B)], by = "feature", all.x = TRUE)
res[, `:=`(d_corr = d / (1 - c_mis), c_mis = c_mis, in_family = feature %in% prim_feats, expected = LOY_LOX_EXPECT[feature], p_holm = NA_real_)]
res[in_family == TRUE, p_holm := p.adjust(p_perm, "holm")]
res[, k1_removed := in_family & !is.na(d_ctrl) & p_ctrl < LOY_H2_ALPHA & abs(d_ctrl) >= LOY_H2_D_SMALL & sign(d_ctrl) == sign(d)]
res[, verdict := NA_character_]
res[in_family == TRUE, verdict := fcase(
  k1_removed, "removed_by_K1",
  p_holm < LOY_H2_ALPHA & sign(d) == expected & abs(d) >= LOY_H2_D_SMALL & (is.na(sign_consistent) | sign_consistent), "repeats",
  p_holm < LOY_H2_ALPHA & sign(d) == expected & abs(d) <  LOY_H2_D_SMALL, "same_direction_small",
  p_holm < LOY_H2_ALPHA & sign(d) != expected, "opposite",
  ci_lo > -LOY_H2_D_SMALL & ci_hi < LOY_H2_D_SMALL & isTRUE(gate_p_pass), "no_difference",
  default = "not_resolved")]
res[, `:=`(world = world, gate_p_pass = gate_p_pass)]
checks[["K1"]] <- data.table(gate = "GATE_LX_K1", value = sum(res$k1_removed), detail = paste(sprintf("%s ctrl d %.3f p %.3f", prim_feats, res[match(prim_feats, feature), d_ctrl], res[match(prim_feats, feature), p_ctrl]), collapse = "; "),
                             pass = !any(res$k1_removed))
checks <- rbindlist(checks, fill = TRUE); checks[, world := world]
fwrite(checks, out("checks")); fwrite(res[order(-in_family, role, feature)], out("contrast")); fwrite(C$res[order(role, feature)], out("control"))
message("[4] checks:"); print(checks[, .(gate, value = round(value, 3), pass, detail = substr(detail, 1, 110))])
message("[5] primary family:"); print(res[in_family == TRUE, .(feature, n_A, n_B, d = round(d, 3), ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi), d_corr = round(d_corr, 3), p_perm, p_holm = round(p_holm, 4),
                                                                  sign_consistent, d_ctrl = round(d_ctrl, 3), p_ctrl, verdict)])
message("    all features, case arm (control d beside it):"); print(res[order(role, feature), .(feature, role, d = round(d, 3), ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi), p_perm, d_ctrl = round(d_ctrl, 3), p_ctrl)])
message("    per fine label, primary features:"); print(dcast(per[feature %in% c("P1_HLA2", "P2_F7", "F7_noMP9", "PAR")], stratum + n_A + n_B ~ feature, value.var = "d")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])

## [6] real world only: sensitivity runs, genome-wide table S1, dose readout ----------------------------------------------------
if (world == "real") {
  pf <- c("P1_HLA2", "P2_F7", "F7_noMP9"); sens <- list()
  ix <- L$d$stratum != "CD14 Mono";      sens[["a_no_CD14_Mono"]]  <- fit_arm(L$d[ix], L$S[ix, pf, drop = FALSE], LOY_H2_N_PERM)
  ix <- L$d$nCount_RNA < 8000;           sens[["b_depth_lt_8000"]] <- fit_arm(L$d[ix], L$S[ix, pf, drop = FALSE], LOY_H2_N_PERM)
  cc <- intersect(c("F4_S", "F4_G2M"), colnames(L$S)); sens[["c_plus_cell_cycle"]] <- fit_arm(L$d, L$S[, pf, drop = FALSE], LOY_H2_N_PERM, L$S[, cc, drop = FALSE])
  sens <- rbindlist(sens, idcol = "run"); fwrite(sens, out("sens")); message("[6] sensitivity:"); print(sens[, .(run, feature, n_A, n_B, d = round(d, 3), ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi), p_perm)])
  wil <- function(d, Xm) { Xd <- Xm[keep_s1, d$cell, drop = FALSE]; ia <- d$group == "A"
    p <- apply(as.matrix(Xd), 1, function(v) suppressWarnings(wilcox.test(v[ia], v[!ia])$p.value))
    data.table(gene = rownames(Xd), mean_A = Matrix::rowMeans(Xd[, ia, drop = FALSE]), mean_B = Matrix::rowMeans(Xd[, !ia, drop = FALSE]), p = p) }
  de_l <- wil(L$d, X); de_c <- wil(C$d, Xc)
  de_l[, `:=`(diff = mean_A - mean_B, chrX = gene %in% chrX, par = gene %in% par)]
  pcn <- de_c$p[!is.na(de_c$p)]; pl <- de_l$p[!is.na(de_l$p)]                          # control arm = empirical null, as 05
  de_l[, p_ctrl_le := vapply(p, function(x) if (is.na(x)) NA_real_ else mean(pcn <= x), 0)]
  de_l[, efdr := pmin(1, p_ctrl_le / pmax(1 / length(pl), vapply(p, function(x) if (is.na(x)) NA_real_ else mean(pl <= x), 0)))]
  de_l <- merge(de_l, de_c[, .(gene, diff_ctrl = mean_A - mean_B, p_ctrl = p)], by = "gene", all.x = TRUE)
  fwrite(de_l[order(p)], paste0(out("de"), ".gz"))
  message(sprintf("    S1: %d genes; %d with empirical FDR < 0.05, of which %d on chrX", nrow(de_l), sum(de_l$efdr < 0.05, na.rm = TRUE), sum(de_l$efdr < 0.05 & de_l$chrX, na.rm = TRUE)))
  dose <- function(d, raw, label) { g <- intersect(unique(c(LOY_LOX_DOSE_GENES, par)), rownames(raw)); ia <- d$group == "A"; lib <- Matrix::colSums(raw); rr <- raw[g, , drop = FALSE]
    data.table(arm = label, gene = g, n_A = sum(ia), n_B = sum(!ia), det_A = Matrix::rowMeans(rr[, ia, drop = FALSE] > 0), det_B = Matrix::rowMeans(rr[, !ia, drop = FALSE] > 0),
               umi10k_A = as.numeric(rr[, ia, drop = FALSE] %*% (1e4 / lib[ia])) / sum(ia), umi10k_B = as.numeric(rr[, !ia, drop = FALSE] %*% (1e4 / lib[!ia])) / sum(!ia)) }
  dz <- rbind(dose(L$d, raw_case, "case"), dose(C$d, raw_ctrl, "control")); dz[, ratio := umi10k_A / umi10k_B]; fwrite(dz, out("dose"))
  message("    dose (ratio = XIST-zero / XIST>=2 UMIs per 10k; the control arm shows what dropout alone does):")
  print(dcast(dz[gene %in% LOY_LOX_DOSE_GENES | (det_A >= 0.05 | det_B >= 0.05)], gene ~ arm, value.var = c("det_A", "det_B", "ratio"))[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])
}
message(sprintf("[done] %s in %.1f min", world, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
