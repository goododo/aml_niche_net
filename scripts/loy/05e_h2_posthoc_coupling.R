# 05e_h2_posthoc_coupling.R
# input : results/tables/loy/04_h2/stage0_cells_<case>.csv.gz, stage0_counts.csv, the case's Seurat object,
#         scripts/loy/gene_sets_h2_v2.tsv (F7_TNFA_NFKB, MP9, chrY_exclude, PAR_exclude)
# output: results/tables/loy/04_h2/posthoc_f7_variants.csv     case x feature x design: d of the -Y label, permutation p
#         results/tables/loy/04_h2/posthoc_coupling.csv        case x cell subset x pair: per-cell partial correlation + random-set null
#         results/tables/loy/04_h2/posthoc_dosage_detection.csv case x gene: raw detection and UMI/10k in -Y vs Y+ (PAR1, X-Y pairs)
# what it does: three post-hoc checks (NOT pre-registered; asked at the 2026-11-26 discussion) on the leukaemic matched set of
#   each H2 case. (1) Is "TNF/NF-kB lower" carried by stress / immediate-early genes? F7 is re-scored without the genes it
#   shares with MP9 (dissociation stress), with only those genes, and on a hand-listed NF-kB core. (2) Are "HLA-II higher"
#   and "NF-kB lower" the same cells? Per-cell partial correlation of the two scores given fine label + depth (+ label),
#   and each effect re-fitted with the other score as a covariate. (3) Raw detection of the PAR1 genes (CSF2RA sits below
#   the 5 % filter of 05) and of the X-Y homologue pairs (KDM6A/UTY ...), straight from the counts.
# Usage : sbatch scripts/loy/submit_loy_h2.sbatch 05e

suppressPackageStartupMessages({ library(data.table); library(here); library(Matrix); library(Seurat); library(SeuratObject) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "loy", "config_loy.R")); source(here::here("scripts", "loy", "00_loy_common.R"))
set.seed(LOY_H2_SEED); N_PERM <- 999L; N_NULL <- 200L; CASES <- c("AML15", "AML13", "PT43"); t0 <- Sys.time()

HLA2_STRUCT <- c("HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DMA","HLA-DMB","HLA-DOA","CD74")   # as 05c
# NF-kB subunits, inhibitors and textbook targets that are in F7 and NOT in MP9; listed by hand on 2026-11-26 after the S1 tables
# were seen, so this set is post-hoc twice over and can only describe.
NFKB_CORE <- c("NFKB1","NFKB2","REL","RELA","RELB","NFKBIE","TNIP1","BIRC2","BIRC3","TRAF1","CD83","ICAM1","CXCL2","CXCL3",
               "SOD2","BCL2A1","TNFAIP2","TNFAIP8","BCL3")
XY_PAIRS  <- c("KDM6A","UTY","KDM5C","KDM5D","DDX3X","DDX3Y","EIF1AX","EIF1AY","ZFX","ZFY","RPS4X","RPS4Y1","USP9X","USP9Y")

score_sets <- function(Xc, sets) {                                           # as 05: z-score each gene over these cells, mean over the set (>= 5 present)
  Xd <- as.matrix(Matrix::t(Xc)); mu <- colMeans(Xd); sdv <- apply(Xd, 2, sd); sdv[sdv == 0] <- 1
  Z <- sweep(sweep(Xd, 2, mu), 2, sdv, "/")
  sc <- lapply(sets, function(g) { g <- intersect(g, colnames(Z)); if (length(g) < 5) return(NULL); rowMeans(Z[, g, drop = FALSE]) })
  do.call(cbind, sc[!vapply(sc, is.null, TRUE)])
}
base_design <- function(d) model.matrix(if (uniqueN(d$fine) > 1) ~ fine + log_depth else ~ log_depth, data = d[, .(fine = factor(fine), log_depth)])
fit_arm <- function(d, S, n_perm, Zextra = NULL) {                           # as 05c: Frisch-Waugh, permutation within fine x depth bin
  g <- as.numeric(d$group == "A"); base <- base_design(d)
  W <- if (!is.null(Zextra)) cbind(base, Zextra) else base
  Q <- qr.Q(qr(W)); resid <- function(v) v - Q %*% crossprod(Q, v)             # qr.Q = orthonormal basis of the design; v - QQ'v = residual
  Sr <- resid(S); gr <- as.numeric(resid(g)); gg <- sum(gr^2); coef <- as.numeric(crossprod(gr, Sr)) / gg
  res <- Sr - outer(gr, coef); dfr <- nrow(S) - ncol(W) - 1; sigma <- sqrt(colSums(res^2) / dfr); d_hat <- coef / sigma
  strata <- d[, paste(fine, cut(log_depth, unique(quantile(log_depth, 0:LOY_H2_N_DEPTH_BINS / LOY_H2_N_DEPTH_BINS)), include.lowest = TRUE))]
  idx <- split(seq_along(g), strata); perm_d <- matrix(NA_real_, n_perm, ncol(S))
  for (p in seq_len(n_perm)) { gp <- g; for (ix in idx) if (length(ix) > 1) gp[ix] <- g[ix][sample.int(length(ix))]
    gpr <- as.numeric(resid(gp)); cp <- as.numeric(crossprod(gpr, Sr)) / sum(gpr^2); perm_d[p, ] <- cp / sqrt(colSums((Sr - outer(gpr, cp))^2) / dfr) }
  p_perm <- (1 + colSums(abs(perm_d) >= matrix(abs(d_hat), n_perm, ncol(S), byrow = TRUE))) / (n_perm + 1)
  data.table(feature = colnames(S), n_A = sum(g == 1), n_B = sum(g == 0), d = d_hat, p_perm = p_perm)
}
# pcor(): correlation of two per-cell scores after both are residualised on fine label + log depth (+ the -Y label when the
# subset holds both groups). Fisher z interval. A positive value can be purely technical (cell quality lifts every gene),
# which is why the same number is computed against random gene sets below.
pcor <- function(d, a, b) {
  W <- base_design(d); if (uniqueN(d$group) > 1) W <- cbind(W, g = as.numeric(d$group == "A"))
  Q <- qr.Q(qr(W)); ra <- as.numeric(a - Q %*% crossprod(Q, a)); rb <- as.numeric(b - Q %*% crossprod(Q, b))
  r <- cor(ra, rb); se <- 1 / sqrt(length(ra) - ncol(W) - 3)
  c(r = r, lo = tanh(atanh(r) - 1.96 * se), hi = tanh(atanh(r) + 1.96 * se))
}

gs <- fread(LOY_H2_GENESETS); chrY <- gs[set == "chrY_exclude", gene]; par <- gs[set == "PAR_exclude", gene]
f7 <- gs[set == "F7_TNFA_NFKB", gene]; mp9 <- gs[set == "MP9", gene]
sets <- list(HLA2 = HLA2_STRUCT, F7_all = f7, F7_noMP9 = setdiff(f7, mp9), F7_MP9shared = intersect(f7, mp9), NFKB_core = NFKB_CORE, MP9 = mp9)
stopifnot(!length(intersect(HLA2_STRUCT, f7)), !length(intersect(NFKB_CORE, mp9)))     # sanity: the two axes share no gene; the core is stress-free
message(sprintf("[0] set sizes: %s", paste(sprintf("%s=%d", names(sets), lengths(sets)), collapse = "  ")))

run_case <- function(cs, s_full) {
  lib <- LOY_H2_CASES[case == cs & role == "contrast"]
  cl <- fread(file.path(LOY_H2_DIR, sprintf("stage0_cells_%s.csv.gz", cs)))
  elig <- fread(file.path(LOY_H2_DIR, "stage0_counts.csv"))[case == cs & eligible == TRUE]
  leukS <- cl[arm == "leukaemic" & fine != "ALL" & matched == TRUE & fine %in% elig[arm == "leukaemic", fine]]
  all_leuk_cells <- cl[arm == "leukaemic" & fine == "ALL", cell]; need <- unique(c(leukS$cell, all_leuk_cells))
  if (lib$dataset == "GSE269669") {
    mats <- lapply(c("avm_aml", "avm_tme"), function(obj) { s <- read_gse269669(obj); hit <- intersect(colnames(s), need)
      m <- if (length(hit)) counts_of(subset(s, cells = hit)) else NULL; rm(s); gc(); m })
    mats <- mats[!vapply(mats, is.null, TRUE)]; gmn <- Reduce(intersect, lapply(mats, rownames)); cnt <- do.call(cbind, lapply(mats, function(m) m[gmn, , drop = FALSE])); rm(mats); gc()
  } else { s <- subset(s_full, cells = intersect(colnames(s_full), need)); cnt <- counts_of(s); rm(s); gc() }
  stopifnot(all(need %in% colnames(cnt)))
  d <- copy(leukS)[order(cell)]; cells <- d$cell
  lognorm <- function(m) { sf <- 1e4 / Matrix::colSums(m); m@x <- log1p(m@x * rep.int(sf, diff(m@p))); m }   # log1p(counts / depth * 1e4)
  det <- Matrix::rowMeans(cnt[, all_leuk_cells, drop = FALSE] > 0)
  keep <- setdiff(names(det)[det >= LOY_H2_MIN_DETECT & !names(det) %in% chrY], par)                    # the gene filter of 05, unchanged
  libsize <- Matrix::colSums(cnt[, cells, drop = FALSE])
  dos_genes <- intersect(unique(c(par, XY_PAIRS)), rownames(cnt)); raw <- cnt[dos_genes, cells, drop = FALSE]   # raw counts, no filter
  X <- lognorm(cnt[, cells, drop = FALSE]); rm(cnt); gc()

  ## (1) F7 variants and (2b) mutual adjustment ------------------------------------------------------------
  S <- score_sets(X[intersect(keep, unlist(sets)), cells, drop = FALSE], sets)
  n_used <- vapply(sets[colnames(S)], function(g) length(intersect(g, keep)), 1L)
  r0 <- fit_arm(d, S, N_PERM)[, design := "baseline"]
  r1 <- fit_arm(d, S[, setdiff(colnames(S), "HLA2"), drop = FALSE], N_PERM, S[, "HLA2", drop = FALSE])[, design := "plus_HLA2"]
  r2 <- fit_arm(d, S[, "HLA2", drop = FALSE], N_PERM, S[, "F7_all", drop = FALSE])[, design := "plus_F7_all"]
  fv <- rbind(r0, r1, r2)[, `:=`(case = cs, n_genes = n_used[feature])]

  ## (2a) per-cell coupling, with a random-set null -----------------------------------------------------------
  pool <- setdiff(keep, unlist(sets)); nullS <- sapply(seq_len(N_NULL), function(i) rowMeans(score_sets(X[sample(pool, 60L), cells, drop = FALSE], list(x = pool))))
  # score_sets() with the pool as the "set" simply averages the z-scores of the 60 sampled genes
  subsets <- list(all = rep(TRUE, nrow(d)), minusY = d$group == "A", Yplus = d$group == "B")
  cp <- rbindlist(lapply(names(subsets), function(nm) { ix <- subsets[[nm]]; dd <- d[ix]
    nl <- apply(nullS[ix, , drop = FALSE], 2, function(z) pcor(dd, S[ix, "HLA2"], z)["r"])       # HLA2 vs each random set, same subset
    rbindlist(lapply(setdiff(colnames(S), "HLA2"), function(f) {
      v <- pcor(dd, S[ix, "HLA2"], S[ix, f])
      data.table(case = cs, subset = nm, n = sum(ix), pair = paste0("HLA2~", f), r = v["r"], lo = v["lo"], hi = v["hi"],
                 null_mean = mean(nl), null_q05 = quantile(nl, 0.05), null_q95 = quantile(nl, 0.95), frac_null_below = mean(nl < v["r"])) })) }))

  ## (3) raw detection of dosage genes --------------------------------------------------------------------------
  a <- d$group == "A"
  dg <- data.table(case = cs, gene = dos_genes, class = fifelse(dos_genes %in% par, "PAR", "XY_pair"),
                   n_A = sum(a), n_B = sum(!a),
                   det_A = Matrix::rowMeans(raw[, a, drop = FALSE] > 0), det_B = Matrix::rowMeans(raw[, !a, drop = FALSE] > 0),
                   umi10k_A = as.numeric(raw[, a, drop = FALSE] %*% (1e4 / libsize[a])) / sum(a),
                   umi10k_B = as.numeric(raw[, !a, drop = FALSE] %*% (1e4 / libsize[!a])) / sum(!a))
  dg[, ratio := umi10k_A / umi10k_B]
  dg[, p_det := mapply(function(xa, xb) fisher.test(matrix(c(round(xa * sum(a)), sum(a) - round(xa * sum(a)), round(xb * sum(!a)), sum(!a) - round(xb * sum(!a))), 2))$p.value, det_A, det_B)]

  message(sprintf("[%s] cells %d (A %d, B %d); genes kept %d; median depth A %.0f / B %.0f", cs, nrow(d), sum(a), sum(!a), length(keep), median(libsize[a]), median(libsize[!a])))
  message(sprintf("[%s] baseline d: %s", cs, paste(sprintf("%s=%.3f(p%.3f,%dg)", r0$feature, r0$d, r0$p_perm, n_used[r0$feature]), collapse = "  ")))
  message(sprintf("[%s] HLA2 d with F7 as covariate %.3f; F7_all d with HLA2 as covariate %.3f", cs, r2$d, r1[feature == "F7_all", d]))
  message(sprintf("[%s] coupling HLA2~F7_noMP9: %s", cs, paste(sprintf("%s r=%.3f [%.3f, %.3f] null %.3f", cp[pair == "HLA2~F7_noMP9", subset], cp[pair == "HLA2~F7_noMP9", r],
                  cp[pair == "HLA2~F7_noMP9", lo], cp[pair == "HLA2~F7_noMP9", hi], cp[pair == "HLA2~F7_noMP9", null_mean]), collapse = "  ")))
  list(fv = fv, cp = cp, dg = dg)
}

all <- list()
lambo <- CASES[LOY_H2_CASES[match(CASES, case), dataset] == "ZENG_Lambo_2023"]
if (length(lambo)) { stem <- LOY_EXT_STUDIES[dataset == "ZENG_Lambo_2023", file_stem]
  s_full <- readRDS(file.path(LOY_ZENG_RDS_DIR, paste0(stem, "_projected_cleaned.rds"))); message("Lambo loaded: ", ncol(s_full), " cells")
  for (cs in lambo) all[[cs]] <- run_case(cs, s_full); rm(s_full); gc() }
for (cs in setdiff(CASES, lambo)) all[[cs]] <- run_case(cs, NULL)

fv <- rbindlist(lapply(all, `[[`, "fv")); cp <- rbindlist(lapply(all, `[[`, "cp")); dg <- rbindlist(lapply(all, `[[`, "dg"))
# sanity check: the baseline F7_all and MP9 d must reproduce 05's contrast tables (same cells, same genes, same model)
ref <- rbindlist(lapply(CASES, function(cs) fread(file.path(LOY_H2_DIR, sprintf("contrast_%s.csv", cs)))[feature %in% c("F7_TNFA_NFKB", "MP9"), .(case, feature = sub("F7_TNFA_NFKB", "F7_all", feature), d_05 = d)]))
chk <- merge(fv[design == "baseline"], ref, by = c("case", "feature")); print(chk[, .(case, feature, d = round(d, 4), d_05 = round(d_05, 4))])
stopifnot(nrow(chk) == 2 * length(CASES), all(abs(chk$d - chk$d_05) < 1e-6))
fwrite(fv, file.path(LOY_H2_DIR, "posthoc_f7_variants.csv")); fwrite(cp, file.path(LOY_H2_DIR, "posthoc_coupling.csv")); fwrite(dg, file.path(LOY_H2_DIR, "posthoc_dosage_detection.csv"))
message("\n=== (1) d by feature, baseline ==="); print(dcast(fv[design == "baseline"], feature + n_genes ~ case, value.var = "d")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])
message("\n=== (1) permutation p, baseline ==="); print(dcast(fv[design == "baseline"], feature ~ case, value.var = "p_perm"))
message("\n=== (2) mutual adjustment ==="); print(fv[(feature == "HLA2" & design != "plus_HLA2") | (feature %in% c("F7_all", "F7_noMP9", "NFKB_core")), .(case, feature, design, d = round(d, 3), p_perm)][order(case, feature, design)])
message("\n=== (2) per-cell partial correlation with HLA2 ==="); print(cp[, .(case, subset, n, pair, r = round(r, 3), lo = round(lo, 3), hi = round(hi, 3), null_mean = round(null_mean, 3), null_q05 = round(null_q05, 3), null_q95 = round(null_q95, 3))])
message("\n=== (3) dosage genes, raw counts (ratio = -Y / Y+ UMI per 10k) ==="); print(dg[det_A >= 0.01 | det_B >= 0.01, .(case, gene, class, det_A = round(det_A, 3), det_B = round(det_B, 3), ratio = round(ratio, 2), p_det = signif(p_det, 2))][order(case, class, gene)], nrows = 200)
message(sprintf("done (%.1f min)", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
