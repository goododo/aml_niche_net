#!/usr/bin/env Rscript
# INPUT : nothing from the project except the 45 HALLMARK covered set SIZES (SCR/01_hallmark_coverage.csv, column
#         inT_noAmb >= 15; sizes only, no gene identities, no expression). Everything else is synthetic.
# OUTPUT: SCR/02_sim_<tag>_raw.csv (one row per set x replicate x method), SCR/02_sim_<tag>_summary.csv, SCR/02_sim_<tag>.log
# what it does: D1 question 2. Synthetic negative-binomial pseudobulk, 15 vs 10, ~7,000 genes in co-expression modules,
#   an optional genome-wide donor-level component that is shifted in AML, and set-level effects of known size.
#   Fits the S5 engine (edgeR QL, robust; filterByExpr + CPM>1 in half; TMM; 6-column design) and runs
#   camera (inter.gene.cor 0.01), camera (estimated), fry (directional + mixed) and preranked fgsea on the signed QL z.
#   usage: Rscript 02_sim_geneset.R <tag> <n_rep> <n_cores>      (tag selects the scenario block at the bottom)
suppressPackageStartupMessages({ library(edgeR); library(limma); library(fgsea); library(data.table); library(parallel) })
Sys.setenv(OMP_PROC_BIND = "false"); Sys.setenv(OMP_NUM_THREADS = "1")
SCR <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D1"
a <- commandArgs(TRUE); TAG <- if (length(a) >= 1) a[1] else "smoke"; NREP <- if (length(a) >= 2) as.integer(a[2]) else 2L
NCORE <- if (length(a) >= 3) as.integer(a[3]) else 1L
SIZES <- fread(file.path(SCR, "01_hallmark_coverage.csv"))[inT_noAmb >= 15, inT_noAmb]   # 45 sizes, 15-181
stopifnot(length(SIZES) == 45L)

## ---- fixed simulation constants (stated in D1 section 2) ----------------------------------------------------
G      <- 7000L
N_AML  <- 15L; N_H <- 10L
BCV_MED <- 0.35       # edgeR User's Guide: typical common BCV 0.4 for human data (SCR/edgeRUsersGuide.txt:1186-1188)
BCV_SDLOG <- 0.30     # gene-to-gene spread of BCV (lognormal)
BG_MODULE <- 50L      # background genes sit in modules of 50 with the same within-module correlation

sim_libs <- function() {   # kept-T-cell numbers: AML median 336 (20-1,208), healthy 617 (123-1,917) (verify_G07.json item 8)
  ca <- pmin(pmax(round(rlnorm(N_AML, log(336), 0.8)), 20), 1208); ch <- pmin(pmax(round(rlnorm(N_H, log(617), 0.6)), 123), 1917)
  c(ca * 3000, ch * 2200)  # UMI per kept T cell ~3,000 AML / ~2,200 healthy (A1-G01 5.6 range 1,684-4,196)
}

## one synthetic world (genes, modules, loadings); compartments share donors and module factors
make_world <- function(p) {
  modules <- list(); start <- 1L
  for (k in seq_along(SIZES)) { modules[[k]] <- start:(start + SIZES[k] - 1L); start <- start + SIZES[k] }
  rest <- start:G; nb <- ceiling(length(rest) / BG_MODULE)
  bg <- split(rest, rep(seq_len(nb), each = BG_MODULE, length.out = length(rest)))
  allmod <- c(modules, bg)
  modid <- integer(G); for (k in seq_along(allmod)) modid[allmod[[k]]] <- k
  log2cpm <- 3 + rgamma(G, shape = 2, scale = 1.3)   # tested genes: filterByExpr at ~1M-UMI libraries keeps log2CPM >~ 3
  bcv <- BCV_MED * exp(rnorm(G, 0, BCV_SDLOG))
  lT <- rnorm(G, 0, p$sd_l)                                                 # genome-wide loadings (natural log per unit D)
  lB <- p$r_TB * lT + sqrt(1 - p$r_TB^2) * rnorm(G, 0, p$sd_l)               # second compartment, correlated loadings
  status <- rep(c(1, 0), c(N_AML, N_H))
  D <- p$delta_D * status + p$tau * rnorm(N_AML + N_H)                       # donor-level factor, shifted in AML
  Fm <- matrix(rnorm(length(allmod) * (N_AML + N_H)), length(allmod))       # module activity per donor (shared by compartments)
  ## set effects (log2), chosen sets and types
  eff_sets <- sample(seq_along(SIZES), p$n_eff)
  types <- rep(p$eff_types$type, p$eff_types$n)[seq_len(p$n_eff)]
  ET <- data.table(set = eff_sets, type = types)
  ET <- merge(ET, p$eff_types[, .(type, beta, pi, where)], by = "type")
  ET[, dir := sample(c(-1, 1), .N, replace = TRUE)]
  list(modules = modules, allmod = allmod, modid = modid, log2cpm = log2cpm, bcv = bcv, lT = lT, lB = lB,
       status = status, D = D, Fm = Fm, ET = ET)
}

gen_counts <- function(w, p, comp = c("T", "B"), keep_samples = NULL) {
  comp <- match.arg(comp)
  n <- length(w$status); lib <- sim_libs(); if (comp == "B") lib <- lib * p$lib_B
  rho <- p$rho
  base <- log(2^w$log2cpm / 1e6)
  l <- if (comp == "T") w$lT else w$lB
  eta <- matrix(base, G, n) + outer(l, w$D) + sqrt(rho) * w$bcv * w$Fm[w$modid, , drop = FALSE]
  for (k in seq_len(nrow(w$ET))) { e <- w$ET[k]
    if (comp == "B" && e$where == "T_only") next
    g <- w$modules[[e$set]]; if (e$pi < 1) g <- g[seq_len(max(1, round(e$pi * length(g))))]
    eta[g, w$status == 1] <- eta[g, w$status == 1] + e$dir * e$beta * log(2) }
  mu <- exp(eta) * matrix(lib, G, n, byrow = TRUE)
  size <- 1 / (w$bcv^2 * (1 - rho))
  cnt <- matrix(rnbinom(G * n, mu = mu, size = size), G, n)
  rownames(cnt) <- paste0("g", seq_len(G)); colnames(cnt) <- paste0(comp, "_s", seq_len(n))
  if (!is.null(keep_samples)) cnt <- cnt[, keep_samples, drop = FALSE]
  cnt
}

s5_design <- function(status) {   # 6 columns as in S5: status + hb + gran + sex + subtype_share (covariates without effect)
  n <- length(status)
  d <- data.frame(status = factor(ifelse(status == 1, "AML", "healthy"), c("healthy", "AML")),
                  hb = rnorm(n), gran = rnorm(n), sex = factor(sample(c("F", "M"), n, replace = TRUE), c("F", "M")),
                  share = runif(n, 0.2, 0.8))
  if (length(unique(d$sex)) < 2) d$sex[1:2] <- c("F", "M")
  model.matrix(~ status + hb + gran + sex + share, d)
}

s5_fit <- function(cnt, X, status) {
  grp <- factor(status)
  y <- DGEList(cnt, group = grp)
  keep <- filterByExpr(y, group = grp) & (rowMeans(cpm(y) > 1) >= 0.5)
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, X, robust = TRUE); f <- glmQLFit(y, X, robust = TRUE); t <- glmQLFTest(f, coef = 2)
  tt <- t$table; z <- sign(tt$logFC) * qnorm(pmax(tt$PValue, 1e-300) / 2, lower.tail = FALSE); names(z) <- rownames(tt)
  list(y = y, tab = tt, z = z)
}

set_tests <- function(fit, X, sets) {
  y <- fit$y; idx <- lapply(sets, function(g) which(rownames(y) %in% g))
  names(idx) <- names(sets); idx <- idx[lengths(idx) >= 15]
  Z <- edgeR:::.zscoreDGE(y = y, design = X, contrast = 2)               # what edgeR's camera/fry.DGEList do internally
  c1 <- camera(Z, idx, X, contrast = 2, inter.gene.cor = 0.01, sort = FALSE)
  c2 <- camera(Z, idx, X, contrast = 2, inter.gene.cor = NA, sort = FALSE)
  fr <- fry(Z, idx, X, contrast = 2, standardize = "none", sort = "none")
  zs <- fit$z; pw <- lapply(idx, function(i) rownames(y)[i])
  fg <- suppressWarnings(fgsea(pw, zs, minSize = 15, maxSize = 500, eps = 0, BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)))
  fg <- fg[match(names(idx), fg$pathway)]
  cpr <- cameraPR(zs, idx, inter.gene.cor = 0.01, sort = FALSE)        # preranked camera on the same signed z (reference)
  sgn <- function(d) ifelse(d == "Up", 1, -1)
  rbind(
    data.table(set = names(idx), m = lengths(idx), method = "camera_0.01", p = c1$PValue, dirn = sgn(c1$Direction), cor_est = NA_real_),
    data.table(set = names(idx), m = lengths(idx), method = "camera_est", p = c2$PValue, dirn = sgn(c2$Direction), cor_est = c2$Correlation),
    data.table(set = names(idx), m = lengths(idx), method = "fry_dir", p = fr$PValue, dirn = sgn(fr$Direction), cor_est = NA_real_),
    data.table(set = names(idx), m = lengths(idx), method = "fry_mixed", p = fr$PValue.Mixed, dirn = NA_real_, cor_est = NA_real_),
    data.table(set = names(idx), m = lengths(idx), method = "fgsea_pre", p = fg$pval, dirn = sign(fg$NES), cor_est = NA_real_),
    data.table(set = names(idx), m = lengths(idx), method = "cameraPR_0.01", p = cpr$PValue, dirn = sgn(cpr$Direction), cor_est = NA_real_))
}

one_rep_T <- function(rep, p) {
  set.seed(p$seed + rep)
  w <- make_world(p)
  cnt <- gen_counts(w, p, "T")
  X <- s5_design(w$status)
  fit <- s5_fit(cnt, X, w$status)
  sets <- setNames(lapply(w$modules, function(i) paste0("g", i)), paste0("S", seq_along(w$modules)))
  r <- set_tests(fit, X, sets)
  r[, k := as.integer(sub("^S", "", set))]
  r <- merge(r, w$ET[, .(k = set, type, beta, pi, where, true_dir = dir)], by = "k", all.x = TRUE)
  r[is.na(type), type := "null"]
  r[, `:=`(rep = rep, n_genes_tested = nrow(fit$tab), n_s5_hits = sum(p.adjust(fit$tab$PValue, "BH") < 0.10 & abs(fit$tab$logFC) > 0.5))]
  ## true genome-wide set mean (log2) from the loadings, for reading the competitive vs self-contained null
  dd <- p$delta_D; gwm <- vapply(seq_along(w$modules), function(kk) mean(w$lT[w$modules[[kk]]]) * dd / log(2), 1)
  r[, set_mean_true_gw := gwm[k]]
  r
}

summarise <- function(R) {
  R[, fdr := p.adjust(p, "BH"), by = .(rep, method)]
  R[, hit05 := p < 0.05]; R[, hitFDR := fdr < 0.05]
  R[, dir_ok := is.na(dirn) | is.na(true_dir) | dirn == true_dir]
  S <- R[, .(n_tests = .N, rate_p05 = mean(hit05), rate_FDR05 = mean(hitFDR),
             rate_FDR05_rightdir = mean(hitFDR & dir_ok)), by = .(method, type)]
  S[order(type, method)]
}

## ---- scenarios -------------------------------------------------------------------------------------------------
EFF <- data.table(type = c("E1_b0.15_all", "E2_b0.30_all", "E3_b0.50_30pct"), n = c(4L, 3L, 3L),
                  beta = c(0.15, 0.30, 0.50), pi = c(1, 1, 0.3), where = "T_only")
base_p <- list(seed = 1000L, sd_l = 0, r_TB = 0.65, delta_D = 0, tau = 0, rho = 0.05, n_eff = 10L, eff_types = EFF, lib_B = 0.25)
SCEN <- list(
  A_noGW_rho05  = modifyList(base_p, list(seed = 1000L)),
  B_GW_rho05    = modifyList(base_p, list(seed = 2000L, sd_l = 0.10, delta_D = 2, tau = 1)),
  C_GW_rho20    = modifyList(base_p, list(seed = 3000L, sd_l = 0.10, delta_D = 2, tau = 1, rho = 0.20)),
  D_noGW_rho20  = modifyList(base_p, list(seed = 4000L, rho = 0.20)),
  smoke         = modifyList(base_p, list(seed = 9000L, sd_l = 0.10, delta_D = 2, tau = 1)))
if (!TAG %in% names(SCEN)) stop("unknown tag ", TAG)
p <- SCEN[[TAG]]
t0 <- proc.time()
RES <- mclapply(seq_len(NREP), one_rep_T, p = p, mc.cores = NCORE, mc.set.seed = FALSE)
bad <- vapply(RES, inherits, TRUE, "try-error"); if (any(bad)) { print(RES[bad][[1]]); stop("replicate failed") }
R <- rbindlist(RES)
tt <- proc.time() - t0
fwrite(R, file.path(SCR, sprintf("02_sim_%s_raw.csv", TAG)))
S <- summarise(R); S[, scenario := TAG]
fwrite(S, file.path(SCR, sprintf("02_sim_%s_summary.csv", TAG)))
sink(file.path(SCR, sprintf("02_sim_%s.log", TAG)))
cat(sprintf("scenario %s | reps %d | cores %d | elapsed %.1f s | cpu self %.1f s + children %.1f s\n", TAG, NREP, NCORE, tt[["elapsed"]],
            tt[["user.self"]] + tt[["sys.self"]], tt[["user.child"]] + tt[["sys.child"]]))
cat(sprintf("genes tested per rep: median %d (range %d-%d); S5-rule hits per rep: median %d (range %d-%d)\n",
            as.integer(median(R[, n_genes_tested[1], by = rep]$V1)), min(R$n_genes_tested), max(R$n_genes_tested),
            as.integer(median(R[, n_s5_hits[1], by = rep]$V1)), min(R$n_s5_hits), max(R$n_s5_hits)))
cat("camera_est estimated inter-gene correlation over all sets: median", round(median(R[method == "camera_est", cor_est]), 3),
    "IQR", paste(round(quantile(R[method == "camera_est", cor_est], c(.25, .75)), 3), collapse = "-"), "\n")
cat("null sets: SD over sets of the true genome-wide set mean (log2):", round(sd(R[method == "camera_0.01" & type == "null", set_mean_true_gw]), 4), "\n")
print(S)
sink()
print(S); cat(readLines(file.path(SCR, sprintf("02_sim_%s.log", TAG)))[1:4], sep = "\n")
