#!/usr/bin/env Rscript
# input : D1's synthetic engine (scratch_D1/02_sim_geneset.R, sourced up to its scenario block); REAL design covariates
#         only (ws_s5_scores.csv, select = sample, status, hb, gran, sex, subtype_share, compartment). No expression read.
# output: review1/L1/04_sim_k5_stacked_{raw.csv,log}
# what it does: runs the exact 3.3(b) model of the draft (donor fixed effects + cT + cT:AML + cT:hb + cT:gran +
#   cT:sexmale + cT:CD8share + cT:Xshare; TMM; estimateDisp(robust) on the stacked design; camera est. on "cT:AML")
#   on SYNTHETIC stacked counts that carry the real donor pairing and covariate values, plants the K5 effects at the
#   count level (y + Poisson(y (2^b - 1)) in AML samples), refits T's and X's own S5 fits for the ratio, and applies
#   the 3.3(b) reading. Scenario "realistic": half of the module noise shared between compartments, X BCV x 1.3,
#   X libraries 0.2 x T. Variant "all20": all 20 sets of a type planted in one refit (B only).
#   usage: Rscript 04_sim_k5_stacked.R <n_worlds> <n_cores>
src <- readLines("/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D1/02_sim_geneset.R")
eval(parse(text = src[1:(grep("^## ---- scenarios", src) - 1L)]))      # D1 engine: make_world, sim_libs, s5_fit, constants
OUT <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/review1/L1"
a <- commandArgs(TRUE); NW <- if (length(a) >= 1) as.integer(a[1]) else 1L; NCORE <- if (length(a) >= 2) as.integer(a[2]) else 1L

## real design covariates (no outcome column read); synthetic donor i = real T row i, AML first (D1's status order)
D <- fread("/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws/ws_s5_scores.csv",
           select = c("sample", "status", "hb", "gran", "sex", "subtype_share", "compartment"))
DT <- D[compartment == "T"][order(-(status == "AML"), sample)]; stopifnot(nrow(DT) == 25, DT[1:15, all(status == "AML")])
own_design <- function(d) model.matrix(~ status + hb + gran + sex + subtype_share,
  transform(d, status = factor(status, c("healthy", "AML")), sex = factor(sex, c("female", "male"))))
XT <- own_design(DT); stopifnot(colnames(XT)[2] == "statusAML")

P <- list(seed = 70300L, sd_l = 0.10, r_TB = 0.65, delta_D = 2, tau = 1, rho = 0.05, n_eff = 0L,
          eff_types = data.table(type = character(), n = integer(), beta = numeric(), pi = numeric(), where = character()),
          lib_B = 0.25, share_F = 0.5, bcv_X = 1.3, lib_X = 0.20)

gen2 <- function(w, p, comp, Fx) {           # D1's gen_counts with partly compartment-specific noise and depth
  n <- length(w$status); lib <- sim_libs(); bcv <- w$bcv; Fm <- w$Fm; l <- w$lT
  if (comp == "X") { lib <- lib * p$lib_X; bcv <- bcv * p$bcv_X; Fm <- sqrt(p$share_F) * w$Fm + sqrt(1 - p$share_F) * Fx; l <- w$lB }
  eta <- matrix(log(2^w$log2cpm / 1e6), G, n) + outer(l, w$D) + sqrt(p$rho) * bcv * Fm[w$modid, , drop = FALSE]
  cnt <- matrix(rnbinom(G * n, mu = exp(eta) * matrix(lib, G, n, byrow = TRUE), size = 1 / (bcv^2 * (1 - p$rho))), G, n)
  dimnames(cnt) <- list(paste0("g", seq_len(G)), DT$sample); cnt
}
plant <- function(cnt, genes, cols, beta) {   # K2b/K5 count-level plant; beta < 0 = binomial thinning
  s <- cnt[genes, cols, drop = FALSE]
  cnt[genes, cols] <- if (beta > 0) s + matrix(rpois(length(s), s * (2^beta - 1)), nrow(s)) else matrix(rbinom(length(s), s, 2^beta), nrow(s))
  cnt }
cam <- function(y, idx, X, contrast) {
  r <- camera(y, idx, X, contrast = contrast, inter.gene.cor = NA, allow.neg.cor = FALSE, sort = FALSE)
  data.table(set = rownames(r), NGenes = r$NGenes, Direction = r$Direction, PValue = r$PValue, FDR = p.adjust(r$PValue, "BH")) }
shift <- function(tab, g, common) { lf <- tab[common, "logFC"]; names(lf) <- common; mean(lf[intersect(g, common)]) - mean(lf[setdiff(common, g)]) }

one_world <- function(wi, Xnames = c("B", "NK"), n_per = 5L, tag = "realistic") {
  set.seed(P$seed + wi); w <- make_world(P); Fx <- matrix(rnorm(length(w$Fm)), nrow(w$Fm))
  cT0 <- gen2(w, P, "T", Fx); cX0 <- gen2(w, P, "X", Fx)
  modsets <- setNames(lapply(w$modules, function(i) paste0("g", i)), paste0("M", seq_along(w$modules)))
  res <- list(); univ <- list()
  ## common universe before planting (the draft draws plant genes from the common universe)
  fT0 <- s5_fit(cT0, XT, as.integer(DT$status == "AML"))
  for (Xn in Xnames) {
    DX <- D[compartment == Xn]; DX <- DX[match(intersect(DT$sample, DX$sample), DX$sample)]   # paired donors, T order
    XX <- own_design(DX); stX <- as.integer(DX$status == "AML")
    fX0 <- s5_fit(cX0[, DX$sample], XX, stX)
    common0 <- intersect(rownames(fT0$y), rownames(fX0$y))
    ## plant sets: n_per disjoint random 50-gene sets (same sets for the three types; T plants identical across types)
    pg <- matrix(sample(common0, 50L * n_per), 50L)
    psets <- setNames(lapply(seq_len(n_per), function(j) pg[, j]), paste0("P", seq_len(n_per)))
    amlT <- DT$sample[DT$status == "AML"]; amlX <- DX$sample[DX$status == "AML"]
    cT <- plant(cT0, as.vector(pg), amlT, 0.6)
    fT <- s5_fit(cT, XT, as.integer(DT$status == "AML"))
    for (ty in c("T_only", "shared", "half")) {
      bX <- c(T_only = 0, shared = 0.6, half = 0.3)[[ty]]
      cX <- if (bX > 0) plant(cX0, as.vector(pg), amlX, bX) else cX0
      fX <- s5_fit(cX[, DX$sample], XX, stX)
      common <- intersect(rownames(fT$y), rownames(fX$y))
      univ[[length(univ) + 1]] <- data.table(world = wi, X = Xn, type = ty, nT0 = nrow(fT0$y), nT = nrow(fT$y), nX0 = nrow(fX0$y), nX = nrow(fX$y),
        common0 = length(common0), common = length(common), planted_lost = sum(!as.vector(pg) %in% common))
      ## stacked model exactly as draft 3.3(b)
      S <- data.table(sample = rep(DX$sample, 2), cT = rep(c(1, 0), each = nrow(DX)), AML = rep(stX, 2), hb = rep(DX$hb, 2), gran = rep(DX$gran, 2),
                      sexmale = rep(as.integer(DX$sex == "male"), 2), CD8share = rep(DT[match(DX$sample, sample), subtype_share], 2), Xshare = rep(DX$subtype_share, 2))
      Xs <- cbind(model.matrix(~ 0 + factor(sample), S), cT = S$cT, `cT:AML` = S$cT * S$AML, `cT:hb` = S$cT * S$hb, `cT:gran` = S$cT * S$gran,
                  `cT:sexmale` = S$cT * S$sexmale, `cT:CD8share` = S$cT * S$CD8share, `cT:Xshare` = S$cT * S$Xshare)
      stopifnot(qr(Xs)$rank == ncol(Xs))
      cnt <- cbind(cT[common, DX$sample], cX[common, DX$sample]); colnames(cnt) <- paste0(rep(c("T_", "X_"), each = nrow(DX)), DX$sample)
      y <- calcNormFactors(DGEList(cnt), method = "TMM")
      t0 <- proc.time()[["elapsed"]]; y <- estimateDisp(y, Xs, robust = TRUE); tdisp <- proc.time()[["elapsed"]] - t0
      allsets <- c(modsets, psets); idx <- lapply(allsets, intersect, common); idx <- idx[lengths(idx) >= 15]
      cI <- cam(y, idx, Xs, "cT:AML")
      cTT <- cam(fT$y, lapply(allsets, intersect, rownames(fT$y))[names(idx)], XT, "statusAML")   # T primary on the same sets
      setnames(cTT, c("NGenes", "Direction", "PValue", "FDR"), paste0("T_", c("NGenes", "Direction", "PValue", "FDR")))
      r <- merge(cI, cTT, by = "set")
      r[, sT := vapply(set, function(s) shift(fT$tab, allsets[[s]], common), 1)]
      r[, sX := vapply(set, function(s) shift(fX$tab, allsets[[s]], common), 1)]
      r[, n_planted := vapply(set, function(s) sum(allsets[[s]] %in% pg), 1L)]
      r[, `:=`(world = wi, X = Xn, type = ty, tag = tag, disp_s = tdisp, common_n = length(common),
               prior_df = if (is.null(y$prior.df)) NA_real_ else median(y$prior.df), common_bcv = sqrt(y$common.dispersion))]
      res[[length(res) + 1]] <- r
    }
  }
  list(res = rbindlist(res), univ = rbindlist(univ))
}

t0 <- proc.time()
RL <- mclapply(seq_len(NW), function(wi) tryCatch(one_world(wi), error = function(e) e), mc.cores = NCORE, mc.set.seed = FALSE)
## variant: all 20 sets of a type planted in a single refit (B only), one world
RA <- tryCatch(one_world(99L, Xnames = "B", n_per = 20L, tag = "all20"), error = function(e) e)
bad <- vapply(c(RL, list(RA)), inherits, TRUE, "error"); if (any(bad)) { print(c(RL, list(RA))[bad][[1]]); stop("world failed") }
R <- rbindlist(c(lapply(RL, `[[`, "res"), list(RA$res))); U <- rbindlist(c(lapply(RL, `[[`, "univ"), list(RA$univ)))
tt <- proc.time() - t0
fwrite(R, file.path(OUT, "04_sim_k5_stacked_raw.csv")); fwrite(U, file.path(OUT, "04_sim_k5_universe.csv"))

R[, rho := sX / sT]
R[, call := fifelse(FDR < 0.05 & Direction == "Up", "larger_T", fifelse(FDR < 0.05 & Direction == "Down", "larger_X",
            fifelse(rho >= 0.5, "shared", "not_resolved")))]
sink(file.path(OUT, "04_sim_k5_stacked.log"), split = TRUE)
cat(sprintf("worlds %d (+1 all20) | elapsed %.0f s | cpu self %.0f + children %.0f s\n", NW, tt[["elapsed"]], tt[["user.self"]] + tt[["sys.self"]], tt[["user.child"]] + tt[["sys.child"]]))
cat("scenario: GW component (delta_D 2, loadings sd 0.10, T-X loading cor 0.65), rho 0.05, module noise shared 50 %, X BCV x1.3, X libs 0.2 x T\n")
cat("estimateDisp on stacked design: median s", round(median(R$disp_s), 1), "| common BCV median", round(median(R$common_bcv), 3), "| prior.df median", round(median(R$prior_df), 1), "\n")
cat("\nuniverse (genes tested; planted genes lost from the common universe after the refit):\n"); print(U[, .(nT0 = median(nT0), nT = median(nT), nX0 = median(nX0), nX = median(nX), common0 = median(common0), common = median(common), planted_lost = sum(planted_lost)), by = .(tag, X, type)])
cat("\n== K5 reading on planted sets (n per cell; share of calls)\n")
K <- R[n_planted >= 45][, .(n = .N, larger_T = sum(call == "larger_T"), shared = sum(call == "shared"), not_res = sum(call == "not_resolved"),
        larger_X = sum(call == "larger_X"), T_flag = sum(T_FDR < 0.05 & T_Direction == "Up"), rho_med = round(median(rho), 2),
        rho_q10 = round(quantile(rho, .1), 2), rho_q90 = round(quantile(rho, .9), 2), I_p_med = signif(median(PValue), 2)), by = .(tag, X, type)]
print(K)
cat("\n== null module sets (no planted gene): interaction FDR < 0.05 rate, and call shares among T-flagged ones\n")
print(R[n_planted == 0, .(n = .N, I_FDR05 = round(mean(FDR < 0.05), 3), I_p05 = round(mean(PValue < 0.05), 3), T_flagged = sum(T_FDR < 0.05)), by = .(tag, X, type)])
cat("\n== K5 bars on these numbers: (i) larger_T >= 14/20; (ii) larger_T <= 2/20 and shared >= 14/20 (per X, realistic tag)\n")
print(K[tag == "realistic", .(X, type, n, larger_T, shared, pass_i = if (type[1] == "T_only") larger_T >= 0.7 * n else NA,
                               pass_ii = if (type[1] == "shared") larger_T <= 0.1 * n & shared >= 0.7 * n else NA), by = .(X, type)][, -(1:2)])
sink()
