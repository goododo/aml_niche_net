#!/usr/bin/env Rscript
# INPUT : the functions of 02_sim_geneset.R (sourced up to its scenario block); HALLMARK covered sizes only.
# OUTPUT: SCR/03_crosscomp_raw.csv, SCR/03_crosscomp_summary.csv, SCR/03_crosscomp.log
# what it does: D1 question 3. Same synthetic world, now with a second compartment ("B": 13 AML + 10 healthy of the
#   same donors, 1/4 of T's library sizes, genome-wide loadings correlated 0.65 with T's, the same donor factor and
#   module activities). Effect sets are either T-only or shared (same size in T and B). Readouts:
#   (1) camera per compartment; naive rule "T-specific = FDR < 0.05 in T and p >= 0.05 in B";
#   (2) one model on stacked T + B pseudobulks, donor as a fixed effect, compartment x status interaction,
#       camera/fry on the interaction coefficient (T effect minus B effect);
#   (3) descriptive set-level ratio: mean B logFC / mean T logFC over the set's genes tested in both.
#   usage: Rscript 03_sim_crosscomp.R <n_rep> <n_cores>
src <- readLines("/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D1/02_sim_geneset.R")
eval(parse(text = src[1:(grep("^## ---- scenarios", src) - 1L)]))
a <- commandArgs(TRUE); NREP <- if (length(a) >= 1) as.integer(a[1]) else 2L; NCORE <- if (length(a) >= 2) as.integer(a[2]) else 1L

EFF3 <- data.table(type = c("T_only", "shared"), n = c(5L, 5L), beta = c(0.30, 0.30), pi = c(1, 1), where = c("T_only", "both"))
P3 <- list(seed = 5000L, sd_l = 0.10, r_TB = 0.65, delta_D = 2, tau = 1, rho = 0.05, n_eff = 10L, eff_types = EFF3, lib_B = 0.25)
B_KEEP <- c(1:13, 16:25)          # B: 13 AML + 10 healthy (S5 eligibility), all among the T donors

cam2 <- function(y, X, coef, idx) {
  Z <- edgeR:::.zscoreDGE(y = y, design = X, contrast = coef)
  c1 <- camera(Z, idx, X, contrast = coef, inter.gene.cor = 0.01, sort = FALSE)
  c2 <- camera(Z, idx, X, contrast = coef, inter.gene.cor = NA, sort = FALSE)
  fr <- fry(Z, idx, X, contrast = coef, standardize = "none", sort = "none")
  sg <- function(d) ifelse(d == "Up", 1, -1)
  data.table(set = names(idx), m = lengths(idx), cam01_p = c1$PValue, cam01_d = sg(c1$Direction),
             camest_p = c2$PValue, camest_d = sg(c2$Direction), fry_p = fr$PValue, fry_d = sg(fr$Direction))
}

one_rep_X <- function(rep) {
  p <- P3; set.seed(p$seed + rep)
  w <- make_world(p)
  cT <- gen_counts(w, p, "T"); cB <- gen_counts(w, p, "B", keep_samples = B_KEEP)
  XT <- s5_design(w$status); XB <- s5_design(w$status[B_KEEP])
  fT <- s5_fit(cT, XT, w$status); fB <- s5_fit(cB, XB, w$status[B_KEEP])
  sets <- setNames(lapply(w$modules, function(i) paste0("g", i)), paste0("S", seq_along(w$modules)))
  both <- intersect(rownames(fT$y), rownames(fB$y))
  idx_of <- function(genes) { ix <- lapply(sets, function(g) which(genes %in% intersect(g, both))); ix[lengths(ix) >= 15] }
  ## (1) per compartment, sets restricted to genes tested in both (so the three readouts see the same genes)
  rT <- cam2(fT$y, XT, 2, idx_of(rownames(fT$y))); rB <- cam2(fB$y, XB, 2, idx_of(rownames(fB$y)))
  ## (2) stacked model with donor fixed effects; interaction = T effect minus B effect
  cnt <- cbind(cT[both, ], cB[both, ])
  donor <- factor(c(seq_len(25), B_KEEP)); compT <- rep(c(1, 0), c(25, length(B_KEEP)))
  aml <- c(w$status, w$status[B_KEEP])
  X <- cbind(model.matrix(~ 0 + donor), compT = compT, compT_AML = compT * aml)
  stopifnot(qr(X)$rank == ncol(X))
  y <- DGEList(cnt); y <- calcNormFactors(y, method = "TMM")
  # simplification for the CPU budget: gene-wise dispersion = mean of the T-only and B-only tagwise estimates
  # (the truth has the same BCV in both compartments); a real analysis would run estimateDisp on the stacked design
  y$tagwise.dispersion <- (fT$y$tagwise.dispersion[match(both, rownames(fT$y))] + fB$y$tagwise.dispersion[match(both, rownames(fB$y))]) / 2
  rI <- cam2(y, X, ncol(X), idx_of(rownames(y)))
  ## (3) descriptive ratio of set-mean logFC (B / T) over genes in both
  lT <- fT$tab[both, "logFC"]; lB <- fB$tab[both, "logFC"]; names(lT) <- names(lB) <- both
  ratio <- vapply(names(idx_of(both)), function(s) { g <- intersect(sets[[s]], both); mean(lB[g]) / mean(lT[g]) }, 1)
  setnames(rT, c("cam01_p", "cam01_d", "camest_p", "camest_d", "fry_p", "fry_d"), paste0("T_", c("cam01_p", "cam01_d", "camest_p", "camest_d", "fry_p", "fry_d")))
  setnames(rB, c("cam01_p", "cam01_d", "camest_p", "camest_d", "fry_p", "fry_d"), paste0("B_", c("cam01_p", "cam01_d", "camest_p", "camest_d", "fry_p", "fry_d")))
  setnames(rI, c("cam01_p", "cam01_d", "camest_p", "camest_d", "fry_p", "fry_d"), paste0("I_", c("cam01_p", "cam01_d", "camest_p", "camest_d", "fry_p", "fry_d")))
  r <- Reduce(function(a, b) merge(a, b[, !"m"], by = "set"), list(rT, rB, rI))
  r[, ratio_BT := ratio[set]]
  r[, k := as.integer(sub("^S", "", set))]
  r <- merge(r, w$ET[, .(k = set, type, true_dir = dir)], by = "k", all.x = TRUE); r[is.na(type), type := "null"]
  r[, `:=`(rep = rep, n_T = nrow(fT$y), n_B = nrow(fB$y), n_both = length(both))]
  r
}

t0 <- proc.time()
RES <- mclapply(seq_len(NREP), one_rep_X, mc.cores = NCORE, mc.set.seed = FALSE)
bad <- vapply(RES, inherits, TRUE, "try-error"); if (any(bad)) { print(RES[bad][[1]]); stop("replicate failed") }
R <- rbindlist(RES); tt <- proc.time() - t0
fwrite(R, file.path(SCR, "03_crosscomp_raw.csv"))
for (pre in c("T", "B", "I")) for (mm in c("cam01", "camest", "fry")) { pc <- paste0(pre, "_", mm, "_p"); R[, (paste0(pre, "_", mm, "_fdr")) := p.adjust(.SD[[1]], "BH"), by = rep, .SDcols = pc] }
dirT <- function(d) is.na(R$true_dir) | d == R$true_dir
S <- R[, .(n = .N,
  T_camest_FDR05 = mean(T_camest_fdr < 0.05), T_cam01_FDR05 = mean(T_cam01_fdr < 0.05),
  B_camest_FDR05 = mean(B_camest_fdr < 0.05), B_cam01_FDR05 = mean(B_cam01_fdr < 0.05),
  naive_Tspecific_camest = mean(T_camest_fdr < 0.05 & B_camest_p >= 0.05),
  naive_Tspecific_cam01 = mean(T_cam01_fdr < 0.05 & B_cam01_p >= 0.05),
  I_camest_FDR05 = mean(I_camest_fdr < 0.05), I_cam01_FDR05 = mean(I_cam01_fdr < 0.05), I_fry_FDR05 = mean(I_fry_fdr < 0.05),
  I_camest_p05 = mean(I_camest_p < 0.05), I_cam01_p05 = mean(I_cam01_p < 0.05),
  ratio_median = median(ratio_BT), ratio_q10 = quantile(ratio_BT, 0.1), ratio_q90 = quantile(ratio_BT, 0.9)), by = type]
## interaction readouts in the direction of the true T effect
S2 <- R[type != "null", .(I_camest_FDR05_rightdir = mean(I_camest_fdr < 0.05 & I_camest_d == true_dir),
                          I_cam01_FDR05_rightdir = mean(I_cam01_fdr < 0.05 & I_cam01_d == true_dir)), by = type]
fwrite(merge(S, S2, by = "type", all.x = TRUE), file.path(SCR, "03_crosscomp_summary.csv"))
sink(file.path(SCR, "03_crosscomp.log"))
cat(sprintf("reps %d | cores %d | elapsed %.1f s | cpu self %.1f + children %.1f s\n", NREP, NCORE, tt[["elapsed"]],
            tt[["user.self"]] + tt[["sys.self"]], tt[["user.child"]] + tt[["sys.child"]]))
cat(sprintf("genes tested: T median %d, B median %d, both median %d\n", as.integer(median(R$n_T)), as.integer(median(R$n_B)), as.integer(median(R$n_both))))
print(merge(S, S2, by = "type", all.x = TRUE))
sink(); cat(readLines(file.path(SCR, "03_crosscomp.log")), sep = "\n")
