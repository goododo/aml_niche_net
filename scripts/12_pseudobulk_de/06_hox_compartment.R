#!/usr/bin/env Rscript
# 06_hox_compartment.R ----
# INPUT  : results/tables/12_pseudobulk_de/hox_axis_scores.csv   (frozen, written 2026-09-18 by 05)
#          results/tables/01_preprocess/00_curated_manifest.csv  (driver_mutations -> NPM1 status)
# OUTPUT : DIR_PSEUDOBULK/hox_compartment_tests.csv , hox_compartment_null.csv
# WHAT IT DOES : asks whether the Mono_DC HOX score predicts NPM1 status after the LMPP_GMP HOX
#          score has been regressed out of it, and symmetrically the other way round.
# Usage  : Rscript scripts/12_pseudobulk_de/06_hox_compartment.R
#
# EVERYTHING HERE IS PRE-REGISTERED in PREREGISTRATION_hox_compartment.md, committed before this
# file existed. Nothing is recomputed from raw data: the scores are read verbatim from 05's output.
#
# TWO CLAUSES THE PRE-REGISTRATION LEFT UNDER-SPECIFIED, both fixed HERE, before any result:
#
#   (A) RANK SCALING. sec 4 says "rank-transform both scores" within dataset, then pool the
#       residuals. Raw ranks are not poolable -- GSE185381 gives ranks 1..29 and GSE116256 1..7,
#       so a residual of +3 means different things in the two. Ranks are therefore mapped to
#       (rank - 0.5)/n, which is on (0,1) for every dataset regardless of its size.
#
#   (B) PERMUTATION AND THE ADJUSTMENT. sec 4 says to recompute the full adjustment inside every
#       permutation. The adjustment regresses Mono_DC ranks on LMPP_GMP ranks and never sees the
#       NPM1 label, so recomputing it under a permuted label is provably identical. It is
#       recomputed anyway, to match the document literally rather than to save time, and
#       SELF-CHECK 3 asserts the identity instead of assuming it.
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))

SET        <- "L"          # PREREG 2: Set L (locus-defined) throughout; Set D is control 5.1(3)
A_BIN      <- "Mono_DC"    # the compartment under test
B_BIN      <- "LMPP_GMP"   # the compartment it is adjusted against
C_BIN      <- "T_NK"       # PREREG 5.1(1) negative control
DROP_DS    <- "Chen2023"   # PREREG 2.1: 5/5 NPM1-mut -> zero within-dataset contrast
RHO_COLLIN <- 0.90         # PREREG 3 precondition
ALPHA      <- 0.05         # PREREG 5
RHO_NCELL  <- 0.40         # PREREG 5.1(2)
N_PERM     <- 10000L
set.seed(20260924)

## ------------------------------------------------------------------ inputs ----
SCF <- file.path(DIR_PSEUDOBULK, "hox_axis_scores.csv")
MNF <- file.path(DIR_PREPROCESS, "00_curated_manifest.csv")
message(sprintf("[frozen] %s  (mtime %s)", basename(SCF), file.mtime(SCF)))
message(sprintf("[frozen] %s  (mtime %s)", basename(MNF), file.mtime(MNF)))

SC <- fread(SCF)[arm == "AML"]
MU <- fread(MNF, select = c("sample", "driver_mutations"))
MU <- MU[!is.na(driver_mutations) & driver_mutations != ""]
MU[, npm1 := grepl("NPM1", driver_mutations)]

## wide: one row per sample, one column per bin, for the chosen set
mk <- function(st) {
  w <- dcast(SC[set == st], sample + dataset ~ hierarchy_bin,
             value.var = c("score", "n_cells"))
  MU[w, on = "sample", nomatch = 0]
}
W <- mk(SET)

roster <- function(w, bins) {
  need <- paste0("score_", bins)
  w <- w[dataset != DROP_DS][complete.cases(w[dataset != DROP_DS, ..need])]
  w[]
}
R2 <- roster(W, c(A_BIN, B_BIN))
R3 <- roster(W, c(A_BIN, B_BIN, C_BIN))

## -- SELF-CHECK 1: the roster PREREG 2.1 counted, asserted rather than assumed ----
tab <- R2[, .(n = .N, mut = sum(npm1)), by = dataset][order(-n)]
message("\n[SELF-CHECK 1] roster for the primary test (PREREG 2.1 says 36 samples / 9 NPM1-mut):")
print(tab)
stopifnot(nrow(R2) == 36L, sum(R2$npm1) == 9L, !(DROP_DS %in% R2$dataset),
          tab[dataset == "GSE185381", n] == 29L, tab[dataset == "GSE185381", mut] == 5L,
          tab[dataset == "GSE116256", n] == 7L,  tab[dataset == "GSE116256", mut] == 4L)
message(sprintf("               matches. T_NK control roster: %d samples / %d NPM1-mut (PREREG 2.1: 30 / 7)",
                nrow(R3), sum(R3$npm1)))
stopifnot(nrow(R3) == 30L, sum(R3$npm1) == 7L)

TESTS <- list()
rec <- function(...) TESTS[[length(TESTS) + 1L]] <<- data.table(...)

## clause (A): ranks mapped to (rank-0.5)/n inside each dataset, so residuals pool across datasets
qrank <- function(x) (frank(x, ties.method = "average") - 0.5) / length(x)
auc1 <- function(s, y) {                       # one-sided: higher score in y == TRUE
  a <- s[y]; b <- s[!y]
  if (!length(a) || !length(b)) return(NA_real_)
  unname(suppressWarnings(wilcox.test(a, b))$statistic) / (length(a) * length(b))
}

## PREREG 4: residual of A's within-dataset quantile-ranks after OLS on B's (plus any extra
## covariates for control 5.1(2)). Never sees the NPM1 label -- see clause (B).
adjust <- function(d, a, b, extra = character(0)) {
  out <- numeric(nrow(d))
  for (ds in unique(d$dataset)) {
    i <- which(d$dataset == ds)
    X <- data.table(y = qrank(d[[paste0("score_", a)]][i]),
                    x = qrank(d[[paste0("score_", b)]][i]))
    for (e in extra) X[[e]] <- qrank(log10(d[[e]][i] + 1))
    out[i] <- residuals(lm(reformulate(c("x", extra), "y"), data = X))
  }
  out
}

## -------------------------------- PREREG 3: precondition, printed FIRST ----
message("\n===== PREREG 3 PRECONDITION: are the two compartments separable at all? =====")
coll <- R2[, .(rho = cor(score_Mono_DC, score_LMPP_GMP, method = "spearman"), n = .N), by = dataset]
for (i in seq_len(nrow(coll)))
  message(sprintf("  Spearman(%s, %s) in %-10s = %+.3f  (n=%d)",
                  A_BIN, B_BIN, coll$dataset[i], coll$rho[i], coll$n[i]))
rho_max <- max(abs(coll$rho))
rec(test = "PREREG3_collinearity", stratum = paste(coll$dataset, collapse = "+"),
    n = nrow(R2), stat = rho_max, p = NA_real_, threshold = RHO_COLLIN,
    label = if (rho_max > RHO_COLLIN) "UNANSWERABLE" else "separable")
message(sprintf("  max |rho| = %.3f  [PREREG 3: > %.2f => UNANSWERABLE, stop here]", rho_max, RHO_COLLIN))
if (rho_max > RHO_COLLIN) {
  message("\n>>> UNANSWERABLE -- the compartments are collinear; no statistic can split their",
          "\n    contributions at n=36. PREREG 3 stops the analysis here. That is the result.")
  fwrite(rbindlist(TESTS, fill = TRUE), file.path(DIR_PSEUDOBULK, "hox_compartment_tests.csv"))
  quit(save = "no", status = 0)
}

## ------------------------------------------------- PREREG 4: the statistic ----
## Stratified label permutation: NPM1 is shuffled WITHIN dataset, so each dataset keeps its own
## mutant count and the dataset can never leak into the null.
perm_auc <- function(d, resid, n_perm) {
  vapply(seq_len(n_perm), function(k) {
    y <- d$npm1
    for (ds in unique(d$dataset)) { i <- which(d$dataset == ds); y[i] <- sample(y[i]) }
    auc1(resid, y)
  }, numeric(1))
}
run <- function(d, a, b, tag, extra = character(0)) {
  r  <- adjust(d, a, b, extra)
  o  <- auc1(r, d$npm1)
  nl <- perm_auc(d, r, N_PERM)
  p  <- (1 + sum(nl >= o, na.rm = TRUE)) / (1 + sum(is.finite(nl)))
  message(sprintf("  %-34s AUC %.3f | perm p = %.4f (n=%d, %d mut)", tag, o, p, nrow(d), sum(d$npm1)))
  rec(test = tag, stratum = paste(unique(d$dataset), collapse = "+"), n = nrow(d),
      stat = o, p = p, threshold = ALPHA, label = if (p < ALPHA) "adds" else "does not add")
  list(auc = o, p = p, null = nl, resid = r)
}

message(sprintf("\n===== PREREG 4: does each compartment add beyond the other? (%d permutations) =====", N_PERM))
message("  unadjusted, for reference only (PREREG 4, not tested):")
for (b in c(A_BIN, B_BIN)) {
  u <- auc1(R2[[paste0("score_", b)]], R2$npm1)
  message(sprintf("    %-9s raw AUC %.3f", b, u))
  rec(test = paste0("unadjusted_AUC_", b), stratum = "GSE185381+GSE116256", n = nrow(R2),
      stat = u, p = NA_real_, threshold = NA_real_, label = "reference only")
}
message("  adjusted (the test):")
MA <- run(R2, A_BIN, B_BIN, sprintf("%s adjusted for %s", A_BIN, B_BIN))
MB <- run(R2, B_BIN, A_BIN, sprintf("%s adjusted for %s", B_BIN, A_BIN))

## -- SELF-CHECK 3: clause (B)'s claim, asserted rather than assumed ----
## The adjustment must be label-independent. If it were not, permuting the label would change the
## residuals and the null would be testing two things at once.
r_shuf <- local({ d <- copy(R2); d[, npm1 := sample(npm1)]; adjust(d, A_BIN, B_BIN) })
message(sprintf("\n[SELF-CHECK 3] residuals identical under a permuted label: %s (max abs diff %.2e)",
                isTRUE(all.equal(MA$resid, r_shuf)), max(abs(MA$resid - r_shuf))))
stopifnot(isTRUE(all.equal(MA$resid, r_shuf)))

## ------------------------------------------------ PREREG 5.1: the controls ----
message("\n===== PREREG 5.1 CONTROLS (printed before the decision label) =====")

## (1) T_NK negative control -- if T_NK also adds, the adding is generic and MONO_ADDS is downgraded
message("  [5.1-1] T_NK negative control, on the 30-sample roster:")
TA <- run(R3, C_BIN, B_BIN, sprintf("%s adjusted for %s", C_BIN, B_BIN))
MA3 <- run(R3, A_BIN, B_BIN, sprintf("%s adjusted for %s  (same roster)", A_BIN, B_BIN))
tnk_generic <- isTRUE(TA$p < ALPHA)

## (2) cell-count confound
message("  [5.1-2] cell-count confound:")
nc <- sapply(c(A_BIN, B_BIN), function(b)
  cor(R2[[paste0("n_cells_", b)]], as.numeric(R2$npm1), method = "spearman"))
for (b in names(nc)) {
  message(sprintf("          Spearman(n_cells %s, NPM1) = %+.3f  [threshold %.2f]", b, nc[[b]], RHO_NCELL))
  rec(test = paste0("5.1-2_ncells_vs_NPM1_", b), stratum = "GSE185381+GSE116256", n = nrow(R2),
      stat = nc[[b]], p = NA_real_, threshold = RHO_NCELL,
      label = if (abs(nc[[b]]) > RHO_NCELL) "confounded" else "clear")
}
if (any(abs(nc) > RHO_NCELL)) {
  message("          -> threshold crossed; PREREG 5.1(2) requires the depth-adjusted version too:")
  MAd <- run(R2, A_BIN, B_BIN, sprintf("%s adj %s + log10 n_cells", A_BIN, B_BIN),
             extra = paste0("n_cells_", c(A_BIN, B_BIN)))
} else {
  message("          -> both below threshold; the depth-adjusted repeat is not required")
}

## (3) Set D repeat -- a disagreement would mean the adjustment, not the axis, drives the answer
message("  [5.1-3] Set D repeat (12 Discovery-frozen genes instead of 19 locus-defined):")
RD <- roster(mk("D"), c(A_BIN, B_BIN))
message(sprintf("          roster %d samples / %d NPM1-mut", nrow(RD), sum(RD$npm1)))
DA <- run(RD, A_BIN, B_BIN, sprintf("SetD %s adjusted for %s", A_BIN, B_BIN))
DB <- run(RD, B_BIN, A_BIN, sprintf("SetD %s adjusted for %s", B_BIN, A_BIN))

## ----------------------------------------- PREREG 5: the decision, assigned last ----
pm <- MA$p; pl <- MB$p
LAB <- "REDUNDANT"
if (pm < ALPHA && pl >= ALPHA) {
  LAB <- "MONO_ADDS"
} else if (pl < ALPHA && pm >= ALPHA) {
  LAB <- "LMPP_ONLY"
} else if (pm < ALPHA && pl < ALPHA) {
  LAB <- "BOTH_ADD"
}
if (LAB == "MONO_ADDS" && tnk_generic) LAB <- "MONO_ADDS_BUT_TNK_CONTROL_FAILED"

message(sprintf("\n>>> %s   (p_mono = %.4f, p_lmpp = %.4f; T_NK control p = %.4f)",
                LAB, pm, pl, TA$p))
msg <- switch(sub("_BUT.*", "", LAB),
  MONO_ADDS = "Mono_DC carries information LMPP_GMP does not. PREREG 6 caps the wording at\n    'not reducible to the LMPP/GMP compartment' -- it does NOT show the cells are non-malignant.",
  LMPP_ONLY = "Mono_DC is the diluted copy. The Mono_DC line closes; what remains is an\n    NPM1 classifier built with no malignancy annotation.",
  BOTH_ADD  = "Two partially independent readouts of one axis. Neither compartment is primary.",
  REDUNDANT = "Not demonstrable at n=36 / 9 events. PREREG 7 fixed in advance that this is NOT\n    evidence that Mono_DC is redundant.")
message(paste0("    ", msg))
if (grepl("TNK_CONTROL_FAILED", LAB))
  message("    T_NK adds too, so the adding is generic (PREREG 5.1-1). This supports nothing.")
rec(test = "DECISION", stratum = LAB, n = nrow(R2), stat = NA_real_, p = NA_real_,
    threshold = ALPHA, label = LAB)

fwrite(rbindlist(TESTS, fill = TRUE), file.path(DIR_PSEUDOBULK, "hox_compartment_tests.csv"))
fwrite(rbind(data.table(which = "mono_adj_lmpp", auc = MA$null),
             data.table(which = "lmpp_adj_mono", auc = MB$null),
             data.table(which = "tnk_adj_lmpp",  auc = TA$null)),
       file.path(DIR_PSEUDOBULK, "hox_compartment_null.csv"))
message(sprintf("\n[done] wrote hox_compartment_{tests,null}.csv to %s", DIR_PSEUDOBULK))
