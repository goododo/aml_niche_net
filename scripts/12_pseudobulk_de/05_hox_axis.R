#!/usr/bin/env Rscript
# 05_hox_axis.R ----
# INPUT  : DIR_PSEUDOBULK/pb_sample_manifest.csv, /depth/*.csv, /de/*.csv, /hitlist_frozen.csv
#          PB_RDS_DIR/<ds>/<sample>__pseudobulk.rds ; 05_ccc/undercall_contamination.csv ;
#          00_project/metadata/GSE116256__samples.tsv
# OUTPUT : DIR_PSEUDOBULK/hox_axis_{scores,genes,tests,null}.csv
# WHAT IT DOES : scores every sample on the HOXA/HOXB-MEIS1 axis inside three hierarchy bins, then
#          asks whether that score is blast burden (vs clinical blast %) and whether it tracks NPM1.
# Usage  : Rscript scripts/12_pseudobulk_de/05_hox_axis.R
#
# EVERYTHING HERE IS PRE-REGISTERED in PREREGISTRATION_hox_axis.md, committed before this file
# existed. This script implements that document and decides nothing about the outcome.
#
# FOUR CLAUSES THE PRE-REGISTRATION LEFT QUALITATIVE. Each is made numeric HERE, i.e. before any
# result exists, and is listed so the choice is auditable rather than retrofitted:
#
#   (A) CELL GATE. sec 4 never fixed one. MIN_CELLS = 30 is used, identical to Discovery (PREREG
#       pseudobulk 5.1) and to 04_validation.R. Without it a 4-cell bin yields a score.
#   (B) sec 7.1 "T_NK tracks blast % about as well as Mono_DC" is read as
#       |rho_T_NK| >= 0.80 * |rho_Mono_DC| in BOTH testable strata.
#   (C) sec 4.3 z-scoring needs within-dataset variance. A gene with sd == 0 inside a dataset cannot
#       be z-scored, so it is dropped for that dataset under the same rule as sec 4.4, and the count
#       is printed.
#   (D) sec 10 orders "controls before tests", but controls 7.1 and 7.5 are defined RELATIVE to the
#       Q1 correlation and cannot precede it. Resolution: the three outcome-independent controls
#       (7.2 depth, 7.3 healthy floor, 7.4 set agreement) run first and can kill before blast % or
#       NPM1 is ever touched; 7.1 and 7.5 are computed together with Q1 and the Q1 decision label is
#       assigned only after they are evaluated. Q2 runs only if nothing has killed the axis. The
#       intent of sec 10 -- a killed axis is never converted into a claim -- is preserved.
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))

MIN_CELLS  <- CCC_MIN_CELLS_PER_OCCUPIED_BIN            # 30 -- clause (A)
BINS       <- c("Mono_DC", "LMPP_GMP", "T_NK")          # PREREG 4.6; B_Plasma has 2 genes, excluded
PRIMARY_BIN<- "Mono_DC"
RHO_PROXY  <- 0.60 ; RHO_FLAT <- 0.30                   # PREREG 5 decision rule
RHO_DEPTH  <- 0.50                                      # PREREG 7.2 kill
RHO_SETAGR <- 0.80                                      # PREREG 7.4 kill
SPEC_FRAC  <- 0.80                                      # clause (B)
N_NULL     <- 1000L                                     # PREREG 7.5
set.seed(20260918)

## ------------------------------------------------------------------ roster ----
MAN <- fread(file.path(DIR_PSEUDOBULK, "pb_sample_manifest.csv"))[in_analysis == TRUE]
message(sprintf("[roster] %d samples | Discovery %d AML + %d healthy | Validation %d AML + %d healthy",
                nrow(MAN), MAN[split == "Discovery" & arm == "AML", .N],
                MAN[split == "Discovery" & arm == "healthy", .N],
                MAN[split == "Validation" & arm == "AML", .N],
                MAN[split == "Validation" & arm == "healthy", .N]))
stopifnot(nrow(MAN) == 67L)                             # PREREG 3, frozen input

DEP <- rbindlist(lapply(list.files(file.path(DIR_PSEUDOBULK, "depth"), full.names = TRUE), fread))
PB  <- lapply(seq_len(nrow(MAN)), function(i)
  readRDS(file.path(PB_RDS_DIR, MAN$dataset[i], paste0(MAN$sample[i], "__pseudobulk.rds"))))
names(PB) <- MAN$sample
stopifnot(!anyDuplicated(MAN$sample))                   # scores are keyed on sample alone

## ------------------------------------------------------- gene sets (PREREG 4.1) ----
## Set L is a GENOMIC LOCUS rule applied to the bin's tested universe, so it contains genes our own
## DE did not pick. Set D is the same rule applied to the frozen hit list. L is primary precisely
## because it is not derived from our result.
univ_of_bin <- function(bin) {
  f <- list.files(file.path(DIR_PSEUDOBULK, "de"), pattern = sprintf("__%s\\.csv$", bin),
                  full.names = TRUE)
  stopifnot(length(f) == 1L)
  u <- fread(f)
  if (!"gene" %in% names(u)) setnames(u, 1L, "gene")
  u$gene
}
is_hox <- function(g) grepl("^HOX[AB]", g) | g %in% c("MEIS1", "PBX3")
HIT <- fread(file.path(DIR_PSEUDOBULK, "hitlist_frozen.csv"))
SETS <- list()
for (b in BINS) {
  u <- univ_of_bin(b)
  SETS[[b]] <- list(L = sort(u[is_hox(u)]),
                    D = sort(HIT[hierarchy_bin == b & is_hox(gene), gene]))
  message(sprintf("[set] %-9s universe %5d | Set L %2d | Set D %2d",
                  b, length(u), length(SETS[[b]]$L), length(SETS[[b]]$D)))
}
stopifnot(length(SETS[[PRIMARY_BIN]]$L) >= 10L)

## ----------------------------------------------- per (dataset,bin) log2 CPM ----
## PREREG 4.2 log2 CPM within the bin, that bin's own column sum as library size, prior count 1.
## PREREG 4.4 a gene absent from ANY sample of a dataset is dropped for that WHOLE dataset, so every
## sample of a dataset is scored on one identical gene list.
cpm_block <- function(ds, bin) {
  keep <- DEP[dataset == ds & hierarchy_bin == bin & n_cells >= MIN_CELLS, .(dataset, sample, n_cells)]
  keep <- keep[MAN, on = c("dataset", "sample"), nomatch = 0]
  if (nrow(keep) < 3L) return(NULL)
  gl <- Reduce(intersect, lapply(keep$sample, function(s) PB[[s]]$genes))
  gl <- intersect(gl, univ_of_bin(bin))
  m <- vapply(keep$sample, function(s) {
    v <- PB[[s]]$counts[gl, bin]; log2(v / sum(v) * 1e6 + 1)
  }, numeric(length(gl)))
  dimnames(m) <- list(gl, keep$sample)
  list(cpm = m, meta = keep)
}

## PREREG 4.3/4.5: z within dataset across ALL samples of that dataset in that bin (healthy anchors
## the low end), then an UNWEIGHTED mean. Weighting by logFC would import the Discovery effect sizes
## and re-open exactly the circularity Set L was chosen to close.
zscore_rows <- function(m) {
  s <- apply(m, 1L, sd)
  m <- m[s > 0, , drop = FALSE]                        # clause (C)
  t(scale(t(m)))
}
score_of <- function(blk, genes) {
  g <- intersect(genes, rownames(blk$cpm))
  if (length(g) < 3L) return(NULL)
  z <- zscore_rows(blk$cpm[g, , drop = FALSE])
  list(score = colMeans(z), genes = rownames(z))
}

SC <- list(); GENELOG <- list()
for (b in BINS) for (ds in unique(MAN$dataset)) {
  blk <- cpm_block(ds, b); if (is.null(blk)) next
  for (st in c("L", "D")) {
    r <- score_of(blk, SETS[[b]][[st]]); if (is.null(r)) next
    SC[[length(SC) + 1L]] <- data.table(dataset = ds, hierarchy_bin = b, set = st,
                                        sample = names(r$score), score = as.numeric(r$score),
                                        n_genes_used = length(r$genes))[
      blk$meta[, .(sample, n_cells)], on = "sample", nomatch = 0]
    GENELOG[[length(GENELOG) + 1L]] <- data.table(dataset = ds, hierarchy_bin = b, set = st,
                                                  gene = r$genes,
                                                  n_dropped = length(SETS[[b]][[st]]) - length(r$genes))
  }
}
SC <- rbindlist(SC); GEN <- rbindlist(GENELOG)
SC <- MAN[, .(sample, arm, split)][SC, on = "sample"]
SC[, lib_size := vapply(seq_len(.N), function(i)
  sum(PB[[sample[i]]]$counts[, hierarchy_bin[i]]), numeric(1))]
fwrite(SC, file.path(DIR_PSEUDOBULK, "hox_axis_scores.csv"))
fwrite(GEN, file.path(DIR_PSEUDOBULK, "hox_axis_genes.csv"))
message(sprintf("[score] %d (sample,bin,set) rows | genes used in %s Set L: %s",
                nrow(SC), PRIMARY_BIN,
                paste(unique(SC[hierarchy_bin == PRIMARY_BIN & set == "L", n_genes_used]), collapse = "/")))

TESTS <- list()
rec <- function(...) TESTS[[length(TESTS) + 1L]] <<- data.table(...)
sp <- function(x, y) { k <- is.finite(x) & is.finite(y)
  if (sum(k) < 4L) return(list(rho = NA_real_, p = NA_real_, n = sum(k)))
  h <- suppressWarnings(cor.test(x[k], y[k], method = "spearman"))
  list(rho = unname(h$estimate), p = h$p.value, n = sum(k)) }

wide <- dcast(SC, sample + dataset + arm + split + hierarchy_bin ~ set,
              value.var = c("score", "n_cells", "lib_size"))

## ===================== CONTROLS 7.2 / 7.3 / 7.4 -- outcome-independent, run first ==============
message("\n===== CONTROLS (outcome-independent; these run before blast % or NPM1 is touched) =====")
KILL <- character(0)

## 7.4 set agreement -- if L and D disagree, the two definitions are not one axis and neither is carried
w <- wide[hierarchy_bin == PRIMARY_BIN]
a <- sp(w$score_L, w$score_D)
message(sprintf("[7.4 set agreement ] Spearman(Set L, Set D) in %s = %+.3f (n=%d)  [kill if < %.2f]",
                PRIMARY_BIN, a$rho, a$n, RHO_SETAGR))
rec(control = "7.4_set_agreement", stratum = PRIMARY_BIN, n = a$n, rho = a$rho, p = a$p,
    threshold = RHO_SETAGR, passed = isTRUE(a$rho >= RHO_SETAGR))
if (!isTRUE(a$rho >= RHO_SETAGR)) KILL <- c(KILL, "7.4 set agreement")

## 7.2 depth -- a score that is really a sequencing-depth readout
for (v in c("lib_size_L", "n_cells_L")) {
  x <- if (v == "lib_size_L") log10(w$lib_size_L) else w$n_cells_L
  d <- sp(w$score_L, x)
  message(sprintf("[7.2 depth        ] Spearman(score, %-11s) = %+.3f (n=%d)  [kill if |rho| > %.2f]",
                  sub("_L$", "", v), d$rho, d$n, RHO_DEPTH))
  rec(control = "7.2_depth", stratum = sub("_L$", "", v), n = d$n, rho = d$rho, p = d$p,
      threshold = RHO_DEPTH, passed = isTRUE(abs(d$rho) <= RHO_DEPTH))
  if (isTRUE(abs(d$rho) > RHO_DEPTH)) KILL <- c(KILL, paste("7.2 depth /", v))
}

## 7.3 healthy floor -- reported as AUC, per dataset (cross-dataset pooling is meaningless after 4.3)
for (ds in unique(w$dataset)) {
  hh <- w[dataset == ds & arm == "healthy", score_L]; aa <- w[dataset == ds & arm == "AML", score_L]
  if (!length(hh) || !length(aa)) next
  u <- suppressWarnings(wilcox.test(aa, hh))
  auc <- unname(u$statistic) / (length(aa) * length(hh))
  message(sprintf("[7.3 healthy floor] %-10s AML %2d vs healthy %2d | AUC = %.3f (p=%.3g)",
                  ds, length(aa), length(hh), auc, u$p.value))
  rec(control = "7.3_healthy_floor", stratum = ds, n = length(aa) + length(hh), rho = auc,
      p = u$p.value, threshold = NA_real_, passed = NA)
}

if (length(KILL)) {
  message(sprintf("\n[KILLED] %s -- PREREG 10.4: sections 5 and 6 are NOT run. That is the result.",
                  paste(KILL, collapse = "; ")))
  fwrite(rbindlist(TESTS, fill = TRUE), file.path(DIR_PSEUDOBULK, "hox_axis_tests.csv"))
  quit(save = "no", status = 0)
}
message("[controls] 7.2 / 7.3 / 7.4 clear -- proceeding to Q1")

## ================================ Q1 + CONTROLS 7.1 / 7.5 ======================================
U <- fread(file.path(TAB_DIR, "05_ccc", "undercall_contamination.csv"))[
  , .(sample, blast_pct_clinical, malignant_frac)]
Q <- U[wide[arm == "AML"], on = "sample", nomatch = 0][is.finite(blast_pct_clinical)]
## n is SAMPLES IN THE PRIMARY BIN, not rows: Q holds one row per (sample,bin) and counting those
## would promote a 3-sample dataset into a testable stratum.
NPB <- Q[hierarchy_bin == PRIMARY_BIN, .N, by = dataset]
STRATA <- NPB[N >= 5L, dataset]
message(sprintf("\n===== Q1: is the axis blast burden? =====\n  strata: %s | tested (n>=5): %s",
                paste(sprintf("%s n=%d", NPB$dataset, NPB$N), collapse = ", "),
                paste(STRATA, collapse = ", ")))
## PREREG 5 quoted n=9 (GSE116256) and n=18 (GSE185381) from the ROSTER. Clause (A)'s 30-cell gate
## is applied on top of that, so the tested n is smaller. The loss is printed, not absorbed.
for (ds in NPB$dataset) {
  roster <- MAN[dataset == ds & arm == "AML", .N]
  havebp <- U[sample %in% MAN[dataset == ds & arm == "AML", sample] & is.finite(blast_pct_clinical), .N]
  message(sprintf("  [gate] %-10s roster AML %2d -> with blast%% %2d -> clearing %d-cell %s gate %2d",
                  ds, roster, havebp, MIN_CELLS, PRIMARY_BIN, NPB[dataset == ds, N]))
}

q1 <- list()
for (ds in Q[, unique(dataset)]) {
  g <- Q[dataset == ds & hierarchy_bin == PRIMARY_BIN]
  r <- sp(g$score_L, g$blast_pct_clinical); rc <- sp(g$malignant_frac, g$blast_pct_clinical)
  tested <- ds %in% STRATA
  message(sprintf("  %-10s n=%2d | HOX axis rho %+.3f (p=%.3f) | inferCNV malignant_frac rho %+.3f (p=%.3f)%s",
                  ds, r$n, r$rho, r$p, rc$rho, rc$p, if (tested) "" else "   [n<5, reported not tested]"))
  ## `passed` means "cleared the threshold", so it must NOT carry stratum testability here -- a
  ## rho of -0.43 is not a pass. Testability is encoded in the control name instead.
  rec(control = if (tested) "Q1_blast" else "Q1_blast_NOT_TESTED_n_below_5",
      stratum = ds, n = r$n, rho = r$rho, p = r$p, threshold = RHO_PROXY,
      passed = if (tested) isTRUE(r$rho >= RHO_PROXY) else NA)
  rec(control = "Q1_blast_inferCNV_comparator", stratum = ds, n = rc$n, rho = rc$rho, p = rc$p,
      threshold = NA_real_, passed = NA)
  if (tested) q1[[ds]] <- r$rho
}
message("  (van Galen signature caller on this same ruler, recorded 2026-09-16: rho = -0.165)")
message(sprintf("  exact two-sided 0.05 boundary: |rho| >= 0.700 at n=9, ~0.47 at n=18"))

## 7.1 bin specificity -- if T_NK tracks blast % about as well, this is ambient RNA / doublets
message("\n[7.1 bin specificity] same correlation in the control bins:")
spec <- list()
for (b in BINS) for (ds in STRATA) {
  g <- Q[dataset == ds & hierarchy_bin == b]
  if (!nrow(g)) next
  r <- sp(g$score_L, g$blast_pct_clinical)
  message(sprintf("   %-10s %-9s rho %+.3f (n=%d)", ds, b, r$rho, r$n))
  rec(control = "7.1_bin_specificity", stratum = paste(ds, b), n = r$n, rho = r$rho, p = r$p,
      threshold = NA_real_, passed = NA)
  spec[[paste(ds, b)]] <- r$rho
}
tnk_fail <- vapply(STRATA, function(ds) {
  a <- spec[[paste(ds, "T_NK")]]; m <- spec[[paste(ds, PRIMARY_BIN)]]
  isTRUE(is.finite(a) && is.finite(m) && abs(a) >= SPEC_FRAC * abs(m))
}, logical(1))
if (length(tnk_fail) && all(tnk_fail)) KILL <- c(KILL, "7.1 bin specificity (T_NK tracks blast% too)")
message(sprintf("   clause (B): T_NK >= %.2f x Mono_DC in every stratum ? %s",
                SPEC_FRAC, if (length(tnk_fail) && all(tnk_fail)) "YES -> KILL" else "no"))

## 7.5 random-set null -- size- and expression-decile-matched. If any 19 genes would do, HOX is decoration.
ds0 <- STRATA[which.max(Q[dataset %in% STRATA, .N, by = dataset][match(STRATA, dataset), N])]
blk <- cpm_block(ds0, PRIMARY_BIN)
obs <- q1[[ds0]]
gl  <- rownames(blk$cpm); mu <- rowMeans(blk$cpm)
dec <- cut(mu, quantile(mu, seq(0, 1, .1)), include.lowest = TRUE, labels = FALSE)
real <- intersect(SETS[[PRIMARY_BIN]]$L, gl)
tab  <- table(dec[match(real, gl)])
smp_ids <- Q[dataset == ds0 & hierarchy_bin == PRIMARY_BIN, sample]
bp      <- Q[dataset == ds0 & hierarchy_bin == PRIMARY_BIN, blast_pct_clinical]
nullr <- vapply(seq_len(N_NULL), function(i) {
  pick <- unlist(lapply(names(tab), function(d)
    sample(gl[dec == as.integer(d)], tab[[d]], replace = FALSE)))
  z <- zscore_rows(blk$cpm[pick, , drop = FALSE])
  s <- colMeans(z)[smp_ids]
  suppressWarnings(cor(s, bp, method = "spearman"))
}, numeric(1))
p95 <- quantile(abs(nullr), .95, na.rm = TRUE)
pemp <- mean(abs(nullr) >= abs(obs), na.rm = TRUE)
message(sprintf("\n[7.5 random null  ] %s, %d matched random sets of %d genes", ds0, N_NULL, length(real)))
message(sprintf("   observed |rho| %.3f | null 95th pct %.3f | empirical p = %.4f  [kill if observed <= p95]",
                abs(obs), p95, pemp))
rec(control = "7.5_random_null", stratum = ds0, n = length(smp_ids), rho = obs, p = pemp,
    threshold = p95, passed = isTRUE(abs(obs) > p95))
fwrite(data.table(dataset = ds0, rho = nullr), file.path(DIR_PSEUDOBULK, "hox_axis_null.csv"))
if (!isTRUE(abs(obs) > p95)) KILL <- c(KILL, "7.5 random-set null")

## Q1 decision -- assigned only now, after 7.1 and 7.5
rh <- unlist(q1)
stopifnot(length(rh) >= 1L, all(is.finite(rh)))   # a stratum that cannot be correlated must not vote
Q1 <- "INDETERMINATE"
if (length(KILL)) {
  Q1 <- "KILLED"
} else if (all(rh >= RHO_PROXY)) {
  Q1 <- "BLAST-PROXY"
} else if (all(abs(rh) < RHO_FLAT)) {
  Q1 <- "NOT-BLAST"
}
message(sprintf("\n>>> Q1 = %s   (rho: %s)", Q1,
                paste(sprintf("%s %+.3f", names(rh), rh), collapse = ", ")))
rec(control = "Q1_DECISION", stratum = Q1, n = length(rh), rho = NA_real_, p = NA_real_,
    threshold = RHO_PROXY, passed = NA)

## ========================================= Q2: NPM1 ===========================================
if (length(KILL)) {
  message(sprintf("[KILLED] %s -- PREREG 10.4: section 6 is NOT run.", paste(KILL, collapse = "; ")))
} else {
  G <- fread(file.path(PROJECT_ROOT, "00_project", "metadata", "GSE116256__samples.tsv"))[, .(sample, mutations)]
  V <- wide[dataset == "GSE116256" & arm == "AML" & hierarchy_bin == PRIMARY_BIN][G, on = "sample", nomatch = 0]
  V[, npm1 := grepl("NPM1", mutations)]
  setorder(V, -score_L)
  message(sprintf("\n===== Q2: NPM1 (%d mut vs %d wt, GSE116256 held-out arm) =====",
                  sum(V$npm1), sum(!V$npm1)))
  message(sprintf("   the table IS the result at n=%d:", nrow(V)))
  for (i in seq_len(nrow(V)))
    message(sprintf("   %-12s score %+6.3f  %-7s  %s", V$sample[i], V$score_L[i],
                    if (V$npm1[i]) "NPM1mut" else "NPM1wt", V$mutations[i]))
  wt <- suppressWarnings(wilcox.test(V[npm1 == TRUE, score_L], V[npm1 == FALSE, score_L],
                                     alternative = "greater", exact = TRUE))
  above <- all(V[npm1 == TRUE, score_L] > median(V[npm1 == FALSE, score_L]))
  Q2 <- if (wt$p.value < 0.05 && above) "SUPPORT" else "NULL"
  minp <- 1 / choose(nrow(V), sum(V$npm1))
  message(sprintf("\n   one-sided Wilcoxon p = %.4f | minimum attainable at this n = 1/%d = %.4f%s | all mutants above wt median: %s",
                  wt$p.value, choose(nrow(V), sum(V$npm1)), minp,
                  if (isTRUE(all.equal(wt$p.value, minp))) "  <- p IS AT THE FLOOR (perfect separation)" else "",
                  above))
  message(sprintf(">>> Q2 = %s%s", Q2,
                  if (Q2 == "SUPPORT") "  -- report as 'NPM1-mut and/or FLT3-ITD' (PREREG 6 confound)"
                  else "  -- a null at 4 vs 5 is NOT evidence of absence (PREREG 6)"))
  rec(control = "Q2_NPM1", stratum = "GSE116256", n = nrow(V), rho = NA_real_, p = wt$p.value,
      threshold = 0.05, passed = (Q2 == "SUPPORT"))
  tp53 <- V[grepl("TP53", mutations), .(sample, score_L, rank = match(sample, V$sample))]
  message(sprintf("   TP53 direction check (n=%d here, no inference): %s", nrow(tp53),
                  paste(sprintf("%s rank %d/%d", tp53$sample, tp53$rank, nrow(V)), collapse = ", ")))

  ## ---- POST-HOC SENSITIVITY, added 2026-09-18 AFTER the block above had been run ----
  ## PREREG 6 registered 4 mutant vs 5 wild-type. The realised test is 4 vs 3: clause (A)'s 30-cell
  ## gate removes AML916-D0 (14 Mono_DC cells) and AML707B-D0 (18), both NPM1-wt. That changes the
  ## minimum attainable one-sided p from 1/C(9,4)=0.0079 to 1/C(7,4)=0.0286, so the gated test can
  ## only ever reach significance by PERFECT separation. This block restores the registered n by
  ## dropping the gate for this dataset alone. It is a move TOWARD the pre-registration, not away
  ## from it, but it was written after seeing the gated result and is labelled accordingly. The two
  ## restored samples are pseudobulks of 14 and 18 cells -- noisier, which is why the gate exists.
  blkA <- local({
    keep <- DEP[dataset == "GSE116256" & hierarchy_bin == PRIMARY_BIN, .(dataset, sample, n_cells)][
      MAN, on = c("dataset", "sample"), nomatch = 0]
    gl <- Reduce(intersect, lapply(keep$sample, function(x) PB[[x]]$genes))
    gl <- intersect(gl, univ_of_bin(PRIMARY_BIN))
    m <- vapply(keep$sample, function(x) { v <- PB[[x]]$counts[gl, PRIMARY_BIN]
      log2(v / sum(v) * 1e6 + 1) }, numeric(length(gl)))
    dimnames(m) <- list(gl, keep$sample); list(cpm = m, meta = keep)
  })
  rA <- score_of(blkA, SETS[[PRIMARY_BIN]]$L)
  VA <- data.table(sample = names(rA$score), score_L = as.numeric(rA$score))[
    MAN[, .(sample, arm)], on = "sample", nomatch = 0][arm == "AML"][G, on = "sample", nomatch = 0]
  VA[, npm1 := grepl("NPM1", mutations)]
  setorder(VA, -score_L)
  message(sprintf("\n   [POST-HOC, gate removed] %d mut vs %d wt; min attainable p = %.4f",
                  sum(VA$npm1), sum(!VA$npm1), 1 / choose(nrow(VA), sum(VA$npm1))))
  for (i in seq_len(nrow(VA)))
    message(sprintf("   %-12s score %+6.3f  %-7s %s", VA$sample[i], VA$score_L[i],
                    if (VA$npm1[i]) "NPM1mut" else "NPM1wt",
                    if (VA$sample[i] %in% c("AML916-D0", "AML707B-D0")) "<- restored, below cell gate" else ""))
  wtA <- suppressWarnings(wilcox.test(VA[npm1 == TRUE, score_L], VA[npm1 == FALSE, score_L],
                                      alternative = "greater", exact = TRUE))
  message(sprintf("   one-sided Wilcoxon p = %.4f | all mutants above wt median: %s",
                  wtA$p.value, all(VA[npm1 == TRUE, score_L] > median(VA[npm1 == FALSE, score_L]))))
  rec(control = "Q2_NPM1_POSTHOC_no_gate", stratum = "GSE116256", n = nrow(VA), rho = NA_real_,
      p = wtA$p.value, threshold = 0.05, passed = NA)

  ## ---- POST-HOC CONTROL, added 2026-09-18: the random-set null for Q2 ----
  ## PREREG 7.5 specified a matched random-set null for Q1 and NOT for Q2. That is a gap in the
  ## pre-registration, and it matters most precisely here: a p sitting at the floor is the weakest
  ## possible significant result, so "would any 19 genes separate NPM1 just as well?" has to be
  ## answered or the HOX identity is decoration. Same matching as 7.5: size and expression decile.
  glA <- rownames(blkA$cpm); muA <- rowMeans(blkA$cpm)
  decA <- cut(muA, quantile(muA, seq(0, 1, .1)), include.lowest = TRUE, labels = FALSE)
  realA <- intersect(SETS[[PRIMARY_BIN]]$L, glA)
  tabA <- table(decA[match(realA, glA)])
  obsA <- wtA$p.value
  nullp <- vapply(seq_len(N_NULL), function(i) {
    pick <- unlist(lapply(names(tabA), function(d)
      sample(glA[decA == as.integer(d)], tabA[[d]], replace = FALSE)))
    sc <- colMeans(zscore_rows(blkA$cpm[pick, , drop = FALSE]))[VA$sample]
    suppressWarnings(wilcox.test(sc[VA$npm1], sc[!VA$npm1],
                                 alternative = "greater", exact = TRUE)$p.value)
  }, numeric(1))
  pe <- mean(nullp <= obsA, na.rm = TRUE)
  message(sprintf("   [POST-HOC null] %d matched random %d-gene sets | %.1f%% separate NPM1 at least as well | empirical p = %.3f",
                  N_NULL, length(realA), 100 * pe, pe))
  rec(control = "Q2_random_null_POSTHOC", stratum = "GSE116256", n = nrow(VA), rho = NA_real_,
      p = pe, threshold = 0.05, passed = isTRUE(pe < 0.05))
}

fwrite(rbindlist(TESTS, fill = TRUE), file.path(DIR_PSEUDOBULK, "hox_axis_tests.csv"))
message(sprintf("\n[done] wrote hox_axis_{scores,genes,tests,null}.csv to %s", DIR_PSEUDOBULK))
