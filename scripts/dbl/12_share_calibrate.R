#!/usr/bin/env Rscript
# 12_share_calibrate.R ----
# Doublet-pairing redesign (share design), calibration stage. Only mode so far: --reuse-smoke = stage 0b (plan:
# reports/2026-11_dbl_pairing_smoke/redesign_2026-11-11/STAGE0B_plan.md; user decision 2026-11-12).
# INPUT  : DBL_TAB_DIR/dbl_typing_injected_devel_<tag>_smoke.csv.gz (stage 04; tags DBL_OE_CAL_TAGS; held-out SYNTHETIC rows only),
#          DBL_LANE_DIR/<lane>.rds + <lane>__labels{,_bmm}.csv.gz (class-A cells only, as stage 10),
#          DBL_TAB_DIR/dbl_nbs_rates_{tune,confirm}_smoke.csv (stage 10 --rates-only; the known answer for M1)
# OUTPUT : DBL_TAB_DIR/dbl_share_rates_smoke.csv   tag x lane x pair type x truth x flag (M1, D, Y_w): n, k
#          DBL_TAB_DIR/dbl_checks_12_smoke.csv     validity checks (the rates are used only if all pass)
# What it does: rebuilds stage 04's held-out synthetic doublets exactly as stage 10 does (same seeds, code copied verbatim),
#   and flags each one with the share design's numerator Y = M1 (stage 10 verbatim) and denominator
#   D = M1 OR [called AND (unambiguous origin naming the LEUK lineage + one X in P, OR NNLS top two groups = {LEUK, X in P},
#   minor >= DBL_SHARE_D_MINOR, gap >= DBL_SHARE_D_GAP)], P = DBL_SHARE_P; plus Y_w = M1 inside the lane's LEUK-share window.
# MASKING: stage 04's tables also hold re-entered real singlets ("singlet:*"); those rows are dropped at read, and the rebuilt
#   injection keeps synthetic doublets only, so no real barcode is ever flagged. No per-barcode file is written.
# Usage  : Rscript scripts/dbl/12_share_calibrate.R --reuse-smoke [--force] [--cores=8]

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here); library(parallel); library(nnls) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
MD5_FROZEN <- c(config_dbl = "91d7f54db724b7d87e9a82a4f3838e5e", config_paths = "b2c7ca9b7d65f2c39774e6ee2cca2523")   # configs as tagged (prereg-dbl-0b-v1)
MD5 <- unname(tools::md5sum(c(here::here("scripts", "dbl", "config_dbl.R"), here::here("scripts", "config", "config_paths.R"))))
if (!identical(MD5, unname(MD5_FROZEN))) stop("config md5 differs from the frozen literals: fix and resubmit (not a rerun)")
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))   # the dated copy that runs
ARGS <- commandArgs(TRUE); stopifnot("--reuse-smoke" %in% ARGS); FORCE <- "--force" %in% ARGS
argv <- function(key, default) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), ARGS, value = TRUE)); if (length(v)) v else default }
CORES <- as.integer(argv("cores", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
OUT <- file.path(DBL_TAB_DIR, c(rates = "dbl_share_rates_smoke.csv", checks = "dbl_checks_12_smoke.csv")); names(OUT) <- c("rates", "checks")
stopifnot(all(startsWith(OUT, DBL_SMOKE_DIR)))                     # outputs never inside the repo
if (any(file.exists(OUT)) && !FORCE) stop("outputs exist (use --force)")
ts  <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-84s expect %-8s got %-8s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **"))
  invisible(ok) }
chk(sprintf("md5 of the running copy %s (recorded)", basename(SELF)), unname(tools::md5sum(SELF)), unname(tools::md5sum(SELF)))
STATES <- DBL_LEUK_STATES; LIN <- c(STATES, "LEUK_other")
TAGS <- DBL_OE_CAL_TAGS                                            # tuning folds 1-2 + fold 3, as stage 11
## -- copied verbatim from 09_recall_variants.R (lines 53-103): stage 04's construction, NNLS with merged LEUK lineage ----
nz <- function(z) fifelse(z == "", NA_character_, z)
lin_origin <- function(origin) vapply(strsplit(ifelse(is.na(origin), "", origin), "+", fixed = TRUE), function(z)
  length(z) == 2 && "T" %in% z && sum(z %in% LIN) == 1, TRUE)

## -- stage 04's construction (verbatim, as in 06/07) ----
sample_cells <- function(pool, n) if (length(pool) >= n) sample(pool, n) else sample(pool, n, replace = TRUE)
build_pairs <- function(src, n, typ, comp_a, comp_b, st = NA_character_) {
  pa <- src[[comp_a]]; pb <- src[[comp_b]]
  if (length(pa) == 0 || length(pb) == 0 || n == 0) return(NULL)
  data.table(type = typ, true_state = st, a = sample_cells(pa, n), b = sample_cells(pb, n), first = comp_a) }
build_leuk_pairs <- function(src, n, typ, comp_a) {
  av <- lengths(src[intersect(STATES, names(src))]); av <- av[av > 0]
  if (!length(av) || length(src[[comp_a]]) == 0 || n == 0) return(NULL)
  st <- sample(names(av), n, replace = TRUE, prob = av)
  data.table(type = typ, true_state = NA_character_, a = sample_cells(src[[comp_a]], n), b = vapply(st, function(s) sample_cells(src[[s]], 1), ""), first = comp_a) }
build_injection <- function(src, run_id) {
  prs <- rbindlist(c(
    lapply(STATES, function(s) build_pairs(src, DBL_SYN_PER_STATE, "T|LEUK", "T", s, s)),
    lapply(c("Tamb", "NK", "B", "Ery", "MonoDC", "LEUK_other", "T"), function(cp) build_pairs(src, DBL_SYN_PER_OTHER, paste0("T|", cp), "T", cp)),
    lapply(c("Tamb", "NK", "B", "Ery", "MonoDC"), function(cp) build_leuk_pairs(src, DBL_SYN_PER_OTHER, paste0(cp, "|LEUK*"), cp)),
    list(build_leuk_pairs(src[intersect(names(src), c(STATES, "LEUK_other"))] |> (\(s) { s[["LEUKx"]] <- unlist(s[intersect(STATES, names(s))]); s })(),
                          DBL_SYN_PER_OTHER, "LEUK*|LEUK*", "LEUKx"))), fill = TRUE)
  comps <- names(src)[lengths(src) > 0]
  per   <- max(1L, floor(DBL_SYN_SINGLETS / length(comps)))
  sng   <- rbindlist(lapply(comps, function(cp) data.table(type = paste0("singlet:", cp), true_state = NA_character_,
                             a = if (length(src[[cp]]) >= per) sample(src[[cp]], per) else src[[cp]], b = NA_character_, first = cp)))
  inj <- rbind(prs, sng, fill = TRUE); inj[, syn := sprintf("SYN_%s_%04d", run_id, .I)]; inj }
## -- NNLS with merged LEUK lineage (stage 06's basis) ----
make_basis <- function(M, ref) {
  keep <- ref[, .N, by = compartment][N >= 5L, compartment]; ref <- ref[compartment %in% keep]
  ok <- !grepl("^MT-|^RP[SL]|^HB[^P]", rownames(M))
  S <- sapply(keep, function(cp) Matrix::rowSums(M[ok, ref[compartment == cp, cell], drop = FALSE]))
  CP <- sweep(S, 2, colSums(S), "/") * 1e6
  top_finite <- function(v, n) head(names(sort(v[is.finite(v)], decreasing = TRUE)), n)
  G <- character(0)
  for (cp in keep) { o <- setdiff(keep, cp); lfc <- log2((CP[, cp] + 1) / (rowMeans(CP[, o, drop = FALSE]) + 1)); lfc[CP[, cp] < 10] <- -Inf
    G <- union(G, top_finite(lfc, 50L)) }
  ls <- intersect(STATES, keep)
  if (length(ls) >= 2) for (i in seq_along(ls)) for (j in seq_along(ls)) if (i < j) {
    d <- abs(log2((CP[, ls[i]] + 1) / (CP[, ls[j]] + 1))); d[pmax(CP[, ls[i]], CP[, ls[j]]) < 10] <- -Inf
    G <- union(G, top_finite(d, 20L)) }
  P <- S[G, , drop = FALSE] + 0.5; P <- sweep(P, 2, colSums(P), "/"); list(G = G, P = P) }
nnls_lin <- function(X, B) {   # merged-lineage call for every column of X
  Y <- as.matrix(X); Y <- sweep(Y, 2, pmax(colSums(Y), 1), "/")
  w <- 1 / sqrt(rowMeans(B$P) + 1e-6); Pw <- B$P * w; cn <- colnames(B$P)
  grp <- ifelse(cn %in% LIN, "LEUK", cn)
  vapply(seq_len(ncol(Y)), function(j) {
    b <- nnls::nnls(Pw, Y[, j] * w)$x; b <- b / max(sum(b), 1e-12)
    gw <- sort(tapply(b, grp, sum), decreasing = TRUE)
    length(gw) >= 2 && setequal(names(gw)[1:2], c("T", "LEUK")) && min(gw[1:2]) >= 0.15 &&
      (gw[2] - (if (length(gw) > 2) gw[3] else 0)) >= 0.05 }, TRUE) }
## -- new for the share design: the D flag's two typing routes ----
lx_origin <- function(origin) vapply(strsplit(ifelse(is.na(origin), "", origin), "+", fixed = TRUE), function(z)   # LEUK + one X in P
  length(z) == 2 && sum(z %in% LIN) == 1 && all(z[!z %in% LIN] %in% DBL_SHARE_P), TRUE)
nnls_lx <- function(X, B) {   # nnls_lin generalised: top two groups are {LEUK, X in P}, with the same minor / gap rule
  Y <- as.matrix(X); Y <- sweep(Y, 2, pmax(colSums(Y), 1), "/")
  w <- 1 / sqrt(rowMeans(B$P) + 1e-6); Pw <- B$P * w; cn <- colnames(B$P)
  grp <- ifelse(cn %in% LIN, "LEUK", cn)
  vapply(seq_len(ncol(Y)), function(j) {
    b <- nnls::nnls(Pw, Y[, j] * w)$x; b <- b / max(sum(b), 1e-12)
    gw <- sort(tapply(b, grp, sum), decreasing = TRUE)
    length(gw) >= 2 && "LEUK" %in% names(gw)[1:2] && all(setdiff(names(gw)[1:2], "LEUK") %in% DBL_SHARE_P) &&
      min(gw[1:2]) >= DBL_SHARE_D_MINOR && (gw[2] - (if (length(gw) > 2) gw[3] else 0)) >= DBL_SHARE_D_GAP }, TRUE) }

## -- [1] per lane: rebuild the held-out synthetic doublets (stage 10's construction), flag Y, D, Y_w ----
FL <- list()
for (L in DBL_SMOKE_LANES) {
  ts("[1] lane ", L)
  x <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds"))); M <- x$counts
  lab  <- fread(file.path(DBL_LANE_DIR, paste0(L, "__labels.csv.gz")))[match(colnames(M), cell)]
  labb <- fread(file.path(DBL_LANE_DIR, paste0(L, "__labels_bmm.csv.gz")))
  A <- lab[class == "A"]
  chk(sprintf("%s: class-A compartments identical in the knn and bmm label files", L), TRUE,
      identical(A$compartment, labb[match(A$cell, cell), compartment]))
  pools <- split(A$cell, A$compartment)
  cpm <- function(cells) { s <- M[, cells, drop = FALSE]; Matrix::rowSums(s) / sum(s) * 1e6 }
  badg <- grepl("^MT-|^RP[SL]|^HB[^P]", rownames(M))
  cT <- cpm(pools[["T"]]); cL <- cpm(unlist(pools[intersect(STATES, names(pools))]))
  GT4 <- rownames(M)[!badg & cT > 5 & cL < 1]; GL4 <- rownames(M)[!badg & cL > 5 & cT < 1]
  leuk_share <- function(Mx) { l <- Matrix::colSums(Mx[GL4, , drop = FALSE]); t <- Matrix::colSums(Mx[GT4, , drop = FALSE]); ifelse(l + t > 0, l / (l + t), NA_real_) }
  set.seed(DBL_SEED_SYN)
  A[, fold := ((sample(.N) - 1L) %% DBL_N_FOLDS) + 1L, by = compartment]
  I04 <- lapply(TAGS, function(TAG) { d <- fread(file.path(DBL_TAB_DIR, paste0("dbl_typing_injected_devel_", TAG, "_smoke.csv.gz")))
    d <- d[lane == L & kind == "heldout" & !startsWith(type, "singlet:")]          # MASKING: real singlets dropped at read
    d[, `:=`(origin = nz(origin), true_state = nz(true_state))] }); names(I04) <- TAGS
  runs <- sort(unique(unlist(lapply(I04, `[[`, "run"))))
  ts("    runs: ", paste(runs, collapse = " "))
  one <- function(rid) {   # stage 10's rebuild (same seed and construction); singlet rows leave before any count is formed
    f <- as.integer(sub("^f([0-9]+)_b.*", "\\1", rid)); b <- as.integer(sub(".*_b([0-9]+)$", "\\1", rid))
    set.seed(DBL_SEED_SYN + f * 10L + b)
    held <- A[fold == f, cell]
    inj <- build_injection(split(A[cell %in% held, cell], A[cell %in% held, compartment]), rid)
    inj <- inj[!startsWith(type, "singlet:")]                       # synthetic doublets only (ids were assigned on the full set)
    Mi <- M[, inj$a, drop = FALSE] + M[, inj$b, drop = FALSE]; colnames(Mi) <- inj$syn
    B <- make_basis(M, A[!(cell %in% held), .(cell, compartment)])
    inj[, nnls := nnls_lin(Mi[B$G, , drop = FALSE], B)]
    inj[, nnls_lx := nnls_lx(Mi[B$G, , drop = FALSE], B)]
    inj[, share_rb := leuk_share(Mi)]
    inj[, tid := Matrix::colSums(Mi[intersect(DBL_T_ID_GENES, rownames(Mi)), , drop = FALSE]) >= DBL_T_ID_MIN_UMI]
    inj[, .(run = rid, syn, type, true_state, nnls, nnls_lx, share_rb, tid)] }
  RB <- mclapply(runs, one, mc.cores = min(CORES, length(runs)))   # forked workers; each run sets its own seed
  stopifnot(!vapply(RB, inherits, TRUE, "try-error"), !vapply(RB, is.null, TRUE)); RB <- rbindlist(RB)   # a killed fork returns NULL
  for (TAG in TAGS) {
    m <- merge(I04[[TAG]], RB, by = c("run", "syn"), suffixes = c("", ".rb"))
    chk(sprintf("%s %s: rebuilt injection == stage 04's synthetic rows", TAG, L), nrow(I04[[TAG]]), nrow(m))
    chk(sprintf("%s %s: rebuilt injection == stage 04's types", TAG, L), TRUE, m[, all(type == type.rb)])
    chk(sprintf("%s %s: rebuilt injection == stage 04's LEUK share (1e-9, NA-aware)", TAG, L), TRUE,
        m[, all((is.na(share) & is.na(share_rb)) | (!is.na(share) & !is.na(share_rb) & abs(share - share_rb) < 1e-9))])
    m[, lin_o := lin_origin(origin)]                                                 # stage 10 verbatim from here ...
    m[, `:=`(V1a = called %in% TRUE & ambiguous %in% FALSE & lin_o, V3 = called %in% TRUE & nnls)]
    m[, V5 := V1a | V3]
    m[, M1 := V5 & tid]                                                              # ... to here (its singlet branch is not needed)
    m[, D := M1 | (called %in% TRUE & ((ambiguous %in% FALSE & lx_origin(origin)) | nnls_lx))]
    chk(sprintf("%s %s: every flagged barcode is one of stage 04's synthetic doublets", TAG, L), TRUE, all(m$syn %in% I04[[TAG]]$syn) && !any(startsWith(m$type, "singlet:")))
    FL[[paste(TAG, L)]] <- m[, .(tag = TAG, lane = L, run, syn, type, true_state, share_rb, M1, D)] }
  rm(x, M, RB); gc(FALSE) }
FL <- rbindlist(FL)
FL[, truth := fifelse(type == "T|LEUK", true_state, fifelse(type == "T|LEUK_other", "LEUK_other", NA_character_))]   # as stage 10
WIN <- FL[type == "T|LEUK", .(lo = quantile(share_rb, DBL_SHARE_WIN_Q[1], na.rm = TRUE, type = 7),
                              hi = quantile(share_rb, DBL_SHARE_WIN_Q[2], na.rm = TRUE, type = 7)), by = lane]
FL <- merge(FL, WIN, by = "lane"); FL[, Y_w := M1 & !is.na(share_rb) & share_rb >= lo & share_rb <= hi]
print(WIN)

## -- [2] checks: M1 reproduces stage 10's stored counts; M1 inside D; synthetic rows only ----
CNT <- melt(FL, id.vars = c("tag", "lane", "type", "truth"), measure.vars = c("M1", "D", "Y_w"), variable.name = "flag", value.name = "v")[
  , .(n = .N, k = sum(v)), by = .(tag, lane, type, truth, flag)]
ref <- rbind(fread(file.path(DBL_TAB_DIR, "dbl_nbs_rates_tune_smoke.csv")), fread(file.path(DBL_TAB_DIR, "dbl_nbs_rates_confirm_smoke.csv")))[
  variant == "M1" & kind == "pair" & tag %in% TAGS, .(tag, lane, type, truth = nz(truth), n, k)]
cmp <- merge(ref, CNT[flag == "M1"], by = c("tag", "lane", "type", "truth"), all = TRUE, suffixes = c(".10", ""))
chk("M1 k/n == stage 10's stored counts, every tag x lane x pair type x truth", TRUE, cmp[, all(!is.na(n.10) & !is.na(n) & n == n.10 & k == k.10)])
chk("M1 is a subset of D for every synthetic doublet", TRUE, FL[, all(!M1 | D)])
print(dcast(CNT[, .(n = sum(n), k = sum(k)), by = .(lane = sub("_.*", "", lane), type, flag)][, rate := round(k / n, 3)], lane + type ~ flag, value.var = "rate"))
CHK <- rbindlist(CHK)
fwrite(CNT, OUT[["rates"]]); fwrite(CHK, OUT[["checks"]])
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
if (!all(CHK$pass)) quit(status = 1)
