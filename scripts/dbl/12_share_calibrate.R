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
# FREEZE 1 (tag prereg-dbl-share-f1): the rules of --calibrate are fixed in scripts/dbl/prereg_share/PREREGISTRATION_SHARE_F1.md
#   section 4 and REDESIGN_share_plan_2026-11-11.md section 4; the implementation details in STAGE2_plan.md (user approval
#   2026-11-12) and config_dbl.R section 11d. The --reuse-smoke block below is the frozen stage 0b code, unchanged.
# MODE --calibrate --lane=<lane>|array (Stage 2 step 2a; one lane per call; array = the SLURM_ARRAY_TASK_ID-th lane to calibrate)
#   INPUT  : DBL_LANE_DIR/<lane>.rds + <lane>__labels.csv.gz (01/02 --share); scripts/dbl/dbl_lanes_share.tsv;
#            DBL_TAB_DIR/dbl_share_stage1_eligibility.csv (Stage 1)
#   OUTPUT : DBL_SHARE_DIR/<lane>__cal[_quick][_stub].rds   synthetic barcodes only: parents, flags (called, origin, NNLS rules,
#            T-identity, LEUK share, Y = M1, D, DT = T-partner flag, Y_w), run facts, the plan, class-A composition, checks; plus ONE
#            aggregate row of the masked real run (threshold, called fraction, cluster-28 called fraction)
#   What it does: 5 folds x NB batches of held-out runs (one fifth of class A removed per run, stage 04's way, folds balanced within
#     donor x compartment). Each run injects within-donor synthetic doublets of the eligible donors at the design allocation
#     (DBL_SHARE_ALLOC; NB keeps every run <= DBL_SHARE_INJ_MAX of the lane) as knownDoublets, knownUse "discard", and flags them.
#     No real singlet is re-entered (strict masking). DT, the T-partner flag (qualifier psi_T; not defined in the freeze, fixed in
#     STAGE2_plan.md by symmetry with D): M1 OR [called AND (unambiguous origin {T, one Z in LEUK or P minus T} OR NNLS top two
#     groups {T, Z}, minor >= DBL_SHARE_D_MINOR, gap >= DBL_SHARE_D_GAP)].
#   --quick: batch 1 of every fold only (same seeds as the full run; Stage 2 test). --stub: a fake caller with scDblFinder's output
#     format, for the masking harness and fast structural tests (output marked _stub; never read as calibration).
#   Seeds: DBL_SEED_SHARE_CAL + 100 x lane index (folds) and + 100 x lane index + run (injection); lane index = position in the
#     committed lane table. A lane with > 99 runs shares seed values with later lanes; harmless (different cell pools), reported.
# Usage  : Rscript scripts/dbl/12_share_calibrate.R --reuse-smoke [--force] [--cores=8]
#          Rscript scripts/dbl/12_share_calibrate.R --calibrate --lane=<lane>|array [--quick] [--stub] [--force] [--cores=16]

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here); library(parallel); library(nnls) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))   # the dated copy that runs
ARGS <- commandArgs(TRUE); FORCE <- "--force" %in% ARGS
MODE <- c("reuse-smoke", "calibrate")[c("--reuse-smoke", "--calibrate") %in% ARGS]; stopifnot(length(MODE) == 1)   # exactly one mode
argv <- function(key, default) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), ARGS, value = TRUE)); if (length(v)) v else default }
CORES <- as.integer(argv("cores", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
if (MODE == "reuse-smoke") {   # stage 0b (frozen): the configs as tagged, fixed outputs
  MD5_FROZEN <- c(config_dbl = "91d7f54db724b7d87e9a82a4f3838e5e", config_paths = "b2c7ca9b7d65f2c39774e6ee2cca2523")   # configs as tagged (prereg-dbl-0b-v1)
  MD5 <- unname(tools::md5sum(c(here::here("scripts", "dbl", "config_dbl.R"), here::here("scripts", "config", "config_paths.R"))))
  if (!identical(MD5, unname(MD5_FROZEN))) stop("config md5 differs from the frozen literals: fix and resubmit (not a rerun)")
  OUT <- file.path(DBL_TAB_DIR, c(rates = "dbl_share_rates_smoke.csv", checks = "dbl_checks_12_smoke.csv")); names(OUT) <- c("rates", "checks")
  stopifnot(all(startsWith(OUT, DBL_SMOKE_DIR)))                     # outputs never inside the repo
  if (any(file.exists(OUT)) && !FORCE) stop("outputs exist (use --force)") }
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
## -- Stage 2: the T-partner flag DT's two typing routes (D with T in LEUK's place; partners Z = LEUK or P minus T) ----
PZ <- setdiff(DBL_SHARE_P, "T")
tx_origin <- function(origin) vapply(strsplit(ifelse(is.na(origin), "", origin), "+", fixed = TRUE), function(z)   # T + one Z
  length(z) == 2 && sum(z == "T") == 1 && all(z[z != "T"] %in% c(LIN, PZ)), TRUE)
nnls_tx <- function(X, B) {   # nnls_lx with T in LEUK's place
  Y <- as.matrix(X); Y <- sweep(Y, 2, pmax(colSums(Y), 1), "/")
  w <- 1 / sqrt(rowMeans(B$P) + 1e-6); Pw <- B$P * w; cn <- colnames(B$P)
  grp <- ifelse(cn %in% LIN, "LEUK", cn)
  vapply(seq_len(ncol(Y)), function(j) {
    b <- nnls::nnls(Pw, Y[, j] * w)$x; b <- b / max(sum(b), 1e-12)
    gw <- sort(tapply(b, grp, sum), decreasing = TRUE)
    length(gw) >= 2 && "T" %in% names(gw)[1:2] && all(setdiff(names(gw)[1:2], "T") %in% c("LEUK", PZ)) &&
      min(gw[1:2]) >= DBL_SHARE_D_MINOR && (gw[2] - (if (length(gw) > 2) gw[3] else 0)) >= DBL_SHARE_D_GAP }, TRUE) }

## == MODE --calibrate: Stage 2 step 2a, one lane (STAGE2_plan.md; user approval 2026-11-12) ====
if (MODE == "calibrate") {
  stopifnot(system2("git", c("-C", here::here(), "rev-parse", "-q", "--verify", "refs/tags/prereg-dbl-share-f1"), stdout = FALSE) == 0)
  QUICK <- "--quick" %in% ARGS; STUB <- "--stub" %in% ARGS
  SCDBL_MODE <- "devel"
  if (!STUB) { .libPaths(c(DBL_DEVEL_LIB, .libPaths()))             # devel 1.27.12 first, as stage 04 --scdbl=devel
    suppressPackageStartupMessages({ library(SingleCellExperiment); library(scDblFinder) })
    chk("scDblFinder version (devel)", "1.27.12", as.character(packageVersion("scDblFinder"))) }
  ## -- copied from 04_type_doublets.R (proc_320, call_scdbl; its MODE renamed SCDBL_MODE) ----
  proc_320 <- function(e, dims = NULL) {   # devel's default processing rebuilt with Bioc 3.20 functions (see 04's header)
    e <- as(e, "dgCMatrix"); if (is.null(dims)) dims <- 20
    stopifnot(ncol(e) <= 50000)              # devel skips normalisation above 50,000 barcodes (misc.R doNorm); not mirrored here
    sf <- Matrix::colSums(e); sf <- sf / mean(sf)
    le <- scuttle::normalizeCounts(e, size.factors = sf, center.size.factors = FALSE)
    pca <- scater::calculatePCA(le, ncomponents = dims, subset_row = seq_len(nrow(le)), ntop = nrow(le), BSPARAM = BiocSingular::IrlbaParam())
    if (is.list(pca)) pca <- pca$x
    row.names(pca) <- colnames(e); pca }
  call_scdbl <- function(Mrun, labs, known, dbr) {   # one scDblFinder call; returns per-barcode results + run facts
    sce  <- SingleCellExperiment(list(counts = Mrun))
    args <- list(sce, clusters = labs, dbr = dbr)
    if (any(known)) { args$knownDoublets <- known; args$knownUse <- "discard" }
    if (SCDBL_MODE == "main") args$adjustSize <- DBL_ADJUST_SIZE      # passed through ... to the artificial-doublet generator
    if (SCDBL_MODE == "devel") args$processing <- proc_320
    warns <- character(0)
    set.seed(DBL_SEED_CALL)
    t0  <- proc.time()[["elapsed"]]
    sce <- withCallingHandlers(do.call(scDblFinder::scDblFinder, args),
             warning = function(cond) { warns <<- c(warns, conditionMessage(cond)); invokeRestart("muffleWarning") },
             message = function(cond) invokeRestart("muffleMessage"))
    thr <- S4Vectors::metadata(sce)$scDblFinder.threshold
    list(res = data.table(cell = colnames(sce), called = sce$scDblFinder.class == "doublet", score = sce$scDblFinder.score,
                          origin = as.character(sce$scDblFinder.mostLikelyOrigin), ambiguous = as.logical(sce$scDblFinder.originAmbiguous)),
         secs = proc.time()[["elapsed"]] - t0, warns = warns, threshold = if (is.null(thr)) NA_real_ else as.numeric(thr)[1]) }
  call_stub <- function(Mrun, labs, known, dbr) {   # --stub: random calls in call_scdbl's format (harness and structural tests only)
    set.seed(DBL_SEED_CALL); n <- ncol(Mrun); u <- sort(unique(labs))
    list(res = data.table(cell = colnames(Mrun), called = runif(n) < 0.3, score = runif(n),
                          origin = vapply(seq_len(n), function(i) paste(sort(u[sample.int(length(u), 2L)]), collapse = "+"), ""),
                          ambiguous = runif(n) < 0.2), secs = 0, warns = character(0), threshold = 0.5) }
  CALL <- if (STUB) call_stub else call_scdbl
  ## -- [c1] lane, donors, weights ----
  LT <- fread(here::here("scripts", "dbl", DBL_SHARE_LANES_FILE), colClasses = list(character = c("patient", "lane")))[role == "analysis"]
  LT[, donor := fifelse(dataset == "GSE185381", patient, lane)]       # stage 01's donor id (as 00_share_eligibility.R)
  EL <- fread(file.path(DBL_TAB_DIR, "dbl_share_stage1_eligibility.csv"), colClasses = list(character = "patient"))[eligible == TRUE, patient]
  LANES_ALL <- unique(LT$lane)                                          # lane index = position in the committed lane table
  WORK <- LANES_ALL[LANES_ALL %in% LT[patient %in% EL, lane]]           # lanes of Stage-1-eligible patients, in table order
  L <- argv("lane", ""); if (identical(L, "array")) L <- WORK[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
  stopifnot(length(L) == 1, !is.na(L), L %in% WORK); LI <- match(L, LANES_ALL)
  NL  <- LT[, .(nl = uniqueN(lane)), by = patient]
  DON <- NL[LT[lane == L & patient %in% EL, .(patient, donor)], on = "patient"][, .(donor, w = fifelse(nl >= DBL_SHARE_HALF_LANES, 0.5, 1))]
  HASHED <- isTRUE(LT[lane == L, hashed[1]])
  OUTF <- file.path(DBL_SHARE_DIR, paste0(L, "__cal", if (QUICK) "_quick", if (STUB) "_stub", ".rds"))
  stopifnot(startsWith(OUTF, DBL_SMOKE_DIR)); if (file.exists(OUTF) && !FORCE) stop("output exists (use --force)")
  dir.create(DBL_SHARE_DIR, showWarnings = FALSE, recursive = TRUE)
  ts(sprintf("[c1] lane %s: index %d of %d (%d to calibrate); donors %s; hashed %s%s%s", L, LI, length(LANES_ALL), length(WORK),
             paste(sprintf("%s (w %.1f)", DON$donor, DON$w), collapse = ", "), HASHED, if (QUICK) "; QUICK" else "", if (STUB) "; STUB CALLER" else ""))
  x <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds"))); M <- x$counts
  lab <- fread(file.path(DBL_LANE_DIR, paste0(L, "__labels.csv.gz")), colClasses = list(character = "donor"))[match(colnames(M), cell)]
  stopifnot(identical(lab$cell, colnames(M)), !anyNA(lab$compartment))
  N <- ncol(M); n_dep <- nrow(x$bc); dbr <- min(DBL_DBR_CAP, n_dep * DBL_DBR_PER_BARCODE)   # M1's dbr rule (stage 04)
  A <- lab[class == "A", .(cell, donor, compartment, bmm_fine)]
  chk(sprintf("%s: every donor to calibrate has class-A T and LEUK cells in the lane", L), TRUE,
      all(vapply(DON$donor, function(d) A[donor == d, any(compartment == "T") && any(compartment %in% LIN)], TRUE)))
  set.seed(DBL_SEED_SHARE_CAL + 100L * LI)
  A[, fold := ((sample(.N) - 1L) %% DBL_N_FOLDS) + 1L, by = .(donor, compartment)]   # balanced folds within donor x compartment
  BIG <- names(which.max(table(lab$compartment)))   # stage 04: label of every injected known barcode (an existing cluster)
  ## ambient marker sets and LEUK share: stage 04 / --reuse-smoke verbatim
  pools <- split(A$cell, A$compartment)
  cpm <- function(cells) { s <- M[, cells, drop = FALSE]; Matrix::rowSums(s) / sum(s) * 1e6 }
  badg <- grepl("^MT-|^RP[SL]|^HB[^P]", rownames(M))
  cT <- cpm(pools[["T"]]); cL <- cpm(unlist(pools[intersect(STATES, names(pools))]))
  GT4 <- rownames(M)[!badg & cT > 5 & cL < 1]; GL4 <- rownames(M)[!badg & cL > 5 & cT < 1]
  leuk_share <- function(Mx) { l <- Matrix::colSums(Mx[GL4, , drop = FALSE]); t <- Matrix::colSums(Mx[GT4, , drop = FALSE]); ifelse(l + t > 0, l / (l + t), NA_real_) }
  chk(sprintf("%s: >= 10 genes in each ambient marker set (T %d, LEUK %d)", L, length(GT4), length(GL4)), TRUE, length(GT4) >= 10 && length(GL4) >= 10)
  pL <- Matrix::rowSums(M[, A[compartment %in% LIN, cell], drop = FALSE]); gpos <- which(pL > 0); pL <- pL[gpos] / sum(pL)   # LEUK pseudo-bulk
  PLT <- A[bmm_fine %in% DBL_SHARE_PS_PLT, cell]; USE_PLT <- length(PLT) >= DBL_SHARE_PS_PLT_MIN
  ts(sprintf("    %d analysis barcodes (%d deposited, dbr %.4f), %d class A; %d '%s' cells -> T+platelet %s", N, n_dep, dbr, nrow(A),
             length(PLT), DBL_SHARE_PS_PLT, if (USE_PLT) "used" else "SKIPPED"))
  ## -- [c2] the plan: items per donor (planned n_full, n after availability), runs, round-robin split ----
  getp <- function(src, cp) if (cp == "LEUK") unlist(src[intersect(LIN, names(src))], use.names = FALSE) else
    if (cp == "any") unlist(src, use.names = FALSE) else if (is.null(src[[cp]])) character(0) else src[[cp]]
  NNP <- rbind(data.table(x1 = PZ, x2 = PZ), as.data.table(t(combn(PZ, 2)))[, .(x1 = V1, x2 = V2)])   # 21 pairs within P minus T
  PLAN <- rbindlist(lapply(seq_len(nrow(DON)), function(i) { d <- DON$donor[i]; al <- round(DON$w[i] * DBL_SHARE_ALLOC)
    srcA <- split(A[donor == d, cell], A[donor == d, compartment]); nA <- function(cp) length(getp(srcA, cp))
    nn <- NNP[mapply(function(a, b) nA(a) > 0 && nA(b) > (a == b), x1, x2)]          # pairs present in the donor's class A
    p <- rbind(
      data.table(grp = "T_LEUK", type = "T|LEUK", ca = "T", cb = LIN, n_full = al[["T_LEUK"]]),
      data.table(grp = "X_LEUK", type = paste0(PZ, "|LEUK"), ca = PZ, cb = "LEUK", n_full = al[["X_LEUK"]]),
      data.table(grp = "LEUK_LEUK", type = "LEUK|LEUK", ca = "LEUK", cb = "LEUK", n_full = al[["LEUK_LEUK"]]),
      data.table(grp = "T_X", type = paste0("T|", DBL_SHARE_P[c(2:7, 1)]), ca = "T", cb = DBL_SHARE_P[c(2:7, 1)], n_full = al[["T_X"]]),
      if (nrow(nn)) data.table(grp = "NN", type = paste0(nn$x1, "|", nn$x2), ca = nn$x1, cb = nn$x2,
                               n_full = diff(round(seq(0, al[["NN"]], length.out = nrow(nn) + 1)))),   # equal split, sums exactly
      data.table(grp = "TRIPLET", type = "triplet", ca = "any", cb = "any", n_full = al[["TRIPLET"]]),
      if (HASHED) data.table(grp = "CROSS", type = "cross", ca = "any", cb = "other donor", n_full = al[["CROSS"]]),
      data.table(grp = "PS_AMB", type = "ps_amb", ca = "T", cb = "ambient", n_full = al[["PS_AMB"]]),
      data.table(grp = "PS_PLT", type = "ps_plt", ca = "T", cb = DBL_SHARE_PS_PLT, n_full = al[["PS_PLT"]]))
    p[, ok := mapply(function(g, a, b) switch(g, TRIPLET = nA("any") >= 3, CROSS = A[donor != d, .N] > 0, PS_AMB = nA("T") > 0,
                                              PS_PLT = USE_PLT && nA("T") > 0, nA(a) > 0 && nA(b) > (a == b)), grp, ca, cb)]
    p[, `:=`(donor = d, n = fifelse(ok, n_full, 0))] }))
  PLAN[, item := .I]
  I_TOT <- sum(PLAN$n); CAP <- floor(DBL_SHARE_INJ_MAX * N)
  NB <- max(DBL_SHARE_MIN_BATCH, ceiling(I_TOT / (DBL_N_FOLDS * CAP))); R <- DBL_N_FOLDS * NB
  u <- rep(PLAN$item, PLAN$n)                                           # one entry per planned synthetic barcode
  K <- unclass(table(factor(u, PLAN$item), factor((seq_along(u) - 1L) %% R + 1L, seq_len(R))))   # item x run: round robin
  RUNS <- if (QUICK) (seq_len(DBL_N_FOLDS) - 1L) * NB + 1L else seq_len(R)   # quick: batch 1 of every fold
  print(PLAN[, .(n_full = sum(n_full), n = sum(n), items_absent = sum(!ok)), by = .(donor, grp)])
  ts(sprintf("[c2] %d synthetic barcodes planned; cap %d per run -> %d batches x %d folds = %d runs (max %d per run); running %d%s",
             I_TOT, CAP, NB, DBL_N_FOLDS, R, max(colSums(K)), length(RUNS), if (LI * 100L + R >= (LI + 1L) * 100L) "; seeds overlap later lanes" else ""))
  chk(sprintf("%s: planned injections per run <= the cap (%d)", L, CAP), TRUE, max(colSums(K)) <= CAP)
  ## -- [c3] one held-out run ----
  build_run <- function(r) {   # synthetic barcodes of run r: rows (parents) + counts; parents are held-out class-A cells of fold f
    f <- (r - 1L) %/% NB + 1L; H <- A[fold == f]
    set.seed(DBL_SEED_SHARE_CAL + 100L * LI + r)
    rows <- list()
    for (i in which(K[, r] > 0)) { p <- PLAN[i]; k <- K[i, r]
      src <- split(H[donor == p$donor, cell], H[donor == p$donor, compartment])
      a <- b <- cc <- NA_character_
      if (p$grp == "TRIPLET") { pa <- getp(src, "any"); if (length(pa) < 3) next
        tr <- matrix(vapply(seq_len(k), function(j) pa[sample.int(length(pa), 3L)], character(3)), nrow = 3)   # 3 distinct cells
        a <- tr[1, ]; b <- tr[2, ]; cc <- tr[3, ]
      } else if (p$grp == "PS_AMB") { pa <- getp(src, "T"); if (!length(pa)) next; a <- sample_cells(pa, k)
      } else {
        pa <- getp(src, p$ca)
        pb <- switch(p$grp, CROSS = H[donor != p$donor, cell], PS_PLT = H[cell %in% PLT, cell], getp(src, p$cb))
        if (p$grp == "CROSS") pa <- getp(src, "any")
        if (!length(pa) || !length(pb) || (p$ca == p$cb && length(pa) < 2)) next   # homotypic needs 2 cells in the fold
        a <- sample_cells(pa, k); b <- sample_cells(pb, k)
        for (j in which(a == b)) { q <- setdiff(pb, a[j]); if (length(q)) b[j] <- q[sample.int(length(q), 1L)] }   # never a cell with itself
      }
      rows[[length(rows) + 1L]] <- data.table(item = i, donor = p$donor, grp = p$grp, type = p$type, a = a, b = b, p3 = cc) }
    inj <- rbindlist(rows); inj[, syn := sprintf("SYN_%d_%05d", r, .I)]
    Mi <- M[, inj$a, drop = FALSE]
    ib <- which(!is.na(inj$b)); if (length(ib)) Mi[, ib] <- Mi[, ib] + M[, inj$b[ib], drop = FALSE]
    ic <- which(!is.na(inj$p3)); if (length(ic)) Mi[, ic] <- Mi[, ic] + M[, inj$p3[ic], drop = FALSE]
    ia <- which(inj$grp == "PS_AMB")
    if (length(ia)) { sz <- round(DBL_SHARE_PS_AMB * Matrix::colSums(M[, inj$a[ia], drop = FALSE]))
      AMB <- vapply(sz, function(s) { v <- numeric(nrow(M)); v[gpos] <- rmultinom(1L, s, pL)[, 1]; v }, numeric(nrow(M)))
      Mi[, ia] <- Mi[, ia] + Matrix::Matrix(AMB, sparse = TRUE) }
    colnames(Mi) <- inj$syn
    list(inj = inj, Mi = Mi, held = H$cell, f = f) }
  run_one <- function(r) {
    bi <- build_run(r); inj <- bi$inj; Mi <- bi$Mi
    keep <- setdiff(colnames(M), bi$held)
    Mrun <- cbind(M[, keep, drop = FALSE], Mi)
    labs <- c(lab$compartment[match(keep, lab$cell)], rep(BIG, ncol(Mi)))
    cr <- CALL(Mrun, labs, c(rep(FALSE, length(keep)), rep(TRUE, ncol(Mi))), dbr)
    out <- merge(inj, cr$res[cell %in% inj$syn], by.x = "syn", by.y = "cell")   # MASKING: synthetic barcodes only from here on
    BS <- make_basis(M, A[!(cell %in% bi$held), .(cell, compartment)])
    Xg <- Mi[BS$G, out$syn, drop = FALSE]
    out[, `:=`(nnls = nnls_lin(Xg, BS), nnls_lx = nnls_lx(Xg, BS), nnls_tx = nnls_tx(Xg, BS))]
    out[, share := leuk_share(Mi[, out$syn, drop = FALSE])]
    out[, tid := Matrix::colSums(Mi[intersect(DBL_T_ID_GENES, rownames(Mi)), out$syn, drop = FALSE]) >= DBL_T_ID_MIN_UMI]
    out[, lin_o := lin_origin(origin)]                                                 # stage 10 verbatim from here ...
    out[, `:=`(V1a = called %in% TRUE & ambiguous %in% FALSE & lin_o, V3 = called %in% TRUE & nnls)]
    out[, V5 := V1a | V3]
    out[, M1 := V5 & tid]                                                              # ... to here
    out[, D  := M1 | (called %in% TRUE & ((ambiguous %in% FALSE & lx_origin(origin)) | nnls_lx))]   # as --reuse-smoke
    out[, DT := M1 | (called %in% TRUE & ((ambiguous %in% FALSE & tx_origin(origin)) | nnls_tx))]   # T-partner flag
    par_ok <- all(na.omit(c(inj$a, inj$b, inj$p3)) %in% bi$held)                       # every parent is a held-out class-A cell
    list(out = out[, run := r], par_ok = par_ok,
         run = data.table(run = r, fold = bi$f, n_real = length(keep), n_injected = ncol(Mi), n_planned = sum(K[, r]), secs = cr$secs,
                          threshold = cr$threshold, n_warn = length(cr$warns), first_warn = substr(c(cr$warns, "")[1], 1, 150))) }
  real_one <- function() {   # MASKED real run: one aggregate row leaves this function, nothing per barcode
    cr <- CALL(M, lab$compartment, rep(FALSE, N), dbr)
    ok <- nrow(cr$res) == N && identical(cr$res$cell, colnames(M)) && all(cr$res$called %in% c(TRUE, FALSE))
    cl28 <- intersect(x$ann[!is.na(atlas_cl) & atlas_cl == DBL_ATLAS_HOMOTYPIC_CL, cell], colnames(M))
    agg <- data.table(threshold = cr$threshold, called_frac = mean(cr$res$called), n_cl28 = length(cl28),
                      called_cl28 = if (length(cl28)) mean(cr$res$called[match(cl28, cr$res$cell)]) else NA_real_,
                      secs = cr$secs, n_warn = length(cr$warns), ok_format = ok)
    rm(cr); list(real = agg) }
  ## -- [c4] all runs + the masked real run, forked ----
  JOBS <- c(0L, RUNS)                                                   # 0 = the masked real run
  ts(sprintf("[c4] %d scDblFinder calls on %d cores", length(JOBS), min(CORES, length(JOBS))))
  RES <- mclapply(JOBS, function(j) tryCatch(if (j == 0L) real_one() else run_one(j), error = function(e) list(error = conditionMessage(e))),
                  mc.cores = min(CORES, length(JOBS)), mc.preschedule = FALSE)
  bad <- vapply(RES, function(z) !is.list(z) || !is.null(z$error), TRUE)   # a killed fork returns NULL or a try-error
  for (j in which(bad)) ts("    ** call ", JOBS[j], " failed: ", if (is.list(RES[[j]])) RES[[j]]$error else paste(as.character(RES[[j]]), collapse = " "))
  chk(sprintf("%s: every scDblFinder call finished", L), 0L, sum(bad))
  REAL <- if (!bad[1]) RES[[1]]$real else data.table()
  OK <- RES[-1][!bad[-1]]
  RUNTAB <- rbindlist(lapply(OK, `[[`, "run")); SYN <- rbindlist(lapply(OK, `[[`, "out"))
  print(RUNTAB[, .(run, fold, n_real, n_injected, n_planned, secs = round(secs), threshold = round(threshold, 3), n_warn)], nrows = 12)
  print(REAL)
  ## -- [c5] window, labels of the parents, checks ----
  WIN <- SYN[type == "T|LEUK", quantile(share, DBL_SHARE_WIN_Q, na.rm = TRUE, type = 7)]
  SYN[, Y_w := M1 & !is.na(share) & share >= WIN[1] & share <= WIN[2]]
  ci <- function(v) A$compartment[match(v, A$cell)]
  SYN[, `:=`(la = ci(a), lb = ci(b), lc = ci(p3), donor_b = A$donor[match(b, A$cell)])]
  pk <- PLAN[SYN$item]                                                  # the planned item of every synthetic barcode
  is_c <- function(lbl, cp) { cp <- rep_len(cp, length(lbl)); fifelse(cp == "LEUK", lbl %in% LIN, fifelse(cp == "any", !is.na(lbl), lbl == cp)) }
  id_ok <- SYN[, fifelse(grp == "PS_AMB", la == "T" & is.na(b), fifelse(grp == "PS_PLT", la == "T" & b %in% PLT,
                 fifelse(grp == "CROSS", !is.na(la) & !is.na(lb) & donor_b != donor,
                 is_c(la, pk$ca) & is_c(lb, pk$cb) & (grp != "TRIPLET" | is_c(lc, "any")))))]
  within <- SYN[!grp %in% c("CROSS", "PS_PLT"), all(A$donor[match(na.omit(c(a, b, p3)), A$cell)] == rep(donor, 3)[!is.na(c(a, b, p3))])]
  if (!STUB) chk(sprintf("%s: devel trains its classifier in every run (no xgboost warning)", L), FALSE,
                 any(RUNTAB[, grepl("invalid function arguments|have been removed", first_warn)]))
  chk(sprintf("%s: injected per run <= cap and <= planned", L), TRUE, RUNTAB[, all(n_injected <= CAP & n_injected <= n_planned)])
  chk(sprintf("%s: every output barcode is synthetic; one row per injected barcode", L), TRUE,
      all(startsWith(SYN$syn, "SYN_")) && nrow(SYN) == sum(RUNTAB$n_injected) && !anyDuplicated(SYN[, paste(run, syn)]))
  chk(sprintf("%s: every parent is a held-out class-A cell of its run's fold", L), TRUE, all(vapply(OK, `[[`, TRUE, "par_ok")))
  chk(sprintf("%s: injection identity (parent labels match the pair type)", L), TRUE, all(id_ok))
  chk(sprintf("%s: within-donor types have parents of one donor; no cell paired with itself", L), TRUE,
      within && SYN[!is.na(b), all(a != b)] && SYN[!is.na(p3), all(p3 != a & p3 != b)])
  chk(sprintf("%s: M1 is a subset of D and of DT for every synthetic barcode", L), TRUE, SYN[, all(!M1 | (D & DT))])
  chk(sprintf("%s: the masked real run returned one aggregate row and nothing else", L), TRUE,
      nrow(REAL) == 1 && identical(names(REAL), c("threshold", "called_frac", "n_cl28", "called_cl28", "secs", "n_warn", "ok_format")) && isTRUE(REAL$ok_format))
  PLAN[, realised := SYN[, .N, by = item][PLAN, on = "item", fifelse(is.na(N), 0L, N)]]
  print(SYN[, .(n = .N, called = round(mean(called), 3), M1 = round(mean(M1), 3), D = round(mean(D), 3), DT = round(mean(DT), 3),
                Y_w = round(mean(Y_w), 3)), by = .(grp, type)][order(grp, type)], nrows = 60)
  ts(sprintf("    window (T|LEUK LEUK share %g-%g): %.3f-%.3f", DBL_SHARE_WIN_Q[1], DBL_SHARE_WIN_Q[2], WIN[1], WIN[2]))
  CK <- rbindlist(CHK)
  saveRDS(list(lane = L, lane_index = LI, dataset = LT[lane == L, dataset[1]], chemistry = LT[lane == L, chemistry[1]], hashed = HASHED,
               quick = QUICK, stub = STUB, scdbl = if (STUB) "stub" else as.character(packageVersion("scDblFinder")),
               md5 = c(self = unname(tools::md5sum(SELF)), config_dbl = unname(tools::md5sum(here::here("scripts", "dbl", "config_dbl.R")))),
               dbr = dbr, N = N, n_dep = n_dep, cap = CAP, batches = NB, runs_planned = R, runs_done = RUNS, donors = DON, plan = PLAN,
               runs = RUNTAB, real = REAL, window = WIN, n_plt = length(PLT), markers = c(T = length(GT4), LEUK = length(GL4)),
               comp_A = A[, .N, by = .(donor, compartment)],
               syn = SYN[, .(run, syn, item, donor, grp, type, a, b, p3, la, lb, lc, donor_b, called, score, origin, ambiguous,
                             nnls, nnls_lx, nnls_tx, tid, share, M1, D, DT, Y_w)], checks = CK), OUTF)
  ts(sprintf("[done] %s; %d checks, %d FAIL", OUTF, nrow(CK), sum(!CK$pass)))
  quit(status = if (all(CK$pass)) 0L else 1L) }

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
