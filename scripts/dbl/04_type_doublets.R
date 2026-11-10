#!/usr/bin/env Rscript
# 04_type_doublets.R ----
# Doublet-pairing SMOKE, stage 4: scDblFinder calling + origin typing at leukemic-state resolution, per lane, measured ONLY on
# barcodes whose truth we built (PROPOSAL_doublet_pairing_external_v1.3 section 4.1; appendix C; handoff v1.1 4.7):
#   held-out runs  one fifth of class A is removed from the lane; synthetic doublets (counts of two held-out cells summed) and
#                  held-out singlets are injected as knownDoublets with knownUse = "discard" (never parents of artificial
#                  doublets, called by the threshold like any barcode) -> confusion matrix C (retention x reported state),
#                  false-positive rates of other pair types and singlets, per-state recall / precision
#   lodo run       (pooled lane) the same with one donor's class-A cells held out = the proposal's held-out patient
#   planted runs   T + LEUK_GMP contact doublets of one donor, NOT flagged, at 0.5 / 1 / 2 % of the lane -> recovery curve
#   real run       the lane as it is. MASKED [D2]: the worker returns aggregates only -- the called fraction of atlas cluster-28
#                  barcodes (homotypic doublets, positive control) and of hashtag cross-donor / same-donor-two-tag barcodes
#                  (stage 03), plus format assertions and run time.
#   Re-entered held-out singlets are REAL cells: since 2026-11-10 (review M-mask F1, user decision) only "reported as any
#   T-LEUK" leaves the worker for them; their origin, score and typed state are blanked inside run_one(). Runs before that
#   date (jobs 3687233/3687239 main, 3687240 devel) wrote singlet states; in devel 2 class-A barcodes of GSM5613785 were typed
#   T-LEUK. Those files are kept and registered as a known exposure (SMOKE_REPORT section 8), never read by the analyst.
# INPUT  : DBL_LANE_DIR/<lane>.rds (01), <lane>__labels.csv.gz (02), <pooled lane>__hto_calls.csv.gz (03)
# OUTPUT : DBL_TAB_DIR/dbl_typing_confusion<SFX>.csv   lane x run kind x true state x reported label: injected T-LEUK doublets
#          DBL_TAB_DIR/dbl_typing_summary<SFX>.csv     lane x state: n, retained, recall, recall given retention, precisions
#          DBL_TAB_DIR/dbl_typing_fp<SFX>.csv          lane x false-positive source: n injected, reported as T-LEUK (any / state)
#          DBL_TAB_DIR/dbl_typing_expected<SFX>.csv    lane x donor x state: model expectation under random co-encapsulation
#                                                      (signal, false positives) from singlet composition + synthetic rates
#          DBL_TAB_DIR/dbl_typing_planted<SFX>.csv     lane x rate: planted, recovered as T-LEUK_GMP / any T-LEUK / called
#          DBL_TAB_DIR/dbl_typing_controls<SFX>.csv    lane: cluster-28 and hashtag controls (aggregate fractions only)
#          DBL_TAB_DIR/dbl_typing_ambient<SFX>.csv     lane x state: LEUK-marker share quantiles of injected T-LEUK doublets
#          DBL_TAB_DIR/dbl_typing_runs<SFX>.csv        one row per scDblFinder call: kind, barcodes, injected, seconds, warnings
#          DBL_TAB_DIR/dbl_typing_losses<SFX>.csv      lane x run kind x true state: where injected T-LEUK doublets are lost
#          DBL_TAB_DIR/dbl_typing_injected<SFX>.csv.gz one row per INJECTED barcode: synthetic doublets with call, origin, share;
#                                                      re-entered singlets with called / reported-as-any-T-LEUK only (masked)
#          DBL_TAB_DIR/dbl_checks_04<SFX>.csv
#          DBL_LANE_DIR/<lane>__calib<SFX>.rds          synthetic-derived calibration + singlet composition, read by stage 05
#          (SFX = "_<mode>[_<tag>]_smoke"; --tag keeps a verification rerun from overwriting earlier files)
# Usage  : Rscript scripts/dbl/04_type_doublets.R [--scdbl=main|devel] [--lanes=single,pooled] [--cores=8] [--tag=name]
#          main  = installed scDblFinder 1.20.2 with adjustSize = 0 (issue #145 workaround) [D4]. Its xgboost classifier never
#                  trains under xgboost 3.2.1.1 (1.20.2 swallows the error), so every 'main' call is the pre-training
#                  heuristic score (cxds + kNN ratio); recorded as an explicit check, not a silent fallback (review R-04 F1)
#          devel = scDblFinder 1.27.12 (commit 1a0831bd) from DBL_DEVEL_LIB, its own defaults [D4], except one call-site
#                  argument: its default processing hands a SparseArray to Bioc 3.20's scater::calculatePCA, which this
#                  SparseArray (1.6.2) cannot divide ("arithmetic ... SparseArray object and a vector of length != 1");
#                  processing = proc_320() rebuilds the same algorithm (centred library-size factors, log2(x/sf + 1),
#                  irlba PCA on all selected features) with scuttle instead of scrapper. Checked on mockDoubletSCE
#                  2026-11-09: runs, trains its classifier, no warnings. No new package.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here); library(parallel) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
ARGS  <- commandArgs(TRUE)
argv  <- function(key, default) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), ARGS, value = TRUE)); if (length(v)) v else default }
MODE  <- argv("scdbl", "main"); stopifnot(MODE %in% c("main", "devel"))
CORES <- as.integer(argv("cores", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
la    <- argv("lanes", ""); LANES <- if (nzchar(la)) DBL_SMOKE_LANES[strsplit(la, ",")[[1]]] else DBL_SMOKE_LANES
stopifnot(!anyNA(LANES))
if (MODE == "devel") .libPaths(c(DBL_DEVEL_LIB, .libPaths()))   # before ANY package that the devel library also holds
suppressPackageStartupMessages({ library(SingleCellExperiment); library(scDblFinder) })
TAG <- argv("tag", "")
SFX <- paste0("_", MODE, if (nzchar(TAG)) paste0("_", TAG), "_smoke")
ts  <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-78s expect %-12s got %-12s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
VER <- c(scDblFinder = as.character(packageVersion("scDblFinder")), xgboost = as.character(packageVersion("xgboost")))
ts(sprintf("[0] mode %s: scDblFinder %s from %s; xgboost %s; %d cores", MODE, VER[1], find.package("scDblFinder"), VER[2], CORES))
chk("scDblFinder version of this mode", if (MODE == "main") "1.20.2" else "1.27.12", VER[["scDblFinder"]])
chk("scDblFinder accepts knownDoublets/knownUse", TRUE, all(c("knownDoublets", "knownUse") %in% names(formals(scDblFinder::scDblFinder))))
STATES <- DBL_LEUK_STATES
WIN <- DBL_AMBIENT_WINDOW

## -- helpers ----
parse_state <- function(origin) {   # "LEUK_GMP+T" -> "LEUK_GMP" when the origin is exactly {T, one LEUK state}; else NA
  p <- strsplit(ifelse(is.na(origin), "", origin), "+", fixed = TRUE)
  vapply(p, function(z) { s <- z[z %in% STATES]; if (length(z) == 2 && "T" %in% z && length(s) == 1) s else NA_character_ }, "") }
proc_320 <- function(e, dims = NULL) {   # devel's default processing rebuilt with Bioc 3.20 functions (see header)
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
  if (MODE == "main") args$adjustSize <- DBL_ADJUST_SIZE      # passed through ... to the artificial-doublet generator
  if (MODE == "devel") args$processing <- proc_320
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

ALL <- list(conf = list(), summ = list(), fp = list(), expd = list(), plant = list(), ctrl = list(), amb = list(), runs = list(),
            lost = list(), inj = list())
for (L in LANES) {
  ts("[1] lane ", L)
  x   <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds")))
  M   <- x$counts
  lab <- fread(file.path(DBL_LANE_DIR, paste0(L, "__labels.csv.gz")))[match(colnames(M), cell)]
  stopifnot(identical(lab$cell, colnames(M)), !anyNA(lab$compartment))
  N     <- ncol(M)
  n_dep <- nrow(x$bc)
  dbr   <- min(DBL_DBR_CAP, n_dep * DBL_DBR_PER_BARCODE)
  ts(sprintf("    %d analysis barcodes (of %d deposited); dbr %.4f", N, n_dep, dbr))
  ## -- [2] singlet pools (class A), ambient marker sets, folds ----
  A <- lab[class == "A"]
  pools <- split(A$cell, A$compartment)
  ts("    class-A pools: ", paste(sprintf("%s %d", names(pools), lengths(pools)), collapse = ", "))
  cpm <- function(cells) { s <- M[, cells, drop = FALSE]; Matrix::rowSums(s) / sum(s) * 1e6 }   # pooled CPM of a cell set
  bad <- grepl("^MT-|^RP[SL]|^HB[^P]", rownames(M))
  cT  <- cpm(pools[["T"]]); cL <- cpm(unlist(pools[intersect(STATES, names(pools))]))
  GT  <- rownames(M)[!bad & cT > 5 & cL < 1]; GL <- rownames(M)[!bad & cL > 5 & cT < 1]   # proposal D: CPM rule
  ts(sprintf("    ambient marker sets: %d T-specific, %d LEUK-specific genes", length(GT), length(GL)))
  chk(sprintf("%s: >= 10 genes in each ambient marker set", L), TRUE, length(GT) >= 10 && length(GL) >= 10)
  leuk_share <- function(Mx) { l <- Matrix::colSums(Mx[GL, , drop = FALSE]); t <- Matrix::colSums(Mx[GT, , drop = FALSE]); ifelse(l + t > 0, l / (l + t), NA_real_) }
  set.seed(DBL_SEED_SYN)
  A[, fold := ((sample(.N) - 1L) %% DBL_N_FOLDS) + 1L, by = compartment]   # random, balanced folds within each compartment
  donors <- sort(unique(A$donor))
  BIG <- names(which.max(table(lab$compartment)))   # label of every injected KNOWN barcode: inert (never a parent), and an
                                                    # existing cluster, so no cluster is made only of discarded barcodes
  ## -- [3] run definitions ----
  RUNS <- list(list(kind = "real", id = "real"))
  for (f in DBL_FOLDS_RUN) for (b in seq_len(DBL_SYN_BATCHES)) RUNS[[length(RUNS) + 1]] <- list(kind = "heldout", id = sprintf("f%d_b%d", f, b), fold = f, batch = b)
  if (length(donors) > 1) {   # leave-one-donor-out: the donor with the most T x LEUK material
    lodo <- A[, .(nT = sum(compartment == "T"), nL = sum(compartment %in% STATES)), by = donor][order(-pmin(nT, nL))][1, donor]
    RUNS[[length(RUNS) + 1]] <- list(kind = "lodo", id = paste0("lodo_", lodo), donor = lodo) }
  for (r in seq_along(DBL_INJECT_RATES)) RUNS[[length(RUNS) + 1]] <- list(kind = "planted", id = sprintf("plant_%g", DBL_INJECT_RATES[r]), rate = DBL_INJECT_RATES[r], ri = r)
  HTO <- if (L == DBL_SMOKE_LANES[["pooled"]]) fread(file.path(DBL_LANE_DIR, paste0(L, "__hto_calls.csv.gz")))[in_counts == TRUE, .(cell, call)] else NULL
  CL28 <- x$ann[atlas_cl == DBL_ATLAS_HOMOTYPIC_CL, cell]
  ## -- [4] one run ----
  sample_cells <- function(pool, n) if (length(pool) >= n) sample(pool, n) else sample(pool, n, replace = TRUE)
  build_pairs <- function(src, n, typ, comp_a, comp_b, st = NA_character_) {   # n doublets: first parent from comp_a, second from comp_b
    pa <- src[[comp_a]]; pb <- src[[comp_b]]
    if (length(pa) == 0 || length(pb) == 0 || n == 0) return(NULL)
    data.table(type = typ, true_state = st, a = sample_cells(pa, n), b = sample_cells(pb, n), first = comp_a) }
  build_leuk_pairs <- function(src, n, typ, comp_a) {   # partner = a LEUK state drawn by availability, then a cell of it
    av <- lengths(src[intersect(STATES, names(src))]); av <- av[av > 0]
    if (!length(av) || length(src[[comp_a]]) == 0 || n == 0) return(NULL)
    st <- sample(names(av), n, replace = TRUE, prob = av)
    data.table(type = typ, true_state = NA_character_, a = sample_cells(src[[comp_a]], n), b = vapply(st, function(s) sample_cells(src[[s]], 1), ""), first = comp_a) }
  run_one <- function(R) {
    t_build <- proc.time()[["elapsed"]]
    if (R$kind == "real") {
      cr <- call_scdbl(M, lab$compartment, rep(FALSE, N), dbr)
      res <- cr$res
      st  <- parse_state(res$origin)
      ok_format <- nrow(res) == N && all(res$cell == colnames(M)) && all(res$called %in% c(TRUE, FALSE)) &&
                   all(is.na(st) | st %in% STATES)
      ctrl <- data.table(n_cl28 = length(CL28), called_cl28 = if (length(CL28)) mean(res[cell %in% CL28, called]) else NA_real_)
      if (!is.null(HTO)) { h <- merge(res[, .(cell, called)], HTO, by = "cell")
        ctrl[, `:=`(n_cross_donor = h[call == "cross_donor", .N], called_cross_donor = h[call == "cross_donor", mean(called)],
                    n_same_two_tags = h[call == "same_donor_two_tags", .N], called_same_two_tags = h[call == "same_donor_two_tags", mean(called)],
                    n_single_A = h[call == "single" & cell %in% A$cell, .N], called_single_A = h[call == "single" & cell %in% A$cell, mean(called)])] }
      run <- data.table(id = "real", kind = "real", n_real = N, n_injected = 0L, secs = cr$secs, threshold = cr$threshold,
                        n_warn = length(cr$warns), first_warn = substr(c(cr$warns, "")[1], 1, 150))
      rm(cr, res, st)   # MASKED [D2]: the return list below is the guarantee -- it holds no per-barcode object
      return(list(kind = "real", id = "real", ok_format = ok_format, ctrl = ctrl, run = run)) }
    if (R$kind %in% c("heldout", "lodo")) {
      set.seed(DBL_SEED_SYN + if (R$kind == "lodo") 999L else R$fold * 10L + R$batch)
      held <- if (R$kind == "lodo") A[donor == R$donor, cell] else A[fold == R$fold, cell]
      src  <- split(A[cell %in% held, cell], A[cell %in% held, compartment])
      prs  <- rbindlist(c(
        lapply(STATES, function(s) build_pairs(src, DBL_SYN_PER_STATE, "T|LEUK", "T", s, s)),
        lapply(c("Tamb", "NK", "B", "Ery", "MonoDC", "LEUK_other", "T"), function(cp) build_pairs(src, DBL_SYN_PER_OTHER, paste0("T|", cp), "T", cp)),
        lapply(c("Tamb", "NK", "B", "Ery", "MonoDC"), function(cp) build_leuk_pairs(src, DBL_SYN_PER_OTHER, paste0(cp, "|LEUK*"), cp)),
        list(build_leuk_pairs(src[intersect(names(src), c(STATES, "LEUK_other"))] |> (\(s) { s[["LEUKx"]] <- unlist(s[intersect(STATES, names(s))]); s })(),
                              DBL_SYN_PER_OTHER, "LEUK*|LEUK*", "LEUKx"))), fill = TRUE)
      comps <- names(src)[lengths(src) > 0]
      per   <- max(1L, floor(DBL_SYN_SINGLETS / length(comps)))
      sng   <- rbindlist(lapply(comps, function(cp) data.table(type = paste0("singlet:", cp), true_state = NA_character_,
                                 a = if (length(src[[cp]]) >= per) sample(src[[cp]], per) else src[[cp]], b = NA_character_, first = cp)))
      inj <- rbind(prs, sng, fill = TRUE)
      inj[, syn := sprintf("SYN_%s_%04d", R$id, .I)]
      Mi <- M[, inj$a, drop = FALSE]
      dbl <- which(!is.na(inj$b)); Mi[, dbl] <- Mi[, dbl] + M[, inj$b[dbl], drop = FALSE]   # synthetic doublet = summed counts
      colnames(Mi) <- inj$syn
      keep <- setdiff(colnames(M), held)
      Mrun <- cbind(M[, keep, drop = FALSE], Mi)
      labs <- c(lab$compartment[match(keep, lab$cell)], rep(BIG, nrow(inj)))
      cr <- call_scdbl(Mrun, labs, c(rep(FALSE, length(keep)), rep(TRUE, nrow(inj))), dbr)
      out <- cr$res[cell %in% inj$syn]                                # injected barcodes only; real ones are dropped here
      out <- merge(inj[, .(syn, type, true_state)], out, by.x = "syn", by.y = "cell")
      out[, share := leuk_share(Mi[, syn, drop = FALSE])]
      out[, st_raw := parse_state(origin)]
      out[, reported := fifelse(called & !is.na(ambiguous) & !ambiguous & !is.na(st_raw), st_raw, "none")]
      out[, reported_win := fifelse(reported != "none" & !is.na(share) & share >= WIN[1] & share <= WIN[2], reported, "none")]
      out[, tleuk_any := reported_win != "none"]
      ## MASKED [D2; user decision 2026-11-10 after review finding M-mask F1]: re-entered singlets are REAL held-out class-A
      ## cells, so only "reported as any T-LEUK" leaves this worker -- never their origin, score or typed state
      sg <- startsWith(out$type, "singlet:")
      out[sg, `:=`(score = NA_real_, origin = NA_character_, ambiguous = NA, st_raw = NA_character_,
                   reported = NA_character_, reported_win = NA_character_)]
      return(list(kind = R$kind, id = R$id, out = out, pools_held = lengths(src),
                  run = data.table(id = R$id, kind = R$kind, n_real = length(keep), n_injected = nrow(inj), secs = cr$secs,
                                   threshold = cr$threshold, n_warn = length(cr$warns), first_warn = substr(c(cr$warns, "")[1], 1, 150)))) }
    if (R$kind == "planted") {
      set.seed(DBL_SEED_PLANT + R$ri)
      n  <- round(R$rate * N)
      av <- A[, .(nT = sum(compartment == "T"), nG = sum(compartment == DBL_PLANT_STATE)), by = donor][pmin(nT, nG) > 0]
      dd <- sample(av$donor, n, replace = TRUE, prob = pmin(av$nT, av$nG))          # same-donor contacts
      pr <- rbindlist(lapply(split(seq_len(n), dd), function(ix) { d <- dd[ix[1]]
        data.table(a = sample_cells(A[donor == d & compartment == "T", cell], length(ix)),
                   b = sample_cells(A[donor == d & compartment == DBL_PLANT_STATE, cell], length(ix))) }))
      pr[, syn := sprintf("PLANT_%s_%04d", R$id, .I)]
      Mi <- M[, pr$a, drop = FALSE] + M[, pr$b, drop = FALSE]; colnames(Mi) <- pr$syn
      keep <- setdiff(colnames(M), unique(c(pr$a, pr$b)))            # the parents became the contact doublet
      Mrun <- cbind(M[, keep, drop = FALSE], Mi)
      labs <- c(lab$compartment[match(keep, lab$cell)], rep("T", nrow(pr)))
      cr <- call_scdbl(Mrun, labs, rep(FALSE, ncol(Mrun)), dbr)
      out <- cr$res[cell %in% pr$syn]
      out[, share := leuk_share(Mi[, cell, drop = FALSE])]
      out[, st_raw := parse_state(origin)]
      out[, reported := fifelse(called & !is.na(ambiguous) & !ambiguous & !is.na(st_raw), st_raw, "none")]
      out[, reported_win := fifelse(reported != "none" & !is.na(share) & share >= WIN[1] & share <= WIN[2], reported, "none")]
      plant <- data.table(rate = R$rate, n_planted = n, n_unique_parents = uniqueN(c(pr$a, pr$b)), called = sum(out$called),
                          rec_state = sum(out$reported_win == DBL_PLANT_STATE), rec_any = sum(out$reported_win != "none"),
                          rec_state_nowin = sum(out$reported == DBL_PLANT_STATE))
      return(list(kind = "planted", id = R$id, plant = plant,
                  run = data.table(id = R$id, kind = "planted", n_real = length(keep), n_injected = n, secs = cr$secs,
                                   threshold = cr$threshold, n_warn = length(cr$warns), first_warn = substr(c(cr$warns, "")[1], 1, 150)))) }
  }
  ts(sprintf("[3] %d scDblFinder runs on %d cores: %s", length(RUNS), min(CORES, length(RUNS)), paste(vapply(RUNS, `[[`, "", "id"), collapse = ", ")))
  RES <- mclapply(RUNS, function(R) tryCatch(run_one(R), error = function(e) list(kind = R$kind, id = R$id, error = conditionMessage(e))),
                  mc.cores = min(CORES, length(RUNS)), mc.preschedule = FALSE)
  RES <- Map(function(r, R) if (is.list(r) && !inherits(r, "try-error")) r else   # a worker killed outright (e.g. out of memory)
               list(kind = R$kind, id = R$id, error = paste("worker died:", paste(as.character(r), collapse = " "))), RES, RUNS)
  err <- Filter(function(r) !is.null(r$error), RES)
  for (r in err) ts("    ** run ", r$id, " failed: ", r$error)
  chk(sprintf("%s: every scDblFinder run finished", L), 0L, length(err))
  RES <- Filter(function(r) is.null(r$error), RES)
  if (!length(RES)) { ts("    ** no run of this lane finished; nothing to summarise"); rm(x, M); gc(FALSE); next }   # empty-result branch
  RUNTAB <- rbindlist(lapply(RES, `[[`, "run"))[, lane := L]
  xgb_fail <- RUNTAB[, grepl("invalid function arguments|have been removed", first_warn)]
  chk(sprintf("%s: %s", L, if (MODE == "main") "main runs on the pre-training heuristic (xgboost 3.2 rejects 1.20.2's call) in every run"
                           else "devel trains its classifier in every run (no xgboost warning)"),
      MODE == "main", if (MODE == "main") all(xgb_fail) else any(xgb_fail))
  ALL$runs[[L]] <- RUNTAB
  print(RUNTAB[, .(id, kind, n_real, n_injected, secs = round(secs, 1), threshold = round(threshold, 3), n_warn)])
  ## -- [5] real run: format assertions + aggregate controls only ----
  real <- Filter(function(r) r$kind == "real", RES)
  if (length(real)) {
    chk(sprintf("%s: real run output has the expected format (values masked)", L), TRUE, real[[1]]$ok_format)
    ALL$ctrl[[L]] <- real[[1]]$ctrl[, lane := L]
    print(real[[1]]$ctrl) }
  ## -- [6] held-out synthetic: confusion, recall, precision, false-positive sources ----
  SYN <- rbindlist(lapply(Filter(function(r) r$kind %in% c("heldout", "lodo"), RES), function(r) r$out[, `:=`(run = r$id, kind = r$kind)]), fill = TRUE)
  if (nrow(SYN)) {
    tl <- SYN[type == "T|LEUK"]
    CONF <- tl[, .N, by = .(kind, true_state, reported = reported_win)][, lane := L]
    ALL$conf[[L]] <- CONF
    summ <- rbindlist(lapply(c("heldout", "lodo"), function(kd) rbindlist(lapply(STATES, function(s) {
      tk <- tl[kind == kd & true_state == s]; rj <- SYN[kind == kd & reported_win == s]
      data.table(kind = kd, state = s, n = nrow(tk), retained = sum(tk$reported_win != "none"), correct = sum(tk$reported_win == s),
                 correct_nowin = sum(tk$reported == s), reported_from_tleuk = rj[type == "T|LEUK", .N], reported_all = nrow(rj)) }))))
    summ <- summ[n > 0 | reported_all > 0]
    ## where injected T-LEUK doublets are lost: not called -> origin ambiguous -> origin not T + one LEUK state -> ambient window
    lost <- tl[, .(n = .N, called = sum(called), not_ambiguous = sum(called & !is.na(ambiguous) & !ambiguous),
                   origin_tleuk_any = sum(called & !is.na(ambiguous) & !ambiguous & !is.na(st_raw)),
                   origin_right_state = sum(reported == true_state), window_pass = sum(reported_win != "none"),
                   origin_has_T = sum(called & grepl("(^|\\+)T(\\+|$)", origin))), by = .(kind, true_state)]
    ALL$lost[[L]] <- lost[, lane := L]
    print(lost)
    ALL$inj[[L]] <- SYN[, .(lane = L, run, kind, syn, type, true_state, called, score, origin, ambiguous, share, reported, reported_win, tleuk_any)]
    summ[, `:=`(recall = correct / n, recall_given_retained = correct / retained, recall_nowin = correct_nowin / n,
                precision_among_tleuk = correct / reported_from_tleuk, precision_injected_mix = correct / reported_all)]
    ALL$summ[[L]] <- summ[, lane := L]
    print(summ[, .(kind, state, n, retained, correct, recall = round(recall, 2), rec_given_ret = round(recall_given_retained, 2),
                   prec_tleuk = round(precision_among_tleuk, 2), prec_mix = round(precision_injected_mix, 2))])
    FP <- SYN[type != "T|LEUK" & kind == "heldout", .(n = .N, reported_any = sum(tleuk_any),
               reported_by_state = if (startsWith(.BY$type, "singlet:")) NA_character_ else   # singlets: masked by state
                 paste(sprintf("%s:%d", STATES, vapply(STATES, function(s) sum(reported_win == s, na.rm = TRUE), 0L)), collapse = ";")), by = type]
    ALL$fp[[L]] <- FP[, lane := L]
    print(FP)
    AMB <- tl[, .(n = .N, q05 = quantile(share, 0.05, na.rm = TRUE), q50 = median(share, na.rm = TRUE), q95 = quantile(share, 0.95, na.rm = TRUE),
                  in_window = mean(share >= WIN[1] & share <= WIN[2], na.rm = TRUE)), by = true_state]
    ALL$amb[[L]] <- AMB[, lane := L]
    ## -- [7] expectation under random co-encapsulation (model only: singlet composition x synthetic rates) ----
    hk <- tl[kind == "heldout"]
    Cmat <- sapply(STATES, function(k) vapply(c(STATES, "none"), function(j) hk[true_state == k, mean(reported_win == j)], 0))   # (K+1) x K
    no_eval <- STATES[colSums(is.nan(Cmat)) > 0]                    # a state without held-out T-LEUK doublets (empty pool)
    if (length(no_eval)) ts("    ** state(s) not evaluable in held-out runs: ", paste(no_eval, collapse = ", "))
    Cmat[is.nan(Cmat)] <- NA_real_
    Cz <- Cmat; Cz[is.na(Cz)] <- 0                                  # expectations: an unevaluable state contributes nothing
    p_pair <- SYN[kind == "heldout" & type != "T|LEUK" & !startsWith(type, "singlet:"),
                  setNames(lapply(STATES, function(s) mean(reported_win == s)), STATES), by = type]
    p_sing <- SYN[kind == "heldout" & startsWith(type, "singlet:"), .(p_any = mean(tleuk_any)),   # pooled over states (masked)
                  by = .(comp = sub("^singlet:", "", type))]
    EXP <- rbindlist(lapply(donors, function(d) {
      fd  <- A[donor == d, .N, by = compartment][, f := N / sum(N)]
      fv  <- setNames(fd$f, fd$compartment); g <- function(cp) if (cp %in% names(fv)) fv[[cp]] else 0
      pd  <- A[, mean(donor == d)]                                        # donor share of the lane's class-A cells
      Nd  <- N * pd; Nd_dbl <- dbr * N * pd^2                             # same-donor doublets only enter a patient's counts
      fL  <- sum(vapply(STATES, g, 0))
      sig <- Nd_dbl * as.vector(Cz[STATES, , drop = FALSE] %*% vapply(STATES, function(k) 2 * g("T") * g(k), 0))
      pi_u <- c("T|Tamb" = 2 * g("T") * g("Tamb"), "T|NK" = 2 * g("T") * g("NK"), "T|B" = 2 * g("T") * g("B"), "T|Ery" = 2 * g("T") * g("Ery"),
                "T|MonoDC" = 2 * g("T") * g("MonoDC"), "T|LEUK_other" = 2 * g("T") * g("LEUK_other"), "T|T" = g("T")^2,
                "Tamb|LEUK*" = 2 * g("Tamb") * fL, "NK|LEUK*" = 2 * g("NK") * fL, "B|LEUK*" = 2 * g("B") * fL, "Ery|LEUK*" = 2 * g("Ery") * fL,
                "MonoDC|LEUK*" = 2 * g("MonoDC") * fL, "LEUK*|LEUK*" = fL^2)
      fp_pair <- vapply(STATES, function(s) sum(vapply(names(pi_u), function(u) { r <- p_pair[type == u][[s]]; if (length(r)) pi_u[[u]] * r else 0 }, 0)), 0) * Nd_dbl
      ## singlet false positives are known only pooled over states (masked): spread them like the pair-type false positives
      w_fp    <- if (sum(fp_pair) > 0) fp_pair / sum(fp_pair) else rep(1 / length(STATES), length(STATES))   # [SMOKE]
      fp_sing <- (Nd - Nd_dbl) * sum(vapply(names(fv), function(cp) { r <- p_sing[comp == cp, p_any]; if (length(r)) fv[[cp]] * r else 0 }, 0)) * w_fp
      data.table(donor = d, state = STATES, f_state = vapply(STATES, g, 0), f_T = g("T"), exp_signal = sig,
                 exp_signal_correct = Nd_dbl * diag(Cz[STATES, , drop = FALSE]) * 2 * g("T") * vapply(STATES, g, 0),
                 exp_fp_pairs = fp_pair, exp_fp_singlets = fp_sing) }))
    EXP[, `:=`(exp_fp = exp_fp_pairs + exp_fp_singlets)][, precision_natural := exp_signal_correct / (exp_signal + exp_fp)]
    ALL$expd[[L]] <- EXP[, lane := L]
    print(EXP[, .(donor, state, f_state = round(f_state, 3), f_T = round(f_T, 3), exp_signal = round(exp_signal, 2), exp_fp = round(exp_fp, 2),
                  precision_natural = round(precision_natural, 2))])
    saveRDS(list(lane = L, mode = MODE, version = VER, dbr = dbr, N = N, n_dep = n_dep, states = STATES, C = Cmat,
                 n_heldout = hk[, .N, by = true_state], p_pair = p_pair, p_sing = p_sing, expected = EXP,
                 comp = A[, .N, by = .(donor, compartment)]), file.path(DBL_LANE_DIR, paste0(L, "__calib", SFX, ".rds")))
  }
  ## -- [8] planted contact doublets: recovery curve ----
  PL <- rbindlist(lapply(Filter(function(r) r$kind == "planted", RES), `[[`, "plant"))
  if (nrow(PL)) {
    setorder(PL, rate); ALL$plant[[L]] <- PL[, lane := L]; print(PL)
    chk(sprintf("%s: planted contacts recovered as T-LEUK_GMP rise with the planting rate (known answer)", L), TRUE, !is.unsorted(PL$rec_state, strictly = TRUE))
    chk(sprintf("%s: >= 1 planted contact recovered at every rate", L), TRUE, all(PL$rec_state > 0)) }
  rm(x, M, RES); gc(FALSE)
}
w <- function(tab, nm) if (length(tab)) fwrite(rbindlist(tab, fill = TRUE), file.path(DBL_TAB_DIR, paste0("dbl_typing_", nm, SFX, ".csv")))
w(ALL$conf, "confusion"); w(ALL$summ, "summary"); w(ALL$fp, "fp"); w(ALL$expd, "expected"); w(ALL$plant, "planted")
w(ALL$ctrl, "controls"); w(ALL$amb, "ambient"); w(ALL$runs, "runs"); w(ALL$lost, "losses")
if (length(ALL$inj)) fwrite(rbindlist(ALL$inj), file.path(DBL_TAB_DIR, paste0("dbl_typing_injected", SFX, ".csv.gz")))   # singlet rows masked
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, paste0("dbl_checks_04", SFX, ".csv")))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
