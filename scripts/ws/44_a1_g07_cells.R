#!/usr/bin/env Rscript
# INPUT  : ws_pb_GSE185381.rds ($cells, $raw$T, $fine); 02_seurat_objects/04_annotated/GSE185381/<sample>.rds (RNA counts layer);
#          results/tables/ws/ ws_s5_scores.csv (T rows: sample, status), ws_s5_imprint_genes.csv, ws_s8_t_imprint_ambient_prone.csv
# OUTPUT : results/tables/ws/ws_a1_g07_{s5b_checks,s5b_patients,s5b_labels}.csv; under /LARGE1/.../ws_a1_g07/: input_md5_s5b.csv, k6b_per_donor.csv,
#          percell_scores_T1000_P370.rds, s5b_references.csv (per-label reference values and fallbacks of the main builds)
# what it does: PREREGISTRATION_A1_G07.md 3.5 (S5b, T only; a dated amendment of S5 4.2) and section 6 rows K6a-K6e and K8.
#   --known  K6a-K6e, K8 on healthy cells only (no AML score is computed)      --score  3.5 for the 15 AML patients, then the cohort reading
#   --smoke  synthetic cells on the real sample IDs and status (no real count, label or depth enters a smoke computation)
# The rule functions (reference, z, f0, statistics, mixture) follow D2's 13_sim.R mk_ref / zf / f0 / pstats and 17_final_rule.R:8-10.
# This script deletes nothing, modifies no input, and never runs or sources another analysis script.

options(warn = 1)   # warnings print where they happen (smoke review)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here); library(parallel) })
.a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% .a
MODES <- c("known", "score"); MODE <- MODES[paste0("--", MODES) %in% .a]
if (length(MODE) != 1L) stop("give exactly one of --known --score [--smoke]")
bad_arg <- setdiff(.a, c(paste0("--", MODES), "--smoke"))   # a mistyped --smoke (e.g. --smok) would otherwise run on real cells
if (length(bad_arg)) stop("unknown argument(s): ", paste(bad_arg, collapse = " "))
stopifnot(packageVersion("UCell") == "2.10.1")   # 3.5(a)
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R"))
RNGkind(sample.kind = "Rejection")                                                    # 3.5(b)
NC <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))); if (is.na(NC) || NC < 1L) NC <- 2L
TABWS <- file.path(TAB_DIR, "ws"); GDIR <- file.path(LARGE1_DIR, "ws_a1_g07"); SDIR <- file.path(GDIR, "smoke")
WDIR <- if (SMOKE) SDIR else GDIR; ODIR <- if (SMOKE) SDIR else TABWS; dir.create(WDIR, recursive = TRUE, showWarnings = FALSE)
SFX <- if (SMOKE) "_smoke" else ""
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
MD5 <- if (length(SELF) && !is.na(SELF) && file.exists(SELF)) unname(tools::md5sum(SELF)) else NA_character_
if (!SMOKE && is.na(MD5)) stop("run the script with Rscript (its md5 is needed for every output)")
if (!SMOKE && file.access(SELF, 2) == 0) stop("a real run must use the dated read-only copy: sbatch scripts/ws/submit_ws_a1_g07.sbatch 44 --", MODE)   # an edit cannot reach a running job
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
wtab <- function(x) { f <- file.path(ODIR, paste0("ws_a1_g07_", x, SFX, ".csv")); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
put  <- function(d, x) { d <- as.data.table(d); if (!ncol(d)) d <- data.table(script_md5 = character(0)) else d[, script_md5 := MD5]; fwrite(d, wtab(x)); invisible(d) }   # an empty table is written as a header
wg   <- function(x) { f <- file.path(WDIR, sub("(\\.[a-z.]+)$", paste0(SFX, "\\1"), x)); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
rg   <- function(x) file.path(WDIR, sub("(\\.[a-z.]+)$", paste0(SFX, "\\1"), x))
if (!SMOKE) { tg <- switch(MODE, known = c(file.path(TABWS, "ws_a1_g07_s5b_checks.csv"), file.path(GDIR, c("k6b_per_donor.csv", "input_md5_s5b.csv"))),   # work files included
                          score = c(file.path(TABWS, paste0("ws_a1_g07_", c("s5b_patients", "s5b_labels"), ".csv")), file.path(GDIR, c("percell_scores_T1000_P370.rds", "s5b_references.csv"))))
  if (any(file.exists(tg))) stop("refusing to overwrite ", paste(tg[file.exists(tg)], collapse = ", "), " (a rerun needs a dated note and the old file renamed)") }
ts(sprintf("[%s%s] %d worker processes", MODE, if (SMOKE) " smoke" else "", NC))

## -- 0. constants -------------------------------------------------------------------------------------------------------------------
TARGET <- 1000L; TARGET_S <- 1500L; MAXRANK <- 1500L; MIN_CELLS <- 50L; FLOOR <- 7                               # 3.5(b), (e)
REF_MIN_CELLS <- 20L; REF_MIN_DON <- 2L; REF_MIN_PER_DON <- 5L                                                       # 3.5(c)
PIG <- seq(0.1, 1, 0.1); DG <- seq(0, 20, 0.5); SG <- seq(-3, 20, 1)                                               # D2 13_sim.R:52
T15 <- paste0("Control", 1:5); PER_LABEL <- c("CD4 Naive", "CD4 Central Memory", "CD8 Effector Memory 1", "CD8 Effector Memory 2")
CD4L <- c("CD4 Naive", "CD4 Central Memory", "CD4 Effector Memory", "CD4 Regulatory")
seed_thin <- function(i, target) SEED + (if (target == TARGET) 7100L else 7300L) + i                                 # 8.3; i = T-row position
seed_k6 <- function(i) SEED + 7200L + i

## -- 1. data ---------------------------------------------------------------------------------------------------------------------------
S5T <- fread(file.path(TABWS, "ws_s5_scores.csv"), select = c("compartment", "sample", "status"))[compartment == "T", .(sample, status)]   # design columns only (2.1)
S25 <- S5T$sample; AMLS <- S5T[status == "AML", sample]; HEAL <- S5T[status == "healthy", sample]; IDX <- setNames(seq_along(S25), S25)
IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))$gene
F402 <- IMP[, .(gene, s = sign(logFC))]; P370 <- F402[!gene %in% AMB]
SETS <- list(P370 = list(up = P370[s > 0, gene], down = P370[s < 0, gene]), F402 = list(up = F402[s > 0, gene], down = F402[s < 0, gene]))
stopifnot(length(S25) == 25L, lengths(SETS$P370) == c(156L, 214L), lengths(SETS$F402) == c(171L, 231L))
MAP <- fread(here::here("scripts", "ws", WS_S4_MAP_FILE), fill = TRUE); TL <- MAP[group4 == "T", bmm_fine]; CD8L <- setdiff(TL, CD4L)
PB <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds")); G <- PB$genes
KC <- PB$cells[kept == TRUE & group4 == "T" & sample %in% S25, .(cell, sample, fine)]
if (!SMOKE) { RAWT <- PB$raw$T[, S25]; FINE0 <- PB$fine[TL] }
rm(PB); invisible(gc(FALSE))
AN <- file.path(LARGE1_DIR, "02_seurat_objects", "04_annotated", "GSE185381")
IN44 <- c(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"), file.path(TABWS, c("ws_s5_scores.csv", "ws_s5_imprint_genes.csv", "ws_s8_t_imprint_ambient_prone.csv")),
          here::here("scripts", "ws", WS_S4_MAP_FILE), here::here("scripts", "config", "config_paths.R"), here::here("scripts", "ws", "config_ws.R"), file.path(AN, paste0(S25, ".rds")))   # K0 for 44 (2.x incl. 2.4)
load_counts <- function() {   # 40:64-83: per sample, the kept T cells of the 04_annotated object, in $cells order
  if (SMOKE) {   # synthetic: random labels and depths; AML patients carry a uniform (odd) or a 20 % (even) up-set shift
    set.seed(SEED + 710L); mu <- setNames(numeric(length(G)), G); on <- sample.int(length(G), 9000L); mu[on] <- exp(rnorm(9000L, 0, 1.6))
    mu[c(SETS$F402$up, SETS$F402$down)] <- exp(rnorm(402L, 1.5, 0.8))
    out <- lapply(S25, function(s) { a <- s %in% AMLS; n <- if (a) 400L else 220L   # 400 AML cells, so the per-label path (3.5(f)) runs in smoke
      lab <- sample(TL, n, TRUE, prob = c(4, 3, 2, 1, 2, 1, 3, 5, 1)[seq_along(TL)])
      umi <- round(exp(rnorm(n, if (a) log(2900) else log(2000), 0.45))); M <- sapply(seq_len(n), function(c) { m <- mu
        if (a && (IDX[[s]] %% 2 == 1 || c <= n / 5)) m[SETS$P370$up] <- m[SETS$P370$up] * (if (IDX[[s]] %% 2 == 1) 1.8 else 5)
        rmultinom(1, umi[c], m / sum(m))[, 1] })
      X <- Matrix(M, sparse = TRUE); dimnames(X) <- list(G, sprintf("%s_c%04d", s, seq_len(n))); list(X = X, lab = lab, ok = TRUE, synthetic = TRUE) })
    names(out) <- S25; return(out) }
  suppressPackageStartupMessages(library(SeuratObject))
  out <- mclapply(S25, function(s) { o <- readRDS(file.path(AN, paste0(s, ".rds"))); X <- SeuratObject::LayerData(o[["RNA"]], layer = "counts"); rm(o)
    k <- KC[sample == s]; if (!identical(rownames(X), G) || !all(k$cell %in% colnames(X))) return(list(X = NULL, lab = NULL, ok = FALSE))
    list(X = as(X[, k$cell, drop = FALSE], "CsparseMatrix"), lab = k$fine, ok = TRUE) }, mc.cores = min(NC, 8L))
  names(out) <- S25; out }
loaded_or_stop <- function(CNT) { bad <- names(CNT)[!vapply(CNT, function(z) is.list(z) && isTRUE(z$ok), TRUE)]   # name the samples that failed to load
  if (length(bad)) stop("per-cell load failed: ", paste(bad, collapse = " ")) }
need_syn <- function(z, s) if (SMOKE && !isTRUE(z$synthetic)) stop("smoke run: real cells of ", s, " would enter a computation")

## -- 2. thinning (3.5(b)), UCell (3.5(a)) -------------------------------------------------------------------------------------------
thin <- function(X, target, drop_below, seed) {   # exact downsampling without replacement, one sample per call, cells in $cells order
  set.seed(seed); u <- Matrix::colSums(X); keep <- if (drop_below) u >= target else rep(TRUE, length(u)); X <- X[, keep, drop = FALSE]
  X <- as(X, "CsparseMatrix"); for (j in seq_len(ncol(X))) { a <- X@p[j] + 1L; b <- X@p[j + 1L]; if (b < a) next; x <- X@x[a:b]; tot <- sum(x)
    if (tot > target) { idx <- rep.int(seq_along(x), x); X@x[a:b] <- tabulate(idx[sample.int(length(idx), target)], nbins = length(x)) } }
  list(X = drop0(X), keep = keep) }
ucell <- function(X, sets) { sc <- UCell::ScoreSignatures_UCell(matrix = X, features = sets, maxRank = MAXRANK, ties.method = "average", BPPARAM = BiocParallel::SerialParam(), chunk.size = 500)   # serial: no socket per call
  sc[, paste0(names(sets), "_UCell"), drop = FALSE] }
detected <- function(X, g) Matrix::colSums(X[intersect(g, rownames(X)), , drop = FALSE] > 0)
pipe_scores <- function(CNT, pipe, sets) {   # per-cell scores of every sample under one pipeline
  rbindlist(mcl(names(CNT), function(s) { need_syn(CNT[[s]], s); X <- CNT[[s]]$X; lab <- CNT[[s]]$lab
    if (pipe == "raw") { Y <- X; kp <- rep(TRUE, ncol(X)) } else { th <- thin(X, if (pipe == "T1000") TARGET else TARGET_S, pipe == "T1000", seed_thin(IDX[[s]], if (pipe == "T1000") TARGET else TARGET_S)); Y <- th$X; kp <- th$keep }
    if (!ncol(Y)) return(data.table())
    sc <- ucell(Y, sets); data.table(sample = s, cell = colnames(Y), label = lab[kp], umi = Matrix::colSums(Y), up = sc[, 1], down = sc[, 2],
                                     det_up = detected(Y, sets$up), det_down = detected(Y, sets$down)) })) }
mcl <- function(X, FUN) { r <- mclapply(X, FUN, mc.cores = NC, mc.preschedule = FALSE)
  bad <- vapply(r, function(z) is.null(z) || inherits(z, "try-error"), TRUE); if (any(bad)) stop(sum(bad), " of ", length(r), " workers failed: ", paste(unique(unlist(lapply(r[bad], as.character))), collapse = " | ")); r }

## -- 3. reference, z, f0, statistics (D2 13_sim.R:44-67; review-2 minimum-reference rule) ------------------------------------------
build_ref <- function(SC, set, donors) {   # SC: per-cell scores of healthy cells; reference from `donors` only
  R <- SC[sample %in% donors]; s <- R[[set]]
  pool_of <- function(k) if (k %in% CD4L) CD4L else CD8L
  par <- rbindlist(lapply(TL, function(k) { nd <- R[label == k, .N, by = sample]; own <- R[label == k, .N] >= REF_MIN_CELLS && sum(nd$N >= REF_MIN_PER_DON) >= REF_MIN_DON
    src <- if (own) k else pool_of(k); e <- R[label %in% src][[set]]; a <- mean(e); sdv <- max(sd(e - a), 1e-6)
    data.table(label = k, a = a, sd = sdv, z95 = unname(quantile((e - a) / sdv, 0.95, type = 7)), z5 = unname(quantile((e - a) / sdv, 0.05, type = 7)),
               source = if (own) "own" else if (k %in% CD4L) "pooled CD4" else "pooled CD8", n_ref = length(e)) }))
  tab <- table(factor(R$label, levels = TL), R$sample)   # 3.5(c) "the code asserts the rule": a second, table-based count must give the same sources
  own2 <- vapply(TL, function(k) sum(tab[k, ]) >= REF_MIN_CELLS && sum(tab[k, ] >= REF_MIN_PER_DON) >= REF_MIN_DON, TRUE)
  if (!identical(unname(own2), par$source == "own")) stop("minimum-reference rule: the two counts disagree")
  zr <- (s - par$a[match(R$label, par$label)]) / par$sd[match(R$label, par$label)]
  d <- density(zr, n = 2048, from = min(zr) - 3, to = max(zr) + 8)
  list(par = par, zr = zr, dx = d$x, dy = pmax(d$y, 1e-8), v = var(zr)) }
fb <- function(r) { x <- r$par[source != "own"]; if (!nrow(x)) "none" else paste(sprintf("%s -> %s (n %d)", x$label, x$source, x$n_ref), collapse = "; ") }   # every fallback, reported
zf <- function(r, s, lab) (s - r$par$a[match(lab, r$par$label)]) / r$par$sd[match(lab, r$par$label)]
f0 <- function(r, x) approx(r$dx, r$dy, x, rule = 2)$y
pstats <- function(r, s, lab, mix = TRUE, down = FALSE) {
  z <- zf(r, s, lab); n <- length(z); mz <- mean(z)
  F95 <- mean(z > r$par$z95[match(lab, r$par$label)]); F5 <- mean(z < r$par$z5[match(lab, r$par$label)]); logVR <- log(var(z) / r$v)
  LR <- NA_real_; pih <- NA_real_
  if (mix && n > 3) { lu <- function(d) sum(log(f0(r, z - d))); g0 <- SG[which.max(sapply(SG, lu))]
    llU <- optimize(lu, c(g0 - 1, g0 + 1), maximum = TRUE)$objective
    F0 <- f0(r, z); FD <- sapply(DG, function(d) f0(r, z - d)); llM <- -Inf; pih <- 1
    for (p in PIG) { ll <- colSums(log((1 - p) * F0 + p * FD)); if (max(ll) > llM + 1e-9) { llM <- max(ll); pih <- p } }
    if (llM - llU < 1e-6) pih <- 1; LR <- max(0, llM - llU) }
  data.table(n = n, mean_z = mz, F95 = F95, F5 = F5, se = sqrt(exp(logVR) / n), LR = LR, pihat = pih) }
call_of <- function(st, hmax, n_min = MIN_CELLS) {   # 3.5(e) up set; 17_final_rule.R:8-10
  if (st$n < n_min) return("too few cells"); if (!(st$mean_z - 1.96 * st$se > hmax)) return("no detectable shift")
  if (st$pihat <= 0.5 && st$LR >= 10) "sub-population" else if (st$pihat >= 0.9 || st$LR < 2) "uniform" else "not separable" }
down_gate <- function(st, hmin) if (st$n < MIN_CELLS) "too few cells" else if (st$mean_z + 1.96 * st$se < hmin) "shift" else "no detectable shift"
lodo_means <- function(SC, set, donors) vapply(donors, function(e) { r <- build_ref(SC, set, setdiff(donors, e)); x <- SC[sample == e]; mean(zf(r, x[[set]], x$label)) }, 1)

## -- 4. --known: K8 and K6 (healthy cells only) -------------------------------------------------------------------------------------
if (MODE == "known") {
  CHK <- list(); chk <- function(id, name, expected, got, kind = "deterministic", ok = isTRUE(all.equal(expected, got))) {
    CHK[[length(CHK) + 1L]] <<- data.table(id = id, kind = kind, check = name, expected = paste(format(expected), collapse = ";"), got = paste(format(got), collapse = ";"),
                                           state = if (kind == "reported") "REPORTED" else if (ok) "PASS" else if (kind == "statistical") "FAIL_statistical" else "FAIL_deterministic")
    message(sprintf("  %-4s %-92s expect %-14s got %-14s %s", id, name, paste(format(expected), collapse = ";"), paste(format(got), collapse = ";"), if (kind == "reported") "[REPORTED]" else if (ok) "PASS" else "** FAIL **")) }
  ts("[K0] input md5s (44's inputs, incl. the 25 per-sample objects of 2.4)")
  md5 <- unname(tools::md5sum(IN44)); M43 <- file.path(WDIR, paste0("input_md5", SFX, ".csv"))
  if (!file.exists(M43)) stop("run 43 --known first: ", M43)
  I43 <- fread(M43); sh <- intersect(IN44, I43$path)
  chk("K0", "all 44 inputs exist (7 files + 25 per-sample objects)", length(IN44), sum(!is.na(md5)))
  chk("K0", "ws_pb_GSE185381.rds md5 as recorded", "6c2fb685c49d4a23ef20a6842aa7aa67", md5[1])
  chk("K0", sprintf("the %d inputs shared with 43 --known have the md5s it recorded", length(sh)), I43$md5[match(sh, I43$path)], md5[match(sh, IN44)])
  fwrite(data.table(path = IN44, md5 = md5, script_md5 = MD5), wg("input_md5_s5b.csv"))
  CNT <- load_counts()
  ts("[K8] per-cell inputs")
  chk("K8", "every sample loaded on the pseudobulk gene axis with all kept T cells", 25L, sum(vapply(CNT, function(z) is.list(z) && isTRUE(z$ok), TRUE)))
  if (sum(vapply(CNT, function(z) is.list(z) && isTRUE(z$ok), TRUE)) < 25L) { fwrite(rbindlist(CHK)[, script_md5 := MD5], wtab("s5b_checks")); loaded_or_stop(CNT) }   # the FAIL row reaches disk
  if (!SMOKE) {
    chk("K8", "kept-T sums equal $raw$T (max |diff| over the 25)", 0, max(vapply(S25, function(s) max(abs(Matrix::rowSums(CNT[[s]]$X) - RAWT[, s])), 1)))
    lab_diff <- function(s) max(vapply(TL, function(k) { j <- CNT[[s]]$lab == k; a <- if (any(j)) Matrix::rowSums(CNT[[s]]$X[, j, drop = FALSE]) else 0   # per gene
      b <- if (!is.null(FINE0[[k]]) && s %in% colnames(FINE0[[k]])) FINE0[[k]][, s] else 0; max(abs(a - b)) }, 1))
    chk("K8", "label sums equal $fine (max |diff| over genes, labels and samples)", 0, max(vapply(S25, lab_diff, 1)))
    chk("K8", "kept T cells in the 25", 14609L, sum(vapply(CNT, function(z) ncol(z$X), 1L)))
  }
  chk("K8", "every P370 and F402 set gene is in the matrix rownames", TRUE, all(unlist(SETS) %in% rownames(CNT[[1]]$X)))
  chk("K8", "labels are exactly the nine T labels, no NA", TRUE, all(unlist(lapply(CNT, `[[`, "lab")) %in% TL) && !anyNA(unlist(lapply(CNT, `[[`, "lab"))))
  ts("[K6a] thinning and UCell engine")
  s1 <- HEAL[1]; a1 <- thin(CNT[[s1]]$X, TARGET, TRUE, seed_thin(IDX[[s1]], TARGET)); a2 <- thin(CNT[[s1]]$X, TARGET, TRUE, seed_thin(IDX[[s1]], TARGET))
  chk("K6a", "same seed twice: identical thinned counts", TRUE, identical(a1$X, a2$X))
  chk("K6a", "every kept cell sums to exactly 1,000; cells below 1,000 dropped", TRUE, all(Matrix::colSums(a1$X) == TARGET) && sum(!a1$keep) == sum(Matrix::colSums(CNT[[s1]]$X) < TARGET))
  chk("K6a", "UCell twice on the same thinned donor: identical scores", TRUE, identical(ucell(a1$X, SETS$P370), ucell(a2$X, SETS$P370)))
  ts("[K6e, K6d] healthy scores under S5b-1000; detection floor")
  SCH <- pipe_scores(CNT[HEAL], "T1000", SETS$P370); SCH <- SCH[sample %in% HEAL]
  sc1 <- ucell(a1$X, SETS$P370); h1 <- SCH[sample == s1]   # a1: s1 thinned alone, in this process; SCH: the same donor inside the forked multi-donor run
  chk("K6a", "a donor thinned alone equals the same donor in the full run (cells, UMI, up, down)", TRUE, identical(h1$cell, colnames(a1$X)) &&
      identical(unname(h1$umi), unname(Matrix::colSums(a1$X))) && identical(unname(sc1[, 1]), unname(h1$up)) && identical(unname(sc1[, 2]), unname(h1$down)))
  med_det <- median(SCH$det_up)
  chk("K6e", "median detected up-set genes per healthy reference cell (S5b-1000)", ">= 7", med_det, kind = "statistical", ok = med_det >= FLOOR)
  chk("K6d", "down-set detected genes per healthy cell (median)", "", median(SCH$det_down), kind = "reported")
  ts("[K6b-K6c] planted shape on real healthy cells (nested leave-one-out)")
  K6 <- rbindlist(mcl(HEAL, function(d) { need_syn(CNT[[d]], d); i <- IDX[[d]]; oth <- setdiff(HEAL, d)
    rU <- build_ref(SCH, "up", oth); hmax <- max(vapply(oth, function(e) { r <- build_ref(SCH, "up", setdiff(oth, e)); x <- SCH[sample == e]; mean(zf(r, x$up, x$label)) }, 1))
    rD <- build_ref(SCH, "down", oth); hmin <- min(vapply(oth, function(e) { r <- build_ref(SCH, "down", setdiff(oth, e)); x <- SCH[sample == e]; mean(zf(r, x$down, x$label)) }, 1))
    X <- CNT[[d]]$X; lab <- CNT[[d]]$lab; up <- SETS$P370$up; umi <- Matrix::colSums(X)
    prof <- lapply(setNames(TL, TL), function(k) { v <- Reduce(`+`, lapply(oth, function(e) { j <- CNT[[e]]$lab == k; if (any(j)) Matrix::rowSums(CNT[[e]]$X[up, j, drop = FALSE]) else numeric(length(up)) }))
      tot <- sum(vapply(oth, function(e) { j <- CNT[[e]]$lab == k; if (any(j)) sum(CNT[[e]]$X[, j]) else 0 }, 1)); if (tot > 0) v / tot else NULL })
    avail <- names(Filter(Negate(is.null), prof))
    for (k in TL) if (is.null(prof[[k]])) { src <- intersect(if (k %in% CD4L) CD4L else CD8L, avail); prof[[k]] <- Reduce(`+`, prof[src]) / length(src) }   # lineage mean
    MU <- sapply(seq_len(ncol(X)), function(c) umi[c] * prof[[lab[c]]])                   # up genes x cells, mu-hat
    set.seed(seed_k6(i)); beta <- 0.6 * pmin(2.5, 0.5 + rexp(length(up), rate = 1 / 0.4))
    Xa <- X; Xa[up, ] <- Xa[up, ] + matrix(rpois(length(MU), MU * (2^beta - 1)), nrow(MU))
    S <- sample.int(ncol(X), round(ncol(X) / 5)); m <- (2^beta - 1) * rowSums(MU) / pmax(rowSums(MU[, S, drop = FALSE]), 1e-12)
    Xb <- X; Xb[up, S] <- Xb[up, S] + matrix(rpois(length(MU[, S]), MU[, S] * m), length(up))
    k6c <- function(Y) log2((sum(Y[up, ]) / sum(Y)) / (sum(X[up, ]) / sum(X)))
    one <- function(Y) { th <- thin(Y, TARGET, TRUE, seed_thin(i, TARGET)); sc <- ucell(th$X, SETS$P370); l <- lab[th$keep]
      st <- pstats(rU, sc[, 1], l); list(call = call_of(st, hmax), gate = st$mean_z - 1.96 * st$se > hmax, dgate = down_gate(pstats(rD, sc[, 2], l, mix = FALSE), hmin), st = st) }
    pa <- one(Xa); pb <- one(Xb); p0 <- one(X)
    data.table(donor = d, n = ncol(X), share_dropped = mean(umi < TARGET), share_thinned = mean(umi > TARGET), hmax = hmax,
               call_a = pa$call, call_b = pb$call, call_unplanted = p0$call, gate_a = pa$gate, gate_b = pb$gate, down_gate_unplanted = p0$dgate,
               k6c_a = k6c(Xa), k6c_b = k6c(Xb), lr_a = pa$st$LR, pi_a = pa$st$pihat, lr_b = pb$st$LR, pi_b = pb$st$pihat,
               fallback_up = fb(rU), fallback_down = fb(rD), tall_pool = d %in% T15) }))
  fwrite(K6[, script_md5 := MD5], wg("k6b_per_donor.csv"))
  sub_ok <- K6[call_b == "sub-population", .N] >= 8L && K6[call_a == "sub-population", .N] <= 1L; uni_ok <- K6[call_a == "uniform", .N] >= 7L
  chk("K6b", "graded: (b) sub-population >= 8/10 and (a) sub-population <= 1/10 -> sub-population calls readable", "TRUE",
      c(K6[call_b == "sub-population", .N], K6[call_a == "sub-population", .N]), kind = "statistical", ok = sub_ok)
  chk("K6b", "graded: (a) uniform >= 7/10 -> uniform calls readable", "TRUE", K6[call_a == "uniform", .N], kind = "statistical", ok = uni_ok)
  chk("K6b", "gate-pass rates of (a) and (b)", "", c(K6[, sum(gate_a)], K6[, sum(gate_b)]), kind = "reported")
  k6c <- median(abs(K6$k6c_a - K6$k6c_b)); chk("K6c", "median over donors of |pseudobulk up-set log2FC (a) - (b)|", "<= 0.05", sprintf("%.4f", k6c), kind = "statistical", ok = k6c <= 0.05)
  chk("K6d", "unplanted held-out donors: no detectable shift (cannot fail by construction; false-shift rate <= 1/10 per held-out donor)", "", K6[call_unplanted == "no detectable shift", .N], kind = "reported")
  chk("K6d", "unplanted held-out donors: down-set gate 'no detectable shift'", "", K6[down_gate_unplanted == "no detectable shift", .N], kind = "reported")
  chk("K6d", "reference fallbacks of the K6 builds (up; down)", "", c(paste(unique(K6$fallback_up), collapse = " | "), paste(unique(K6$fallback_down), collapse = " | ")), kind = "reported")
  for (s in S25) { u <- Matrix::colSums(CNT[[s]]$X)   # depth only (UMI per cell; no score), printed after every K6 state (dated note)
    chk("K6d", sprintf("%s%s: share dropped (< 1,000 UMI); share thinned (S5b-1000); share thinned (S5b-1500, > 1,500); median UMI raw", s, if (s %in% T15) " (T-ALL lane)" else ""), "",
        c(sprintf("%.3f", mean(u < TARGET)), sprintf("%.3f", mean(u > TARGET)), sprintf("%.3f", mean(u > TARGET_S)), median(u)), kind = "reported") }
  C <- rbindlist(CHK)[, script_md5 := MD5]; fwrite(C, wtab("s5b_checks"))
  ts(sprintf("[known] %d rows; deterministic FAIL %d; statistical FAIL %d", nrow(C), C[state == "FAIL_deterministic", .N], C[state == "FAIL_statistical", .N])); quit(status = 0)
}

## -- 5. --score: 3.5 ----------------------------------------------------------------------------------------------------------------
if (MODE == "score") {
  C <- fread(file.path(ODIR, paste0("ws_a1_g07_s5b_checks", SFX, ".csv")))
  if (!SMOKE && !identical(unique(C$script_md5), MD5)) stop("the S5b checks were written by another version of this script; rerun --known")
  if (any(C$state == "FAIL_deterministic")) stop("a deterministic S5b known answer failed")
  I44 <- fread(rg("input_md5_s5b.csv"))   # 44's inputs must be unchanged since --known
  if (!SMOKE && !identical(unique(I44$script_md5), MD5)) stop("input_md5_s5b.csv was written by another version of this script; rerun --known")
  if (!identical(I44$path, IN44) || !identical(unname(tools::md5sum(IN44)), I44$md5)) stop("44's inputs changed since --known")
  SUB_OK <- identical(C[id == "K6b" & grepl("sub-population calls readable", check), state], "PASS")
  UNI_OK <- identical(C[id == "K6b" & grepl("uniform calls readable", check), state], "PASS")
  K6C_OK <- identical(C[id == "K6c", state], "PASS"); FLOOR_OK <- identical(C[id == "K6e", state], "PASS")
  SUB_R <- SUB_OK && K6C_OK; UNI_R <- UNI_OK && K6C_OK   # a K6c failure leaves both shape calls unreadable (dated note)
  GATEAB <- as.numeric(strsplit(C[id == "K6b" & grepl("gate-pass", check), got], ";")[[1]])
  CNT <- load_counts(); loaded_or_stop(CNT); OUT <- list(); LAB <- list(); REFS <- list()
  read_call <- function(cl, descriptive) {   # a call as the K6 states allow it to be read (3.5(e)2, section 6); F402 reuses the P370 states (3.5(a))
    out <- cl; shp <- out %in% c("uniform", "sub-population", "not separable")
    if (!FLOOR_OK) out[shp] <- "shape not readable (detection floor)" else {
      if (!SUB_R) out[out == "sub-population"] <- "sub-population (not readable; counted as not separable)"
      if (!UNI_R) out[out == "uniform"] <- "uniform or not separable" }
    if (descriptive) out[nzchar(out)] <- paste(out[nzchar(out)], "(descriptive)"); out }
  run_view <- function(pipe, setname, ref_donors, view) {
    raw <- pipe == "raw"; desc <- view %in% c("sensitivity S5b-1500", "sensitivity: reference without Control1-5")
    SC <- pipe_scores(CNT, pipe, SETS[[setname]])
    rU <- build_ref(SC, "up", ref_donors); rD <- build_ref(SC, "down", ref_donors)
    REFS[[view]] <<- rbind(copy(rU$par)[, ref := "up"], copy(rD$par)[, ref := "down"])[, `:=`(view = view, pipeline = pipe, sets = setname)]   # 3.5(c) reference values, main build
    ts(sprintf("[score] %s: reference fallbacks, up: %s; down: %s", view, fb(rU), fb(rD)))
    hmax <- max(lodo_means(SC, "up", ref_donors)); hmin <- min(lodo_means(SC, "down", ref_donors))
    P <- rbindlist(lapply(c(AMLS, ref_donors), function(s) { x <- SC[sample == s]; heal <- s %in% ref_donors
      ru <- if (heal) build_ref(SC, "up", setdiff(ref_donors, s)) else rU; rd <- if (heal) build_ref(SC, "down", setdiff(ref_donors, s)) else rD
      su <- pstats(ru, x$up, x$label, mix = !raw); sd_ <- pstats(rd, x$down, x$label, mix = FALSE)   # S5b-raw: the registered summaries only (3.5(d))
      data.table(view = view, pipeline = pipe, sets = setname, sample = s, status = if (heal) "healthy (leave-one-out)" else "AML", n_cells = nrow(x),
                 mean_ucell_up = mean(x$up), mean_ucell_down = mean(x$down), mean_z_up = su$mean_z, mean_z_down = sd_$mean_z, F95_up = su$F95, F95_down = sd_$F95, F5_down = sd_$F5,
                 LR = su$LR, pihat = su$pihat, hmax = hmax, hmin = hmin, call_up = if (heal || raw) "" else call_of(su, hmax), gate_down = if (heal || raw) "" else down_gate(sd_, hmin),
                 det_up_median = median(x$det_up), det_down_median = median(x$det_down), ref_fallback_up = fb(ru), ref_fallback_down = fb(rd)) }))
    P[, call_up_read := if (raw) "" else read_call(call_up, desc)]
    P[, depth_warning := if (raw) "raw counts: AML T cells are deeper (median UMI 2,952 vs 2,000); in D2's synthetic study depth alone shifted null AML patients by +0.45 healthy SD (up set; 1.4)" else ""]
    if (pipe == "T1000" && setname == "P370" && view == "primary") LAB[[view]] <<- rbindlist(lapply(AMLS, function(s) rbindlist(lapply(PER_LABEL, function(k) {
      x <- SC[sample == s & label == k]; if (nrow(x) < MIN_CELLS) return(NULL); st <- pstats(rU, x$up, x$label)
      data.table(sample = s, label = k, n_cells = nrow(x), mean_z = st$mean_z, LR = st$LR, pihat = st$pihat) }))))
    if (view == "primary" && pipe == "T1000" && setname == "P370") { setattr(SC, "script_md5", MD5); saveRDS(SC, wg("percell_scores_T1000_P370.rds")) }
    P }
  OUT$prim <- run_view("T1000", "P370", HEAL, "primary")
  OUT$raw <- run_view("raw", "P370", HEAL, "registered S5b-raw")
  OUT$s1500 <- run_view("T1500", "P370", HEAL, "sensitivity S5b-1500")
  OUT$noT15 <- run_view("T1000", "P370", setdiff(HEAL, T15), "sensitivity: reference without Control1-5")
  OUT$f402 <- run_view("T1000", "F402", HEAL, "sensitivity: F402")
  PT <- rbindlist(OUT)
  fwrite(rbindlist(REFS)[, script_md5 := MD5], wg("s5b_references.csv"))
  ## cohort reading (3.5(e)5; 7.5), primary view
  A <- PT[view == "primary" & status == "AML"]; elig <- A[n_cells >= MIN_CELLS]; SHAPES <- c("uniform", "sub-population", "not separable")
  raw_calls <- elig$call_up
  ns_raw <- if (FLOOR_OK) sum(raw_calls == "sub-population") else 0L; nn_raw <- if (FLOOR_OK) sum(raw_calls == "not separable") else 0L   # counted before any conversion (7.5)
  calls <- raw_calls; if (!FLOOR_OK) calls[calls %in% SHAPES] <- "shape not readable"
  if (!SUB_R) calls[calls == "sub-population"] <- "not separable"   # reported, but counted as not separable (3.5(e)5)
  sh <- calls[calls %in% SHAPES]; npass <- length(sh)
  nu <- sum(sh == "uniform"); ns <- sum(sh == "sub-population"); nn <- sum(sh == "not separable"); nz <- sum(calls == "no detectable shift")
  COH <- if (npass + sum(calls == "shape not readable") < 6L || (!SUB_R && !UNI_R) || !FLOOR_OK) "not readable" else
    if (nu >= 2 / 3 * npass) { if (length(GATEAB) == 2 && GATEAB[2] < 0.8 * GATEAB[1]) sprintf("mostly uniform among patients with a detectable shift (the gate favours uniform shifts here: K6b (a) %d/10, (b) %d/10)", GATEAB[1], GATEAB[2]) else "mostly uniform" } else
    if (SUB_R && ns >= 2 / 3 * npass) "mostly sub-population" else if (SUB_R && nu >= 2 && ns >= 2) "mixed" else if (nn >= npass / 2) "mostly not separable" else "no majority pattern"
  if (!UNI_R) COH <- sub("uniform", "uniform or not separable", COH)
  sent <- if (!FLOOR_OK) sprintf("Shape not readable (detection floor); %d of %d showed no detectable shift. Cohort reading: %s.", nz, nrow(elig), COH) else
    sprintf(paste("With depth-equalised per-cell scores (cells thinned to 1,000 UMI), among the %d AML patients with \u2265 50 eligible T cells, %d were %s, %d %s, %d not separable",
                  "and %d showed no detectable shift (patient-level summaries). Cohort reading: %s."),
            nrow(elig), nu, if (UNI_R) "uniform" else "uniform or not separable", ns_raw, if (SUB_R) "sub-population" else "sub-population (not readable; counted as not separable)", nn_raw, nz, COH)
  sent <- paste(Filter(nzchar, c(sent, sprintf("%d of %d showed a detectable down-set shift (descriptive).", sum(elig$gate_down == "shift"), nrow(elig)),
                                 if (grepl("uniform|mixed", COH)) "A uniform shift is what a change in all cells would give, and also what handling, ambient RNA or a residual depth difference would give." else "",
                                 if (grepl("sub-population|mixed", COH)) "A minority-cell call can be composition below the label level, a clonal expansion, label error, cells rich in ambient RNA, or doublets." else "",
                                 "\"Not separable\" is not evidence for either shape.")), collapse = " ")
  PT[, `:=`(cohort_reading = "", cohort_sentence = "")][view == "primary", `:=`(cohort_reading = COH, cohort_sentence = sent)]   # the cohort reading belongs to the primary view
  PT[, `:=`(sub_readable = SUB_R, uniform_readable = UNI_R, floor_ok = FLOOR_OK, k6c_ok = K6C_OK,
            tall_pool = fifelse(sample %in% T15, "T-ALL lane (P02/P03, captured with T-ALL/ETP-ALL samples; soup not removed)", ""))]   # Control1-5 flagged in every table (2.1)
  ts("[score] by-status detection (printed after K6b and K6e have states)")
  print(PT[view == "primary", .(det_up = median(det_up_median), det_down = median(det_down_median)), by = status])
  if (SMOKE) PT[, status := "smoke"]   # 8.5: no real status next to a smoke statistic (set after the cohort code used it)
  put(PT, "s5b_patients"); put(if (length(LAB)) rbindlist(LAB) else data.table(), "s5b_labels")
  ts(sprintf("[score] %s", sent)); quit(status = 0)
}
