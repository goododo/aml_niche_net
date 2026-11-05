#!/usr/bin/env Rscript
# INPUT  : the inputs of PREREGISTRATION_A1_G03.md section 2.2 (tag prereg-a1-g03-v1), each md5-checked by K0:
#          ws_pb_GSE185381.rds ($raw$T, $raw$NK, $cells); results/tables/ws/ ws_s5_scores.csv, ws_s5_imprint_genes.csv,
#          ws_s5_de_T.csv, ws_s8_t_imprint_ambient_prone.csv, ws_a1_g02_s5_perm_draws.csv, ws_a1_batch_covariates.csv,
#          decontx/GSE185381/ws_s1_ambient_blacklist.csv; the September pseudobulks /LARGE1/.../12_pseudobulk_de/<ds>/<sample>__pseudobulk.rds
#          and depth tables (rung samples only, opened by explicit name); pb_sample_manifest.csv; the K3 targets (secondary bins)
# OUTPUT : results/tables/ws/ws_a1_g03_{checks,ladder,extend,sens,p370,read,summary}.csv; draws, fits and null counts under
#          /LARGE1/.../ws_a1_g03/; GATES.md block A1-G03 (--combine). Smoke runs write only to /LARGE1/.../ws_a1_g03/smoke/ (*_smoke).
# what it does: PREREGISTRATION_A1_G03.md, nothing else.
#   --known         section 6's known answers (K0-K7) and the stored permutation draws; computes no rung count beyond R0 (K1)
#   --ladder-obs    observed fits of R0-R6, 5.1, 5.2, 5.8, K6       --ladder-perm --task=N   one chunk of null draws (SLURM array)
#   --ladder-state  rung states, extension plan (3.2), 5.7           --extend   E1, E2, 5.5, 5.6     --sens   section 4
#   --read          5.3 and 5.4                                       --combine  section 7, GATES block
#   --smoke         with any mode: synthetic counts and covariates on the real sample and library IDs; outputs end in _smoke
#   usage: sbatch scripts/ws/submit_ws_a1_g03.sbatch --known   (the wrapper runs a dated read-only copy of this file)
# This script deletes nothing, modifies no input, and never runs or sources another analysis script: the functions it needs are
#   copied below with their source lines. PAUMTZ and PAXMIJ are never read (G04 2.3).

suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma); library(statmod); library(Matrix); library(here); library(parallel) })
.a <- commandArgs(TRUE); SMOKE_B <- "--smoke-b" %in% .a; SMOKE <- "--smoke" %in% .a || SMOKE_B   # --smoke-b: planted signal on the S5 side only (R5, R6 and E lose it)
MODES <- c("known", "ladder-obs", "ladder-perm", "ladder-state", "extend", "sens", "read", "combine")
MODE <- MODES[paste0("--", MODES) %in% .a]
if (length(MODE) != 1L) stop("give exactly one of ", paste0("--", MODES, collapse = " "), " [--smoke]")
TASK <- suppressWarnings(as.integer(sub("^--task=", "", grep("^--task=", .a, value = TRUE)[1])))
ROUND <- if ("--round=2" %in% .a) 2L else 1L
source(here::here("scripts", "config", "config_paths.R"))   # SEED and the roots; its md5 is recorded by K0 (config_ws.R is not needed)
NC <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))); if (is.na(NC) || NC < 1L) NC <- 2L
TABWS <- file.path(TAB_DIR, "ws"); SEPT_TAB <- file.path(TAB_DIR, "12_pseudobulk_de")
GDIR <- file.path(LARGE1_DIR, "ws_a1_g03"); SDIR <- file.path(GDIR, if (SMOKE_B) "smoke_b" else "smoke")
WDIR <- if (SMOKE) SDIR else GDIR                     # draws, fits, null chunks
ODIR <- if (SMOKE) SDIR else TABWS                    # tables: smoke tables never go under results/ (the 23:30 cron)
dir.create(file.path(WDIR, "null"), recursive = TRUE, showWarnings = FALSE)
SFX <- if (SMOKE) "_smoke" else ""
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])   # path of the script copy being run
MD5 <- if (length(SELF) && !is.na(SELF) && file.exists(SELF)) unname(tools::md5sum(SELF)) else NA_character_
if (!SMOKE && is.na(MD5)) stop("run the script with Rscript (its md5 is needed for every output)")
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
## write guards (41's pattern): the only places this script writes; a real (non-smoke) file is never overwritten; tables carry the md5
wtab <- function(x) { f <- file.path(ODIR, paste0("ws_a1_g03_", x, SFX, ".csv")); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
put  <- function(d, x) { d <- as.data.table(d); d[, script_md5 := MD5]; fwrite(d, wtab(x)); invisible(d) }
rtab <- function(x) fread(file.path(ODIR, paste0("ws_a1_g03_", x, SFX, ".csv")))
wg   <- function(x) { f <- file.path(WDIR, sub("(\\.[a-z]+)$", paste0(SFX, "\\1"), x)); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
rg   <- function(x) file.path(WDIR, sub("(\\.[a-z]+)$", paste0(SFX, "\\1"), x))
ts(sprintf("[%s%s] %d worker processes", MODE, if (SMOKE) " smoke" else "", NC))
if (!SMOKE) {   # every fixed target of this mode is checked before any work (41's pattern); --ladder-perm checks its chunk file itself
  tg <- switch(MODE, known = c(file.path(TABWS, "ws_a1_g03_checks.csv"), file.path(GDIR, c(paste0("draws_", c("free25", "lib25", "free24", "lib24", "ext25", "ext24", "enum25", "enum24"), ".csv"), "input_md5.csv"))),
               `ladder-obs` = c(file.path(GDIR, "obs_fits.rds"), file.path(TABWS, c("ws_a1_g03_checks_obs.csv", "ws_a1_g03_p370.csv"))),
               `ladder-state` = if (ROUND == 1L) file.path(GDIR, c("ext_plan.csv", "ladder_round1.csv")) else file.path(TABWS, "ws_a1_g03_ladder.csv"),
               extend = file.path(TABWS, "ws_a1_g03_extend.csv"), sens = file.path(TABWS, "ws_a1_g03_sens.csv"),
               read = file.path(TABWS, "ws_a1_g03_read.csv"), combine = file.path(TABWS, "ws_a1_g03_summary.csv"), character(0))
  if (length(tg) && any(file.exists(tg))) stop("refusing to overwrite ", paste(tg[file.exists(tg)], collapse = ", "), " (section 6: a rerun needs a dated note and the old file renamed)")
  if (MODE == "combine" && any(grepl("A1-G03 start", readLines(here::here("scripts", "ws", "GATES.md"))))) stop("GATES.md already holds an A1-G03 block")
}

## -- 0. constants (section numbers refer to the pre-registration) ----------------------------------------------------------
Q_S5 <- 0.10; L_S5 <- 0.5          # S5's rule: BH q < 0.10 and |logFC| > 0.5 (strict)
Q_SEP <- 0.05; L_SEP <- 1.0        # September's rule: q < 0.05 and |logFC| >= 1 (02:67-69)
ALPHA <- 0.05; N_DRAW <- if (SMOKE) 40L else 1000L; N_EXT <- if (SMOKE) 120L else 10000L; CAP_CPUH <- 300
OFF <- c(free25 = 601L, lib25 = 602L, free24 = 603L, lib24 = 604L, ext25 = 605L, ext24 = 606L)   # SEED + offset (3.2); +607 unused (p_lib is enumerated)
CLOSED <- c("PAUMTZ", "PAXMIJ")                                                                     # never read (2.3)
CHILD <- c("AML001", "AML003", "AML005", "AML006", "AML009", "AML012", "AML022", "AML025", "AML026", "AML028", "AML038", "AML043",
           "AML048", "AML051", "AML052", "AML055", "AML056", "AML073")                                  # 2.3
CHEN <- c("AML103", "AML162", "AML24", "AML304", "AML320", "NBM1", "NBM2", "NBM3", "NBM4")            # 2.3 (manifest order)
T15 <- paste0("Control", 1:5)                                                                       # P02 + P03 (4.4)
K2B_SAMPLES <- c("AML0048", "Control0004")                                                          # K2b re-aggregation (declared here)
stopifnot(length(CHILD) == 18L, !any(CLOSED %in% c(CHILD, CHEN)))

## known-answer bookkeeping (41's pattern; section 6's states plus "reported")
CHK <- list()
chk <- function(id, scope, name, expected, got, tol = NULL) {
  ok <- if (is.null(tol)) isTRUE(all.equal(expected, got)) else length(got) == length(expected) && isTRUE(all(abs(expected - got) <= tol))
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = scope, kind = "deterministic", check = name, expected = paste(format(expected), collapse = ";"),
                                         got = paste(format(got), collapse = ";"), state = if (ok) "PASS" else "FAIL_deterministic")
  message(sprintf("  %-4s %-84s expect %-14s got %-14s %s", id, name, paste(format(expected), collapse = ";"), paste(format(got), collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
rep_ <- function(id, scope, name, value, state = "REPORTED") {
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = scope, kind = "reported", check = name, expected = "", got = paste(value, collapse = ";"), state = state)
  message(sprintf("  %-4s %-84s %s [%s]", id, name, paste(value, collapse = ";"), state)) }
checks <- function() { f <- file.path(ODIR, paste0("ws_a1_g03_checks", SFX, ".csv")); if (!file.exists(f)) stop("run --known first: ", f); fread(f) }
gate <- function(sc) {   # a mode runs only when every deterministic check of its scope passed (section 6), under this script version
  C <- checks(); stopifnot(nrow(C[kind == "deterministic"]) > 0)   # sc, not scope: inside C[...] "scope" is the column
  if (!SMOKE && !identical(unique(C$script_md5), MD5)) stop("the checks table was written by another version of this script; rerun --known (section 6)")
  all(C[grepl(sc, scope) & kind == "deterministic", state] == "PASS") }
syn <- function(M) { attr(M, "synthetic") <- TRUE; M }           # marks a smoke matrix; need_syn() refuses anything else in a smoke run
need_syn <- function(M, what) if (SMOKE && !isTRUE(attr(M, "synthetic"))) stop("smoke run: a real ", what, " matrix would enter a computation")
mcl <- function(X, FUN) { r <- mclapply(X, FUN, mc.cores = NC)   # every worker must return; a dead worker stops the mode (code review 1)
  bad <- vapply(r, function(z) is.null(z) || inherits(z, "try-error"), TRUE); if (any(bad)) stop(sum(bad), " of ", length(r), " parallel workers failed"); r }
recheck_inputs <- function() {   # every mode after --known re-checks the md5s that --known recorded (code review 1)
  f <- rg("input_md5.csv"); if (!file.exists(f)) stop("run --known first: ", f); I <- fread(f)
  now <- unname(tools::md5sum(I$path)); bad <- I$path[is.na(now) | now != I$md5]
  if (length(bad)) stop("inputs changed since --known: ", paste(basename(bad), collapse = ", ")) }
cols <- function(M, j) { X <- M[, j, drop = FALSE]; if (isTRUE(attr(M, "synthetic"))) syn(X) else X }   # subsetting drops attributes; keep the smoke mark
ANNO_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "03_annotation_reconciled")   # config_hierarchy.R:51

## -- 1. S5 side: blacklist, the 25, P370 (27:28-49; 37:45-80) ------------------------------------------------------------------
BL <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_s1_ambient_blacklist.csv"))
GO <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
PBW <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"))
G <- PBW$genes
BLACK <- unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene],
                  "XIST", GO[chr == "chrY", gene], grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))   # 27:47-49
S5T <- fread(file.path(TABWS, "ws_s5_scores.csv"))[compartment == "T"]          # file order = 27's design order (AML first)
COV <- fread(file.path(TABWS, "ws_a1_batch_covariates.csv"))                    # library_id (September format) and period for the 25
IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))$gene
F402 <- IMP[, .(gene, s = sign(logFC))]; P370 <- F402[!gene %in% AMB]
S25 <- S5T$sample; S24 <- setdiff(S25, "AML0310")
stopifnot(length(S25) == 25L, length(S24) == 24L, nrow(F402) == 402L, nrow(P370) == 370L, all(S25 %in% COV$sample))
D25 <- S5T[, .(sample, status = factor(status, c("healthy", "AML")), hb, gran, sex = factor(sex), subtype_share)]
D25[, `:=`(library_id = COV$library_id[match(sample, COV$sample)], period = COV$period[match(sample, COV$sample)], dataset = "GSE185381")]
if (!SMOKE) { RAW_T <- PBW$raw$T[, S25]; RAW_NK <- PBW$raw$NK[, S25]; PBW_CELLS <- PBW$cells[sample %in% S25] }
rm(PBW); invisible(gc())                                                        # only the 25 columns are kept (2.3: closed columns dropped at load)

## -- 2. September side: manifest, pinned universe, T_NK columns of the rung samples (02:82-108; 01:132-141) -------------------------
.perm_seed <- function(label) {   # 03:74-78 verbatim (PREREG 8.4 of September)
  M <- 2147483647; h <- 0
  for (ch in utf8ToInt(label)) h <- (h * 131 + ch) %% M
  as.integer((h + SEED) %% M)
}
MAN <- fread(file.path(SEPT_TAB, "pb_sample_manifest.csv"))[in_analysis == TRUE]
RUNG_SAMPLES <- MAN[(dataset == "GSE185381" & sample %in% c(S24, CHILD)) | (dataset == "Chen2023" & sample %in% CHEN)]
stopifnot(nrow(RUNG_SAMPLES) == 51L, !any(RUNG_SAMPLES$sample %in% CLOSED))
sept_file <- function(ds, s) { stopifnot(!s %in% CLOSED); file.path(PB_RDS_DIR, ds, paste0(s, "__pseudobulk.rds")) }
dep_file  <- function(ds, s) { stopifnot(!s %in% CLOSED); file.path(SEPT_TAB, "depth", paste0(ds, "__", s, "__pb_depth.csv")) }
# the pinned universe (decision 9.6): 02:82-89 with the first Discovery sample of each dataset, named (AML103, AML001), minus sex genes
GO4 <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, col.names = c("gene", "chr", "start", "end"))
SEXG <- union(GO4[chr %in% c("chrY", "Y"), unique(gene)], "XIST")
FIRST <- MAN[split == "Discovery", .SD[1], by = dataset][, .(dataset, sample)]
stopifnot(identical(FIRST$sample, c("AML103", "AML001")))
GENES <- list(Chen2023 = readRDS(sept_file("Chen2023", "AML103"))$genes, GSE185381 = readRDS(sept_file("GSE185381", "AML001"))$genes)
UNIV <- setdiff(Reduce(intersect, GENES), SEXG)
# T_NK columns (own axis) and depth rows of the 51 rung samples; nothing else of these files is kept
load_sept <- function() {
  vec <- list(); dep <- list(); ax_ok <- logical(0)
  for (i in seq_len(nrow(RUNG_SAMPLES))) {
    ds <- RUNG_SAMPLES$dataset[i]; s <- RUNG_SAMPLES$sample[i]; o <- readRDS(sept_file(ds, s))
    ax_ok[s] <- identical(o$genes, GENES[[ds]]) && identical(rownames(o$counts), o$genes)
    vec[[s]] <- setNames(as.numeric(o$counts[, "T_NK"]), o$genes); rm(o)
    dep[[s]] <- fread(dep_file(ds, s))[hierarchy_bin == "T_NK"]
  }
  list(vec = vec, dep = rbindlist(dep), ax_ok = ax_ok) }

## -- 3. engines: one function for every rung (27:70-87 edgeR; 02:145-171 voom/blocked; 03:149-160 frozen) ---------------------------
fml <- function(terms) as.formula(paste("~", paste(terms, collapse = " + ")))
MAIN <- c("status", "hb", "gran", "sex", "subtype_share")
# spec: engine in edgeR / voom / voom_block; univ in S5 / SEP; terms. d: sample, status (healthy first), covariates, library_id, dataset
fit_rung <- function(cnt, d, engine, univ, terms, keep_obj = FALSE, U = UNIV) {
  need_syn(cnt, "count")
  if ("sex" %in% terms && uniqueN(d$sex) < 2) terms <- setdiff(terms, "sex")          # 37:52's guard
  d <- as.data.frame(d); X <- model.matrix(fml(terms), data = d)
  if (!"statusAML" %in% colnames(X) || qr(X)$rank < ncol(X)) stop("design not estimable")
  if (nrow(X) - ncol(X) < 2) stop("residual df < 2: not estimable")               # edgeR segfaults at 1 residual df (37:56)
  if (univ == "S5") {                                                              # 27:73-75: blacklist first, then the two filters
    M <- cnt[!rownames(cnt) %in% BLACK, d$sample, drop = FALSE]; y <- DGEList(M, group = d$status)
    keep <- filterByExpr(y, group = d$status) & (rowMeans(cpm(y) > 1) >= 0.5)
  } else {                                                                         # 02:153: the pinned base, then filterByExpr(design)
    M <- cnt[U, d$sample, drop = FALSE]; y <- DGEList(M); keep <- filterByExpr(y, X)
  }
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  cons <- NA_real_; v <- NULL; blk <- NULL
  if (engine == "edgeR") {                                                         # 27:78-80
    y <- estimateDisp(y, X, robust = TRUE); f <- glmQLFit(y, X, robust = TRUE); t <- glmQLFTest(f, coef = "statusAML")
    tt <- as.data.table(topTags(t, n = Inf, sort.by = "none")$table, keep.rownames = "gene")[, .(gene, logFC, P = PValue, q = FDR)]
  } else {
    if (engine == "voom") { v <- voom(y, X); fit <- lmFit(v, X) }                  # 02:166 (unblocked branch)
    else {                                                                         # 02:157-163 (two passes)
      blk <- factor(d$library_id)
      v  <- voom(y, X); c1 <- duplicateCorrelation(v, X, block = blk)
      v  <- voom(y, X, block = blk, correlation = c1$consensus)
      c2 <- duplicateCorrelation(v, X, block = blk); cons <- c2$consensus
      fit <- lmFit(v, X, block = blk, correlation = cons) }
    fit <- eBayes(fit)
    tt <- as.data.table(topTable(fit, coef = "statusAML", number = Inf, sort.by = "none"), keep.rownames = "gene")[, .(gene, logFC, P = P.Value, q = adj.P.Val)]
  }
  out <- list(tt = tt, n_genes = sum(keep), cons = cons)
  if (keep_obj) out$obj <- list(v = v, d = d, terms = terms, blk = blk, cons = cons, engine = engine)
  out }
# frozen refit (03:149-160): the observed voom object, gene set and consensus are reused; only lmFit/eBayes rerun under new labels
refit_frozen <- function(O, status_new, E = NULL) {
  d2 <- O$d; d2$status <- factor(status_new, levels = c("healthy", "AML")); X <- model.matrix(fml(O$terms), data = d2)
  if (!"statusAML" %in% colnames(X) || qr(X)$rank < ncol(X)) return(NULL)
  v <- O$v; if (!is.null(E)) v$E <- E
  fit <- if (O$engine == "voom_block") lmFit(v, X, block = O$blk, correlation = O$cons) else lmFit(v, X)
  as.data.table(topTable(eBayes(fit), coef = "statusAML", number = Inf, sort.by = "none"), keep.rownames = "gene")[, .(gene, logFC, P = P.Value, q = adj.P.Val)] }
# the four counts of 3.2 (restricted = outside S5's blacklist; q from the whole universe, then the genes dropped)
counts_of <- function(tt) {
  if (is.null(tt)) return(c(c_S5 = NA_integer_, c_Sep = NA_integer_, c_S5r = NA_integer_, c_Sepr = NA_integer_))
  s5 <- tt$q < Q_S5 & abs(tt$logFC) > L_S5; sp <- tt$q < Q_SEP & abs(tt$logFC) >= L_SEP; r <- !tt$gene %in% BLACK
  c(c_S5 = sum(s5, na.rm = TRUE), c_Sep = sum(sp, na.rm = TRUE), c_S5r = sum(s5 & r, na.rm = TRUE), c_Sepr = sum(sp & r, na.rm = TRUE)) }

## -- 4. permutation draws (3.2): stored as sample-keyed labellings, "A" = AML, "H" = healthy -------------------------------------------
lab_chr <- function(st) ifelse(as.character(st) == "AML", "A", "H")
draw_free <- function(d, n, seed) { set.seed(seed); st <- lab_chr(d$status)
  M <- t(vapply(seq_len(n), function(i) sample(st), character(nrow(d)))); colnames(M) <- d$sample; M }
lib_units <- function(d) {   # 03:116-117: arm-pure libraries are units; mixed libraries keep their real labels
  dd <- data.table(sample = d$sample, lab = lab_chr(d$status), lib = d$library_id)
  LP <- dd[, .(n = .N, n_arm = uniqueN(lab), lab1 = lab[1]), by = lib][, pure := n_arm == 1L]
  list(dd = dd, LP = LP) }
draw_lib <- function(d, n, seed) {   # 03:118-127 for one dataset: shuffle the unit labels of the pure libraries
  U <- lib_units(d); P <- U$LP[pure == TRUE]; set.seed(seed)
  M <- t(vapply(seq_len(n), function(i) { a <- U$dd$lab; nl <- sample(P$lab1); for (k in seq_len(nrow(P))) a[U$dd$lib == P$lib[k]] <- nl[k]; a }, character(nrow(d))))
  colnames(M) <- d$sample; M }
enum_lib <- function(d) {   # the exact library space (review 2, finding 6): every choice of which pure units are healthy
  U <- lib_units(d); P <- U$LP[pure == TRUE]; k <- sum(P$lab1 == "H"); cmb <- combn(nrow(P), k)
  M <- t(apply(cmb, 2, function(j) { a <- U$dd$lab; for (u in seq_len(nrow(P))) a[U$dd$lib == P$lib[u]] <- if (u %in% j) "H" else "A"; a }))
  colnames(M) <- d$sample; M }
n_units <- function(d) {   # 03:168-177 for the sample and library nulls (per dataset product)
  dd <- as.data.table(d)[, .(sample, lab = lab_chr(status), lib = library_id, dataset)]
  samp <- prod(dd[, .(u = choose(.N, sum(lab == "H"))), by = dataset]$u)
  LP <- dd[, .(n_arm = uniqueN(lab), lab1 = lab[1]), by = .(dataset, lib)][n_arm == 1L]
  lib <- prod(LP[, .(u = choose(.N, sum(lab1 == "H"))), by = dataset]$u)
  c(sample = samp, library = lib) }
status_from <- function(M, i, d) factor(ifelse(M[i, d$sample] == "A", "AML", "healthy"), levels = c("healthy", "AML"))   # joined by sample ID (K7)

## -- 5. smoke data: synthetic counts and covariates on the real IDs (no real count or covariate enters a smoke computation) ----------
if (SMOKE) {
  set.seed(SEED + 699L); SMS <- sample.int(.Machine$integer.max, 20L)   # synthetic-only seeds sit above +607
  set.seed(SMS[1]); ng <- length(G); mu <- rep(0, ng); on <- sample.int(ng, 12000L); mu[on] <- exp(rnorm(12000L, 0, 2))
  AMLS <- D25[status == "AML", sample]; PLANT <- on[1:300]                       # a planted shift (x 2.5) in the AML columns, smoke only
  synth <- function(ids, depth, plant = TRUE) { M <- vapply(ids, function(s) rnbinom(ng, size = 8, mu = mu * depth / sum(mu) * ifelse(plant & seq_len(ng) %in% PLANT & s %in% AMLS, 2.5, 1)), numeric(ng))
    dimnames(M) <- list(G, ids); syn(M) }
  RAW_T <- synth(S25, 2e6); RAW_NK <- synth(S25, 4e5)
  set.seed(SMS[2]); D25[, `:=`(hb = runif(.N, 0, 0.02), gran = runif(.N, 0, 0.01), sex = factor(sample(c("female", "male"), .N, TRUE), levels = c("female", "male")),
                              subtype_share = runif(.N, 0.2, 0.6), period = sample(period))]   # real level names; the real capture period is shuffled
}

## September's permutation procedure (03), used by K3, E1, E2, 4.5 and the R6 part of the extension sentence
lib_purity3 <- function(md) md[, .(n = .N, n_arm = uniqueN(arm), arm1 = as.character(arm)[1]), by = .(dataset, library_id)][, pure := n_arm == 1L][]   # 03:116-117
perm_library3 <- function(md, LP) {   # 03:118-127 verbatim
  a <- as.character(md$arm)
  for (d in unique(md$dataset)) { P <- LP[dataset == d & pure == TRUE]; if (nrow(P) < 2L) next
    newlab <- sample(P$arm1); for (k in seq_len(nrow(P))) a[md$dataset == d & md$library_id == P$library_id[k]] <- newlab[k] }
  factor(a, levels = levels(md$arm)) }
perm_sample3x <- function(md) { a <- as.character(md$arm); for (d in unique(md$dataset)) { i <- which(md$dataset == d); a[i] <- sample(a[i]) }; factor(a, levels = levels(md$arm)) }   # 03:108-112
ext_design <- function(rung) {   # E1 / E2 in manifest order (September's roster order, so .perm_seed gives September's procedure; code review 1)
  MAN[(dataset == "GSE185381" & sample %in% c(S24, CHILD)) | (rung == "E2" & dataset == "Chen2023" & sample %in% CHEN)][
    , .(sample, status = factor(arm, c("healthy", "AML")), library_id, dataset)] }
ext_counts <- function(d) {   # counts on the pinned universe; the stored September T_NK columns (smoke: synthetic)
  if (SMOKE) { set.seed(SMS[4]); M <- vapply(d$sample, function(s) rnbinom(length(UNIV), size = 8, mu = mu[match(UNIV, G)] * 1.5e6 / sum(mu) *
                ifelse(!SMOKE_B & s %in% AMLS & UNIV %in% G[PLANT], 2.5, 1)), numeric(length(UNIV))); dimnames(M) <- list(UNIV, d$sample); return(syn(M)) }
  M <- vapply(seq_len(nrow(d)), function(i) { o <- readRDS(sept_file(d$dataset[i], d$sample[i])); as.numeric(o$counts[UNIV, "T_NK"]) }, numeric(length(UNIV)))
  dimnames(M) <- list(UNIV, d$sample); M }
dep_rows <- function(d) rbindlist(lapply(seq_len(nrow(d)), function(i) fread(dep_file(d$dataset[i], d$sample[i]))[hierarchy_bin == "T_NK"]))
reportable <- function(base, dep) {   # 02:190-197: a hit in the baseline fit and in the fit with log10(median_umi)
  hb <- base[q < Q_SEP & abs(logFC) >= L_SEP, gene]; hd <- dep[q < Q_SEP & abs(logFC) >= L_SEP, gene]; h <- intersect(hb, hd)
  c(c_Sep_rep = length(h), c_Sepr_rep = sum(!h %in% BLACK)) }
marker_gate_T <- function(cnt, d) {   # 02:112-139, T_NK branch (amended rule), on the pinned universe
  need_syn(cnt, "count")
  cp <- edgeR::cpm(cnt[intersect(UNIV, rownames(cnt)), d$sample, drop = FALSE]); mk <- c("CD3E", "CD3D")
  rbindlist(lapply(c("healthy", "AML"), function(a) { mu1 <- rowMeans(cp[mk, d$status == a, drop = FALSE]); mu2 <- rowMeans(cp[mk, d$status != a, drop = FALSE])
    data.table(arm = a, marker = mk, mean_cpm = mu1, cpm_other_arm = mu2, pass = mu1 >= 100 & mu2 >= 100 & abs(log2(mu1 / mu2)) <= 1) }))[, arm_pass := any(pass), by = arm][] }
planted <- function(O, rung) {   # 03:202-224 with the rung's frozen voom object
  set.seed(.perm_seed(paste0("planted|g03_", rung, "|T_NK")))
  mu1 <- rowMeans(O$v$E); ter <- cut(rank(mu1, ties.method = "first"), 3, labels = FALSE)
  pg <- sample(rownames(O$v$E)[ter == 2L], min(200L, sum(ter == 2L))); aml <- O$d$status == "AML"
  cur <- rbindlist(lapply(c(0.5, 1, 1.5, 2, 3), function(L) { E2 <- O$v$E; E2[pg, aml] <- E2[pg, aml] + L
    tt <- refit_frozen(O, O$d$status, E = E2)
    data.table(rung = rung, planted_lfc = L, n_planted = length(pg), recovery_full = tt[gene %in% pg & q < Q_SEP & abs(logFC) >= L_SEP, .N] / length(pg),
               recovery_qonly = tt[gene %in% pg & q < Q_SEP, .N] / length(pg)) }))
  ok <- cur[recovery_full >= 0.80]; cur[, LFC80 := if (nrow(ok)) min(ok$planted_lfc) else NA_real_][] }
p_rows <- function(rung, obs, Z, nl, scheme) rbindlist(lapply(names(obs), function(cn) { z <- Z[, cn]; n <- sum(!is.na(z)); x <- sum(z >= obs[[cn]], na.rm = TRUE); s <- cp_state(x, n)
  data.table(rung = rung, count = cn, null = nl, scheme = scheme, obs = obs[[cn]], x = x, n = n, p = (1 + x) / (1 + n), lo = s$lo, hi = s$hi, state = s$state) }))
frozen_null <- function(O, md, nl, label) {   # 1,000 draws of September's procedure: labels drawn in sequence from one seed (03:180), refits in parallel
  set.seed(.perm_seed(label)); LP <- lib_purity3(md)
  A <- lapply(seq_len(N_DRAW), function(k) if (nl == "free") perm_sample3x(md) else perm_library3(md, LP))
  do.call(rbind, mcl(A, function(a) counts_of(tryCatch(refit_frozen(O, a), error = function(e) NULL)))) }

## 5.4's expected-calls function (the addendum's method: recovery curve anchored at (0, 0), linear, flat above 3; K9 checks it)
PW <- fread(file.path(SEPT_TAB, "power.csv"))[hierarchy_bin == "T_NK"]; setorder(PW, planted_lfc)
rec <- function(v) approx(c(0, PW$planted_lfc), c(0, PW$recovery_full), xout = v, rule = 2)$y
L370 <- abs(IMP[gene %in% P370$gene, logFC]); stopifnot(length(L370) == 370L)
exp_calls <- function(f) sum(rec(abs(f) * L370))   # |f| when f < 0: the curve is symmetric in sign (5.4)

## -- 6. --known (section 6) ------------------------------------------------------------------------------------------------------------
if (MODE == "known") {
  ts("[K0] md5 of the declared inputs")
  DECL <- data.table(path = c(file.path(SEPT_TAB, c("de/primary__T_NK.csv", "perm_null.csv", "perm_summary.csv", "hitlist_frozen.csv", "de_summary.csv",
                                                   "power.csv", "marker_gate.csv", "validation/T_NK.csv", "pb_sample_manifest.csv")),
                              file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds")),
                     md5 = c("73e73469d6a01e4a29dca1d3d2c6210a", "1abdc4bb449d3f3a8fd2482f2492ebb1", "dc03aae846d3c2af24a77f8fe3307698",
                             "8e3d082471ea33878abb00960dd736c6", "f4f0fc82ddad0179a4b16958cab0b507", "f23833718500d7861119cf73ba87ef63",
                             "f2c105102f425490b9a6fc97f7ddbc32", "ba2b152140d80916221940d7276ff72d", "2e36ce07c3883942d2334db3b26b3c09",
                             "6c2fb685c49d4a23ef20a6842aa7aa67"))
  got0 <- unname(tools::md5sum(DECL$path))
  chk("K0", "L;E", sprintf("md5 of %d declared files (inventory A section 5; ws_pb)", nrow(DECL)), nrow(DECL), sum(!is.na(got0) & got0 == DECL$md5))
  OTHER <- c(file.path(TABWS, c("ws_s5_scores.csv", "ws_s5_imprint_genes.csv", "ws_s5_de_T.csv", "ws_s8_t_imprint_ambient_prone.csv",
                                "ws_a1_g02_s5_perm_draws.csv", "ws_a1_batch_covariates.csv", "decontx/GSE185381/ws_s1_ambient_blacklist.csv")),
             file.path(SEPT_TAB, c("secondary_exhaustive_perm.csv", "de/secondary__LMPP_GMP.csv", "de/secondary__HSC_MPP.csv")),
             file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), here::here("scripts", "config", "config_paths.R"),
             file.path(TABWS, c("ws_a1_g04_primary.csv", "ws_a1_g04_children.csv")),                                   # 5.4
             file.path(QC_RDS_DIR, "GSE185381", paste0(K2B_SAMPLES, ".rds")),                                           # K2b
             file.path(ANNO_DIR, "GSE185381", paste0(S24, "__anno_percell.csv")),                                       # K2c
             mapply(sept_file, RUNG_SAMPLES$dataset, RUNG_SAMPLES$sample), mapply(dep_file, RUNG_SAMPLES$dataset, RUNG_SAMPLES$sample))
  m_oth <- unname(tools::md5sum(OTHER))
  chk("K0", "L;E", sprintf("the other %d inputs exist (md5 recorded for later modes)", length(OTHER)), length(OTHER), sum(!is.na(m_oth)))
  MD5IN <- data.table(path = c(DECL$path, OTHER), md5 = c(got0, m_oth))

  ts("[K1-K1b] R0 reproduces S5 (real data only)")
  if (!SMOKE) {
    r0 <- fit_rung(RAW_T, D25, "edgeR", "S5", MAIN)
    h0 <- r0$tt[q < Q_S5 & abs(logFC) > L_S5]
    chk("K1", "L", "R0: genes tested / hits / up / down", c(7379, 402, 171, 231), c(r0$n_genes, nrow(h0), sum(h0$logFC > 0), sum(h0$logFC < 0)))
    chk("K1", "L", "R0: the 402 hit genes equal the frozen S5 list (and so P370)", TRUE, setequal(h0$gene, F402$gene))
    D0 <- fread(file.path(TABWS, "ws_s5_de_T.csv")); jj <- match(r0$tt$gene, D0$gene)
    chk("K1", "L", "R0: max |logFC - ws_s5_de_T logFC|", 0, max(abs(r0$tt$logFC - D0$logFC[jj])), tol = 1e-6)
    SD <- fread(file.path(TABWS, "ws_a1_g02_s5_perm_draws.csv"))                 # positions = T rows of ws_s5_scores.csv in file order (K1b)
    stopifnot(nrow(SD) == 200L, all(nchar(SD$labels) == 25L))
    nh <- unlist(mcl(seq_len(nrow(SD)), function(i) {
      d <- copy(D25); d[, status := factor(ifelse(strsplit(SD$labels[i], "")[[1]] == "1", "AML", "healthy"), levels = c("healthy", "AML"))]
      tryCatch(counts_of(fit_rung(RAW_T, d, "edgeR", "S5", MAIN)$tt)[["c_S5"]], error = function(e) NA_integer_) }))
    chk("K1b", "L", "the edgeR pipeline on S5's 200 stored labellings reproduces every n_hits", 200L, sum(nh == SD$n_hits, na.rm = TRUE))
    chk("K1b", "L", "... including the four draws >= 402 (sorted n_hits)", c(1615, 1166, 707, 563), sort(nh[SD$n_hits >= 402], decreasing = TRUE))
  } else rep_("K1", "L", "skipped in a smoke run (S5's real counts)", "skipped")

  ts("[K2-K2d] September inputs, gene axes, cell sets, library IDs")
  if (!SMOKE) {
    SP <- load_sept()
    tot <- vapply(names(SP$vec), function(s) sum(SP$vec[[s]]), 1)
    dd <- SP$dep[match(names(SP$vec), sample)]
    chk("K2", "L;E", "51 rung samples: stored T_NK column sum equals the depth table's total_counts", 51L, sum(tot == dd$total_counts))
    chk("K2", "L;E", "51 rung samples: n_cells >= 30 and passes_30 agree and are TRUE", 51L, sum(dd$n_cells >= 30 & dd$passes_30 == TRUE))
    chk("K2", "L;E", "the pinned universe rebuilt from gene names (AML103, AML001; minus sex genes)", 21843L, length(UNIV))
    chk("K2", "L;E", "every rung file's gene axis equals its dataset's first-sample axis", 51L, sum(SP$ax_ok))
    gs <- RUNG_SAMPLES[dataset == "GSE185381", sample]
    chk("K2b", "L", "ws_pb's gene axis equals every GSE185381 September file's axis (42)", TRUE, identical(G, GENES$GSE185381))
    for (s in K2B_SAMPLES) {
      cells <- PBW_CELLS[sample == s & kept == TRUE & group4 %in% c("T", "NK"), cell]
      o <- readRDS(file.path(QC_RDS_DIR, "GSE185381", paste0(s, ".rds"))); X <- SeuratObject::LayerData(o, assay = "RNA", layer = "counts"); rm(o)
      stopifnot(identical(rownames(X), G), all(cells %in% colnames(X)))
      re <- as.numeric(Matrix::rowSums(X[, cells, drop = FALSE])); rm(X)
      chk("K2b", "L", sprintf("%s: R3's cells re-aggregated from 01_per_sample_qc equal $raw$T + $raw$NK (identical)", s), TRUE,
          identical(re, as.numeric(RAW_T[, s] + RAW_NK[, s]))) }
    k2c <- rbindlist(lapply(S24, function(s) {
      a <- fread(file.path(ANNO_DIR, "GSE185381", paste0(s, "__anno_percell.csv")),
                 select = c("cell", "hierarchy_bin", "in_ccc_graph", "high_error", "bmm_broad"))
      a[!nzchar(trimws(hierarchy_bin)), hierarchy_bin := NA_character_]; a[, `:=`(in_ccc_graph = as.logical(in_ccc_graph), high_error = as.logical(high_error))]
      c5 <- a[in_ccc_graph == TRUE & high_error == FALSE & hierarchy_bin == "T_NK" & !(bmm_broad %in% "Early Lymphoid"), cell]   # 01:139-141
      c3 <- PBW_CELLS[sample == s & kept == TRUE & group4 %in% c("T", "NK"), cell]
      data.table(sample = s, n_R3 = length(c3), n_R5 = length(c5), n_both = length(intersect(c3, c5))) }))
    rep_("K2c", "L", "R5's cell numbers equal the depth tables' n_cells (of 24; reported: September's join need only reach 0.999)", sum(k2c$n_R5 == dd$n_cells[match(k2c$sample, dd$sample)]))
    rep_("K2c", "L", "R3 / R5 / shared cells, summed over the 24", k2c[, c(sum(n_R3), sum(n_R5), sum(n_both))])
    rep_("K2c", "L", "per-sample shared share of R3 (min, median)", round(k2c[, c(min(n_both / n_R3), median(n_both / n_R3))], 3))
  } else rep_("K2", "L;E", "skipped in a smoke run (real files)", "skipped")
  lm <- MAN[dataset == "GSE185381" & sample %in% S24]
  chk("K2d", "L", "on the 24, library_id of ws_a1_batch_covariates equals pb_sample_manifest's", 24L,
      sum(D25$library_id[match(lm$sample, D25$sample)] == lm$library_id))

  ts("[K3-K3c] the production engine (fit_rung, refit_frozen, frozen_null) on September's Chen2023 secondary bins (real data only)")
  if (!SMOKE) {
    CH <- MAN[split == "Discovery" & dataset == "Chen2023"]                              # 03's roster, in manifest order
    DEP3 <- rbindlist(lapply(CH$sample, function(s) fread(dep_file("Chen2023", s))))
    PB3 <- lapply(CH$sample, function(s) readRDS(sept_file("Chen2023", s))$counts); names(PB3) <- CH$sample
    SEC <- c("LMPP_GMP", "HSC_MPP")                                                       # only the secondary rows are kept (code review 1)
    PN <- fread(file.path(SEPT_TAB, "perm_null.csv"))[hierarchy_bin %in% SEC]; DS <- fread(file.path(SEPT_TAB, "de_summary.csv"))[hierarchy_bin %in% SEC]
    for (bin in SEC) {
      keep <- DEP3[hierarchy_bin == bin & n_cells >= 30, .(dataset, sample)][CH, on = c("dataset", "sample"), nomatch = 0]   # 03:96-98
      d <- keep[, .(sample, status = factor(arm, c("healthy", "AML")), library_id, dataset)]
      cnt <- vapply(d$sample, function(s) PB3[[s]][UNIV, bin], numeric(length(UNIV))); dimnames(cnt) <- list(UNIV, d$sample)
      f <- fit_rung(cnt, d, "voom", "SEP", "status", keep_obj = TRUE); obs <- counts_of(f$tt)[["c_Sep"]]
      ref <- fread(file.path(SEPT_TAB, "de", paste0("secondary__", bin, ".csv")))
      chk("K3", "L;E", sprintf("%s: n (AML v healthy), genes tested, baseline hits = de_summary", bin),
          as.numeric(DS[hierarchy_bin == bin, c(n_AML, n_healthy, n_genes_tested, hits_baseline)]), as.numeric(c(sum(d$status == "AML"), sum(d$status == "healthy"), f$n_genes, obs)))
      chk("K3", "L;E", sprintf("%s: max |logFC - de/secondary__%s.csv|", bin, bin), 0, max(abs(f$tt$logFC - ref$logFC[match(f$tt$gene, ref$gene)])), tol = 1e-8)
      Z <- frozen_null(f$obj, data.table(arm = d$status, dataset = d$dataset, library_id = d$library_id), "free", paste0("perm|sample_within_dataset|secondary|", bin))
      chk("K3", "L;E", sprintf("%s: frozen_null's 1,000 sample-null draws reproduce perm_null.csv's n_hits", bin), 1000L,
          sum(Z[, "c_Sep"] == PN[hierarchy_bin == bin & null == "sample_within_dataset"][order(perm), n_hits], na.rm = TRUE))
      if (bin == "LMPP_GMP") {   # K3c: 07:140-159, the full refit over the 126 arrangements, through fit_rung
        EX <- fread(file.path(SEPT_TAB, "secondary_exhaustive_perm.csv"))[hierarchy_bin == bin & scheme == "full_refit"]
        n <- nrow(d); k <- sum(d$status == "healthy"); cmb <- combn(n, k)
        hv <- unlist(mcl(seq_len(ncol(cmb)), function(j) { a <- rep("AML", n); a[cmb[, j]] <- "healthy"; dj <- copy(d)[, status := factor(a, c("healthy", "AML"))]
          counts_of(fit_rung(cnt, dj, "voom", "SEP", "status")$tt)[["c_Sep"]] }))
        chk("K3c", "L", "LMPP_GMP full refit (fit_rung) over 126 labellings: observed_hits, n_ge_observed, null_max", as.numeric(EX[, c(observed_hits, n_ge_observed, null_max)]),
            as.numeric(c(obs, sum(hv >= obs), max(hv)))) }
    }
    rm(PB3)
  } else rep_("K3", "L;E", "skipped in a smoke run (September's real files)", "skipped")

  ts("[K4] the blocked fit and the frozen refit equal September's code (synthetic counts, shared libraries)")
  set.seed(SEED + 698L)
  ns <- 16L; Y <- matrix(rnbinom(3000L * ns, size = 6, mu = rep(exp(rnorm(3000L, 3, 1.5)), ns)), 3000L, ns, dimnames = list(c(UNIV[1:2990], paste0("synG", 1:10)), paste0("s", 1:ns)))
  mdk <- data.table(dataset = factor(rep(c("dA", "dB"), each = 8)), arm = relevel(factor(rep(rep(c("AML", "healthy"), each = 4), 2)), ref = "healthy"),
                    library_id = factor(c("L1", "L1", "L2", "L3", "L4", "L4", "L5", "L6", "L7", "L7", "L8", "L9", "L10", "L11", "L11", "L12")))
  fit_bin_02 <- function(B) {   # 02:145-171, primary tier, verbatim apart from the return value
    md <- B$meta; design <- model.matrix(as.formula("~ dataset + arm"), data = md)
    y <- DGEList(B$counts); keep <- filterByExpr(y, design); y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM")
    v  <- voom(y, design); c1 <- duplicateCorrelation(v, design, block = md$library_id)
    v  <- voom(y, design, block = md$library_id, correlation = c1$consensus); c2 <- duplicateCorrelation(v, design, block = md$library_id)
    fit <- eBayes(lmFit(v, design, block = md$library_id, correlation = c2$consensus))
    list(tt = as.data.table(topTable(fit, coef = "armAML", number = Inf), keep.rownames = "gene"), cons = c2$consensus, v = v, design = design, md = md) }
  Yk <- Y[intersect(rownames(Y), UNIV), ]
  ref02 <- fit_bin_02(list(counts = Yk, meta = mdk))
  dk <- data.table(sample = colnames(Y), status = factor(as.character(mdk$arm), levels = c("healthy", "AML")), dataset = mdk$dataset, library_id = as.character(mdk$library_id))
  mine <- fit_rung(syn(Yk), dk, "voom_block", "SEP", c("dataset", "status"), keep_obj = TRUE, U = rownames(Yk))
  jj <- match(ref02$tt$gene, mine$tt$gene)
  chk("K4", "L;E", "blocked fit: max |logFC diff|, max |P diff|, |consensus diff| vs 02's fit_bin", c(0, 0, 0),
      c(max(abs(ref02$tt$logFC - mine$tt$logFC[jj])), max(abs(ref02$tt$P.Value - mine$tt$P[jj])), abs(ref02$cons - mine$cons)), tol = c(1e-10, 1e-12, 1e-12))
  set.seed(SEED + 697L); ap <- sample(as.character(mdk$arm))
  refit_03 <- function(O, md, arm_new) {   # 03:149-160 verbatim (primary tier)
    md2 <- copy(md)[, arm := arm_new]; design <- model.matrix(as.formula("~ dataset + arm"), data = md2)
    fit <- eBayes(lmFit(O$v, design, block = md$library_id, correlation = O$cons))
    as.data.table(topTable(fit, coef = "armAML", number = Inf), keep.rownames = "gene") }
  r03 <- refit_03(ref02, mdk, relevel(factor(ap), ref = "healthy")); rme <- refit_frozen(mine$obj, ap)
  chk("K4", "L;E", "frozen refit: max |logFC diff| vs 03's refit_hits under one permuted labelling", 0, max(abs(r03$logFC - rme$logFC[match(r03$gene, rme$gene)])), tol = 1e-10)

  ts("[K5, K7] permutable units; the stored draws")
  D24 <- D25[sample %in% S24]
  u25 <- n_units(D25); u24 <- n_units(D24)
  chk("K5", "L", "25: sample / library units", c(3268760, 8008), unname(u25))
  chk("K5", "L", "24: sample / library units", c(1961256, 8008), unname(u24))
  mE <- MAN[sample %in% c(S24, CHILD, CHEN)][, .(sample, status = factor(arm, c("healthy", "AML")), library_id, dataset)]
  chk("K5", "E", "E1 (42): sample / library units", c(choose(42, 10), 3876), unname(n_units(mE[dataset == "GSE185381"])))
  chk("K5", "E", "E2 (51): sample / library units", c(choose(42, 10) * choose(9, 4), 488376), unname(n_units(mE)))
  DR <- list(free25 = draw_free(D25, N_DRAW, SEED + OFF[["free25"]]), lib25 = draw_lib(D25, N_DRAW, SEED + OFF[["lib25"]]),
             free24 = draw_free(D24, N_DRAW, SEED + OFF[["free24"]]), lib24 = draw_lib(D24, N_DRAW, SEED + OFF[["lib24"]]),
             ext25 = draw_free(D25, N_EXT - N_DRAW, SEED + OFF[["ext25"]]), ext24 = draw_free(D24, N_EXT - N_DRAW, SEED + OFF[["ext24"]]),
             enum25 = enum_lib(D25), enum24 = enum_lib(D24))
  chk("K5", "L", "enumerated library spaces (25, 24)", c(8008, 8008), c(nrow(DR$enum25), nrow(DR$enum24)))
  if (SMOKE) { DR$enum25 <- DR$enum25[1:60, ]; DR$enum24 <- DR$enum24[1:60, ] }   # smoke: 60 labellings, so the exact path runs in minutes
  chk("K7", "L", "every free draw keeps the arm sizes (15/10 and 14/10)", TRUE,
      all(rowSums(DR$free25 == "A") == 15L, rowSums(DR$ext25 == "A") == 15L, rowSums(DR$free24 == "A") == 14L, rowSums(DR$ext24 == "A") == 14L))
  fixed_ok <- function(M, d) { U <- lib_units(d); mix <- U$LP[pure == FALSE, lib]; i <- d$library_id %in% mix; all(M[, d$sample[i]] == matrix(lab_chr(d$status[i]), nrow(M), sum(i), byrow = TRUE)) }
  chk("K7", "L", "library draws hold the mixed library (P05) at its real labels, whole libraries move together", TRUE,
      fixed_ok(DR$lib25, D25) && fixed_ok(DR$lib24, D24) && fixed_ok(DR$enum24, D24) &&
      all(vapply(lib_units(D24)$LP[pure == TRUE, lib], function(l) all(apply(DR$lib24[, D24[library_id == l, sample], drop = FALSE], 1, uniqueN) == 1L), TRUE)))
  dperm <- D24[sample(.N)]                                                        # any design order: the labels are joined by sample ID
  chk("K7", "L", "labels joined by sample ID are order-free (draw 1, shuffled design order)", TRUE,
      identical(as.character(status_from(DR$free24, 1, dperm)), unname(ifelse(DR$free24[1, dperm$sample] == "A", "AML", "healthy"))))
  for (nm in names(DR)) fwrite(data.table(draw = seq_len(nrow(DR[[nm]])), DR[[nm]]), wg(paste0("draws_", nm, ".csv")))
  rep_("K7", "L", "draw sets stored (rows)", paste(names(DR), vapply(DR, nrow, 1), sep = "=", collapse = " "))

  ts("[K8] seconds per refit (one draw per engine; timing only, no count is written)")
  d1 <- copy(D24); d1[, status := status_from(DR$free24, 1, d1)]
  cntT <- cols(RAW_T, S24); cntTN <- cols(RAW_T, S24) + cols(RAW_NK, S24); if (SMOKE) cntTN <- syn(cntTN)
  tm <- function(expr) as.numeric(system.time(expr)[["elapsed"]])
  t_edger <- tm(fit_rung(cntT, d1, "edgeR", "S5", "status")); t_voom <- tm(fit_rung(cntTN, d1, "voom", "S5", "status"))
  t_blk <- tm(ob <- fit_rung(cntTN, d1, "voom_block", "S5", "status", keep_obj = TRUE)); t_frz <- tm(refit_frozen(ob$obj, d1$status))
  s24 <- if (SMOKE) { set.seed(SMS[3]); synth(S24, 1.5e6) } else vapply(S24, function(s) { o <- readRDS(sept_file("GSE185381", s)); as.numeric(o$counts[, "T_NK"]) }, numeric(length(G)))
  if (!SMOKE) dimnames(s24) <- list(G, S24)
  t_blk6 <- tm(fit_rung(s24, d1, "voom_block", "SEP", "status"))                    # R6's universe tests more genes (code review 1)
  rep_("K8", "L", "seconds per refit: edgeR / voom / blocked (two passes) / frozen / blocked on R6's universe", round(c(t_edger, t_voom, t_blk, t_frz, t_blk6), 2))

  ts("[K9] 5.4's method reproduces the addendum's prediction rows (power.csv and S5's own sizes only)")
  addm <- vapply(c(1, 1 / 1.43, 0.5, 14 / 39, 14 / 39 / 1.43), exp_calls, 1)
  chk("K9", "R", "expected calls: all / shrunk / halved / adults only / adults only shrunk", c(155.9, 63.3, 16.8, 3.8, 1.2), round(addm, 1), tol = rep(0.05, 5))
  CK <- rbindlist(CHK); put(CK, "checks")
  fwrite(MD5IN, wg("input_md5.csv"))
  ts(sprintf("known answers: %d rows; deterministic FAIL %d", nrow(CK), sum(CK$state == "FAIL_deterministic")))
}


## -- 7. the ladder (section 3): rung table, data per rung, chunked null draws -------------------------------------------------------
RUNG <- data.table(id = c("R0", "R1", "R2", "R3", "R4a", "R4b", "R5", "R6", "R0f", "R1f"),
                   set = c("25", rep("24", 7), "25", "24"), cells = c("T", "T", "T", "TNK", "TNK", "TNK", "SEP", "SEP", "T", "T"),
                   terms = c("MAIN", "MAIN", rep("status", 6), "MAIN", "MAIN"), engine = c(rep("edgeR", 4), "voom", rep("voom_block", 3), "edgeR", "edgeR"),
                   univ = c(rep("S5", 7), "SEP", "S5", "S5"), frozen = c(rep(FALSE, 5), TRUE, FALSE, TRUE, FALSE, FALSE),
                   class = c("anchor", "validity", "specification", "specification", "specification", "validity", "specification", "specification", "sens 4.2", "sens 4.2"))
setkey(RUNG, id)
LADDER <- c("R0", "R1", "R2", "R3", "R4a", "R4b", "R5", "R6")
D24 <- D25[sample %in% S24]                                                       # S5 order; every rung joins labels by sample ID
CHUNK <- if (SMOKE) c(edgeR = 10L, voom = 20L, voom_block = 10L) else c(edgeR = 250L, voom = 500L, voom_block = 50L)
SEPT24 <- NULL
rung_data <- function(id) {
  r <- RUNG[id]; d <- copy(if (r$set == "25") D25 else D24)
  if (id %in% c("R0f", "R1f")) { stopifnot(identical(levels(d$sex), c("female", "male")), d[sample == "Control4", sex] == "male")   # 4.2 (code review 1, B1)
    d[sample == "Control4", sex := factor("female", levels = levels(sex))] }
  cnt <- switch(r$cells, T = cols(RAW_T, d$sample), TNK = { x <- cols(RAW_T, d$sample) + cols(RAW_NK, d$sample); if (SMOKE) syn(x) else x },
                SEP = { if (is.null(SEPT24)) SEPT24 <<- sept24(); SEPT24 })
  list(cnt = cnt, d = d, terms = if (r$terms == "MAIN") MAIN else "status") }
sept24 <- function() {   # R5/R6 counts: the stored September T_NK columns of the 24 (GSE185381 axis = G, K2b)
  if (SMOKE) { set.seed(SMS[3]); return(synth(S24, 1.5e6, plant = !SMOKE_B)) }
  M <- vapply(S24, function(s) { o <- readRDS(sept_file("GSE185381", s)); stopifnot(identical(o$genes, G)); as.numeric(o$counts[, "T_NK"]) }, numeric(length(G)))
  dimnames(M) <- list(G, S24); M }
draw_set <- function(id, null, round) {   # which stored draw set a (rung, null, round) reads
  s <- RUNG[id]$set; if (round == 1L) paste0(null, s) else if (null == "free") paste0("ext", s) else paste0("enum", s) }
read_draws <- function(nm) { x <- fread(rg(paste0("draws_", nm, ".csv"))); M <- as.matrix(x[, -1]); rownames(M) <- x$draw; M }
lab_key <- function(M) apply(M, 1, paste, collapse = "")
tasks_round1 <- function() rbindlist(lapply(RUNG$id, function(id) rbindlist(lapply(c("free", "lib"), function(nl) {
  st <- seq(1L, N_DRAW, by = CHUNK[[RUNG[id]$engine]]); data.table(round = 1L, rung = id, null = nl, from = st, to = pmin(st + CHUNK[[RUNG[id]$engine]] - 1L, N_DRAW)) }))))
tasks_round2 <- function() {
  P <- fread(rg("ext_plan.csv"))[funded == TRUE]
  if (!nrow(P)) return(data.table(round = integer(0), rung = character(0), null = character(0), from = integer(0), to = integer(0)))
  rbindlist(lapply(seq_len(nrow(P)), function(i) { n <- P$n_rows[i]; ch <- CHUNK[[RUNG[P$rung[i]]$engine]]; st <- seq(1L, n, by = ch)
    data.table(round = 2L, rung = P$rung[i], null = P$null[i], from = st, to = pmin(st + ch - 1L, n)) })) }
ext_rows <- function(id, null) {   # round 2: the extension draws, or the enumerated labellings not already fitted in round 1
  M <- read_draws(draw_set(id, null, 2L)); if (null == "free") return(M)
  M[!lab_key(M) %in% lab_key(read_draws(draw_set(id, null, 1L))[, colnames(M), drop = FALSE]), , drop = FALSE] }
chunk_file <- function(round, id, null, from) file.path(WDIR, "null", sprintf("r%d_%s_%s_%05d%s.csv", round, id, null, from, SFX))

gate_pass <- function(sc) { C <- checks(); all(C[grepl(sc, scope) & kind == "deterministic", state] == "PASS") }
need_tables <- function(extra = character(0)) {   # the tables a gate outcome requires before --read / --combine run (code review 1)
  L <- gate_pass("L"); E <- gate_pass("E")
  req <- c(if (L) c("checks_obs", "p370", "ladder"), if (E) "extend", if (L && E) "sens", extra)
  miss <- req[!file.exists(file.path(ODIR, paste0("ws_a1_g03_", req, SFX, ".csv")))]
  if (length(miss)) stop("run these modes first: ", paste(miss, collapse = ", "))
  md5s <- unique(unlist(lapply(c("checks", req), function(x) unique(rtab(x)$script_md5))))
  if (!SMOKE && !identical(md5s, MD5)) stop("tables from different script versions: ", paste(md5s, collapse = " "))
  list(L = L, E = E) }
load_obs <- function() { O <- readRDS(rg("obs_fits.rds")); if (!SMOKE && !identical(attr(O, "script_md5"), MD5)) stop("obs_fits.rds was written by another script version"); O }
if (MODE == "ladder-obs") {
  if (!gate("L")) stop("section 6: a scope-L deterministic known answer has not passed; the ladder is 'not decisive'")
  recheck_inputs()
  ts("[3.1] observed fits of R0-R6 (and 4.2's R0f, R1f)")
  OBS <- list()
  for (id in RUNG$id) { x <- rung_data(id); r <- RUNG[id]
    f <- fit_rung(x$cnt, x$d, r$engine, r$univ, x$terms, keep_obj = r$frozen)
    OBS[[id]] <- list(tt = f$tt, n_genes = f$n_genes, cons = f$cons, counts = counts_of(f$tt), obj = f$obj)
    ts(sprintf("    %-4s genes %5d | consensus %s | c_S5 %d c_Sep %d c_S5' %d c_Sep' %d", id, f$n_genes, if (is.na(f$cons)) "-" else sprintf("%.3f", f$cons),
               OBS[[id]]$counts[[1]], OBS[[id]]$counts[[2]], OBS[[id]]$counts[[3]], OBS[[id]]$counts[[4]])) }
  ts("[K6] each rung changes the per-gene result")
  K6 <- rbindlist(lapply(2:length(LADDER), function(k) { a <- OBS[[LADDER[k - 1]]]$tt; b <- OBS[[LADDER[k]]]$tt
    same_set <- setequal(a$gene, b$gene); dif <- if (same_set) max(abs(a$logFC - b$logFC[match(a$gene, b$gene)])) else Inf
    data.table(id = "K6", scope = "L", check = sprintf("%s -> %s: per-gene logFC differs", LADDER[k - 1], LADDER[k]), got = dif, state = if (dif > 1e-6) "PASS" else "FAIL_deterministic") }))
  print(K6[, .(check, got, state)])
  ts("[5.1, 5.2] P370 and AC007952.4 per rung (5.1 at R0-R6 is circular)")
  MAINFIT <- fread(file.path(TABWS, "ws_s5_de_T.csv"))[, .(gene, lfc0 = logFC)]; if (SMOKE) MAINFIT[, lfc0 := 1]   # smoke: no real S5 logFC
  P51 <- rbindlist(lapply(RUNG$id, function(id) { tt <- OBS[[id]]$tt; x <- merge(merge(P370, tt, by = "gene"), MAINFIT, by = "gene")
    data.table(section = "5.1", rung = id, n_genes = OBS[[id]]$n_genes, consensus = OBS[[id]]$cons, n_tested = nrow(x), hits_S5 = x[q < Q_S5 & abs(logFC) > L_S5, .N], hits_Sep = x[q < Q_SEP & abs(logFC) >= L_SEP, .N],
               sign_share = mean(sign(x$logFC) == x$s), attenuation = median(x$s * x$logFC / abs(x$lfc0)), circular = TRUE) }))
  P52 <- rbindlist(lapply(RUNG$id, function(id) { tt <- OBS[[id]]$tt[gene == "AC007952.4"]
    data.table(section = "5.2", rung = id, in_universe = nrow(tt) == 1L, logFC = if (nrow(tt)) tt$logFC else NA_real_, q = if (nrow(tt)) tt$q else NA_real_) }))
  ts("[5.8] references for R1 (leave one other AML out of R0) and R2 (one covariate dropped from R1)")
  x0 <- rung_data("R0"); x1 <- rung_data("R1")
  ref1 <- rbindlist(mcl(setdiff(D25[status == "AML", sample], "AML0310"), function(s) { d <- x0$d[sample != s]
    data.table(section = "5.8", rung = "R1 ref", dropped = s, c_S5 = tryCatch(counts_of(fit_rung(x0$cnt, d, "edgeR", "S5", MAIN)$tt)[["c_S5"]], error = function(e) NA_integer_)) }))
  ref2 <- rbindlist(lapply(c("hb", "gran", "sex", "subtype_share"), function(v)
    data.table(section = "5.8", rung = "R2 ref", dropped = v, c_S5 = counts_of(fit_rung(x1$cnt, x1$d, "edgeR", "S5", setdiff(MAIN, v))$tt)[["c_S5"]])))
  attr(OBS, "script_md5") <- MD5; saveRDS(OBS, wg("obs_fits.rds"))
  put(K6[, .(id, scope, kind = "deterministic", check, expected = "> 1e-6", got = as.character(got), state)], "checks_obs")
  put(rbindlist(list(P51, P52, ref1, ref2), fill = TRUE), "p370")
  ts(sprintf("ladder-obs: K6 FAIL %d; tasks for --ladder-perm round 1: %d", sum(K6$state != "PASS"), nrow(tasks_round1())))
}

if (MODE == "ladder-perm") {
  if (!gate("L")) stop("section 6: a scope-L deterministic known answer has not passed")
  recheck_inputs()
  TK <- if (ROUND == 1L) tasks_round1() else tasks_round2()
  if (is.na(TASK) || TASK < 1L || TASK > nrow(TK)) stop(sprintf("--task=N with N in 1..%d for round %d", nrow(TK), ROUND))
  tk <- TK[TASK]; id <- tk$rung; r <- RUNG[id]; out <- chunk_file(ROUND, id, tk$null, tk$from)
  if (file.exists(out) && !SMOKE) stop("refusing to overwrite ", out)
  ts(sprintf("[3.2] round %d task %d/%d: %s, %s null, rows %d-%d", ROUND, TASK, nrow(TK), id, tk$null, tk$from, tk$to))
  M <- if (ROUND == 1L) read_draws(draw_set(id, tk$null, 1L)) else ext_rows(id, tk$null)
  x <- rung_data(id); OBS <- load_obs(); O <- OBS[[id]]$obj
  res <- rbindlist(mcl(tk$from:tk$to, function(i) {
    st <- status_from(M, i, x$d); d <- copy(x$d); d[, status := st]
    full <- tryCatch(counts_of(fit_rung(x$cnt, d, r$engine, r$univ, x$terms)$tt), error = function(e) counts_of(NULL))
    row <- data.table(row = i, lab = paste(M[i, x$d$sample], collapse = ""), t(full))   # "key" would be read as data.table()'s key argument
    if (r$frozen) { fz <- counts_of(tryCatch(refit_frozen(O, st), error = function(e) NULL)); row <- cbind(row, t(setNames(fz, paste0(names(fz), "_frozen")))) }
    row }))
  if (nrow(res) != tk$to - tk$from + 1L) stop("chunk incomplete: ", nrow(res), " rows")
  res[, script_md5 := MD5]; fwrite(res, out)
  ts(sprintf("written %s (%d rows, %d failed fits)", out, nrow(res), sum(is.na(res$c_S5))))
}

## rung states (3.2): Clopper-Pearson for Monte Carlo draws; exact for the enumerated library space
cp_state <- function(x, n) { if (n < 1L) return(list(lo = NA, hi = NA, state = "not computed"))
  ci <- binom.test(x, n)$conf.int; list(lo = ci[1], hi = ci[2], state = if (ci[2] <= ALPHA) "held" else if (ci[1] > ALPHA) "lost" else "borderline") }
worse <- function(s) { o <- c("held", "borderline", "lost", "not computed"); s[!s %in% o] <- "not computed"; if (!length(s)) "not computed" else o[max(match(s, o))] }   # an unknown state is the worst
deciding <- function(id) if (id == "R6") "c_S5r" else "c_S5"
read_chunks <- function(round) { fs <- list.files(file.path(WDIR, "null"), pattern = sprintf("^r%d_[A-Za-z0-9]+_(free|lib)_[0-9]{5}%s\\.csv$", round, SFX), full.names = TRUE)
  if (!length(fs)) return(NULL)
  rbindlist(lapply(fs, function(f) { m <- regmatches(basename(f), regexec("^r(\\d)_([A-Za-z0-9]+)_(free|lib)_", basename(f)))[[1]]
    cbind(data.table(round = as.integer(m[2]), rung = m[3], null = m[4]), fread(f)) }), fill = TRUE) }

if (MODE == "ladder-state") {
  if (!gate("L")) stop("section 6: a scope-L deterministic known answer has not passed"); recheck_inputs()
  OBS <- load_obs(); N1 <- read_chunks(1L)
  if (!SMOKE && !all(N1$script_md5 == MD5)) stop("round-1 chunks from another script version")
  exp1 <- tasks_round1()[, .(n = sum(to - from + 1L)), by = .(rung, null)]
  got1 <- N1[, .N, by = .(rung, null)]; miss <- merge(exp1, got1, all.x = TRUE)[is.na(N) | N != n]
  if (nrow(miss)) stop("round 1 is incomplete: ", paste(miss$rung, miss$null, collapse = "; "))
  N2 <- if (ROUND == 2L) read_chunks(2L) else NULL
  if (!is.null(N2) && !SMOKE && !all(N2$script_md5 == MD5)) stop("round-2 chunks from another script version")
  K8 <- checks()[id == "K8", as.numeric(strsplit(got, ";")[[1]])]; names(K8) <- c("edgeR", "voom", "voom_block", "frozen", "voom_block6")
  sec_per <- function(id) K8[[if (id == "R6") "voom_block6" else RUNG[id]$engine]]
  # K7 (cross-rung): every rung that shares a draw set assigns the same labels to each draw (code review 1)
  K7B <- N1[, .(n_labels = uniqueN(lab)), by = .(set = vapply(rung, function(r) RUNG[r]$set, ""), null, row)]
  K7B <- data.table(id = "K7b", check = "each draw carries one labelling across the rungs that share its draw set", got = sum(K7B$n_labels != 1L),
                    state = if (all(K7B$n_labels == 1L)) "PASS" else "FAIL_deterministic")
  print(K7B)
  pstate <- function(id, nl, cnt, scheme = "full") {   # one (rung, count, null) at the current round; inside N1[...] "null" is the column, nl the argument
    col <- if (scheme == "full") cnt else paste0(cnt, "_frozen"); obs <- OBS[[id]]$counts[[cnt]]
    z <- N1[rung == id & null == nl, ..col][[1]]; ext <- FALSE; exact <- FALSE
    if (ROUND == 2L && scheme == "full" && cnt == deciding(id) && !is.null(N2) && nrow(N2[rung == id & null == nl])) {
      z2 <- N2[rung == id & null == nl, ..col][[1]]
      if (nl == "free") { z <- c(z, z2); ext <- TRUE }
      else { E <- read_draws(draw_set(id, "lib", 2L)); k1 <- N1[rung == id & null == "lib"]; k2 <- N2[rung == id & null == "lib"]
             ix <- match(lab_key(E[, if (RUNG[id]$set == "25") S25 else S24, drop = FALSE]), c(k1$lab, k2$lab))
             if (anyNA(ix)) stop("enumeration incomplete for ", id); z <- c(k1[[col]], k2[[col]])[ix]; exact <- TRUE } }
    n <- sum(!is.na(z)); x <- sum(z >= obs, na.rm = TRUE); nf <- sum(is.na(z))   # failed fits are excluded from the valid draws and counted
    if (exact) { p <- x / n; return(data.table(rung = id, count = cnt, null = nl, scheme, obs, x, n, n_failed = nf, p, lo = p, hi = p, state = if (p <= ALPHA) "held" else "lost", draws = "exact enumeration")) }
    s <- cp_state(x, n)
    data.table(rung = id, count = cnt, null = nl, scheme, obs, x, n, n_failed = nf, p = (1 + x) / (1 + n), lo = s$lo, hi = s$hi, state = s$state, draws = if (ext) sprintf("%d + %d", N_DRAW, N_EXT - N_DRAW) else as.character(N_DRAW)) }
  ROWS <- rbindlist(lapply(RUNG$id, function(id) rbindlist(lapply(c("free", "lib"), function(nl)
    rbindlist(c(lapply(c("c_S5", "c_Sep", "c_S5r", "c_Sepr"), function(cn) pstate(id, nl, cn)),
                if (RUNG[id]$frozen) lapply(c("c_S5", "c_Sep", "c_S5r", "c_Sepr"), function(cn) pstate(id, nl, cn, "frozen"))))))))
  ST <- ROWS[scheme == "full" & count == vapply(rung, deciding, ""), .(state = worse(state), free = state[null == "free"], lib = state[null == "lib"]), by = rung]
  if (ROUND == 1L) {
    ts("[3.2] round-1 states and the extension plan")
    pri <- LADDER   # 3.2: rungs up to and including the first lost rung, then the later rungs in order = rung order
    force <- if (SMOKE) c("R1 free", "R1 lib", "R4a free") else character(0)   # smoke only: exercise both extension paths
    PL <- rbindlist(lapply(pri, function(id) { s <- ST[rung == id]; out <- list()
      if ((s$free == "borderline" && s$lib != "lost") || paste(id, "free") %in% force) out$f <- data.table(rung = id, null = "free", n_rows = N_EXT - N_DRAW)
      if ((s$lib == "borderline" && s$free != "lost") || paste(id, "lib") %in% force) out$l <- data.table(rung = id, null = "lib", n_rows = nrow(ext_rows(id, "lib")))
      rbindlist(out) }))
    if (nrow(PL)) { PL[, cpuh := n_rows * sec_per(rung) / 3600, by = seq_len(nrow(PL))]; left <- CAP_CPUH; PL[, funded := FALSE]
      for (i in seq_len(nrow(PL))) if (PL$cpuh[i] <= left) { PL$funded[i] <- TRUE; left <- left - PL$cpuh[i] } }
    else PL <- data.table(rung = character(0), null = character(0), n_rows = integer(0), cpuh = numeric(0), funded = logical(0))
    print(ST); print(PL)
    fwrite(PL, wg("ext_plan.csv")); fwrite(rbind(ROWS, K7B, fill = TRUE), wg("ladder_round1.csv"))
    ts(sprintf("round 2: %d funded extensions, %d array tasks (submit --ladder-perm --round=2, then --ladder-state --round=2)", sum(PL$funded), nrow(tasks_round2())))
  } else {
    ts("[3.2, 5.7] final rung states")
    PL <- fread(rg("ext_plan.csv")); exp2 <- tasks_round2()[, .(n = sum(to - from + 1L)), by = .(rung, null)]
    if (nrow(exp2)) { got2 <- if (is.null(N2)) data.table(rung = character(0), null = character(0), N = integer(0)) else N2[, .N, by = .(rung, null)]
      miss <- merge(exp2, got2, all.x = TRUE)[is.na(N) | N != n]; if (nrow(miss)) stop("round 2 is incomplete: ", paste(miss$rung, miss$null, collapse = "; ")) }
    ST[, unfunded := rung %in% PL[funded == FALSE, rung]]
    # 5.7: the draws at or above the observed deciding count (round 1, both nulls)
    lab <- D25[, .(sample, true = lab_chr(status), lib = library_id, period)]
    E57 <- rbindlist(lapply(LADDER, function(id) { cn <- deciding(id); obs <- OBS[[id]]$counts[[cn]]; s <- if (RUNG[id]$set == "25") S25 else S24
      rbindlist(lapply(c("free", "lib"), function(nl) { z <- N1[rung == id & null == nl & get(cn) >= obs]
        if (!nrow(z)) return(data.table(rung = id, null = nl, n_exceed = 0L))
        L <- do.call(rbind, strsplit(z$lab, "")); colnames(L) <- s; tl <- lab[match(s, sample)]
        data.table(rung = id, null = nl, n_exceed = nrow(z), true_AML_labelled_AML = paste(rowSums(L[, tl$true == "A", drop = FALSE] == "A"), collapse = ";"),
                   P01_labelled_AML = paste(rowSums(L[, tl$lib == "GSE185381__lib__P01", drop = FALSE] == "A"), collapse = ";"),
                   P02_labelled_AML = paste(rowSums(L[, tl$lib == "GSE185381__lib__P02", drop = FALSE] == "A"), collapse = ";"),
                   P15_labelled_AML = paste(rowSums(L[, tl$lib == "GSE185381__lib__P15", drop = FALSE] == "A"), collapse = ";"),
                   period_of_AML_labels = paste(apply(L, 1, function(v) paste(names(table(tl$period[v == "A"])), table(tl$period[v == "A"]), sep = ":", collapse = "/")), collapse = ";")) })) }), fill = TRUE)
    E57[, note := "round-1 draws (the first 1,000 of each null)"]
    put(rbindlist(list(ROWS[, row_type := "p"], ST[, row_type := "rung state"], E57[, row_type := "5.7"], K7B[, row_type := "check"]), fill = TRUE), "ladder")
    print(ST)
  }
}

## -- 8. extension rungs and September-procedure readouts (3.3, 5.5, 5.6; depth arm at R5, R6, E1, E2) ---------------------------------

if (MODE == "extend") {   # self-contained (code review 1): needs scope E only, not the ladder's outputs
  if (!gate("E")) stop("section 6: a scope-E deterministic known answer has not passed; the extension readout is 'not run'")
  recheck_inputs(); OUT <- list(); PL <- list(); MG <- list(); DA <- list(); P370E <- list()
  MF <- fread(file.path(TABWS, "ws_s5_de_T.csv"))[, .(gene, lfc0 = logFC)]; if (SMOKE) MF[, lfc0 := 1]   # smoke: no real S5 logFC
  p370_row <- function(rung, f, n) { xx <- merge(merge(P370, f$tt, by = "gene"), MF, by = "gene")
    data.table(section = "5.1", rung = rung, n_tested = nrow(xx), hits_S5 = xx[q < Q_S5 & abs(logFC) > L_S5, .N], hits_Sep = xx[q < Q_SEP & abs(logFC) >= L_SEP, .N],
               sign_share = mean(sign(xx$logFC) == xx$s), attenuation = median(xx$s * xx$logFC / abs(xx$lfc0)),
               circular = if (rung == "R6") "yes" else "partly (the 24 are shared)", AC007952_logFC = f$tt[gene == "AC007952.4", logFC][1],
               AC007952_q = f$tt[gene == "AC007952.4", q][1], n_samples = n, n_genes = f$n_genes, consensus = f$cons) }
  if (SMOKE) set.seed(SMS[5])
  ts("[3.3, 4.1, 5.5, 5.6, depth arm] R5 and R6 (the 24): own fits; R6's frozen nulls on the stored draws feed the extension sentence")
  x5 <- rung_data("R5"); d5 <- x5$d; dr <- if (SMOKE) data.table(sample = d5$sample, median_umi = 10^runif(nrow(d5), 3.5, 4.2)) else dep_rows(d5)
  d5[, log_depth := log10(dr$median_umi[match(sample, dr$sample)])]
  for (id in c("R5", "R6")) { r <- RUNG[id]
    f <- fit_rung(x5$cnt, x5$d, r$engine, r$univ, "status", keep_obj = (id == "R6"))
    fd <- fit_rung(x5$cnt, d5, r$engine, r$univ, c("log_depth", "status"))$tt
    DA[[id]] <- data.table(rung = id, t(reportable(f$tt, fd)))
    if (id == "R6") { obs <- counts_of(f$tt)
      for (nl in c("free", "lib")) { M <- read_draws(paste0(nl, "24"))
        Z <- do.call(rbind, mcl(seq_len(nrow(M)), function(i) counts_of(tryCatch(refit_frozen(f$obj, status_from(M, i, x5$d)), error = function(e) NULL))))
        OUT[[paste("R6", nl)]] <- p_rows("R6", obs, Z, nl, "frozen (September)") }
      PL$R6 <- planted(f$obj, "R6"); P370E$R6 <- p370_row("R6", f, 24L) } }
  MG$R5R6 <- cbind(rung = "R5 = R6 (same counts)", marker_gate_T(x5$cnt, x5$d))
  for (rung in c("E1", "E2")) {
    ts(sprintf("[3.3] %s", rung))
    d <- ext_design(rung); cnt <- ext_counts(d); terms <- if (rung == "E2") c("dataset", "status") else "status"
    f <- fit_rung(cnt, d, "voom_block", "SEP", terms, keep_obj = TRUE); obs <- counts_of(f$tt)
    dr <- if (SMOKE) data.table(sample = d$sample, median_umi = 10^runif(nrow(d), 3.5, 4.2)) else dep_rows(d); d2 <- copy(d)[, log_depth := log10(dr$median_umi[match(sample, dr$sample)])]
    fd <- fit_rung(cnt, d2, "voom_block", "SEP", c(setdiff(terms, "status"), "log_depth", "status"))$tt
    DA[[rung]] <- data.table(rung = rung, t(reportable(f$tt, fd)))
    md <- data.table(arm = d$status, dataset = d$dataset, library_id = d$library_id)
    for (nl in c("free", "lib")) {
      lab <- sprintf("perm|%s|g03_%s|T_NK", if (nl == "free") "sample_within_dataset" else "library", rung)
      Z <- frozen_null(f$obj, md, nl, lab); OUT[[paste(rung, nl)]] <- p_rows(rung, obs, Z, nl, "frozen (September)") }
    PL[[rung]] <- planted(f$obj, rung); MG[[rung]] <- cbind(rung = rung, marker_gate_T(cnt, d)); P370E[[rung]] <- p370_row(rung, f, nrow(d))
    ts(sprintf("    %s: n %d, genes %d, consensus %.3f | c_Sep %d (c_Sep' %d) | c_S5 %d (c_S5' %d)", rung, nrow(d), f$n_genes, f$cons, obs[["c_Sep"]], obs[["c_Sepr"]], obs[["c_S5"]], obs[["c_S5r"]]))
  }
  put(rbindlist(list(rbindlist(OUT)[, row_type := "p"], rbindlist(DA, fill = TRUE)[, row_type := "depth arm"], rbindlist(PL)[, row_type := "5.5 planted"],
                     rbindlist(MG)[, row_type := "5.6 marker gate"], rbindlist(P370E)[, row_type := "5.1/5.2"]), fill = TRUE), "extend")
}

## -- 9. sensitivities (section 4) --------------------------------------------------------------------------------------------------------
if (MODE == "sens") {
  if (!gate("L") || !gate("E")) stop("section 6: --sens needs scope L (4.1, 4.2, 4.4 R6) and scope E (4.3, 4.4 E1, 4.5 E1)")
  recheck_inputs(); LAD <- rtab("ladder"); S <- list()
  ts("[4.1, 4.2] from the ladder table (frozen scheme at R4b and R6; Control4 coded female at R0 and R1)")
  S$s41 <- LAD[row_type == "p" & scheme == "frozen"][, section := "4.1"]
  S$s42 <- LAD[row_type == "p" & rung %in% c("R0f", "R1f") & count == "c_S5"][, section := "4.2"]
  ts("[4.3, 4.4] E1 without AML073; R6 and E1 without Control1-5 (counts, no p)")
  x6 <- rung_data("R6"); dE <- ext_design("E1"); cE <- ext_counts(dE)
  one <- function(lbl, cnt, d) data.table(section = lbl, t(counts_of(fit_rung(cnt, d, "voom_block", "SEP", "status")$tt)))
  S$s43 <- one("4.3 E1 without AML073", cE, dE[sample != "AML073"])
  S$s44 <- rbind(one("4.4 R6 without Control1-5", x6$cnt, x6$d[!sample %in% T15]), one("4.4 E1 without Control1-5", cE, dE[!sample %in% T15]))
  ts("[4.5] R6 and E1 without the library block (frozen scheme, the same draws)")
  for (rung in c("R6", "E1")) {
    if (rung == "R6") { cnt <- x6$cnt; d <- x6$d } else { cnt <- cE; d <- dE }
    f <- fit_rung(cnt, d, "voom", "SEP", "status", keep_obj = TRUE); obs <- counts_of(f$tt)
    for (nl in c("free", "lib")) {
      Z <- if (rung == "R6") { M <- read_draws(paste0(nl, "24")); do.call(rbind, mcl(seq_len(nrow(M)), function(i) counts_of(tryCatch(refit_frozen(f$obj, status_from(M, i, d)), error = function(e) NULL)))) }
           else frozen_null(f$obj, data.table(arm = d$status, dataset = d$dataset, library_id = d$library_id), nl, sprintf("perm|%s|g03_E1|T_NK", if (nl == "free") "sample_within_dataset" else "library"))
      S[[paste("s45", rung, nl)]] <- p_rows(rung, obs, Z, nl, "frozen, no block")[, section := "4.5"] } }
  put(rbindlist(S, fill = TRUE), "sens")
}

## -- 10. --read: 5.3 the stored September table; 5.4 the prediction recomputed with G04's estimates ------------------------------------
if (MODE == "read") {   # the September T_NK table is opened only after the ladder, extension and sensitivity tables exist (code review 1)
  invisible(gate("")); recheck_inputs(); GP <- need_tables(); R <- list()
  f53 <- file.path(SEPT_TAB, "de", "primary__T_NK.csv")
  m53 <- unname(tools::md5sum(f53)); stopifnot(identical(m53, "73e73469d6a01e4a29dca1d3d2c6210a"))
  ts("[5.3] de/primary__T_NK.csv opened (md5 checked)")
  if (!SMOKE) { T53 <- fread(f53); x <- merge(P370, T53, by = "gene")
    R$r53 <- data.table(section = "5.3", item = c("P370 genes in the table", "of 9,281 tested", "hits (baseline rule)", "reportable", "share with S5's sign", "median sign_S5 x logFC", "median signed z"),
                        value = c(nrow(x), nrow(T53), sum(x$hit_baseline), sum(x$reportable), mean(sign(x$logFC) == x$s), median(x$s * x$logFC),
                                  median(x$s * sign(x$logFC) * qnorm(x$P.Value / 2, lower.tail = FALSE))), note = "descriptive; circular (24 of 53 samples shared)") }
  ts("[5.4] the addendum's prediction, then G04's estimates")
  r <- 1.43
  G4P <- fread(file.path(TABWS, "ws_a1_g04_primary.csv"))[row == "summary"]; G4C <- fread(file.path(TABWS, "ws_a1_g04_children.csv"))[section == "3.4 primary" & item == "phi"]
  rho <- c(point = G4P$rho, L = G4P$rho_L, U = G4P$rho_U); phi <- c(point = G4C$value, L = G4C$lo, U = G4C$hi)
  grid <- rbind(data.table(corner = "point", rho = rho[["point"]], phi = phi[["point"]]),
                CJ(rho = rho[c("L", "U")], phi = phi[c("L", "U")])[, corner := "interval corner"], fill = TRUE)
  grid[, `:=`(f_unshrunk = (14 + (5 * rho + 20 * phi) / r) / 39, f_shrunk = (14 + 5 * rho + 20 * phi) / (39 * r))]
  grid[, `:=`(expected_unshrunk = vapply(f_unshrunk, exp_calls, 1), expected_shrunk = vapply(f_shrunk, exp_calls, 1), observed = 1L)]
  R$r54 <- grid[, section := "5.4"]
  R$ab <- fread(file.path(TABWS, "ws_a1_g04_children.csv"))[section == "3.4 primary" & item %in% c("A_b", "phi"), .(section = "G04(b) handed over", item, value, lo, hi)]
  put(rbindlist(R, fill = TRUE), "read")
  print(grid)
}

## -- 11. --combine: section 7 (first match wins), the fixed sentences, the GATES block -------------------------------------------------
if (MODE == "combine") {   # stops (writes nothing) when a table its gate outcome requires is missing (code review 1)
  invisible(gate("")); GP <- need_tables("read"); C <- checks()
  LAD <- if (GP$L) rtab("ladder"); CO <- if (GP$L) rtab("checks_obs"); P3 <- if (GP$L) rtab("p370"); EXT <- if (GP$E) rtab("extend"); RDT <- rtab("read")
  L_ok <- GP$L && all(CO$state == "PASS") && all(LAD[row_type == "check", state] == "PASS")          # K6 and K7b join scope L
  E_ok <- GP$E
  f3 <- function(p) formatC(p, format = "f", digits = 3)
  pq <- function(rg_, cn, sch = "full", tab = LAD) { z <- tab[row_type == "p" & rung == rg_ & count == cn & scheme == sch]
    sprintf("p_free %s, p_lib %s", f3(z[null == "free", p]), f3(z[null == "lib", p])) }
  obs_of <- function(rg_, cn, sch = "full", tab = LAD) tab[row_type == "p" & rung == rg_ & count == cn & scheme == sch & null == "free", obs]
  AXW <- c(held = "beyond its null", borderline = "borderline", lost = "inside its null", `not computed` = "not computed")   # 7's words for the second axis
  FACTOR <- c(R1 = "removing AML0310", R2 = "S5's covariates (hb, gran, sex, CD8 share)", R3 = "T cells alone, without NK cells", R4a = "S5's edgeR quasi-likelihood engine",
              R4b = "the library block", R5 = "the S5 cell definition: marker concordance, the fine-label system and the in_ccc_graph filter", R6 = "S5's gene universe")
  LIST <- "not set by G03"; AXIS <- "not read"
  if (!L_ok) { CAT <- "not decisive"; SENT <- "A scope-L known answer failed; the ladder is not read."
  } else {
    ST <- LAD[row_type == "rung state"]; st <- setNames(ST$state, ST$rung); unf <- setNames(as.logical(ST$unfunded), ST$rung)
    AXIS <- AXW[[worse(LAD[row_type == "p" & rung == "R6" & count == "c_Sepr" & scheme == "full", state])]]
    nulls_of <- function(k) { z <- ST[rung == k]; paste(c(if (z$free == "lost") "sample null", if (z$lib == "lost") "library null"), collapse = " and ") }
    bl_note <- function(k) if (isTRUE(unf[[k]])) "(1,000 draws; extension not funded under the cap)" else "(Monte Carlo interval straddles 0.05 after 10,000 draws)"
    bad <- LADDER[!st[LADDER] %in% c("held", "borderline", "lost")]
    lost <- LADDER[-1][st[LADDER[-1]] == "lost"]; border <- LADDER[st[LADDER] == "borderline"]; VALID_LOSS <- FALSE; LINES <- character(0)
    if (length(bad)) { CAT <- "not decisive"; SENT <- sprintf("Rung state not computed at %s; the ladder is not read.", paste(bad, collapse = ", "))
    } else if (st[["R0"]] == "lost") {
      CAT <- "S5's count is not beyond its own null"; LIST <- "exploratory"
      z0 <- LAD[row_type == "p" & rung == "R0" & count == "c_S5" & scheme == "full"]
      SENT <- sprintf("S5's own count is not beyond its permutation null (%s; draws: sample null %s, library null %s). The per-gene list is reported as exploratory.",
                      pq("R0", "c_S5"), z0[null == "free", draws], z0[null == "lib", draws])
    } else if (length(lost)) {
      for (k in lost) { prev <- LADDER[match(k, LADDER) - 1L]
        if (st[[prev]] != "held") { LINES <- c(LINES, sprintf("%s: lost (%s), not attributable to %s, because %s is %s.", k, nulls_of(k), FACTOR[[k]], prev, st[[prev]])); next }
        if (k == "R1") { r1 <- P3[section == "5.8" & rung == "R1 ref", range(c_S5, na.rm = TRUE)]
          LINES <- c(LINES, sprintf("R1 (validity): S5's per-gene count does not survive removing AML0310 (lost under the %s; 5.8: removing any other single AML gives c_S5 from %d to %d).", nulls_of(k), r1[1], r1[2])); VALID_LOSS <- TRUE }
        else if (k == "R4b") { LINES <- c(LINES, sprintf("R4b (validity): S5's per-gene count does not survive the library-level correlation of a blocked model (lost under the %s).", nulls_of(k))); VALID_LOSS <- TRUE }
        else LINES <- c(LINES, sprintf("%s (specification): S5's per-gene count needs %s. Without it the count falls inside its permutation null (rung table; lost under the %s).", k, FACTOR[[k]], nulls_of(k))) }
      i1 <- match(lost[1], LADDER); after <- if (i1 < length(LADDER)) LADDER[(i1 + 1L):length(LADDER)] else character(0)
      again <- after[st[after] == "held"]
      if (length(again)) LINES <- c(LINES, sprintf("Held again at %s.", paste(again, collapse = ", ")))
      if (length(border)) LINES <- c(LINES, sprintf("Borderline: %s.", paste(sprintf("%s %s", border, vapply(border, bl_note, "")), collapse = "; ")))
      CAT <- sprintf("lost at %s", paste(lost, collapse = ", "))
      LIST <- if (VALID_LOSS) "exploratory" else "stands with the qualifier(s) above"
      SENT <- paste(c(LINES, if (VALID_LOSS) "The per-gene list is reported as exploratory. The sample-level separation (A1-G02, A1-BATCH G2) is a different statistic and stands as reported there."
                              else "The per-gene list is reported with the qualifier(s) above."), collapse = " ")
    } else if (length(border)) {
      CAT <- "borderline"; LIST <- "stands with the borderline qualifier"
      SENT <- sprintf("S5's per-gene count is borderline at %s; no rung loses it. The per-gene list carries this qualifier.", paste(sprintf("%s %s", border, vapply(border, bl_note, "")), collapse = "; "))
    } else if (all(st[LADDER] == "held")) {
      CAT <- "held through R6"; LIST <- "stands"
      z6 <- worse(LAD[row_type == "p" & rung == "R6" & count == "c_Sepr" & scheme == "full", state])
      SENT <- sprintf(paste("On the 24 adults that S5 and the September fit share, S5's rule keeps %d hits beyond both permutation nulls (%s) through the September cell definition,",
                            "engine with library block, and gene universe. Under the September rule the count is %d (%s; %s). The two analyses count different effects: S5 counts",
                            "|logFC| > 0.5 at FDR 0.10, September |logFC| >= 1 at q 0.05. The September T_NK result stays reported as a null for |logFC| >= 1.5 per gene."),
                      obs_of("R6", "c_S5r"), pq("R6", "c_S5r"), obs_of("R6", "c_Sepr"), pq("R6", "c_Sepr"), AXW[[z6]])
    } else stop("unreachable rung-state pattern")
  }
  EXTS <- if (!E_ok) "Extension readout not run (a scope-E known answer failed)." else {
    e <- function(r, cn) EXT[row_type == "p" & rung == r & count == cn]
    pe <- function(r, cn) { z <- e(r, cn); sprintf("p_free %s, p_lib %s", f3(z[null == "free", p]), f3(z[null == "lib", p])) }
    oc <- function(r) sprintf("%d (%d)", e(r, "c_Sep")[null == "free", obs], e(r, "c_Sepr")[null == "free", obs])
    h <- EXT[row_type == "5.1/5.2"]; hp <- function(r) h[rung == r, hits_Sep]
    sprintf(paste("On the 24 adults the September specification finds %s hits (%s). Adding the 18 children gives %s (%s); adding Chen2023 gives %s (%s).",
                  "The September fit on 53 found 1, inside its nulls. G03 does not test which samples account for any difference.",
                  "P370 genes that are hits under the September rule: R6 %d, E1 %d, E2 %d."),
            oc("R6"), pe("R6", "c_Sep"), oc("E1"), pe("E1", "c_Sep"), oc("E2"), pe("E2", "c_Sep"), hp("R6"), hp("E1"), hp("E2")) }
  K9ok <- all(C[id == "K9", state] == "PASS")
  PRED <- if (!K9ok) "5.4 not reported (K9 failed: the method does not reproduce the addendum)." else { g <- RDT[section == "5.4"]; pt <- g[corner == "point"]; cn <- g[corner != "point"]
    sprintf(paste("With G04's estimates, the addendum's method expects %.1f (unshrunk) and %.1f (shrunk) September calls among the 370, ranging from %.1f to %.1f over the",
                  "interval corners; the fit called 1. These are rough upper estimates from a planted-curve approximation. They are not a test, and they are not read as explaining the September result."),
            pt$expected_unshrunk, pt$expected_shrunk, min(cn[, c(expected_unshrunk, expected_shrunk)]), max(cn[, c(expected_unshrunk, expected_shrunk)])) }
  z <- RDT[section == "G04(b) handed over"]
  AB <- sprintf("A_b %.3f [%.3f, %.3f]; phi %.3f [%.3f, %.3f]", z[item == "A_b", value], z[item == "A_b", lo], z[item == "A_b", hi], z[item == "phi", value], z[item == "phi", lo], z[item == "phi", hi])
  RT <- if (!L_ok) "" else { s <- LAD[row_type == "rung state"][match(LADDER, rung)]; zz0 <- LAD[row_type == "p" & scheme == "full"]
    paste(vapply(LADDER, function(k) { cn <- deciding(k); zz <- zz0[rung == k & count == cn]
      sprintf("| %s | %s | %d | %s (%s) | %s (%s) | %d | %s |", k, RUNG[k]$class, zz[null == "free", obs], f3(zz[null == "free", p]), zz[null == "free", state],
              f3(zz[null == "lib", p]), zz[null == "lib", state], zz0[rung == k & count == "c_Sep" & null == "free", obs], s[rung == k, state]) }, ""), collapse = "\n") }
  SUM <- data.table(item = c("ladder category", "second axis (R6, September rule, c_Sep')", "per-gene list status", "ladder sentence", "extension sentence", "prediction sentence",
                             "G04(b) handed over", "scope-L known answers (incl. K6, K7b)", "scope-E known answers"),
                    value = c(CAT, AXIS, LIST, SENT, EXTS, PRED, AB, if (L_ok) "all PASS" else "FAIL", if (E_ok) "all PASS" else "FAIL"))
  put(SUM, "summary")
  BLOCK <- c("<!-- A1-G03 start -->", sprintf("## S5's T imprint vs the September T_NK null, one factor at a time (A1-G03, prereg-a1-g03-v1) -- %s", format(Sys.Date())), "",
             "Script: scripts/ws/42_a1_g03.R. Ladder R0-R6 on S5's samples (R1-R6: the 24 shared with the September fit); E1 +18 children, E2 +9 Chen2023 (G04 5.7).", "",
             if (L_ok) c("| rung | class | deciding count | p_free (state) | p_lib (state) | c_Sep | rung state |", "|---|---|---|---|---|---|---|", RT, ""),
             sprintf("**Verdict: A1-G03 ladder = %s; R6 under the September rule: %s; per-gene list: %s.** G03 has no pass/fail gate of its own.", CAT, AXIS, LIST), "",
             paste("Ladder:", SENT), "", paste("Extension:", EXTS), "", paste("Prediction:", PRED), "", paste("G04(b), handed over:", AB), "<!-- A1-G03 end -->")
  if (SMOKE) writeLines(BLOCK, rg("gates_block.md")) else cat(paste0("\n", paste(BLOCK, collapse = "\n"), "\n"), file = here::here("scripts", "ws", "GATES.md"), append = TRUE)
  ts("[combine] done"); writeLines(BLOCK)
}
