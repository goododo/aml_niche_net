#!/usr/bin/env Rscript
# INPUT  : PREREGISTRATION_A1_G09.md (tag prereg-a1-g09-v1). /LARGE1/.../ws_pseudobulk/<ds>/ws_pb_<ds>.rds for GSE120221 (Oetjen),
#          GSE185381, GSE227690 (Griffin), GSE276360, GSE116256, GSE289435, GSE227903; per-sample objects in 02_seurat_objects/04_annotated/<ds>/;
#          results/tables/ws/ ws_s4_pseudobulk_meta_<ds>.csv, ws_s5_imprint_genes.csv, ws_s8_t_imprint_ambient_prone.csv, ws_s8_tests.csv,
#          ws_s8_scores_{rep,rep2}.csv, ws_s5_scores.csv, ws_a1_g02_fold_genes.csv, ws_a1_g01_axis.csv, ws_s9_scores_GSE120221.csv (count columns only)
# OUTPUT : results/tables/ws/ws_a1_g09_{checks,tests,balance,draws8,summary}.csv; per-donor work tables under /LARGE1/.../ws_a1_g09/;
#          GATES.md block A1-G09 (and, only under 7.2(iii), a dated erratum note after block A1-S8). Smoke runs write only to /LARGE1/.../ws_a1_g09/smoke/.
# what it does: PREREGISTRATION_A1_G09.md sections 2-7. S8's T2 (up- and down-half means move oppositely across samples) and T1 (axis genes
#   co-vary within a sample's T cells) on healthy donors (a) and, in the AML cohorts, against per-half and loading-matched random sets (b).
#   Nulls: N0 = S8's own (head/tail split of a bin-ordered list), N1 = per-half on S8's strata, N1-L = per-half on N2's universe, N2 = per-half x L.
#   --known  K0-K7 (section 6)      --run  3.1-3.3, section 4, K8 (balance; before any statistic)      --combine  section 7, GATES block
# Functions copied from 32_s8_replication.R (S8, unchanged there) and 39_a1_g01.R (G01 5.5 loading), with line references.
# This script deletes nothing and modifies none of its inputs.

options(warn = 1)
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(Seurat); library(SeuratObject); library(here); library(parallel) })
.a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% .a
MODES <- c("known", "run", "combine"); MODE <- MODES[paste0("--", MODES) %in% .a]
if (length(MODE) != 1L) stop("give exactly one of ", paste0("--", MODES, collapse = " "), " [--smoke]")
bad_arg <- setdiff(.a, c(paste0("--", MODES), "--smoke"))   # a mistyped --smoke would otherwise run on real data
if (length(bad_arg)) stop("unknown argument(s): ", paste(bad_arg, collapse = " "))
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R"))
NC <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))); if (is.na(NC) || NC < 1L) NC <- 2L
NC <- min(NC, 16L)   # section 8: 16 workers (the node can report 32 CPUs for c=16, which halves the memory per worker)
TABWS <- file.path(TAB_DIR, "ws"); METADIR <- file.path(RAW_DIR, "metadata"); GATES <- here::here("scripts", "ws", "GATES.md")
GDIR <- file.path(LARGE1_DIR, "ws_a1_g09"); SDIR <- file.path(GDIR, "smoke")
WDIR <- if (SMOKE) SDIR else GDIR; ODIR <- if (SMOKE) SDIR else TABWS; dir.create(WDIR, recursive = TRUE, showWarnings = FALSE)
SFX <- if (SMOKE) "_smoke" else ""
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
MD5 <- if (length(SELF) && !is.na(SELF) && file.exists(SELF)) unname(tools::md5sum(SELF)) else NA_character_
if (!SMOKE && is.na(MD5)) stop("run the script with Rscript (its md5 is needed for every output)")
if (!SMOKE && file.access(SELF, 2) == 0) stop("a real run must use the dated read-only copy: sbatch scripts/ws/submit_ws_a1_g09.sbatch --", MODE)
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
wtab <- function(x) { f <- file.path(ODIR, paste0("ws_a1_g09_", x, SFX, ".csv")); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
put  <- function(d, x) { d <- as.data.table(d); if (!ncol(d)) d <- data.table(script_md5 = character(0)) else d[, script_md5 := MD5]; fwrite(d, wtab(x)); invisible(d) }
rtab <- function(x) { f <- file.path(ODIR, paste0("ws_a1_g09_", x, SFX, ".csv")); if (!file.exists(f)) stop("missing ", f, " (run its mode first)"); fread(f) }
wg   <- function(x) { f <- file.path(WDIR, sub("(\\.[a-z.]+)$", paste0(SFX, "\\1"), x)); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
if (!SMOKE) { tg <- switch(MODE, known = file.path(TABWS, "ws_a1_g09_checks.csv"),
                           run = c(file.path(TABWS, paste0("ws_a1_g09_", c("tests", "balance", "draws8"), ".csv")), file.path(GDIR, c("donor_values.csv", "naive_share_spread.csv", "k0b.csv"))),
                           combine = file.path(TABWS, "ws_a1_g09_summary.csv"))
  if (any(file.exists(tg))) stop("refusing to overwrite ", paste(tg[file.exists(tg)], collapse = ", "), " (a rerun needs a dated note and the old file renamed)")
  if (MODE == "combine" && any(grepl("A1-G09 start", readLines(GATES)))) stop("GATES.md already holds an A1-G09 block") }
ts(sprintf("[%s%s] %d worker processes", MODE, if (SMOKE) " smoke" else "", NC))

## -- 0. constants (section numbers refer to the pre-registration) ---------------------------------------------------------------
N_RAND <- if (SMOKE) 100L else 1000L; N_RAND_CELL <- if (SMOKE) 20L else 200L; CAP <- 500L; MINC <- 20L; DET_MIN <- 0.10; MIN_HALF <- 20L   # S8 (32:23)
N_D8 <- if (SMOKE) 6L else 500L; N_D8_NULL <- if (SMOKE) 40L else 200L                                                                    # 3.1
N_K3 <- if (SMOKE) 10L else 100L; N_K4 <- if (SMOKE) 5L else 20L; N_K5 <- if (SMOKE) 5L else 20L; N_K7 <- if (SMOKE) 5L else 20L       # section 6
OFF <- c(cap = 901L, d8 = 902L, n0 = 903L, n1 = 904L, n2 = 905L, worlds = 906L, k7 = 907L, n1l = 908L, smoke = 909L)                    # 8.3 (n1l, smoke: dated note)
T15 <- paste0("Control", 1:5)
NAIVE <- c("CD4 Naive", "CD8 Naive"); EM <- c("CD8 Effector Memory 1", "CD8 Effector Memory 2", "CD4 Effector Memory")   # 39:119
MERGE_RATIO <- 3; BAL_MAX <- 0.1; K_L_T2 <- 20L; K_L_T1 <- 10L; K_D_T1 <- 10L   # L ventiles (T2), L and detection deciles (T1 N1-L/N2); dated pre-run note (smoke balance)
hash_int <- function(x) { v <- utf8ToInt(x); as.integer(sum(as.numeric(v) * seq_along(v)) %% 99991) }   # deterministic, fork-safe label hash
bseed <- function(key, label) SEED + OFF[[key]] * 100000L + hash_int(label)   # one seed per (stream, cohort/gene set/null/donor) block (8.3)

## known-answer bookkeeping (G07's three states)
CHK <- list()
chk <- function(id, name, expected, got, tol = NULL) { ok <- if (is.null(tol)) isTRUE(all.equal(expected, got)) else length(got) == length(expected) && isTRUE(all(abs(expected - got) <= tol))
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, kind = "deterministic", check = name, expected = paste(format(expected), collapse = ";"), got = paste(format(got), collapse = ";"), state = if (ok) "PASS" else "FAIL_deterministic")
  message(sprintf("  %-4s %-96s %s", id, name, if (ok) "PASS" else paste("** FAIL ** got", paste(format(got), collapse = ";")))) }
chk_stat <- function(id, name, bar, got, ok) { CHK[[length(CHK) + 1L]] <<- data.table(id = id, kind = "statistical", check = name, expected = bar, got = paste(got, collapse = ";"), state = if (isTRUE(ok)) "PASS" else "FAIL_statistical")
  message(sprintf("  %-4s %-96s bar %-12s got %-10s %s", id, name, bar, paste(got, collapse = ";"), if (isTRUE(ok)) "PASS" else "** FAIL (statistical) **")) }
rep_ <- function(id, name, value) { CHK[[length(CHK) + 1L]] <<- data.table(id = id, kind = "reported", check = name, expected = "", got = paste(value, collapse = ";"), state = "REPORTED")
  message(sprintf("  %-4s %-96s %s [REPORTED]", id, name, paste(value, collapse = ";"))) }
na_chk <- function(id, name, v) chk(id, sprintf("%s: every world was computed (no stratum FAIL or gene shortfall)", name), 0L, sum(is.na(v)))   # "not computed" is never "did not pass"
mcl <- function(X, FUN) { r <- mclapply(X, FUN, mc.cores = NC, mc.preschedule = FALSE)
  bad <- vapply(r, function(z) inherits(z, "try-error") || is.null(z), TRUE); if (any(bad)) stop(sum(bad), " of ", length(r), " workers failed: ", paste(unique(unlist(lapply(r[bad], as.character))), collapse = " | ")); r }
syn <- function(M) { attr(M, "synthetic") <- TRUE; M }
cols <- function(M, j) { X <- M[, j, drop = FALSE]; if (isTRUE(attr(M, "synthetic"))) syn(X) else X }   # column subset that keeps the smoke mark
need_syn <- function(M, what) if (SMOKE && !isTRUE(attr(M, "synthetic"))) stop("smoke run: a real ", what, " matrix would enter a computation")

## -- 1. gene sets and the S8 blacklist (32:36-44) --------------------------------------------------------------------------------
IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))
GS <- list(primary = IMP[!gene %in% AMB$gene], full = IMP)
stopifnot(GS$primary[, sprintf("%d/%d/%d", .N, sum(direction == "up"), sum(direction == "down"))] == "370/156/214",
          GS$full[, sprintf("%d/%d/%d", .N, sum(direction == "up"), sum(direction == "down"))] == "402/171/231")
UPDN <- lapply(GS, function(g) list(up = g[direction == "up", gene], dn = g[direction == "down", gene]))
BL <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_s1_ambient_blacklist.csv"))
GO <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
black_for <- function(G) unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene], "XIST", GO[chr == "chrY", gene],
                                  grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))

## -- 2. S8's helpers, verbatim (32:47-61) -----------------------------------------------------------------------------------------
logz <- function(M) { need_syn(M, "pseudobulk")
  M <- M[!rownames(M) %in% black_for(rownames(M)), , drop = FALSE]; y <- DGEList(M)
  keep <- rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  L <- cpm(y, log = TRUE, prior.count = 1); Z <- t(scale(t(L)))
  U <- data.table(gene = rownames(L), mean_lcpm = rowMeans(L), det = rowMeans(M[rownames(L), , drop = FALSE] > 0))
  U[, bin := paste(cut(mean_lcpm, quantile(mean_lcpm, 0:10 / 10), include.lowest = TRUE), cut(det, c(-Inf, 0.5, 0.8, 0.95, Inf)))]
  list(Z = Z, U = U) }
matched <- function(U, genes) { tb <- table(U[gene %in% genes, bin])
  unlist(lapply(names(tb), function(b) { pool <- U[bin == b & !gene %in% genes, gene]; if (!length(pool)) pool <- U[!gene %in% genes, gene]; sample(pool, tb[[b]], replace = length(pool) < tb[[b]]) })) }
t2_stat <- function(Z, up, dn) { up <- intersect(up, rownames(Z)); dn <- intersect(dn, rownames(Z)); if (length(up) < 5 || length(dn) < 5 || ncol(Z) < 4) return(NA_real_)
  cor(colMeans(Z[up, , drop = FALSE]), colMeans(Z[dn, , drop = FALSE]), method = "spearman") }
split_random <- function(U, up, dn) { g <- matched(U, c(up, dn)); list(up = head(g, length(intersect(up, U$gene))), dn = tail(g, length(intersect(dn, U$gene)))) }
coherence <- function(C, up, dn) { if (length(up) < 2 || length(dn) < 1) return(NA_real_); A <- C[up, up]; mean(A[upper.tri(A)]) - mean(C[up, dn]) }
norm_cells <- function(Xc) { need_syn(Xc, "per-cell")   # 32:81, per-cell test: genes detected in >= 10 % of the capped cells, log-normalised to 1e4
  det <- Matrix::rowMeans(Xc > 0); keepg <- names(det)[det >= DET_MIN]; Xn <- log1p(t(t(as.matrix(Xc[keepg, , drop = FALSE])) / colSums(Xc)) * 1e4)
  list(X = Xn, det = det[keepg]) }

## -- 3. G01's loading (39:75-76, 91-94, 536-540), verbatim ------------------------------------------------------------------------
P185 <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds")); G185 <- P185$genes
BLACK_G01 <- unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene],
                      "XIST", GO[chr == "chrY", gene], grep("^MT-", G185, value = TRUE), grep("^RP[LS][0-9]", G185, value = TRUE)))
logz_g01 <- function(M) { M <- M[!rownames(M) %in% BLACK_G01, , drop = FALSE]; y <- DGEList(M)
  keep <- rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE]); cpm(y, log = TRUE, prior.count = 1) }
L_of <- function(fine, cols = NULL) {   # cols: library columns to keep (leave-one-donor-out, or a donor subset); NULL = every column (G01)
  pool <- function(ks) Reduce(`+`, lapply(ks, function(k) { M <- fine[[k]]; if (!is.null(cols)) M <- M[, intersect(colnames(M), cols), drop = FALSE]; rowSums(M) }))
  Mo <- cbind(naive = pool(NAIVE), em = pool(EM)); if (all(vapply(fine[c(NAIVE, EM)], function(M) isTRUE(attr(M, "synthetic")), TRUE))) Mo <- syn(Mo)
  need_syn(Mo, "loading"); Lo <- logz_g01(Mo)
  cp <- cpm(calcNormFactors(DGEList(Mo[rownames(Lo), ])))
  log2((cp[, "em"] + 1) / (cp[, "naive"] + 1)) }

## -- 4. per-half nulls: strata, merges, draws (section 3) -------------------------------------------------------------------------
qbin <- function(x, k) cut(rank(x, ties.method = "first"), k, labels = FALSE)            # rank quantiles (G01 39:544 idiom), integer codes
make_strata <- function(U, up, dn, dims, Lv = NULL, k_l = 10L) {
  # U: data.table(gene, e, d) with integer codes (missing dims = 1); returns groups after the 3 x need merge rule, or an error label
  U <- copy(U); if (!"e" %in% names(U)) U[, e := 1L]; if (!"d" %in% names(U)) U[, d := 1L]
  U[, l := if ("l" %in% dims) qbin(Lv[gene], k_l) else 1L]
  imp <- c(up, dn); U[, `:=`(isu = gene %in% up, isd = gene %in% dn)]
  cells <- unique(U[, .(e, d, l)]); cells[, grp := seq_len(.N)]
  U <- merge(U, cells, by = c("e", "d", "l")); med_l <- stats::median(unique(U$l)); n_merge <- 0L; done <- integer(0)
  repeat {
    gs <- U[, .(need = sum(isu | isd), pool = sum(!(isu | isd))), by = grp]; bad <- gs[need > 0 & pool < MERGE_RATIO * need & !grp %in% done]
    if (!nrow(bad)) break
    g <- bad$grp[1]; mem <- unique(U[grp == g, .(e, d, l)]); cand <- integer(0)
    for (dm in dims) {   # merge order: detection, then expression (T2), then L toward its median
      nb <- if (dm == "l") { sg <- sign(med_l - mem$l); m1 <- copy(mem)[sg != 0][, l := l + sg[sg != 0]]; m0 <- mem[sg == 0]
        rbind(m1, copy(m0)[, l := l - 1L], copy(m0)[, l := l + 1L]) } else
        unique(rbindlist(lapply(c(-1L, 1L), function(st) { m <- copy(mem); m[, (dm) := get(dm) + st]; m })))
      cand <- unique(U[nb, on = c("e", "d", "l"), nomatch = 0][grp != g, grp]); if (length(cand)) break }
    if (!length(cand)) { if (bad$pool[1] >= bad$need[1]) { done <- c(done, g); next }   # no neighbour left: section 3 fails only when pool < need
      return(list(fail = sprintf("stratum %d: pool %d < need %d after every merge", g, bad$pool[1], bad$need[1]))) }
    tgt <- gs[grp %in% cand][order(pool)][1, grp]; U[grp == tgt, grp := g]; n_merge <- n_merge + 1L }
  pools <- split(U[!(isu | isd), gene], U[!(isu | isd), grp])
  list(nu = U[isu == TRUE, .N, by = grp], nd = U[isd == TRUE, .N, by = grp], pools = pools, n_merge = n_merge, grp_of = setNames(U$grp, U$gene)) }
draw_halves <- function(S) {   # matched_x: strata from each target half, pool excludes every imprint gene (and the drawn up half for the down half)
  gu <- unlist(lapply(seq_len(nrow(S$nu)), function(i) { p <- S$pools[[as.character(S$nu$grp[i])]]; p[sample.int(length(p), S$nu$N[i])] }), use.names = FALSE)
  gd <- unlist(lapply(seq_len(nrow(S$nd)), function(i) { p <- setdiff(S$pools[[as.character(S$nd$grp[i])]], gu); p[sample.int(length(p), S$nd$N[i])] }), use.names = FALSE)
  list(up = gu, dn = gd) }
draws_of <- function(S, n, seed) { set.seed(seed); lapply(seq_len(n), function(i) draw_halves(S)) }
balance <- function(S, D, up, dn, Lv = NULL, det = NULL) {   # K8 on merged strata
  cnt_ok <- all(vapply(D, function(r) { a <- table(S$grp_of[r$up]); b <- setNames(S$nu$N, S$nu$grp); c2 <- table(S$grp_of[r$dn]); d2 <- setNames(S$nd$N, S$nd$grp)
    identical(as.integer(a[as.character(names(b))]), as.integer(b)) && identical(as.integer(c2[as.character(names(d2))]), as.integer(d2)) && !anyNA(c(r$up, r$dn)) }, TRUE))
  jac <- if (length(D) > 1) mean(vapply(2:min(length(D), 50), function(i) length(intersect(D[[i]]$up, D[[i - 1]]$up)) / length(union(D[[i]]$up, D[[i - 1]]$up)), 1)) else NA_real_
  bl <- function(v) if (is.null(v)) NA_real_ else { s <- sd(v[names(S$grp_of)], na.rm = TRUE); max(mean(vapply(D, function(r) abs(mean(v[r$up]) - mean(v[up])), 1)), mean(vapply(D, function(r) abs(mean(v[r$dn]) - mean(v[dn])), 1))) / s }   # SD over the matched universe the strata were built on
  list(counts_ok = cnt_ok, jaccard = jac, bal_L = bl(Lv), bal_det = bl(det), n_merge = S$n_merge) }

## -- 5. cohorts: donors, pseudobulks, per-cell matrices (S9 33:110-126 pooling for Oetjen; Griffin donors from the sample name) -------
COH <- list(oetjen = list(ds = "GSE120221", roles = "control"), g185h = list(ds = "GSE185381", roles = "control"), griffin = list(ds = "GSE227690", roles = "control"),
            gse276360 = list(ds = "GSE276360", roles = "control"), seqwell_h = list(ds = "GSE116256", roles = "control"),
            rep = list(ds = "GSE289435", roles = "AML_Dx"), rep2 = list(ds = "GSE227903", roles = "AML_Dx"))
donor_map <- function(ds, libs) {
  if (ds == "GSE120221") { don <- fread(file.path(METADIR, ds, "donor_ages_GSE120221_v1.0.csv"), colClasses = "character")
    d <- don$donor[match(sub("_.*$", "", libs), don$gsm)]; a <- as.numeric(don$age[match(sub("_.*$", "", libs), don$gsm)]); return(data.table(sample = libs, donor = d, age = a)) }
  if (ds == "GSE227690") return(data.table(sample = libs, donor = sub("^GSM[0-9]+_(BM[0-9]+)_.*$", "\\1", libs), age = NA_real_))
  data.table(sample = libs, donor = libs, age = NA_real_) }
load_cohort <- function(key) { ds <- COH[[key]]$ds
  P <- if (ds == "GSE185381") P185 else readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", ds, paste0("ws_pb_", ds, ".rds")))
  META <- fread(file.path(TABWS, paste0("ws_s4_pseudobulk_meta_", ds, ".csv")))
  libs <- unique(META[role %in% COH[[key]]$roles & paediatric == FALSE, sample])
  if (key == "rep2") stopifnot(length(libs) == 9L)                                                     # K6 (asserted again in --known)
  mp <- donor_map(ds, libs); stopifnot(!anyNA(mp$donor))
  kc <- P$cells[kept == TRUE & group4 == "T" & sample %in% libs]
  nT <- kc[, .N, by = sample]; mp[, n_T := nT$N[match(sample, nT$sample)]][is.na(n_T), n_T := 0L]
  DON <- mp[, .(n_lib = .N, libraries = paste(sort(sample), collapse = ";"), n_T = sum(n_T), age = age[order(sample)][1]), by = donor][n_T >= MINC][order(donor)]
  ns <- kc[, .(naive = sum(fine %in% NAIVE), em = sum(fine %in% EM), n = .N), by = sample]; ns <- merge(mp, ns, by = "sample")[, .(naive_share = sum(naive) / sum(n), em_share = sum(em) / sum(n)), by = donor]
  DON <- merge(DON, ns, by = "donor", all.x = TRUE)
  PBT <- if (SMOKE) NULL else sapply(DON$donor, function(d) rowSums(P$raw$T[, mp[donor == d, sample], drop = FALSE]))
  if (!is.null(PBT)) colnames(PBT) <- DON$donor
  list(key = key, ds = ds, genes = P$genes, mp = mp, DON = DON, PBT = PBT, cells = kc[, .(cell, sample, fine)],
       fine = if (SMOKE || key %in% c("rep", "rep2")) NULL else P$fine) }   # K6: no AML per-label matrix (and no MRD/relapse column) is kept
read_cells <- function(co, d, seed, lib_only = NULL) {   # raw counts of a donor's kept T cells (pooled over libraries), capped at 500 (one subsample per donor)
  libs <- if (is.null(lib_only)) co$mp[donor == d, sample] else lib_only
  AN <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated", co$ds)
  X <- do.call(cbind, lapply(libs, function(sm) { o <- readRDS(file.path(AN, paste0(sm, ".rds"))); M <- LayerData(o, assay = "RNA", layer = "counts"); rm(o)
    cc <- co$cells[sample == sm]; M <- M[, cc$cell, drop = FALSE]; stopifnot(identical(rownames(M), co$genes)); M }))
  lab <- co$cells[match(colnames(X), cell), fine]
  j <- if (ncol(X) > CAP) { set.seed(seed); sample.int(ncol(X), CAP) } else seq_len(ncol(X))
  list(X = X[, j, drop = FALSE], lab = lab[j], n_all = ncol(X)) }

## -- 6. smoke worlds: synthetic counts on the real IDs and gene axes (no real count enters a smoke computation) ------------------------
if (SMOKE) {
  set.seed(SEED + OFF[["smoke"]])
  GALL <- sort(unique(unlist(lapply(COH, function(x) readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", x$ds, paste0("ws_pb_", x$ds, ".rds")))$genes))))   # gene names only
  MU0 <- setNames(exp(rnorm(length(GALL), -1, 2)), GALL); MU0[sample.int(length(GALL), length(GALL) %/% 3)] <- 0   # one biology shared by every cohort
  EMW0 <- setNames(exp(rnorm(length(GALL), 0, 0.7)), GALL)
  u0 <- intersect(UPDN$primary$up, GALL); d0 <- intersect(UPDN$primary$dn, GALL); MU0[c(u0, d0)] <- exp(rnorm(length(c(u0, d0)), 1.5, 0.7)); EMW0[u0] <- EMW0[u0] * 2; EMW0[d0] <- EMW0[d0] * 0.5
  smoke_world <- function(co, t2_struct, t1_struct, comp) {   # t2_struct: donor factor on the axis halves; t1_struct: cell factor; comp: composition-driven world
    G <- co$genes; ng <- length(G); mu <- MU0[G]; emw <- EMW0[G]                                                   # EM vs naive differs genome-wide; up genes lean EM, down genes naive
    up <- intersect(UPDN$primary$up, G); dn <- intersect(UPDN$primary$dn, G)
    nd <- nrow(co$DON); f <- rnorm(nd); s <- runif(nd, 0.2, 0.8)
    PB <- sapply(seq_len(nd), function(i) { m <- mu * 2e5 / sum(mu)
      if (t2_struct) { m[up] <- m[up] * 2^(0.5 * f[i]); m[dn] <- m[dn] * 2^(-0.5 * f[i]) }
      if (comp) m <- m * (s[i] * 1 + (1 - s[i]) * emw) / (0.5 + 0.5 * emw)
      rnbinom(ng, size = 20, mu = m) }); dimnames(PB) <- list(G, co$DON$donor)
    co$PBT <- syn(PB); co$DON[, naive_share := s]
    co$smoke_cells <- function(d, seed) { set.seed(seed); n <- min(150L, co$DON[donor == d, n_T]); g <- rnorm(n); isem <- runif(n) < 0.5
      M <- sapply(seq_len(n), function(c) { m <- mu * 1500 / sum(mu)
        if (t1_struct) { m[up] <- m[up] * 2^(0.8 * g[c]); m[dn] <- m[dn] * 2^(-0.8 * g[c]) }
        if (comp) m <- m * (if (isem[c]) emw else 1)
        rpois(ng, m) }); X <- Matrix(M, sparse = TRUE); dimnames(X) <- list(G, sprintf("%s_c%03d", d, seq_len(n)))
      list(X = syn(X), lab = ifelse(isem, "CD8 Effector Memory 2", "CD4 Naive"), n_all = n) }
    co$fine <- lapply(setNames(c(NAIVE, EM), c(NAIVE, EM)), function(k) { w <- if (k %in% EM) emw else rep(1, ng)
      M <- sapply(co$mp$sample, function(sm) rnbinom(ng, size = 20, mu = mu * w * 5e4 / sum(mu * w))); dimnames(M) <- list(G, co$mp$sample); syn(M) })   # one column per library, as in the real objects
    co }
}
get_cells <- function(co, d, seed) if (SMOKE) co$smoke_cells(d, seed) else read_cells(co, d, seed)
SPEC <- list(oetjen = c(FALSE, FALSE, TRUE), g185h = c(FALSE, FALSE, FALSE), griffin = c(FALSE, FALSE, FALSE), gse276360 = c(FALSE, FALSE, FALSE),
             seqwell_h = c(FALSE, FALSE, FALSE), rep = c(TRUE, TRUE, FALSE), rep2 = c(FALSE, FALSE, FALSE))   # smoke: (T2 structure, T1 structure, composition)
CO <- lapply(setNames(names(COH), names(COH)), function(k) { co <- load_cohort(k); if (SMOKE) co <- smoke_world(co, SPEC[[k]][1], SPEC[[k]][2], SPEC[[k]][3]); co })
rm(P185); invisible(gc(FALSE))
LD <- list(aml = L_of(CO$oetjen$fine))                                                                         # all-Oetjen L (AML runs)
noT15 <- setdiff(CO$g185h$DON$donor, T15)
LD$t2_oetjen <- L_of(CO$g185h$fine, cols = CO$g185h$mp[donor %in% noT15, sample])                              # 9.7 / 9.10
LD$t2_oetjen_all10 <- L_of(CO$g185h$fine, cols = CO$g185h$mp$sample)                                           # sensitivity 4.2
lodo_L <- function(d) L_of(CO$oetjen$fine, cols = setdiff(unlist(lapply(CO$oetjen$fine, colnames)), CO$oetjen$mp[donor == d, sample]))

## -- 7. the two tests under the four nulls -------------------------------------------------------------------------------------
t2_run <- function(PB, up, dn, Lv, label, nulls = c("N0", "N1", "N1L", "N2"), nr = N_RAND, keep_draws = FALSE) {
  TZ <- logz(PB); Z <- TZ$Z; U <- copy(TZ$U); out <- list(); cov <- mean(c(up, dn) %in% U$gene)
  U[, `:=`(e = cut(mean_lcpm, quantile(mean_lcpm, 0:10 / 10), include.lowest = TRUE, labels = FALSE), d = cut(det, c(-Inf, 0.5, 0.8, 0.95, Inf), labels = FALSE))]
  u1 <- intersect(up, U$gene); d1 <- intersect(dn, U$gene); obs <- t2_stat(Z, u1, d1)
  U2 <- U[gene %in% names(Lv)]; u2 <- intersect(up, U2$gene); d2 <- intersect(dn, U2$gene); obs2 <- t2_stat(Z, u2, d2)
  one <- function(nl, S, uu, dd, ob, Lb = NULL) { if (!is.null(S$fail)) return(data.table(null = nl, statistic = ob, P5 = NA_real_, pass = NA, fail = S$fail))
    D <- draws_of(S, nr, bseed(if (nl == "N1L") "n1l" else tolower(nl), paste(label, nl))); v <- vapply(D, function(r) t2_stat(Z, r$up, r$dn), 1)
    b <- balance(S, D, uu, dd, Lv = Lb); if (keep_draws) out$draws[[nl]] <<- v
    data.table(null = nl, statistic = ob, P5 = unname(quantile(v, 0.05)), P95 = unname(quantile(v, 0.95)), pass = isTRUE(ob < 0 && ob < quantile(v, 0.05)),
               na_draws = sum(is.na(v)), counts_ok = b$counts_ok, bal_L = b$bal_L, jaccard = b$jaccard, n_merge = b$n_merge, fail = "") }
  R <- list()
  if ("N0" %in% nulls) { set.seed(bseed("n0", paste(label, "N0"))); v0 <- replicate(nr, { r <- split_random(U, u1, d1); t2_stat(Z, r$up, r$dn) })
    R$N0 <- data.table(null = "N0", statistic = obs, P5 = unname(quantile(v0, 0.05, na.rm = TRUE)), P95 = unname(quantile(v0, 0.95, na.rm = TRUE)), pass = isTRUE(obs < 0 && obs < quantile(v0, 0.05, na.rm = TRUE)), na_draws = sum(is.na(v0)), fail = "") }
  if ("N1" %in% nulls) R$N1 <- one("N1", make_strata(U, u1, d1, c("d", "e")), u1, d1, obs)
  if ("N1L" %in% nulls) R$N1L <- one("N1L", make_strata(U2, u2, d2, c("d", "e")), u2, d2, obs2)
  if ("N2" %in% nulls) R$N2 <- one("N2", make_strata(U2, u2, d2, c("d", "e", "l"), Lv = Lv, k_l = K_L_T2), u2, d2, obs2, Lb = Lv)
  res <- rbindlist(R, fill = TRUE)[, `:=`(test = "T2", n = ncol(Z), coverage = cov, coverage_ok = cov >= 0.8, n_up = fifelse(null %in% c("N0", "N1"), length(u1), length(u2)), n_dn = fifelse(null %in% c("N0", "N1"), length(d1), length(d2)))]
  list(res = res, Z = Z, draws = out$draws) }
t1_donor <- function(cl, up, dn, Lv, label, nulls = c("N0", "N1", "N1L", "N2"), nr = N_RAND_CELL) {
  u <- intersect(up, rownames(cl$X)); d <- intersect(dn, rownames(cl$X))
  if (length(u) < MIN_HALF || length(d) < MIN_HALF) return(data.table(null = nulls, obs = NA_real_, p95 = NA_real_, above = NA, n_up = length(u), n_dn = length(d), skipped = TRUE))
  C <- cor(t(cl$X[c(u, d), ])); obs <- coherence(C, u, d); R <- list()
  if ("N0" %in% nulls) { Ud <- data.table(gene = names(cl$det), det = cl$det); Ud[, bin := cut(det, unique(quantile(det, 0:10 / 10)), include.lowest = TRUE)]
    set.seed(bseed("n0", paste(label, "N0"))); v <- replicate(nr, { g <- matched(Ud, c(u, d)); gu <- head(g, length(u)); gd <- tail(g, length(d)); coherence(cor(t(cl$X[c(gu, gd), ])), gu, gd) })
    R$N0 <- data.table(null = "N0", obs = obs, p95 = unname(quantile(v, 0.95, na.rm = TRUE)), na_draws = sum(is.na(v))) }
  U1 <- data.table(gene = names(cl$det), d = cut(cl$det, unique(quantile(cl$det, 0:10 / 10)), include.lowest = TRUE, labels = FALSE))
  g2 <- setdiff(intersect(names(cl$det), names(Lv)), black_for(names(cl$det))); u2 <- intersect(u, g2); d2 <- intersect(d, g2)
  U2 <- data.table(gene = g2, d = qbin(cl$det[g2], K_D_T1)); obs2 <- if (length(u2) >= MIN_HALF && length(d2) >= MIN_HALF) coherence(C[c(u2, d2), c(u2, d2)], u2, d2) else NA_real_
  one <- function(nl, S, uu, dd, ob, Lb = NULL) { if (!is.null(S$fail)) return(data.table(null = nl, obs = ob, p95 = NA_real_, fail = S$fail))
    D <- draws_of(S, nr, bseed(if (nl == "N1L") "n1l" else tolower(nl), paste(label, nl)))
    v <- vapply(D, function(r) coherence(cor(t(cl$X[c(r$up, r$dn), ])), r$up, r$dn), 1); b <- balance(S, D, uu, dd, Lv = Lb, det = if (nl == "N2") cl$det else NULL)
    data.table(null = nl, obs = ob, p95 = unname(quantile(v, 0.95, na.rm = TRUE)), na_draws = sum(is.na(v)), counts_ok = b$counts_ok, bal_L = b$bal_L, bal_det = b$bal_det, jaccard = b$jaccard, n_merge = b$n_merge, fail = "") }
  if ("N1" %in% nulls) R$N1 <- one("N1", make_strata(U1, u, d, "d"), u, d, obs)
  if (any(c("N1L", "N2") %in% nulls) && (length(u2) < MIN_HALF || length(d2) < MIN_HALF)) { for (nl in intersect(c("N1L", "N2"), nulls)) R[[nl]] <- data.table(null = nl, obs = NA_real_, p95 = NA_real_, fail = "fewer than 20 genes per half on N2's universe") } else {
    if ("N1L" %in% nulls) R$N1L <- one("N1L", make_strata(U2, u2, d2, "d"), u2, d2, obs2)
    if ("N2" %in% nulls) R$N2 <- one("N2", make_strata(U2, u2, d2, c("d", "l"), Lv = Lv, k_l = K_L_T1), u2, d2, obs2, Lb = Lv) }
  rbindlist(R, fill = TRUE)[, `:=`(above = obs > p95, n_up = length(u), n_dn = length(d), skipped = FALSE)] }
t1_pass <- function(x) { x <- x[!is.na(above)]; n_ok <- nrow(x); n_above <- sum(x$above); list(n_ok = n_ok, n_above = n_above, pass = isTRUE(n_ok >= 1 && n_above >= min(6L, n_ok) && n_above / n_ok >= 0.75)) }   # 32:109
t1_cohort <- function(co, gsn, Lfun, label, nulls) {
  rbindlist(mcl(co$DON$donor, function(d) { j <- match(d, co$DON$donor); cl0 <- get_cells(co, d, SEED + OFF[["cap"]] + 1000L * match(co$key, names(COH)) + j)
    cl <- norm_cells(cl0$X); r <- t1_donor(cl, UPDN[[gsn]]$up, UPDN[[gsn]]$dn, Lfun(d), paste(label, d), nulls); r[, `:=`(donor = d, n_cells = ncol(cl$X), n_all = cl0$n_all)] })) }

## -- 8. --known: K0-K7 (section 6) --------------------------------------------------------------------------------------------------
t1_obs <- function(cl, up, dn) { u <- intersect(up, rownames(cl$X)); d <- intersect(dn, rownames(cl$X)); if (length(u) < MIN_HALF || length(d) < MIN_HALF) return(NA_real_)
  coherence(cor(t(cl$X[c(u, d), ])), u, d) }
perm_counts <- function(M, seed) { set.seed(seed); lib <- colSums(M); C <- t(t(M) / lib)   # CPM-scale permutation across columns per gene, rescaled to each column's library size
  Cp <- t(apply(C, 1, function(v) v[sample.int(length(v))])); out <- t(t(Cp) * lib); dimnames(out) <- dimnames(M); if (isTRUE(attr(M, "synthetic"))) syn(out) else out }
match_L_set <- function(Ugenes, Lv, up, dn, k, seed) {   # fresh test set: P370's per-L-quantile counts, drawn from genes outside P370
  g <- intersect(Ugenes, names(Lv)); q <- setNames(qbin(Lv[g], k), g); pool <- setdiff(g, c(UPDN$primary$up, UPDN$primary$dn)); set.seed(seed)
  pick <- function(h, excl) unlist(lapply(split(intersect(h, g), q[intersect(h, g)]), function(x) { b <- q[x[1]]; p <- setdiff(pool[q[pool] == b], excl); p[sample.int(length(p), min(length(p), length(x)))] }), use.names = FALSE)
  su <- pick(up, character(0)); list(up = su, dn = pick(dn, su)) }
if (MODE == "known") {
  up <- UPDN$primary$up; dn <- UPDN$primary$dn
  chk("K0b", "R build and collation locale of S8's run (R 4.4.3, en_US.UTF-8; table() orders S8's bin labels by locale)", c("4.4.3", "en_US.UTF-8"),
      c(paste(R.version$major, R.version$minor, sep = "."), Sys.getlocale("LC_COLLATE")))
  ts("[K0] S8's observed AML statistics"); S8T <- fread(file.path(TABWS, "ws_s8_tests.csv"))
  for (k in c("rep", "rep2")) { co <- CO[[k]]
    if (SMOKE) { rep_("K0", sprintf("%s: smoke world (S8's stored values not applicable)", k), "n/a"); next }
    o2 <- t2_stat(logz(co$PBT)$Z, up, dn); st <- S8T[set == k & gene_set == "primary" & test == "T2_opposite_movement", statistic]
    chk("K0", sprintf("%s: T2 rho equals ws_s8_tests.csv", k), st, o2, tol = 1e-12)
    SC8 <- fread(file.path(TABWS, sprintf("ws_s8_scores_%s.csv", k))); small <- co$DON[n_T <= CAP, donor]
    o1 <- vapply(small, function(d) t1_obs(norm_cells(read_cells(co, d, 1L)$X), up, dn), 1)
    chk("K0", sprintf("%s: T1 coherence of the %d patients with <= 500 T cells equals ws_s8_scores", k, length(small)), SC8[match(small, sample), t1_coherence], unname(o1), tol = 1e-12) }
  ts("[K1] rosters")
  S9 <- fread(file.path(TABWS, "ws_s9_scores_GSE120221.csv"), select = c("level", "id", "n_libraries", "n_kept_T"))[level == "donor"]
  D1 <- CO$oetjen$DON
  chk("K1", "Oetjen: 20 donors; libraries per donor and kept T cells equal S9's", TRUE, nrow(D1) == 20L && identical(as.integer(D1$n_lib), as.integer(S9[match(D1$donor, id), n_libraries])) && identical(as.integer(D1$n_T), as.integer(S9[match(D1$donor, id), n_kept_T])))
  chk("K1", "Griffin: donors BM1-BM5 with 1250 / 820 / 349 / 2080 / 533 kept T cells", c(1250L, 820L, 349L, 2080L, 533L), as.integer(CO$griffin$DON[match(paste0("BM", 1:5), donor), n_T]))
  ts("[K2] the loading is G01's")
  if (SMOKE) rep_("K2", "smoke world (G01's stored values not applicable)", "n/a") else {
    P185 <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"))
    S5 <- fread(file.path(TABWS, "ws_s5_scores.csv"))[compartment == "T"]; SAMP <- S5$sample; STAT <- as.character(S5$status)
    FG <- fread(file.path(TABWS, "ws_a1_g02_fold_genes.csv")); uni <- unique(FG$gene)
    Mo <- cbind(naive = Reduce(`+`, lapply(NAIVE, function(k) rowSums(CO$oetjen$fine[[k]]))), em = Reduce(`+`, lapply(EM, function(k) rowSums(CO$oetjen$fine[[k]])))); Lo <- logz_g01(Mo)
    LOAD <- LD$aml; EG <- log2(cpm(DGEList(matrix(rowSums(Mo[rownames(Lo), ]), dimnames = list(rownames(Lo), "x"))))[, 1] + 1)   # 39:539
    AXS <- fread(file.path(TABWS, "ws_a1_g01_axis.csv"))
    for (k in c("CD4 Naive", "CD8 Effector Memory 2")) {
      CN <- matrix(0, length(G185), length(SAMP), dimnames = list(G185, SAMP)); h <- intersect(colnames(P185$fine[[k]]), SAMP); CN[, h] <- P185$fine[[k]][, h, drop = FALSE]
      nc <- P185$cells[kept == TRUE & sample %in% SAMP & fine == k, .N, by = sample]; ncv <- setNames(nc$N[match(SAMP, nc$sample)], SAMP); ncv[is.na(ncv)] <- 0
      E <- SAMP[ncv >= WS_S4_MIN_CELLS]; sE <- STAT[match(E, SAMP)]; Xk <- logz_g01(CN[, E]); g <- intersect(rownames(Xk), names(LOAD))
      dec <- cut(rank(EG[g], ties.method = "first"), 10, labels = FALSE)
      dg <- rowMeans(Xk[g, sE == "AML", drop = FALSE]) - rowMeans(Xk[g, sE == "healthy", drop = FALSE])
      sta <- mean(vapply(1:10, function(q) suppressWarnings(cor(dg[dec == q], LOAD[g][dec == q], method = "spearman")), 1), na.rm = TRUE)
      gu <- intersect(uni, g); ur <- cor(dg[gu], LOAD[gu], method = "spearman"); a <- AXS[label == k]
      chk("K2", sprintf("%s: G01 5.5 statistic and union rho reproduced", k), c(a$statistic, a$union_rho_number_only), c(sta, ur), tol = 1e-12)
      chk("K2", sprintf("%s: n_genes and n_union", k), c(a$n_genes, a$n_union), c(length(g), length(gu))) }
    rm(P185) }
  allcols <- unlist(lapply(CO$oetjen$fine, colnames))
  chk("K2b", "leave-one-donor-out L with no donor left out equals the all-Oetjen L", TRUE, isTRUE(all.equal(L_of(CO$oetjen$fine, cols = unique(allcols)), LD$aml, tolerance = 0)))
  chk("K2b", "leave-one-donor-out drops all of a donor's libraries (C 3, S 4)", c(3L, 4L), vapply(c("C", "S"), function(d) length(intersect(unique(allcols), CO$oetjen$mp[donor == d, sample])), 1L) |> unname())
  ts("[K3-K5] T2 worlds on Oetjen pseudobulks")
  PB <- CO$oetjen$PBT; Lt2 <- LD$t2_oetjen
  k3 <- unlist(mcl(seq_len(N_K3), function(i) { M <- perm_counts(PB, bseed("worlds", paste("K3", i))); t2_run(M, up, dn, Lt2, paste("K3", i), nulls = "N1")$res$pass }))
  na_chk("K3", "null world", k3)
  chk_stat("K3", sprintf("null world: T2 (N1) passes in at most %d of %d permutations", ceiling(0.1 * N_K3), N_K3), sprintf("<= %d", ceiling(0.1 * N_K3)), sum(k3, na.rm = TRUE), sum(k3, na.rm = TRUE) <= ceiling(0.1 * N_K3))
  plant <- function(M, b, seed) { set.seed(seed); f <- rnorm(ncol(M)); u <- intersect(up, rownames(M)); d <- intersect(dn, rownames(M))
    M[u, ] <- t(t(M[u, ]) * 2^(b * f)); M[d, ] <- t(t(M[d, ]) * 2^(-b * f)); M }
  k4 <- lapply(c(0.3, 0.1, 0.05), function(b) unlist(mcl(seq_len(N_K4), function(i) { M <- plant(perm_counts(PB, bseed("worlds", paste("K4p", b, i))), b, bseed("worlds", paste("K4f", b, i))); t2_run(M, up, dn, Lt2, paste("K4", b, i), nulls = "N1")$res$pass })))
  na_chk("K4", "planted worlds (0.3, 0.1, 0.05)", unlist(k4))
  chk_stat("K4", sprintf("machinery: 0.3 log2 plant, T2 (N1) passes in at least %d of %d", ceiling(0.8 * N_K4), N_K4), sprintf(">= %d", ceiling(0.8 * N_K4)), sum(k4[[1]], na.rm = TRUE), sum(k4[[1]], na.rm = TRUE) >= ceiling(0.8 * N_K4))
  rep_("K4", "plants at 0.1 and 0.05 log2: T2 (N1) passes", c(sum(k4[[2]], na.rm = TRUE), sum(k4[[3]], na.rm = TRUE)))
  prof <- function(fine, cols) { pl <- function(ks) Reduce(`+`, lapply(ks, function(k) { M <- fine[[k]]; rowSums(M[, intersect(colnames(M), cols), drop = FALSE]) })); cbind(naive = pl(NAIVE), em = pl(EM)) }
  PR185 <- prof(CO$g185h$fine, CO$g185h$mp[donor %in% noT15, sample]); PROE <- prof(CO$oetjen$fine, unique(allcols))
  comp_world <- function(M, s_d, PRF) { cp <- t(t(PRF) / colSums(PRF)); g <- intersect(rownames(M), rownames(cp)); sb <- mean(s_d)
    for (j in seq_len(ncol(M))) { r <- (s_d[j] * cp[g, "naive"] + (1 - s_d[j]) * cp[g, "em"]) / (sb * cp[g, "naive"] + (1 - sb) * cp[g, "em"]); r[!is.finite(r)] <- 1; M[g, j] <- M[g, j] * r }; M }
  sr <- range(CO$oetjen$DON$naive_share)
  k5 <- rbindlist(mcl(seq_len(N_K5), function(i) { set.seed(bseed("worlds", paste("K5s", i))); s_d <- runif(ncol(PB), sr[1], sr[2])
    out <- list()
    for (src in c("g185", "oetjen")) { M <- comp_world(perm_counts(PB, bseed("worlds", paste("K5p", i))), s_d, if (src == "g185") PR185 else PROE)
      Ug <- logz(M)$U$gene; ts1 <- match_L_set(Ug, LD$aml, up, dn, 10L, bseed("worlds", paste("K5t", i)))
      r1 <- t2_run(M, ts1$up, ts1$dn, LD$aml, paste("K5", src, i), nulls = c("N1L", "N2"))$res
      out[[src]] <- data.table(world = i, src = src, set = "P370-matched", N1L = r1[null == "N1L", pass], N2 = r1[null == "N2", pass]) }
    gL <- intersect(Ug, names(LD$aml)); o <- order(LD$aml[gL]); ex <- list(up = tail(gL[o], 50), dn = head(gL[o], 50))
    M <- comp_world(perm_counts(PB, bseed("worlds", paste("K5p", i))), s_d, PR185); r2 <- t2_run(M, ex$up, ex$dn, LD$aml, paste("K5x", i), nulls = c("N1L", "N2"))$res
    rbind(rbindlist(out), data.table(world = i, src = "g185", set = "100-extreme", N1L = r2[null == "N1L", pass], N2 = r2[null == "N2", pass])) }))
  a5 <- k5[src == "g185" & set == "P370-matched"]; na_chk("K5", "composition worlds (gating set)", c(a5$N1L, a5$N2))
  chk_stat("K5", sprintf("composition world, T2: N1-L passes in at least %d of %d", ceiling(0.8 * N_K5), N_K5), sprintf(">= %d", ceiling(0.8 * N_K5)), sum(a5$N1L, na.rm = TRUE), sum(a5$N1L, na.rm = TRUE) >= ceiling(0.8 * N_K5))
  chk_stat("K5", sprintf("composition world, T2: N2 passes in at most %d of %d", floor(0.2 * N_K5), N_K5), sprintf("<= %d", floor(0.2 * N_K5)), sum(a5$N2, na.rm = TRUE), sum(a5$N2, na.rm = TRUE) <= floor(0.2 * N_K5))
  rep_("K5", "reported: same-source profiles (N1-L; N2) and the 100-extreme set (N1-L; N2); NA = not computed", c(k5[src == "oetjen", sum(N1L)], k5[src == "oetjen", sum(N2)], k5[set == "100-extreme", sum(N1L)], k5[set == "100-extreme", sum(N2)]))
  ts("[K5b, K7] T1 worlds in Oetjen donor K")
  coK <- CO$oetjen; clK <- get_cells(coK, "K", SEED + OFF[["cap"]] + 1000L * match("oetjen", names(COH)) + match("K", coK$DON$donor)); LK <- lodo_L("K")
  XK <- as.matrix(clK$X[Matrix::rowSums(clK$X) > 0, , drop = FALSE]); labK <- clK$lab   # all-zero genes dropped (memory; they stay zero and never pass DET_MIN)
  k5b <- rbindlist(mcl(seq_len(N_K5), function(i) { set.seed(bseed("worlds", paste("K5b", i)))
    Xp <- t(apply(XK, 1, function(v) v[sample.int(length(v))])); dimnames(Xp) <- dimnames(XK)
    w <- ifelse(labK %in% EM, 1, ifelse(labK %in% NAIVE, 0, NA_real_)); wb <- mean(w, na.rm = TRUE); w[is.na(w)] <- wb
    pn <- rowSums(XK[, labK %in% NAIVE, drop = FALSE]); pe <- rowSums(XK[, labK %in% EM, drop = FALSE]); pn <- pn / sum(pn); pe <- pe / sum(pe)
    R <- outer(pe, w) + outer(pn, 1 - w); R <- R / (wb * pe + (1 - wb) * pn); R[!is.finite(R)] <- 1; R <- R / pmax(apply(R, 1, max), 1e-12)
    Xt <- matrix(rbinom(length(Xp), as.vector(Xp), as.vector(R)), nrow(Xp), dimnames = dimnames(Xp)); Xt <- Matrix(Xt, sparse = TRUE); if (isTRUE(attr(clK$X, "synthetic"))) Xt <- syn(Xt)
    cl <- norm_cells(Xt); g2 <- setdiff(intersect(names(cl$det), names(LK)), black_for(names(cl$det)))
    tsb <- match_L_set(g2, LK, up, dn, 5L, bseed("worlds", paste("K5bt", i)))
    r <- t1_donor(cl, tsb$up, tsb$dn, LK, paste("K5b", i), nulls = c("N1L", "N2")); data.table(world = i, N1L = r[null == "N1L", above][1], N2 = r[null == "N2", above][1]) }))   # NA = not computed
  na_chk("K5b", "T1 composition worlds", c(k5b$N1L, k5b$N2))
  chk_stat("K5b", sprintf("composition world, T1: N1-L above P95 in at least %d of %d", ceiling(0.8 * N_K5), N_K5), sprintf(">= %d", ceiling(0.8 * N_K5)), sum(k5b$N1L, na.rm = TRUE), sum(k5b$N1L, na.rm = TRUE) >= ceiling(0.8 * N_K5))
  chk_stat("K5b", sprintf("composition world, T1: N2 above P95 in at most %d of %d", floor(0.2 * N_K5), N_K5), sprintf("<= %d", floor(0.2 * N_K5)), sum(k5b$N2, na.rm = TRUE), sum(k5b$N2, na.rm = TRUE) <= floor(0.2 * N_K5))
  ts("[K6] GSE227903 AML_Dx subset")
  m27 <- fread(file.path(TABWS, "ws_s4_pseudobulk_meta_GSE227903.csv")); dx <- unique(m27[role == "AML_Dx" & paediatric == FALSE, sample])
  chk("K6", "GSE227903: exactly the 9 AML_Dx samples, 7 eligible, z-scores over those 7 only", "9/7/7", sprintf("%d/%d/%d", length(intersect(CO$rep2$mp$sample, dx)) * as.integer(setequal(CO$rep2$mp$sample, dx)), nrow(CO$rep2$DON), ncol(CO$rep2$PBT)))
  k7 <- unlist(mcl(seq_len(N_K7), function(i) { cl <- norm_cells(clK$X); set.seed(bseed("k7", paste("K7", i))); X <- t(apply(cl$X, 1, function(v) v[sample.int(length(v))])); dimnames(X) <- dimnames(cl$X)
    t1_donor(list(X = X, det = cl$det), up, dn, LK, paste("K7", i), nulls = "N1")[null == "N1", above][1] }))   # NA = not computed
  na_chk("K7", "T1 null worlds", k7)
  chk_stat("K7", sprintf("null world for T1 (N1): above P95 in at most %d of %d", ceiling(0.15 * N_K7), N_K7), sprintf("<= %d", ceiling(0.15 * N_K7)), sum(k7, na.rm = TRUE), sum(k7, na.rm = TRUE) <= ceiling(0.15 * N_K7))
  C <- rbindlist(CHK)[, script_md5 := MD5]; fwrite(C, wtab("checks"))
  ts(sprintf("[known] %d rows; deterministic FAIL %d; statistical FAIL %d", nrow(C), C[state == "FAIL_deterministic", .N], C[state == "FAIL_statistical", .N])); quit(status = 0)
}

## -- 9. --run: 3.1-3.3, section 4, K8 (no statistic is printed before K8 has a state) ---------------------------------------------------
if (MODE == "run") {
  CK <- fread(file.path(ODIR, paste0("ws_a1_g09_checks", SFX, ".csv")))
  if (!SMOKE && !identical(unique(CK$script_md5), MD5)) stop("the checks table was written by another version of this script; rerun --known")
  if (any(CK$state == "FAIL_deterministic")) stop("a deterministic known answer failed: ", paste(CK[state == "FAIL_deterministic", id], collapse = ", "))
  ROWS <- list(); DONOR <- list()
  addT2 <- function(r, cohort, gsn, view) ROWS[[length(ROWS) + 1L]] <<- r$res[, `:=`(cohort = cohort, gene_set = gsn, view = view)]
  addT1 <- function(x, cohort, gsn, view) { x <- copy(x); for (cn in c("counts_ok", "na_draws", "bal_L", "bal_det", "n_merge", "fail")) if (!cn %in% names(x)) x[, (cn) := NA]
    DONOR[[length(DONOR) + 1L]] <<- copy(x)[, `:=`(cohort = cohort, gene_set = gsn, view = view)]
    ROWS[[length(ROWS) + 1L]] <<- x[, { p <- t1_pass(.SD); .(test = "T1", n = p$n_ok, statistic = p$n_above / max(1, p$n_ok), n_above = p$n_above, pass = p$pass,
      med_obs = median(obs, na.rm = TRUE), med_p95 = median(p95, na.rm = TRUE),   # 3.3
      counts_ok = all(counts_ok %in% c(TRUE, NA)), na_draws = sum(na_draws, na.rm = TRUE), bal_L = max(c(bal_L, -Inf), na.rm = TRUE), bal_det = max(c(bal_det, -Inf), na.rm = TRUE),
      n_merge = sum(n_merge, na.rm = TRUE), fail = paste(unique(na.omit(fail[nzchar(fail)])), collapse = "; ")) }, by = null][, `:=`(cohort = cohort, gene_set = gsn, view = view)] }
  for (gsn in c("primary", "full")) { up <- UPDN[[gsn]]$up; dn <- UPDN[[gsn]]$dn; nl <- if (gsn == "primary") c("N0", "N1", "N1L", "N2") else c("N1", "N2")
    ts(sprintf("[run] %s: Oetjen (a)", gsn))
    r <- t2_run(CO$oetjen$PBT, up, dn, LD$t2_oetjen, paste("oetjen", gsn), nulls = nl); addT2(r, "oetjen", gsn, "primary"); if (gsn == "primary") ZOE <- r$Z
    addT1(t1_cohort(CO$oetjen, gsn, lodo_L, paste("oetjen", gsn), nl), "oetjen", gsn, "primary")
    ts(sprintf("[run] %s: AML cohorts (b)", gsn))
    for (k in c("rep", "rep2")) { addT2(t2_run(CO[[k]]$PBT, up, dn, LD$aml, paste(k, gsn), nulls = nl), k, gsn, "primary")
      addT1(t1_cohort(CO[[k]], gsn, function(d) LD$aml, paste(k, gsn), nl), k, gsn, "primary") } }
  up <- UPDN$primary$up; dn <- UPDN$primary$dn
  ts("[run] sensitivities and reported cohorts")
  addT2(t2_run(CO$oetjen$PBT, up, dn, LD$t2_oetjen_all10, "oetjen L185all", nulls = "N2"), "oetjen", "primary", "4.2 L from all 10 GSE185381 healthy donors")
  addT2(t2_run(CO$oetjen$PBT, up, dn, LD$aml, "oetjen Loetjen", nulls = "N2"), "oetjen", "primary", "4.2 L from all 20 Oetjen donors (circular)")
  g185 <- t1_cohort(CO$g185h, "primary", function(d) LD$aml, "g185h", c("N0", "N1"))
  addT2(t2_run(CO$g185h$PBT, up, dn, LD$aml, "g185h", nulls = c("N0", "N1")), "g185h", "primary", "reported: all 10"); addT1(g185, "g185h", "primary", "reported: all 10")
  addT2(t2_run(cols(CO$g185h$PBT, noT15), up, dn, LD$aml, "g185h noT15", nulls = c("N0", "N1")), "g185h", "primary", "reported: without Control1-5")
  addT1(g185[donor %in% noT15], "g185h", "primary", "reported: without Control1-5")
  for (k in c("griffin", "gse276360", "seqwell_h")) addT1(t1_cohort(CO[[k]], "primary", function(d) LD$aml, k, c("N0", "N1")), k, "primary", "supporting: T1 only")
  ll <- rbindlist(lapply(list(c("oetjen", "C"), c("oetjen", "S"), c("griffin", paste0("BM", 1:5))), function(v) { co <- CO[[v[1]]]
    rbindlist(lapply(v[-1], function(d) { lib <- co$mp[donor == d][order(-n_T)][1, sample]
      cl0 <- if (SMOKE) co$smoke_cells(d, 1L) else read_cells(co, d, SEED + OFF[["cap"]] + 7L, lib_only = lib)
      t1_donor(norm_cells(cl0$X), up, dn, LD$aml, paste("largestlib", v[1], d), nulls = "N1")[, `:=`(donor = d, cohort = v[1], library = lib)] })) }))
  DONOR[[length(DONOR) + 1L]] <- ll[, `:=`(gene_set = "primary", view = "4.4 largest library only")]
  ts("[run] 8-of-20 draws (descriptive)")
  set.seed(SEED + OFF[["d8"]]); DD <- lapply(seq_len(N_D8), function(i) sort(sample(CO$oetjen$DON$donor, 8L)))
  t1flags <- rbindlist(DONOR, fill = TRUE)[cohort == "oetjen" & gene_set == "primary" & view == "primary" & null == "N1", .(donor, above)]
  D8 <- rbindlist(mcl(seq_len(N_D8), function(i) { dd <- DD[[i]]; r <- t2_run(cols(CO$oetjen$PBT, dd), up, dn, LD$t2_oetjen, paste("d8", i), nulls = "N1", nr = N_D8_NULL)$res
    data.table(draw = i, donors = paste(dd, collapse = ";"), T2_rho = r$statistic, T2_P5 = r$P5, T2_pass = r$pass, T1_pass = t1_pass(t1flags[donor %in% dd])$pass) }))
  TT <- rbindlist(ROWS, fill = TRUE); DV <- rbindlist(DONOR, fill = TRUE)
  ## K8: counts and NA (deterministic), balance (statistical, per run and per donor) -- checked before anything is written or printed
  k8_det <- TT[null %in% c("N1", "N1L", "N2") & (counts_ok %in% FALSE | (!is.na(na_draws) & na_draws > 0) | grepl("^stratum|; stratum", fail))]
  BAL <- TT[null == "N2", .(cohort, gene_set, view, test, bal_L, bal_det,   # a balance that could not be computed (NA, -Inf) is not readable; T2 has no detection term
                            readable = is.finite(bal_L) & bal_L < BAL_MAX & (test == "T2" | (is.finite(bal_det) & bal_det < BAL_MAX)))]
  put(BAL[, k8_deterministic_ok := nrow(k8_det) == 0L], "balance")
  if (nrow(k8_det)) stop("K8 deterministic FAIL (per-stratum counts, NA draws or a stratum that could not be filled): ", paste(unique(paste(k8_det$cohort, k8_det$test, k8_det$null, k8_det$fail)), collapse = " | "))
  ## naive-share spread, age against the halves, the published S8 null (K0b)
  SPR <- rbindlist(lapply(names(CO), function(k) CO[[k]]$DON[, .(cohort = k, n = .N, naive_sd = sd(naive_share), naive_min = min(naive_share), naive_max = max(naive_share))]))
  ag <- CO$oetjen$DON[match(colnames(ZOE), donor), age]; uz <- colMeans(ZOE[intersect(up, rownames(ZOE)), , drop = FALSE]); dz <- colMeans(ZOE[intersect(dn, rownames(ZOE)), , drop = FALSE])
  AGE <- data.table(cohort = "oetjen", test = "age", null = c("up half", "down half"), statistic = c(cor(uz, ag, method = "spearman"), cor(dz, ag, method = "spearman")), view = "3.3 descriptive", gene_set = "primary")
  S8T <- fread(file.path(TABWS, "ws_s8_tests.csv"))
  K0B <- rbindlist(lapply(c("rep", "rep2"), function(k) { tl <- c(T2 = "T2_opposite_movement", T1 = "T1_within_patient_coherence")   # smoke: S8's published values against synthetic re-runs (exercises the code and 7.2(iv) only)
    rbindlist(lapply(names(tl), function(tt) { s8 <- S8T[set == k & gene_set == "primary" & test == tl[[tt]]]; rr <- TT[cohort == k & gene_set == "primary" & view == "primary" & test == tt & null == "N0"]
      stopifnot(nrow(s8) == 1L, nrow(rr) == 1L)
      data.table(cohort = k, test = tt, published_pass = s8$pass, rerun_pass = rr$pass,   # 7.2(iv): "(published p, re-run p)" = P5 for T2, k of n for T1
                 published_value = if (tt == "T2") sprintf("P5 %.2f", s8$null_P5) else sub("^([0-9]+ of [0-9]+).*$", "\\1", s8$detail),
                 rerun_value = if (tt == "T2") sprintf("P5 %.2f", rr$P5) else sprintf("%d of %d", as.integer(rr$n_above), as.integer(rr$n)),
                 r_version = R.version.string, collate = Sys.getlocale("LC_COLLATE")) })) }))
  if (nrow(K0B)) K0B[, equal := published_pass == rerun_pass]
  put(rbind(TT, AGE, fill = TRUE), "tests"); put(D8, "draws8"); fwrite(DV[, script_md5 := MD5], wg("donor_values.csv"))
  fwrite(SPR[, script_md5 := MD5], wg("naive_share_spread.csv")); fwrite(K0B[, script_md5 := MD5], wg("k0b.csv"))
  ts(sprintf("[run] done: %d test rows; K8 balance readable %d of %d", nrow(TT), BAL[readable == TRUE, .N], nrow(BAL))); quit(status = 0)
}

## -- 10. --combine: section 7 ----------------------------------------------------------------------------------------------------------
if (MODE == "combine") {
  CK <- rtab("checks"); TT <- rtab("tests"); BAL <- rtab("balance"); D8 <- rtab("draws8")
  SPR <- fread(file.path(WDIR, paste0("naive_share_spread", SFX, ".csv"))); K0B <- fread(file.path(WDIR, paste0("k0b", SFX, ".csv")))
  if (!SMOKE && !identical(unique(c(CK$script_md5, TT$script_md5, BAL$script_md5, D8$script_md5, SPR$script_md5, K0B$script_md5)), MD5)) stop("tables from different script versions")
  okK <- function(ids) all(CK[id %in% ids & kind == "statistical", state] == "PASS")
  f2 <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 2))
  P <- TT[gene_set == "primary" & view == "primary"]
  g <- function(co, te, nl, col) P[cohort == co & test == te & null == nl][[col]][1]
  narrower <- isTRUE(SPR[cohort == "oetjen", naive_sd] < min(SPR[cohort %in% c("rep", "rep2"), naive_sd]))
  CTX <- "Structure is covariance, not level: this does not change the A1-G02 numbers. It is not evidence for or against composition below the label level (A1-G01 11.5). The direction stays shown in one study (A1-G04)."
  ## (b) per AML cohort and test. Precedence (dated note 11.4): (iii) first; then N1-L not computable or failing; then N2 not computable;
  ## then not readable (K5/K5b or K8 balance, the failing check named); then (i) / (ii)
  CN <- c(rep = "GSE289435", rep2 = "GSE227903")
  comp <- function(co, te, nl) { r <- P[cohort == co & test == te & null == nl]   # computed: one row, no fail text, a statistic (T2) or >= 1 donor (T1)
    if (nrow(r) != 1L || nzchar(paste(na.omit(r$fail), collapse = ""))) return(FALSE)
    if (te == "T2") is.finite(r$statistic) && is.finite(r$P5) else isTRUE(r$n >= 1) }
  why_nc <- function(co, te, nl) { f <- paste(na.omit(P[cohort == co & test == te & null == nl, fail]), collapse = ""); if (nzchar(f)) f else "too few genes for the test" }
  B <- rbindlist(lapply(c("rep", "rep2"), function(co) rbindlist(lapply(c("T2", "T1"), function(te) {
    cname <- CN[[co]]; kid <- if (te == "T2") "K5" else "K5b"
    p1 <- isTRUE(g(co, te, "N1", "pass")); p1l <- isTRUE(g(co, te, "N1L", "pass")); p2 <- isTRUE(g(co, te, "N2", "pass"))
    balok <- isTRUE(BAL[cohort == co & test == te & gene_set == "primary" & view == "primary", readable][1])
    why <- c(if (!okK(kid)) kid, if (!balok) "K8 balance")
    case <- if (!p1) "iii" else if (!comp(co, te, "N1L")) "not computable" else if (!p1l) "n1l_fail" else if (!comp(co, te, "N2")) "not computable" else
      if (length(why)) "not readable" else if (p2) "i" else "ii"
    s <- switch(case,
      iii = sprintf("Against random sets matched separately for each half, S8's %s in %s does not pass; S8's published pass used a null that split the random halves by bin order.", te, cname),
      `not computable` = sprintf("Whether S8's %s in %s goes beyond the naive-to-effector axis is not readable: the loading-matched null could not be computed (%s).", te, cname,
                                 why_nc(co, te, if (!comp(co, te, "N1L")) "N1L" else "N2")),
      n1l_fail = sprintf("S8's %s in %s does not pass on the loading-matched universe even before matching on the loading (N1-L); no (b) reading is made.", te, cname),
      `not readable` = sprintf("Whether S8's %s in %s goes beyond the naive-to-effector axis is not readable: the loading-matched null failed its known answer (%s).", te, cname, paste(why, collapse = ", ")),
      i = if (te == "T2") sprintf("S8's T2 in %s exceeds random gene sets matched separately for each half on expression, detection and the naive-to-effector loading.", cname) else
            sprintf("S8's T1 in %s exceeds random gene sets matched separately for each half on detection and the naive-to-effector loading.", cname),
      ii = sprintf("S8's %s in %s is consistent with the naive-to-effector differentiation axis.", te, cname))
    data.table(cohort = co, test = te, case = case, sentence = paste(s, "This concerns the covariance of S8's structure, not the within-label shift of A1-G01 (5.5 did not fire).")) }))))
  ## (a) Oetjen. Only readable passed tests count; the (ii) sentence names only the readable tests (dated note 11.4)
  t2p <- isTRUE(g("oetjen", "T2", "N1", "pass")); t1p <- isTRUE(g("oetjen", "T1", "N1", "pass")); ra <- c(T2 = okK(c("K3", "K4")), T1 = okK("K7"))
  r2 <- f2(g("oetjen", "T2", "N1", "statistic")); p5 <- f2(g("oetjen", "T2", "N1", "P5")); k1 <- as.integer(g("oetjen", "T1", "N1", "n_above")); n1 <- as.integer(g("oetjen", "T1", "N1", "n"))
  txt <- c(T2 = sprintf("T2: the up-half and down-half means move oppositely across donors, rho %s, random P5 %s", r2, p5),
           T1 = sprintf("T1: they co-vary within donors, %d of %d donors above their random P95", k1, n1))
  val <- c(T2 = sprintf("T2 rho %s, random P5 %s", r2, p5), T1 = sprintf("T1 %d of %d donors", k1, n1))
  pz <- c(if (t2p && ra[["T2"]]) "T2", if (t1p && ra[["T1"]]) "T1"); rd <- names(ra)[ra]
  A <- if (!all(ra)) sprintf("Whether the healthy Oetjen donors show S8's %s is not readable: its known answer failed (%s). ", paste(names(ra)[!ra], collapse = " and "),
                             paste(c(if (!ra[["T2"]]) "K3/K4", if (!ra[["T1"]]) "K7"), collapse = ", ")) else ""
  if (length(pz)) { hit <- intersect(B[cohort == "rep" & case == "iii", test], pz)
    A <- paste0(A, sprintf("Among the T cells of one healthy marrow cohort (Oetjen 2018, 20 donors, 10x 3' v2, one protocol; T labels in-sample to BoneMarrowMap), the axis genes show the structure that S8 tested (%s). ", paste(txt[pz], collapse = "; ")),
                if (length(hit)) sprintf("This structure is therefore not specific to AML; S8's AML pass for %s did not survive the per-half null (7.2(iii)).", paste(hit, collapse = " and ")) else
                  "S8 therefore shows that the axis's structure is reproducible on 10x Chromium cohorts, not that it is specific to AML.")
  } else if (length(rd)) {
    A <- paste0(A, sprintf("Among 20 healthy Oetjen donors (one cohort, one protocol), this structure was not detected (%s). This is not evidence that the structure is specific to AML: depth, cell numbers, protocol and age differ from the AML cohorts%s, and the power to detect %s among healthy donors was not established.",
                           paste(val[rd], collapse = "; "), if (narrower) ", and the spread of T-subset composition is narrower among them" else "", paste(rd, collapse = " or "))) }
  A <- paste(A, CTX)
  joint <- if (!length(pz) && all(ra) && any(B$case == "ii") && narrower) "The AML structure is consistent with the naive-to-effector axis, so its absence among healthy donors may reflect a narrower spread of T-subset composition." else ""
  K0S <- if (nrow(K0B) && "equal" %in% names(K0B) && any(!K0B$equal)) paste(K0B[equal == FALSE, sprintf("Re-running S8's own null changed its pass flag for %s, %s (published %s, re-run %s); Monte-Carlo error of the null.", test, CN[cohort], published_value, rerun_value)], collapse = " ") else ""
  COV <- unique(TT[test == "T2" & gene_set == "primary" & !is.na(coverage), .(cohort, view, coverage)])   # 3.3: S8 requires >= 0.8 (reported)
  covtxt <- if (COV[coverage < 0.8, .N]) paste("P370 coverage below S8's 0.8 in:", COV[coverage < 0.8, paste(sprintf("%s (%s) %.2f", cohort, view, coverage), collapse = "; ")]) else
    sprintf("P370 coverage >= 0.8 in every T2 run (lowest %.2f)", min(COV$coverage))
  CTX74 <- paste("The 32 ambient-prone genes stay excluded from every statement; F402 appears in tables only, and no P370/F402 difference is attributed to the 32.",
                 "S8 is shown on 10x Chromium cohorts; not reproduced on Seq-Well GSE116256; its T3 is uninformative.",
                 "Oetjen is one cohort and one protocol (10x 3' v2), with labels in-sample to BoneMarrowMap. GSE185381's Control1-5 are T-ALL lanes and Control4 is female by expression.",
                 "A1-G01 category: \"beyond T-subset composition, detected inside the main labels\", with its 11.5 sentence: \"In GSE185381 (15 adult AML, 10 healthy donors), the mix of BoneMarrowMap T subtypes accounts for about 5% of the cross-fitted T-imprint gap (C 0.045, 95% CI -0.100 to 0.191; C_H 0.051, -0.107 to 0.193; C_A 0.039, -0.104 to 0.196), and a shift is detected inside the four main T labels (S_w 0.664, 95% CI 0.497 to 0.772; a sign check, not a share of the gap). The intervals are descriptive, with the fold lists held fixed. With gamma-delta and MAIT cells as separate labels the composition share is 0.057 (the flags are lower bounds).\" (and A1-G01's section-0 limits).",
                 "\"Samples were taken at diagnosis; acute illness, infection, transfusion and medication are unrecorded, and no non-AML patient comparator was run, so 'AML-associated' includes being acutely ill. The commercial-versus-biobank source of the arms is not separated from disease.\" (A1-G07 7.7) CMV status is unrecorded (CMV expands TEMRA-like cells inside the CD8 EM labels and rises with age).",
                 "S9 bounds the age share (0.102, upper bound 0.463).")
  SUM <- data.table(item = c("(a) healthy donors (7.1)", paste0("(b) ", unname(CN[B$cohort]), " ", B$test, " (7.2)"), "joint (7.2)", "K0b (7.2(iv))", "known answers", "coverage (3.3)", "context (7.4)", "8-of-20 draws (descriptive; never compared with the AML passes)"),
                    value = c(A, B$sentence, joint, K0S, sprintf("deterministic FAIL %d; statistical FAIL %s", CK[state == "FAIL_deterministic", .N], paste(unique(CK[state == "FAIL_statistical", id]), collapse = ", ")),
                              covtxt, CTX74, sprintf("T2 passed in %d of %d draws; T1 in %d of %d", D8[, sum(T2_pass)], nrow(D8), D8[, sum(T1_pass)], nrow(D8))))
  put(SUM, "summary")
  BLOCK <- c("<!-- A1-G09 start -->", sprintf("## Is S8's structure specific to AML? (A1-G09 stage 1, prereg-a1-g09-v1) -- %s", format(Sys.Date())), "",
             sprintf("Script: scripts/ws/45_a1_g09.R (md5 %s). P370; Oetjen GSE120221 (20 healthy donors) decides (a); GSE289435 and GSE227903 AML_Dx for (b).", MD5), "",
             "| item | reading |", "|---|---|", SUM[nzchar(value), sprintf("| %s | %s |", item, gsub("\\|", "/", value))], "",
             "**Reading: stage 1 of G09 (PREREGISTRATION_A1_G09.md section 7; its never-written list is section 7.3). Stage 2 (power, Seq-Well) is a separate amendment.**",
             "<!-- A1-G09 end -->")
  ## 7.2(iii): dated erratum after the "Reading" paragraph of block A1-S8, outside its markers; lists all four N1 results
  n1txt <- function(co, te) { pa <- if (isTRUE(g(co, te, "N1", "pass"))) "pass" else "fail"
    if (te == "T2") sprintf("T2, %s, rho %s, P5 %s: %s", CN[[co]], f2(g(co, "T2", "N1", "statistic")), f2(g(co, "T2", "N1", "P5")), pa) else
      sprintf("T1, %s, %d of %d: %s", CN[[co]], as.integer(g(co, "T1", "N1", "n_above")), as.integer(g(co, "T1", "N1", "n")), pa) }
  ERRN <- if (any(B$case == "iii")) { gate <- isTRUE(g("rep", "T2", "N1", "pass")) && isTRUE(g("rep", "T1", "N1", "pass"))
    sprintf("Erratum (A1-G09, %s): S8's random-set null drew one gene per imprint gene and split the list by bin order, so the random halves were not matched separately. Against a per-half null (A1-G09 N1): %s. GATE_S8 under the per-half null: %s (T2 and T1 on GSE289435).",
            format(Sys.Date()), paste(c(n1txt("rep", "T2"), n1txt("rep", "T1"), n1txt("rep2", "T2"), n1txt("rep2", "T1")), collapse = "; "), if (gate) "PASS" else "FAIL") } else NULL
  insert_erratum <- function(GL) { e8 <- grep("<!-- A1-S8 end -->", GL, fixed = TRUE); r8 <- grep("^Reading of the A1-S8 block", GL)
    stopifnot(length(e8) == 1L, length(r8) == 1L, r8 > e8, !any(grepl("A1-S8 erratum (A1-G09)", GL, fixed = TRUE)))
    nb <- which(seq_along(GL) > r8 & !nzchar(trimws(GL)))[1]; stopifnot(!is.na(nb))   # first blank line after the Reading paragraph
    append(GL, c("", "<!-- A1-S8 erratum (A1-G09) -->", ERRN), after = nb - 1L) }
  GL <- readLines(GATES); if (!is.null(ERRN)) GL <- insert_erratum(GL)
  if (SMOKE) { writeLines(c(BLOCK, if (!is.null(ERRN)) c("", "[erratum note after A1-S8]", ERRN)), file.path(WDIR, "gates_block_smoke.md"))
    writeLines(c(GL, "", BLOCK), file.path(WDIR, "GATES_copy_smoke.md"))   # smoke: the real insert path, on a copy outside the repo
  } else writeLines(c(GL, "", BLOCK), GATES)
  ts("[combine] done"); writeLines(c(BLOCK, if (!is.null(ERRN)) ERRN))
}
