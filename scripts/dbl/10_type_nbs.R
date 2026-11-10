#!/usr/bin/env Rscript
# 10_type_nbs.R ----
# Doublet-pairing SMOKE, stage 10 (user decision 2026-11-10: "try a Neighbor-seq-style pair classifier"; the marker idea enters
# as a T-identity gate). Lineage level as in 08/09: LEUK = five states + LEUK_other; detection = devel scDblFinder (stage 04).
# CLASSIFIER (Neighbor-seq defaults, see config_dbl.R section 9): per lane and held-out fold, the lineage types of the lane's
#   NON-held-out class-A cells (LEUK, T, Tamb, NK, B, MonoDC, Ery, Other; >= 5 cells) -> top-50 positive markers per type ->
#   100 artificial barcodes per class (every pair type incl. homotypic = summed raw counts of one cell of each type; plus 100
#   cells of each type) -> stratified 80/20 split -> xgboost random forest (200 trees, depth 20, subsample 0.632, colsample
#   sqrt(p)/p by tree and by node, one round) on raw marker counts -> each test barcode gets its most probable class.
#   Deviation: singlet classes sample with replacement only for types with < 100 cells (Neighbor-seq switches all of them).
# Pre-listed variants (no others are tried):
#   N1   classifier alone: most probable class is the T+LEUK pair                                 (no scDblFinder at all)
#   N2   devel scDblFinder called the barcode AND N1
#   M1   stage 09's best (V5 = V1a OR V3) AND T-identity UMI (CD3D/E/G, TRAC, TRBC1/2) >= 1 in the barcode  (marker gate)
#   N2M  N2 AND the same T-identity gate
#   V5   reference only, never selectable: must reproduce stage 09's tuning numbers exactly (validity check)
#   each with stage 04 runs made with B/C labels from 02 (knn) or 02b (bmm); the classifier itself sees class-A cells only.
# Natural recall / precision per qualifying patient: stage 08/09's random co-encapsulation model (block copied verbatim).
# SELECTION (tuning folds 1-2; fixed 2026-11-10 before any stage-10 number was computed): the (variant, labels) among N1, N2,
#   M1, N2M x {knn, bmm} with the largest min over qualifying patients of min(natural recall, natural precision); ties -> N1 <
#   N2 < M1 < N2M, knn before bmm. Frozen to dbl_nbs_choice_smoke.csv. If it is below 0.7 on the tuning folds, fold 3 is NOT
#   run: verdict "the pair classifier does not reach the bar either (already on the tuning folds); stop".
# READING RULE (confirmation, only if the bar was reached): stage 04 --folds=3 --heldout-only --labels=<chosen> --tag=f3_<labels>,
#   then --mode=confirm. Passes if natural recall >= 0.7 and natural precision >= 0.7 for every qualifying patient in both
#   lanes -> "lineage-level typing is good enough to design the O/E test"; otherwise "stop (no further tuning)".
# VALIDITY (the result is not read if any check fails): rebuilt injection == stage 04's rows/types/LEUK share; V5 reproduces
#   stage 09 (tune mode); training cells disjoint from the held-out fold; shuffled-label classifier (same marker genes, class-A
#   lineage labels permuted before the artificial barcodes are built) recalls <= 0.10 of the held-out synthetic T-LEUK-lineage
#   doublets; every barcode given to the classifier is synthetic or a re-entered singlet (no real B/C barcode is typed).
# MASKING: re-entered singlets are real cells -- for them only pooled per-compartment counts of each variant's flag leave this
#   script. The classifier is never applied to class B/C barcodes. The confusion table covers synthetic doublets only.
# INPUT  : DBL_TAB_DIR/dbl_typing_injected_devel_<tag>_smoke.csv.gz, DBL_LANE_DIR/<lane>__calib_devel_<tag>_smoke.rds (04),
#          DBL_LANE_DIR/<lane>.rds, <lane>__labels.csv.gz, <lane>__labels_bmm.csv.gz; tune: DBL_TAB_DIR/dbl_recall_tune_smoke.csv (09)
# OUTPUT : tune:    DBL_TAB_DIR/dbl_nbs_tune_smoke.csv, dbl_nbs_choice_smoke.csv
#          confirm: DBL_TAB_DIR/dbl_nbs_confirm_smoke.csv, dbl_nbs_verdict_smoke.csv
#          both:    dbl_nbs_classifier_<mode>_smoke.csv (per lane x fold x world: genes, classes, own-test accuracy, held-out
#                   T-LEUK recall), dbl_nbs_confusion_<mode>_smoke.csv (synthetic doublet type -> predicted class), dbl_checks_10_<mode>_smoke.csv
# Usage  : Rscript scripts/dbl/10_type_nbs.R --mode=tune --runs=tune2_knn:knn,tune2_bmm:bmm
#          Rscript scripts/dbl/10_type_nbs.R --mode=confirm --runs=f3_<labels>:<labels>

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here); library(parallel); library(nnls)
  library(Seurat); library(SeuratObject); library(xgboost) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
ARGS <- commandArgs(TRUE)
argv <- function(key, default) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), ARGS, value = TRUE)); if (length(v)) v else default }
MODE <- argv("mode", ""); stopifnot(MODE %in% c("tune", "confirm"))
RUNSPEC <- rbindlist(lapply(strsplit(strsplit(argv("runs", ""), ",")[[1]], ":"), function(z) data.table(tag = z[1], labels = z[2])))
stopifnot(nrow(RUNSPEC) > 0, all(RUNSPEC$labels %in% c("knn", "bmm")))
CORES <- as.integer(argv("cores", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
ts  <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-84s expect %-8s got %-8s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
STATES <- DBL_LEUK_STATES; LIN <- c(STATES, "LEUK_other"); WIN <- DBL_AMBIENT_WINDOW
NEWV <- c("N1", "N2", "M1", "N2M")                 # selectable, in tie-break order
VARS <- c("V5", NEWV)                              # V5 = stage 09 reference
TL   <- "LEUK|T"                                   # the classifier's T+LEUK pair class (pair names are sorted)
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

## -- Neighbor-seq core, re-implemented (prep_cell_mat / artificial_multiplets / multiplet_rf / xgpred) ----
lineage_of <- function(comp) fifelse(comp %in% LIN, "LEUK", comp)   # LEUK = five states + LEUK_other
nbs_markers <- function(M, cells, types) {          # top-N positive markers per type on log-normalized counts
  so <- NormalizeData(CreateSeuratObject(M[, cells]), verbose = FALSE)
  Idents(so) <- factor(types)
  mk <- as.data.table(FindAllMarkers(so, logfc.threshold = 1, only.pos = TRUE, min.pct = 0.25, max.cells.per.ident = 200, verbose = FALSE))
  unique(mk[order(cluster, -avg_log2FC)][, head(.SD, DBL_NBS_TOPN), by = cluster]$gene) }
nbs_artificial <- function(M, G, cells, types, n) { # n summed-count barcodes per pair class (incl. homotypic) + n cells per type
  idx <- split(cells, types); u <- sort(names(idx)); pr <- CJ(a = u, b = u)[a <= b]   # CJ = all combinations; a <= b keeps each pair once
  draw <- function(tp, k, rep = TRUE) idx[[tp]][sample.int(length(idx[[tp]]), k, replace = rep)]
  pa <- unlist(lapply(pr$a, draw, k = n)); pb <- unlist(lapply(pr$b, draw, k = n))
  sg <- unlist(lapply(u, function(tp) draw(tp, n, rep = length(idx[[tp]]) < n)))
  X <- cbind(M[G, pa, drop = FALSE] + M[G, pb, drop = FALSE], M[G, sg, drop = FALSE])
  list(X = t(as.matrix(X)), y = c(rep(paste(pr$a, pr$b, sep = "|"), each = n), rep(u, each = n))) }
nbs_train <- function(X, y, seed) {                 # stratified 80/20 split, xgboost in random-forest mode (Neighbor-seq settings)
  set.seed(seed)
  tr <- unlist(lapply(split(seq_along(y), y), function(i) i[sample.int(length(i), round(DBL_NBS_TRAIN_FRAC * length(i)))]))
  lv <- sort(unique(y)); p <- ncol(X)
  fit <- xgb.train(params = list(objective = "multi:softprob", eval_metric = "mlogloss", num_class = length(lv), max_depth = 20,
                                 num_parallel_tree = 200, subsample = 0.632, colsample_bytree = sqrt(p) / p,
                                 colsample_bynode = sqrt(p) / p, nthread = CORES, seed = seed),
                   data = xgb.DMatrix(X[tr, , drop = FALSE], label = match(y[tr], lv) - 1L), nrounds = 1)
  m <- list(fit = fit, lv = lv, G = colnames(X))
  hit <- nbs_predict(m, X[-tr, , drop = FALSE]) == y[-tr]
  c(m, list(acc = mean(hit), acc_TL = if (TL %in% lv) mean(hit[y[-tr] == TL]) else NA_real_)) }
nbs_predict <- function(m, Xq) {                    # most probable class of every row
  pr <- predict(m$fit, Xq[, m$G, drop = FALSE])
  if (is.null(dim(pr))) pr <- matrix(pr, ncol = length(m$lv), byrow = TRUE)   # older xgboost returns a flat vector
  m$lv[max.col(pr, ties.method = "first")] }

## -- [1] per lane: rebuild stage 04's held-out injections, train per fold, flag every variant per run tag ----
FLAGS <- list(); CAL <- list(); CLS <- list(); CONF <- list()
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
  I04 <- lapply(RUNSPEC$tag, function(TAG) { d <- fread(file.path(DBL_TAB_DIR, paste0("dbl_typing_injected_devel_", TAG, "_smoke.csv.gz")))[lane == L & kind == "heldout"]
    d[, `:=`(origin = nz(origin), true_state = nz(true_state))]; d }); names(I04) <- RUNSPEC$tag
  runs <- sort(unique(I04[[1]]$run))
  chk(sprintf("%s: every run tag holds the same held-out runs", L), TRUE, all(vapply(I04, function(d) identical(sort(unique(d$run)), runs), TRUE)))
  folds <- sort(unique(as.integer(sub("^f([0-9]+)_b.*", "\\1", runs))))
  ts("    runs: ", paste(runs, collapse = " "))
  ## marker genes per fold (training cells only)
  TRN <- list(); GEN <- list()
  for (f in folds) {
    trn <- A[fold != f][, lin := lineage_of(compartment)]
    trn <- trn[lin %in% trn[, .N, by = lin][N >= DBL_NBS_MIN_CELLS, lin]]
    chk(sprintf("%s fold %d: training cells disjoint from the held-out fold", L, f), 0L, length(intersect(trn$cell, A[fold == f, cell])))
    set.seed(DBL_SEED_NBS + f * 10L)
    GEN[[f]] <- nbs_markers(M, trn$cell, trn$lin); TRN[[f]] <- trn
    ts(sprintf("    fold %d: %d training cells, %d types (%s), %d marker genes", f, nrow(trn), uniqueN(trn$lin),
               paste(trn[, .N, by = lin][order(-N), sprintf("%s %d", lin, N)], collapse = ", "), length(GEN[[f]])))
  }
  ## rebuild every held-out run (verbatim construction as 09) + NNLS + LEUK share + T-identity gate; classifier features kept
  one <- function(rid) {
    f <- as.integer(sub("^f([0-9]+)_b.*", "\\1", rid)); b <- as.integer(sub(".*_b([0-9]+)$", "\\1", rid))
    set.seed(DBL_SEED_SYN + f * 10L + b)
    held <- A[fold == f, cell]
    inj <- build_injection(split(A[cell %in% held, cell], A[cell %in% held, compartment]), rid)
    Mi <- M[, inj$a, drop = FALSE]; dbl <- which(!is.na(inj$b)); Mi[, dbl] <- Mi[, dbl] + M[, inj$b[dbl], drop = FALSE]; colnames(Mi) <- inj$syn
    B <- make_basis(M, A[!(cell %in% held), .(cell, compartment)])
    inj[, nnls := nnls_lin(Mi[B$G, , drop = FALSE], B)]
    inj[, share_rb := leuk_share(Mi)]
    inj[, tid := Matrix::colSums(Mi[intersect(DBL_T_ID_GENES, rownames(Mi)), , drop = FALSE]) >= DBL_T_ID_MIN_UMI]
    list(dt = inj[, .(run = rid, fold = f, syn, type, true_state, nnls, share_rb, tid)], X = t(as.matrix(Mi[GEN[[f]], , drop = FALSE]))) }
  RB <- mclapply(runs, one, mc.cores = min(CORES, length(runs)))   # forks before any xgboost call (OpenMP is not fork-safe)
  stopifnot(!vapply(RB, inherits, TRUE, "try-error"))
  chk(sprintf("%s: every barcode given to the classifier is synthetic or re-entered (SYN_ ids)", L), TRUE,
      all(vapply(RB, function(r) all(startsWith(rownames(r$X), "SYN_")), TRUE)))
  ## train per fold, real and shuffled-label worlds, then predict every run
  MOD <- list()
  for (f in folds) for (w in c("real", "null")) {
    trn <- TRN[[f]]; sd <- DBL_SEED_NBS + f * 10L + (w == "null")
    set.seed(sd); ty <- if (w == "null") sample(trn$lin) else trn$lin
    AM <- nbs_artificial(M, GEN[[f]], trn$cell, ty, DBL_NBS_N_AM)
    MOD[[paste(f, w)]] <- nbs_train(AM$X, AM$y, sd)
    ts(sprintf("    fold %d %-4s: %d classes, own 20 %% test accuracy %.3f, own T+LEUK recall %.3f", f, w,
               length(MOD[[paste(f, w)]]$lv), MOD[[paste(f, w)]]$acc, MOD[[paste(f, w)]]$acc_TL)) }
  RBD <- rbindlist(lapply(RB, function(r) { f <- r$dt$fold[1]
    r$dt[, `:=`(nbs = nbs_predict(MOD[[paste(f, "real")]], r$X), nbs_null = nbs_predict(MOD[[paste(f, "null")]], r$X))] }))
  ## classifier diagnostics and confusion: synthetic doublets only (singlets are real cells)
  syn_t <- RBD[type %in% c("T|LEUK", "T|LEUK_other")]
  CLS[[L]] <- rbindlist(lapply(names(MOD), function(k) { f <- as.integer(sub(" .*", "", k)); w <- sub("^[0-9]+ ", "", k)
    data.table(lane = L, fold = f, world = w, n_genes = length(MOD[[k]]$G), n_classes = length(MOD[[k]]$lv), own_test_acc = MOD[[k]]$acc,
               own_test_TL_recall = MOD[[k]]$acc_TL,
               heldout_TL_recall = syn_t[fold == f, mean((if (w == "real") nbs else nbs_null) == TL)], n_heldout_TL = syn_t[fold == f, .N]) }))
  print(CLS[[L]][, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])
  for (r in which(CLS[[L]]$world == "null"))
    chk(sprintf("%s fold %d: shuffled-label classifier recalls <= %.2f of held-out T-LEUK", L, CLS[[L]]$fold[r], DBL_NBS_NULL_MAX), TRUE,
        CLS[[L]]$heldout_TL_recall[r] <= DBL_NBS_NULL_MAX)
  CONF[[L]] <- RBD[!startsWith(type, "singlet:"), .N, by = .(type, predicted = nbs)][, frac := N / sum(N), by = type][, lane := L]
  ## per run tag: merge with stage 04, identity checks, variant flags
  for (k in seq_len(nrow(RUNSPEC))) {
    TAG <- RUNSPEC$tag[k]; LB <- RUNSPEC$labels[k]; inj04 <- I04[[TAG]]
    m <- merge(inj04, RBD, by = c("run", "syn"), suffixes = c("", ".rb"))
    chk(sprintf("%s %s: rebuilt injection == stage 04's rows", TAG, L), nrow(inj04), nrow(m))
    chk(sprintf("%s %s: rebuilt injection == stage 04's types", TAG, L), TRUE, m[, all(type == type.rb)])
    chk(sprintf("%s %s: rebuilt injection == stage 04's LEUK share (1e-9, NA-aware)", TAG, L), TRUE,
        m[, all((is.na(share) & is.na(share_rb)) | (!is.na(share) & !is.na(share_rb) & abs(share - share_rb) < 1e-9))])
    chk(sprintf("%s %s: stage 04 carries the exact singlet no-window flag", TAG, L), TRUE, "tlin_nowin_any" %in% names(m))
    sg <- startsWith(m$type, "singlet:")
    m[, lin_o := lin_origin(origin)]
    m[, `:=`(V1a = called %in% TRUE & ambiguous %in% FALSE & lin_o, V3 = called %in% TRUE & nnls)]
    m[sg, V1a := tlin_nowin_any %in% TRUE]                              # singlet origins are masked in 04: its pooled-safe flag
    m[, V5 := V1a | V3]
    m[, `:=`(N1 = nbs == TL, N2 = called %in% TRUE & nbs == TL)]
    m[, `:=`(M1 = V5 & tid, N2M = N2 & tid)]
    FLAGS[[paste(TAG, L)]] <- m[, c("run", "syn", "type", "true_state", "share", "called", VARS), with = FALSE][, `:=`(tag = TAG, labels = LB, lane = L)]
    CAL[[paste(TAG, L)]] <- readRDS(file.path(DBL_LANE_DIR, paste0(L, "__calib_devel_", TAG, "_smoke.rds")))[c("dbr", "N", "comp")]
  }
  rm(x, M, RB, RBD); gc(FALSE)
}
FL <- rbindlist(FLAGS, fill = TRUE)
FL[, truth := fifelse(type == "T|LEUK", true_state, fifelse(type == "T|LEUK_other", "LEUK_other", NA_character_))]
## MASKED singlets: pooled per compartment right away; the per-barcode singlet rows are dropped
SGP <- melt(FL[startsWith(type, "singlet:")], id.vars = c("tag", "labels", "lane", "type"), measure.vars = VARS, variable.name = "variant", value.name = "flag")[
  , .(n = .N, k = sum(flag)), by = .(tag, labels, lane, variant, comp = sub("^singlet:", "", type))]
FL <- FL[!startsWith(type, "singlet:")]
## -- [2] copied verbatim from 09_recall_variants.R (lines 168-194): natural recall / precision ----
## -- [2] natural recall / precision per qualifying patient (stage 08's model) ----
RATES <- melt(FL, id.vars = c("tag", "labels", "lane", "type", "truth"), measure.vars = VARS, variable.name = "variant", value.name = "flag")[
  , .(rate = mean(flag)), by = .(tag, labels, lane, variant, type, truth)]
NAT <- list()
for (key in names(CAL)) { cb <- CAL[[key]]; TAG <- sub(" .*", "", key); L <- sub("^[^ ]+ ", "", key)
  comp <- cb$comp; N_tot <- cb$N; dbr <- cb$dbr
  qual <- comp[, .(nT = sum(N[compartment == "T"]), nL = sum(N[compartment %in% STATES])), by = donor][nT >= 30 & nL >= 50, donor]
  for (v in VARS) for (d in qual) {
    cd <- comp[donor == d][, f := N / sum(N)]; g <- function(cp) cd[compartment == cp, if (.N) f else 0]
    pd <- comp[, sum(N[donor == d]) / sum(N)]; Nd_dbl <- dbr * N_tot * pd^2; Nd_sing <- N_tot * pd - Nd_dbl
    R <- RATES[tag == TAG & lane == L & variant == v]
    rk <- function(kk) { r <- if (kk %in% STATES) R[type == "T|LEUK" & truth == kk, rate] else R[type == "T|LEUK_other", rate]; if (length(r)) r else NA_real_ }
    ru <- function(u) { r <- R[type == u, rate]; if (length(r)) r else 0 }
    fL <- sum(vapply(LIN, g, 0)); sig_k <- vapply(LIN, function(kk) 2 * g("T") * g(kk) * rk(kk), 0)
    fp_pairs <- 2 * g("T") * g("Tamb") * ru("T|Tamb") + 2 * g("T") * g("NK") * ru("T|NK") + 2 * g("T") * g("B") * ru("T|B") +
                2 * g("T") * g("Ery") * ru("T|Ery") + 2 * g("T") * g("MonoDC") * ru("T|MonoDC") + g("T")^2 * ru("T|T") +
                2 * g("Tamb") * fL * ru("Tamb|LEUK*") + 2 * g("NK") * fL * ru("NK|LEUK*") + 2 * g("B") * fL * ru("B|LEUK*") +
                2 * g("Ery") * fL * ru("Ery|LEUK*") + 2 * g("MonoDC") * fL * ru("MonoDC|LEUK*") + fL^2 * ru("LEUK*|LEUK*")
    sp <- SGP[tag == TAG & lane == L & variant == v][, r := k / n]
    fp_sing <- sum(vapply(cd$compartment, function(cp) { r <- sp[comp == cp, r]; if (length(r)) g(cp) * r else 0 }, 0))
    w <- vapply(LIN, g, 0); ok <- !is.na(sig_k)
    signal <- Nd_dbl * sum(sig_k, na.rm = TRUE); fp <- Nd_dbl * fp_pairs + Nd_sing * fp_sing
    NAT[[length(NAT) + 1]] <- data.table(tag = TAG, labels = RUNSPEC[tag == TAG, labels], lane = L, variant = v, donor = d,
      recall_nat = sum(w[ok] * vapply(LIN[ok], rk, 0)) / sum(w[ok]), precision_nat = signal / (signal + fp),
      exp_signal = signal, exp_fp = fp) } }
NAT <- rbindlist(NAT)
NAT[, score := pmin(recall_nat, precision_nat)]
## -- [3] validity, then selection (tune) or the one-shot confirmation ----
CLS <- rbindlist(CLS); CONF <- rbindlist(CONF)
fwrite(CLS,  file.path(DBL_TAB_DIR, paste0("dbl_nbs_classifier_", MODE, "_smoke.csv")))
fwrite(CONF, file.path(DBL_TAB_DIR, paste0("dbl_nbs_confusion_", MODE, "_smoke.csv")))
if (MODE == "tune") {                                                  # known answer: the copied pipeline reproduces stage 09
  t9 <- fread(file.path(DBL_TAB_DIR, "dbl_recall_tune_smoke.csv"))[variant == "V5"]
  v5 <- merge(NAT[variant == "V5"], t9, by = c("tag", "lane", "donor"), suffixes = c("", ".09"))
  chk("V5 reproduces stage 09: same patient rows", nrow(t9), nrow(v5))
  chk("V5 reproduces stage 09: natural recall and precision (1e-9)", TRUE,
      v5[, all(abs(recall_nat - recall_nat.09) < 1e-9 & abs(precision_nat - precision_nat.09) < 1e-9)])
}
VALID <- all(rbindlist(CHK)$pass)
CAND <- NAT[variant %in% NEWV, .(lanes = uniqueN(lane), patients = .N, worst = min(score), worst_recall = min(recall_nat),
                                 worst_precision = min(precision_nat)), by = .(labels, variant)]
CAND[, `:=`(vrank = match(variant, NEWV), lrank = match(labels, c("knn", "bmm")))]
CAND[is.na(worst), worst := -Inf]                                       # nothing flagged -> never chosen
setorder(CAND, -worst, vrank, lrank)
print(NAT[, .(lane = sub("_.*", "", lane), labels, variant, donor, recall = round(recall_nat, 2), precision = round(precision_nat, 2))][order(labels, variant, lane, donor)])
print(CAND[, .(labels, variant, patients, worst = round(worst, 3), worst_recall = round(worst_recall, 3), worst_precision = round(worst_precision, 3))])
if (!VALID) ts("** a validity check FAILED: the numbers above are not read (see the checks table) **")
if (MODE == "tune") {
  best <- CAND[lanes == length(DBL_SMOKE_LANES)][1]
  fwrite(NAT, file.path(DBL_TAB_DIR, "dbl_nbs_tune_smoke.csv"))
  fwrite(best[, .(labels, variant, worst, valid = VALID, reaches_bar = VALID & worst >= 0.7, chosen_on = "tuning folds 1-2")],
         file.path(DBL_TAB_DIR, "dbl_nbs_choice_smoke.csv"))
  ts(sprintf("[3] CHOSEN on the tuning folds: %s labels, %s (worst min(recall, precision) %.3f) -> %s", best$labels, best$variant, best$worst,
             if (!VALID) "INVALID run: nothing is read" else if (best$worst >= 0.7) "run fold 3 once with these labels"
             else "below 0.7 already on the tuning folds: the pair classifier does not reach the bar either; fold 3 is NOT run; stop"))
} else {
  ch <- fread(file.path(DBL_TAB_DIR, "dbl_nbs_choice_smoke.csv"))
  if (!isTRUE(ch$reaches_bar)) stop("the frozen choice did not reach 0.7 on the tuning folds: fold 3 must not be run")
  if (!all(RUNSPEC$labels == ch$labels)) stop("--runs labels differ from the frozen choice (", ch$labels, ")")
  res <- NAT[labels == ch$labels & variant == ch$variant]
  pass <- VALID && res[, uniqueN(lane) == length(DBL_SMOKE_LANES) && all(recall_nat >= 0.7 & precision_nat >= 0.7)]
  verdict <- if (!VALID) "INVALID run: not read" else if (pass) "lineage-level typing is good enough to design the O/E test" else "stop (no further tuning)"
  print(res[, .(lane = sub("_.*", "", lane), donor, recall = round(recall_nat, 3), precision = round(precision_nat, 3))])
  ts("[3] VERDICT on fold 3 (", ch$labels, " labels, ", ch$variant, "): ", verdict)
  fwrite(NAT, file.path(DBL_TAB_DIR, "dbl_nbs_confirm_smoke.csv"))
  fwrite(data.table(labels = ch$labels, variant = ch$variant, verdict = verdict), file.path(DBL_TAB_DIR, "dbl_nbs_verdict_smoke.csv"))
}
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, paste0("dbl_checks_10_", MODE, "_smoke.csv")))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
