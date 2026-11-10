#!/usr/bin/env Rscript
# 09_recall_variants.R ----
# Doublet-pairing SMOKE, stage 9 (user decision 2026-11-10: "try to raise recall"): lineage-level T-LEUK typing variants, chosen
# on TUNING folds 1-2 and confirmed ONCE on the never-used fold 3. LEUK (lineage) = five states + LEUK_other; detection = devel
# scDblFinder (stage 04 --scdbl=devel --heldout-only) in every variant.
# Pre-listed variants (no others are tried):
#   V0  baseline       scDblFinder origin {T, LEUK-lineage}, not ambiguous, ambient window 0.25-0.80 (stage 08's scdbl_devel)
#   V1a no window      as V0 without the window                                                        (recall lever 1)
#   V1b lane window    as V0 with a per-lane window = 5th-95th percentile of the LEUK share of the tuning folds' synthetic
#                      T-LEUK-lineage doublets that V1a reports (frozen at tuning, reused unchanged at confirmation)
#   V3  NNLS merged    called by devel; NNLS against compartment profiles of non-held-out class-A cells, LEUK-lineage
#                      components summed into one group; the two heaviest groups are T and LEUK, the minor >= 15 %, and
#                      the second group exceeds the third by >= 0.05                                    (recall lever 3)
#   V4  V1a AND V3     both agree (precision-leaning)
#   V5  V1a OR V3      either (recall-leaning)
#   each with B/C labels from stage 02 (knn) or from 02b (bmm projection)                             (recall lever 2)
# Natural recall / precision per qualifying patient: stage 08's random co-encapsulation model (composition, dbr, N from the
# run's calibration; unsimulated partners borrow the nearest simulated type).
# SELECTION (tuning, fixed 2026-11-10 before any tuning result was read): the (variant, labels) with the largest
#   min over qualifying patients of min(natural recall, natural precision); ties -> fewer components (V0 < V1a < V1b < V3 <
#   V4 < V5; knn before bmm). Frozen to dbl_recall_choice_smoke.csv. If no candidate reaches 0.7 on the tuning folds, fold 3
#   is NOT run: verdict "recall could not be raised enough (already on the tuning folds); stop".
# READING RULE (confirmation): the chosen candidate passes if, on fold 3, natural recall >= 0.7 and natural precision >= 0.7 for
#   every qualifying patient in both lanes -> "recall raised: lineage-level typing is good enough to design the O/E test";
#   otherwise "recall could not be raised enough; stop (no further tuning)". Only this result counts.
# MASKING: re-entered singlets are real cells -- for them only pooled per-compartment counts of a variant's lineage flag leave
#   this script. Their NNLS call is computed here and immediately pooled. Singlet flags need stage 04's tlin_nowin_any: the
#   first tuning runs (jobs 3687425/3687426) predate it, so their V1a/V1b/V4/V5 singlet false positives were lower bounds
#   (review wf_7d4cf7f8-170 F2); the tuning was rerun with exact flags (tags tune2_knn / tune2_bmm, same seeds).
# INPUT  : DBL_TAB_DIR/dbl_typing_injected_devel_<tag>_smoke.csv.gz, DBL_LANE_DIR/<lane>__calib_devel_<tag>_smoke.rds (04),
#          DBL_LANE_DIR/<lane>.rds, <lane>__labels{,_bmm}.csv.gz
# OUTPUT : tune:    DBL_TAB_DIR/dbl_recall_tune_smoke.csv, dbl_recall_choice_smoke.csv
#          confirm: DBL_TAB_DIR/dbl_recall_confirm_smoke.csv, dbl_recall_verdict_smoke.csv;  both: dbl_checks_09_<mode>_smoke.csv
# Usage  : Rscript scripts/dbl/09_recall_variants.R --mode=tune --runs=tune_knn:knn,tune_bmm:bmm
#          Rscript scripts/dbl/09_recall_variants.R --mode=confirm --runs=f3_<labels>:<labels>

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here); library(parallel); library(nnls) })
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
VARS <- c("V0", "V1a", "V1b", "V3", "V4", "V5")
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

## -- [1] per (run tag, lane): rebuild, check identity, compute variant flags ----
FLAGS <- list(); SGL <- list(); CAL <- list()
for (k in seq_len(nrow(RUNSPEC))) for (L in DBL_SMOKE_LANES) {
  TAG <- RUNSPEC$tag[k]; LB <- RUNSPEC$labels[k]
  ts(sprintf("[1] %s (%s labels), lane %s", TAG, LB, L))
  inj04 <- fread(file.path(DBL_TAB_DIR, paste0("dbl_typing_injected_devel_", TAG, "_smoke.csv.gz")))[lane == L]
  inj04[, `:=`(origin = nz(origin), true_state = nz(true_state))]
  x   <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds"))); M <- x$counts
  lab <- fread(file.path(DBL_LANE_DIR, paste0(L, if (LB == "bmm") "__labels_bmm.csv.gz" else "__labels.csv.gz")))[match(colnames(M), cell)]
  A <- lab[class == "A"]
  pools <- split(A$cell, A$compartment)
  cpm <- function(cells) { s <- M[, cells, drop = FALSE]; Matrix::rowSums(s) / sum(s) * 1e6 }
  badg <- grepl("^MT-|^RP[SL]|^HB[^P]", rownames(M))
  cT <- cpm(pools[["T"]]); cL <- cpm(unlist(pools[intersect(STATES, names(pools))]))
  GT4 <- rownames(M)[!badg & cT > 5 & cL < 1]; GL4 <- rownames(M)[!badg & cL > 5 & cT < 1]
  leuk_share <- function(Mx) { l <- Matrix::colSums(Mx[GL4, , drop = FALSE]); t <- Matrix::colSums(Mx[GT4, , drop = FALSE]); ifelse(l + t > 0, l / (l + t), NA_real_) }
  set.seed(DBL_SEED_SYN)
  A[, fold := ((sample(.N) - 1L) %% DBL_N_FOLDS) + 1L, by = compartment]
  runs <- unique(inj04[kind == "heldout", run])
  one <- function(rid) {
    f <- as.integer(sub("^f([0-9]+)_b.*", "\\1", rid)); b <- as.integer(sub(".*_b([0-9]+)$", "\\1", rid))
    set.seed(DBL_SEED_SYN + f * 10L + b)
    held <- A[fold == f, cell]
    inj <- build_injection(split(A[cell %in% held, cell], A[cell %in% held, compartment]), rid)
    Mi <- M[, inj$a, drop = FALSE]; dbl <- which(!is.na(inj$b)); Mi[, dbl] <- Mi[, dbl] + M[, inj$b[dbl], drop = FALSE]; colnames(Mi) <- inj$syn
    B <- make_basis(M, A[!(cell %in% held), .(cell, compartment)])
    inj[, nnls := nnls_lin(Mi[B$G, , drop = FALSE], B)]
    inj[, share_rb := leuk_share(Mi)]
    inj[, .(run = rid, syn, type, true_state, nnls, share_rb)] }
  RB <- rbindlist(mclapply(runs, one, mc.cores = min(CORES, length(runs))))
  m <- merge(inj04[kind == "heldout"], RB, by = c("run", "syn"), suffixes = c("", ".rb"))
  chk(sprintf("%s %s: rebuilt injection == stage 04's rows", TAG, L), inj04[kind == "heldout", .N], nrow(m))
  chk(sprintf("%s %s: rebuilt injection == stage 04's types", TAG, L), TRUE, m[, all(type == type.rb)])
  chk(sprintf("%s %s: rebuilt injection == stage 04's LEUK share (1e-9, NA-aware)", TAG, L), TRUE,
      m[, all((is.na(share) & is.na(share_rb)) | (!is.na(share) & !is.na(share_rb) & abs(share - share_rb) < 1e-9))])
  sg <- startsWith(m$type, "singlet:")
  m[, lin_o := lin_origin(origin)]
  m[, `:=`(V0 = called %in% TRUE & ambiguous %in% FALSE & lin_o & !is.na(share) & share >= WIN[1] & share <= WIN[2],
           V1a = called %in% TRUE & ambiguous %in% FALSE & lin_o, V3 = called %in% TRUE & nnls)]
  ## singlets: origin is masked in stage 04 -> pooled flags only (V1a lower bound at tuning when the no-window flag is absent)
  has_nw <- "tlin_nowin_any" %in% names(m)
  m[sg, `:=`(V0 = tlin_any %in% TRUE, V1a = if (has_nw) tlin_nowin_any %in% TRUE else tlin_any %in% TRUE)]
  FLAGS[[paste(TAG, L)]] <- m[, c("run", "syn", "type", "true_state", "share", "called", VARS[VARS %in% names(m)]), with = FALSE][
    , `:=`(tag = TAG, labels = LB, lane = L, exact_singlets = has_nw)]
  CAL[[paste(TAG, L)]] <- readRDS(file.path(DBL_LANE_DIR, paste0(L, "__calib_devel_", TAG, "_smoke.rds")))[c("dbr", "N", "comp")]
  rm(x, M); gc(FALSE)
}
FL <- rbindlist(FLAGS, fill = TRUE)
FL[, truth := fifelse(type == "T|LEUK", true_state, fifelse(type == "T|LEUK_other", "LEUK_other", NA_character_))]
## V1b window: frozen at tuning per lane from V1a-reported synthetic T-LEUK-lineage doublets of the tuning folds
WINF <- file.path(DBL_TAB_DIR, "dbl_recall_window_smoke.csv")
if (MODE == "tune") {
  WL <- FL[!startsWith(type, "singlet:") & !is.na(truth) & V1a == TRUE, .(lo = quantile(share, 0.05, na.rm = TRUE), hi = quantile(share, 0.95, na.rm = TRUE)), by = .(labels, lane)]
  fwrite(WL, WINF); print(WL)
} else WL <- fread(WINF)
FL <- merge(FL, WL, by = c("labels", "lane"))
FL[!startsWith(type, "singlet:"), V1b := V1a & !is.na(share) & share >= lo & share <= hi]
FL[startsWith(type, "singlet:"), V1b := V1a & !is.na(share) & share >= lo & share <= hi]   # frozen per-lane window (review F1)
FL[, `:=`(V4 = V1a & V3, V5 = V1a | V3)]
## MASKED singlets: pooled per compartment right away; the per-barcode singlet rows are dropped
SGP <- melt(FL[startsWith(type, "singlet:")], id.vars = c("tag", "labels", "lane", "type"), measure.vars = VARS, variable.name = "variant", value.name = "flag")[
  , .(n = .N, k = sum(flag)), by = .(tag, labels, lane, variant, comp = sub("^singlet:", "", type))]
FL <- FL[!startsWith(type, "singlet:")]
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
CAND <- NAT[, .(lanes = uniqueN(lane), patients = .N, worst = min(score), worst_recall = min(recall_nat), worst_precision = min(precision_nat)), by = .(labels, variant)]
CAND[, `:=`(vrank = match(variant, VARS), lrank = match(labels, c("knn", "bmm")))]
CAND[is.na(worst), worst := -Inf]                                       # nothing flagged -> never chosen
setorder(CAND, -worst, vrank, lrank)
print(NAT[, .(lane = sub("_.*", "", lane), labels, variant, donor, recall = round(recall_nat, 2), precision = round(precision_nat, 2))][order(labels, variant, lane, donor)])
print(CAND[, .(labels, variant, patients, worst = round(worst, 3), worst_recall = round(worst_recall, 3), worst_precision = round(worst_precision, 3))])
if (MODE == "tune") {
  best <- CAND[lanes == length(DBL_SMOKE_LANES)][1]
  fwrite(NAT, file.path(DBL_TAB_DIR, "dbl_recall_tune_smoke.csv"))
  fwrite(best[, .(labels, variant, worst, reaches_bar = worst >= 0.7, chosen_on = "tuning folds 1-2")], file.path(DBL_TAB_DIR, "dbl_recall_choice_smoke.csv"))
  ts(sprintf("[3] CHOSEN on the tuning folds: %s labels, %s (worst min(recall, precision) %.3f) -> %s", best$labels, best$variant, best$worst,
             if (best$worst >= 0.7) "run fold 3 once with these labels" else "below 0.7 already on the tuning folds: fold 3 is NOT run; stop"))
} else {
  ch <- fread(file.path(DBL_TAB_DIR, "dbl_recall_choice_smoke.csv"))
  if (!isTRUE(ch$reaches_bar)) stop("the frozen choice did not reach 0.7 on the tuning folds: fold 3 must not be run")
  if (!all(RUNSPEC$labels == ch$labels)) stop("--runs labels differ from the frozen choice (", ch$labels, ")")
  res <- NAT[labels == ch$labels & variant == ch$variant]
  pass <- res[, uniqueN(lane) == length(DBL_SMOKE_LANES) && all(recall_nat >= 0.7 & precision_nat >= 0.7)]
  verdict <- if (pass) "recall raised: lineage-level typing is good enough to design the O/E test" else "recall could not be raised enough; stop (no further tuning)"
  print(res[, .(lane = sub("_.*", "", lane), donor, recall = round(recall_nat, 3), precision = round(precision_nat, 3))])
  ts("[3] VERDICT on fold 3 (", ch$labels, " labels, ", ch$variant, "): ", verdict)
  fwrite(NAT, file.path(DBL_TAB_DIR, "dbl_recall_confirm_smoke.csv"))
  fwrite(data.table(labels = ch$labels, variant = ch$variant, verdict = verdict), file.path(DBL_TAB_DIR, "dbl_recall_verdict_smoke.csv"))
}
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, paste0("dbl_checks_09_", MODE, "_smoke.csv")))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
