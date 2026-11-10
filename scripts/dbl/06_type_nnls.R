#!/usr/bin/env Rscript
# 06_type_nnls.R ----
# Doublet-pairing SMOKE, stage 6 (probe; user decision 2026-11-10): can a DECOMPOSITION type T-LEUK doublets by leukemic state
# where scDblFinder's origin vote could not (stage 04)? Each synthetic doublet is decomposed against compartment profiles
# built from class-A cells that were NOT held out:
#   nnls   (primary)   weighted NNLS on gene proportions; the pair = the two heaviest components (proposal 4.2 / appendix D)
#   mnlik  (cross-check) two-component multinomial likelihood over every compartment pair, mixing fraction optimised
# Doublets: (1) exactly the held-out and leave-one-donor-out synthetic doublets of stage 04 (rebuilt by replaying its seeds and
# checked against its injected table), so both methods are scored on identical barcodes; (2) pooled lane only, new within-donor
# doublets of the donors that pass the singlet floors, typed against that donor's own profiles (the proposal's patient-specific
# references). Scopes: 'all' = every synthetic doublet (pure typing); 'called_main' / 'called_devel' = only the doublets that
# scDblFinder called in stage 04 (detect, then type).
# MASKING [D2 + decision 2026-11-10]: no real barcode is decomposed. Re-entered held-out singlets are real cells: for them only
#   "decomposed as any T-LEUK" is kept. Stage 04's injected tables are read WITHOUT their singlet rows.
# READING RULE (fixed 2026-11-10 before any run; not changed after results):
#   primary = held-out runs, lane-level profiles, nnls, scope called_devel; recall = correctly typed / all injected T-LEUK(k)
#   (detection loss included); precision = correctly typed / all synthetic doublets typed as that state (pair types included,
#   singlets excluded). A merge level (K5..K2, appendix-B order) passes in a lane if every state at that level with >= 30
#   injected doublets has recall >= 0.7 and precision >= 0.7. NNLS is "worth writing into the proposal" only if, in BOTH lanes,
#   some level with >= 2 states passes; otherwise "state typing does not work in these two lanes". Everything else is secondary.
# INPUT  : DBL_LANE_DIR/<lane>.rds (01), <lane>__labels.csv.gz (02), DBL_TAB_DIR/dbl_typing_injected_{main,devel}_smoke.csv.gz (04)
# OUTPUT : DBL_TAB_DIR/dbl_nnls_summary_smoke.csv   lane x kind x refs x scope x method x level x state: n, correct, recall, precision
#          DBL_TAB_DIR/dbl_nnls_levels_smoke.csv    lane x kind x refs x scope x method x level: states evaluated, pass (the decision)
#          DBL_TAB_DIR/dbl_nnls_fp_smoke.csv        false-positive sources: n, decomposed as any T-LEUK (singlets pooled only)
#          DBL_TAB_DIR/dbl_nnls_injected_smoke.csv.gz  synthetic barcodes: pair, weights, minor fraction (singlet rows masked)
#          DBL_TAB_DIR/dbl_nnls_genes_smoke.csv     genes used per run;  DBL_TAB_DIR/dbl_checks_06_smoke.csv
# Usage  : Rscript scripts/dbl/06_type_nnls.R [--cores=8]

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here); library(parallel); library(nnls) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
ARGS  <- commandArgs(TRUE)
argv  <- function(key, default) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), ARGS, value = TRUE)); if (length(v)) v else default }
CORES <- as.integer(argv("cores", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
ts  <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-82s expect %-10s got %-10s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
STATES <- DBL_LEUK_STATES
NNLS_TOP_ONE   <- 50L     # [SMOKE] genes per compartment, one-versus-rest log fold change of pseudobulk CPM (CPM >= 10)
NNLS_TOP_PAIR  <- 20L     # [SMOKE] genes per pair of LEUK states, |log fold change| (state-versus-state discrimination)
NNLS_MIN_REF   <- 5L      # [SMOKE] a compartment enters the basis with >= 5 reference cells
NNLS_MIN_MINOR <- 0.15    # [SMOKE] the second component must carry >= 15 % of the weight (else: a singlet plus noise)
NNLS_MARGIN    <- 0.05    # [SMOKE] nnls: second minus third weight >= 0.05; mnlik: best minus second-best pair log-lik >= 2
MNLIK_MARGIN   <- 2
LEVELS <- list(K5 = setNames(STATES, STATES),
               K4 = setNames(c("HSC+MPP", "HSC+MPP", "LMPP", "GMP", "ProMono"), STATES),
               K3 = setNames(c("prim", "prim", "prim", "GMP", "ProMono"), STATES),
               K2 = setNames(c("prim+GMP", "prim+GMP", "prim+GMP", "prim+GMP", "ProMono"), STATES))
bad_gene <- function(g) grepl("^MT-|^RP[SL]|^HB[^P]", g)

## -- helpers: the decomposition ----
make_basis <- function(M, ref) {   # ref: data.table(cell, compartment) -> list(G = genes, P = genes x compartments (proportions))
  keep <- ref[, .N, by = compartment][N >= NNLS_MIN_REF, compartment]
  ref <- ref[compartment %in% keep]
  ok <- !bad_gene(rownames(M))
  S <- sapply(keep, function(cp) Matrix::rowSums(M[ok, ref[compartment == cp, cell], drop = FALSE]))   # genes x comps, raw sums
  CP <- sweep(S, 2, colSums(S), "/") * 1e6
  top_finite <- function(v, n) head(names(sort(v[is.finite(v)], decreasing = TRUE)), n)   # genes failing the CPM floor are -Inf
  G <- character(0)
  for (cp in keep) { o <- setdiff(keep, cp); lfc <- log2((CP[, cp] + 1) / (rowMeans(CP[, o, drop = FALSE]) + 1)); lfc[CP[, cp] < 10] <- -Inf
    G <- union(G, top_finite(lfc, NNLS_TOP_ONE)) }
  ls <- intersect(STATES, keep)
  if (length(ls) >= 2) for (i in seq_along(ls)) for (j in seq_along(ls)) if (i < j) {
    d <- abs(log2((CP[, ls[i]] + 1) / (CP[, ls[j]] + 1))); d[pmax(CP[, ls[i]], CP[, ls[j]]) < 10] <- -Inf
    G <- union(G, top_finite(d, NNLS_TOP_PAIR)) }
  P <- S[G, , drop = FALSE] + 0.5; P <- sweep(P, 2, colSums(P), "/")   # pseudo-count 0.5 keeps the likelihood finite
  list(G = G, P = P) }
type_nnls <- function(X, B) {   # X: genes(B$G) x barcodes counts -> one row per barcode
  Y <- sweep(as.matrix(X), 2, pmax(colSums(as.matrix(X)), 1), "/")
  w <- 1 / sqrt(rowMeans(B$P) + 1e-6)                              # approximately Poisson: variance ~ mean
  Pw <- B$P * w
  rbindlist(lapply(seq_len(ncol(Y)), function(j) {
    b <- nnls::nnls(Pw, Y[, j] * w)$x; b <- b / max(sum(b), 1e-12); o <- order(b, decreasing = TRUE)
    data.table(c1 = colnames(B$P)[o[1]], c2 = colnames(B$P)[o[2]], w1 = b[o[1]], w2 = b[o[2]], w3 = if (length(o) > 2) b[o[3]] else 0) })) }
type_mnlik <- function(X, B) {   # best two-component multinomial mixture over all compartment pairs
  X <- as.matrix(X); lp <- colnames(B$P); pr <- combn(lp, 2, simplify = FALSE)
  rbindlist(lapply(seq_len(ncol(X)), function(j) {
    x <- X[, j]
    fits <- vapply(pr, function(p) { p1 <- B$P[, p[1]]; p2 <- B$P[, p[2]]
      o <- optimize(function(a) sum(x * log(a * p1 + (1 - a) * p2)), c(0.001, 0.999), maximum = TRUE); c(o$objective, o$maximum) }, numeric(2))
    o <- order(fits[1, ], decreasing = TRUE); a <- fits[2, o[1]]
    data.table(c1 = pr[[o[1]]][1], c2 = pr[[o[1]]][2], w1 = a, w2 = 1 - a, margin = fits[1, o[1]] - fits[1, o[2]]) })) }
call_state <- function(d, method) {   # -> reported state if the pair is exactly {T, one LEUK state} and passes the margins
  st <- ifelse(d$c1 == "T" & d$c2 %in% STATES, d$c2, ifelse(d$c2 == "T" & d$c1 %in% STATES, d$c1, NA_character_))
  minor <- pmin(d$w1, d$w2)
  ok <- !is.na(st) & minor >= NNLS_MIN_MINOR & (if (method == "nnls") (d$w2 - d$w3) >= NNLS_MARGIN else d$margin >= MNLIK_MARGIN)
  ifelse(ok, st, "none") }

## -- the synthetic doublets of stage 04, rebuilt (code copied verbatim from 04_type_doublets.R, same seeds) ----
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
build_injection <- function(src, run_id, with_singlets = TRUE) {   # identical to stage 04's held-out construction
  prs <- rbindlist(c(
    lapply(STATES, function(s) build_pairs(src, DBL_SYN_PER_STATE, "T|LEUK", "T", s, s)),
    lapply(c("Tamb", "NK", "B", "Ery", "MonoDC", "LEUK_other", "T"), function(cp) build_pairs(src, DBL_SYN_PER_OTHER, paste0("T|", cp), "T", cp)),
    lapply(c("Tamb", "NK", "B", "Ery", "MonoDC"), function(cp) build_leuk_pairs(src, DBL_SYN_PER_OTHER, paste0(cp, "|LEUK*"), cp)),
    list(build_leuk_pairs(src[intersect(names(src), c(STATES, "LEUK_other"))] |> (\(s) { s[["LEUKx"]] <- unlist(s[intersect(STATES, names(s))]); s })(),
                          DBL_SYN_PER_OTHER, "LEUK*|LEUK*", "LEUKx"))), fill = TRUE)
  inj <- prs
  if (with_singlets) {
    comps <- names(src)[lengths(src) > 0]
    per   <- max(1L, floor(DBL_SYN_SINGLETS / length(comps)))
    sng   <- rbindlist(lapply(comps, function(cp) data.table(type = paste0("singlet:", cp), true_state = NA_character_,
                               a = if (length(src[[cp]]) >= per) sample(src[[cp]], per) else src[[cp]], b = NA_character_, first = cp)))
    inj <- rbind(prs, sng, fill = TRUE) }
  inj[, syn := sprintf("SYN_%s_%04d", run_id, .I)]
  inj }
counts_of <- function(M, inj) {   # synthetic doublet = summed counts of its two parents; a singlet row = its own counts
  Mi <- M[, inj$a, drop = FALSE]; dbl <- which(!is.na(inj$b)); Mi[, dbl] <- Mi[, dbl] + M[, inj$b[dbl], drop = FALSE]
  colnames(Mi) <- inj$syn; Mi }

## -- known answer on the functions: noise-free profile mixtures ----
## (run per lane below, on the lane's own basis)

INJ04 <- rbindlist(lapply(c("main", "devel"), function(m) {
  x <- fread(file.path(DBL_TAB_DIR, paste0("dbl_typing_injected_", m, "_smoke.csv.gz")))
  x <- x[!startsWith(type, "singlet:")]                  # masked: never load stage 04's singlet rows past this line
  x[true_state == "", true_state := NA_character_]       # fwrite wrote NA as "": without this the identity check returns NA
  x[, .(lane, run, syn, type, true_state, share, called, mode = m)] }))   # (job 3687415 hit it; NA-aware recheck: identical)
RES <- list(); FPS <- list(); GEN <- list(); INJ <- list()
for (L in DBL_SMOKE_LANES) {
  ts("[1] lane ", L)
  x   <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds"))); M <- x$counts
  lab <- fread(file.path(DBL_LANE_DIR, paste0(L, "__labels.csv.gz")))[match(colnames(M), cell)]
  stopifnot(identical(lab$cell, colnames(M)))
  A <- lab[class == "A"]
  ## stage 04's ambient marker sets, needed only for the identity check on 'share'
  pools <- split(A$cell, A$compartment)
  cpm <- function(cells) { s <- M[, cells, drop = FALSE]; Matrix::rowSums(s) / sum(s) * 1e6 }
  bad <- grepl("^MT-|^RP[SL]|^HB[^P]", rownames(M))
  cT  <- cpm(pools[["T"]]); cL <- cpm(unlist(pools[intersect(STATES, names(pools))]))
  GT  <- rownames(M)[!bad & cT > 5 & cL < 1]; GL <- rownames(M)[!bad & cL > 5 & cT < 1]
  leuk_share <- function(Mx) { l <- Matrix::colSums(Mx[GL, , drop = FALSE]); t <- Matrix::colSums(Mx[GT, , drop = FALSE]); ifelse(l + t > 0, l / (l + t), NA_real_) }
  set.seed(DBL_SEED_SYN)
  A[, fold := ((sample(.N) - 1L) %% DBL_N_FOLDS) + 1L, by = compartment]   # stage 04's folds
  donors <- sort(unique(A$donor))
  ## run list: stage-04 held-out (+ lodo) runs, then within-donor runs (pooled lane)
  RUNS <- list()
  for (f in DBL_FOLDS_RUN) for (b in seq_len(DBL_SYN_BATCHES)) RUNS[[length(RUNS) + 1]] <- list(kind = "heldout", id = sprintf("f%d_b%d", f, b), fold = f, batch = b, refs = "lane")
  if (length(donors) > 1) {
    lodo <- A[, .(nT = sum(compartment == "T"), nL = sum(compartment %in% STATES)), by = donor][order(-pmin(nT, nL))][1, donor]
    RUNS[[length(RUNS) + 1]] <- list(kind = "lodo", id = paste0("lodo_", lodo), donor = lodo, refs = "lane")
    floors <- A[, .(nT = sum(compartment == "T"), nL = sum(compartment %in% STATES)), by = donor][nT >= 30 & nL >= 50, donor]
    for (di in seq_along(floors)) for (f in DBL_FOLDS_RUN)
      RUNS[[length(RUNS) + 1]] <- list(kind = "perpatient", id = sprintf("pp_%s_f%d", floors[di], f), donor = floors[di], fold = f, di = di, refs = "patient") }
  one_run <- function(R) {
    if (R$kind == "perpatient") {
      set.seed(DBL_SEED_SYN + 5000L + R$di * 100L + R$fold * 10L)
      held <- A[donor == R$donor & fold == R$fold, cell]; refc <- A[donor == R$donor & !(cell %in% held)]
      inj  <- build_injection(split(A[cell %in% held, cell], A[cell %in% held, compartment]), R$id, with_singlets = FALSE)
    } else {
      set.seed(DBL_SEED_SYN + if (R$kind == "lodo") 999L else R$fold * 10L + R$batch)
      held <- if (R$kind == "lodo") A[donor == R$donor, cell] else A[fold == R$fold, cell]
      refc <- A[!(cell %in% held)]
      inj  <- build_injection(split(A[cell %in% held, cell], A[cell %in% held, compartment]), R$id, with_singlets = TRUE) }
    Mi <- counts_of(M, inj)
    inj[, share := leuk_share(Mi)]
    B  <- make_basis(M, refc[, .(cell, compartment)])
    X  <- Mi[B$G, , drop = FALSE]
    ## known answer: noise-free mixtures 0.35 T + 0.65 state are typed as that state (both methods)
    ka <- intersect(STATES, colnames(B$P))
    if ("T" %in% colnames(B$P) && length(ka)) {
      Mix <- sapply(ka, function(k) round(3000 * (0.35 * B$P[, "T"] + 0.65 * B$P[, k])))
      ka_n <- mean(call_state(type_nnls(Mix, B), "nnls") == ka); ka_m <- mean(call_state(type_mnlik(Mix, B), "mnlik") == ka)
    } else { ka_n <- NA_real_; ka_m <- NA_real_ }
    ## negative test: shuffle the LEUK state labels among the LEUK reference cells -> state typing should collapse
    rp <- copy(refc[, .(cell, compartment)]); il <- which(rp$compartment %in% STATES)
    set.seed(DBL_SEED_SYN + 7000L); rp[il, compartment := sample(compartment)]
    Bp <- make_basis(M, rp)
    tl <- which(inj$type == "T|LEUK")
    neg <- if (length(tl)) mean(call_state(type_nnls(Mi[Bp$G, tl, drop = FALSE], Bp), "nnls") == inj$true_state[tl]) else NA_real_
    dn <- type_nnls(X, B); dm <- type_mnlik(X, B)
    out <- rbind(cbind(inj[, .(syn, type, true_state, share)], dn[, .(c1, c2, w1, w2, w3, margin = NA_real_)], method = "nnls"),
                 cbind(inj[, .(syn, type, true_state, share)], dm[, .(c1, c2, w1, w2, w3 = NA_real_, margin)], method = "mnlik"))
    out[, reported := call_state(.SD, .BY$method), by = method]      # .SD = this method's rows only
    out[, tleuk_any := reported != "none"]
    sg <- startsWith(out$type, "singlet:")                          # MASKED: real held-out cells keep only tleuk_any
    out[sg, `:=`(c1 = NA_character_, c2 = NA_character_, w1 = NA_real_, w2 = NA_real_, w3 = NA_real_, margin = NA_real_, reported = NA_character_)]
    out[, `:=`(run = R$id, kind = R$kind, refs = R$refs)]
    obs <- out[method == "nnls" & type == "T|LEUK", mean(reported == true_state)]   # same scope as neg: all T|LEUK, K5
    list(out = out, genes = data.table(run = R$id, gene = B$G), basis = colnames(B$P),
         check = data.table(run = R$id, kind = R$kind, ka_nnls = ka_n, ka_mnlik = ka_m, obs_recall = obs, neg_recall = neg,
                            n_genes = length(B$G), basis = paste(colnames(B$P), collapse = ","))) }
  ts(sprintf("[2] %d runs on %d cores: %s", length(RUNS), min(CORES, length(RUNS)), paste(vapply(RUNS, `[[`, "", "id"), collapse = ", ")))
  RR <- mclapply(RUNS, function(R) tryCatch(one_run(R), error = function(e) list(error = conditionMessage(e), id = R$id)),
                 mc.cores = min(CORES, length(RUNS)), mc.preschedule = FALSE)
  RR <- Map(function(r, R) if (is.list(r) && !inherits(r, "try-error")) r else list(error = paste("worker died:", as.character(r)), id = R$id), RR, RUNS)
  for (r in Filter(function(r) !is.null(r$error), RR)) ts("    ** run ", r$id, " failed: ", r$error)
  chk(sprintf("%s: every decomposition run finished", L), 0L, length(Filter(function(r) !is.null(r$error), RR)))
  RR <- Filter(function(r) is.null(r$error), RR)
  if (!length(RR)) next
  OUT <- rbindlist(lapply(RR, `[[`, "out"))[, lane := L]
  CK  <- rbindlist(lapply(RR, `[[`, "check")); print(CK)
  chk(sprintf("%s: known answer -- noise-free mixtures typed correctly (nnls, every run)", L), TRUE, CK[, all(ka_nnls == 1, na.rm = TRUE) && any(!is.na(ka_nnls))])
  chk(sprintf("%s: known answer -- noise-free mixtures typed correctly (mnlik, every run)", L), TRUE, CK[, all(ka_mnlik == 1, na.rm = TRUE) && any(!is.na(ka_mnlik))])
  ng <- CK[kind == "heldout" & obs_recall >= 0.3]
  if (nrow(ng)) chk(sprintf("%s: negative test -- shuffled state profiles lose >= 0.1 recall in every informative held-out run", L), TRUE,
                    ng[, all(neg_recall <= obs_recall - 0.1)])
  else ts("    negative test not informative: no held-out run reaches recall 0.3 (scope all, K5)")
  ## identity with stage 04: the same synthetic doublets (type, true state, LEUK share), non-singlet rows
  idc <- merge(unique(OUT[kind != "perpatient" & !startsWith(type, "singlet:"), .(run, syn, type, true_state, share)]),
               unique(INJ04[lane == L, .(run, syn, type4 = type, true_state4 = true_state, share4 = share)]), by = c("run", "syn"))
  chk(sprintf("%s: rebuilt doublets == stage 04's (same rows)", L), INJ04[lane == L & mode == "main", .N], nrow(idc))
  chk(sprintf("%s: rebuilt doublets == stage 04's (type, true state, share to 1e-9)", L), TRUE,
      idc[, all(type == type4 & (true_state == true_state4 | (is.na(true_state) & is.na(true_state4))) & (abs(share - share4) < 1e-9 | (is.na(share) & is.na(share4))))])
  ## scopes
  called <- dcast(INJ04[lane == L], run + syn ~ mode, value.var = "called")
  OUT <- merge(OUT, called, by = c("run", "syn"), all.x = TRUE)
  setnames(OUT, c("main", "devel"), c("called_main", "called_devel"), skip_absent = TRUE)
  INJ[[L]] <- OUT
  GEN[[L]] <- rbindlist(lapply(RR, `[[`, "genes"))[, lane := L]
  ## metrics per scope x method x level
  for (sc in c("all", "called_main", "called_devel")) for (mt in c("nnls", "mnlik")) for (kd in unique(OUT$kind)) for (lv in names(LEVELS)) {
    if (sc != "all" && kd == "perpatient") next                     # no scDblFinder call for the new within-donor doublets
    D <- OUT[method == mt & kind == kd & !startsWith(type, "singlet:")]
    if (!nrow(D)) next
    hit <- if (sc == "all") rep(TRUE, nrow(D)) else (D[[sc]] %in% TRUE)   # detection: an uncalled doublet is never typed
    mp <- LEVELS[[lv]]
    D[, `:=`(tm = unname(mp[true_state]), rm = fifelse(hit & reported != "none", unname(mp[reported]), "none"))]
    tl <- D[type == "T|LEUK"]
    r <- tl[, .(n = .N, correct = sum(rm == tm), retained = sum(rm != "none")), by = .(merged = tm)]
    p <- D[rm != "none", .(typed = .N, typed_from_tleuk = sum(type == "T|LEUK"), typed_correct = sum(type == "T|LEUK" & rm == tm, na.rm = TRUE)), by = .(merged = rm)]
    s <- merge(r, p, by = "merged", all = TRUE)
    s[, `:=`(recall = correct / n, recall_given_retained = correct / retained, precision = typed_correct / typed,
             precision_among_tleuk = typed_correct / typed_from_tleuk)]
    RES[[length(RES) + 1]] <- s[, `:=`(lane = L, kind = kd, refs = if (kd == "perpatient") "patient" else "lane", scope = sc, method = mt, level = lv)] }
  ## false-positive sources (held-out runs): singlets pooled only
  FPS[[L]] <- OUT[kind == "heldout" & type != "T|LEUK", .(n = .N, typed_any_tleuk = sum(tleuk_any),
                 by_state = if (startsWith(.BY$type, "singlet:")) NA_character_ else
                   paste(sprintf("%s:%d", sub("LEUK_", "", STATES), vapply(STATES, function(s) sum(reported == s, na.rm = TRUE), 0L)), collapse = ";")),
                 by = .(method, type)][, lane := L]
  rm(x, M); gc(FALSE)
}
SUM <- rbindlist(RES, fill = TRUE)
## the decision table (reading rule in the header)
LEV <- SUM[n >= 30, .(states = .N, pass = all(recall >= 0.7 & !is.na(precision) & precision >= 0.7)), by = .(lane, kind, refs, scope, method, level)]
LEV[, multi := states >= 2]
PRIM <- LEV[kind == "heldout" & refs == "lane" & method == "nnls" & scope == "called_devel"]
dec_lane <- PRIM[, .(any_multi_pass = any(pass & multi)), by = lane]
verdict <- if (nrow(dec_lane) == length(DBL_SMOKE_LANES) && all(dec_lane$any_multi_pass)) "NNLS is worth writing into the proposal" else
  "state typing does not work in these two lanes"
ts("[3] PRIMARY (held-out, lane profiles, nnls, called by devel):")
print(SUM[kind == "heldout" & refs == "lane" & method == "nnls" & scope == "called_devel" & level == "K5",
          .(lane = sub("_.*", "", lane), state = merged, n, correct, recall = round(recall, 2), rec_given_ret = round(recall_given_retained, 2), precision = round(precision, 2))])
print(PRIM[, .(lane = sub("_.*", "", lane), level, states, pass)])
ts("[3] VERDICT by the pre-registered reading rule: ", verdict)
print(LEV[, .(any_multi_pass = any(pass & multi)), by = .(lane = sub("_.*", "", lane), kind, refs, scope, method)][order(lane, kind, refs, scope, method)])
fwrite(SUM, file.path(DBL_TAB_DIR, "dbl_nnls_summary_smoke.csv")); fwrite(LEV, file.path(DBL_TAB_DIR, "dbl_nnls_levels_smoke.csv"))
fwrite(rbindlist(FPS, fill = TRUE), file.path(DBL_TAB_DIR, "dbl_nnls_fp_smoke.csv")); fwrite(rbindlist(GEN), file.path(DBL_TAB_DIR, "dbl_nnls_genes_smoke.csv"))
fwrite(rbindlist(INJ, fill = TRUE), file.path(DBL_TAB_DIR, "dbl_nnls_injected_smoke.csv.gz"))
fwrite(data.table(verdict = verdict, rule = "held-out, lane profiles, nnls, called_devel; level passes if every state with n >= 30 has recall and precision >= 0.7; both lanes need a passing level with >= 2 states"),
       file.path(DBL_TAB_DIR, "dbl_nnls_verdict_smoke.csv"))
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, "dbl_checks_06_smoke.csv"))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
