#!/usr/bin/env Rscript
# 07_partner_position.R ----
# Doublet-pairing SMOKE, stage 7 (probe; user decision 2026-11-10 after stage 06 failed): instead of five discrete states, estimate
# WHERE on the BoneMarrowMap pseudotime axis the leukemic half of a T-LEUK doublet sits, and measure how accurately, on
# synthetic doublets whose leukemic parent (and so its true pseudotime) is known.
#   knn   (primary)     genes that leukemic cells express and T cells do not (pooled CPM T < 1, LEUK > 5), log-normalised, PCA
#                       on the reference leukemic singlets; position = median pseudotime of the k = 15 nearest reference singlets
#   nnls  (cross-check) weighted NNLS against T / Tamb / NK and leukemic profiles in 8 pseudotime-quantile bins; position =
#                       weight-averaged bin centre
# References: class-A cells of the five leukemic states (and T/Tamb/NK for nnls) that were NOT held out. The same estimator is
# applied to held-out leukemic SINGLETS too, so pairings and availability are measured with the same ruler (singlet positions
# are singlet information, allowed under D2). Doublets: stage 04's held-out and leave-one-donor-out synthetic doublets
# (rebuilt with its seeds; identity checked) and, in the pooled lane, stage 06's within-donor doublets with per-patient
# references. Scopes as in stage 06: all / called_main / called_devel.
# MASKING [D2]: no real doublet is positioned or typed; no re-entered singlet is used. Stage 04's tables are read without
#   their singlet rows.
# READING RULE (fixed 2026-11-10 before any run; user chose the threshold rho >= 0.7):
#   primary = held-out runs, lane references, knn, T|LEUK synthetic doublets called by devel scDblFinder. A lane passes if the
#   Spearman rho between estimated and true pseudotime is >= 0.7 AND |median(estimated - true)| <= 0.1 x IQR(true). Verdict
#   "the continuous axis is worth designing a pairing test around" only if BOTH lanes pass; otherwise "the continuous axis
#   can't be measured accurately in these two lanes either". Everything else is secondary.
# POST-RUN NOTE (2026-11-10, after job 3687419; code unchanged): the negative test's bar |rho| < 0.1 was mis-specified. Synthetic
#   doublets of one state share few distinct parents (about 8 per run for HSC/MPP in GSM5613785), so under shuffling each state's
#   doublets get one shared random median and the shuffled rho has a wide null. Read it as observed-versus-shuffled instead
#   (GSM5613785 0.72 vs 0.14-0.31; GSM5613744 0.36 vs -0.24-0.13). The verdict does not depend on it.
# INPUT  : DBL_LANE_DIR/<lane>.rds, <lane>__labels.csv.gz; DBL_TAB_DIR/dbl_typing_injected_{main,devel}_smoke.csv.gz (04);
#          LARGE1/02_seurat_objects/03_bmm_projected/GSE185381/<donor>__bmm_percell.csv (predicted_Pseudotime)
# OUTPUT : DBL_TAB_DIR/dbl_position_{summary,singlets,detection,fp,verdict}_smoke.csv, dbl_position_injected_smoke.csv.gz,
#          dbl_checks_07_smoke.csv
# Usage  : Rscript scripts/dbl/07_partner_position.R [--cores=8]

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
  message(sprintf("  %-84s expect %-8s got %-8s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
STATES <- DBL_LEUK_STATES
POS_K       <- 15L   # [SMOKE] neighbours (primary); k = 5 reported as a sensitivity
POS_K_SENS  <- 5L
POS_NPC     <- 20L   # [SMOKE] PCs on the leukemic-gene space
POS_NBIN    <- 8L    # [SMOKE] pseudotime-quantile bins for the nnls cross-check
bad_gene <- function(g) grepl("^MT-|^RP[SL]|^HB[^P]", g)
nz <- function(z) fifelse(z == "", NA_character_, z)

## -- stage 04's synthetic doublets, rebuilt (construction code identical to 04 / 06) ----
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
build_pairs_only <- function(src, run_id) {   # the pair part of stage 04's held-out construction (singlets come after: RNG unaffected)
  prs <- rbindlist(c(
    lapply(STATES, function(s) build_pairs(src, DBL_SYN_PER_STATE, "T|LEUK", "T", s, s)),
    lapply(c("Tamb", "NK", "B", "Ery", "MonoDC", "LEUK_other", "T"), function(cp) build_pairs(src, DBL_SYN_PER_OTHER, paste0("T|", cp), "T", cp)),
    lapply(c("Tamb", "NK", "B", "Ery", "MonoDC"), function(cp) build_leuk_pairs(src, DBL_SYN_PER_OTHER, paste0(cp, "|LEUK*"), cp)),
    list(build_leuk_pairs(src[intersect(names(src), c(STATES, "LEUK_other"))] |> (\(s) { s[["LEUKx"]] <- unlist(s[intersect(STATES, names(s))]); s })(),
                          DBL_SYN_PER_OTHER, "LEUK*|LEUK*", "LEUKx"))), fill = TRUE)
  prs[, syn := sprintf("SYN_%s_%04d", run_id, .I)]
  prs }

INJ04 <- rbindlist(lapply(c("main", "devel"), function(m) {
  x <- fread(file.path(DBL_TAB_DIR, paste0("dbl_typing_injected_", m, "_smoke.csv.gz")))
  x <- x[!startsWith(type, "singlet:")]                  # masked: stage 04's singlet rows are never loaded past this line
  x[, .(lane, run, syn, type, true_state = nz(true_state), share, called, mode = m)] }))
BMM_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "03_bmm_projected", DBL_DATASET)
SUMM <- list(); SING <- list(); DET <- list(); FPT <- list(); INJ <- list()
for (L in DBL_SMOKE_LANES) {
  ts("[1] lane ", L)
  x   <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds"))); M <- x$counts
  lab <- fread(file.path(DBL_LANE_DIR, paste0(L, "__labels.csv.gz")))[match(colnames(M), cell)]
  stopifnot(identical(lab$cell, colnames(M)))
  A <- lab[class == "A"]
  pt <- rbindlist(lapply(sort(unique(A$donor)), function(s) fread(file.path(BMM_DIR, paste0(s, "__bmm_percell.csv")), select = c("cell", "predicted_Pseudotime"))))
  A <- merge(A, pt, by = "cell", all.x = TRUE, sort = FALSE)[match(lab[class == "A", cell], cell)]   # keep stage 04's row order
  setnames(A, "predicted_Pseudotime", "pseudotime")
  chk(sprintf("%s: every class-A cell of the five states has a pseudotime", L), 0L, A[compartment %in% STATES & is.na(pseudotime), .N])
  PT <- setNames(A$pseudotime, A$cell)
  ## stage 04's ambient sets (identity check on 'share' only)
  pools <- split(A$cell, A$compartment)
  cpm <- function(cells) { s <- M[, cells, drop = FALSE]; Matrix::rowSums(s) / sum(s) * 1e6 }
  badg <- bad_gene(rownames(M))
  cT  <- cpm(pools[["T"]]); cL <- cpm(unlist(pools[intersect(STATES, names(pools))]))
  GT4 <- rownames(M)[!badg & cT > 5 & cL < 1]; GL4 <- rownames(M)[!badg & cL > 5 & cT < 1]
  leuk_share <- function(Mx) { l <- Matrix::colSums(Mx[GL4, , drop = FALSE]); t <- Matrix::colSums(Mx[GT4, , drop = FALSE]); ifelse(l + t > 0, l / (l + t), NA_real_) }
  set.seed(DBL_SEED_SYN)
  A[, fold := ((sample(.N) - 1L) %% DBL_N_FOLDS) + 1L, by = compartment]
  donors <- sort(unique(A$donor))
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
      set.seed(DBL_SEED_SYN + 5000L + R$di * 100L + R$fold * 10L)       # stage 06's within-donor doublets
      held <- A[donor == R$donor & fold == R$fold, cell]; refc <- A[donor == R$donor & !(cell %in% held)]
    } else {
      set.seed(DBL_SEED_SYN + if (R$kind == "lodo") 999L else R$fold * 10L + R$batch)
      held <- if (R$kind == "lodo") A[donor == R$donor, cell] else A[fold == R$fold, cell]
      refc <- A[!(cell %in% held)] }
    prs <- build_pairs_only(split(A[cell %in% held, cell], A[cell %in% held, compartment]), R$id)
    Mi  <- M[, prs$a, drop = FALSE] + M[, prs$b, drop = FALSE]; colnames(Mi) <- prs$syn
    prs[, `:=`(share = leuk_share(Mi), true_pt = PT[b])]               # b = the leukemic parent of every '|LEUK' type
    ## the leukemic-gene space, from reference cells only
    rL <- refc[compartment %in% STATES]; rT <- refc[compartment == "T"]
    cTr <- cpm(rT$cell); cLr <- cpm(rL$cell)
    GL <- rownames(M)[!badg & cLr > 5 & cTr < 1]
    lnorm <- function(Mx) { Mx <- Mx[GL, , drop = FALSE]; s <- pmax(Matrix::colSums(Mx), 1); log1p(t(t(as.matrix(Mx)) / s) * 1e4) }
    Xr <- lnorm(M[, rL$cell, drop = FALSE]); mu <- rowMeans(Xr); sdv <- apply(Xr, 1, sd); sdv[sdv == 0] <- 1
    pc <- prcomp(t((Xr - mu) / sdv), rank. = POS_NPC, center = FALSE, scale. = FALSE)
    proj <- function(Mx) t((lnorm(Mx) - mu) / sdv) %*% pc$rotation
    knn_pos <- function(Q, k, ptref = rL$pseudotime, exclude_self = FALSE) {
      nn <- RANN::nn2(data = pc$x, query = Q, k = k + exclude_self)$nn.idx
      if (exclude_self) nn <- nn[, -1, drop = FALSE]
      apply(matrix(ptref[nn], nrow = nrow(nn)), 1, median) }
    Qd <- proj(Mi)
    prs[, `:=`(est_knn = knn_pos(Qd, POS_K), est_knn5 = knn_pos(Qd, POS_K_SENS))]
    ## negative test: reference pseudotimes shuffled -> the estimate must lose its relation to the truth
    set.seed(DBL_SEED_SYN + 7100L)
    ti  <- which(prs$type == "T|LEUK")                                 # indices taken OUTSIDE data.table: inside j, 'type' would be the subset's
    neg <- cor(knn_pos(Qd[ti, , drop = FALSE], POS_K, ptref = sample(rL$pseudotime)), prs$true_pt[ti], method = "spearman")
    ## same ruler: held-out leukemic singlets and leave-one-out references
    hs <- A[cell %in% held & compartment %in% STATES]
    sing <- data.table(run = R$id, kind = R$kind, n_heldout_singlets = nrow(hs),
      rho_heldout_singlets = if (nrow(hs) >= 10) cor(knn_pos(proj(M[, hs$cell, drop = FALSE]), POS_K), hs$pseudotime, method = "spearman") else NA_real_,
      rho_ref_loo = cor(knn_pos(pc$x, POS_K, exclude_self = TRUE), rL$pseudotime, method = "spearman"), n_ref = nrow(rL), n_genes = length(GL), neg_rho = neg)
    ## nnls cross-check: T / Tamb / NK + leukemic pseudotime bins
    br <- unique(quantile(rL$pseudotime, seq(0, 1, length.out = POS_NBIN + 1)))
    rL[, bin := cut(pseudotime, br, include.lowest = TRUE, labels = FALSE)]
    centres <- rL[, .(centre = median(pseudotime)), by = bin][order(bin)]
    GT <- rownames(M)[!badg & cTr > 5 & cLr < 1]; G <- union(GL, GT)
    prof <- function(cells) { s <- Matrix::rowSums(M[G, cells, drop = FALSE]) + 0.5; s / sum(s) }
    P <- cbind(sapply(intersect(c("T", "Tamb", "NK"), refc[, .N, by = compartment][N >= 5, compartment]), function(cp) prof(refc[compartment == cp, cell])),
               sapply(centres$bin, function(bb) prof(rL[bin == bb, cell])))
    colnames(P)[(ncol(P) - nrow(centres) + 1):ncol(P)] <- paste0("bin", centres$bin)
    w <- 1 / sqrt(rowMeans(P) + 1e-6); Pw <- P * w
    Y <- as.matrix(Mi[G, , drop = FALSE]); Y <- sweep(Y, 2, pmax(colSums(Y), 1), "/")
    bins <- grep("^bin", colnames(P))
    prs[, est_nnls := apply(Y, 2, function(y) { b <- nnls::nnls(Pw, y * w)$x[bins]; if (sum(b) > 0) sum(b * centres$centre) / sum(b) else NA_real_ })]
    prs[, `:=`(run = R$id, kind = R$kind, refs = R$refs, iqr_ref = IQR(rL$pseudotime))]
    list(out = prs[, .(run, kind, refs, syn, type, true_state, share, true_pt, est_knn, est_knn5, est_nnls, iqr_ref)], sing = sing) }
  ts(sprintf("[2] %d runs on %d cores", length(RUNS), min(CORES, length(RUNS))))
  RR <- mclapply(RUNS, function(R) tryCatch(one_run(R), error = function(e) list(error = conditionMessage(e), id = R$id)),
                 mc.cores = min(CORES, length(RUNS)), mc.preschedule = FALSE)
  RR <- Map(function(r, R) if (is.list(r) && !inherits(r, "try-error")) r else list(error = paste("worker died:", as.character(r)), id = R$id), RR, RUNS)
  for (r in Filter(function(r) !is.null(r$error), RR)) ts("    ** run ", r$id, " failed: ", r$error)
  chk(sprintf("%s: every position run finished", L), 0L, length(Filter(function(r) !is.null(r$error), RR)))
  RR <- Filter(function(r) is.null(r$error), RR); if (!length(RR)) next
  OUT <- rbindlist(lapply(RR, `[[`, "out"))[, lane := L]
  SG  <- rbindlist(lapply(RR, `[[`, "sing"))[, lane := L]; SING[[L]] <- SG; print(SG)
  ## identity with stage 04 (non-perpatient rows)
  idc <- merge(OUT[kind != "perpatient", .(run, syn, type, true_state, share)], INJ04[lane == L & mode == "main", .(run, syn, type4 = type, ts4 = true_state, share4 = share)], by = c("run", "syn"))
  chk(sprintf("%s: rebuilt doublets == stage 04's (rows)", L), INJ04[lane == L & mode == "main", .N], nrow(idc))
  chk(sprintf("%s: rebuilt doublets == stage 04's (type, true state, share, NA-aware)", L), TRUE,
      idc[, all(type == type4 & ((is.na(true_state) & is.na(ts4)) | (!is.na(true_state) & !is.na(ts4) & true_state == ts4)) &
                ((is.na(share) & is.na(share4)) | (!is.na(share) & !is.na(share4) & abs(share - share4) < 1e-9)))])
  chk(sprintf("%s: negative test -- shuffled reference pseudotimes give |rho| < 0.1 (every run)", L), TRUE, SG[, all(abs(neg_rho) < 0.1)])
  called <- dcast(INJ04[lane == L], run + syn ~ mode, value.var = "called")
  OUT <- merge(OUT, called, by = c("run", "syn"), all.x = TRUE)
  setnames(OUT, c("main", "devel"), c("called_main", "called_devel"), skip_absent = TRUE)
  INJ[[L]] <- OUT
  ## known answer: the estimated position rises from HSC-like to promonocyte-like (medians of T|LEUK doublets, held-out runs)
  med <- OUT[kind == "heldout" & type == "T|LEUK", .(m = median(est_knn)), by = true_state][match(STATES, true_state), m]
  chk(sprintf("%s: knn position rises with the state order (HSC/MPP <= LMPP < GMP < ProMono medians)", L), TRUE,
      isTRUE(max(med[1:2], na.rm = TRUE) <= med[3] && med[3] < med[4] && med[4] < med[5]))
  ## accuracy
  for (sc in c("all", "called_main", "called_devel")) for (kd in unique(OUT$kind)) for (mt in c("est_knn", "est_knn5", "est_nnls")) {
    if (sc != "all" && kd == "perpatient") next
    D <- OUT[kind == kd & type == "T|LEUK" & !is.na(get(mt))]
    if (sc != "all") D <- D[get(sc) %in% TRUE]
    if (nrow(D) < 10) next
    iq <- IQR(D$true_pt)
    SUMM[[length(SUMM) + 1]] <- D[, .(lane = L, kind = kd, refs = refs[1], scope = sc, method = sub("^est_", "", mt), n = .N,
      rho = cor(get(mt), true_pt, method = "spearman"), mae_iqr = median(abs(get(mt) - true_pt)) / iq,
      bias_iqr = median(get(mt) - true_pt) / iq, slope_rank = unname(coef(lm(rank(get(mt)) ~ rank(true_pt)))[2]))] }
  ## detection rate along the true position (stage 04 calls), T|LEUK held-out doublets
  DET[[L]] <- melt(OUT[kind == "heldout" & type == "T|LEUK", .(true_pt, called_main, called_devel)], id.vars = "true_pt", variable.name = "mode", value.name = "called")[
    , .(n = .N, called = mean(called %in% TRUE)), by = .(mode, quintile = cut(true_pt, unique(quantile(true_pt, 0:5 / 5)), include.lowest = TRUE, labels = FALSE))][, lane := L][order(mode, quintile)]
  ## false-positive sources: how well their leukemic partner is positioned (held-out, knn)
  FPT[[L]] <- OUT[kind == "heldout" & grepl("\\|LEUK\\*$", type) & !startsWith(type, "LEUK*"), .(n = .N, rho = cor(est_knn, true_pt, method = "spearman")), by = type][, lane := L]
  rm(x, M); gc(FALSE)
}
SUM <- rbindlist(SUMM, fill = TRUE)
SUM[, pass := rho >= 0.7 & abs(bias_iqr) <= 0.1]
PRIM <- SUM[kind == "heldout" & refs == "lane" & method == "knn" & scope == "called_devel"]
verdict <- if (nrow(PRIM) == length(DBL_SMOKE_LANES) && all(PRIM$pass)) "the continuous axis is worth designing a pairing test around" else
  "the continuous axis can't be measured accurately in these two lanes either"
ts("[3] PRIMARY (held-out, lane references, knn, called by devel):")
print(PRIM[, .(lane = sub("_.*", "", lane), n, rho = round(rho, 3), mae_iqr = round(mae_iqr, 3), bias_iqr = round(bias_iqr, 3), pass)])
ts("[3] VERDICT by the pre-registered reading rule: ", verdict)
print(SUM[, .(lane = sub("_.*", "", lane), kind, refs, scope, method, n, rho = round(rho, 2), mae_iqr = round(mae_iqr, 2), bias_iqr = round(bias_iqr, 2))][order(lane, kind, scope, method)])
fwrite(SUM, file.path(DBL_TAB_DIR, "dbl_position_summary_smoke.csv")); fwrite(rbindlist(SING), file.path(DBL_TAB_DIR, "dbl_position_singlets_smoke.csv"))
fwrite(rbindlist(DET), file.path(DBL_TAB_DIR, "dbl_position_detection_smoke.csv")); fwrite(rbindlist(FPT), file.path(DBL_TAB_DIR, "dbl_position_fp_smoke.csv"))
fwrite(rbindlist(INJ, fill = TRUE), file.path(DBL_TAB_DIR, "dbl_position_injected_smoke.csv.gz"))
fwrite(data.table(verdict = verdict, rule = "held-out, lane references, knn (k 15), called by devel; lane passes if Spearman rho >= 0.7 and |median bias| <= 0.1 x IQR(true); both lanes must pass"),
       file.path(DBL_TAB_DIR, "dbl_position_verdict_smoke.csv"))
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, "dbl_checks_07_smoke.csv"))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
