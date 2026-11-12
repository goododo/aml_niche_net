#!/usr/bin/env Rscript
# 03_call_hashtags.R ----
# Doublet-pairing SMOKE, stage 3: hashtag (HTO) calls for the pooled lane -- which barcodes carry cells of two donors
# (cross-donor doublets: the proposal's technical-background calibration set, section 4.1 "Background calibration").
# Runs before the doublet typing (stage 04) so that 04 can measure the caller on these barcodes in memory and never write a
# per-barcode real doublet call. (The plan of 2026-11-09 numbered hashtags 04 and typing 03; swapped for that reason.)
# INPUT  : DBL_LANE_DIR/<lane>.rds (stage 01: hto = hashtag rows of the raw GEO matrix x every deposited barcode)
# OUTPUT : DBL_LANE_DIR/<lane>__hto_calls.csv.gz   cell, class, donor (author), hto_total, n_pos, donors_pos, call
#          DBL_TAB_DIR/dbl_hto_map_smoke.csv         hashtag -> donor, background threshold, own-donor top share (known answer)
#          DBL_TAB_DIR/dbl_hto_summary_smoke.csv     class x call barcode counts (hashtag level only; no doublet typing)
#          DBL_TAB_DIR/dbl_checks_03_smoke.csv
# Usage  : Rscript scripts/dbl/03_call_hashtags.R [--force]
# Calls: a hashtag is positive above the DBL_HTO_BG_QUANTILE quantile of its counts among author-assigned barcodes of the
#   OTHER donors (and >= DBL_HTO_MIN_POS UMI). negative = no positive hashtag; single = one donor; same_donor_two_tags = two
#   hashtags of one donor (donors split over two tags: a random within-donor doublet); cross_donor = >= 2 donors.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
FORCE <- "--force" %in% commandArgs(TRUE)
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-78s expect %-12s got %-12s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }

## == MODE --share: Stage 2 step 2c (reports/.../redesign_2026-11-11/STAGE2_plan.md; user approval 2026-11-12) ====
# Hashtag level only (untyped; allowed after freeze 1). Lanes: hashed analysis pools of Stage-1-eligible patients or K1 pools.
# OUTPUT : DBL_LANE_DIR/<lane>__hto_calls_share.csv.gz   per barcode: class, author donor, in_counts, nc, n_pos, donors_pos, call
#          DBL_TAB_DIR/dbl_hto_{map,calls,sensitivity,pairs,lambda,escape}_share.csv, dbl_checks_03_share.csv
# Tag map: a tag whose name (HTO- prefix, .k suffix and punctuation removed) is a deposited donor maps to it; any other tag maps to
#   the author donor of >= DBL_HTO_MAP_MIN of the author-assigned barcodes whose top tag it is (>= DBL_SHARE_HTO_MIN_TOP of them);
#   else it is an orphan tag (cells of a sample not in the deposit), kept as its own pseudo-donor. Thresholds and calls: stage 03.
# Escape and detection: DBL_SHARE_HTO_MC random class-A cell pairs of different donors, hashtag rows summed (the redesign's
#   construction), called with the lane's thresholds; tabulated by donor x lineage compartment of both cells. lambda-hat = called
#   cross / (in-QC barcodes x (1 - sum q^2) x P(detected)), q = single-call donor shares; also on doublet-like calls (0a rule).
if ("--share" %in% commandArgs(TRUE)) {
  stopifnot(system2("git", c("-C", here::here(), "rev-parse", "-q", "--verify", "refs/tags/prereg-dbl-share-f1"), stdout = FALSE) == 0)
  LIN <- c(DBL_LEUK_STATES, "LEUK_other")
  LT <- fread(here::here("scripts", "dbl", DBL_SHARE_LANES_FILE), colClasses = list(character = c("patient", "lane")))[role == "analysis"]
  EL <- fread(file.path(DBL_TAB_DIR, "dbl_share_stage1_eligibility.csv"), colClasses = list(character = "patient"))[eligible == TRUE, patient]
  LANES_ALL <- unique(LT$lane); HL <- LANES_ALL[LANES_ALL %in% LT[hashed == TRUE & (patient %in% EL | k1_lane == TRUE), lane]]
  OUT <- setNames(file.path(DBL_TAB_DIR, paste0("dbl_", c("hto_map", "hto_calls", "hto_sensitivity", "hto_pairs", "hto_lambda", "hto_escape", "hto_owntop", "checks_03"), "_share.csv")),
                  c("map", "calls", "sens", "pairs", "lambda", "escape", "owntop", "checks"))
  if (any(file.exists(OUT)) && !FORCE) stop("outputs exist (use --force)")
  norm <- function(z) toupper(gsub("[^A-Za-z0-9]", "", z))
  R <- list(map = list(), calls = list(), sens = list(), pairs = list(), lambda = list(), escape = list(), owntop = list())
  for (L in HL) {
    ts("[s] ", L)
    x <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds")))
    bc <- x$bc[, .(cell, donor, demuxed, class, in_counts, nc)]; H <- x$hto; rm(x); gc(FALSE)   # hashtags + 6 barcode columns only
    stopifnot(!is.null(H), setequal(colnames(H), bc$cell)); H <- as.matrix(H[, bc$cell, drop = FALSE])
    labA <- fread(file.path(DBL_LANE_DIR, paste0(L, "__labels.csv.gz")), select = c("cell", "class", "compartment"))[class == "A"]
    tags <- rownames(H); donors <- sort(unique(bc[demuxed == TRUE, donor])); asg <- which(bc$demuxed)
    topi <- max.col(t(H[, asg, drop = FALSE]), ties.method = "first")   # top tag of every author-assigned barcode
    MAP <- rbindlist(lapply(seq_along(tags), function(i) {
      tb <- sort(table(bc$donor[asg][topi == i]), decreasing = TRUE)
      nm <- donors[norm(donors) == norm(sub("\\.[0-9]+$", "", sub("^HTO-", "", tags[i])))]
      maj <- if (length(tb)) names(tb)[1] else NA_character_; pur <- if (length(tb)) tb[[1]] / sum(tb) else NA_real_
      d <- if (length(nm) == 1) nm else if (isTRUE(pur >= DBL_HTO_MAP_MIN) && sum(tb) >= DBL_SHARE_HTO_MIN_TOP) maj else paste0("orphan:", tags[i])
      data.table(lane = L, hashtag = tags[i], donor = d, rule = if (length(nm) == 1) "name" else if (startsWith(d, "orphan:")) "orphan" else "majority",
                 n_top = sum(tb), majority = maj, purity = pur) }))
    tdon <- MAP$donor
    MAP[, threshold := vapply(seq_along(tags), function(i) {           # stage 03: 99.5 % of the tag on other donors' author barcodes
      bg <- as.numeric(H[i, asg[bc$donor[asg] != tdon[i]]]); max(quantile(bg, DBL_HTO_BG_QUANTILE), DBL_HTO_MIN_POS) }, 0)]
    ## stage 03's calls (lines 50-55, with this tag map)
    callf <- function(Hm) { P <- Hm > MAP$threshold                  # threshold recycles down each column: row i vs tag i
      dp <- vapply(seq_len(ncol(P)), function(j) paste(sort(unique(tdon[P[, j]])), collapse = ";"), "")
      nd <- lengths(strsplit(dp, ";", fixed = TRUE)); np <- colSums(P)
      data.table(n_pos = np, donors_pos = dp, call = fifelse(np == 0, "negative", fifelse(nd >= 2, "cross_donor", fifelse(np >= 2, "same_donor_two_tags", "single")))) }
    calls <- cbind(bc[, .(cell, class, donor, in_counts, nc)], callf(H))
    ## own-top known answer on threshold-scaled counts (amended 2026-11-12 after GSM5613760/61 failed on raw counts: HTO-AML043 has
    ##   ~5x the others' background, so raw maxima favour it; the calls name the author's donor in 97-99 % of every donor's class-A
    ##   barcodes). Raw shares are reported next to the scaled ones.
    tops <- tdon[max.col(t(H[, asg, drop = FALSE] / MAP$threshold), ties.method = "first")]
    own_top <- vapply(donors, function(d) mean(tdon[topi][bc$donor[asg] == d] == d), 0); own_top_s <- vapply(donors, function(d) mean(tops[bc$donor[asg] == d] == d), 0)
    R$owntop[[L]] <- data.table(lane = L, donor = donors, own_top_raw = own_top, own_top_scaled = own_top_s)
    chk(sprintf("%s: every deposited donor's tags are the top scaled tag of >= %.0f %% of its author barcodes", L, 100 * DBL_HTO_MAP_MIN), TRUE, all(own_top_s >= DBL_HTO_MAP_MIN))
    chk(sprintf("%s: class-A single calls name the author's donor (>= 0.95)", L), TRUE, calls[class == "A" & call == "single", mean(donors_pos == donor)] >= 0.95)
    fa <- calls[class == "A", mean(call == "cross_donor")]; fc <- calls[class == "C" & in_counts == TRUE, mean(call == "cross_donor")]
    chk(sprintf("%s: class-A cross_donor share < 0.05 (%.3f)", L, fa), TRUE, fa < 0.05)
    if (calls[class == "C" & in_counts == TRUE, .N] > 0) chk(sprintf("%s: in-QC class C is more often cross_donor than class A", L), TRUE, fc > fa)
    if (L == DBL_SMOKE_LANES[["pooled"]]) { old <- fread(file.path(DBL_LANE_DIR, paste0(L, "__hto_calls.csv.gz")))   # stage 03 of 2026-11-09
      m <- calls[match(old$cell, cell)]
      chk(sprintf("%s: regression: calls, positive donors and n_pos == stage 03's stored smoke calls", L), TRUE,
          nrow(old) == nrow(calls) && identical(m$call, old$call) && identical(m$donors_pos, fifelse(is.na(old$donors_pos), "", old$donors_pos)) &&
          all(m$n_pos == old$n_pos)) }
    fwrite(calls[, .(cell, class, donor, in_counts, nc, n_pos, donors_pos, call)], file.path(DBL_LANE_DIR, paste0(L, "__hto_calls_share.csv.gz")))
    ## sensitivity: class-A cells keep their own donor's tag (single or two own tags), by lineage compartment
    cA <- merge(calls[class == "A"], labA[, .(cell, comp = fifelse(compartment %in% LIN, "LEUK", compartment))], by = "cell")
    R$sens[[L]] <- cA[, .(n = .N, own_tag = mean(call %in% c("single", "same_donor_two_tags") & donors_pos == donor),
                          negative = mean(call == "negative"), cross = mean(call == "cross_donor")), by = .(donor, comp)][, lane := L]
    ## donor-pair counts (in-QC exactly-two-donor cross calls; all and doublet-like) and propensities
    t_dl <- unname(quantile(calls[class == "A" & call == "single", nc], DBL_0A_DL_Q, type = 7))
    q <- calls[in_counts == TRUE & call == "single", .N, by = .(d = donors_pos)][, q := N / sum(N)]
    DP <- if (nrow(q) >= 2) as.data.table(t(combn(sort(q$d), 2)))[, .(d1 = V1, d2 = V2)] else data.table(d1 = character(), d2 = character())
    cr <- calls[in_counts == TRUE & call == "cross_donor" & lengths(strsplit(donors_pos, ";")) == 2]
    DP[, `:=`(n_all = cr[, sum(donors_pos == paste(d1, d2, sep = ";"))], n_dl = cr[nc >= t_dl, sum(donors_pos == paste(d1, d2, sep = ";"))]), by = .(d1, d2)]
    DP[, `:=`(q1 = q$q[match(d1, q$d)], q2 = q$q[match(d2, q$d)])]
    if (nrow(q) >= 3) { Z <- sapply(sort(q$d), function(d) as.numeric(DP$d1 == d | DP$d2 == d))   # pair x donor incidence
      for (v in c("all", "dl")) { y <- DP[[paste0("n_", v)]]; fit <- glm.fit(Z, y, family = poisson())
        p <- exp(fit$coefficients)
        DP[, paste0("fit_", v) := fit$fitted.values]
        R$pairs_prop[[paste(L, v)]] <- data.table(lane = L, set = v, donor = sort(q$d), propensity = p / sum(p, na.rm = TRUE), q = q$q[match(sort(q$d), q$d)],
                                                  x2_df = sum((y - fit$fitted.values)^2 / pmax(fit$fitted.values, 1e-9)) / max(nrow(DP) - ncol(Z), 1)) } }
    R$pairs[[L]] <- DP[, lane := L]
    ## Monte Carlo cross pairs with summed hashtag rows: detection, escapes, doublet-like share
    set.seed(DBL_SEED_SHARE_HTO + match(L, LANES_ALL))
    ia <- match(cA$cell, bc$cell); a <- sample.int(nrow(cA), DBL_SHARE_HTO_MC, TRUE); b <- sample.int(nrow(cA), DBL_SHARE_HTO_MC, TRUE)
    k <- cA$donor[a] != cA$donor[b]; a <- a[k]; b <- b[k]
    sc <- callf(H[, ia[a], drop = FALSE] + H[, ia[b], drop = FALSE])
    MC <- data.table(da = cA$donor[a], ca = cA$comp[a], db = cA$donor[b], cb = cA$comp[b], call = sc$call, dp = sc$donors_pos,
                     dl = cA$nc[a] + cA$nc[b] >= t_dl)
    sw <- MC$da > MC$db; MC[sw, `:=`(da = db, ca = cb, db = da, cb = ca)]   # unordered donor pair: da < db
    R$escape[[L]] <- MC[, .(n = .N, p_cross = mean(call == "cross_donor"), p_single_a = mean(call == "single" & dp == da),
                            p_single_b = mean(call == "single" & dp == db), p_negative = mean(call == "negative"), p_dl = mean(dl)), by = .(da, ca, db, cb)][, lane := L]
    N_in <- calls[in_counts == TRUE, .N]; den <- N_in * (1 - sum(q$q^2))
    n_cr <- calls[in_counts == TRUE & call == "cross_donor", .N]; n_cr_dl <- calls[in_counts == TRUE & call == "cross_donor" & nc >= t_dl, .N]
    p_det <- MC[, mean(call == "cross_donor")]; p_det_dl <- MC[, mean(call == "cross_donor" & dl)]
    R$lambda[[L]] <- data.table(lane = L, n_in_qc = N_in, donors_pos = nrow(q), sum_q2 = sum(q$q^2), n_cross = n_cr, n_cross_dl = n_cr_dl, t_dl = t_dl,
                                mc_pairs = nrow(MC), p_detect = p_det, p_detect_dl = p_det_dl, lambda_raw = n_cr / den, lambda_det = n_cr / (den * p_det),
                                lambda_dl = n_cr_dl / (den * p_det_dl), dbr_input = min(DBL_DBR_CAP, nrow(bc) * DBL_DBR_PER_BARCODE))
    R$map[[L]] <- MAP; R$calls[[L]] <- calls[, .N, by = .(class, in_counts, call)][, lane := L]
    print(MAP[, .(hashtag, donor, rule, n_top, purity = round(purity, 3), threshold = round(threshold, 1))]); print(R$lambda[[L]])
  }
  CK <- rbindlist(CHK); fwrite(CK, OUT[["checks"]])
  if (all(CK$pass)) for (nm in setdiff(names(OUT), "checks")) fwrite(rbindlist(R[[nm]], fill = TRUE), OUT[[nm]])
  if (all(CK$pass) && length(R$pairs_prop)) fwrite(rbindlist(R$pairs_prop), file.path(DBL_TAB_DIR, "dbl_hto_propensity_share.csv"))
  ts(sprintf("[done] %d lanes, %d checks, %d FAIL%s", length(HL), nrow(CK), sum(!CK$pass), if (all(CK$pass)) "" else ": NOT READ, only the checks are written"))
  quit(status = if (all(CK$pass)) 0L else 1L) }

L   <- DBL_SMOKE_LANES[["pooled"]]
out <- file.path(DBL_LANE_DIR, paste0(L, "__hto_calls.csv.gz"))
if (file.exists(out) && !FORCE) { ts("[skip] ", out, " exists (use --force)"); quit(save = "no") }
x  <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds")))
H  <- x$hto; bc <- x$bc
stopifnot(!is.null(H), setequal(colnames(H), bc$cell))
H  <- H[, bc$cell, drop = FALSE]
## -- [1] hashtag -> donor from the hashtag names, verified on author assignments ----
tags  <- rownames(H)
tdon  <- sub("\\.[0-9]+$", "", sub("^HTO-", "", tags))          # "HTO-AML0160.2" -> "AML0160"
donors <- sort(unique(bc[demuxed == TRUE, donor]))
chk("every hashtag names a deposited donor of the lane", TRUE, all(tdon %in% donors))
chk("every deposited donor has a hashtag", TRUE, all(donors %in% tdon))
asg <- bc[, which(demuxed == TRUE)]                               # author-assigned barcodes (classes A and B)
top <- tdon[max.col(t(as.matrix(H[, asg])), ties.method = "first")]
MAP <- rbindlist(lapply(seq_along(tags), function(i) {
  bg <- as.numeric(H[i, asg[bc$donor[asg] != tdon[i]]])           # counts of this hashtag on other donors' barcodes
  data.table(hashtag = tags[i], donor = tdon[i], threshold = max(quantile(bg, DBL_HTO_BG_QUANTILE), DBL_HTO_MIN_POS),
             own_top_share = mean(top[bc$donor[asg] == tdon[i]] == tdon[i])) }))
print(MAP)
for (d in donors) chk(sprintf("known answer: %s's hashtag is the top hashtag of >= %.0f %% of its assigned barcodes", d, 100 * DBL_HTO_MAP_MIN),
                      TRUE, MAP[donor == d, own_top_share[1]] >= DBL_HTO_MAP_MIN)
## -- [2] per-barcode calls ----
POS <- as.matrix(H) > MAP$threshold                               # hashtag x barcode, positive?
calls <- data.table(cell = bc$cell, class = bc$class, donor = bc$donor, in_counts = bc$in_counts,
                    hto_total = as.numeric(Matrix::colSums(H)), n_pos = colSums(POS))
calls[, donors_pos := vapply(seq_len(ncol(POS)), function(j) paste(sort(unique(tdon[POS[, j]])), collapse = ";"), "")]
calls[, n_donors := lengths(strsplit(donors_pos, ";", fixed = TRUE))]
calls[, call := fifelse(n_pos == 0, "negative", fifelse(n_donors >= 2, "cross_donor", fifelse(n_pos >= 2, "same_donor_two_tags", "single")))]
## -- [3] known answers on the author-assigned barcodes ----
own <- calls[class == "A" & call == "single", mean(donors_pos == donor)]
chk("known answer: class-A 'single' calls name the author's donor (>= 0.95)", TRUE, own >= 0.95)
fa <- calls[class == "A", mean(call == "cross_donor")]; fc <- calls[class == "C" & in_counts == TRUE, mean(call == "cross_donor")]
ts(sprintf("    cross_donor share: class A %.3f, class C (line-QC passing) %.3f", fa, fc))
chk("known answer: class-A barcodes are rarely cross_donor (< 0.05)", TRUE, fa < 0.05)
chk("known answer: line-QC-passing class C is more often cross_donor than class A", TRUE, fc > fa)
SUM <- calls[, .N, by = .(class, in_counts, call)][order(class, -in_counts, call)]
print(dcast(SUM, class + in_counts ~ call, value.var = "N", fill = 0L))
fwrite(calls[, .(cell, class, donor, in_counts, hto_total, n_pos, donors_pos, call)], out); ts("[write] ", out)
fwrite(MAP, file.path(DBL_TAB_DIR, "dbl_hto_map_smoke.csv"))
fwrite(SUM, file.path(DBL_TAB_DIR, "dbl_hto_summary_smoke.csv"))
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, "dbl_checks_03_smoke.csv"))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
