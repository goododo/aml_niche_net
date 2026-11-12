#!/usr/bin/env Rscript
# 03b_hto_desk.R ----
# Doublet-pairing redesign, stage 0a: do the P01 hashtag doublets behave as random doublets? Untyped, hashtag level only.
# INPUT  : DBL_LANE_DIR/<lane>.rds of the 3 P01 lanes (stage 01). Only $hto and 6 columns of $bc are kept; $ann, $counts never.
#          DBL_LANE_DIR/<count-1>__hto_calls.csv.gz, DBL_TAB_DIR/dbl_hto_{map,summary}_smoke.csv   stage 03 of 2026-11-09 (regression)
#          DBL_TAB_DIR/dbl_checks_01_p01.csv; TAB_DIR/01_preprocess/03_qc_report__GSE185381.csv (n_raw, n_doublet per donor)
# OUTPUT : DBL_TAB_DIR/dbl_checks_03b_p01.csv      validity checks, md5s, unit tests                          (always)
#          DBL_TAB_DIR/dbl_hto_desk_known_p01.csv   known-answer worlds on synthetic pseudo-lanes              (always)
#          DBL_TAB_DIR/dbl_hto_desk_{map,pairs,gates}_p01.csv   per-lane thresholds + calls md5 / tag-pair counts and
#                                                   fits / gates, reported items, verdict                       (only after (c))
# Usage  : Rscript scripts/dbl/03b_hto_desk.R [--selftest | --dryrun=null | --dryrun=excess] [--force]
#   --selftest       (0), (a) md5 + count-1 regression, (b) on count-1, gating worlds W1-W5 on count-1 class-A singles only
#                    (3 pseudo-lanes of DBL_0A_SELFTEST_N); prints only, writes nothing, never forms a real tag pair.
#   --dryrun=WORLD   the whole run (b)-(d) on 3 fake lanes = count-1 class-A singles (class A) + DBL_0A_SELFTEST_N synthetic
#                    doublets (class C); WORLD null must give "K1 usable", excess (two-tag doublets x2) must give BN-1b. (a)'s
#                    data checks are skipped; outputs go to DBL_SMOKE_DIR/dryrun_0a_<WORLD>/. Selftest and both dryruns run before the tag.
#   --force          only under the rerun rule; earlier outputs are moved to DBL_TAB_DIR/superseded_<job>/, never deleted.
# What it does: recomputes stage 03's hashtag calls on the three P01 lanes, proves the machinery on known-answer worlds, and
#   only then computes G1-G3 on the real tag-pair counts and maps them to the pre-written verdict.
#
# ==== FROZEN RULES: tag prereg-dbl-0a-v1 (plan: reports/2026-11_dbl_pairing_smoke/redesign_2026-11-11/STAGE0A_plan_v2.md) ====
# RUN ORDER. (0) decide() and gates() unit tests. (a) md5 of config_dbl.R and config_paths.R == MD5_FROZEN, before any data is
#   read (mismatch: stop, fix, resubmit; not a rerun); count-1 regression (calls recomputed in memory: cell order, call,
#   donors_pos, class, donor, in_counts identical; n_pos, hto_total identical as integers; thresholds, own-top shares and the
#   class x in_counts x call counts all.equal to the stored files); dbl_checks_01_p01.csv all PASS; per P01 donor over the 3
#   lanes, B_union == the QC report's n_doublet and demuxed == n_raw. (b) per lane: stage 03's known answers (own-top share >=
#   0.80; class-A single calls name the author's donor >= 0.95; class-A cross_donor share < 0.05; in-QC class C cross share >
#   class A); >= 50 class-A single calls per tag; every tag in >= 1 in-QC exactly-two-tag cross call. A count-2/3 lane failing
#   (b) is excluded from (c) and (d) (df 12 per remaining lane; both excluded -> count-1 alone); count-1 failing (b) counts as
#   (a). (c) gating worlds W1-W5 on the remaining lanes. (d) only then G1-G3 on the real counts; a non-converged GLM or a
#   fitted-vs-observed tag margin > 1e-6 x lane total in (d) counts as (c). Failure at (a) or (c): one rerun only for a
#   demonstrated code bug (minimal reproduction showing the code departs from this header; same seeds; dated amendment; old
#   outputs moved to superseded_<job>/), otherwise BN-1 "the cross-donor check could not be validated". Until (c) passes only
#   the checks and known tables are written (booleans, thresholds, synthetic results; never a real tag-pair count).
# CALLS. Stage 03 lines 35-55 verbatim: a tag is positive iff count > max(quantile(bg, 0.995, type 7), 5); bg = its counts on
#   author-assigned barcodes of other donors (any class, in or out of line QC). Sets: in_counts & exactly two positive tags;
#   cross = tags of two donors (19 pairs), two-tag = {AML0160.1, .2} or {AML3762.1, .2}; >= 3-tag calls enter no gate; nCount = nc.
# G1. Per lane glm.fit(Z, y, poisson()) on the 19 cross counts, zeros kept, Z = tag incidence (rank 7). G1 = sum X2 / (12 x
#   lanes), X2 = sum((y - mu)^2 / mu). Fail iff G1 > 2. Its 7 free tag effects absorb loading, hashtag sensitivity and generic
#   per-donor stickiness, so no gate tests generic stickiness.
# G2. Median nc of in-QC exactly-two-tag cross calls / median nc of class-A single calls, both pooled over lanes. Fail iff < 1.4.
# G3. t_l = q75 (type 7) of the lane's class-A single nc (lane-pooled by design); doublet-like iff nc >= t_l. pi_jk = exact share
#   of all pairs of the lane's class-A singles of tags j, k with nc_a + nc_b >= t_l. Refit G1's model on doublet-like cross
#   counts with offset log pi; P = exp(eta_d.1 + eta_d.2) pi_d; O = doublet-like two-tag count; G3 = sum O / sum P over lanes
#   and both split donors, unrounded. Fail iff G3 < 0.67 or G3 > 1.5; equality passes every bar.
# WORLDS (200 replicates each; a pseudo-cohort is 3 pseudo-lanes, one per passing lane, cycling over the passing lanes if one
#   was excluded, so every bar is judged at the 3-lane size). Pseudo-lane of lane l: N_l = round(n_l / r_l) doublets of two
#   class-A single calls drawn uniformly with replacement (same-tag pairs kept); n_l = in-QC exactly-two-tag cross count, r_l = share of 1e5 random pairs called
#   that way. HTO one-ambient: H_a, plus max(H_b[t_b] - m(t_b), 0) on b's tag, m(t) = mean count of t on class-A singles of
#   other donors; nc = nc_a + nc_b; the real thresholds held fixed; reference (t_l, pi, G2 denominator) unmodified.
#   Gating, bars >= 0.95: W1 random: all gates pass, and abs(mean log G3) <= 0.05. W2 + Poisson(n_d) extra (d.1, d.2) doublets
#   (n_d = the replicate's (d.1, d.2) draws): G3 > 1.5. W3 + Poisson(0.75 n_d) single calls of split donor d (uniform over both
#   aliquots) given their sibling tag: raw O/P >= 1.5 and G3 in band. W4 + Poisson(n_AB) extra AML0160-AML3133 doublets: G1 > 2.
#   W5 each exactly-two-tag cross call's nc replaced w.p. 0.7 by a random single's: G2 < 1.4.
#   Reported only, in this order (seed + 10 + k): excess x1.5, x1.75; sub-pool x1.3, then x1.5, then x2 over the 10 sorted
#   donor pairs (AML0160|AML3133 x2 is W4, so 29 rows); borrowing proportional to nc; W1 and W2 with the plain HTO sum.
#   Seeds: DBL_SEED_0A before sizing; + i before Wi; + 10 + k before reported world k; + 99 before the bootstrap; + 60 dry-run lanes.
# REPORTED, NEVER DECIDING. Per lane x split donor O, P, O/P raw and doublet-like; G3 on count-1 alone and on count-2+3 (blind
#   replicate), each with a CI; CI = 2,000 Poisson redraws of each lane's 21 doublet-like counts, refit with the same offsets,
#   type-7 2.5/97.5 % (failed refits counted, dropped); G3 at q90; G3 with pi at kappa 0.8; G3 with pi from class-A + in-QC
#   B_mad singles; G2 on all in-QC cross_donor calls; >= 3-tag counts; kappa_hat, f_hat (plan section 6); class of two-tag
#   calls; nc medians by call; r_low deciles; pooled standardized residual per cross pair; leave-one-donor-pair-out G3 range
#   (a pass whose range crosses 0.67 or 1.5 is "pass, model-sensitive"; verdict unchanged); restricted X2/df; fitted tag
#   shares beside class-A single shares.
# VERDICT. All pass: "K1 usable" -> plan 0b. G1: BN-1 "cross-donor doublets do not follow random pairing after pooling".
#   G2: BN-1 "cross-donor calls do not carry doublet-level RNA" (+ median ratio, kappa_hat, f_hat; with f_hat <= 0.25 real
#   doublets carry less RNA than summed synthetics, so the summed-synthetic calibration K1 checks would not transfer).
#   G3: BN-1b "same-donor doublets of the P01 pool depart from random formation after the sibling aliquots met (O/P, CI);
#   aggregation or a tag artefact cannot be separated from pairing, so a share result could not be read as pairing". The
#   typed two-tag share stays a reported negative control in Stage 4. These hold whatever .1/.2 later turn out to mark.
# MEANING. Gates pool the three lanes (user, 2026-11-11; amends redesign section 7's "count-1"). What .1/.2 mark is not public;
#   most likely two separately hashed tubes of one specimen that met at pooling: G3 tests association formed after they met,
#   and a pass does not exclude aggregation before pooling. The class-A reference was cleaned by the authors (HTODemux
#   multi-tag barcodes and their recoverDoublets neighbours removed). The verdict holds for this calling rule only (0.995 /
#   5 UMI, per-lane thresholds, name-based tag mapping); 03b replaces the planned '03 --desk'; the md5 of the per-lane
#   "cell call" lines lets a later generalised 03 prove identical calls.
# MASKING. A lane object is reduced to $hto and bc[, cell, donor, demuxed, class, in_counts, nc] at load; $ann of count-2/3 is
#   not read before freeze 1; no label file, no author metadata, no composition by tag or call in any output.
# EXPOSED BEFORE THE FREEZE (nobody has computed anything else about P01 cross or two-tag barcodes): count-1 class x call table
#   (in-QC cross 2,945 = A 115 + B_mad 16 + B_union 16 + C 2,798; two-tag 244 = 8 + 4 + 1 + 231; all classes cross 5,874,
#   two-tag 602, negative 1,275), thresholds 5/5/32/29.275/124.275/14/46, own-top shares; ratios from it (>= 1.75x with
#   class-A loads, 2.2x, 1.32x with author-assigned loads, purity-corrected 1.27-1.71x, 0.122 out of QC); scDblFinder called
#   fractions (cross 0.43/0.39, two-tag 0.25/0.24, class-A single 0.0145/0.0047); author-assigned donor counts and shares of
#   count-1..3; count-2/3 barcode classes (01); metadata structure, HTO names; review computations on count-1 class-A singles
#   and synthetic doublets only. No hashtag call existed for count-2/3 before this run.
# ================================================================================================================================

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
MD5_FROZEN <- c(config_dbl = "d933a1b5767b9ccc5316548c31a91f9a", config_paths = "b2c7ca9b7d65f2c39774e6ee2cca2523")
ARGS <- commandArgs(TRUE); FORCE <- "--force" %in% ARGS; SELFTEST <- "--selftest" %in% ARGS
DRY  <- sub("^--dryrun=", "", grep("^--dryrun=", ARGS, value = TRUE))     # "null" or "excess" in a dry run, else empty
stopifnot(length(DRY) == 0 || (length(DRY) == 1 && DRY %in% c("null", "excess") && !SELFTEST))
REAL <- !SELFTEST && length(DRY) == 0
ts  <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {   # record and print one check; returns its pass so callers can collect them
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-92s expect %-8s got %-8s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **"))
  ok }
OUT <- setNames(file.path(DBL_TAB_DIR, paste0(c("dbl_checks_03b", "dbl_hto_desk_known", "dbl_hto_desk_map", "dbl_hto_desk_pairs",
                                                "dbl_hto_desk_gates"), "_p01.csv")), c("checks", "known", "map", "pairs", "gates"))
if (length(DRY)) { OUT[] <- file.path(DBL_SMOKE_DIR, paste0("dryrun_0a_", DRY), basename(OUT)); dir.create(dirname(OUT[[1]]), showWarnings = FALSE) }
stopifnot(all(startsWith(OUT, DBL_SMOKE_DIR)), !startsWith(DBL_SMOKE_DIR, here::here()))   # never write inside the public repo
if (REAL && any(file.exists(OUT))) {
  if (!FORCE) stop("03b outputs exist; a rerun is allowed only under the rerun rule (header), with --force")
  sup <- file.path(DBL_TAB_DIR, paste0("superseded_", Sys.getenv("SLURM_JOB_ID", format(Sys.time(), "%Y%m%d_%H%M%S"))))
  dir.create(sup, showWarnings = FALSE); for (f in OUT[file.exists(OUT)]) stopifnot(file.rename(f, file.path(sup, basename(f))))
  ts("[moved] earlier 03b outputs to ", sup) }

## -- design and functions (nothing here reads data) ----
TAGS <- c("HTO-AML0160.1", "HTO-AML0160.2", "HTO-AML3133", "HTO-AML3762.1", "HTO-AML3762.2", "HTO-AML0612", "HTO-AML0310")
TDON <- sub("\\.[0-9]+$", "", sub("^HTO-", "", TAGS)); NT <- length(TAGS)
PR <- t(combn(NT, 2)); CP <- PR[TDON[PR[, 1]] != TDON[PR[, 2]], ]; TP <- PR[TDON[PR[, 1]] == TDON[PR[, 2]], , drop = FALSE]
pcode <- function(j, k) (j - 1L) * NT + k                        # one integer per tag pair j < k, so tabulate() counts pairs
CPc <- pcode(CP[, 1], CP[, 2]); TPc <- pcode(TP[, 1], TP[, 2])
Z <- matrix(0, nrow(CP), NT); Z[cbind(seq_len(nrow(CP)), CP[, 1])] <- 1; Z[cbind(seq_len(nrow(CP)), CP[, 2])] <- 1
SPLIT <- lapply(seq_len(nrow(TP)), function(i) TP[i, ])          # (d.1, d.2) of AML0160 and of AML3762
DPAIRS <- combn(sort(unique(TDON)), 2, simplify = FALSE)          # the 10 donor pairs
stopifnot(nrow(CP) == 19, nrow(TP) == 2, qr(Z)$rank == NT)
two_codes <- function(P) {   # pair code of every column of a logical tag x barcode matrix that has exactly two positive tags
  w <- which(colSums(P) == 2L)
  rr <- matrix((which(P[, w, drop = FALSE]) - 1L) %% NT + 1L, ncol = 2, byrow = TRUE)   # which() walks column by column
  list(w = w, code = pcode(rr[, 1], rr[, 2])) }
fitP <- function(y, off = numeric(length(y)), Zm = Z) {   # log mu = eta_j + eta_k (+ offset); ok = converged and margins match
  f <- suppressWarnings(glm.fit(Zm, y, family = poisson(), offset = off)); mu <- f$fitted.values
  list(eta = f$coefficients, mu = mu, X2 = sum((y - mu)^2 / mu), ok = isTRUE(f$converged) && all(is.finite(f$coefficients)) &&
       max(abs(crossprod(Zm, y - mu))) <= 1e-6 * max(sum(y), 1)) }
pi_tab <- function(ncs, t, kap = 1) {   # exact P(kap * (nc_a + nc_b) >= t) over ALL pairs of singles of tags j and k
  f <- function(j, k) { a <- ncs[[j]]; b <- ncs[[k]]               # b is sorted: findInterval counts b < t / kap - a
    sum(length(b) - findInterval(t / kap - a, b, left.open = TRUE)) / (length(a) * length(b)) }
  list(cross = mapply(f, CP[, 1], CP[, 2]), two = mapply(f, TP[, 1], TP[, 2])) }
gates <- function(EL, RL, TL = lapply(RL, `[[`, "t"), PIL = lapply(RL, `[[`, "pi")) {
  # EL: per lane list(code, nc) of in-QC exactly-two-tag calls; RL: per-lane class-A single reference. Same code for worlds and data.
  X2 <- 0; O <- 0; P <- 0; Or <- 0; Pr <- 0; ok <- TRUE; per <- list()
  for (l in seq_along(EL)) {
    y  <- tabulate(EL[[l]]$code, NT * NT); yd <- tabulate(EL[[l]]$code[EL[[l]]$nc >= TL[[l]]], NT * NT)
    f  <- fitP(y[CPc]); g <- fitP(yd[CPc], log(PIL[[l]]$cross))
    pr <- exp(f$eta[TP[, 1]] + f$eta[TP[, 2]]); pd <- exp(g$eta[TP[, 1]] + g$eta[TP[, 2]]) * PIL[[l]]$two
    X2 <- X2 + f$X2; Or <- Or + sum(y[TPc]); Pr <- Pr + sum(pr); O <- O + sum(yd[TPc]); P <- P + sum(pd); ok <- ok && f$ok && g$ok
    per[[l]] <- list(y = y, yd = yd, f = f, g = g, pr = pr, pd = pd) }
  ncx <- unlist(lapply(EL, function(e) e$nc[e$code %in% CPc]))
  G <- c(G1 = X2 / (12 * length(EL)), G2 = median(ncx) / median(unlist(lapply(RL, `[[`, "nc"))), G3 = O / P, raw = Or / Pr,
         df = 12 * length(EL), ok = ok)
  attr(G, "per") <- per; G }
CLAUSE <- c(G1 = "BN-1 (G1): cross-donor doublets do not follow random pairing after pooling",
            G2 = "BN-1 (G2): cross-donor calls do not carry doublet-level RNA",
            G3 = "BN-1b (G3): same-donor doublets of the P01 pool depart from random formation after the sibling aliquots met")
decide <- function(G, valid = TRUE) {   # the pre-written verdict; equality passes every bar
  if (!valid) return("BN-1: the cross-donor check could not be validated")
  f <- c(G1 = G[["G1"]] > DBL_0A_G1_MAX, G2 = G[["G2"]] < DBL_0A_G2_MIN, G3 = G[["G3"]] < DBL_0A_G3_BAND[1] || G[["G3"]] > DBL_0A_G3_BAND[2])
  if (!any(f)) "K1 usable" else paste(CLAUSE[f], collapse = " | ") }
call_lane <- function(X) {   # stage 03 lines 35-55 verbatim (bar the X$ unpacking and the returned list)
  H <- X$H; bc <- X$bc
  tags  <- rownames(H)
  tdon  <- sub("\\.[0-9]+$", "", sub("^HTO-", "", tags))          # "HTO-AML0160.2" -> "AML0160"
  donors <- sort(unique(bc[demuxed == TRUE, donor]))
  asg <- bc[, which(demuxed == TRUE)]                               # author-assigned barcodes (classes A and B)
  top <- tdon[max.col(t(as.matrix(H[, asg])), ties.method = "first")]
  MAP <- rbindlist(lapply(seq_along(tags), function(i) {
    bg <- as.numeric(H[i, asg[bc$donor[asg] != tdon[i]]])           # counts of this hashtag on other donors' barcodes
    data.table(hashtag = tags[i], donor = tdon[i], threshold = max(quantile(bg, DBL_HTO_BG_QUANTILE), DBL_HTO_MIN_POS),
               own_top_share = mean(top[bc$donor[asg] == tdon[i]] == tdon[i])) }))
  POS <- as.matrix(H) > MAP$threshold                               # hashtag x barcode, positive?
  calls <- data.table(cell = bc$cell, class = bc$class, donor = bc$donor, in_counts = bc$in_counts,
                      hto_total = as.numeric(Matrix::colSums(H)), n_pos = colSums(POS))
  calls[, donors_pos := vapply(seq_len(ncol(POS)), function(j) paste(sort(unique(tdon[POS[, j]])), collapse = ";"), "")]
  calls[, n_donors := lengths(strsplit(donors_pos, ";", fixed = TRUE))]
  calls[, call := fifelse(n_pos == 0, "negative", fifelse(n_donors >= 2, "cross_donor", fifelse(n_pos >= 2, "same_donor_two_tags", "single")))]
  c(X, list(donors = donors, MAP = MAP, POS = POS, calls = calls)) }
read_lane <- function(L) {   # masking by construction: keep the hashtags and 6 barcode columns; ann and counts dropped at once
  x <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds")))
  bc <- x$bc[, .(cell, donor, demuxed, class, in_counts, nc)]; H <- x$hto; rm(x); gc(FALSE)
  stopifnot(!is.null(H), setequal(colnames(H), bc$cell), identical(rownames(H), TAGS), !anyNA(bc[in_counts == TRUE, nc]))
  call_lane(list(L = L, bc = bc, H = H[, bc$cell, drop = FALSE])) }   # hto columns come in raw GEO order: reorder to bc
ref_lane <- function(X) {   # the lane's class-A single calls: reference for t_l, pi, G2 and the worlds
  i <- which(X$calls$class == "A" & X$calls$n_pos == 1); H <- as.matrix(X$H[, i, drop = FALSE])
  tag <- apply(X$POS[, i, drop = FALSE], 2, which); nc <- X$bc$nc[i]; don <- TDON[tag]   # the one positive tag; its donor
  by <- split(seq_along(tag), factor(tag, levels = seq_len(NT))); ncs <- lapply(by, function(ix) sort(nc[ix]))
  t <- unname(quantile(nc, DBL_0A_DL_Q, type = 7))
  list(H = H, tag = tag, nc = nc, don = don, by = by, ncs = ncs, t = t, pi = pi_tab(ncs, t), thr = X$MAP$threshold,
       m = vapply(seq_len(NT), function(l) mean(H[l, don != TDON[l]]), 0),            # background of tag l (one-ambient)
       own = vapply(seq_len(NT), function(l) median(H[l, tag == l]), 0)) }             # median own-tag count (r_low)
real_E <- function(X) {   # in-QC barcodes with exactly two positive tags: pair code, nCount, row in bc
  tc <- two_codes(X$POS[, X$calls$in_counts, drop = FALSE]); w <- which(X$calls$in_counts)[tc$w]
  list(code = tc$code, nc = X$bc$nc[w], w = w) }
syn_H <- function(R, a, b, sum_hto = FALSE) {   # synthetic doublets of singles a and b: HTO (one-ambient unless sum_hto) and nCount
  H <- R$H[, a, drop = FALSE]; tbb <- R$tag[b]; ix <- cbind(tbb, seq_along(b))
  if (sum_hto) H <- H + R$H[, b, drop = FALSE] else H[ix] <- H[ix] + pmax(R$H[cbind(tbb, b)] - R$m[tbb], 0)   # b adds only its own tag
  list(H = H, nc = R$nc[a] + R$nc[b]) }
pseudo_lane <- function(R, N, world = "random", mult = 1, sub = NULL, sum_hto = FALSE, by_nc = FALSE, fcont = 0) {
  n <- length(R$nc); a <- sample.int(n, N, TRUE); b <- sample.int(n, N, TRUE); ta <- R$tag[a]; tb <- R$tag[b]
  npair <- function(s) sum((ta == s[1] & tb == s[2]) | (ta == s[2] & tb == s[1]))   # this replicate's draws of tag pair s
  draw  <- function(ix, k) ix[sample.int(length(ix), k, TRUE)]                       # k uniform draws from the index set ix
  if (world == "excess") for (s in SPLIT) { k <- rpois(1, (mult - 1) * npair(s)); a <- c(a, draw(R$by[[s[1]]], k)); b <- c(b, draw(R$by[[s[2]]], k)) }
  if (world == "subpool") { k <- rpois(1, (mult - 1) * sum((R$don[a] == sub[1] & R$don[b] == sub[2]) | (R$don[a] == sub[2] & R$don[b] == sub[1])))
    a <- c(a, draw(which(R$don == sub[1]), k)); b <- c(b, draw(which(R$don == sub[2]), k)) }
  D <- syn_H(R, a, b, sum_hto); P <- D$H > R$thr; nc <- D$nc      # thr recycles down each column: row i vs threshold i
  if (world == "borrow") for (s in SPLIT) {                        # singles of split donor d carrying their sibling tag
    k <- rpois(1, (mult - 1) * npair(s)); pool <- c(R$by[[s[1]]], R$by[[s[2]]])
    v <- if (by_nc) pool[sample.int(length(pool), k, TRUE, prob = R$nc[pool])] else draw(pool, k)
    Pv <- R$H[, v, drop = FALSE] > R$thr; Pv[cbind(ifelse(R$tag[v] == s[1], s[2], s[1]), seq_along(v))] <- TRUE
    P <- cbind(P, Pv); nc <- c(nc, R$nc[v]) }
  tc <- two_codes(P); ncw <- nc[tc$w]
  if (fcont > 0) { ic <- which(tc$code %in% CPc); r <- ic[runif(length(ic)) < fcont]; ncw[r] <- R$nc[sample.int(n, length(r), TRUE)] }
  list(code = tc$code, nc = ncw) }
fake_lanes <- function(world) {   # --dryrun: 3 fake lanes = count-1 class-A singles (class A) + synthetic doublets (class C)
  R1 <- ref_lane(read_lane(DBL_P01_LANES[["p01_1"]])); n <- length(R1$nc); set.seed(DBL_SEED_0A + 60)
  lapply(1:3, function(k) { a <- sample.int(n, DBL_0A_SELFTEST_N, TRUE); b <- sample.int(n, DBL_0A_SELFTEST_N, TRUE)
    if (world == "excess") { ta <- R1$tag[a]; tb <- R1$tag[b]           # the known answer: two-tag doublets x2
      for (s in SPLIT) { m <- rpois(1, sum((ta == s[1] & tb == s[2]) | (ta == s[2] & tb == s[1])))
        a <- c(a, R1$by[[s[1]]][sample.int(length(R1$by[[s[1]]]), m, TRUE)]); b <- c(b, R1$by[[s[2]]][sample.int(length(R1$by[[s[2]]]), m, TRUE)]) } }
    D <- syn_H(R1, a, b); H <- cbind(R1$H, D$H); colnames(H) <- c(paste0("s", seq_len(n)), paste0("d", seq_along(a)))
    bc <- data.table(cell = colnames(H), donor = c(R1$don, rep("UNDEMUX", length(a))), demuxed = rep(c(TRUE, FALSE), c(n, length(a))),
                     class = rep(c("A", "C"), c(n, length(a))), in_counts = TRUE, nc = c(R1$nc, D$nc))
    call_lane(list(L = sprintf("dryrun_%s_%d", world, k), bc = bc, H = H)) }) }
run_world <- function(RL, N, seed, a) { set.seed(seed)   # DBL_0A_N_REP pseudo-cohorts of one pseudo-lane per lane of RL
  t(vapply(seq_len(DBL_0A_N_REP), function(r) gates(Map(function(R, n) do.call(pseudo_lane, c(list(R = R, N = n), a)), RL, N), RL)[c("G1", "G2", "G3", "raw", "ok")], numeric(5))) }
summ <- function(M) { ok <- M[, "ok"] == 1; inb <- M[, "G3"] >= DBL_0A_G3_BAND[1] & M[, "G3"] <= DBL_0A_G3_BAND[2]
  c(p_all_pass = mean(ok & M[, "G1"] <= DBL_0A_G1_MAX & M[, "G2"] >= DBL_0A_G2_MIN & inb), mean_log_G3 = mean(log(M[, "G3"])),
    p_G1_fail = mean(ok & M[, "G1"] > DBL_0A_G1_MAX), p_G2_fail = mean(ok & M[, "G2"] < DBL_0A_G2_MIN), p_G3_high = mean(ok & M[, "G3"] > DBL_0A_G3_BAND[2]),
    p_G3_in = mean(ok & inb), p_raw_ge_1.5 = mean(ok & M[, "raw"] >= 1.5), mean_G3 = mean(M[, "G3"]), mean_raw = mean(M[, "raw"]), p_fit_ok = mean(ok)) }
GW <- list(W1 = list(), W2 = list(world = "excess", mult = DBL_0A_MULT[["excess"]]), W3 = list(world = "borrow", mult = DBL_0A_MULT[["borrow"]]),
           W4 = list(world = "subpool", mult = DBL_0A_MULT[["subpool"]], sub = DBL_0A_SUBPOOL), W5 = list(fcont = DBL_0A_CONTAM))
BAR <- list(W1 = function(s) s[["p_all_pass"]] >= DBL_0A_BAR && abs(s[["mean_log_G3"]]) <= DBL_0A_BIAS_MAX,
            W2 = function(s) s[["p_G3_high"]] >= DBL_0A_BAR, W3 = function(s) s[["p_raw_ge_1.5"]] >= DBL_0A_BAR && s[["p_G3_in"]] >= DBL_0A_BAR,
            W4 = function(s) s[["p_G1_fail"]] >= DBL_0A_BAR, W5 = function(s) s[["p_G2_fail"]] >= DBL_0A_BAR)
RW <- c(list(excess_x1.5 = list(world = "excess", mult = 1.5), excess_x1.75 = list(world = "excess", mult = 1.75)),
        unlist(lapply(c(1.3, 1.5, 2), function(m) setNames(lapply(DPAIRS, function(p) list(world = "subpool", mult = m, sub = p)),
                                                          sprintf("subpool_%s_x%s", vapply(DPAIRS, paste, "", collapse = "|"), m))), recursive = FALSE),
        list(borrow_by_nc = list(world = "borrow", mult = DBL_0A_MULT[["borrow"]], by_nc = TRUE), W1_plain_sum = list(sum_hto = TRUE),
             W2_plain_sum = list(world = "excess", mult = DBL_0A_MULT[["excess"]], sum_hto = TRUE)))
RW <- RW[names(RW) != sprintf("subpool_%s_x2", paste(DBL_0A_SUBPOOL, collapse = "|"))]   # that row is the gating W4

## -- [0] unit tests (no data) ----
ts("[0] unit tests")
for (u in list(list(1, 2, 1, "K1 usable"), list(2.01, 2, 1, "BN-1 (G1)"), list(1, 1.39, 1, "BN-1 (G2)"), list(1, 2, 1.51, "BN-1b (G3)"),
               list(1, 2, 0.66, "BN-1b (G3)"), list(2, 1.4, 0.67, "K1 usable"), list(1, 2, 1.5, "K1 usable")))
  chk(sprintf("(0) decide(G1 %s, G2 %s, G3 %s) starts with '%s'", u[[1]], u[[2]], u[[3]], u[[4]]), TRUE,
      startsWith(decide(c(G1 = u[[1]], G2 = u[[2]], G3 = u[[3]])), u[[4]]))
invisible(chk("(0) decide() of an invalid run", TRUE, startsWith(decide(c(G1 = 1, G2 = 2, G3 = 1), valid = FALSE), "BN-1: the cross-donor check could not")))
toyP <- matrix(FALSE, NT, 3); toyP[c(1, 3), 1] <- TRUE; toyP[c(4, 5), 2] <- TRUE; toyP[c(2, 6, 7), 3] <- TRUE
invisible(chk("(0) pair codes of a toy tag matrix", c(pcode(1L, 3L), pcode(4L, 5L)), two_codes(toyP)$code))
fakeE <- list(code = rep(c(CPc, TPc), each = 40L), nc = rep(c(4000L, 9000L), length.out = 21L * 40L))
fakeR <- list(nc = c(3000, 5000, 7000), t = 6000, pi = list(cross = rep(0.9, 19), two = c(0.95, 0.9)))
invisible(chk("(0) gates() on two lanes (one excluded) uses df 24", 24, unname(gates(list(fakeE, fakeE), list(fakeR, fakeR))[["df"]])))
if (!all(vapply(CHK, `[[`, TRUE, "pass"))) stop("unit tests failed: a pre-freeze bug")

## -- [a] provenance and inputs ----
ts("[a] provenance, count-1 regression, inputs")
MD5 <- unname(tools::md5sum(c(here::here("scripts", "dbl", "config_dbl.R"), here::here("scripts", "config", "config_paths.R"))))
if (!chk("(a) md5 of config_dbl.R and config_paths.R == the frozen literals", unname(MD5_FROZEN), MD5))
  stop("config md5 differs from the frozen literals: fix and resubmit (not a rerun); no data was read")
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
invisible(chk(sprintf("(a) md5 of the running copy %s (recorded)", basename(self)), unname(tools::md5sum(self)), unname(tools::md5sum(self))))
okA <- c(rule = chk("(a) calling rule is the 0.995 quantile / 5 UMI", c(0.995, 5), c(DBL_HTO_BG_QUANTILE, DBL_HTO_MIN_POS)))
if (REAL) { c01 <- fread(file.path(DBL_TAB_DIR, "dbl_checks_01_p01.csv"))
  okA["c01"] <- chk("(a) stage 01 checks on count-2/3 all PASS", TRUE, nrow(c01) > 0 && all(c01$pass)) }
XL <- if (length(DRY)) fake_lanes(DRY) else lapply(if (SELFTEST) DBL_P01_LANES["p01_1"] else DBL_P01_LANES, function(L) { ts("    lane ", L); read_lane(L) })
for (X in XL) { ts("    thresholds ", X$L); print(X$MAP) }
if (!length(DRY)) {   # count-1 regression against stage 03's stored files (skipped in a dry run, whose lanes are fake)
  X1 <- XL[[1]]; st <- fread(file.path(DBL_LANE_DIR, paste0(X1$L, "__hto_calls.csv.gz")))
  mp <- fread(file.path(DBL_TAB_DIR, "dbl_hto_map_smoke.csv")); sm <- fread(file.path(DBL_TAB_DIR, "dbl_hto_summary_smoke.csv"))
  okA["reg_calls"] <- chk("(a) count-1 regression: cell order, call, donors_pos, class, donor, in_counts identical", TRUE,
    all(vapply(c("cell", "call", "donors_pos", "class", "donor", "in_counts"), function(v) identical(st[[v]], X1$calls[[v]]), TRUE)))
  okA["reg_int"] <- chk("(a) count-1 regression: n_pos and hto_total identical as integers", TRUE,
    identical(as.integer(st$n_pos), as.integer(X1$calls$n_pos)) && identical(as.integer(st$hto_total), as.integer(X1$calls$hto_total)))
  okA["reg_map"] <- chk("(a) count-1 regression: thresholds and own-top shares all.equal to the stored map", TRUE,
    identical(mp$hashtag, X1$MAP$hashtag) && isTRUE(all.equal(mp$threshold, X1$MAP$threshold)) && isTRUE(all.equal(mp$own_top_share, X1$MAP$own_top_share)))
  okA["reg_sum"] <- chk("(a) count-1 regression: class x in_counts x call counts all.equal to the stored summary", TRUE,
    isTRUE(all.equal(as.data.frame(sm), as.data.frame(X1$calls[, .N, by = .(class, in_counts, call)][order(class, -in_counts, call)]), check.attributes = FALSE)))
  rm(st, mp, sm)
}
if (REAL) { QCR <- fread(file.path(TAB_DIR, "01_preprocess", "03_qc_report__GSE185381.csv")); bcs <- rbindlist(lapply(XL, `[[`, "bc"))
  okA["sums"] <- all(vapply(sort(unique(TDON)), function(s) chk(sprintf("(a) %s over the 3 lanes: B_union == n_doublet, demuxed == n_raw", s),
    as.numeric(QCR[Sample == s, c(n_doublet, n_raw)]), as.numeric(c(bcs[class == "B_union" & donor == s, .N], bcs[demuxed == TRUE & donor == s, .N]))), TRUE)) }

## -- [b] per lane: references and stage 03's known answers ----
ts("[b] per-lane known answers")
okB <- logical(length(XL))
for (l in seq_along(XL)) {   # booleans only before (c): no share or count of a real call is printed here
  X <- XL[[l]]; cl <- X$calls; L <- X$L; R <- XL[[l]]$R <- ref_lane(X)
  ok <- c(chk(sprintf("(b) %s: every hashtag names a deposited donor; every donor has a hashtag", L), TRUE, setequal(TDON, X$donors)),
          vapply(X$donors, function(d) chk(sprintf("(b) %s: %s's hashtag is the top hashtag of >= %.0f %% of its barcodes", L, d, 100 * DBL_HTO_MAP_MIN),
                                           TRUE, X$MAP[donor == d, own_top_share[1]] >= DBL_HTO_MAP_MIN), TRUE))
  fa <- cl[class == "A", mean(call == "cross_donor")]; fc <- cl[class == "C" & in_counts == TRUE, mean(call == "cross_donor")]
  ok <- c(ok, chk(sprintf("(b) %s: class-A 'single' calls name the author's donor (>= 0.95)", L), TRUE, cl[class == "A" & call == "single", mean(donors_pos == donor)] >= 0.95),
          chk(sprintf("(b) %s: class-A barcodes are rarely cross_donor (< 0.05)", L), TRUE, fa < 0.05),
          chk(sprintf("(b) %s: line-QC-passing class C is more often cross_donor than class A", L), TRUE, fc > fa),
          chk(sprintf("(b) %s: >= %d class-A single calls per tag", L, DBL_0A_MIN_SINGLES), TRUE, all(lengths(R$by) >= DBL_0A_MIN_SINGLES)))
  if (!SELFTEST) { e <- XL[[l]]$E <- real_E(X); cc <- e$code[e$code %in% CPc]   # in memory only until (d)
    ok <- c(ok, chk(sprintf("(b) %s: every tag is in >= 1 in-QC exactly-two-tag cross call", L), TRUE, all(tabulate(c(CP[match(cc, CPc), ]), NT) > 0))) }
  okB[l] <- all(ok) }
okA["count1_b"] <- okB[[1]]                                        # count-1 failing (b) counts as (a)
keep <- if (SELFTEST) 1L else c(1L, which(okB)[which(okB) > 1])     # a failing count-2/3 lane is excluded, not a close
if (length(keep) < length(XL)) ts("    excluded lane(s): ", paste(vapply(XL[-keep], `[[`, "", "L"), collapse = ", "))
XL <- XL[keep]; RL <- if (SELFTEST) rep(list(XL[[1]]$R), 3) else lapply(XL, `[[`, "R")

## -- [c] known-answer worlds on synthetic pseudo-lanes ----
KN <- list(); okC <- FALSE
if (all(okA)) {
  ts("[c] known-answer worlds")
  set.seed(DBL_SEED_0A)
  rl <- vapply(RL, function(R) { e <- pseudo_lane(R, DBL_0A_N_SIZE); sum(e$code %in% CPc) / DBL_0A_N_SIZE }, 0)
  N  <- if (SELFTEST) rep(DBL_0A_SELFTEST_N, 3) else round(vapply(XL, function(X) sum(X$E$code %in% CPc), 0) / rl)
  iw <- rep_len(seq_along(RL), 3); RLw <- RL[iw]; Nw <- N[iw]     # always 3 pseudo-lanes, cycling over the passing lanes
  ts(sprintf("    share of random pairs called exactly-two-tag cross r = %s; pseudo-lane sizes N = %s", paste(round(rl, 4), collapse = ";"), paste(Nw, collapse = ";")))
  krow <- function(nm, gating, off, s, met) as.data.table(c(list(world = nm, gating = gating, seed = DBL_SEED_0A + off, N = paste(Nw, collapse = ";")), as.list(s), list(bar_met = met)))
  okW <- logical(length(GW))
  for (i in seq_along(GW)) { s <- summ(run_world(RLw, Nw, DBL_SEED_0A + i, GW[[i]])); okW[i] <- BAR[[i]](s)
    KN[[length(KN) + 1]] <- krow(names(GW)[i], TRUE, i, s, okW[i]); chk(sprintf("(c) gating world %s meets its bar", names(GW)[i]), TRUE, okW[i]) }
  okC <- all(okW)
  if (!SELFTEST) for (k in seq_along(RW)) KN[[length(KN) + 1]] <- krow(names(RW)[k], FALSE, 10 + k, summ(run_world(RLw, Nw, DBL_SEED_0A + 10 + k, RW[[k]])), NA)
  KN <- rbindlist(KN); num <- setdiff(names(KN)[vapply(KN, is.numeric, TRUE)], "seed")
  print(KN[, c(list(world = world, bar_met = bar_met), lapply(.SD, round, 3)), .SDcols = num]) }
VALID <- all(okA) && okC
invisible(chk("VALIDITY: (a), (b) for count-1, and (c) all pass", TRUE, VALID))
if (SELFTEST) { ts(sprintf("[selftest done] %d checks, %d FAIL; nothing written", length(CHK), sum(!vapply(CHK, `[[`, TRUE, "pass")))); quit(save = "no") }
if (!VALID) { fwrite(rbindlist(CHK), OUT[["checks"]]); if (length(KN)) fwrite(KN, OUT[["known"]])
  ts("[INVALID] (a) or (c) failed: no gate computed. Rerun rule (header): one rerun only for a demonstrated code bug, else BN-1"); quit(save = "no") }

## -- [d] the real gates, then the reported items ----
ts("[d] G1-G3 on the real tag-pair counts")
EL <- lapply(XL, `[[`, "E"); G <- gates(EL, RL); per <- attr(G, "per"); nl <- length(XL); LN <- vapply(XL, `[[`, "", "L")
if (!chk("(d) real fits converged and fitted tag margins equal the observed ones", TRUE, G[["ok"]] == 1)) {
  fwrite(rbindlist(CHK), OUT[["checks"]]); fwrite(KN, OUT[["known"]]); ts("[INVALID] a real fit failed: counts as (c); rerun rule"); quit(save = "no") }
verdict <- decide(G)
set.seed(DBL_SEED_0A + 99)                                         # G3's CI: Poisson redraws of the 21 doublet-like counts
B <- t(vapply(seq_len(DBL_0A_N_BOOT), function(r) { o <- p <- numeric(nl); ok <- TRUE
  for (l in seq_len(nl)) { g <- fitP(rpois(19, per[[l]]$yd[CPc]), log(RL[[l]]$pi$cross)); ok <- ok && g$ok
    o[l] <- sum(rpois(2, per[[l]]$yd[TPc])); p[l] <- sum(exp(g$eta[TP[, 1]] + g$eta[TP[, 2]]) * RL[[l]]$pi$two) }
  c(ok, o, p) }, numeric(1 + 2 * nl)))
Bk <- B[B[, 1] == 1, , drop = FALSE]
ci <- function(ix) unname(quantile(rowSums(Bk[, 1 + ix, drop = FALSE]) / rowSums(Bk[, 1 + nl + ix, drop = FALSE]), c(0.025, 0.975), type = 7))
loo <- vapply(DPAIRS, function(dp) { rows <- which(!((TDON[CP[, 1]] == dp[1] & TDON[CP[, 2]] == dp[2]) | (TDON[CP[, 1]] == dp[2] & TDON[CP[, 2]] == dp[1])))
  if (qr(Z[rows, ])$rank < NT) return(NA_real_); O <- 0; P <- 0
  for (l in seq_len(nl)) { g <- fitP(per[[l]]$yd[CPc][rows], log(RL[[l]]$pi$cross[rows]), Z[rows, ])
    O <- O + sum(per[[l]]$yd[TPc]); P <- P + sum(exp(g$eta[TP[, 1]] + g$eta[TP[, 2]]) * RL[[l]]$pi$two) }
  O / P }, 0)
invisible(chk("(d) every leave-one-donor-pair-out design has rank 7", TRUE, !anyNA(loo)))
mixF <- function(x, RLs, W) { num <- 0                             # exact CDF at x of nc_a + nc_b over cross pairs, weights W
  for (l in seq_along(RLs)) for (p in which(W[[l]] > 0)) { a <- RLs[[l]]$ncs[[CP[p, 1]]]; b <- RLs[[l]]$ncs[[CP[p, 2]]]
    num <- num + W[[l]][p] * sum(findInterval(x - a, b)) / (length(a) * length(b)) }   # findInterval counts b <= x - a
  num / sum(unlist(W)) }
mixQ <- function(q, RLs, W) { lo <- 0; hi <- 2 * max(unlist(lapply(RLs, `[[`, "nc")))   # bisection for the q-quantile
  for (it in 1:60) { mid <- (lo + hi) / 2; if (mixF(mid, RLs, W) < q) lo <- mid else hi <- mid }; hi }
Wc <- lapply(per, function(p) p$y[CPc]); ncx <- lapply(EL, function(e) e$nc[e$code %in% CPc])
kap <- unname(quantile(unlist(ncx), 0.75, type = 7)) / mixQ(0.75, RL, Wc); kapc <- min(max(kap, 0.6), 1)
q25 <- vapply(RL, function(R) unname(quantile(R$nc, 0.25, type = 7)), 0)
Fx <- vapply(seq_len(nl), function(l) mean(ncx[[l]] <= q25[l]), 0); Fs <- vapply(seq_len(nl), function(l) mixF(q25[l] / kapc, RL[l], Wc[l]), 0)
fhat <- sum(lengths(ncx) * (Fx - Fs)) / sum(lengths(ncx) * (0.25 - Fs))
sens <- verdict == "K1 usable" && (min(loo) < DBL_0A_G3_BAND[1] || max(loo) > DBL_0A_G3_BAND[2])
note <- paste0(verdict, if (sens) "; pass, model-sensitive" else "",
               if (G[["G2"]] < DBL_0A_G2_MIN) sprintf("; G2 %.3f, kappa_hat %.3f, f_hat %.3f", G[["G2"]], kap, fhat) else "",
               if (grepl("G3", verdict)) sprintf("; G3 %.3f [%.3f, %.3f]", G[["G3"]], ci(seq_len(nl))[1], ci(seq_len(nl))[2]) else "")
CORE <- data.table(item = c("G1", "G2", "G3", "G3_ci_lo", "G3_ci_hi", "raw_two_tag_OP", "df", "verdict"), subset = "pooled",
                   value = c(unname(G[c("G1", "G2", "G3")]), ci(seq_len(nl)), unname(G[c("raw", "df")]), NA), note = c(rep("", 7), note))
fwrite(rbindlist(CHK), OUT[["checks"]]); fwrite(KN, OUT[["known"]]); fwrite(CORE, OUT[["gates"]])   # verdict on disk before the reported items
ts("[write] checks, known, core gates and verdict; the reported items follow")
ext_ncs <- function(X) { i <- which((X$calls$class == "A" | (X$calls$class == "B_mad" & X$calls$in_counts)) & X$calls$n_pos == 1)
  lapply(split(X$bc$nc[i], factor(apply(X$POS[, i, drop = FALSE], 2, which), levels = seq_len(NT))), sort) }
t90 <- lapply(RL, function(R) unname(quantile(R$nc, 0.9, type = 7)))
G3v <- c(q90 = gates(EL, RL, t90, Map(function(R, t) pi_tab(R$ncs, t), RL, t90))[["G3"]],
         kappa_0.8 = gates(EL, RL, PIL = lapply(RL, function(R) pi_tab(R$ncs, R$t, 0.8)))[["G3"]],
         ref_A_plus_Bmad = gates(EL, RL, PIL = Map(function(X, R) pi_tab(ext_ncs(X), R$t), XL, RL))[["G3"]])
GT <- rbindlist(list(
  data.table(item = c("G1", "G2", "G3", "raw_two_tag_OP", "df", "lanes"), subset = "pooled", value = c(unname(G[c("G1", "G2", "G3", "raw", "df")]), nl), note = c("", "", "", "", "", paste(LN, collapse = ";"))),
  data.table(item = c("G3_ci_lo", "G3_ci_hi"), subset = "pooled", value = ci(seq_len(nl)), note = sprintf("%d/%d bootstrap refits failed", nrow(B) - nrow(Bk), nrow(B))),
  data.table(item = "G3", subset = "count-1", value = unname(gates(EL[1], RL[1])[["G3"]]), note = sprintf("CI %.3f-%.3f", ci(1)[1], ci(1)[2])),
  if (nl > 1) data.table(item = "G3", subset = paste(paste(sub("^.*-count-", "count-", LN[-1]), collapse = "+"), "(blind replicate)"), value = unname(gates(EL[-1], RL[-1])[["G3"]]), note = sprintf("CI %.3f-%.3f", ci(2:nl)[1], ci(2:nl)[2])),
  data.table(item = paste0("G3_", names(G3v)), subset = "pooled", value = unname(G3v), note = "reported variant"),
  data.table(item = c("G3_loo_min", "G3_loo_max"), subset = "pooled", value = range(loo), note = paste(sprintf("%s %.3f", vapply(DPAIRS, paste, "", collapse = "|"), loo), collapse = "; ")),
  data.table(item = c("G2_all_cross_donor", "kappa_hat", "f_hat", "restricted_X2_per_df"), subset = "pooled",
             value = c(median(unlist(lapply(XL, function(X) X$bc$nc[X$calls$in_counts & X$calls$call == "cross_donor"]))) / median(unlist(lapply(RL, `[[`, "nc"))),
                       kap, fhat, sum(vapply(per, function(p) p$g$X2, 0)) / (12 * nl)), note = ""),
  data.table(item = "n_ge3_tags_inQC", subset = LN, value = vapply(XL, function(X) sum(X$calls$in_counts & X$calls$n_pos >= 3), 0), note = ""),
  rbindlist(lapply(seq_len(nl), function(l) data.table(item = c("O_raw", "P_raw", "O_doublet_like", "P_doublet_like"), subset = paste(LN[l], rep(TDON[TP[, 1]], each = 4)),
    value = c(rbind(per[[l]]$y[TPc], per[[l]]$pr, per[[l]]$yd[TPc], per[[l]]$pd)), note = ""))),
  rbindlist(lapply(seq_len(nl), function(l) { X <- XL[[l]]
    data.table(nc = X$bc$nc, call = X$calls$call)[X$calls$in_counts][, .(item = "median_nc_by_call", value = as.numeric(median(nc)), note = paste("N", .N)), by = .(subset = paste(LN[l], call))] })),
  rbindlist(lapply(seq_len(nl), function(l) { X <- XL[[l]]; w <- X$E$w[X$E$code %in% TPc]
    data.table(cls = X$calls$class[w])[, .N, by = cls][, .(item = "two_tag_calls_by_class", subset = paste(LN[l], cls), value = N, note = "")] })),
  rbindlist(lapply(seq_len(nl), function(l) { e <- exp(per[[l]]$f$eta)
    data.table(item = "fitted_tag_share", subset = paste(LN[l], TAGS), value = e / sum(e), note = sprintf("class-A single share %.3f", tabulate(RL[[l]]$tag, NT) / length(RL[[l]]$tag))) })),
  data.table(item = "verdict", subset = "pooled", value = NA_real_, note = note)), fill = TRUE)
RLOW <- rbindlist(lapply(seq_len(nl), function(l) { X <- XL[[l]]; E <- X$E; R <- RL[[l]]
  j <- (E$code - 1L) %/% NT + 1L; k <- (E$code - 1L) %% NT + 1L; h <- as.matrix(X$H[, E$w, drop = FALSE])   # decode the pair code
  data.table(kind = ifelse(E$code %in% TPc, "two_tag", "cross"), doublet_like = E$nc >= R$t,
             r_low = pmin(h[cbind(j, seq_along(j))] / R$own[j], h[cbind(k, seq_along(k))] / R$own[k])) }))
GT <- rbind(GT, RLOW[, .(item = "r_low_deciles", value = quantile(r_low, seq(0.1, 0.9, 0.1), type = 7), note = paste0("q", 1:9, "0")),
                       by = .(subset = paste(kind, ifelse(doublet_like, "doublet-like", "not doublet-like")))], fill = TRUE)
Ysum <- Reduce(`+`, lapply(per, function(p) p$y[CPc])); Msum <- Reduce(`+`, lapply(per, function(p) p$f$mu))
PT <- rbind(rbindlist(lapply(seq_len(nl), function(l) { p <- per[[l]]; R <- RL[[l]]
         rbind(data.table(lane = LN[l], tag_j = TAGS[CP[, 1]], tag_k = TAGS[CP[, 2]], kind = "cross", y_all = p$y[CPc], fit_all = p$f$mu,
                          y_doublet_like = p$yd[CPc], fit_doublet_like = p$g$mu, pi = R$pi$cross),
               data.table(lane = LN[l], tag_j = TAGS[TP[, 1]], tag_k = TAGS[TP[, 2]], kind = "two_tag", y_all = p$y[TPc], fit_all = p$pr,
                          y_doublet_like = p$yd[TPc], fit_doublet_like = p$pd, pi = R$pi$two)) })),
            data.table(lane = "pooled", tag_j = TAGS[CP[, 1]], tag_k = TAGS[CP[, 2]], kind = "cross", y_all = Ysum, fit_all = Msum, std_resid = (Ysum - Msum) / sqrt(Msum)), fill = TRUE)
MAPT <- rbindlist(lapply(seq_len(nl), function(l) { X <- XL[[l]]; R <- RL[[l]]; tf <- tempfile(); writeLines(paste(X$calls$cell, X$calls$call), tf)
  m5 <- unname(tools::md5sum(tf)); unlink(tf)                      # md5 of the lane's "cell call" lines, in bc order
  X$MAP[, .(lane = X$L, hashtag, donor, threshold, own_top_share, n_single_A = tabulate(R$tag, NT), bg_mean = R$m, own_median = R$own, calls_md5 = m5)] }))
print(GT[item %in% c("G1", "G2", "G3", "G3_ci_lo", "G3_ci_hi", "kappa_hat", "f_hat", "G3_loo_min", "G3_loo_max") | subset == "count-1" | grepl("blind", subset),
         .(item, subset, value = signif(value, 4))])
print(PT[lane == "pooled", .(tag_j, tag_k, y_all, fit_all = round(fit_all, 1), std_resid = round(std_resid, 2))])
if (length(DRY)) invisible(chk(sprintf("dryrun %s: the verdict is the known answer", DRY), TRUE, startsWith(verdict, c(null = "K1 usable", excess = "BN-1b")[[DRY]])))
TABS <- list(checks = rbindlist(CHK), known = KN, map = MAPT, pairs = PT, gates = GT)
for (nm in names(TABS)) fwrite(TABS[[nm]], OUT[[nm]])
ts("[write] ", paste(basename(OUT), collapse = ", "))
ts(sprintf("[done] %d checks, %d FAIL%s. VERDICT: %s", length(CHK), sum(!vapply(CHK, `[[`, TRUE, "pass")),
           if (length(keep) < 3) " (incl. the (b) rows of the excluded lane(s))" else "", note))
