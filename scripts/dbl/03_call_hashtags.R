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
