#!/usr/bin/env Rscript
# 21_verify_infercnv_reference.R ----
# BATCH 2 of the pipeline audit: how the inferCNV reference is assembled (00_infercnv_common.R).
#
# An error here is invisible and total: every burden value, every malignancy call, and the healthy
# FPR that JUDGES the reference change all come out of this one matrix. The checks below are the
# invariants that, if violated, would still produce a complete and plausible-looking cohort.
#
# The load-bearing one is leave-one-out. The dataset-matched reference is judged by the healthy
# donor FPR; if a donor's own cells sit in the reference used to score it, that number measures
# the reference against itself and would improve no matter what.
#
#   Rscript scripts/99_admin/21_verify_infercnv_reference.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "02_malignancy", "00_infercnv_common.R"))

FAIL <- 0L; N <- 0L; SKIPPED <- 0L
ok  <- function(m) { N <<- N + 1L; cat(sprintf("  [PASS] %s\n", m)) }
bad <- function(m) { N <<- N + 1L; FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", m)) }
chk <- function(cond, m, d = "") if (isTRUE(cond)) ok(m) else bad(paste0(m, if (nzchar(d)) paste0(" -- ", d)))

# A CHECK THAT CANNOT FAIL MUST NOT REPORT PASS. `reachable` states the precondition under which the
# assertion could actually be violated; when it does not hold the check is SKIPPED and says why.
# Concretely: 2.4 asserts that no bin is covered by BOTH the matched and the BMM source, but
# INFERCNV_MATCHED_REF is FALSE, .assemble_foreign_ref() never enters the matched branch, and the
# foreign block comes back with 0 matched cells (verified: Chen2023 -> 2250 BMM, 0 matched). With one
# source present the conflict is impossible, so the old [PASS] meant "unreachable", not "verified".
skipped <- function(m, why) { SKIPPED <<- SKIPPED + 1L; cat(sprintf("  [SKIP] %s -- %s\n", m, why)) }
chk_reachable <- function(reachable, cond, m, why, d = "") {
  if (!isTRUE(reachable)) skipped(m, why) else chk(cond, m, d)
}

DS <- c("Chen2023", "GSE116256", "GSE185381")   # the datasets with an in-dataset healthy donor

cat("\n=========== 2.1 matched pools are well formed ===========\n")
pools <- list()
for (ds in DS) {
  m <- get_matched_ref(ds); pools[[ds]] <- m
  chk(!is.null(m), sprintf("%s pool builds", ds)); if (is.null(m)) next
  chk(ncol(m$counts) == length(m$bin) && ncol(m$counts) == length(m$donor),
      sprintf("%s: counts/bin/donor lengths agree", ds))
  chk(all(m$bin %in% HIERARCHY_BINS), sprintf("%s: pool bins are legal", ds))
  chk(!anyDuplicated(colnames(m$counts)), sprintf("%s: no duplicate reference cell ids", ds))
  chk(all(grepl("^MATCHREF_", colnames(m$counts))),
      sprintf("%s: reference cells are namespaced (cannot collide with sample cells)", ds))
  chk(all(m$donor %in% m$donors), sprintf("%s: donor vector consistent with donor list", ds))
}

cat("\n=========== 2.2 LEAVE-ONE-OUT (the honesty constraint) ===========\n")
for (ds in DS) {
  m <- pools[[ds]]; if (is.null(m)) next
  for (d in m$donors) {
    blk <- .matched_ref_lineage_block(m, drop_bins = c("T_NK", "B_Plasma"), exclude_donor = d)
    n_self <- if (is.null(blk)) 0L else sum(m$donor[blk$idx] == d)
    chk(n_self == 0L, sprintf("%s :: %s excluded from its own reference", ds, d),
        sprintf("%d own cells leaked", n_self))
  }
}

cat("\n=========== 2.3 a matched block never SHRINKS the reference ===========\n")
ext <- get_external_ref(); bmm_sizes <- table(ext$bin)
for (ds in DS) {
  m <- pools[[ds]]; if (is.null(m)) next
  blk <- .matched_ref_lineage_block(m, drop_bins = c("T_NK", "B_Plasma"))
  if (is.null(blk)) { cat(sprintf("  [note] %s: no bin qualifies -> all BMM\n", ds)); next }
  got <- table(sub("^reference_matched__", "", blk$groups))
  for (b in names(got)) {
    want <- INFERCNV_MATCHED_REF_PER_BIN[[b]]
    chk(as.integer(got[[b]]) == want, sprintf("%s/%s: block is the configured size", ds, b),
        sprintf("got %d want %d", as.integer(got[[b]]), want))
    if (b %in% names(bmm_sizes))
      chk(as.integer(got[[b]]) >= as.integer(bmm_sizes[[b]]),
          sprintf("%s/%s: not smaller than the BMM block it replaces", ds, b),
          sprintf("matched %d vs BMM %d", as.integer(got[[b]]), as.integer(bmm_sizes[[b]])))
  }
}

cat("\n=========== 2.4 matched and BMM never cover the same bin ===========\n")
# If both sources contributed groups for one bin, inferCNV's per-group [min,max] interval would be
# widened by an extra group nobody intended -- silently making the reference more permissive.
fake_counts <- function(genes, n) {
  m <- matrix(0, nrow = length(genes), ncol = n,
              dimnames = list(genes, paste0("obs", seq_len(n)))); m
}
for (ds in DS) {
  m <- pools[[ds]]; if (is.null(m)) next
  s_cnt <- fake_counts(rownames(m$counts), 5L)
  fb <- .assemble_foreign_ref(s_cnt, ds, sid = NA_character_, drop_bins = c("T_NK", "B_Plasma"),
                              ext_ref_fn = get_external_ref, matched_ref_fn = function(x) m)
  chk(!is.null(fb), sprintf("%s: foreign block assembles", ds)); if (is.null(fb)) next
  bin_of_group <- sub("^reference_(matched|external)__", "", fb$groups)
  src_of_group  <- ifelse(grepl("^reference_matched__", fb$groups), "matched", "bmm")
  both <- names(which(tapply(src_of_group, bin_of_group, function(v) uniqueN(v) > 1)))
  n_matched_cells <- sum(src_of_group == "matched")
  chk_reachable(n_matched_cells > 0, length(both) == 0,
                sprintf("%s: no bin is covered by both sources", ds),
                sprintf(paste("INFERCNV_MATCHED_REF=%s and the foreign block has %d matched cells,",
                              "so a two-source conflict cannot arise: this assertion is UNREACHABLE,",
                              "not satisfied"), INFERCNV_MATCHED_REF, n_matched_cells),
                paste(both, collapse = ","))
  chk(nrow(fb$counts) == length(fb$genes), sprintf("%s: counts rows == shared gene list", ds))
  chk(ncol(fb$counts) == length(fb$groups), sprintf("%s: one group label per reference cell", ds))
  chk(identical(rownames(fb$counts), fb$genes), sprintf("%s: gene order aligned", ds))
}

cat("\n=========== 2.5 group names survive the downstream rules ===========\n")
# WAS A TAUTOLOGY. The previous version built a literal group vector, ran grepl("reference", .) on it,
# and asserted grepl's answer. It never touched 41_infercnv_to_percell.R, so the rule could be changed
# there and all three checks would still pass. Replaced with the same invariant read off 41's ACTUAL
# output against its ACTUAL input, over every sample that has both.
BURDEN_DIR <- file.path(LARGE1_DIR, "05_cnv_snv", "infercnv_burden")
pc_files <- list.files(INFERCNV_ROOT, pattern = "__infercnv_percell\\.csv$",
                       recursive = TRUE, full.names = TRUE)
n_pairs <- 0L; bad_rows <- character(0); bad_ref <- character(0); bad_ext <- character(0)
for (f in pc_files) {
  smp <- sub("__infercnv_percell\\.csv$", "", basename(f))
  ds  <- basename(dirname(dirname(f)))
  bf  <- file.path(BURDEN_DIR, ds, paste0(smp, "_infercnv_burden.csv"))
  if (!file.exists(bf)) next
  b <- fread(bf, select = c("cell", "group")); pc <- fread(f, select = c("cell", "is_ref"))
  n_pairs <- n_pairs + 1L
  n_obs <- sum(b$group == "observation"); n_rn <- sum(b$group == "reference_normal")
  if (nrow(pc) != n_obs + n_rn) bad_rows <- c(bad_rows, sprintf("%s/%s(%d!=%d)", ds, smp, nrow(pc), n_obs + n_rn))
  if (sum(as.logical(pc$is_ref)) != n_rn)  bad_ref  <- c(bad_ref,  sprintf("%s/%s", ds, smp))
  ext <- b$cell[grepl("^reference_external__", b$group)]
  if (any(ext %in% pc$cell)) bad_ext <- c(bad_ext, sprintf("%s/%s", ds, smp))
}
chk(n_pairs > 100, sprintf("there are per-cell/burden pairs to check at all (%d)", n_pairs))
chk(length(bad_rows) == 0,
    sprintf("per-cell rows == observation + reference_normal, in all %d samples", n_pairs),
    paste(head(bad_rows, 3), collapse = " "))
chk(length(bad_ref) == 0,
    sprintf("is_ref TRUE count == reference_normal count, in all %d samples", n_pairs),
    paste(head(bad_ref, 3), collapse = " "))
chk(length(bad_ext) == 0,
    sprintf("no reference_external__ cell reaches the per-cell table, in all %d samples", n_pairs),
    paste(head(bad_ext, 3), collapse = " "))

# NEGATIVE TEST: the invariant above must FIRE when it is violated, or it proves nothing.
.b <- data.table(cell = c("c1","c2","c3"), group = c("observation","reference_normal","reference_external__HSC_MPP"))
.pc_good <- data.table(cell = c("c1","c2"), is_ref = c(FALSE, TRUE))
.pc_bad  <- data.table(cell = c("c1","c2","c3"), is_ref = c(FALSE, TRUE, TRUE))   # external leaked in
.rule <- function(b, pc) {
  n_obs <- sum(b$group == "observation"); n_rn <- sum(b$group == "reference_normal")
  ext <- b$cell[grepl("^reference_external__", b$group)]
  nrow(pc) == n_obs + n_rn && sum(as.logical(pc$is_ref)) == n_rn && !any(ext %in% pc$cell)
}
chk(isTRUE(.rule(.b, .pc_good)), "NEGATIVE TEST: the rule accepts a correct per-cell table")
chk(isFALSE(.rule(.b, .pc_bad)), "NEGATIVE TEST: the rule REJECTS an external cell that leaked through")

cat(sprintf("\n=========== BATCH 2: %d checks, %d failed, %d skipped ===========\n", N, FAIL, SKIPPED))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 2 PASS\n")
