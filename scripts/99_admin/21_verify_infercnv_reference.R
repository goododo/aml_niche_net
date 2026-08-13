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

FAIL <- 0L; N <- 0L
ok  <- function(m) { N <<- N + 1L; cat(sprintf("  [PASS] %s\n", m)) }
bad <- function(m) { N <<- N + 1L; FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", m)) }
chk <- function(cond, m, d = "") if (isTRUE(cond)) ok(m) else bad(paste0(m, if (nzchar(d)) paste0(" -- ", d)))

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
  chk(length(both) == 0, sprintf("%s: no bin is covered by both sources", ds),
      paste(both, collapse = ","))
  chk(nrow(fb$counts) == length(fb$genes), sprintf("%s: counts rows == shared gene list", ds))
  chk(ncol(fb$counts) == length(fb$groups), sprintf("%s: one group label per reference cell", ds))
  chk(identical(rownames(fb$counts), fb$genes), sprintf("%s: gene order aligned", ds))
}

cat("\n=========== 2.5 group names survive the downstream rules ===========\n")
# 41 keeps a cell as the sample's own iff it is NOT a reference group other than reference_normal,
# and thresholds on every group matching "reference". Both must hold for the new matched labels.
grp <- c("observation", "reference_normal", "reference_external__HSC_MPP", "reference_matched__HSC_MPP")
is_ref <- grepl("reference", grp, ignore.case = TRUE)
chk(identical(is_ref, c(FALSE, TRUE, TRUE, TRUE)), "matched groups count as reference for the threshold")
dropped <- is_ref & grp != "reference_normal"
chk(identical(dropped, c(FALSE, FALSE, TRUE, TRUE)), "matched groups are dropped from per-cell output")
chk(identical(grp[!dropped], c("observation", "reference_normal")),
    "only the sample's own cells survive to the per-cell table")

cat(sprintf("\n=========== BATCH 2: %d checks, %d failed ===========\n", N, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 2 PASS\n")
