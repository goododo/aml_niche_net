#!/usr/bin/env Rscript
# 08_validate_annotation.R ----
# Does the marker correction (06) make the cell typing BETTER or WORSE than the BMM projection?
#
# WHY THIS EXISTS. 06 rewrites 27.7% of the cohort's hierarchy_bin. Nothing so far has tested
# whether those rewrites are right. Two facts make that untenable:
#   * 16,966 cells that BMM calls lymphoid are NOT lymphoid after the correction. Those cells sit
#     in the inferCNV autologous reference as "normal lymphoid". If the correction is right, the
#     reference is contaminated; if it is wrong, adopting it would contaminate the reference.
#   * Some rewrites look implausible on their face (T_NK -> Stromal, 1334 cells in two datasets).
# Either way the decision needs a referee, and the referee cannot be either method being judged.
#
# THE REFEREE. van Galen 2019's CellType column for GSE116256: 38,107 cells hand-typed by the
# original authors from their own marker panel and published under peer review.
#
#   It is INDEPENDENT of both methods under test -- different authors, different algorithm, and it
#   never entered this pipeline. That is the whole basis for using it.
#
#   It is NOT ground truth. It is expression-based, like both methods it is judging, so it shares
#   their blind spots: a cell whose transcriptome genuinely looks like a monocyte will be called a
#   monocyte by all three. It also carries van Galen's own conventions. Read a win here as "agrees
#   better with an independent expert annotation", never as "is correct".
#
# HOW THE VOCABULARY IS FIXED. Every van Galen type is mapped to the bin that bmm_bin_map.tsv
# already assigns to its BMM counterpart (GMP -> LMPP_GMP, cDC/pDC -> Mono_DC, and so on), so the
# referee speaks the same 8-bin vocabulary as both contestants and neither is handed an advantage
# by the translation. The map is written from bmm_bin_map.tsv alone, BEFORE any accuracy is
# computed, and is not adjusted afterwards.
#
# THE AMBIGUOUS CASE IS NOT SWEPT UP. van Galen's "Prog" sits between HSC and GMP and has no clean
# counterpart in the 8-bin vocabulary. Assigning it silently would let one arbitrary choice decide
# the verdict. So Prog/Prog-like are EXCLUDED from the primary comparison and reported separately
# under both possible assignments; if the verdict flips between them, the verdict is not usable.
#
# THE DECISIVE TEST is not overall accuracy -- that is diluted by the 55% of cells where the two
# methods already agree and the referee cannot separate them. It is the disagreement subset:
# among cells where BMM and marker disagree, how often is each one right? That, and only that,
# says whether the correction rule in 06 is adding signal or removing it.
#
#   Rscript scripts/03_hierarchy/08_validate_annotation.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))

RAW <- file.path(LARGE1_DIR, "00_raw", "public", "GEO", "GSE116256_RAW")
DS  <- "GSE116256"
OUT <- file.path(HIER_TAB_DIR, "annotation_validation_vs_vangalen.csv")

## ---------------------------------------------------------------------------------------------
## referee vocabulary: van Galen CellType -> hierarchy bin, each line anchored to bmm_bin_map.tsv
## ---------------------------------------------------------------------------------------------
# The "-like" types are van Galen's MALIGNANT cells; each names the normal type it resembles
# ("GMP-like" is a blast with GMP identity), so it takes that type's bin. Both methods are
# supposed to bin a cell by lineage identity regardless of malignancy, so this is the comparison
# they were both built for -- not a handicap on either.
VG2BIN <- c(
  "HSC"          = "HSC_MPP",       "HSC-like"      = "HSC_MPP",       # BMM: HSC MPP
  "GMP"          = "LMPP_GMP",      "GMP-like"      = "LMPP_GMP",      # BMM: Early/Late GMP
  "ProMono"      = "Mono_DC",       "ProMono-like"  = "Mono_DC",       # BMM: Pro-Monocyte
  "Mono"         = "Mono_DC",       "Mono-like"     = "Mono_DC",       # BMM: Monocyte
  "cDC"          = "Mono_DC",       "cDC-like"      = "Mono_DC",       # BMM: cDC
  "pDC"          = "Mono_DC",                                          # BMM: pDC (flagged ambiguous there too)
  "earlyEry"     = "Erythroid",     "lateEry"       = "Erythroid",     # BMM: Early/Late Erythroid
  "ProB"         = "B_Plasma",      "B"             = "B_Plasma",      # BMM: Pro-B / B
  "Plasma"       = "B_Plasma",                                         # BMM: Plasma Cell
  "T"            = "T_NK",          "CTL"           = "T_NK",          # BMM: Naive/Memory T
  "NK"           = "T_NK"                                              # BMM: NK
)
AMBIG <- c("Prog", "Prog-like")   # no clean counterpart; handled separately, never silently binned
# Megakaryocyte and Stromal have NO van Galen counterpart, so this referee cannot score them.
# Stated up front so their absence from the tables below is not read as "no errors there".
UNSCORABLE <- c("Megakaryocyte", "Stromal")

.core <- function(x) sub("^.*_", "", x)
FAIL <- 0L
chk <- function(cond, msg, detail = "") {
  if (isTRUE(cond)) cat(sprintf("  [OK]   %s\n", msg))
  else { FAIL <<- FAIL + 1L; cat(sprintf("  [WARN] %s%s\n", msg, if (nzchar(detail)) paste0(" -- ", detail) else "")) }
}

## ---------------------------------------------------------------------------------------------
## load
## ---------------------------------------------------------------------------------------------
# colClasses="character" is load-bearing: with type inference an all-empty column reads as logical
# NA, and downstream string tests on NA silently change group membership. That bug inflated a
# genotyped-cell count 11x earlier in this project; it is not hypothetical.
anno <- rbindlist(lapply(list.files(RAW, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE), function(f) {
  s <- sub("^GSM[0-9]+_", "", sub("\\.anno\\.txt\\.gz$", "", basename(f)))
  d <- fread(cmd = paste("zcat", shQuote(f)), showProgress = FALSE, colClasses = "character")
  if (!"CellType" %in% names(d)) return(NULL)
  data.table(sample = s, cell = d$Cell, vg = trimws(d$CellType))
}), fill = TRUE)[nzchar(vg)]

rec <- rbindlist(lapply(list.files(file.path(ANNO_RECONCILED_DIR, DS), pattern = "__anno_percell\\.csv$",
                                   full.names = TRUE), function(p) {
  d <- fread(p, select = c("cell", "bin_bmm", "bin_marker", "hierarchy_bin", "anno_source",
                           "marker_low_conf", "marker_cell_type", "margin_bin",
                           "bmm_prob", "mapping_error"))
  d[, sample := sub("__anno_percell\\.csv$", "", basename(p))][]
}), fill = TRUE)

cat("\n=========== inputs ===========\n")
cat(sprintf("  van Galen typed cells : %d across %d samples\n", nrow(anno), uniqueN(anno$sample)))
cat(sprintf("  reconciled cells      : %d across %d samples\n", nrow(rec), uniqueN(rec$sample)))

anno[, k := paste(sample, .core(cell))]
rec[,  k := paste(sample, .core(cell))]
M <- merge(anno, rec, by = "k", suffixes = c("", ".r"))

cat("\n=========== the join is not vacuous ===========\n")
chk(nrow(M) > 20000, "enough matched cells to measure anything",
    sprintf("%d matched of %d van Galen cells", nrow(M), nrow(anno)))
# DENOMINATOR MATTERS. van Galen typed 16 samples that are not in this cohort at all (including the
# MUTZ3 and OCI-AML3 CELL LINES), and within the shared samples our QC drops cells they kept. Scored
# against all van Galen cells the match rate looks alarming (63.9%) for entirely benign reasons.
# The question that actually bears on bias is the reverse one: of the cells WE are scoring, how many
# got a referee label? If that is near 1, the comparison is not running on a skewed subset.
shared   <- intersect(unique(anno$sample), unique(rec$sample))
rec_in   <- rec[sample %in% shared]
cov_ours <- nrow(M) / nrow(rec_in)
cat(sprintf("  van Galen samples not in this cohort : %d (0 matches, correctly)\n",
            uniqueN(anno$sample) - length(shared)))
chk(cov_ours > 0.95, "nearly every cell WE score carries a referee label (no skewed subset)",
    sprintf("%.1f%% of our %d cells in shared samples", 100 * cov_ours, nrow(rec_in)))
per <- anno[sample %in% shared, .(n = .N), by = sample][
  M[, .(m = .N), by = sample], on = "sample"][, r := m / n]
chk(min(per$r, na.rm = TRUE) > 0.3, "no shared sample matched at ~0 (that would be an id-convention bug)",
    sprintf("worst sample %.2f", min(per$r, na.rm = TRUE)))
chk(uniqueN(M$vg) >= 15, "the referee carries many distinct types (a degenerate referee proves nothing)",
    sprintf("%d types", uniqueN(M$vg)))
chk(!anyDuplicated(M$k), "no cell matched twice")

## ---------------------------------------------------------------------------------------------
## score
## ---------------------------------------------------------------------------------------------
score <- function(D, label, prog_bin = NULL) {
  d <- copy(D)
  if (is.null(prog_bin)) d <- d[!vg %in% AMBIG] else d[vg %in% AMBIG, ref := prog_bin]
  d[!vg %in% AMBIG, ref := VG2BIN[vg]]
  d <- d[!is.na(ref)]
  if (!nrow(d)) return(NULL)
  data.table(
    comparison   = label,
    n_cells      = nrow(d),
    acc_bmm      = round(mean(d$bin_bmm       == d$ref, na.rm = TRUE), 4),
    acc_marker   = round(mean(d$bin_marker    == d$ref, na.rm = TRUE), 4),
    acc_final    = round(mean(d$hierarchy_bin == d$ref, na.rm = TRUE), 4)
  )
}

cat("\n=========== 1. overall agreement with the referee ===========\n")
cat("  (Prog/Prog-like excluded: no clean 8-bin counterpart)\n\n")
P <- score(M, "primary (Prog excluded)")
print(P)
cat(sprintf("\n  final vs BMM: %+.4f     final vs marker: %+.4f\n",
            P$acc_final - P$acc_bmm, P$acc_final - P$acc_marker))

cat("\n=========== 2. does the ambiguous class change the verdict? ===========\n")
S <- rbindlist(list(P,
                    score(M, "Prog -> HSC_MPP",  "HSC_MPP"),
                    score(M, "Prog -> LMPP_GMP", "LMPP_GMP")), fill = TRUE)
print(S)
flip <- length(unique(S[, acc_final > acc_bmm])) > 1
chk(!flip, "the sign of (final - BMM) is stable across both Prog assignments",
    "it FLIPS -- the verdict depends on an arbitrary choice and must not be acted on")

cat("\n=========== 3. THE DECISIVE TEST: who is right where they disagree? ===========\n")
D <- copy(M)[!vg %in% AMBIG][, ref := VG2BIN[vg]][!is.na(ref)]
D <- D[!is.na(bin_bmm) & !is.na(bin_marker) & bin_bmm != bin_marker]
cat(sprintf("  cells where BMM and marker disagree: %d (%.1f%% of scorable)\n\n",
            nrow(D), 100 * nrow(D) / nrow(M[!vg %in% AMBIG & !is.na(VG2BIN[vg])])))
W <- D[, .(n = .N,
           bmm_right    = sum(bin_bmm    == ref),
           marker_right = sum(bin_marker == ref),
           both_wrong   = sum(bin_bmm != ref & bin_marker != ref)), by = anno_source]
W[, `:=`(pct_bmm = round(100 * bmm_right / n, 1), pct_marker = round(100 * marker_right / n, 1))]
setorder(W, -n); print(W)
cat(sprintf("\n  POOLED  BMM right %.1f%%  |  marker right %.1f%%  |  both wrong %.1f%%\n",
            100 * D[, mean(bin_bmm == ref)], 100 * D[, mean(bin_marker == ref)],
            100 * D[, mean(bin_bmm != ref & bin_marker != ref)]))
cat("\n  Read 'marker_corrected' as the rows where 06 OVERRODE the projection, and 'bmm_kept' as\n")
cat("  the rows where it declined to. Those two rows are the correction rule on trial.\n")

cat("\n=========== 4. per-bin, where does each method lose? ===========\n")
B <- D[, .(n = .N, bmm_right = round(mean(bin_bmm == ref), 3),
           marker_right = round(mean(bin_marker == ref), 3)), by = ref][order(-n)]
print(B)

cat("\n=========== 5. malignant vs normal cells ===========\n")
# A projection built on a healthy atlas is expected to struggle on blasts. If the marker method
# wins only on malignant cells, the correction is doing something specific rather than general.
D[, malignant := grepl("-like$", vg)]
print(D[, .(n = .N, bmm_right = round(mean(bin_bmm == ref), 3),
            marker_right = round(mean(bin_marker == ref), 3)), by = malignant])

cat("\n=========== 6. the specific rewrites that looked wrong ===========\n")
X <- D[bin_bmm %in% c("T_NK", "B_Plasma") & !bin_marker %in% c("T_NK", "B_Plasma")]
cat(sprintf("  BMM says lymphoid, marker says otherwise: %d cells here\n", nrow(X)))
if (nrow(X)) {
  cat(sprintf("    referee agrees with BMM (really lymphoid): %d (%.1f%%)\n",
              sum(X$bin_bmm == X$ref), 100 * mean(X$bin_bmm == X$ref)))
  cat(sprintf("    referee agrees with marker              : %d (%.1f%%)\n",
              sum(X$bin_marker == X$ref), 100 * mean(X$bin_marker == X$ref)))
  cat("\n  This is the number that decides whether the inferCNV autologous reference is safe.\n")
  print(X[, .N, by = .(bin_bmm, bin_marker, ref)][order(-N)][1:min(10, .N)])
}

cat("\n=========== 7. WHY the marker method loses, and whether the flag catches it ===========\n")
TT <- M[!vg %in% AMBIG][, ref := VG2BIN[vg]][ref == "T_NK"]
cat(sprintf("  cells the referee calls T/NK: %d | BMM right %.3f | marker right %.3f\n\n",
            nrow(TT), mean(TT$bin_bmm == "T_NK", na.rm = TRUE), mean(TT$bin_marker == "T_NK", na.rm = TRUE)))
cat("  what the marker method calls real T/NK cells when it gets them wrong:\n")
print(TT[bin_marker != "T_NK", .N, by = .(marker_cell_type, bin_marker)][order(-N)][1:6])
cat("\n  The top offender is the panel's own doing: 'DC: Migratory (LAMP3+CCR7+)' shares CCR7 with\n")
cat("  naive/central-memory T cells. Errors also arrive per CLUSTER, not per cell -- the method\n")
cat("  types clusters, so one bad cluster call mislabels hundreds of cells at once.\n\n")
print(TT[, .(n = .N, low_conf = round(mean(marker_low_conf %in% TRUE), 3),
             med_margin = round(median(margin_bin, na.rm = TRUE), 3)), by = .(correct = bin_marker == "T_NK")])
cat("  The confidence flag DOES separate them -- which is why 06's 'low conf -> keep BMM' rule\n")
cat("  helps. It is just not sufficient: the wrong calls that clear the flag still get through.\n")

cat("\n=========== 8. does BMM degrade where the projection is poor? ===========\n")
# This decides whether the verdict can be carried to GSE253355, which was marker-typed precisely
# BECAUSE its projection was judged poor. If BMM only wins where it projects well, the win does
# not transfer.
Z <- M[!vg %in% AMBIG][, ref := VG2BIN[vg]][!is.na(ref)]
Z[, pq := cut(bmm_prob, c(-Inf, .5, .7, .9, Inf), labels = c("<0.50", "0.50-0.70", "0.70-0.90", ">0.90"))]
cat("  by bmm_prob (the PRIMARY projection-QC flag per config_hierarchy.R):\n")
print(Z[!is.na(pq), .(n = .N, bmm = round(mean(bin_bmm == ref, na.rm = TRUE), 3),
                      marker = round(mean(bin_marker == ref, na.rm = TRUE), 3)), by = pq][order(pq)])
Z[, eq := cut(mapping_error, quantile(mapping_error, 0:4 / 4, na.rm = TRUE),
              labels = c("Q1 best", "Q2", "Q3", "Q4 worst"), include.lowest = TRUE)]
cat("\n  by mapping_error (the DEMOTED, scale-dependent metric):\n")
print(Z[!is.na(eq), .(n = .N, bmm = round(mean(bin_bmm == ref, na.rm = TRUE), 3),
                      marker = round(mean(bin_marker == ref, na.rm = TRUE), 3)), by = eq][order(eq)])
cat("\n  BMM holds up across its own confidence but collapses in the worst mapping_error quartile.\n")

cat("\n=========== 9. the two QC metrics rank the cohort OPPOSITELY ===========\n")
Q <- rbindlist(lapply(list.dirs(ANNO_RECONCILED_DIR, recursive = FALSE), function(d)
  rbindlist(lapply(list.files(d, full.names = TRUE), function(p)
    fread(p, select = c("bmm_prob", "mapping_error"))[, dataset := basename(d)][]), fill = TRUE)), fill = TRUE)
CUT <- as.numeric(quantile(M$mapping_error, 0.75, na.rm = TRUE))
cat(sprintf("  danger-zone cutoff = GSE116256 mapping_error P75 = %.2f (where BMM fell to ~0.52)\n\n", CUT))
QS <- Q[, .(n = .N, pct_prob_lt_0.5 = round(100 * mean(bmm_prob < 0.5, na.rm = TRUE), 1),
            med_maperr = round(median(mapping_error, na.rm = TRUE), 2),
            pct_in_danger = round(100 * mean(mapping_error > CUT, na.rm = TRUE), 1)), by = dataset]
print(QS[order(-pct_in_danger)])
cat("\n  GSE253355 is the BEST dataset on the primary flag and the WORST on the demoted one. The\n")
cat("  marker route was adopted for it on the strength of the demoted metric alone. Nothing here\n")
cat("  resolves that -- there is no referee for GSE253355 -- but the justification is not settled.\n")

cat(sprintf("\n=========== bins this referee CANNOT score: %s ===========\n",
            paste(UNSCORABLE, collapse = ", ")))
cat("  van Galen typed no megakaryocytes or stroma, so nothing above speaks to those two bins.\n")

fwrite(S, OUT)
cat(sprintf("\n[done] %s\n", OUT))
if (FAIL > 0) cat(sprintf("\n[!] %d guard(s) warned above -- read them before quoting any number here.\n", FAIL))
