#!/usr/bin/env Rscript
# 99_genotype_concordance.R ----
# The only ORTHOGONAL ground truth in this project. GSE116256 is the one dataset whose authors
# published per-cell targeted genotyping, so it is the only place the malignancy calls can be
# checked against DNA rather than against another expression-based caller.
#
# Result as of 2026-08-18 (see also 96 for the healthy-donor arm). BOTH DENOMINATORS ARE GIVEN:
# the autologous inferCNV reference cells are ASSIGNED malignant=0 and cannot be positive, and
# because the reference set IS the lymphoid population van Galen genotypes as wt-only, the
# contamination is asymmetric (381 of 1685 wt-only vs 10 of 628 mut+) and pushes the flag rate
# DOWN -- i.e. the padded number flatters the caller. 0.484/0.3015 were reported here until
# 2026-08-18 without that qualification.
#   sensitivity on genotyped mut+ cells .................. 0.484  (evaluable cells only: 0.492)
#   flag rate on wt-only cells ........................... 0.302  (evaluable cells only: 0.390)
#   flag rate on known-normal non-reference lymphoid ..... 0.639   (i.e. inverted)
#   pooled AUC (mut+ vs known-normal lymphoid) ........... 0.394
#   AUC within nGene quartiles ........................... 0.56 / 0.45 / 0.41 / 0.44
# The pooled inversion is largely a DEPTH artefact -- cor(burden, nUMI) = -0.41, and mut+ cells
# are far deeper (median 1364 genes) than the lymphoid comparison set (785). Stratifying by
# depth removes the inversion but leaves AUC ~= 0.5: on this dataset the burden carries no
# usable discriminative signal, which is why neither a different threshold nor a larger
# reference block can rescue it.
#
# CAUTION when reading a rerun: fread() must be given colClasses = "character". An all-empty
# MutTranscripts column otherwise reads as logical NA, and nzchar(NA) is TRUE, which silently
# promotes every ungenotyped cell into the mut+ set.
#
#   Rscript scripts/02_malignancy/99_genotype_concordance.R
#
# ---------------------------------------------------------------------------------------------
# Orthogonal anchoring of the malignancy calls against van Galen 2019 single-cell targeted
# GENOTYPING (GSE116256 *.anno.txt.gz, MutTranscripts / WtTranscripts).
#
# This is DNA-level evidence produced by targeted PCR amplification of known AML driver
# mutations, entirely independent of the expression-based CNV inference it is used to judge.
#
# ASYMMETRY OF THE TRUTH SET -- read before using the numbers:
#   mut+ (>=1 mutant transcript)  = STRONG positive. The cell carries the driver lesion.
#   wt-only (wt transcripts, no mutant) = WEAK negative. Most AML driver mutations here are
#     heterozygous, so a malignant cell carries a wt allele too; detecting only wt is fully
#     compatible with allelic dropout. Treat the wt-only flag rate as an UPPER BOUND on the
#     false-positive rate, not as a specificity measurement.
# Sensitivity on mut+ cells is therefore the number this analysis can actually establish.

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
RAW <- "/LARGE1/gr10634/gaozy/aml_niche_net/00_raw/public/GEO/GSE116256_RAW"

tasks <- fread("infercnv_tasks.tsv", header = FALSE, col.names = c("dataset", "sample"))
ours  <- tasks[dataset == "GSE116256"]$sample
message(sprintf("[0] %d GSE116256 samples in the cohort", length(ours)))

geno <- rbindlist(lapply(list.files(RAW, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE), function(f) {
  s <- sub("^GSM[0-9]+_", "", sub("\\.anno\\.txt\\.gz$", "", basename(f)))
  if (!s %in% ours) return(NULL)
  d <- tryCatch(fread(cmd = paste("zcat", shQuote(f)), showProgress = FALSE, colClasses = "character"),
                error = function(e) NULL)
  if (is.null(d) || !all(c("Cell", "MutTranscripts", "WtTranscripts") %in% names(d))) return(NULL)
  mut <- nzchar(trimws(d$MutTranscripts)); wt <- nzchar(trimws(d$WtTranscripts))
  data.table(sample = s, cell_anno = d$Cell,
             celltype_vg = if ("CellType" %in% names(d)) d$CellType else NA_character_,
             truth = fifelse(mut, "mut_pos", fifelse(wt, "wt_only", NA_character_)))[!is.na(truth)]
}), fill = TRUE)
message(sprintf("[1] genotyped cells in cohort samples: %d (mut+ %d, wt-only %d)",
                nrow(geno), sum(geno$truth == "mut_pos"), sum(geno$truth == "wt_only")))

# our per-cell calls
calls <- rbindlist(lapply(unique(geno$sample), function(s) {
  f <- file.path(INFERCNV_ROOT, "GSE116256", s, paste0(s, "__infercnv_percell.csv"))
  if (!file.exists(f)) return(NULL)
  hdr <- names(fread(f, nrows = 0L))
  d <- fread(f, select = if ("is_ref" %in% hdr) c("cell", "malignant", "is_ref") else c("cell", "malignant"))
  if (!"is_ref" %in% names(d)) d[, is_ref := NA]
  d[, sample := s][]
}), fill = TRUE)
message(sprintf("[2] our per-cell calls loaded for %d samples, %d cells",
                uniqueN(calls$sample), nrow(calls)))

# barcode join: van Galen writes "<sample>_<barcode>"; strip any sample prefix on both sides
.core <- function(x) sub("^.*_", "", x)
geno[,  k := paste(sample, .core(cell_anno))]
calls[, k := paste(sample, .core(cell))]
m <- merge(geno, calls[, .(k, malignant, is_ref = as.logical(is_ref))], by = "k")
message(sprintf("[3] joined %d of %d genotyped cells (%.1f%%)",
                nrow(m), nrow(geno), 100 * nrow(m) / nrow(geno)))
if (!nrow(m)) stop("no barcode overlap -- check id conventions")

cat("\n================ inferCNV vs single-cell genotyping (GSE116256) ================\n")
tab <- m[, .(n = .N, called_malignant = sum(malignant == 1, na.rm = TRUE)), by = truth]
tab[, rate := round(called_malignant / n, 4)]
print(tab)

sens <- m[truth == "mut_pos", mean(malignant == 1, na.rm = TRUE)]
fp   <- m[truth == "wt_only", mean(malignant == 1, na.rm = TRUE)]
cat(sprintf("\n  SENSITIVITY on genotypically malignant cells : %.4f  (%d/%d)\n",
            sens, m[truth == "mut_pos", sum(malignant == 1, na.rm = TRUE)], m[truth == "mut_pos", .N]))
cat(sprintf("  flag rate on wt-only cells (UPPER bound FPR) : %.4f  (%d/%d)\n",
            fp, m[truth == "wt_only", sum(malignant == 1, na.rm = TRUE)], m[truth == "wt_only", .N]))

# THE DENOMINATORS ABOVE ARE PADDED, AND ASYMMETRICALLY SO. 41 ASSIGNS malignant=0 to autologous
# reference cells; the reference set IS the lymphoid population van Galen genotypes as wt-only, so
# contamination lands almost entirely on the wt-only side and pushes the flag rate DOWN -- i.e. it
# flatters the caller. Report the evaluable-cell version beside it.
if (m[!is.na(is_ref), .N] > 0) {
  e <- m[!(is_ref %in% TRUE)]
  cat(sprintf("\n  -- excluding autologous reference cells (%d of %d dropped; mut+ %d, wt-only %d) --\n",
              m[is_ref %in% TRUE, .N], nrow(m),
              m[is_ref %in% TRUE & truth == "mut_pos", .N], m[is_ref %in% TRUE & truth == "wt_only", .N]))
  cat(sprintf("  SENSITIVITY, evaluable cells only            : %.4f  (%d/%d)\n",
              e[truth == "mut_pos", mean(malignant == 1, na.rm = TRUE)],
              e[truth == "mut_pos", sum(malignant == 1, na.rm = TRUE)], e[truth == "mut_pos", .N]))
  cat(sprintf("  wt-only flag rate, evaluable cells only      : %.4f  (%d/%d)\n",
              e[truth == "wt_only", mean(malignant == 1, na.rm = TRUE)],
              e[truth == "wt_only", sum(malignant == 1, na.rm = TRUE)], e[truth == "wt_only", .N]))
  cat("  The second pair is the one comparable across callers; the first is what the pipeline acts on.\n")
} else cat("\n  [!] no is_ref column -- rerun 41; the rates above are over padded denominators.\n")

cat("\n-- per sample (samples with >=10 mut+ cells) --\n")
ps <- m[, .(mut_pos = sum(truth == "mut_pos"),
            sens = round(mean(malignant[truth == "mut_pos"] == 1, na.rm = TRUE), 3),
            wt_only = sum(truth == "wt_only"),
            wt_flag = round(mean(malignant[truth == "wt_only"] == 1, na.rm = TRUE), 3)), by = sample]
print(ps[mut_pos >= 10][order(-mut_pos)])

cat("\n-- sensitivity by van Galen cell type (mut+ cells only) --\n")
print(m[truth == "mut_pos", .(n = .N, sens = round(mean(malignant == 1, na.rm = TRUE), 3)),
        by = celltype_vg][order(-n)][1:12])
