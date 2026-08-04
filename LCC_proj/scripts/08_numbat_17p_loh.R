# 08_numbat_17p_loh.R ----
# Extract 17p LOH / deletion from the Numbat ALLELE-based segment consensus. This is the only
# evidence arm in the project that is orthogonal to expression-CNV: LOH is inferred from phased
# allele counts, not from smoothed expression, so it fails in different ways from inferCNV/CopyKAT.
#
# WHY 17p LOH SPECIFICALLY: biallelic TP53 inactivation in AML/MDS most often pairs one missense
# hit with loss of the other allele across 17p -- copy-neutral LOH or deletional. An expression
# caller sees the deletional form weakly and the copy-neutral form not at all; an allele caller
# sees both. Observed LLRs on the real calls here run 16-209, i.e. not marginal.
#
# INPUT  : LCC_NUMBAT_ROOT/<ds>/<sample>/numbat/segs_consensus_<N>.tsv   (highest N wins)
#          LCC_TAB_DIR/03_sample_manifest.csv via load_sample_meta()
# OUTPUT : LCC_TAB_DIR/08_numbat_17p.csv         per-sample 17p LOH call + evidence detail
#          LCC_TAB_DIR/08_numbat_segments17.csv  every non-neutral chr17 segment (audit trail)
# Usage  : Rscript LCC_proj/scripts/08_numbat_17p_loh.R
#
# COVERAGE LIMIT, stated up front: Numbat needs FASTQ/BAM, which only GSE289435 and GSE227903 have.
# This arm therefore covers 28 samples, not the cohort. It cannot be the primary caller -- it is the
# arm that tells us how much to trust the expression-only calls where it overlaps them.

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("LCC_proj", "scripts", "config_lcc.R"))

## -- Step 1. locate the winning iteration per sample ----
message("[1] locating Numbat segment consensus files under ", LCC_NUMBAT_ROOT)
fs <- list.files(LCC_NUMBAT_ROOT, pattern = "^segs_consensus_[0-9]+\\.tsv$",
                 recursive = TRUE, full.names = TRUE)
if (!length(fs)) stop("[08] no segs_consensus_*.tsv under ", LCC_NUMBAT_ROOT)
idx <- data.table(f = fs)[, `:=`(dir = dirname(f),
                                 iter = as.integer(sub(".*_([0-9]+)\\.tsv$", "\\1", f)))]
# [CODING_STANDARDS §9] "Numbat 迭代: 最高号 segs_consensus_{N}.tsv 为准"
setorder(idx, dir, -iter)
idx <- idx[, .SD[1], by = dir]
message("    ", nrow(idx), " samples with a segment consensus")

## -- Step 2. per-sample chr17 readout ----
message("[2] reading chr17 segments")
read_one <- function(f, dirpath) {
  parts <- strsplit(dirpath, "/", fixed = TRUE)[[1]]
  smp <- parts[length(parts) - 1]; ds <- parts[length(parts) - 2]   # <ds>/<sample>/numbat
  s <- fread(f, sep = "\t", showProgress = FALSE)
  s <- s[CHROM == 17]
  if (!nrow(s)) return(NULL)
  # Numbat's own first column is ALSO called `sample` and is INTEGER (a clone index, not a sample
  # name). Assigning our character sample id into it silently coerces to NA -- rename it first.
  if ("sample" %in% names(s)) setnames(s, "sample", "numbat_clone_idx")
  s[, `:=`(dataset = ds, sample = smp,
           overlaps_tp53 = seg_start <= LCC_TP53_END & seg_end >= LCC_TP53_START,
           # bp of this segment that fall on 17p (p arm ends at the centromere start)
           len_on_17p = pmax(0, pmin(seg_end, LCC_CHR17_CEN_START) - pmax(seg_start, 1L)))]
  s[]
}
segs <- rbindlist(lapply(seq_len(nrow(idx)), function(i) read_one(idx$f[i], idx$dir[i])), fill = TRUE)
fwrite_safe(segs[cnv_state != "neu"], file.path(LCC_TAB_DIR, "08_numbat_segments17.csv"))

## -- Step 3. per-sample call ----
message("[3] calling 17p LOH")
loh <- segs[cnv_state %in% LCC_NUMBAT_LOH_STATES & (is.na(LLR) | LLR >= LCC_NUMBAT_MIN_LLR)]
call <- segs[, .(n_seg_chr17 = .N, n_seg_nonneutral = sum(cnv_state != "neu")), by = .(dataset, sample)]
# UNION of covered bp, not the sum: Numbat can emit overlapping LOH and del segments over the same
# span, and summing them produced coverage fractions above 1.
union_bp <- function(start, end) {
  o <- order(start); a <- pmax(start[o], 1L); b <- pmin(end[o], LCC_CHR17_CEN_START)
  keep <- b > a; a <- a[keep]; b <- b[keep]
  if (!length(a)) return(0)
  tot <- 0; cs <- a[1]; ce <- b[1]
  for (i in seq_along(a)[-1]) {
    if (a[i] <= ce) { ce <- max(ce, b[i]) } else { tot <- tot + (ce - cs); cs <- a[i]; ce <- b[i] }
  }
  tot + (ce - cs)
}
agg <- loh[, .(arm17p_frac_loh = union_bp(seg_start, seg_end) / LCC_CHR17_CEN_START,
               arm17p_states   = paste(sort(unique(cnv_state)), collapse = "+"),
               best_llr        = max(LLR, na.rm = TRUE)), by = .(dataset, sample)]
tp53 <- segs[overlaps_tp53 == TRUE,
             .(tp53_seg_state = paste(sort(unique(cnv_state)), collapse = "+"),
               tp53_seg_llr   = max(LLR, na.rm = TRUE),
               tp53_p_loh     = max(p_loh, na.rm = TRUE)), by = .(dataset, sample)]
call <- merge(call, agg,  by = c("dataset", "sample"), all.x = TRUE)
call <- merge(call, tp53, by = c("dataset", "sample"), all.x = TRUE)
for (j in c("arm17p_frac_loh")) call[is.na(get(j)), (j) := 0]
call[is.na(tp53_seg_state), `:=`(tp53_seg_state = "neu", tp53_p_loh = NA_real_)]
# A sample is called only if a LOH/deletion segment covers at least half of 17p. The TP53-segment
# state is reported alongside but is NOT required: segment boundaries fall inside the TP53 gene in
# some samples, so requiring the gene-overlapping segment itself to be LOH would drop real calls.
call[, numbat_17p_loh := arm17p_frac_loh >= LCC_NUMBAT_MIN_ARM_FRAC]
man <- load_sample_meta(include_stroma_ref = TRUE)
call[man, `:=`(timepoint = i.timepoint, uid_patient = i.uid_patient,
               in_selection = i.in_selection), on = c("dataset", "sample")]
setorder(call, -arm17p_frac_loh)
fwrite_safe(call, file.path(LCC_TAB_DIR, "08_numbat_17p.csv"))

message("[4] result")
print(call[numbat_17p_loh == TRUE,
           .(dataset, sample, timepoint, arm17p_frac_loh = round(arm17p_frac_loh, 3),
             arm17p_states, tp53_seg_state, best_llr = round(best_llr, 1))])
# healthy donors are the negative control for this arm too, wherever Numbat covers any
nh <- call[timepoint == "Healthy"]
message(sprintf("    positives: %d / %d Numbat-covered samples; healthy negative control: %d / %d called",
                call[numbat_17p_loh == TRUE, .N], nrow(call),
                nh[numbat_17p_loh == TRUE, .N], nrow(nh)))
message("[done] 08_numbat_17p_loh")
