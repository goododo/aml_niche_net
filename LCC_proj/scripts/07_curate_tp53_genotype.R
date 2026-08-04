# 07_curate_tp53_genotype.R ----
# Build the TP53 GENOTYPE table -- true mutation calls, not CNV proxies -- and use it to measure how
# well the 17p-loss CNV rule from 02 actually performs.
#
# WHY THIS SCRIPT EXISTS: the project inventory recorded exactly ONE TP53-mutant sample in the whole
# cohort (Petti2019 809653, TP53 E286G). That was an undercount. van Galen 2019 (GSE116256) deposited
# PER-CELL single-cell genotyping in its GEO `*.anno.txt.gz` files -- columns MutTranscripts /
# WtTranscripts carry strings like "TP53.R273L/12" -- and those files have been sitting in the raw
# deposit all along. Parsing them yields THREE genotype-confirmed TP53-mutant patients with
# single-cell resolution, which is a different class of evidence from anything inferCNV can give.
#
# LANGUAGE DISCIPLINE [CODING_STANDARDS §9]: only calls sourced here may be labelled "TP53-mutant".
# Anything derived from inferCNV/CopyKAT stays "17p-loss" / "TP53-locus deleted". The two are joined
# in this script precisely so the difference stays visible.
#
# INPUT  : RAW/public/GEO/GSE116256_RAW/*.anno.txt.gz     (van Galen per-cell genotyping)
#          LCC_PANEL_DIR/tp53_genotype_literature.tsv     (hand-curated from papers; provenance per row)
#          LCC_TAB_DIR/02_sample_cnv_proxy.csv            (the 17p-loss CNV calls to be validated)
#          LCC_TAB_DIR/03_sample_manifest.csv via load_sample_meta()
# OUTPUT : LCC_TAB_DIR/07_tp53_genotype_percell.csv       per-cell TP53 mut/wt calls (van Galen only)
#          LCC_TAB_DIR/07_tp53_genotype_sample.csv        per-sample TP53 status + evidence tier
#          LCC_TAB_DIR/07_proxy_vs_genotype.csv           confusion matrix: 17p rule vs genotype truth
# Usage  : Rscript LCC_proj/scripts/07_curate_tp53_genotype.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("LCC_proj", "scripts", "config_lcc.R"))

VG_RAW <- file.path(GEO_RAW_DIR, "GSE116256_RAW")

## -- Step 1. parse van Galen per-cell genotyping ----
message("[1] parsing van Galen per-cell genotyping from ", VG_RAW)
# Format: MutTranscripts / WtTranscripts hold "GENE.VARIANT/UMI" entries, multiple entries separated
# by "; ". A cell is TP53-mutant if any MutTranscripts entry starts with "TP53.", TP53-WT if any
# WtTranscripts entry does and no mutant entry does, and uninformative if TP53 was not covered.
anno <- list.files(VG_RAW, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE)
if (!length(anno)) stop("[07] no van Galen anno files under ", VG_RAW)

parse_anno <- function(f) {
  s <- sub("\\.anno\\.txt\\.gz$", "", sub("^GSM[0-9]+_", "", basename(f)))
  d <- fread(f, sep = "\t", colClasses = "character", showProgress = FALSE)
  if (!all(c("Cell", "MutTranscripts", "WtTranscripts") %in% names(d))) return(NULL)
  has_tp53 <- function(x) grepl("(^|[;, ])TP53\\.", x)
  # a "wt" entry is spelled TP53.wt; a mutant entry is TP53.<aa-change>
  d[, `:=`(tp53_mut_entry = has_tp53(MutTranscripts) & grepl("TP53\\.(?!wt)", MutTranscripts, perl = TRUE),
           tp53_wt_entry  = grepl("TP53\\.wt", WtTranscripts))]
  d[, tp53_cell_call := fifelse(tp53_mut_entry, "mutant",
                         fifelse(tp53_wt_entry, "wildtype", NA_character_))]
  out <- d[, .(sample = s, cell = Cell, cell_type = if ("CellType" %in% names(d)) CellType else NA_character_,
               mut_transcripts = MutTranscripts, wt_transcripts = WtTranscripts, tp53_cell_call)]
  out[!is.na(tp53_cell_call) | grepl("TP53", mut_transcripts) | grepl("TP53", wt_transcripts)]
}
pc <- rbindlist(lapply(anno, parse_anno), fill = TRUE)
pc[, dataset := "GSE116256"]
# the exact variant string, for the provenance column
pc[, tp53_variant := {
  v <- regmatches(mut_transcripts, regexpr("TP53\\.[A-Za-z0-9*_]+", mut_transcripts))
  ifelse(lengths(regmatches(mut_transcripts, gregexpr("TP53\\.[A-Za-z0-9*_]+", mut_transcripts))) > 0,
         vapply(regmatches(mut_transcripts, gregexpr("TP53\\.[A-Za-z0-9*_]+", mut_transcripts)),
                function(x) paste(unique(x), collapse = "+"), ""), NA_character_)
}]
fwrite_safe(pc, file.path(LCC_TAB_DIR, "07_tp53_genotype_percell.csv"))
message("    ", nrow(pc), " cells carry a TP53 genotype read across ",
        uniqueN(pc[!is.na(tp53_cell_call)]$sample), " samples")

## -- Step 2. roll per-cell genotype up to per-sample status ----
message("[2] per-sample TP53 status from single-cell genotyping")
vg <- pc[, .(n_cells_genotyped = sum(!is.na(tp53_cell_call)),
             n_cells_tp53_mut  = sum(tp53_cell_call == "mutant", na.rm = TRUE),
             n_cells_tp53_wt   = sum(tp53_cell_call == "wildtype", na.rm = TRUE),
             tp53_variants     = paste(sort(unique(na.omit(tp53_variant))), collapse = "+")),
         by = .(dataset, sample)]
vg <- vg[n_cells_genotyped > 0 | nzchar(tp53_variants)]
vg[, tp53_status := fifelse(n_cells_tp53_mut > 0, "mutant",
                     fifelse(n_cells_tp53_wt > 0, "wildtype", "uninformative"))]
vg[, `:=`(evidence = "single-cell genotyping (van Galen 2019 GEO anno)", evidence_tier = "A_genotype")]

## -- Step 3. literature-curated calls (hand table; one row per sample or patient) ----
lit_f <- file.path(LCC_PANEL_DIR, "tp53_genotype_literature.tsv")
lit <- if (file.exists(lit_f)) fread_commented(lit_f) else
  data.table(dataset = character(), sample = character(), patient = character(),
             tp53_status = character(), tp53_variants = character(),
             evidence = character(), evidence_tier = character(), citation = character())
message("[3] literature table: ", nrow(lit), " curated rows")

geno <- rbind(vg[, .(dataset, sample, tp53_status, tp53_variants, n_cells_genotyped,
                     n_cells_tp53_mut, n_cells_tp53_wt, evidence, evidence_tier)],
              lit[, .(dataset, sample, tp53_status, tp53_variants, n_cells_genotyped = NA_integer_,
                      n_cells_tp53_mut = NA_integer_, n_cells_tp53_wt = NA_integer_,
                      evidence, evidence_tier)], fill = TRUE)
man <- load_sample_meta(include_stroma_ref = TRUE)
geno[man, `:=`(timepoint = i.timepoint, uid_patient = i.uid_patient,
               in_selection = i.in_selection), on = c("dataset", "sample")]
# patient-level status: a mutation found at ANY timepoint applies to the patient
geno[!is.na(uid_patient), tp53_status_patient :=
       fifelse(any(tp53_status == "mutant"), "mutant",
        fifelse(any(tp53_status == "wildtype"), "wildtype", "uninformative")), by = uid_patient]
setorder(geno, -tp53_status, dataset, sample)
fwrite_safe(geno, file.path(LCC_TAB_DIR, "07_tp53_genotype_sample.csv"))
message("    genotype-confirmed TP53-mutant samples: ", geno[tp53_status == "mutant", .N],
        " across ", uniqueN(geno[tp53_status == "mutant"]$uid_patient), " patients")
print(geno[tp53_status == "mutant",
           .(dataset, sample, timepoint, tp53_variants, n_cells_tp53_mut, n_cells_tp53_wt, in_selection)])

## -- Step 4. the 17p CNV rule against genotype truth ----
message("[4] 17p-loss CNV rule vs genotype truth")
# prefer the all-datasets 17p calls: the genotype truth set is mostly GSE116256, outside the selection
f <- file.path(LCC_TAB_DIR, "02_sample_cnv_proxy_all_datasets.csv")
if (!file.exists(f)) f <- file.path(LCC_TAB_DIR, "02_sample_cnv_proxy.csv")
if (!file.exists(f)) {
  message("    [skip] 02 output absent -- run 02_define_cnv_proxy.R (use --all_datasets to cover GSE116256)")
} else {
  cnv <- fread(f)[, .(dataset, sample, p17_pos, ck_max_arms)]
  cm <- merge(geno[tp53_status %in% c("mutant", "wildtype"),
                   .(dataset, sample, tp53_status, tp53_variants, evidence_tier)],
              cnv, by = c("dataset", "sample"))
  if (!nrow(cm)) {
    message("    [skip] no sample has BOTH a genotype call and a 17p call. ",
            "GSE116256 is outside LCC_CORE_DATASETS, so 02 must be re-run with --all_datasets.")
  } else {
    tab <- cm[, .N, by = .(tp53_status, p17_pos)]
    tp <- sum(cm$tp53_status == "mutant"   & cm$p17_pos); fn <- sum(cm$tp53_status == "mutant"   & !cm$p17_pos)
    fp <- sum(cm$tp53_status == "wildtype" & cm$p17_pos); tn <- sum(cm$tp53_status == "wildtype" & !cm$p17_pos)
    perf <- data.table(n_evaluable = nrow(cm), true_pos = tp, false_neg = fn, false_pos = fp, true_neg = tn,
                       sensitivity = if (tp + fn) tp / (tp + fn) else NA_real_,
                       specificity = if (fp + tn) tn / (fp + tn) else NA_real_,
                       ppv         = if (tp + fp) tp / (tp + fp) else NA_real_)
    fwrite_safe(cbind(perf, cm[, .(detail = paste(sprintf("%s:%s/17p=%s", sample, tp53_status, p17_pos),
                                                  collapse = "; "))]),
                file.path(LCC_TAB_DIR, "07_proxy_vs_genotype.csv"))
    print(tab); print(perf)
  }
}
message("[done] 07_curate_tp53_genotype")
