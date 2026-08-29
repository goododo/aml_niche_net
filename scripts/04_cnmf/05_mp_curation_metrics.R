#!/usr/bin/env Rscript
# 05_mp_curation_metrics.R ----
# INPUT  : CNMF_RES_DIR/malignant/mp_usage/<ds>/<sample>__mp_usage.csv           (malignant cells only)
#          CNMF_RES_DIR/malignant/mp_usage_all_bins/<ds>/<sample>__mp_usage.csv  (all cells in the 7 CCC bins)
#          CNMF_TAB_DIR/malignant_mp_annotation.csv   (per MP: class, max_jaccard, sampleCoverage, silhouette)
#          CNMF_TAB_DIR/malignant_metaprograms.tsv    (MP -> genes; comma-delimited despite the extension)
#          CNMF_SCRIPT_DIR/mp_labels.tsv              (the hand-typed curation, the decision of record)
# OUTPUT : CNMF_TAB_DIR/mp_curation_metrics.csv       (one row per MP x cell_set, every cited number computed)
# WHAT   : Computes, from data, each quantity the hand curation cites, so the curation can be checked
#          instead of taken on trust -- and computes it on BOTH cell sets, because they disagree.
#
# WHY THIS EXISTS. scripts/04_cnmf/mp_labels.tsv decides which meta-programs are trustworthy, and it
# is hand-typed: no script writes it, and its `note` column carries numbers ("91% Mono_DC bin") that
# nothing recomputes. Two consequences. First, the headline "3 of 10 meta-programs are robust tumour
# programs" is not reproducible from code -- the code-derived class in malignant_mp_annotation.csv
# says 5 (MP1, MP5, MP6, MP8, MP10), and the curation overrode it on MP1 and MP8. Second, and worse,
# the numbers it cites were computed on the MALIGNANT-ONLY scoring, which 03_score_metaprograms.R now
# also produces for all cells in the CCC bins -- and the two do not agree.
#
# THE CIRCULARITY THIS IS BUILT TO EXPOSE. On the malignant-only cell set, MP5/MP6/MP10 sit at
# 93/96/97% Mono_DC, which is what "robust_tumor, monocytic" rests on. But inferCNV calls 24.5% of
# Mono_DC cells malignant and only 1.6% of HSC_MPP cells, so the malignant compartment is itself
# Mono_DC-biased. Measuring bin specificity inside it partly restates where the caller finds
# malignant cells, rather than where the program lives. Reporting both cell sets side by side is not
# indecision: it is the honest form of a question that has not been answered yet.
#
# THIS SCRIPT DOES NOT DECIDE. It computes and it flags disagreement. Curation judgements that are
# not mechanisable ("this is technical immediate-early-gene noise") stay in mp_labels.tsv, which
# remains the decision of record -- but every number that justifies one becomes checkable.
#
# Usage : Rscript scripts/04_cnmf/05_mp_curation_metrics.R [--support_min_frac 0.01]

suppressPackageStartupMessages({ library(data.table); library(optparse); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_cnmf.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--support_min_frac", type = "double", default = 0.01,
              help = "a dataset/sample SUPPORTS an MP if the MP is top_MP for at least this fraction of its cells")
)))

## -- Step 1. the two cell sets ----
CELL_SETS <- list(
  malignant    = file.path(CNMF_RES_DIR, "malignant", "mp_usage"),
  all_ccc_bins = file.path(CNMF_RES_DIR, "malignant", "mp_usage_all_bins")
)

read_usage <- function(dir) {
  f <- list.files(dir, pattern = "__mp_usage\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(NULL)
  rbindlist(lapply(f, fread, select = c("cell", "dataset", "sample", "hierarchy_bin", "top_MP")), fill = TRUE)
}

## -- Step 2. the metrics, computed per MP within one cell set ----
# HSI (hierarchy specificity index): of the cells an MP WINS (is their top_MP), what share fall in
# its single most common bin. This is the "91% Mono_DC" the curation hand-computed. It is a property
# of the winning population, not of the score distribution, which is why an MP that never wins has
# no HSI rather than a zero.
metrics_for <- function(U, cell_set) {
  if (is.null(U)) return(NULL)
  won <- U[!is.na(top_MP)]
  bin <- won[, .N, by = .(top_MP, hierarchy_bin)]
  bin[, share := N / sum(N), by = top_MP]
  hsi <- bin[bin[, .I[which.max(share)], by = top_MP]$V1,
             .(top_MP, HSI = round(share, 4), top_bin = hierarchy_bin)]

  # SUPPORT COUNTED BY DATASET, NOT ONLY BY SAMPLE. sampleCoverage in 02's annotation counts samples,
  # so one large cohort can carry a program on its own -- GSE185381 contributes 52 of 212 samples.
  # Counting independent datasets is the protection the two-level reproducibility design exists for
  # (see P8 in scripts/DECISIONS_pending.md); this is the cheap half of it, not SRRS itself.
  tot_ds <- won[, .(n_ds = .N), by = dataset]
  by_ds  <- won[, .N, by = .(top_MP, dataset)][tot_ds, on = "dataset"][, frac := N / n_ds]
  sup_ds <- by_ds[frac >= opt$support_min_frac, .(n_datasets_supporting = .N), by = top_MP]

  tot_s <- won[, .(n_s = .N), by = .(dataset, sample)]
  by_s  <- won[, .N, by = .(top_MP, dataset, sample)][tot_s, on = c("dataset","sample")][, frac := N / n_s]
  sup_s <- by_s[frac >= opt$support_min_frac, .(n_samples_supporting = .N), by = top_MP]

  n_won <- won[, .(n_cells_won = .N, frac_cells_won = .N / nrow(won)), by = top_MP]

  M <- Reduce(function(a, b) merge(a, b, by = "top_MP", all = TRUE), list(hsi, sup_ds, sup_s, n_won))
  setnames(M, "top_MP", "MP")
  M[, `:=`(cell_set = cell_set,
           n_cells_total = nrow(won),
           n_datasets_total = uniqueN(won$dataset),
           n_samples_total  = uniqueN(won[, paste(dataset, sample)]))]
  M[]
}

message("[1] reading both cell sets")
RES <- rbindlist(lapply(names(CELL_SETS), function(cs) {
  U <- read_usage(CELL_SETS[[cs]])
  if (is.null(U)) { message("    [", cs, "] no files at ", CELL_SETS[[cs]], " -- skipped"); return(NULL) }
  message(sprintf("    [%s] %d cells, %d samples, %d datasets",
                  cs, nrow(U), uniqueN(U[, paste(dataset, sample)]), uniqueN(U$dataset)))
  metrics_for(U, cs)
}), fill = TRUE)
stopifnot(nrow(RES) > 0)

## -- Step 3. join the static per-MP quantities ----
ann  <- fread(file.path(CNMF_TAB_DIR, "malignant_mp_annotation.csv"))
setnames(ann, "malignant_MP", "MP")
# .tsv by extension, comma-delimited by content. fread sniffs the separator, but naming it here means
# a future switch to real tabs fails loudly instead of producing one column.
sig  <- fread(file.path(CNMF_TAB_DIR, "malignant_metaprograms.tsv"), sep = ",")
stopifnot(all(c("MP", "gene") %in% names(sig)))
ngene <- sig[, .(n_genes = uniqueN(gene)), by = MP]

# FULL GRID FIRST, THEN JOIN. An MP declared by 02 but never the top program for a single cell is
# invisible to every metric above and would simply be absent -- MP1 is the standing example: 1 gene,
# silhouette -0.031, and it wins nowhere in either cell set. Building the complete MP x cell_set grid
# and left-joining onto it means "never wins" shows up as a row with n_cells_won = 0, not as a
# missing row, and the static per-MP columns are attached exactly once for every MP.
all_mp <- sort(unique(c(ann$MP, sig$MP)))
grid <- CJ(MP = all_mp, cell_set = unique(RES$cell_set), sorted = FALSE)
RES <- merge(grid, RES, by = c("MP", "cell_set"), all.x = TRUE)
RES[is.na(n_cells_won), `:=`(n_cells_won = 0L, frac_cells_won = 0)]
RES[, ever_top_MP := n_cells_won > 0]

RES <- merge(RES, ann[, .(MP, code_class = class, max_jaccard, sampleCoverage, silhouette)],
             by = "MP", all.x = TRUE)
RES <- merge(RES, ngene, by = "MP", all.x = TRUE)
# Support counts are 0, not NA, for an MP that never wins: it is supported by no dataset, which is a
# measurement, whereas NA would read as "not measured".
for (cc in c("n_datasets_supporting", "n_samples_supporting"))
  RES[is.na(get(cc)), (cc) := 0L]
setorder(RES, cell_set, -HSI, na.last = TRUE)

## -- Step 4. the hand curation, and where it disagrees ----
lab_f <- file.path(CNMF_SCRIPT_DIR, "mp_labels.tsv")
lab <- if (file.exists(lab_f)) fread(lab_f) else NULL
if (!is.null(lab)) {
  RES <- merge(RES, lab[, .(MP, hand_confidence = confidence, hand_label = biological_label)],
               by = "MP", all.x = TRUE)
}

out <- file.path(CNMF_TAB_DIR, "mp_curation_metrics.csv")
fwrite_safe(RES, out)

message("\n[2] PER-MP METRICS, both cell sets side by side")
cols <- c("MP","cell_set","HSI","top_bin","n_cells_won","n_datasets_supporting","n_samples_supporting",
          "ever_top_MP","n_genes","silhouette","sampleCoverage","code_class","hand_confidence")
cols <- intersect(cols, names(RES))
print(RES[order(MP, cell_set), ..cols])

## -- Step 5. the two things a reader must not miss ----
if (uniqueN(RES$cell_set) == 2) {
  # Keyed on MP ALONE, not MP + top_bin: MP7 and MP9 change which bin they dominate between the two
  # cell sets, and keying on the bin would split them into two half-empty rows and hide exactly that.
  W <- dcast(RES[!is.na(HSI)], MP ~ cell_set, value.var = c("HSI", "top_bin"))
  if (all(c("HSI_malignant", "HSI_all_ccc_bins") %in% names(W))) {
    setnames(W, c("HSI_malignant","HSI_all_ccc_bins","top_bin_malignant","top_bin_all_ccc_bins"),
             c("malignant","all_ccc_bins","bin_mal","bin_all"))
    W[, delta := round(all_ccc_bins - malignant, 3)]
    W[, bin_changed := bin_mal != bin_all]
    message("\n[3] HSI ON THE TWO CELL SETS. A large negative delta means the program's bin specificity")
    message("    was substantially a property of WHERE THE CALLER FINDS MALIGNANT CELLS, not of the")
    message("    program: inferCNV calls 24.5% of Mono_DC malignant and 1.6% of HSC_MPP.")
    print(W[order(delta)])
  }
}

if (!is.null(lab)) {
  D <- RES[cell_set == "all_ccc_bins" & !is.na(hand_confidence)]
  D[, code_says_tumor := code_class == "tumor_specific"]
  D[, hand_says_tumor := hand_confidence == "robust_tumor"]
  dis <- D[code_says_tumor != hand_says_tumor,
           .(MP, code_class, hand_confidence, HSI, ever_top_MP, n_genes, silhouette)]
  message("\n[4] WHERE THE HAND CURATION OVERRODE THE CODE-DERIVED CLASS")
  if (nrow(dis)) {
    print(dis)
    message("    Each override needs a computed column that supports it, or it is an unrecorded judgement.")
  } else message("    (none -- the two agree on every MP)")
  message(sprintf("\n    code-derived tumor_specific: %d | hand robust_tumor: %d",
                  sum(D$code_says_tumor, na.rm = TRUE), sum(D$hand_says_tumor, na.rm = TRUE)))
}

message("\n[done] wrote ", out)
