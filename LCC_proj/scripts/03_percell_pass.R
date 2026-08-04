# 03_percell_pass.R ----
# ONE pass over a sample's QC object that emits everything downstream needs from raw expression:
# gene detection rates, macrophage-vs-monocyte calls inside Mono_DC, three-level pseudobulk, and
# per-cell UCell pathway scores. Reading each RDS once (not four times) is the whole point.
#
# INPUT  : QC_RDS_DIR/<ds>/<sample>.rds                        (counts-only Seurat v5)
#          LCC_BMM_DIR/<ds>/<sample>__bmm_percell.csv          (hierarchy_bin, bmm_broad)
#          DIR_MALIGNANCY/<ds>/<sample>__consensus_percell.csv (malignant 0/1; absent for 90 samples)
#          LCC_PANEL_DIR/{fibrosis_ecm_panel,pathway_sets}.tsv + msigdb GMT
# OUTPUT : LCC_TAB_DIR/03_detect/<ds>__<sample>__detect.csv    per gene x stratum detection
#          LCC_TAB_DIR/03_myeloid/<ds>__<sample>__myeloid.csv  macrophage/monocyte counts (the CCC gate)
#          LCC_PERCELL_DIR/<ds>/<sample>__lcc_percell.csv.gz   per-cell bin/malignant/scores
#          LCC_PB_DIR/<ds>/<sample>__pseudobulk.rds            list(all, malignant, nonmalignant, by_bin)
# Usage  : Rscript LCC_proj/scripts/03_percell_pass.R --row=$SLURM_ARRAY_TASK_ID
#          Rscript LCC_proj/scripts/03_percell_pass.R --make_manifest      (build the row index first)
#
# WHY UCell and not AddModuleScore: UCell is rank-based per cell, so a score is comparable across
# cells of very different depth and across the 13 platforms in this cohort. AddModuleScore's
# control-bin subtraction is dataset-relative and would confound the study covariate we must model.
#
# WHY pseudobulk is raw summed counts, not averaged normalized values: the group comparison in 06
# uses limma-voom, which needs library sizes to compute precision weights.

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(here); library(Matrix)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--row", type = "integer", default = NA_integer_, help = "1-based row of the manifest"),
  make_option("--make_manifest", action = "store_true", default = FALSE),
  make_option("--workers", type = "integer", default = 1L, help = "UCell cores; match the sbatch alloc"),
  make_option("--force", action = "store_true", default = FALSE)
)))

MANIFEST  <- file.path(LCC_TAB_DIR, "03_sample_manifest.csv")
DETECT_DIR<- file.path(LCC_TAB_DIR, "03_detect")
MYELO_DIR <- file.path(LCC_TAB_DIR, "03_myeloid")

## -- Step 0. manifest ----
if (opt$make_manifest || !file.exists(MANIFEST)) {
  message("[0] building sample manifest from ", QC_RDS_DIR)
  ds <- list.dirs(QC_RDS_DIR, recursive = FALSE)
  man <- rbindlist(lapply(ds, function(d) {
    f <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
    if (!length(f)) return(NULL)
    data.table(dataset = basename(d), sample = sub("\\.rds$", "", basename(f)), rds = f)
  }))
  rm_ <- load_role_manifest()
  man[rm_, `:=`(timepoint = i.Timepoint, uid_patient = i.uid_patient, study = i.Study),
      on = c(dataset = "dataset", sample = "Sample")]
  setorder(man, dataset, sample)
  fwrite_safe(man, MANIFEST)
  message("    ", nrow(man), " samples -> ", MANIFEST)
  if (opt$make_manifest) quit(save = "no")
}
man <- fread(MANIFEST)
if (is.na(opt$row)) stop("[03] --row is required (or use --make_manifest)")
if (opt$row < 1 || opt$row > nrow(man)) stop("[03] --row out of range 1..", nrow(man))
DS <- man$dataset[opt$row]; SM <- man$sample[opt$row]

f_detect <- file.path(DETECT_DIR, sprintf("%s__%s__detect.csv", DS, SM))
f_myelo  <- file.path(MYELO_DIR,  sprintf("%s__%s__myeloid.csv", DS, SM))
f_cell   <- file.path(LCC_PERCELL_DIR, DS, sprintf("%s__lcc_percell.csv.gz", SM))
f_pb     <- file.path(LCC_PB_DIR, DS, sprintf("%s__pseudobulk.rds", SM))
if (!opt$force && all(file.exists(f_detect, f_myelo, f_cell, f_pb))) {
  message("[skip] ", DS, "/", SM, " -- all four outputs present"); quit(save = "no")
}
message("[1] ", DS, " / ", SM)

## -- Step 2. load counts + per-cell annotation ----
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS(man$rds[opt$row])
cts <- get_counts(obj)
rm(obj); invisible(gc())
message("    counts: ", nrow(cts), " genes x ", ncol(cts), " cells")

ann <- data.table(cell = colnames(cts))
f_bmm <- file.path(LCC_BMM_DIR, DS, paste0(SM, "__bmm_percell.csv"))
if (file.exists(f_bmm)) {
  b <- fread(f_bmm, select = c("cell", "bmm_broad", "hierarchy_bin", "high_error"))
  ann[b, `:=`(bmm_broad = i.bmm_broad, hierarchy_bin = i.hierarchy_bin, high_error = i.high_error), on = "cell"]
}
f_con <- file.path(DIR_MALIGNANCY, DS, paste0(SM, "__consensus_percell.csv"))
if (file.exists(f_con)) {
  cn <- fread(f_con, select = c("cell", "malignant"))
  ann[cn, malignant := i.malignant, on = "cell"]
}
if (!"hierarchy_bin" %in% names(ann)) ann[, hierarchy_bin := NA_character_]
if (!"malignant"     %in% names(ann)) ann[, malignant := NA_integer_]
ann[is.na(hierarchy_bin) | hierarchy_bin == "", hierarchy_bin := "unassigned"]

## -- Step 3. normalize once (log1p CP10K), reused by detection means and UCell ----
lib <- Matrix::colSums(cts)
norm <- cts
norm@x <- norm@x / rep.int(pmax(lib, 1), diff(norm@p))   # column-wise CP1 on a dgCMatrix
norm@x <- log1p(norm@x * 1e4)

## -- Step 4. panel gene detection, per stratum ----
message("[2] gene detection")
panel <- load_gene_panel()
want  <- unique(c(panel$gene, panel$alias[nzchar(panel$alias)], LCC_NECTIN_GENES))
have  <- intersect(want, rownames(cts))
missing_from_ref <- setdiff(want, rownames(cts))

detect_block <- function(cells, stratum, level) {
  if (!length(cells) || !length(have)) return(NULL)
  sc <- cts[have, cells, drop = FALSE]; sn <- norm[have, cells, drop = FALSE]
  data.table(gene = have, stratum = stratum, level = level, n_cells = length(cells),
             n_nonzero  = Matrix::rowSums(sc > 0),
             mean_lognorm = Matrix::rowMeans(sn),
             sum_counts = Matrix::rowSums(sc))
}
blocks <- list(detect_block(ann$cell, "all", "sample"))
for (b in setdiff(unique(ann$hierarchy_bin), NA)) blocks <- c(blocks, list(detect_block(ann[hierarchy_bin == b]$cell, b, "bin")))
for (mv in c(0, 1)) blocks <- c(blocks, list(detect_block(ann[malignant == mv]$cell,
                                                          c("nonmalignant", "malignant")[mv + 1], "malignancy")))
det <- rbindlist(blocks[!vapply(blocks, is.null, TRUE)])
det[, `:=`(pct_nonzero = round(100 * n_nonzero / n_cells, 4), dataset = DS, sample = SM)]
# genes absent from the reference are reported explicitly -- a zero detection rate and a missing
# gene are different facts and must not collapse into the same row.
if (length(missing_from_ref))
  det <- rbind(det, data.table(gene = missing_from_ref, stratum = "all", level = "sample",
                               n_cells = ncol(cts), n_nonzero = NA_integer_, mean_lognorm = NA_real_,
                               sum_counts = NA_real_, pct_nonzero = NA_real_, dataset = DS, sample = SM),
               fill = TRUE)
fwrite_safe(det, f_detect)

## -- Step 5. macrophage vs monocyte inside Mono_DC (the CCC node-vocabulary gate) ----
message("[3] myeloid scoring")
score_set <- function(genes) {
  g <- intersect(genes, rownames(norm))
  if (!length(g)) return(rep(NA_real_, ncol(norm)))
  Matrix::colMeans(norm[g, , drop = FALSE])
}
ann[, `:=`(mac_score = score_set(LCC_MAC_MARKERS), mono_score = score_set(LCC_MONO_MARKERS))]
# Relative call inside Mono_DC only. A within-sample median split would force a 50/50 answer, so use
# the absolute contrast instead and report the score distributions so the gate can be re-cut later.
ann[, myeloid_class := NA_character_]
ann[hierarchy_bin == "Mono_DC" & !is.na(mac_score) & !is.na(mono_score),
    myeloid_class := fifelse(mac_score > mono_score, "macrophage_like",
                      fifelse(mono_score > mac_score, "monocyte_like", "ambiguous"))]
myelo <- data.table(dataset = DS, sample = SM,
                    n_cells_total = ncol(cts),
                    n_mono_dc = ann[hierarchy_bin == "Mono_DC", .N],
                    n_macrophage_like = ann[myeloid_class == "macrophage_like", .N],
                    n_monocyte_like   = ann[myeloid_class == "monocyte_like", .N],
                    mac_score_p50 = as.numeric(median(ann[hierarchy_bin == "Mono_DC"]$mac_score, na.rm = TRUE)),
                    mac_score_p90 = as.numeric(quantile(ann[hierarchy_bin == "Mono_DC"]$mac_score, .9, na.rm = TRUE)),
                    n_mac_markers_present = length(intersect(LCC_MAC_MARKERS, rownames(cts))))
fwrite_safe(myelo, f_myelo)

## -- Step 6. UCell pathway scores ----
message("[4] UCell")
suppressPackageStartupMessages(library(UCell))
psets <- fread_commented(LCC_PATHWAY_TSV)
sets  <- read_gmt_subset(psets$set_name)
set.seed(SEED)
uc <- UCell::ScoreSignatures_UCell(norm, features = sets, ncores = opt$workers, maxRank = 1500)
ann <- cbind(ann, as.data.table(uc))
fwrite_safe(ann, f_cell)

## -- Step 7. three-level pseudobulk (raw summed counts) ----
message("[5] pseudobulk")
pb <- function(cells) if (length(cells)) Matrix::rowSums(cts[, cells, drop = FALSE]) else NULL
by_bin <- lapply(setdiff(unique(ann$hierarchy_bin), NA), function(b) pb(ann[hierarchy_bin == b]$cell))
names(by_bin) <- setdiff(unique(ann$hierarchy_bin), NA)
res <- list(dataset = DS, sample = SM, genes = rownames(cts),
            n_cells = c(all = ncol(cts),
                        malignant = ann[malignant == 1L, .N], nonmalignant = ann[malignant == 0L, .N],
                        vapply(by_bin, function(x) 0L, 0L)),
            n_cells_bin = ann[, .N, by = hierarchy_bin],
            all = pb(ann$cell), malignant = pb(ann[malignant == 1L]$cell),
            nonmalignant = pb(ann[malignant == 0L]$cell), by_bin = by_bin)
dir.create(dirname(f_pb), recursive = TRUE, showWarnings = FALSE)
saveRDS(res, f_pb, compress = TRUE)
message("[done] ", DS, "/", SM)
