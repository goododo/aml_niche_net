#!/usr/bin/env Rscript
# 10_stromal_denovo.R ----
# Find stromal cells DE NOVO from raw counts, bypassing the BoneMarrowMap projection label.
#
# WHY THIS SCRIPT EXISTS: the projected `Stromal` bin is not trustworthy in either direction.
#   FALSE POSITIVES -- GSE227903 1886_R has 528 cells labelled Stromal with DCN 0%, LUM 0.2%,
#     CXCL12 0.2%. GSE185991 M95/M93 the same. BoneMarrowMap's reference is haematopoietic-centric,
#     so low-confidence cells land in the nearest reference bin and that bin is often `Stromal`.
#   FALSE NEGATIVES -- the converse is equally possible: a real MSC has no good home in a
#     haematopoietic reference and can be forced into an HSPC bin. A per-gene table scan found
#     DCN+/LUM+/CXCL12+ cells in samples whose projection reports no stroma at all.
# Neither failure can be settled by per-gene counts, because a per-gene count cannot tell one cell
# co-expressing DCN+LUM+COL1A1 from three cells each carrying one transcript of ambient RNA. This
# script therefore works PER CELL and calls stroma only on module CO-EXPRESSION.
#
# THE CALIBRATION, which is the point: CD34/CD117-SORTED libraries (GSE185991, GSE147989,
# Chen2023 *_CD34) contain no stroma BY CONSTRUCTION -- the sort removed it. Any de-novo call in a
# sorted library is a false positive. So the sorted libraries fix the threshold, exactly as the
# healthy donors fixed the 17p arm-fraction threshold in 02. Stroma-enriched libraries
# (Chen2023 *_Niche_Immune, GSE253355) are the positive control that the threshold must not kill.
#
# INPUT  : LCC_TAB_DIR/04_detection_by_sample.csv     (the cheap all-sample screen)
#          QC_RDS_DIR/<ds>/<sample>.rds               (raw counts; screened samples only)
#          LCC_PERCELL_DIR/<ds>/<sample>__lcc_percell.csv.gz  (hierarchy_bin, for the cross-tab)
# OUTPUT : LCC_TAB_DIR/10_stromal_screen.csv          every sample, screened in or out, with counts
#          LCC_TAB_DIR/10_stromal_calibration.csv     threshold grid x sorted-library FPR
#          LCC_TAB_DIR/10_stromal_denovo_sample.csv   per-sample de-novo stromal counts + subtype
#          LCC_TAB_DIR/10_stromal_vs_projection.csv   de-novo call vs BoneMarrowMap `Stromal` bin
# Usage  : Rscript LCC_proj/scripts/10_stromal_denovo.R [--min_screen 20]

# Deliberately does NOT load Seurat. The user library carries Seurat 5.5.1, which needs
# promises >= 1.5.0, while the env has 1.3.3 -- loading it aborts. 03 hit the same wall and solved it
# the same way: utils.R::get_counts() reaches the counts through SeuratObject alone, which is all
# that deserialising a QC object requires.
suppressPackageStartupMessages({
  library(data.table); library(here); library(Matrix); library(optparse)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
# SORTED_LIBRARY_DATASETS is the main line's authoritative list of CD34/CD117-sorted libraries and it
# is what makes the calibration in Step 3 valid. Sourced, never copied: a sorted library that gets
# added upstream must not silently stay in the whole_MNC class here.
source(here::here("scripts", "config", "config_malignancy.R"))
set.seed(SEED)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--min_screen", type = "integer", default = 20L,
              help = "a sample enters the per-cell stage only if some stromal marker is detected in >= this many cells [20]")
)))

## -- Step 0. cell-identity modules ----
# Chosen so that each module needs several genes to fire. Single-gene calls are what produced the
# artefacts this script exists to remove, so no module is allowed to be decided by one gene.
# MOD / SCREEN_GENES / the subtype rule now live in stromal_modules.R so that 16_dsp_handoff.R can
# reproduce exactly these per-cell calls when it builds the SpatialDecon stromal profiles. The
# definitions are unchanged -- this was a move, not an edit.
source(here::here("LCC_proj", "scripts", "stromal_modules.R"))

## -- Step 1. cheap all-sample screen from the 04 table ----
message("[1] screening all samples from 04_detection_by_sample.csv")
det <- fread(file.path(LCC_TAB_DIR, "04_detection_by_sample.csv"))
scr <- dcast(det[stratum == "all" & gene %in% SCREEN_GENES],
             dataset + sample + n_cells ~ gene, value.var = "n_nonzero", fill = 0L)
for (g in SCREEN_GENES) if (!g %in% names(scr)) scr[, (g) := 0L]
scr[, screen_max := pmax(DCN, LUM, COL1A1, CXCL12)]

man <- load_sample_meta(include_stroma_ref = TRUE)
scr[man, `:=`(timepoint = i.timepoint, uid_patient = i.uid_patient), on = c("dataset", "sample")]
# library type drives the whole interpretation, so it is an explicit column, not a comment
scr[, library_type := fifelse(dataset %in% SORTED_LIBRARY_DATASETS | grepl("_CD34$", sample), "sorted_CD34",
                       fifelse(grepl("(?i)niche", sample) | dataset %in% LCC_STROMA_REF_DATASETS,
                               "stroma_enriched", "whole_MNC"))]
scr[, screened_in := screen_max >= opt$min_screen]
setorder(scr, -screen_max)
fwrite_safe(scr, file.path(LCC_TAB_DIR, "10_stromal_screen.csv"))
# NO SILENT CAP: state exactly what the screen drops and on what evidence.
message(sprintf("    screened IN %d / %d samples (>= %d cells positive for DCN/LUM/COL1A1/CXCL12)",
                sum(scr$screened_in), nrow(scr), opt$min_screen))
message(sprintf("    screened OUT %d samples; their max marker+ count is %d cells, i.e. below any",
                sum(!scr$screened_in), max(scr[screened_in == FALSE]$screen_max)))
message("    usable stromal population -- they are excluded from the per-cell stage, not from the report")

## -- Step 2. per-cell module detection on the screened-in samples + both controls ----
# Controls are forced in even if the screen would drop them: the sorted libraries ARE the FPR
# measurement, so silently omitting a sorted library with a low screen count would flatter the
# threshold.
todo <- rbind(scr[screened_in == TRUE], scr[library_type == "sorted_CD34"])[!duplicated(paste(dataset, sample))]
message("[2] per-cell pass over ", nrow(todo), " samples (screened-in + all sorted controls)")

cell_calls <- function(ds, sm) {
  f <- file.path(QC_RDS_DIR, ds, paste0(sm, ".rds"))
  if (!file.exists(f)) return(NULL)
  cnt <- get_counts(readRDS(f))   # utils.R; same accessor 03 used, so the cell sets match exactly
  gn <- rownames(cnt)
  # detection only (count > 0): no normalisation needed, and a binary readout is far more robust to
  # depth differences across the 13 studies than any expression cutoff would be.
  nd <- lapply(MOD, function(gs) {
    idx <- match(intersect(gs, gn), gn)
    if (!length(idx)) return(rep(0L, ncol(cnt)))
    as.integer(Matrix::colSums(cnt[idx, , drop = FALSE] > 0))
  })
  d <- as.data.table(nd)
  d[, `:=`(dataset = ds, sample = sm, cell = colnames(cnt),
           n_genes = as.integer(Matrix::colSums(cnt > 0)))]
  d[]
}

pc <- rbindlist(lapply(seq_len(nrow(todo)), function(i) {
  if (i %% 10 == 0) message("    ", i, "/", nrow(todo))
  tryCatch(cell_calls(todo$dataset[i], todo$sample[i]), error = function(e) {
    message("    [warn] ", todo$sample[i], ": ", conditionMessage(e)); NULL })
}), fill = TRUE)
if (!nrow(pc)) stop("[10] no per-cell data produced")
pc[todo, library_type := i.library_type, on = c("dataset", "sample")]

## -- Step 3. calibrate the call on the sorted libraries ----
message("[3] calibrating: FPR measured on CD34-sorted libraries, sensitivity on stroma-enriched")
grid <- CJ(min_fibro = 2:5, max_haem = 0:2)
cal <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  mf <- grid$min_fibro[i]; mh <- grid$max_haem[i]
  pc[, hit := fibro_msc >= mf & haematopoietic <= mh]
  s <- pc[, .(n_cells = .N, n_hit = sum(hit)), by = library_type]
  data.table(min_fibro = mf, max_haem = mh,
             fpr_sorted   = s[library_type == "sorted_CD34",     n_hit / n_cells],
             rate_enrich  = s[library_type == "stroma_enriched", n_hit / n_cells],
             rate_wholeMNC= s[library_type == "whole_MNC",       n_hit / n_cells],
             n_hit_sorted = s[library_type == "sorted_CD34", n_hit])
}))
cal[, fp_per_10k_sorted := 1e4 * fpr_sorted]
# The operating point: the LOOSEST rule (most sensitive) whose sorted-library false-positive rate
# stays under 1 in 10,000 cells. Stated as a rule before looking at which samples it favours.
cal[, acceptable := fpr_sorted < 1e-4]
setorder(cal, -acceptable, -rate_enrich)
fwrite_safe(cal, file.path(LCC_TAB_DIR, "10_stromal_calibration.csv"))
print(cal[, .(min_fibro, max_haem, fp_per_10k_sorted = round(fp_per_10k_sorted, 2),
              rate_enrich = round(rate_enrich, 4), rate_wholeMNC = round(rate_wholeMNC, 5), acceptable)])
op <- cal[acceptable == TRUE][1]
if (!nrow(op) || is.na(op$min_fibro)) { op <- cal[order(fpr_sorted)][1]
  message("    [!] no grid point reaches FPR < 1e-4; falling back to the strictest available") }
MIN_FIBRO <- op$min_fibro; MAX_HAEM <- op$max_haem
message(sprintf("    operating point: fibro_msc >= %d & haematopoietic <= %d  (sorted FPR %.2e)",
                MIN_FIBRO, MAX_HAEM, op$fpr_sorted))

## -- Step 4. apply, and subtype the calls ----
message("[4] calling de-novo stromal cells and subtyping them")
pc[, stromal_denovo := fibro_msc >= MIN_FIBRO & haematopoietic <= MAX_HAEM]
# Subtype order matters: adipocyte and osteolineage are tested BEFORE the MSC default, because both
# arise FROM the MSC and so still carry the fibro_msc genes that would otherwise capture them.
pc[stromal_denovo == TRUE, stromal_subtype := fifelse(adipocyte >= 2, "adipocyte",
                                              fifelse(osteolineage >= 2, "osteolineage",
                                               fifelse(chondrocyte >= 2, "chondrocyte",
                                                fifelse(endothelial >= 2 & fibro_msc < 4, "endothelial",
                                                 fifelse(pericyte >= 3 & fibro_msc < 4, "pericyte", "MSC_fibroblast")))))]
# adipogenic PRIMING inside the MSC pool -- the readout droplet data can actually support, as opposed
# to the mature-adipocyte count above which the platform destroys.
pc[, adipo_primed := stromal_subtype == "MSC_fibroblast" & adipo_prime >= 3]
smp <- pc[, .(n_cells = .N, n_stromal_denovo = sum(stromal_denovo),
              n_MSC_fibroblast = sum(stromal_subtype == "MSC_fibroblast", na.rm = TRUE),
              n_adipocyte      = sum(stromal_subtype == "adipocyte",      na.rm = TRUE),
              n_MSC_adipo_primed = sum(adipo_primed, na.rm = TRUE),
              n_endothelial    = sum(stromal_subtype == "endothelial",    na.rm = TRUE),
              n_pericyte       = sum(stromal_subtype == "pericyte",       na.rm = TRUE),
              n_osteolineage   = sum(stromal_subtype == "osteolineage",   na.rm = TRUE),
              n_chondrocyte    = sum(stromal_subtype == "chondrocyte",    na.rm = TRUE)),
          by = .(dataset, sample, library_type)]
smp[, frac_stromal := n_stromal_denovo / n_cells]
grp <- fread(file.path(LCC_TAB_DIR, "09_tp53_groups.csv"))[, .(dataset, sample, timepoint, uid_patient,
                                                               tp53_tier, tp53_group)]
smp <- merge(smp, grp, by = c("dataset", "sample"), all.x = TRUE)
setorder(smp, -n_MSC_fibroblast)
fwrite_safe(smp, file.path(LCC_TAB_DIR, "10_stromal_denovo_sample.csv"))

## -- Step 5. de-novo call vs the projection label ----
message("[5] cross-tabulating against the BoneMarrowMap projection")
xt <- rbindlist(lapply(seq_len(nrow(todo)), function(i) {
  f <- file.path(LCC_PERCELL_DIR, todo$dataset[i], paste0(todo$sample[i], "__lcc_percell.csv.gz"))
  if (!file.exists(f)) return(NULL)
  b <- fread(f, select = c("cell", "hierarchy_bin"), showProgress = FALSE)
  b[, `:=`(dataset = todo$dataset[i], sample = todo$sample[i])]
  b[]
}), fill = TRUE)
if (nrow(xt)) {
  pc[xt, hierarchy_bin := i.hierarchy_bin, on = c("dataset", "sample", "cell")]
  cross <- pc[!is.na(hierarchy_bin),
              .(n = .N), by = .(proj_stromal = hierarchy_bin == "Stromal", stromal_denovo)]
  agree <- pc[!is.na(hierarchy_bin), .(n_cells = .N,
                proj_stromal = sum(hierarchy_bin == "Stromal"),
                denovo_stromal = sum(stromal_denovo),
                both = sum(hierarchy_bin == "Stromal" & stromal_denovo),
                proj_only = sum(hierarchy_bin == "Stromal" & !stromal_denovo),
                denovo_only = sum(hierarchy_bin != "Stromal" & stromal_denovo)),
              by = .(dataset, sample, library_type)]
  setorder(agree, -denovo_only)
  fwrite_safe(agree, file.path(LCC_TAB_DIR, "10_stromal_vs_projection.csv"))
  cat("\n-- cell-level agreement, whole cohort --\n"); print(cross)
}

message("[6] summary")
cat("\n-- niche composition by library type: what the 'stromal' compartment is made of --\n")
comp <- smp[, lapply(.SD, sum), by = library_type,
            .SDcols = c("n_cells", "n_stromal_denovo", "n_MSC_fibroblast", "n_adipocyte",
                        "n_MSC_adipo_primed", "n_endothelial", "n_pericyte", "n_osteolineage",
                        "n_chondrocyte")]
comp[, pct_stromal := round(100 * n_stromal_denovo / n_cells, 4)]
print(comp)
fwrite_safe(comp, file.path(LCC_TAB_DIR, "10_niche_composition.csv"))
cat("\n-- samples with a usable de-novo MSC/fibroblast population (>= 30 cells) --\n")
print(smp[n_MSC_fibroblast >= 30, .(dataset, sample, timepoint, library_type, n_cells,
                                    n_MSC_fibroblast, n_endothelial, tp53_tier)])
cat("\n-- TP53-aberrant (tier A/B) samples: any de-novo stroma at all? --\n")
print(smp[tp53_tier %in% c("A_genotype", "B_allele_loh"),
          .(dataset, sample, timepoint, tp53_tier, n_cells, n_stromal_denovo, n_MSC_fibroblast)])
cat("\n-- projection FALSE POSITIVES: labelled Stromal but not called de novo --\n")
if (exists("agree")) print(agree[proj_only >= 20, .(dataset, sample, library_type, proj_stromal,
                                                    both, proj_only, denovo_only)][order(-proj_only)])
message("[done] 10_stromal_denovo")
