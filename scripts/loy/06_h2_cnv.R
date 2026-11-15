#!/usr/bin/env Rscript
# 06_h2_cnv.R ----
# input : results/tables/loy/04_h2/stage0_cells_<case>.csv.gz (04: leukaemic cells with their -Y call and fine label), the library's
#         per-cell chrY table (LOY_PERCELL_DIR; its normal T / NK / B cells become the autologous reference), raw counts (Zeng pack /
#         GSE269669 objects), INFERCNV_GENE_ORDER, LOY_H2_CYTOBAND (UCSC hg38), the BoneMarrowMap external reference (cache)
# output: results/tables/loy/04_h2/
#           cnv_regions_<case>.csv    per region: mean score in -Y / Y+ / reference cells, carrier frequencies, ratio, permutation p, verdict
#           cnv_cells_<case>.csv.gz   per cell: group, fine label, one column per region score
#           cnv_gate_<case>.csv       GATE_H2_CNV (PT43A decides on it; reported for every case) and the chrY machinery check
#         inferCNV run directory: INFERCNV_ROOT/loy_h2/<case>/ (run.final.infercnv_obj kept; re-used unless --force)
# what it does: PREREGISTRATION_H2.md section 6 / amendment 7. One inferCNV run per case: observation = every leukaemic cell with a
#   suspected_LOY or Y_retained call (ONE group, so the clustering never sees the label); reference = the library's own T / NK / B cells
#   (>= 20 per group) plus the project's BoneMarrowMap lineage groups for the lineages those cannot cover (the 02_malignancy
#   autologous route, via .assemble_foreign_ref); project cutoff and denoising; HMM off, analysis_mode "samples", no plots
#   (region means need the residual matrix only). Region score per cell = mean residual over the region's genes; regions = every
#   chromosome arm plus the case's pre-specified cytoband ranges. A cell carries a region's CNV when its score lies beyond
#   LOY_H2_CNV_CARRIER_SD reference SDs in the direction of the Y+ leukaemic mean; a leukaemic CNV has carrier frequency >=
#   LOY_H2_CNV_MIN_FREQ among Y+ cells; GATE_H2_CNV passes when one leukaemic CNV has a -Y / Y+ carrier ratio inside LOY_H2_CNV_RATIO.
# Usage : sbatch scripts/loy/submit_loy_h2.sbatch 06 --case=AML15 [--force]
# WHY: the PT35B pattern (a -Y clone of normal haematopoiesis next to a Y+ leukaemia) can only be told from a -Y leukaemic branch by
#   the lesions the two halves share; for AML15 the diagnostic karyotype predicts one region that separates and one that does not.
suppressPackageStartupMessages({ library(data.table); library(here); library(Matrix) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "loy", "config_loy.R"))
source(here::here("scripts", "loy", "00_loy_common.R"))
source(here::here("scripts", "02_malignancy", "00_infercnv_common.R"))     # sources config_malignancy / hierarchy / qc / utils, loads infercnv
args  <- commandArgs(trailingOnly = TRUE)
argv  <- function(k, default) { v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE)); if (length(v)) v[1] else default }
case  <- argv("case", NA_character_); force <- any(args == "--force")
stopifnot(case %in% LOY_H2_CASES$case)
lib   <- LOY_H2_CASES[case == get("case", envir = globalenv()) & role == "contrast"]; stopifnot(nrow(lib) == 1)
out_dir <- file.path(INFERCNV_ROOT, "loy_h2", case); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(LOY_H2_SEED); t0 <- Sys.time()
message(sprintf("[0] case %s | library %s | %s", case, lib$Sample, format(Sys.time())))

## [1] observation cells (04 handoff) and the autologous reference cells (per-cell table) --------------------
cl  <- fread(file.path(LOY_H2_DIR, sprintf("stage0_cells_%s.csv.gz", case)))
obs <- cl[arm == "leukaemic" & fine == "ALL", .(cell, group, fine = fine_depositor, nCount_RNA)]      # every leukaemic cell with a call
obs[cl[arm == "leukaemic" & fine != "ALL"], on = "cell", fine := i.fine]                             # the Stage-0 fine label (BMM for PT43)
pc  <- fread(file.path(LOY_PERCELL_DIR, paste0(lib$dataset, "__y_percell.csv.gz")))
pc[, Sample := fifelse(is.na(Sample) | Sample == "", as.character(object), as.character(Sample))]
pc  <- pc[Sample == lib$Sample]
pc[, comp := loy_comp_of(hierarchy_bin, bmm_fine)]
pc[, leuk := loy_h2_is_leuk(dataset, object, anno_malignant)]
ref_auto <- pc[!leuk & comp %in% c("T", "NK", "B_Plasma") & !cell %in% obs$cell, .(cell, group = paste0("reference_", sub("_Plasma", "", comp)))]
ref_auto <- ref_auto[, if (.N >= INFERCNV_REF_MIN_PER_GROUP) .SD, by = group]
auto_ok  <- nrow(ref_auto) >= LOY_H2_CNV_MIN_REF
message(sprintf("[1] observation %d cells (A %d, B %d) | autologous normal cells %d in %d groups -> %s",
                nrow(obs), sum(obs$group == "A"), sum(obs$group == "B"), nrow(ref_auto), uniqueN(ref_auto$group),
                if (auto_ok) "autologous lymphoid + external lineage groups" else "external lineage groups only (too few autologous normals)"))
if (!auto_ok) ref_auto <- ref_auto[0]

## [2] counts ------------------------------------------------------------------------------------------------
need <- c(obs$cell, ref_auto$cell)
if (lib$dataset == "GSE269669") {
  mats <- lapply(c("avm_aml", "avm_tme"), function(o) { s <- read_gse269669(o); s <- subset(s, cells = intersect(colnames(s), need)); m <- counts_of(s); rm(s); gc(); m })
  g <- Reduce(intersect, lapply(mats, rownames)); cnt <- do.call(cbind, lapply(mats, function(m) m[g, , drop = FALSE])); rm(mats)
} else {
  f <- file.path(LOY_ZENG_RDS_DIR, paste0(LOY_EXT_STUDIES[dataset == lib$dataset, file_stem], "_projected_cleaned.rds")); message("    reading ", basename(f))
  s <- readRDS(f); s <- subset(s, cells = intersect(colnames(s), need)); cnt <- counts_of(s); rm(s); invisible(gc())
}
stopifnot(all(need %in% colnames(cnt))); cnt <- cnt[, need, drop = FALSE]
message(sprintf("[2] counts %d genes x %d cells (%.1f min)", nrow(cnt), ncol(cnt), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## [3] reference assembly (02_malignancy autologous route) and the inferCNV run ---------------------------------
final_rds <- file.path(out_dir, "run.final.infercnv_obj")
if (file.exists(final_rds) && !force) {
  message("[3] re-using ", final_rds); obj <- readRDS(final_rds)
} else {
  drop <- if (auto_ok) c("T_NK", "B_Plasma") else character(0)
  fb <- .assemble_foreign_ref(cnt, lib$dataset, lib$Sample, drop_bins = drop, ext_ref_fn = get_external_ref, matched_ref_fn = get_matched_ref)
  stopifnot(!is.null(fb))
  counts <- cbind(cnt[fb$genes, , drop = FALSE], fb$counts)
  anno   <- data.table(cell = colnames(counts), group = c(rep("observation", nrow(obs)), ref_auto$group, fb$groups))
  rg     <- setdiff(unique(anno$group), "observation")
  message(sprintf("    reference groups: %s", paste(anno[group != "observation", .N, by = group][, paste0(group, "=", N)], collapse = " ")))
  anno_file <- file.path(out_dir, "annotations.tsv"); fwrite(anno, anno_file, sep = "\t", col.names = FALSE)
  # chr_exclude keeps chrX and chrY (package default drops them): chrY is this arm's machinery check. Only chrM is dropped.
  obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = counts, annotations_file = anno_file, gene_order_file = INFERCNV_GENE_ORDER,
                                        ref_group_names = rg, chr_exclude = "chrM")
  obj <- infercnv::run(obj, cutoff = INFERCNV_CUTOFF, out_dir = out_dir, cluster_by_groups = TRUE, analysis_mode = "samples",
                       denoise = INFERCNV_DENOISE, HMM = FALSE, num_threads = INFERCNV_THREADS, no_plot = TRUE, save_rds = TRUE)
  rm(counts); invisible(gc())
}
E  <- obj@expr.data - 1                                   # genes x cells, residual relative to the reference interval (centred at 1)
go <- as.data.table(obj@gene_order, keep.rownames = "gene"); setnames(go, tolower(names(go)))
if (!"stop" %in% names(go) && "end" %in% names(go)) setnames(go, "end", "stop")
go[, chr := as.character(chr)]
message(sprintf("[3] residual matrix %d genes x %d cells (%.1f min)", nrow(E), ncol(E), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## [4] regions: every arm from the cytoband table + the case's pre-specified bands --------------------------------
cb <- fread(LOY_H2_CYTOBAND, header = FALSE, col.names = c("chr", "start", "end", "band", "stain"))
cb <- cb[chr %in% paste0("chr", c(1:22, "X", "Y"))]
arms <- cb[, .(p_start = 0, p_end = min(start[stain == "acen"]), q_start = max(end[stain == "acen"]), q_end = max(end)), by = chr]
regions <- rbind(arms[, .(region = paste0(chr, "p"), chr, start = p_start, end = p_end)], arms[, .(region = paste0(chr, "q"), chr, start = q_start, end = q_end)])
pre <- LOY_H2_CNV_REGIONS[case == get("case", envir = globalenv())]
if (nrow(pre)) {
  pre_rng <- pre[, {
    b1 <- cb[chr == .BY$chr & band == band_from]; b2 <- cb[chr == .BY$chr & band == band_to]
    stopifnot(nrow(b1) == 1, nrow(b2) == 1)
    .(region = region, start = min(b1$start, b2$start), end = max(b1$end, b2$end), expectation = expectation) }, by = chr]
  regions <- rbind(regions, pre_rng[, .(region, chr, start, end)])
}
reg_genes <- lapply(seq_len(nrow(regions)), function(i) go[chr == regions$chr[i] & start >= regions$start[i] & stop <= regions$end[i], gene])
regions[, n_genes := lengths(reg_genes)]
keep_r <- regions$n_genes >= LOY_H2_CNV_MIN_GENES
message(sprintf("[4] %d regions scored (%d dropped with < %d genes: %s)", sum(keep_r), sum(!keep_r), LOY_H2_CNV_MIN_GENES, paste(regions$region[!keep_r], collapse = ",")))
regions <- regions[keep_r]; reg_genes <- reg_genes[keep_r]
# chrY as a whole, with a 5-gene floor: only a handful of chrY genes pass inferCNV's expression cutoff, yet this row is the arm's
# machinery check (prereg amendment 7), so it is scored separately from the 20-gene rule that governs every other region
y_genes <- go[chr == "chrY", gene]
if (length(y_genes) >= 5) { regions <- rbind(regions, data.table(region = "chrY_all", chr = "chrY", start = 0, end = max(go[chr == "chrY", stop]), n_genes = length(y_genes))); reg_genes <- c(reg_genes, list(y_genes)); message(sprintf("    chrY_all added with %d genes", length(y_genes))) } else message(sprintf("    chrY_all NOT scored: only %d chrY genes in the residual matrix", length(y_genes)))
S <- sapply(reg_genes, function(g) colMeans(E[g, , drop = FALSE])); colnames(S) <- regions$region      # cells x regions
cell_tab <- data.table(cell = rownames(S))
cell_tab[obs, on = "cell", `:=`(group = i.group, fine = i.fine)]
cell_tab[is.na(group), group := "reference"]
fwrite(cbind(cell_tab, as.data.table(S)), file.path(LOY_H2_DIR, sprintf("cnv_cells_%s.csv.gz", case)))

## [5] per-region statistics ---------------------------------------------------------------------------------------
isA <- cell_tab$group == "A"; isB <- cell_tab$group == "B"; isR <- cell_tab$group == "reference"
mean_ref <- colMeans(S[isR, , drop = FALSE]); sd_ref <- apply(S[isR, , drop = FALSE], 2, sd)
mean_A <- colMeans(S[isA, , drop = FALSE]); mean_B <- colMeans(S[isB, , drop = FALSE])
dir_leuk <- sign(mean_B - mean_ref)                                                             # direction of the Y+ leukaemic deviation
thr <- mean_ref + dir_leuk * LOY_H2_CNV_CARRIER_SD * sd_ref
# carrier: score beyond the threshold in the leukaemic direction (multiplying by the direction turns both cases into ">")
carrier <- sweep(S, 2, dir_leuk, "*") > matrix(dir_leuk * thr, nrow(S), ncol(S), byrow = TRUE)
freq_A <- colMeans(carrier[isA, , drop = FALSE]); freq_B <- colMeans(carrier[isB, , drop = FALSE]); freq_R <- colMeans(carrier[isR, , drop = FALSE])
# permutation of the -Y label among leukaemic cells, within fine label, for the difference of region means
ob <- which(isA | isB); grp <- isA[ob]; strata <- split(seq_along(ob), cell_tab$fine[ob])
delta <- mean_A - mean_B
perm  <- matrix(NA_real_, LOY_H2_CNV_N_PERM, ncol(S))
Sob   <- S[ob, , drop = FALSE]
for (p in seq_len(LOY_H2_CNV_N_PERM)) {
  gp <- grp; for (ix in strata) if (length(ix) > 1) gp[ix] <- grp[ix][sample.int(length(ix))]
  perm[p, ] <- colMeans(Sob[gp, , drop = FALSE]) - colMeans(Sob[!gp, , drop = FALSE])
}
p_perm <- (1 + colSums(abs(perm) >= matrix(abs(delta), nrow(perm), ncol(perm), byrow = TRUE))) / (LOY_H2_CNV_N_PERM + 1)
res <- data.table(case = case, region = regions$region, chr = regions$chr, n_genes = regions$n_genes,
                  mean_ref = mean_ref, sd_ref = sd_ref, mean_A_minusY = mean_A, mean_B_Yplus = mean_B, delta_A_minus_B = delta,
                  p_perm = p_perm, p_holm = p.adjust(p_perm, "holm"),
                  freq_ref = freq_R, freq_A_minusY = freq_A, freq_B_Yplus = freq_B, ratio_A_over_B = freq_A / freq_B)
res[, leukaemic_cnv := freq_B_Yplus >= LOY_H2_CNV_MIN_FREQ & !chr %in% c("chrX", "chrY")]   # set after construction: a data.table() argument cannot see its siblings (job 3689465)
res[, cosegregation := fifelse(!leukaemic_cnv, "not_a_leukaemic_cnv",
                       fifelse(ratio_A_over_B >= LOY_H2_CNV_RATIO[1] & ratio_A_over_B <= LOY_H2_CNV_RATIO[2], "shared_by_both_halves",
                       fifelse(ratio_A_over_B < LOY_H2_CNV_RATIO[1], "Yplus_only_or_enriched", "minusY_enriched")))]
if (nrow(pre)) res <- merge(res, pre[, .(region, expectation)], by = "region", all.x = TRUE)
fwrite(res[order(chr, region)], file.path(LOY_H2_DIR, sprintf("cnv_regions_%s.csv", case)))
message("[5] regions with a leukaemic CNV (freq in Y+ cells >= 30 %) or a significant -Y vs Y+ difference:")
print(res[leukaemic_cnv | p_holm < 0.05, .(region, n_genes, mean_ref = round(mean_ref, 3), mean_A = round(mean_A_minusY, 3), mean_B = round(mean_B_Yplus, 3),
                                           p_holm = signif(p_holm, 3), freq_A = round(freq_A_minusY, 2), freq_B = round(freq_B_Yplus, 2), ratio = round(ratio_A_over_B, 2), cosegregation)])
if (nrow(pre)) { message("    pre-specified regions:"); print(res[region %in% pre$region, .(region, mean_A = round(mean_A_minusY, 3), mean_B = round(mean_B_Yplus, 3), p_perm, freq_A = round(freq_A_minusY, 2), freq_B = round(freq_B_Yplus, 2), ratio = round(ratio_A_over_B, 2), cosegregation, expectation)]) }

## [6] gates -------------------------------------------------------------------------------------------------------
yrow <- res[chr == "chrY"]
gate <- data.table(case = case,
  gate   = c("CNV_machinery_chrY", "GATE_H2_CNV"),
  value  = c(if (nrow(yrow)) min(yrow$delta_A_minus_B) else NA_real_, sum(res$cosegregation == "shared_by_both_halves")),
  detail = c(if (nrow(yrow)) paste(yrow[, sprintf("%s delta %.3f p %.4f", region, delta_A_minus_B, p_perm)], collapse = " | ") else "no chrY region with >= 20 genes",
             sprintf("leukaemic CNVs: %d; shared: %s; Y+-only: %s; -Y-enriched: %s", sum(res$leukaemic_cnv),
                     paste(res[cosegregation == "shared_by_both_halves", region], collapse = ","), paste(res[cosegregation == "Yplus_only_or_enriched", region], collapse = ","),
                     paste(res[cosegregation == "minusY_enriched", region], collapse = ","))),
  pass   = c(nrow(yrow) > 0 && all(yrow$delta_A_minus_B < 0 & yrow$p_perm < 0.001), any(res$cosegregation == "shared_by_both_halves")),
  decides = c(FALSE, case == "PT43"))
fwrite(gate, file.path(LOY_H2_DIR, sprintf("cnv_gate_%s.csv", case))); print(gate[, .(gate, value = round(value, 3), pass, decides, detail = substr(detail, 1, 100))])
message(sprintf("[done] %s in %.1f min", case, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
