#!/usr/bin/env Rscript
# 16_dsp_handoff.R ----
# Build the two things the DSP group (Kao lab, "TP53 Kao discussion") can use directly:
#
#   (1) A SpatialDecon cell-profile matrix built from THIS project's AML marrow scRNA-seq compendium.
#       Their GeoMx AOIs are segmented into three masks -- CD45-bright / CD68+ / Leftover -- and
#       slide 6 of their deck states the purpose of H0 as separating "this compartment got bigger"
#       from "expression changed inside the compartment". Three area fractions cannot do that.
#       Deconvolution can, and it needs a profile matrix. The one shipped with SpatialDecon (safeTME)
#       is a normal-tissue immune reference: it has no AML blast column, no marrow stromal columns,
#       and no separation of macrophage from monocyte. All three gaps are exactly where their
#       questions live, so a bespoke AML-marrow reference is worth building.
#
#   (2) A small pre-specified gene-set file (GMT + a supporting table) derived from OUR results, on a
#       DIFFERENT platform and a DIFFERENT cohort, BEFORE seeing their per-gene data. Slide 19 of
#       their deck sets the right rule -- "不用看結果後再換 gene set", prefer external predefined
#       collections. An externally derived set is the legitimate way to satisfy that rule with a
#       bespoke set. And our own P4 is the evidence that the generic collections they picked are the
#       wrong instrument here: across 500 valid control sets the best of 14 MSigDB sets (3 TGF-beta,
#       3 ECM, EMT, hypoxia, 3 inflammatory) reached p<0.05 in 5.8% -- chance -- 8 of 14 never did,
#       and all three TGF-beta sets contain ZERO of our significant genes. A 200-gene rank score
#       cannot register a 5-gene shift.
#
# INPUT  : LCC_PERCELL_DIR/<ds>/<sample>__lcc_percell.csv.gz   cell -> bmm_broad / hierarchy_bin /
#                                                              malignant / myeloid_class
#          QC_RDS_DIR/<ds>/<sample>.rds                        raw counts
#          LCC_TAB_DIR/10_stromal_calibration.csv              the de-novo stromal operating point
#          LCC_TAB_DIR/11_gene_results.csv                     per-gene paired statistics, for (2)
# OUTPUT : LCC_DSP_DIR/16_profile_matrix_AMLmarrow.tsv         genes x cell types, CP10K units
#          LCC_DSP_DIR/16_profile_celltype_meta.csv            cells/samples/RNA content per column
#          LCC_DSP_DIR/16_profile_qc_markers.csv               marker specificity self-check
#          LCC_DSP_DIR/16_genesets_LCC.gmt                     (2), GMT for fry/roast/camera/GSEA
#          LCC_DSP_DIR/16_genesets_LCC_evidence.csv            every gene with the stat it came from
# Usage  : Rscript LCC_proj/scripts/16_dsp_handoff.R [--max_cells_per_type 3000] [--min_cells 25]
#                                                    [--min_samples 3]
#
# Deliberately does NOT load Seurat: the user library carries Seurat 5.5.1, which needs
# promises >= 1.5.0 while the env has 1.3.3, so loading it aborts. 03 and 10 hit the same wall and
# solved it the same way -- utils.R::get_counts() reaches the counts through SeuratObject alone.
suppressPackageStartupMessages({
  library(data.table); library(here); library(Matrix); library(optparse)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
source(here::here("LCC_proj", "scripts", "stromal_modules.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
set.seed(SEED)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--max_cells_per_type", type = "integer", default = 3000L,
              help = "cap cells per (sample, cell type) before averaging [3000]"),
  make_option("--min_cells", type = "integer", default = 25L,
              help = "a sample contributes to a cell type only with >= this many cells [25]"),
  make_option("--min_samples", type = "integer", default = 3L,
              help = "a cell type is kept only if >= this many samples contribute [3]"),
  make_option("--gene_presence", type = "double", default = 0.8,
              help = "keep genes present in the reference of >= this fraction of contributing samples [0.8]")
)))

LCC_DSP_DIR <- file.path(LCC_DIR, "results", "dsp_handoff")
dir.create(LCC_DSP_DIR, recursive = TRUE, showWarnings = FALSE)

## ============================================================================================
## PART 1 -- SpatialDecon profile matrix
## ============================================================================================

## -- Step 1. cell-type vocabulary ----
# Granularity is chosen to match the questions their deck actually asks, not to be maximal. Three
# AOI masks cannot support twenty columns, and every extra column costs collinearity in the
# deconvolution. Each column below earns its place:
#
#   LMPP_GMP      -- AZU1/PRTN3/ELANE/CTSG/MPO/DEFA1B are the azurophilic-granule programme of
#                    promyelocytes/GMP, and they are the top "higher in nonfibrotic" genes in ALL
#                    FOUR of their contrasts (H1, CD45, CD68, Leftover). A fibrosis effect would not
#                    be identical across three different segmentation masks; a case-composition
#                    difference would. This column tests that directly.
#   Macrophage    -- separated from Monocyte on purpose: their CD68+ mask is a macrophage mask, and
#                    C1QC/LGMN are the only lineage-appropriate genes to reach any of their top-10
#                    lists (CD68+ segment, FDR 2.8e-7 / 1.4e-7). Our C1QB is the single most robust
#                    gene in our whole panel (7/9 pairs, p=0.008, 90% of control sets).
#   Plasma        -- separated from B on purpose: IGHG1/IGKC are top "higher in fibrosis" genes in
#                    all four of their contrasts.
#   AML_blast     -- safeTME has no such column, and their tissue is majority blast.
#   MSC_fibroblast / Osteolineage / Adipocyte / Endothelial / Pericyte
#                 -- the compartment that makes collagen, and the one their Leftover mask dilutes.
#                    De-novo calls, not projection labels; see stromal_modules.R.
HAEM_RULES <- list(
  HSC_MPP       = list(col = "hierarchy_bin", vals = "HSC_MPP"),
  LMPP_GMP      = list(col = "hierarchy_bin", vals = "LMPP_GMP"),
  # NAMED Erythroid_MEP, NOT Erythroid, because that is what it contains. A first build called it
  # Erythroid and the QC caught the lie: HBB 0.41 CP10K in this column (lower than the ambient 1.05
  # in AML_blast), GYPA / ALAS2 / SLC4A1 flat zero in EVERY column. There is not one erythroblast in
  # the whole compendium -- droplet prep plus the QC gene-count filter removes them, so the bin holds
  # MEP and early progenitors only. Consequence for the DSP group, and it is a real one: marrow
  # trephine AOIs DO contain erythroid islands and this reference cannot represent them, so that
  # signal will be pushed onto whichever column is nearest. Flagged in the meta table and README.
  Erythroid_MEP = list(col = "hierarchy_bin", vals = "Erythroid"),
  Megakaryocyte = list(col = "hierarchy_bin", vals = "Megakaryocyte"),
  T_NK          = list(col = "hierarchy_bin", vals = "T_NK"),
  Monocyte      = list(col = "bmm_broad", vals = c("Monocyte", "Pro-Monocyte")),
  DC            = list(col = "bmm_broad", vals = c("cDC", "pDC")),
  B             = list(col = "bmm_broad", vals = c("B", "Pre-B", "Pro-B")),
  Plasma        = list(col = "bmm_broad", vals = "Plasma Cell")
)
# Chondrocyte is deliberately folded into MSC_fibroblast rather than shipped as its own column. The
# first build shipped it and the QC showed it is not chondrocyte: ACAN 0.6 and COL2A1 0.9 CP10K (the
# defining markers, essentially absent) against DCN 59.5, COL3A1 40.5, LUM 29.3 -- i.e. the column
# was simply the deepest-sequenced fibroblasts, median 18,011 UMI/cell vs 5,470 for MSC_fibroblast.
# That is the known failure mode of a detection-count subtype rule: more UMIs -> more modules fire ->
# the cell is pulled out of the MSC default into a rarer subtype. Those 555 cells are real stroma, so
# they are kept, just under the label the marker evidence supports.
# Adipo_MSC likewise: ADIPOQ reaches only 1.8 CP10K, so this is the adipogenic/adipo-CAR MSC that
# droplet data CAN see, not the mature adipocyte it destroys (see stromal_modules.R).
# label -> the 10_stromal_denovo subtype(s) it collects. A list, not a named vector, because
# MSC_fibroblast legitimately collects two subtypes and `x[[name]]` on a vector with duplicate names
# silently returns only the first.
STROMAL_TYPES <- list(MSC_fibroblast = c("MSC_fibroblast", "chondrocyte"),
                      Osteolineage   = "osteolineage",
                      Adipo_MSC      = "adipocyte",
                      Endothelial    = "endothelial",
                      Pericyte       = "pericyte")

op <- stromal_operating_point(file.path(LCC_TAB_DIR, "10_stromal_calibration.csv"))
message(sprintf("[0] stromal operating point read back from 10: fibro_msc >= %d & haematopoietic <= %d (sorted FPR %.2e)",
                op$min_fibro, op$max_haem, op$fpr_sorted))

man <- load_sample_meta(include_stroma_ref = TRUE)
man[, library_type := fifelse(dataset %in% SORTED_LIBRARY_DATASETS | grepl("_CD34$", sample), "sorted_CD34",
                       fifelse(grepl("(?i)niche", sample) | dataset %in% LCC_STROMA_REF_DATASETS,
                               "stroma_enriched", "whole_MNC"))]
man[, pc_file  := file.path(LCC_PERCELL_DIR, dataset, paste0(sample, "__lcc_percell.csv.gz"))]
man[, rds_file := file.path(QC_RDS_DIR, dataset, paste0(sample, ".rds"))]
todo <- man[file.exists(pc_file) & file.exists(rds_file)]
message("[1] ", nrow(todo), " samples have both a per-cell table and a counts object")
print(todo[, .N, by = library_type])

## -- Step 2. per-sample, per-cell-type mean profile ----
# Two normalisation decisions, both of which matter for a 13-study compendium:
#   (a) each cell is scaled to CP10K BEFORE averaging, so a deeply sequenced study does not dominate
#       a shallow one. The profile is therefore "expression per average cell", not "per read".
#   (b) the cross-sample aggregate (step 3) is a MEDIAN of per-sample means, not a pooled mean.
#       This is the same lesson the main analysis learned the hard way: a cell-weighted pool lets one
#       15,148-cell sample set the answer on its own, and it disagreed in sign with the sample-level
#       measure for 51 of 145 panel genes.
# Column sums are therefore equal across cell types, i.e. the matrix is in RNA-composition units.
# Cell-abundance units need the RNA-content scaler written out in step 4; see the README.
profile_one <- function(i) {
  ds <- todo$dataset[i]; sm <- todo$sample[i]
  pc <- fread(todo$pc_file[i], select = c("cell", "bmm_broad", "hierarchy_bin",
                                          "malignant", "myeloid_class"), showProgress = FALSE)
  cnt <- get_counts(readRDS(todo$rds_file[i]))
  pc <- pc[match(colnames(cnt), cell)]                 # align to the counts object, drop the rest
  keep <- !is.na(pc$cell)
  if (sum(keep) < opt$min_cells) return(NULL)
  cnt <- cnt[, keep, drop = FALSE]; pc <- pc[keep]

  # de-novo stromal call, identical rule to 10 (shared file, shared operating point)
  mods <- stromal_module_counts(cnt)
  is_str <- mods$fibro_msc >= op$min_fibro & mods$haematopoietic <= op$max_haem
  sub    <- rep(NA_character_, nrow(mods))
  if (any(is_str)) sub[is_str] <- stromal_subtype_of(mods[is_str])

  # label precedence: stroma first (it is the rarest and the most specifically defined), then
  # macrophage, then malignant blast, then the normal haematopoietic bins. A cell that passes the
  # stromal co-expression rule is never also counted as a blast.
  lab <- rep(NA_character_, nrow(pc))
  for (nm in names(STROMAL_TYPES)) lab[is.na(lab) & !is.na(sub) & sub %in% STROMAL_TYPES[[nm]]] <- nm
  lab[is.na(lab) & pc$myeloid_class == "macrophage_like"] <- "Macrophage"
  lab[is.na(lab) & pc$malignant %in% c(1L, TRUE, "1", "TRUE")] <- "AML_blast"
  for (nm in names(HAEM_RULES)) {
    r <- HAEM_RULES[[nm]]
    lab[is.na(lab) & pc[[r$col]] %in% r$vals] <- nm
  }

  cs <- Matrix::colSums(cnt)
  ok <- !is.na(lab) & cs > 0
  if (!any(ok)) return(NULL)

  # map this sample's gene symbols onto the shared vocabulary, extending it with anything new. The
  # 13 studies were quantified against different references, so the symbol sets genuinely differ --
  # which is exactly why gene presence has to be tracked rather than assumed.
  gn  <- rownames(cnt)
  hit <- match(gn, GVOCAB)
  if (anyNA(hit)) {
    new <- unique(gn[is.na(hit)])
    GVOCAB <<- c(GVOCAB, new)
    hit <- match(gn, GVOCAB)
  }
  gidx <- hit

  out <- list()
  for (nm in unique(lab[ok])) {
    idx <- which(ok & lab == nm)
    if (length(idx) < opt$min_cells) next
    if (length(idx) > opt$max_cells_per_type) idx <- sample(idx, opt$max_cells_per_type)
    sub_cnt <- cnt[, idx, drop = FALSE]
    # CP10K per cell, then mean across cells: scale columns by 1e4/total, then row means
    cp <- sub_cnt %*% Matrix::Diagonal(x = 1e4 / cs[idx])
    v  <- as.numeric(Matrix::rowMeans(cp))
    # MEMORY: keep only the nonzero entries. A first version stored the full dense gene x celltype
    # profile with `gene` as a character column and was OOM-killed at sample 180 of 220 -- roughly
    # 55M rows carrying four character vectors. Genes are carried as an integer index into a shared
    # vocabulary and zeros are dropped; the zeros are restored exactly at aggregation time (step 3),
    # so the median is unchanged, not approximated.
    nz <- which(v > 0)
    if (!length(nz)) next
    out[[nm]] <- list(g = gidx[nz], v = v[nz], ct = nm,
                      n_cells = length(idx), median_umi = as.numeric(median(cs[idx])))
  }
  if (!length(out)) return(NULL)
  list(vals = rbindlist(lapply(out, function(o)
         data.table(g = o$g, v = o$v, ct = o$ct))),
       ref_genes = gidx,
       cts = rbindlist(lapply(out, function(o)
         data.table(celltype = o$ct, dataset = ds, sample = sm,
                    n_cells = o$n_cells, median_umi = o$median_umi))))
}

message("[2] building per-sample profiles")
GVOCAB <- character(0)                     # shared gene vocabulary; profile_one() indexes into it
# gene presence per celltype: how many of the samples that profiled celltype k carried gene g at all.
# This is what makes the zero-restoration in step 3 exact -- a gene missing from a sample's reference
# is "not measured there", which is not the same as "measured as zero".
PRESENT <- list()
GPRESENT <- integer(0)                     # over ALL samples, irrespective of celltype
grow <- function(x) if (length(x) < length(GVOCAB)) c(x, integer(length(GVOCAB) - length(x))) else x
res <- vector("list", nrow(todo)); cts <- vector("list", nrow(todo))
for (i in seq_len(nrow(todo))) {
  if (i %% 10 == 0) message("    ", i, "/", nrow(todo))
  r <- tryCatch(profile_one(i), error = function(e) {
    message("    [warn] ", todo$sample[i], ": ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  res[[i]] <- r$vals; cts[[i]] <- r$cts
  GPRESENT <- grow(GPRESENT)
  GPRESENT[r$ref_genes] <- GPRESENT[r$ref_genes] + 1L
  for (k in unique(r$cts$celltype)) {
    PRESENT[[k]] <- grow(if (is.null(PRESENT[[k]])) integer(0) else PRESENT[[k]])
    PRESENT[[k]][r$ref_genes] <- PRESENT[[k]][r$ref_genes] + 1L
  }
}
prof <- rbindlist(res); sample_n <- unique(rbindlist(cts)); rm(res, cts); invisible(gc())
if (!nrow(prof)) stop("[16] no profiles produced")

## -- Step 3. aggregate across samples ----
message("[3] aggregating across samples (median of per-sample means)")
ct_keep <- sample_n[, .(n_samples = .N), by = celltype][n_samples >= opt$min_samples, celltype]
dropped <- setdiff(unique(sample_n$celltype), ct_keep)
if (length(dropped))
  message("    dropped for < ", opt$min_samples, " contributing samples: ", paste(dropped, collapse = ", "))
prof <- prof[ct %in% ct_keep]

# A gene absent from a sample's reference is "not measured there", not zero. Restrict the universe to
# genes carried by most contributing samples rather than silently imputing zeros.
n_tot <- nrow(unique(sample_n[, .(dataset, sample)]))
GPRESENT <- grow(GPRESENT)
genes_keep_idx <- which(GPRESENT / n_tot >= opt$gene_presence)
message(sprintf("    gene universe: %d of %d symbols present in >= %.0f%% of %d samples",
                length(genes_keep_idx), length(GVOCAB), 100 * opt$gene_presence, n_tot))
prof <- prof[g %in% genes_keep_idx]

# Median of per-sample means, with the dropped zeros restored. All stored values are > 0, so after
# sorting within (gene, celltype) the padded vector is [pad zeros, sorted values] and the median can
# be read off by position -- no need to materialise the zeros.
pad <- rbindlist(lapply(ct_keep, function(k)
  data.table(g = genes_keep_idx, ct = k, n_meas = grow(PRESENT[[k]])[genes_keep_idx])))
prof <- merge(prof, pad, by = c("g", "ct"), all.x = TRUE)
setorder(prof, g, ct, v)
agg <- prof[, {
  k <- .N; p <- max(n_meas[1] - k, 0L); n <- k + p
  at <- function(ii) if (ii <= p) 0 else v[ii - p]
  .(value = if (n %% 2L == 1L) at((n + 1L) / 2L) else (at(n / 2L) + at(n / 2L + 1L)) / 2)
}, by = .(g, ct)]
agg[, `:=`(gene = GVOCAB[g], celltype = ct)]
X <- dcast(agg, gene ~ celltype, value.var = "value", fill = 0)
setcolorder(X, c("gene", sort(setdiff(names(X), "gene"))))
# rescale each column to sum 1e4 so that columns are strictly comparable after the median step
for (cc in setdiff(names(X), "gene")) {
  s <- sum(X[[cc]], na.rm = TRUE)
  if (s > 0) set(X, j = cc, value = X[[cc]] * 1e4 / s)
}
fwrite(X, file.path(LCC_DSP_DIR, "16_profile_matrix_AMLmarrow.tsv"), sep = "\t")
message("    wrote profile matrix: ", nrow(X), " genes x ", ncol(X) - 1L, " cell types")

## -- Step 4. per-column metadata, including the RNA-content scaler ----
meta <- sample_n[celltype %in% ct_keep,
                 .(n_samples = .N, n_datasets = uniqueN(dataset), n_cells_total = sum(n_cells),
                   median_umi_per_cell = median(median_umi)), by = celltype]
# SpatialDecon returns abundance in the units of the profile matrix. With equal column sums that is
# RNA abundance. Multiplying column k by rna_content_rel converts to approximate CELL abundance --
# which matters here because an erythroid cell and a megakaryocyte differ several-fold in RNA content
# and their AOIs contain both.
meta[, rna_content_rel := median_umi_per_cell / median(median_umi_per_cell)]

# Per-column honesty flags. A reference matrix that ships without these invites the user to trust
# every column equally, and these columns are NOT equally trustworthy.
meta[, flags := ""]
add_flag <- function(cond, tag) meta[cond, flags := paste0(flags, ifelse(nzchar(flags), ";", ""), tag)]
# marrow stroma exists only in physically enriched libraries, so several stromal columns necessarily
# come from one study; their profile is that study's, and cannot be checked against another.
add_flag(meta$n_datasets == 1L, "single_dataset")
add_flag(meta$n_cells_total < 1000L, "few_cells")
add_flag(meta$celltype == "Erythroid_MEP", "no_erythroblasts_see_README")
add_flag(meta$celltype == "Adipo_MSC", "adipogenic_MSC_not_mature_adipocyte")
setorder(meta, -n_cells_total)
fwrite(meta, file.path(LCC_DSP_DIR, "16_profile_celltype_meta.csv"))
print(meta)

## -- Step 5. QC: does the reference put the canonical markers where they belong? ----
# A profile matrix that fails this is not usable, and the failure would be invisible downstream.
# These are also precisely the genes that drive their DE lists, so this table doubles as the
# evidence for the composition argument.
QC <- list(
  granule_GMP  = c("AZU1", "PRTN3", "ELANE", "CTSG", "MPO", "DEFA1B", "PRG2"),
  # The SECOND maturity gap, same root cause as the erythroid one below and equally worth recording:
  # droplet prep plus the QC gene-count filter loses mature granulocytes (very low RNA content, high
  # ambient). Promyelocyte/GMP genes survive -- MPO 13.1 and ELANE 4.4 CP10K in LMPP_GMP -- but
  # everything from myelocyte onward does not: DEFA1B/CEACAM8 exactly 0, DEFA3/DEFA4/LTF/LCN2/CAMP/
  # FCGR3B/CXCR2 all under 0.2 CP10K and topping essentially random columns, i.e. noise.
  # This is not academic for the DSP group: DEFA1B is a top-4 "higher in nonfibrotic" gene in two of
  # their four contrasts, and this reference cannot attribute it to any cell type. It CAN attribute
  # AZU1/PRTN3/ELANE/CTSG/MPO, which is most of that signal.
  granulocyte_mature = c("DEFA3", "DEFA4", "LTF", "LCN2", "CAMP", "CEACAM8", "FCGR3B", "CXCR2"),
  macrophage   = c("C1QA", "C1QB", "C1QC", "LGMN", "MRC1", "LYVE1", "VSIG4", "CD163"),
  monocyte     = c("FCN1", "VCAN", "S100A12", "S100A8"),
  plasma_Ig    = c("IGHG1", "IGKC", "MZB1", "JCHAIN"),
  # ALAS2/GYPA/SLC4A1 are in this list precisely so that the shipped QC file RECORDS the gap rather
  # than hiding it: they are flat zero in every column, which is the evidence that the compendium
  # contains no erythroblasts. See the Erythroid_MEP comment above.
  erythroid    = c("HBB", "HBA1", "HBA2", "AHSP", "ALAS2", "GYPA", "SLC4A1"),
  megakaryo    = c("PF4", "PPBP", "GP9", "ITGA2B", "VWF"),
  T_NK         = c("CD3D", "CD3E", "IL7R", "GNLY", "NKG7"),
  HSPC         = c("CD34", "AVP", "CRHBP", "SPINK2"),
  fibro_MSC    = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "CXCL12"),
  osteo        = c("BGLAP", "IBSP", "SP7", "RUNX2"),
  endothelial  = c("PECAM1", "CDH5", "VWF", "EMCN"),
  ribosome     = c("RPL13", "RPS23", "RPL30", "RPS11")     # the "should be everywhere" control
)
Xm <- as.matrix(X[, -1]); rownames(Xm) <- X$gene
qc <- rbindlist(lapply(names(QC), function(k) {
  g <- intersect(QC[[k]], rownames(Xm))
  if (!length(g)) return(NULL)
  rbindlist(lapply(g, function(gg) {
    v <- Xm[gg, ]
    data.table(marker_set = k, gene = gg,
               top_celltype = names(v)[which.max(v)],
               top_value = round(max(v), 2),
               share_of_total = round(max(v) / sum(v), 3))
  }))
}))
fwrite(qc, file.path(LCC_DSP_DIR, "16_profile_qc_markers.csv"))
cat("\n-- marker specificity check (each gene's highest-expressing column) --\n")
print(qc, nrows = 100)

## -- Step 5b. collinearity, because it decides how many columns are actually usable ----
# Deconvolution cannot separate columns it cannot tell apart. Their AOIs come from three masks, so
# the practical question is not "is the profile correct" but "are these columns distinguishable".
# Shipping this table means they can see for themselves which estimates to read jointly.
Lg <- log1p(Xm)
cr <- cor(Lg)
ij <- which(upper.tri(cr), arr.ind = TRUE)
coll <- data.table(celltype_a = colnames(cr)[ij[, 1]], celltype_b = colnames(cr)[ij[, 2]],
                   r_log_profile = round(cr[ij], 3))
setorder(coll, -r_log_profile)
fwrite(coll, file.path(LCC_DSP_DIR, "16_profile_collinearity.csv"))
cat("\n-- most collinear column pairs (read these jointly, not as independent estimates) --\n")
print(head(coll, 12))

## ============================================================================================
## PART 2 -- pre-specified gene sets
## ============================================================================================
message("\n[6] writing the pre-specified gene sets")
# Curated from OUR paired sample-level results (11_gene_results.csv, design A, 9 matched pairs),
# on scRNA-seq, in a different cohort, before seeing their per-gene tables. Each set is small and
# DIRECTIONAL by construction: the hypothesis in every case is "higher in the fibrotic / TP53-mutant
# arm". A 7-gene directional set is far better powered at n=3 cases than a 1000-gene matrisome.
#
# Membership is fixed here rather than thresholded from the data on purpose -- a set whose membership
# moves when the p-value cutoff moves is not pre-specified.
SETS <- list(
  LCC_TP53_MACROPHAGE_C1Q = c("C1QA", "C1QB", "C1QC", "MAF", "LYVE1", "MRC1", "VSIG4", "CD163"),
  LCC_COLLAGEN_CROSSLINK  = c("PLOD2", "LOX", "LOXL1", "LOXL2", "MMP14", "TIMP1", "SERPINH1"),
  LCC_PROFIBROTIC_EFFECTOR= c("SPP1", "OSM", "PDGFB", "PDGFRB", "IL11", "TGFB3"),
  LCC_HYPOXIA_GLYCOLYSIS  = c("SLC2A1", "PGK1", "LDHA", "NDRG1", "ADM", "BNIP3"),
  LCC_MEGAKARYOCYTE_AXIS  = c("GP9", "ITGA2B", "GATA1", "PF4", "PPBP", "VWF", "MPL", "NFE2"),
  # the negative control, and it is the point of including it: if their pipeline reports this set as
  # moving while the four above do not, the result is a normalisation artefact, not fibrosis. We find
  # these flat or LOWER in the mutant arm (COL1A1 4/9, COL1A2 2/9, COL3A1 2/9, DCN 2/9, POSTN 1/9)
  # and their own four top-20 lists contain no collagen gene at all.
  LCC_COLLAGEN_STRUCTURAL_NEGCTRL = c("COL1A1", "COL1A2", "COL3A1", "COL6A3", "DCN", "LUM", "POSTN")
)
gmt <- vapply(names(SETS), function(k)
  paste(c(k, "LCC_proj scRNA TP53-mut vs WT, 9 matched pairs", SETS[[k]]), collapse = "\t"), "")
writeLines(gmt, file.path(LCC_DSP_DIR, "16_genesets_LCC.gmt"))

# Every gene travels with the statistic it came from, so nothing has to be taken on trust.
gr <- fread(file.path(LCC_TAB_DIR, "11_gene_results.csv"))[stratum == "all"]
rb <- fread(file.path(LCC_TAB_DIR, "12_robustness_genes.csv"))[design == "A_sample_strict" & kind == "gene"]
ev <- rbindlist(lapply(names(SETS), function(k) data.table(gene_set = k, gene = SETS[[k]])))
ev[gr, `:=`(n_pairs = i.n_pairs, frac_pairs_mut_higher = i.frac_pairs_mut_higher,
            p_sample_higher = i.p_sample_higher, median_delta_log2cpm = i.median_delta_log2cpm,
            panel_category = i.category), on = "gene"]
ev[rb, robustness_share := i.frac_p05, on = c(gene = "feature")]
ev[, pairs_higher := ifelse(is.na(frac_pairs_mut_higher), NA_character_,
                            paste0(round(frac_pairs_mut_higher * n_pairs), "/", n_pairs))]
ev[, in_our_panel := !is.na(n_pairs)]
setcolorder(ev, c("gene_set", "gene", "in_our_panel", "pairs_higher", "p_sample_higher",
                  "robustness_share", "median_delta_log2cpm", "panel_category"))
fwrite(ev[, -c("frac_pairs_mut_higher")], file.path(LCC_DSP_DIR, "16_genesets_LCC_evidence.csv"))
cat("\n-- gene sets with the evidence behind each member --\n")
print(ev[, .(gene_set, gene, pairs_higher, p = round(p_sample_higher, 3),
             robust = robustness_share)], nrows = 100)

message("\n[done] outputs in ", LCC_DSP_DIR)
