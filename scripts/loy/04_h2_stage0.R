#!/usr/bin/env Rscript
# 04_h2_stage0.R ----
# input : per-cell chrY tables of the H2 datasets (LOY_PERCELL_DIR/<ds>__y_percell.csv.gz, written by 01/03), the frozen
#         call tables of 02 (loy_reference_rates.csv, loy_external_calibration.csv, loy_compartment_calls.csv,
#         loy_patient_table.csv), the GSE269669 author objects and the Zeng pack Lambo object (re-projection), the
#         BoneMarrowMap scaffold, the gene-set sources named in PREREGISTRATION_H2.md section 5
# output: results/tables/loy/04_h2/
#           stage0_k0_reproduction.csv       K0: this script's LEUK counts equal 02's loy_compartment_calls.csv for every H2 library
#           stage0_counts.csv                per contrast library x fine label: suspected / retained / matched counts, depth medians, eligible
#           stage0_cells_<case>.csv.gz       per cell of the contrast library: call, group, matched, match_id (the handoff to 05)
#           stage0_pt43_bmm_percell.csv.gz   PT43A/B cells: this project's BMM fine / broad label next to the authors' Azimuth l2 label
#           stage0_k3_label_independence.csv GATE_H2_K3: fine labels with vs without chrY genes in the projection
#           stage0_gates.csv                 GATE_H2_0 and GATE_H2_K3 verdicts
#           h2_denominator.csv               male AML patients with >= 2 readable leukaemic timepoints, -Y class per timepoint
#         scripts/loy/gene_sets_h2_v1.tsv    the frozen feature panel (md5 printed; it goes into the tag message)
# what it does: Stage 0 of the H2 pre-registration -- counts and labels only; no expression value is read for the contrast.
#   [1] loads the per-cell tables of the H2 libraries; [2] re-derives the frozen three-state call and checks it against 02 (K0);
#   [3] leukaemic cells = the depositors' call, technical-zero control = the library's own T cells; [4] GATE_H2_K3 and the PT43
#   re-projection with this project's BoneMarrowMap pipeline (one label system for all cases; Azimuth l2 kept as secondary);
#   [5] depth matching of suspected_LOY to Y_retained cells within fine label (greedy nearest neighbour on log10 nCount, up to
#   LOY_H2_MATCH_RATIO per cell, caliper LOY_H2_CALIPER) and the eligibility floor (GATE_H2_0); [6] freezes the feature panel;
#   [7] the denominator table.
# Usage : sbatch scripts/loy/submit_loy_h2.sbatch 04 [--no_pt43] [--k3=auto|yes|no]
#   --no_pt43  skip the GSE269669 re-projection (Azimuth labels are then used for PT43; for a quick re-run of the counts)
#   --k3=      auto (default): re-project Lambo subsamples only if the scaffold's variable genes contain a chrY gene
# WHY: the prereg (sections 2, 4, 7) fixes the eligibility floor, the control and the gates before 05 reads expression;
#   every number here is a count or a label, which the prereg allows before the tag prereg-loy-h2-v1.
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "loy", "config_loy.R"))
source(here::here("scripts", "loy", "00_loy_common.R"))
args    <- commandArgs(trailingOnly = TRUE)
no_pt43 <- any(args == "--no_pt43")
k3_mode <- sub("^--k3=", "", grep("^--k3=", args, value = TRUE)); if (!length(k3_mode)) k3_mode <- "auto"
stopifnot(k3_mode %in% c("auto", "yes", "no"))
dir.create(LOY_H2_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(LOY_H2_SEED)
t0 <- Sys.time()

## [1] per-cell tables of the H2 libraries -----------------------------------------------------------
cases <- copy(LOY_H2_CASES)
dsl   <- unique(cases$dataset)
cells <- rbindlist(lapply(dsl, function(ds) {
  f <- file.path(LOY_PERCELL_DIR, paste0(ds, "__y_percell.csv.gz"))
  if (!file.exists(f)) { message("[1] WARNING missing per-cell file, dataset dropped: ", ds); return(NULL) }
  d <- fread(f)
  d[, Sample := fifelse(is.na(Sample) | Sample == "", as.character(object), as.character(Sample))]   # as in 02 [1]
  d[Sample %in% cases[dataset == ds, Sample]]
}), fill = TRUE)
ycols <- setdiff(grep("^y_[A-Z0-9]+$", names(cells), value = TRUE), "y_umi")
cells[, n_ygenes := rowSums(as.matrix(.SD) > 0, na.rm = TRUE), .SDcols = ycols]      # genes with >= 1 UMI, as in 02
cells[, `:=`(nCount_RNA = as.numeric(nCount_RNA), y_umi = as.numeric(y_umi), log_depth = log10(as.numeric(nCount_RNA)))]
cells <- merge(cells, cases[, .(dataset, Sample, case, role, timepoint)], by = c("dataset", "Sample"))
message(sprintf("[1] %d cells, %d libraries, %d cases (%.1f min)", nrow(cells), uniqueN(cells$Sample), uniqueN(cells$case),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(cells[, .(cells = .N, objects = paste(unique(object), collapse = "+")), by = .(case, Sample, role)][order(case, role)])
dup <- cells[, .(n_dup = sum(duplicated(cell))), by = dataset]; print(dup)
stopifnot(all(dup$n_dup == 0))                                        # a barcode in two author objects would be counted twice

## [2] frozen three-state call, re-derived; K0 = it reproduces 02 -------------------------------------
cells[, chem := LOY_EXT_ALL$chem_lit[match(dataset, LOY_EXT_ALL$dataset)]]        # external datasets: literature chemistry (7b)
cells[, comp := loy_comp_of(hierarchy_bin, bmm_fine)]
rates <- fread(file.path(LOY_CALL_DIR, "loy_reference_rates.csv"))
calib <- fread(file.path(LOY_CALL_DIR, "loy_external_calibration.csv"))
cells <- merge(cells, rates[, .(chem, comp, r)], by = c("chem", "comp"), all.x = TRUE)
cells <- merge(cells, unique(rates[, .(chem, r_chem)]), by = "chem", all.x = TRUE)
cells[is.na(r), r := r_chem]
cells <- merge(cells, calib[, .(dataset, factor, calib_status)], by = "dataset", all.x = TRUE)
stopifnot(!anyNA(cells$factor), !anyNA(cells$r))
cells[, lambda := nCount_RNA * r * factor]
cells[, call := loy_cell_call(y_umi, n_ygenes, lambda)]
k0 <- cells[comp == "LEUK", .(n = .N, n_susp = sum(call == "suspected_LOY"), n_ret = sum(call == "Y_retained")), by = .(dataset, Sample)]
ref02 <- fread(file.path(LOY_CALL_DIR, "loy_compartment_calls.csv"))[comp == "LEUK", .(dataset, Sample, n_02 = n, n_susp_02 = n_susp, n_ret_02 = n_ret)]
k0 <- merge(k0, ref02, by = c("dataset", "Sample"), all.x = TRUE)
k0[, pass := n == n_02 & n_susp == n_susp_02 & n_ret == n_ret_02]
fwrite(k0[order(dataset, Sample)], file.path(LOY_H2_DIR, "stage0_k0_reproduction.csv"))
message("[2] K0 reproduction of 02's LEUK counts:"); print(k0)
if (!all(k0$pass, na.rm = TRUE)) stop("[K0] re-derived LEUK counts differ from 02 -- the call is not the frozen one; stop")

## [3] leukaemic cells and the technical-zero control -----------------------------------------------
cells[, leuk := loy_h2_is_leuk(dataset, object, anno_malignant)]
cells[, fine_depositor := fifelse(is.na(bmm_fine) | bmm_fine == "", "unlabelled", as.character(bmm_fine))]
cells[, fine := fine_depositor]
cells[, ctrl := comp == "T" & !leuk]                                   # T cells: Y-zero there is dropout, not biology
message("[3] leukaemic / control cells per library:"); print(cells[, .(n = .N, leuk = sum(leuk), ctrl_T = sum(ctrl), leuk_susp = sum(leuk & call == "suspected_LOY"),
                                                                       leuk_ret = sum(leuk & call == "Y_retained")), by = .(case, Sample, role)][order(case, role)])

## [4] GATE_H2_K3 (labels independent of chrY) and the PT43 re-projection ------------------------------
gene_order <- fread(LOY_H2_GENE_ORDER, header = FALSE, col.names = c("gene", "chr", "start", "end"))
chrY_genes <- union(LOY_YGENES, gene_order[chr == "chrY", gene])
message(sprintf("[4] chrY exclusion list: %d genes (%d in the gene order, %d LOY count genes)", length(chrY_genes), sum(gene_order$chr == "chrY"), length(LOY_YGENES)))
k3 <- list()
need_proj <- (!no_pt43 && "GSE269669" %in% cells$dataset) || k3_mode == "yes"
if (need_proj || k3_mode != "no") {
  suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(BoneMarrowMap) })
  source(here::here("scripts", "config", "config_hierarchy.R"))
  if (!exists("BMM_KNN_K")) BMM_KNN_K <- 30L
  ref <- readRDS(BMM_SCAFFOLD_RDS)
  vg  <- ref$vargenes
  vg  <- if (is.data.frame(vg)) as.character(vg[[intersect(c("symbol", "gene", names(vg)), names(vg))[1]]]) else as.character(vg)
  vg_y <- intersect(vg, chrY_genes)
  message(sprintf("    scaffold variable genes: %d; chrY among them: %d (%s)", length(vg), length(vg_y), paste(vg_y, collapse = ",")))
  k3[["scaffold_vargenes"]] <- data.table(check = "scaffold_vargenes_chrY", library = "BoneMarrowMap", n = length(vg), value = length(vg_y),
                                          detail = paste(vg_y, collapse = ","), pass = length(vg_y) == 0)
  pick_col <- function(df, exact, pattern) { for (e in exact) if (e %in% names(df)) return(e); h <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE); if (length(h)) h[1] else NA_character_ }
  project <- function(cnt, md, tag) {
    q <- map_Query(exp_query = cnt, metadata_query = md, ref_obj = ref, do_normalize = TRUE, do_umap = FALSE, verbose = FALSE)
    q <- calculate_MappingError(q, reference = ref, MAD_threshold = MAD_THRESHOLD_BMM)
    q <- predict_CellTypes(query_obj = q, ref_obj = ref, ref_label = BMM_REF_LABEL, k = BMM_KNN_K, include_broad = TRUE)
    qm <- q@meta.data
    out <- data.table(cell = rownames(qm),
                      bmm_fine  = as.character(qm[[pick_col(qm, c("predicted_CellType", "final_CellType"), "^predicted_CellType$")]]),
                      bmm_broad = as.character(qm[[pick_col(qm, c("predicted_CellType_Broad", "predicted_Broad"), "[Bb]road")]]),
                      bmm_prob  = as.numeric(qm[[pick_col(qm, c("predicted_CellType_prob", "predicted_CellType_Broad_prob"), "predicted.*prob|_prob$")]]),
                      mapping_error_QC = as.character(qm[[pick_col(qm, c("mapping_error_QC", "mapping_QC"), "mapping.*QC")]]))
    setnames(out, c("bmm_fine", "bmm_broad", "bmm_prob", "mapping_error_QC"), paste0(c("bmm_fine", "bmm_broad", "bmm_prob", "mapping_error_QC"), tag))
    out
  }
  bin_map <- fread(LOY_BMM_BIN_MAP)
  # 4a. PT43: project once with every gene (the label used downstream) and once without chrY genes (K3)
  if (!no_pt43 && "GSE269669" %in% cells$dataset) {
    pt_cells <- cells[dataset == "GSE269669", cell]; mats <- list()
    for (obj in c("avm_aml", "avm_tme")) {
      s <- read_gse269669(obj); keep <- intersect(colnames(s), pt_cells)
      message(sprintf("    %s: %d of %d cells belong to the H2 libraries", obj, length(keep), ncol(s)))
      s <- subset(s, cells = keep)
      mats[[obj]] <- counts_of(s)
      rm(s); invisible(gc())
    }
    g   <- Reduce(intersect, lapply(mats, rownames))
    cnt <- do.call(cbind, lapply(mats, function(m) m[g, , drop = FALSE])); rm(mats)
    # metadata rownames must equal the cell names exactly: rbind() of a NAMED list prefixes row names ("avm_aml.<cell>"),
    # which made map_Query fail with "Cannot add new cells" (job 3689457); build it from the matrix columns instead
    md  <- data.frame(cell = colnames(cnt), row.names = colnames(cnt), stringsAsFactors = FALSE)
    md$object <- cells$object[match(md$cell, cells$cell)]
    message(sprintf("    PT43 count matrix: %d genes x %d cells", nrow(cnt), ncol(cnt)))
    pa <- project(cnt, md, "")
    pb <- project(cnt[!rownames(cnt) %in% chrY_genes, , drop = FALSE], md, "_noY")
    pt <- merge(pa, pb, by = "cell"); pt <- merge(pt, cells[dataset == "GSE269669", .(cell, Sample, object, leuk, azimuth_l2 = fine_depositor)], by = "cell")
    pt <- merge(pt, bin_map[, .(bmm_broad = CellType_Broad, hierarchy_bin_bmm = hierarchy_bin)], by = "bmm_broad", all.x = TRUE)
    fwrite(pt, file.path(LOY_H2_DIR, "stage0_pt43_bmm_percell.csv.gz"))
    for (lib in unique(pt$Sample)) for (sub in c("all", "leuk")) {
      d <- if (sub == "all") pt[Sample == lib] else pt[Sample == lib & leuk == TRUE]
      k3[[paste(lib, sub)]] <- data.table(check = paste0("reprojection_", sub), library = lib, n = nrow(d), value = mean(d$bmm_fine == d$bmm_fine_noY),
                                          detail = sprintf("broad agreement %.3f", mean(d$bmm_broad == d$bmm_broad_noY)), pass = mean(d$bmm_fine == d$bmm_fine_noY) >= LOY_H2_K3_MIN_AGREE)
    }
    message("    PT43 BMM fine label vs Azimuth l2 (leukaemic cells, PT43A):")
    print(pt[Sample == "PT43A" & leuk == TRUE, .N, by = .(bmm_fine, azimuth_l2)][order(-N)][1:15])
    # downstream fine label for PT43 = this project's BMM label (same system as the Lambo pack); Azimuth l2 stays as fine_depositor
    cells[pt, on = "cell", fine := i.bmm_fine]
    rm(cnt, pa, pb); invisible(gc())
  }
  # 4b. Lambo: by construction if the scaffold's variable genes hold no chrY gene; otherwise re-project a subsample of each contrast library
  if (k3_mode == "yes" || (k3_mode == "auto" && length(vg_y) > 0)) {
    f <- file.path(LOY_ZENG_RDS_DIR, paste0(LOY_EXT_STUDIES[study == "Lambo_2023", file_stem], "_projected_cleaned.rds"))
    message(sprintf("    Lambo re-projection needed; reading %s (%.1f GB)", basename(f), file.size(f) / 1e9))
    s <- readRDS(f)
    for (lib in cases[dataset == "ZENG_Lambo_2023" & role == "contrast", Sample]) {
      cc <- intersect(colnames(s), cells[Sample == lib, cell]); cc <- sample(cc, min(LOY_H2_K3_SUBSAMPLE, length(cc)))
      sub <- subset(s, cells = cc); cnt <- counts_of(sub); rm(sub)
      md  <- data.frame(cell = cc, row.names = cc, stringsAsFactors = FALSE)
      pa  <- project(cnt, md, ""); pb <- project(cnt[!rownames(cnt) %in% chrY_genes, , drop = FALSE], md, "_noY")
      d   <- merge(merge(pa, pb, by = "cell"), cells[Sample == lib, .(cell, fine_depositor)], by = "cell")
      k3[[paste(lib, "sub")]] <- data.table(check = "reprojection_subsample", library = lib, n = nrow(d), value = mean(d$bmm_fine == d$bmm_fine_noY),
                                            detail = sprintf("agreement with the pack label %.3f", mean(d$bmm_fine == d$fine_depositor)),
                                            pass = mean(d$bmm_fine == d$bmm_fine_noY) >= LOY_H2_K3_MIN_AGREE)
      message(sprintf("    %s: %d cells, fine agreement with/without chrY %.3f, with the pack label %.3f", lib, nrow(d), mean(d$bmm_fine == d$bmm_fine_noY), mean(d$bmm_fine == d$fine_depositor)))
    }
    rm(s); invisible(gc())
  } else {
    for (lib in cases[dataset == "ZENG_Lambo_2023" & role == "contrast", Sample])
      k3[[paste(lib, "constr")]] <- data.table(check = "by_construction", library = lib, n = NA_integer_, value = NA_real_,
                                               detail = "scaffold variable genes contain no chrY gene", pass = length(vg_y) == 0)
  }
  rm(ref); invisible(gc())
}
k3 <- rbindlist(k3, fill = TRUE)
if (nrow(k3)) { fwrite(k3, file.path(LOY_H2_DIR, "stage0_k3_label_independence.csv")); message("[4] K3:"); print(k3) }

## [5] depth matching within fine label; eligibility ----------------------------------------------------
# match_depth(): greedy nearest neighbour on log10 depth, without replacement, up to `ratio` partners within `caliper`.
# findInterval() gives the position of the nearest lower value in the sorted vector; the loop then walks outward.
match_depth <- function(ld_s, ld_r, ratio, caliper) {
  ord <- order(ld_r); rs <- ld_r[ord]; used <- logical(length(rs)); out <- vector("list", length(ld_s))
  for (i in seq_along(ld_s)) {
    pos <- findInterval(ld_s[i], rs); lo <- pos; hi <- pos + 1L; picked <- integer()
    while (length(picked) < ratio) {
      while (lo >= 1L && used[lo]) lo <- lo - 1L
      while (hi <= length(rs) && used[hi]) hi <- hi + 1L
      dlo <- if (lo >= 1L) ld_s[i] - rs[lo] else Inf
      dhi <- if (hi <= length(rs)) rs[hi] - ld_s[i] else Inf
      if (min(dlo, dhi) > caliper) break
      if (dlo <= dhi) { picked <- c(picked, lo); used[lo] <- TRUE; lo <- lo - 1L } else { picked <- c(picked, hi); used[hi] <- TRUE; hi <- hi + 1L }
    }
    out[[i]] <- ord[picked]
  }
  out
}
match_stratum <- function(d) {        # d: one library x stratum; group A = the Y-zero side, B = Y_retained
  d <- copy(d); d[, `:=`(matched = FALSE, match_id = NA_integer_)]
  a <- which(d$group == "A"); b <- which(d$group == "B")
  if (!length(a) || !length(b)) return(d)
  a <- a[sample.int(length(a))]                        # random processing order: no cell monopolises the nearest partners
  res <- match_depth(d$log_depth[a], d$log_depth[b], LOY_H2_MATCH_RATIO, LOY_H2_CALIPER)
  for (i in seq_along(a)) if (length(res[[i]])) set(d, i = c(a[i], b[res[[i]]]), j = c("matched", "match_id"), value = list(TRUE, i))
  d
}
summarise_stratum <- function(d) d[, .(n_A = sum(group == "A"), n_B = sum(group == "B"),
                                        matched_A = sum(matched & group == "A"), matched_B = sum(matched & group == "B"),
                                        med_depth_A = as.numeric(median(nCount_RNA[group == "A"])), med_depth_B = as.numeric(median(nCount_RNA[group == "B"])),
                                        med_depth_A_m = as.numeric(median(nCount_RNA[matched & group == "A"])), med_depth_B_m = as.numeric(median(nCount_RNA[matched & group == "B"])))]
counts <- list(); handoff <- list()
for (lib in cases[role == "contrast", Sample]) {
  cs <- cases[Sample == lib, case]
  # leukaemic contrast: A = suspected_LOY, B = Y_retained, uncertain dropped; strata = fine label, plus the pooled "ALL"
  L <- cells[Sample == lib & leuk & call != "uncertain", .(case = cs, Sample, cell, object, fine, fine_depositor, comp, nCount_RNA, log_depth, lambda, y_umi, n_ygenes, call,
                                                           arm = "leukaemic", group = fifelse(call == "suspected_LOY", "A", "B"))]
  # technical-zero control: the library's T cells, A = Y-zero at any depth, B = Y_retained (the pure-dropout version of the same split)
  Ct <- cells[Sample == lib & ctrl & (y_umi == 0 | call == "Y_retained"), .(case = cs, Sample, cell, object, fine, fine_depositor, comp, nCount_RNA, log_depth, lambda, y_umi, n_ygenes, call,
                                                                           arm = "control_T", group = fifelse(y_umi == 0, "A", "B"))]
  both <- rbind(L, Ct)
  res  <- rbindlist(lapply(split(both, by = c("arm", "fine"), drop = TRUE), match_stratum))
  pooled <- rbindlist(lapply(split(both, by = "arm", drop = TRUE), match_stratum)); pooled[, fine := "ALL"]
  res_all <- rbind(res, pooled)
  sm <- res_all[, summarise_stratum(.SD), by = .(case, Sample, arm, fine)]
  sm[, eligible := matched_A >= LOY_H2_MIN_CELLS & matched_B >= LOY_H2_MIN_CELLS]
  counts[[lib]] <- sm; handoff[[lib]] <- res_all
  message(sprintf("[5] %s (%s): strata eligible leukaemic %d, control pooled A/B matched %d/%d", lib, cs, sum(sm[arm == "leukaemic" & fine != "ALL", eligible]),
                  sm[arm == "control_T" & fine == "ALL", matched_A], sm[arm == "control_T" & fine == "ALL", matched_B]))
  print(sm[arm == "leukaemic"][order(-n_A)][, .(fine, n_A, n_B, matched_A, matched_B, med_depth_A, med_depth_B, med_depth_A_m, med_depth_B_m, eligible)])
  fwrite(res_all, file.path(LOY_H2_DIR, sprintf("stage0_cells_%s.csv.gz", cs)))
}
counts <- rbindlist(counts); fwrite(counts, file.path(LOY_H2_DIR, "stage0_counts.csv"))
# sanity: matched halves sit at the same depth (medians within 5 %), and no cell carries two groups
chk <- counts[eligible == TRUE, .(case, Sample, arm, fine, depth_ratio = med_depth_A_m / med_depth_B_m)]
message("[5] depth ratio A/B among eligible matched strata (should be within 0.95-1.05):"); print(chk)
stopifnot(all(abs(log(chk$depth_ratio)) <= log(1.05)))
elig_cases <- counts[arm == "leukaemic" & fine != "ALL" & eligible == TRUE, unique(case)]
ctrl_ok    <- counts[arm == "control_T" & fine == "ALL" & eligible == TRUE, unique(case)]

## [6] feature panel, frozen ---------------------------------------------------------------------------
st <- fread(text = grep("^#", readLines(LOY_H2_STEMNESS_TSV), value = TRUE, invert = TRUE))
gs <- list()
for (nm in names(LOY_H2_STEMNESS)) gs[[nm]] <- st[signature == LOY_H2_STEMNESS[[nm]], .(set = nm, role = "primary", gene, weight = as.numeric(weight),
                                                                                           source = "scripts/以前06_hierarchy/stemness_signatures.tsv (as frozen in PREREGISTRATION_CLONAL_ENV.md)")]
cc <- Seurat::cc.genes.updated.2019
gs[["F4_S"]]   <- data.table(set = "F4_S",   role = "primary", gene = cc$s.genes,   weight = NA_real_, source = "Seurat cc.genes.updated.2019 (Tirosh 2016)")
gs[["F4_G2M"]] <- data.table(set = "F4_G2M", role = "primary", gene = cc$g2m.genes, weight = NA_real_, source = "Seurat cc.genes.updated.2019 (Tirosh 2016)")
stopifnot(file.exists(LOY_H2_GMT))
gmt <- strsplit(readLines(LOY_H2_GMT), "\t"); names(gmt) <- vapply(gmt, `[`, "", 1L); gmt <- lapply(gmt, function(x) x[-(1:2)])
for (nm in names(LOY_H2_HALLMARK)) { stopifnot(LOY_H2_HALLMARK[[nm]] %in% names(gmt))
  gs[[nm]] <- data.table(set = nm, role = "primary", gene = gmt[[LOY_H2_HALLMARK[[nm]]]], weight = NA_real_, source = paste0(LOY_H2_HALLMARK[[nm]], " ", basename(LOY_H2_GMT))) }
mp <- fread(LOY_H2_MP_TSV)
for (m in LOY_H2_MPS) gs[[m]] <- mp[MP == m, .(set = m, role = fifelse(m == "MP9", "artefact_control", "primary"), gene, weight = as.numeric(weight),
                                               source = "this project cNMF results/tables/04_cnmf/malignant_metaprograms.tsv (2026-07-17)")]
gs[["F12_chrY"]]   <- data.table(set = "F12_chrY", role = "machinery", gene = LOY_YGENES, weight = NA_real_, source = "config_loy.R LOY_YGENES (the counted chrY genes)")
gs[["chrY_excl"]]  <- data.table(set = "chrY_exclude", role = "exclude", gene = chrY_genes, weight = NA_real_, source = "gencode_GRCh38_gene_order.txt chrY + LOY_YGENES")
gs <- rbindlist(gs)
leak <- gs[role %in% c("primary", "artefact_control") & gene %in% chrY_genes]
if (nrow(leak)) { message("[6] chrY genes removed from feature sets: "); print(leak); gs <- gs[!(role %in% c("primary", "artefact_control") & gene %in% chrY_genes)] }
sizes <- gs[, .N, by = .(set, role)]; print(sizes); stopifnot(all(sizes$N >= 5))
fwrite(gs, LOY_H2_GENESETS, sep = "\t")
md5 <- tools::md5sum(LOY_H2_GENESETS)
message(sprintf("[6] %s written: %d sets, %d rows, md5 %s", basename(LOY_H2_GENESETS), uniqueN(gs$set), nrow(gs), md5))

## [7] denominator: male AML patients with >= 2 readable leukaemic timepoints -------------------------------
pt   <- fread(file.path(LOY_CALL_DIR, "loy_patient_table.csv"))
cc02 <- fread(file.path(LOY_CALL_DIR, "loy_compartment_calls.csv"))[comp == "LEUK" & readable == TRUE & !is.na(f_hat)]
pack <- unique(fread(file.path(LOY_ZENG_RDS_DIR, "scAML_21cohorts_Annotations_byCell.csv"), select = c("SampleID_Unique", "PatientID", "Timepoint")))
patient_key <- function(ds, s) {                     # one rule per dataset; the pack carries its own PatientID
  fifelse(grepl("^ZENG_", ds), pack$PatientID[match(s, pack$SampleID_Unique)],
  fifelse(ds %in% c("GSE269669", "GSE239721"), sub("[A-Z]$", "", s),
  fifelse(ds == "GSE116256", sub("-D.*$", "", s),
  fifelse(ds == "GSE147989", substr(s, 1, 3),
  fifelse(ds == "GSE201966", sub("M.*$", "", sub("^GSM[0-9]+_", "", s)),
  fifelse(ds == "GSE207356", sub("_.*$", "", sub("^GSM[0-9]+_", "", s)),
  fifelse(ds == "GSE227903", sub("_.*$", "", s), s)))))))
}
den <- merge(cc02[, .(dataset, Sample, n_LEUK = n, f_hat_LEUK = f_hat)],
             pt[sex_group == "M" & is_aml == TRUE, .(dataset, Sample, leuk_class, timepoint_man, age_group, AgeGroup)], by = c("dataset", "Sample"))
den[, patient := patient_key(dataset, Sample)]
den[, timepoint := fifelse(!is.na(timepoint_man) & timepoint_man != "", timepoint_man, pack$Timepoint[match(Sample, pack$SampleID_Unique)])]
den[, age := fifelse(!is.na(age_group) & age_group != "", age_group, AgeGroup)]
dup_studies <- c("ZENG_Ennis_2023", "ZENG_Zeng_2025", "ZENG_vanGalen_2019", "ZENG_Petti_2019")   # pack copies of on-disk datasets
den <- den[!dataset %in% dup_studies]
den <- den[, if (.N >= 2) .(age = age[1], n_readable_leuk_timepoints = .N, libraries = paste(Sample, collapse = ";"), timepoints = paste(timepoint, collapse = ";"),
                            f_hat_LEUK = paste(sprintf("%.3f", f_hat_LEUK), collapse = ";"), n_LEUK = paste(n_LEUK, collapse = ";"),
                            leuk_class = paste(leuk_class, collapse = ";"),
                            status_change = any(leuk_class %in% c("clonal_LOY", "subclonal_LOY")) & uniqueN(leuk_class) > 1), by = .(dataset, patient)]
fwrite(den[order(dataset, patient)], file.path(LOY_H2_DIR, "h2_denominator.csv"))
message(sprintf("[7] denominator: %d patients with >= 2 readable leukaemic timepoints (paediatric %d); %d change -Y class between timepoints",
                nrow(den), sum(den$age %in% c("Pediatric", "infants", "teens", "pediatric_NA", "children")), sum(den$status_change)))
print(den[status_change == TRUE, .(dataset, patient, timepoints, f_hat_LEUK, leuk_class)])

## [8] gates -----------------------------------------------------------------------------------------
gates <- data.table(gate = c("GATE_H2_0", "GATE_H2_K3"),
                    value = c(length(elig_cases), if (nrow(k3)) mean(k3$pass) else NA_real_),
                    detail = c(sprintf("eligible cases: %s | control ok: %s", paste(elig_cases, collapse = ","), paste(ctrl_ok, collapse = ",")),
                               if (nrow(k3)) paste(k3[, sprintf("%s:%s=%s", library, check, ifelse(pass, "pass", "FAIL"))], collapse = " | ") else "not run"),
                    pass = c(length(elig_cases) >= 2, if (nrow(k3)) all(k3$pass) else NA))
fwrite(gates, file.path(LOY_H2_DIR, "stage0_gates.csv")); print(gates)
message(sprintf("[done] Stage 0 in %.1f min; gene-set md5 %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), md5))
