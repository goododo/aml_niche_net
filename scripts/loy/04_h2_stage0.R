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
#         scripts/loy/gene_sets_h2_v2.tsv    the frozen feature panel (md5 printed; it goes into the tag message); v1 lacked the PAR exclusion
# what it does: Stage 0 of the H2 pre-registration -- counts and labels only; no expression value is read for the contrast.
#   [1] loads the per-cell tables of the H2 libraries; [2] re-derives the frozen three-state call and checks it against 02 (K0);
#   [3] leukaemic cells = the depositors' call, technical-zero control = the library's own T cells; [4] GATE_H2_K3 and the PT43
#   re-projection with this project's BoneMarrowMap pipeline (one label system for all cases; Azimuth l2 kept as secondary);
#   [5] depth matching of suspected_LOY to Y_retained cells within fine label (greedy nearest neighbour on log10 nCount, up to
#   LOY_H2_MATCH_RATIO per cell, caliper LOY_H2_CALIPER) and the eligibility floor (GATE_H2_0); [6] freezes the feature panel;
#   [7] the denominator table.
# Usage : sbatch scripts/loy/submit_loy_h2.sbatch 04 [--no_pt43] [--k3=auto|yes|no] [--force]
#   --force    redo the PT43 / Lambo projections even when stage0_pt43_bmm_percell.csv.gz / stage0_k3_lambo_subsample.csv.gz exist
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
force   <- any(args == "--force")                                     # redo the projections even when their cached outputs exist
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
# technical-zero control (prereg section 2 / amendment 3, 8): non-leukaemic cells of compartments whose sample-level f_hat is
# <= LOY_H2_CONTROL_MAX_F in 02's table (T, NK, B, erythroid, normal Mono/DC, megakaryocyte); Y-zero there is dropout, not biology.
cc02 <- fread(file.path(LOY_CALL_DIR, "loy_compartment_calls.csv"))
ctrl_ok_comps <- cc02[readable == TRUE & !is.na(f_hat) & f_hat <= LOY_H2_CONTROL_MAX_F &
                      comp %in% c("T", "NK", "B_Plasma", "Erythroid", "Mono_DC", "Megakaryocyte"), paste(Sample, comp)]
cells[, ctrl := !leuk & paste(Sample, comp) %in% ctrl_ok_comps]
message("[3] leukaemic / control cells per library:"); print(cells[, .(n = .N, leuk = sum(leuk), ctrl = sum(ctrl), ctrl_y0 = sum(ctrl & y_umi == 0), leuk_susp = sum(leuk & call == "suspected_LOY"),
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
  # informational row: a chrY gene among the variable genes only decides whether the Lambo subsample is re-projected (pass = NA)
  k3[["scaffold_vargenes"]] <- data.table(check = "scaffold_vargenes_chrY", library = "BoneMarrowMap", n = length(vg), value = length(vg_y),
                                          detail = paste0(paste(vg_y, collapse = ","), if (length(vg_y)) " -> Lambo subsample re-projected" else " -> labels chrY-free by construction"), pass = NA)
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
  pt_file <- file.path(LOY_H2_DIR, "stage0_pt43_bmm_percell.csv.gz")
  if (!no_pt43 && "GSE269669" %in% cells$dataset && file.exists(pt_file) && !force) {
    message("    PT43 projection cached: ", basename(pt_file)); pt <- fread(pt_file)
  } else if (!no_pt43 && "GSE269669" %in% cells$dataset) {
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
    fwrite(pt, pt_file)
    rm(cnt, pa, pb); invisible(gc())
  }
  if (exists("pt")) {
    # agreement over cells that received a label in both runs (mapping-QC failures carry NA and are counted separately)
    for (lib in unique(pt$Sample)) for (sub in c("all", "leuk")) {
      d <- if (sub == "all") pt[Sample == lib] else pt[Sample == lib & leuk == TRUE]
      ok <- !is.na(d$bmm_fine) & !is.na(d$bmm_fine_noY); agree <- mean(d$bmm_fine[ok] == d$bmm_fine_noY[ok])
      k3[[paste(lib, sub)]] <- data.table(check = paste0("reprojection_", sub), library = lib, n = sum(ok), value = agree,
                                          detail = sprintf("broad agreement %.4f; %d mapping-QC-fail cells excluded", mean(d$bmm_broad[ok] == d$bmm_broad_noY[ok]), sum(!ok)),
                                          pass = agree >= LOY_H2_K3_MIN_AGREE)
    }
    message("    PT43 BMM fine label vs Azimuth l2 (leukaemic cells, PT43A):")
    print(pt[Sample == "PT43A" & leuk == TRUE, .N, by = .(bmm_fine, azimuth_l2)][order(-N)][1:15])
    # downstream fine label for PT43 = this project's BMM label (same system as the Lambo pack); Azimuth l2 stays as fine_depositor
    cells[pt[!is.na(bmm_fine)], on = "cell", fine := i.bmm_fine]          # QC-fail cells keep the depositor label
  }
  # 4b. Lambo: by construction if the scaffold's variable genes hold no chrY gene; otherwise re-project a subsample of each contrast library
  lam_file <- file.path(LOY_H2_DIR, "stage0_k3_lambo_subsample.csv.gz")
  if (k3_mode == "yes" || (k3_mode == "auto" && length(vg_y) > 0)) {
    if (file.exists(lam_file) && !force) { message("    Lambo K3 subsample cached: ", basename(lam_file)); lam <- fread(lam_file) } else {
      f <- file.path(LOY_ZENG_RDS_DIR, paste0(LOY_EXT_STUDIES[study == "Lambo_2023", file_stem], "_projected_cleaned.rds"))
      message(sprintf("    Lambo re-projection needed; reading %s (%.1f GB)", basename(f), file.size(f) / 1e9))
      s <- readRDS(f); lam <- list()
      for (lib in cases[dataset == "ZENG_Lambo_2023" & role == "contrast", Sample]) {
        cc <- intersect(colnames(s), cells[Sample == lib, cell]); cc <- sample(cc, min(LOY_H2_K3_SUBSAMPLE, length(cc)))
        sub <- subset(s, cells = cc); cnt <- counts_of(sub); rm(sub)
        md  <- data.frame(cell = cc, row.names = cc, stringsAsFactors = FALSE)
        pa  <- project(cnt, md, ""); pb <- project(cnt[!rownames(cnt) %in% chrY_genes, , drop = FALSE], md, "_noY")
        lam[[lib]] <- merge(merge(pa, pb, by = "cell"), cells[Sample == lib, .(cell, Sample, fine_depositor)], by = "cell")
      }
      rm(s); invisible(gc()); lam <- rbindlist(lam); fwrite(lam, lam_file)
    }
    for (lib in unique(lam$Sample)) {
      d <- lam[Sample == lib]; ok <- !is.na(d$bmm_fine) & !is.na(d$bmm_fine_noY); agree <- mean(d$bmm_fine[ok] == d$bmm_fine_noY[ok])
      k3[[paste(lib, "sub")]] <- data.table(check = "reprojection_subsample", library = lib, n = sum(ok), value = agree,
                                            detail = sprintf("agreement with the pack label %.3f; %d mapping-QC-fail cells excluded", mean(d$bmm_fine[ok] == d$fine_depositor[ok]), sum(!ok)),
                                            pass = agree >= LOY_H2_K3_MIN_AGREE)
      message(sprintf("    %s: %d cells, fine agreement with/without chrY %.4f, with the pack label %.3f", lib, sum(ok), agree, mean(d$bmm_fine[ok] == d$fine_depositor[ok])))
    }
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
counts <- list(); handoff <- list(); ctrl_note <- list()
cols <- function(d, cs, arm_name, grpA) d[, .(case = cs, Sample, cell, object, fine, fine_depositor, comp, nCount_RNA, log_depth, lambda, y_umi, n_ygenes, call, arm = arm_name, group = fifelse(grpA, "A", "B"))]
# control matching: strata = library x compartment, plus one pooled stratum per library ("ALL"; depth is library-specific, so pooling never crosses libraries)
match_control <- function(Ct) {
  if (!nrow(Ct)) return(Ct[, `:=`(matched = logical(0), match_id = integer(0))])
  Ct <- copy(Ct)[, fine := comp]
  per <- rbindlist(lapply(split(Ct, by = c("Sample", "fine"), drop = TRUE), match_stratum))
  pooled <- rbindlist(lapply(split(Ct, by = "Sample", drop = TRUE), match_stratum)); pooled[, fine := "ALL"]
  rbind(per, pooled)
}
for (lib in cases[role == "contrast", Sample]) {
  cs <- cases[Sample == lib, case]
  # leukaemic contrast: A = suspected_LOY, B = Y_retained, uncertain dropped; strata = fine label, plus the pooled "ALL"
  L <- cols(cells[Sample == lib & leuk & call != "uncertain"], cs, "leukaemic", cells[Sample == lib & leuk & call != "uncertain", call == "suspected_LOY"])
  resL <- rbindlist(lapply(split(L, by = "fine", drop = TRUE), match_stratum))
  pooledL <- match_stratum(L); pooledL[, fine := "ALL"]
  # technical-zero control: A = Y-zero at any depth, B = Y_retained, in the control compartments of the contrast library;
  # when fewer than LOY_H2_MIN_CELLS Y-zero control cells find a match there, the patient's other libraries join (amendment 8)
  sel <- cells$ctrl & (cells$y_umi == 0 | cells$call == "Y_retained")
  Ct  <- cols(cells[Sample == lib & sel], cs, "control", cells[Sample == lib & sel, y_umi == 0])
  resC <- match_control(Ct); src <- "contrast library"
  if (resC[fine == "ALL", sum(matched & group == "A")] < LOY_H2_MIN_CELLS) {
    Ct <- cols(cells[case == cs & sel], cs, "control", cells[case == cs & sel, y_umi == 0]); resC <- match_control(Ct); src <- "all libraries of the patient"
  }
  ctrl_note[[lib]] <- data.table(case = cs, control_source = src, control_libraries = paste(unique(resC$Sample), collapse = ";"),
                                 matched_A = resC[fine == "ALL", sum(matched & group == "A")], matched_B = resC[fine == "ALL", sum(matched & group == "B")])
  res_all <- rbind(resL, pooledL, resC)
  sm <- res_all[, summarise_stratum(.SD), by = .(case, Sample, arm, fine)]
  sm[, eligible := matched_A >= LOY_H2_MIN_CELLS & matched_B >= LOY_H2_MIN_CELLS]
  counts[[lib]] <- sm; handoff[[lib]] <- res_all
  message(sprintf("[5] %s (%s): eligible leukaemic strata %d | control (%s): pooled A/B matched %d/%d", lib, cs, sum(sm[arm == "leukaemic" & fine != "ALL", eligible]),
                  src, ctrl_note[[lib]]$matched_A, ctrl_note[[lib]]$matched_B))
  print(sm[arm == "leukaemic"][order(-n_A)][, .(fine, n_A, n_B, matched_A, matched_B, med_depth_A, med_depth_B, med_depth_A_m, med_depth_B_m, eligible)])
  print(sm[arm == "control"][order(Sample, -n_A)][, .(Sample, fine, n_A, n_B, matched_A, matched_B, med_depth_A_m, med_depth_B_m)])
  fwrite(res_all, file.path(LOY_H2_DIR, sprintf("stage0_cells_%s.csv.gz", cs)))
}
counts <- rbindlist(counts); fwrite(counts, file.path(LOY_H2_DIR, "stage0_counts.csv"))
ctrl_note <- rbindlist(ctrl_note); fwrite(ctrl_note, file.path(LOY_H2_DIR, "stage0_control_source.csv")); print(ctrl_note)
# sanity: matched halves sit at the same depth. Greedy 1:3 matching inside a 0.1 log10 caliper leaves median gaps of a few per cent
# (job 3689458: 0.91-1.07); the hard stop is half the caliper, and log depth stays a covariate in 05 (prereg section 4).
chk <- counts[eligible == TRUE, .(case, Sample, arm, fine, depth_ratio = med_depth_A_m / med_depth_B_m)]
message(sprintf("[5] depth ratio A/B among eligible matched strata (hard stop beyond 10^%.2f = %.2f-%.2f):", LOY_H2_CALIPER / 2, 10^(-LOY_H2_CALIPER / 2), 10^(LOY_H2_CALIPER / 2))); print(chk)
stopifnot(all(abs(log10(chk$depth_ratio)) <= LOY_H2_CALIPER / 2))
elig_cases <- counts[arm == "leukaemic" & fine != "ALL" & eligible == TRUE, unique(case)]
ctrl_ok    <- ctrl_note[matched_A >= LOY_H2_MIN_CELLS & matched_B >= LOY_H2_MIN_CELLS, case]

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
# v2: pseudoautosomal genes sit on chrX in the gene order but have one copy in a -Y cell (SLC25A6 and CD99 were the top two S1
# genes in every case of the v1 run); excluded from feature sets, kept in S1 as the dosage known answer (05 flags them)
par_genes <- gene_order[chr == "chrX" & ((start >= LOY_H2_PAR_GRCH38$PAR1[1] & end <= LOY_H2_PAR_GRCH38$PAR1[2]) |
                                        (start >= LOY_H2_PAR_GRCH38$PAR2[1] & end <= LOY_H2_PAR_GRCH38$PAR2[2])), gene]
message(sprintf("[6] PAR genes (GRCh38 chrX PAR1/PAR2): %d: %s", length(par_genes), paste(par_genes, collapse = ",")))
gs[["par_excl"]]   <- data.table(set = "PAR_exclude", role = "exclude", gene = par_genes, weight = NA_real_, source = "gencode_GRCh38_gene_order.txt chrX inside PAR1/PAR2 (GRCh38 coordinates)")
gs <- rbindlist(gs)
excl_genes <- union(chrY_genes, par_genes)
leak <- gs[role %in% c("primary", "artefact_control") & gene %in% excl_genes]
if (nrow(leak)) { message("[6] chrY / PAR genes removed from feature sets: "); print(leak); gs <- gs[!(role %in% c("primary", "artefact_control") & gene %in% excl_genes)] }
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
# raw-matrix series and ScPCA (config 7e/7f): patient and timepoint come from the 03b/03c scan summaries (2026-11-15)
scan_sum <- rbindlist(lapply(c(LOY_GEO_MTX$dataset, LOY_SCPCA_DATASET), function(ds) {
  f <- file.path(LOY_SCAN_DIR, paste0(ds, "__y_sample_summary.csv"))
  if (file.exists(f)) fread(f, select = c("dataset", "Sample", "Patient_ID", "Timepoint"), colClasses = "character") }), fill = TRUE)
if (nrow(scan_sum)) {
  den[scan_sum, on = .(dataset, Sample), `:=`(patient = i.Patient_ID, timepoint_scan = i.Timepoint)]   # update-join: only matching rows change
  den[is.na(timepoint) & !is.na(timepoint_scan), timepoint := timepoint_scan]; den[, timepoint_scan := NULL]
}
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
k3g <- k3[!is.na(pass)]                                   # the gate is the label-agreement rows only
gates <- data.table(gate = c("GATE_H2_0", "GATE_H2_K3"),
                    value = c(length(elig_cases), if (nrow(k3g)) min(k3g$value) else NA_real_),
                    detail = c(sprintf("eligible cases: %s | control ok: %s", paste(elig_cases, collapse = ","), paste(ctrl_ok, collapse = ",")),
                               if (nrow(k3g)) paste(k3g[, sprintf("%s:%s=%.4f", library, check, value)], collapse = " | ") else "not run"),
                    pass = c(length(elig_cases) >= 2, if (nrow(k3g)) all(k3g$pass) else NA))
fwrite(gates, file.path(LOY_H2_DIR, "stage0_gates.csv")); print(gates)
message(sprintf("[done] Stage 0 in %.1f min; gene-set md5 %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), md5))
