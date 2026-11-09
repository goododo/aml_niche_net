#!/usr/bin/env Rscript
# 02_label_barcodes.R ----
# Doublet-pairing SMOKE, stage 2: one compartment / state label and a cycling flag for every analysis barcode of a lane.
# Class A takes its projection label through dbl_label_map.tsv (+ the strict-T rule, proposal 3.4); classes B and C, which
# have no label anywhere, get one by kNN transfer from the same lane's class-A cells (handoff v1.1 4.4.4).
# INPUT  : DBL_LANE_DIR/<lane>.rds (stage 01); scripts/dbl/dbl_label_map.tsv
# OUTPUT : DBL_LANE_DIR/<lane>__labels.csv.gz  cell, class, donor, bmm_fine, compartment, source, vote, t_id_umi, phase, cycling
#          DBL_TAB_DIR/dbl_class_a_composition_smoke.csv   lane x donor x compartment counts of CLASS-A cells only (= singlet f)
#          DBL_TAB_DIR/dbl_transfer_check_smoke.csv        held-out class-A transfer accuracy per compartment (known answer)
#          DBL_TAB_DIR/dbl_checks_02_smoke.csv
# Usage  : Rscript scripts/dbl/02_label_barcodes.R [--lanes=single,pooled] [--force]
# Masking [D2]: nothing here types a doublet. Labels of B/C barcodes are written for stage 03 but their composition is not
#   printed or tabulated (B/C barcodes are where removed doublets sit).

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here); library(Seurat) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
ARGS  <- commandArgs(TRUE)
FORCE <- "--force" %in% ARGS
la    <- sub("^--lanes=", "", grep("^--lanes=", ARGS, value = TRUE))
LANES <- if (length(la)) DBL_SMOKE_LANES[strsplit(la, ",")[[1]]] else DBL_SMOKE_LANES
stopifnot(!anyNA(LANES))
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-78s expect %-12s got %-12s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
MAP <- fread(here::here("scripts", "dbl", DBL_LABEL_MAP_FILE))
stopifnot(!anyDuplicated(MAP$bmm_fine))

knn_vote <- function(ref_pc, ref_lab, qry_pc, k) {   # majority label among k nearest reference cells, and its vote share
  nn <- RANN::nn2(data = ref_pc, query = qry_pc, k = k)$nn.idx
  lab <- matrix(ref_lab[nn], nrow = nrow(nn))
  top <- apply(lab, 1, function(r) { t <- table(r); c(names(t)[which.max(t)], max(t) / k) })
  data.table(label = top[1, ], vote = as.numeric(top[2, ])) }

COMP <- list(); TRANS <- list()
for (L in LANES) {
  out <- file.path(DBL_LANE_DIR, paste0(L, "__labels.csv.gz"))
  if (file.exists(out) && !FORCE) { ts("[skip] ", out, " exists (use --force)"); next }
  ts("[1] lane ", L)
  x  <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds")))
  M  <- x$counts
  bc <- x$bc[in_counts == TRUE][match(colnames(M), cell)]
  stopifnot(identical(bc$cell, colnames(M)))
  ## -- [2] class A: projection label -> compartment; strict-T rule  [P 3.4, D5] ----
  A <- merge(bc[class == "A", .(cell)], x$ann, by = "cell", all.x = TRUE)
  A <- merge(A, MAP[, .(bmm_fine, compartment)], by = "bmm_fine", all.x = TRUE)
  chk(sprintf("%s: every non-missing class-A bmm_fine is in the label map", L), 0L, A[!is.na(bmm_fine) & is.na(compartment), .N])
  A[is.na(bmm_fine), compartment := unname(DBL_BIN_FALLBACK[hierarchy_bin])]   # projection failed: bin-level fallback
  ts(sprintf("    %d class-A cells without a projection label -> bin fallback", A[is.na(bmm_fine), .N]))
  chk(sprintf("%s: every class-A cell without bmm_fine has a bin fallback", L), 0L, A[is.na(compartment), .N])
  tid <- intersect(DBL_T_ID_GENES, rownames(M))
  chk(sprintf("%s: all T-identity genes present", L), length(DBL_T_ID_GENES), length(tid))
  tid_umi <- Matrix::colSums(M[tid, , drop = FALSE]); names(tid_umi) <- colnames(M)
  A[, t_id_umi := tid_umi[cell]]
  A[compartment %in% DBL_LEUK_STATES & !(bmm_prob >= DBL_MIN_BMM_PROB), compartment := "LEUK_other"]   # state unsure
  A[, strict := compartment == "T_cand" & bmm_prob >= DBL_MIN_BMM_PROB & t_id_umi >= DBL_T_ID_MIN_UMI &
               (!DBL_T_REQUIRE_CONCORDANT | anno_source == "concordant") & (is.na(atlas_cl) | atlas_cl != DBL_ATLAS_EXCLUDE_T_CL)]
  A[, compartment := fifelse(compartment == "T_cand", fifelse(strict, "T", "Tamb"), compartment)]
  ts(sprintf("    class A: %d cells; T_cand -> strict T %d, Tamb %d; atlas cluster %d among class A: %d", nrow(A),
             A[compartment == "T", .N], A[strict == FALSE & bmm_fine %in% MAP[compartment == "T_cand", bmm_fine], .N],
             DBL_ATLAS_HOMOTYPIC_CL, A[atlas_cl == DBL_ATLAS_HOMOTYPIC_CL, .N]))
  ## -- [3] normalise, HVGs and PCA on class A; project B/C; kNN transfer  [H 4.4.4] ----
  so <- CreateSeuratObject(counts = M, min.cells = 0, min.features = 0)
  so <- NormalizeData(so, verbose = FALSE)
  soA <- subset(so, cells = A$cell)
  soA <- FindVariableFeatures(soA, nfeatures = DBL_TRANSFER_NHVG, verbose = FALSE)
  hv  <- VariableFeatures(soA)
  soA <- ScaleData(soA, features = hv, verbose = FALSE)
  set.seed(DBL_SEED_TRANSFER)
  soA <- RunPCA(soA, features = hv, npcs = DBL_TRANSFER_NPC, verbose = FALSE)
  load <- Loadings(soA, "pca")[hv, ]                           # genes x PCs
  D   <- LayerData(so, layer = "data")[hv, , drop = FALSE]
  DA  <- D[, A$cell, drop = FALSE]
  mu  <- Matrix::rowMeans(DA); sdv <- sqrt(pmax(Matrix::rowMeans(DA^2) - mu^2, 0) * ncol(DA) / (ncol(DA) - 1))
  proj <- function(cells) { Z <- (as.matrix(D[, cells, drop = FALSE]) - mu) / sdv; Z[Z > 10] <- 10; t(Z) %*% load }   # ScaleData's clip at 10
  pcA <- Embeddings(soA, "pca")[A$cell, ]
  chk(sprintf("%s: re-projected class-A PCs reproduce RunPCA embeddings", L), TRUE,
      isTRUE(all.equal(unname(proj(head(A$cell, 200))), unname(pcA[head(A$cell, 200), ]), tolerance = 1e-4)))
  ## known answer for the transfer: hold out 20 % of class A, transfer from the other 80 %  (labels of singlets only)
  set.seed(DBL_SEED_TRANSFER + 1L)
  ho <- sample(nrow(A), round(0.2 * nrow(A)))
  tr <- knn_vote(pcA[-ho, ], A$compartment[-ho], pcA[ho, , drop = FALSE], DBL_TRANSFER_K)
  tr[, truth := A$compartment[ho]][, assigned := fifelse(vote >= DBL_TRANSFER_MIN_VOTE, label, "Unassigned")]
  TRANS[[L]] <- tr[, .(n = .N, correct = sum(assigned == truth), unassigned = sum(assigned == "Unassigned")), by = truth][
    , `:=`(lane = L, accuracy = correct / n)][order(truth)]
  acc <- tr[, mean(assigned == truth)]
  lin <- function(z) fifelse(z %in% c(DBL_LEUK_STATES, "LEUK_other"), "LEUK", z)   # lineage level: the five states + LEUK_other merged
  acc_lin <- tr[assigned != "Unassigned", mean(lin(assigned) == lin(truth))]
  ts(sprintf("    held-out class-A transfer accuracy %.3f (unassigned %.3f); among assigned, lineage level %.3f, state level %.3f",
             acc, tr[, mean(assigned == "Unassigned")], acc_lin, tr[assigned != "Unassigned", mean(assigned == truth)]))
  TRANS[[L]][, lineage_accuracy_assigned := acc_lin]
  chk(sprintf("%s: held-out class-A transfer accuracy >= 0.80 (known answer)", L), TRUE, acc >= 0.80)
  ## transfer to B and C
  BC <- bc[class != "A", cell]
  tb <- knn_vote(pcA, A$compartment, proj(BC), DBL_TRANSFER_K)
  tb[, cell := BC][, compartment := fifelse(vote >= DBL_TRANSFER_MIN_VOTE, label, "Unassigned")]
  tb[compartment == "T" & tid_umi[cell] < DBL_T_ID_MIN_UMI, compartment := "Tamb"]   # strict T needs T-identity UMI here too
  ts(sprintf("    B/C transfer: %d barcodes, median vote %.2f, %.3f below %.1f (-> Unassigned)", nrow(tb), median(tb$vote),
             mean(tb$vote < DBL_TRANSFER_MIN_VOTE), DBL_TRANSFER_MIN_VOTE))
  ## -- [4] cycling flag: Seurat CellCycleScoring per lane  [P 3.4, D5] ----
  set.seed(DBL_SEED_TRANSFER + 2L)
  so <- CellCycleScoring(so, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
  cc <- data.table(cell = colnames(so), phase = so$Phase, s_score = so$S.Score, g2m_score = so$G2M.Score)
  ## known answer: BMM's cycling labels are enriched for S/G2M among class-A cells
  kk <- merge(A[, .(cell, bmm_fine)], cc, by = "cell")
  cyc_lab <- c("GMP-Cycle", "Cycling Progenitor", "T Proliferating", "NK Proliferating", "Pro-B Cycling", "Pre-pDC Cycling")
  r_in <- kk[bmm_fine %in% cyc_lab, mean(phase != "G1")]; r_out <- kk[!bmm_fine %in% cyc_lab, mean(phase != "G1")]
  ts(sprintf("    S/G2M share: BMM cycling labels %.2f (n %d) vs other class-A cells %.2f", r_in, kk[bmm_fine %in% cyc_lab, .N], r_out))
  chk(sprintf("%s: S/G2M share higher in BMM cycling labels than in other class-A cells (known answer)", L), TRUE, isTRUE(r_in > r_out))
  ## -- [5] assemble and write ----
  lab <- rbind(A[, .(cell, bmm_fine, compartment, source = "projection", vote = NA_real_)],
               tb[, .(cell, bmm_fine = NA_character_, compartment, source = "knn_transfer", vote)])
  lab <- merge(bc[, .(cell, class, donor)], lab, by = "cell")
  lab <- merge(lab, cc, by = "cell")
  lab[, t_id_umi := tid_umi[cell]][, cycling := phase != "G1"]
  chk(sprintf("%s: every analysis barcode has one label", L), ncol(M), nrow(lab))
  chk(sprintf("%s: no missing compartment", L), 0L, lab[is.na(compartment), .N])
  fwrite(lab, out); ts("    [write] ", out)
  COMP[[L]] <- A[, .N, by = .(donor, compartment)][, lane := L][order(donor, compartment)]
  print(dcast(COMP[[L]], compartment ~ donor, value.var = "N", fill = 0L))
  rm(x, M, so, soA, D, DA); gc(FALSE)
}
if (length(COMP))  fwrite(rbindlist(COMP),  file.path(DBL_TAB_DIR, "dbl_class_a_composition_smoke.csv"))
if (length(TRANS)) fwrite(rbindlist(TRANS), file.path(DBL_TAB_DIR, "dbl_transfer_check_smoke.csv"))
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, "dbl_checks_02_smoke.csv"))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
