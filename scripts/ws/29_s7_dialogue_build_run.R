#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds          (S4 kept cells and their compartment)
#          /LARGE1/.../02_seurat_objects/04_annotated/GSE185381/<sample>.rds  (raw counts)
#          results/tables/ws/ws_s4_pseudobulk_meta_GSE185381.csv             (status, hb/granule shares, paediatric)
#          results/tables/ws/decontx/GSE185381/ws_decontx_index.csv.gz        (hierarchy_bin: progenitor fraction)
#          scripts/ws/PREREGISTRATION_S7_dialogue.md (tag prereg-s7-v1), config_ws.R section 12
# OUTPUT : /LARGE1/.../aml_niche_net/ws_dialogue/rA_<run>.rds        DIALOGUE cell.type objects per run (R1, R2)
#          /LARGE1/.../aml_niche_net/ws_dialogue/R_<run>.rds         DIALOGUE.run results (R1, R2, planted)
#          /LARGE1/.../aml_niche_net/ws_dialogue/<run>/              DIALOGUE's own results.dir
#          results/tables/ws/ws_s7_build.csv                          cells and samples per compartment and run; timings
#          results/tables/ws/ws_s7_build_checks.csv
# what it does: PREREGISTRATION_S7_dialogue.md sections 2 and 4.1. Builds one DIALOGUE cell.type per compartment
#   (T, NK, B, Mono) from the S4 kept cells: Seurat LogNormalize, 2,000 variable genes, 30 PCs within the
#   compartment, cellQ = scaled log UMI, metadata with the covariates; at most 500 cells per sample and
#   compartment. Runs DIALOGUE twice (R1 adult AML; R2 adult AML + healthy with status as covariate) and once more
#   on R1's inputs with a planted multicellular program (known answer; planted BEFORE the PCA, prereg section 8). Evaluation of the criteria, the shuffled
#   null and the stability loop are in 30_s7_dialogue_evaluate.R, which reads the objects saved here.
#   usage: Rscript scripts/ws/29_s7_dialogue_build_run.R [--runs=R1,R2,planted] [--npcs=30] [--tag=]
#          tool checks (prereg section 8, item 7): --runs=plantedR2 plants on R2's inputs; --npcs changes the PCs per compartment;
#          --tag suffixes every output (rA_/R_/planted_truth/results dir/build tables) so the pre-registered objects are not touched
#          v2 (prereg section 9.2): --build-only saves the cell.type objects without running DIALOGUE; --run-only reads them back and runs DIALOGUE
# This script deletes nothing and does not modify its inputs.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject); library(DIALOGUE) })
.a <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
RUNS <- strsplit(.get("runs", "R1,R2,planted"), ",")[[1]]; NPCS <- as.integer(.get("npcs", "30")); TAG <- .get("tag", "")
BUILD_ONLY <- "--build-only" %in% .a; RUN_ONLY <- "--run-only" %in% .a; stopifnot(!(BUILD_ONLY && RUN_ONLY))
suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
ws_patch_dialogue()   # DIALOGUE 1.0's MCP cell-type assignment has a row-order bug (config_ws.R section 13; prereg section 8, item 5)
set.seed(SEED)
options(error = function() {                       # DIALOGUE's own errors are terse; print where they come from
  calls <- sys.calls(); n <- length(calls); message("--- TRACEBACK ---")
  for (i in seq(max(1, n - 15), n)) message(i, ": ", substr(paste(deparse(calls[[i]]), collapse = " "), 1, 220)); q(status = 1) })
DATASET <- WS_S4_DEFINE; TABWS <- file.path(TAB_DIR, "ws")
AN   <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated", DATASET)
OUTD <- file.path(LARGE1_DIR, "ws_dialogue"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
CAP <- 500L; ABN <- 15L; K <- 3L; N_PLANT_GENES <- 100L; PLANT_SD <- 0.5
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-64s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- 1. cells, samples, covariates ---------------------------------------------------------------------------
ts("[1] inputs")
PB   <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", DATASET, paste0("ws_pb_", DATASET, ".rds")))
META <- fread(file.path(TABWS, paste0("ws_s4_pseudobulk_meta_", DATASET, ".csv")))
META[, status := fifelse(role == "control", "healthy", fifelse(role == "AML_Dx" & paediatric == FALSE, "AML", NA_character_))]
SM <- unique(META[!is.na(status), .(sample, status, hb_share_nonery, gran_share_lymph)])
IX <- fread(file.path(TABWS, "decontx", DATASET, "ws_decontx_index.csv.gz"), select = c("sample", "hierarchy_bin", "bmm_fine", "high_error"))
IX[, bin := fifelse(bmm_fine %in% c("CLP", "MLP", "MLP-II", "Pre-ProB"), "EarlyLymphoid", hierarchy_bin)]
PF <- IX[!is.na(bin) & bin != "" & !(high_error %in% 1L), .(prog_frac = mean(bin %in% c("HSC_MPP", "LMPP_GMP", "EarlyLymphoid"))), by = sample]
SM <- merge(SM, PF, by = "sample"); chk("covariates for every adult AML / healthy sample", 30L, nrow(SM))
cells <- PB$cells[kept == TRUE & sample %in% SM$sample, .(cell, sample, group4)]
cells <- cells[, if (.N > CAP) .SD[sample.int(.N, CAP)] else .SD, by = .(sample, group4)]      # cap per sample x compartment
cells <- merge(cells, SM, by = "sample")
BUILD <- cells[, .(cells = .N, samples = uniqueN(sample), samples_ge15 = sum(table(sample) >= ABN)), by = .(group4)]
print(BUILD)

## -- 2. one cell.type per compartment, from raw counts ---------------------------------------------------------
ts("[2] building cell.type objects")
samples_needed <- unique(cells$sample); CT <- list(); RAWC <- list()
mats <- list()
if (!RUN_ONLY) for (sm in samples_needed) { o <- readRDS(file.path(AN, paste0(sm, ".rds"))); X <- LayerData(o, assay = "RNA", layer = "counts")
  mats[[sm]] <- X[, cells[sample == sm, cell], drop = FALSE]; rm(o); gc(FALSE) }
if (!RUN_ONLY) { G <- rownames(mats[[1]]); stopifnot(all(vapply(mats, function(m) identical(rownames(m), G), TRUE))) }
PLANTED <- list()
make_ct <- function(g4, keep_samples, name = g4, plant = NULL) {
  cc <- cells[group4 == g4 & sample %in% keep_samples]
  X <- do.call(cbind, lapply(unique(cc$sample), function(sm) mats[[sm]][, cc[sample == sm, cell], drop = FALSE]))[, cc$cell]
  X <- X[Matrix::rowMeans(X > 0) >= 0.01, , drop = FALSE]                       # genes detected in >= 1% of the compartment's cells (dense tpm downstream)
  so <- CreateSeuratObject(X, min.cells = 0, min.features = 0); so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
  if (!is.null(plant)) {   # planted MCP (prereg sections 4.1 and 8): the sample-level factor is added to the log-expression of 100 genes BEFORE
    # ScaleData/RunPCA, so that the PCs DIALOGUE step I works on carry it (the first version added it to tpm only); planted genes become eligible features
    dat <- as.matrix(LayerData(so, assay = "RNA", layer = "data")); expr <- rownames(dat)[rowMeans(dat > 0) > 0.2]   # dense: sparse row sub-assignment took ~5 min per compartment
    pg <- sample(expr, N_PLANT_GENES); sdv <- apply(dat[pg, ], 1, sd)
    dat[pg, ] <- dat[pg, ] + outer(PLANT_SD * sdv, plant$f[cc$sample])
    LayerData(so, assay = "RNA", layer = "data") <- as(dat, "dgCMatrix"); VariableFeatures(so) <- union(VariableFeatures(so), pg); rm(dat)
    PLANTED[[name]] <<- list(genes = pg, sdv = sdv)
  }
  so <- ScaleData(so, verbose = FALSE); so <- RunPCA(so, npcs = NPCS, verbose = FALSE)
  if (!is.null(plant)) PLANTED[[name]]$dir <<- colSums(Loadings(so, "pca")[pg, ] * sdv)   # the planted shift's direction in PC space (planted genes' loadings weighted by their shift)
  tpm <- as.matrix(LayerData(so, assay = "RNA", layer = "data")); Xp <- Embeddings(so, "pca")
  md <- data.frame(samples = cc$sample, status = cc$status, hb_share_nonery = cc$hb_share_nonery, gran_share_lymph = cc$gran_share_lymph,
                   prog_frac = cc$prog_frac, nFeatures = so$nFeature_RNA, row.names = cc$cell)
  cellQ <- as.vector(scale(log(so$nCount_RNA)))
  r <- make.cell.type(name = name, tpm = tpm, samples = cc$sample, X = Xp, metadata = md, cellQ = cellQ)
  rm(so, X); gc(FALSE); r
}
aml <- SM[status == "AML", sample]; all30 <- SM$sample
run_one <- function(run, keep_samples, covar, plant = FALSE) {
  pl <- NULL; fA <- file.path(OUTD, paste0("rA_", run, TAG, ".rds"))
  if (RUN_ONLY) { stopifnot(!plant, file.exists(fA)); rA <- readRDS(fA); ts("    [", run, TAG, "] cell.type objects read back (--run-only)") } else {
  if (plant) { set.seed(SEED + 1L); f <- setNames(rnorm(length(keep_samples)), keep_samples); pl <- list(f = f); PLANTED <<- list() }   # one factor per sample, same sign in every compartment; own seed so the draw does not depend on the runs before it (added 2026-10-20; the 2026-10-19 planted objects were drawn without it)
  rA <- lapply(setNames(WS_S4_GROUPS, WS_S4_GROUPS), make_ct, keep_samples = keep_samples, plant = pl)
  if (plant) {
    saveRDS(list(factor = f, genes = lapply(PLANTED, `[[`, "genes"), dir = lapply(PLANTED, `[[`, "dir")), file.path(OUTD, paste0("planted_truth", TAG, ".rds")))
    chk("planted: 100 genes recorded in each of the 4 compartments", TRUE, length(PLANTED) == 4 && all(vapply(PLANTED, function(p) length(p$genes), 1L) == N_PLANT_GENES))
    # is the planted factor in the PCs DIALOGUE step I works on? (a) sample means of the PC coordinate along the planted direction vs f;
    # (b) the best single sample-mean PC vs f; each against the P95 of 500 random re-assignments of f to samples (chance for (b): median ~0.55, P95 ~0.7;
    # the first, tpm-only version scored 0.38-0.61 on (b), i.e. chance)
    vis <- rbindlist(lapply(names(rA), function(ct) { r <- rA[[ct]]; A <- rowsum(r@X, r@samples) / as.vector(table(r@samples)[sort(unique(r@samples))]); fs <- f[rownames(A)]
      sc <- as.vector(A %*% PLANTED[[ct]]$dir); nul <- replicate(500, { fp <- sample(fs); c(abs(cor(sc, fp)), max(abs(cor(A, fp)))) })
      data.table(compartment = ct, n_samples = nrow(A), dir_cor = abs(cor(sc, fs)), dir_p95 = unname(quantile(nul[1, ], .95)), pc_max_cor = max(abs(cor(A, fs))), pc_p95 = unname(quantile(nul[2, ], .95))) }))
    print(vis[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
    chk("planted: factor visible along the planted direction in PC space, every compartment (cor above P95 of 500 random f)", TRUE, all(vis$dir_cor > vis$dir_p95))
    chk("planted: factor visible in the best single sample-mean PC, every compartment (above P95 of 500 random f)", TRUE, all(vis$pc_max_cor > vis$pc_p95))
  }
  saveRDS(rA, fA) }
  if (BUILD_ONLY) { ts("    [", run, TAG, "] built and saved; DIALOGUE not run (--build-only)")
    return(data.table(run = paste0(run, TAG), npcs = NPCS, minutes = NA_real_, n_samples = length(keep_samples), message = "build only", run_error = "")) }
  rd <- file.path(OUTD, paste0(run, TAG)); dir.create(rd, showWarnings = FALSE)
  t0 <- Sys.time(); err <- ""
  # DIALOGUE.run ignores param$plot.flag and plots by default (its plots failed on these inputs after the results were
  # saved), so plot.flag = FALSE is passed to DIALOGUE.run itself; the fallback below reads results.dir if it still errors
  R <- tryCatch(DIALOGUE.run(rA = rA, main = run, plot.flag = FALSE, param = DLG.get.param(k = K, results.dir = paste0(rd, "/"), conf = "cellQ", covar = covar,
                                                                        abn.c = ABN, seed1 = SEED, plot.flag = FALSE)),
                error = function(e) { err <<- conditionMessage(e); f <- file.path(rd, paste0("DLG.full.output_", run, ".rds"))
                  if (!file.exists(f)) stop("DIALOGUE failed before writing results: ", err); ts("    DIALOGUE.run error after its results were written (", err, "); reading ", basename(f)); readRDS(f) })
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins")); saveRDS(R, file.path(OUTD, paste0("R_", run, TAG, ".rds")))
  ts(sprintf("    %s done in %.1f min | message: %s | class: %s | names: %s", run, el, if (!is.null(R$message)) R$message else "-", class(R)[1], paste(head(names(R), 20), collapse = ",")))
  data.table(run = paste0(run, TAG), npcs = NPCS, minutes = round(el, 1), n_samples = length(keep_samples), message = if (!is.null(R$message)) R$message else "", run_error = err)
}
TIM <- list()
if ("R1" %in% RUNS)      TIM[["R1"]] <- run_one("R1", aml, covar = WS_S7_COVAR$R1)          # covariates frozen in config_ws.R section 14
if ("R2" %in% RUNS)      TIM[["R2"]] <- run_one("R2", all30, covar = WS_S7_COVAR$R2)
if ("planted" %in% RUNS) TIM[["planted"]] <- run_one("planted", aml, covar = WS_S7_COVAR$R1, plant = TRUE)
if ("plantedR2" %in% RUNS) TIM[["plantedR2"]] <- run_one("plantedR2", all30, covar = WS_S7_COVAR$R2, plant = TRUE)   # tool check: planted on R2's 23-sample input
BT <- cbind(rbindlist(TIM), BUILD[1, .(note = paste(BUILD[, sprintf("%s %d cells/%d samples", group4, cells, samples)], collapse = "; "))])
bf <- file.path(TABWS, "ws_s7_build.csv"); if (file.exists(bf)) BT <- rbind(fread(bf)[!run %in% paste0(RUNS, TAG)], BT, fill = TRUE)   # a partial rerun keeps the other runs' rows
fwrite(BT, bf)
C <- rbindlist(CHK); fwrite(C, file.path(TABWS, paste0("ws_s7_build_checks", TAG, ".csv")))
ts(sprintf("[done] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
