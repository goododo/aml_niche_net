#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds          (S4: per-cell kept flag and group4)
#          results/tables/ws/ws_s4_pseudobulk_meta_GSE185381.csv           (eligibility, status, covariates)
#          /LARGE1/.../02_seurat_objects/04_annotated/GSE185381/<sample>.rds  (raw counts of the eligible samples)
#          results/tables/ws/ws_s5_imprint_genes.csv                        (T imprint, for the reported overlap only)
#          PREREGISTRATION_S5_imprint.md section 12 (settings fixed before this run)
# OUTPUT : /LARGE1/.../aml_niche_net/ws_dialogue/<run>/            DIALOGUE's own results directory and the full R object (rds)
#          results/tables/ws/ws_s7_mcp_genes.csv      run, MCP, compartment, direction, gene
#          results/tables/ws/ws_s7_pref.csv           run, MCP, compartment pair, cross-compartment R and p
#          results/tables/ws/ws_s7_sample_scores.csv  run, MCP, compartment, sample, mean cell score, status, age group
#          results/tables/ws/ws_s7_stability.csv      Jaccard of top-50 genes, 10 runs on 80% of samples (lite version)
#          results/tables/ws/ws_s7_permutation.csv    cross-compartment R under 10 sample-label permutations (lite version)
#          results/tables/ws/ws_s7_reading.csv        the section-12 reading rule per MCP
#          results/tables/ws/ws_s7_checks.csv
# what it does: S7 of v2.0, as fixed in PREREGISTRATION_S5_imprint.md section 12. Builds one DIALOGUE cell.type per
#   compartment from the S4 kept cells of GSE185381 adults and healthy donors (log1p CP10k, 2,000 variable genes,
#   15 PCs, cellQ = log10 UMIs, covariates cellQ + the two raw ambient shares), runs DIALOGUE (k = 3) on the samples
#   shared by all compartments of the run, and evaluates the three-part reading rule: every-pair positive R with
#   p < 0.05, above the permutation 95th percentile, and stable under patient resampling. Two runs: primary
#   (T, NK, B, Mono) and secondary (T, Mono). Association with status and age group is reported, never tested.
#   usage: Rscript scripts/ws/29_s7_dialogue.R [--run=primary|secondary] [--smoke]   (--smoke: 3 resamples / 3 permutations)
# This script deletes nothing and does not modify its inputs.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject); library(DIALOGUE) })
.a <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
RUNS <- strsplit(.get("run", "primary,secondary"), ",")[[1]]; SMOKE <- "--smoke" %in% .a
suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
set.seed(SEED)
DATASET <- WS_S4_DEFINE; TABWS <- file.path(TAB_DIR, "ws"); AN <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated", DATASET)
OUTD <- file.path(LARGE1_DIR, "ws_dialogue"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
PB   <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", DATASET, paste0("ws_pb_", DATASET, ".rds")))
META <- fread(file.path(TABWS, paste0("ws_s4_pseudobulk_meta_", DATASET, ".csv")))
IMP  <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))
RUNDEF <- list(primary = c("T", "NK", "B", "Mono"), secondary = c("T", "Mono"))
N_HVG <- 2000L; N_PC <- 15L; K <- 3L; N_RES <- if (SMOKE) 3L else 10L; N_PERM <- if (SMOKE) 3L else 10L; TOPN <- 50L
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-64s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- 1. samples, cells, per-compartment inputs (built once, reused by both runs) --------------------------------
META[, status := fifelse(role == "control", "healthy", fifelse(role == "AML_Dx" & paediatric == FALSE, "AML", NA_character_))]
SM <- unique(META[!is.na(status), .(sample, status, age_group, hb = hb_share_nonery, gran = gran_share_lymph)])
ELIG <- dcast(META[!is.na(status), .(sample, group4, eligible)], sample ~ group4, value.var = "eligible")
cells <- PB$cells[kept == TRUE & sample %in% SM$sample, .(cell, sample, group4, fine)]
ts(sprintf("[1] %d samples (%d AML, %d healthy), %d kept cells", nrow(SM), SM[status == "AML", .N], SM[status == "healthy", .N], nrow(cells)))
CT <- list(); QC <- list()
for (sm in SM$sample) {                                   # one pass over the objects: log-normalised HVG-ready data per compartment
  o <- readRDS(file.path(AN, paste0(sm, ".rds"))); X <- LayerData(o, assay = "RNA", layer = "counts"); cc <- cells[sample == sm, cell]
  X <- X[, intersect(colnames(X), cc), drop = FALSE]; QC[[sm]] <- data.table(cell = colnames(X), umi = Matrix::colSums(X))
  for (g4 in WS_S4_GROUPS) { j <- cells[sample == sm & group4 == g4 & cell %in% colnames(X), cell]; if (length(j)) CT[[g4]][[sm]] <- X[, j, drop = FALSE] }
  rm(o, X); gc(FALSE)
}
QC <- rbindlist(QC)
build_ct <- function(g4, samples) {                        # a DIALOGUE cell.type for one compartment on the given samples
  M <- do.call(cbind, CT[[g4]][samples]); so <- CreateSeuratObject(M, min.cells = 0, min.features = 0)
  so <- NormalizeData(so, verbose = FALSE); so <- FindVariableFeatures(so, nfeatures = N_HVG, verbose = FALSE)
  so <- ScaleData(so, features = VariableFeatures(so), verbose = FALSE); so <- RunPCA(so, features = VariableFeatures(so), npcs = N_PC, verbose = FALSE)
  smp <- cells[match(colnames(so), cell), sample]; q <- log10(QC[match(colnames(so), cell), umi])
  md <- SM[match(smp, sample), .(status, age_group, hb, gran)]; md <- as.data.frame(md); rownames(md) <- colnames(so)
  tpm <- as.matrix(LayerData(so, assay = "RNA", layer = "data")[VariableFeatures(so), ])
  make.cell.type(name = g4, tpm = tpm, samples = smp, X = Embeddings(so, "pca"), metadata = md, cellQ = q)
}
shared_samples <- function(comps) { e <- ELIG[, c("sample", comps), with = FALSE]; e[rowSums(e[, -1] == TRUE) == length(comps), sample] }

## -- 2. one DIALOGUE run, with extraction -----------------------------------------------------------------------
sig_table <- function(R, run) {                            # sig: list named "MCPk.up.<ct>" / "MCPk.down.<ct>" (DIALOGUE 1.0)
  s <- R$sig; if (is.null(s) || !length(s)) return(NULL)
  rbindlist(lapply(names(s), function(n) { p <- regmatches(n, regexec("^(MCP[0-9]+)\\.(up|down)\\.(.+)$", n))[[1]]
    if (length(p) < 4) return(NULL); data.table(run = run, mcp = p[2], direction = p[3], compartment = p[4], gene = s[[n]], rank = seq_along(s[[n]])) }))
}
pref_table <- function(R, run) {
  p <- R$pref; if (is.null(p)) return(NULL)
  out <- tryCatch({ if (is.data.frame(p) || is.matrix(p)) as.data.table(as.data.frame(p), keep.rownames = "key") else rbindlist(lapply(names(p), function(n) data.table(key = n, as.data.table(as.data.frame(p[[n]]), keep.rownames = "pair"))), fill = TRUE) }, error = function(e) NULL)
  if (!is.null(out)) out[, run := run]; out
}
run_dialogue <- function(name, comps, samples, full = TRUE, seed = SEED, relabel = NULL) {
  rA <- lapply(comps, function(g4) { ct <- build_ct(g4, samples); if (!is.null(relabel)) ct@samples <- relabel[[g4]][ct@samples]; ct })
  names(rA) <- comps
  dir.create(file.path(OUTD, name), showWarnings = FALSE, recursive = TRUE)
  param <- DLG.get.param(k = K, results.dir = file.path(OUTD, name, ""), seed1 = seed, conf = "cellQ", covar = c("cellQ", "hb", "gran"),
                         n.genes = 200, plot.flag = FALSE, full.version = full)
  R <- DIALOGUE.run(rA = rA, main = name, param = param, plot.flag = FALSE)
  list(R = R, rA = rA)
}

## -- 3. primary and secondary runs ---------------------------------------------------------------------------------
GENES <- list(); PREF <- list(); SCORES <- list(); STAB <- list(); PERM <- list(); READ <- list()
for (run in RUNS) {
  comps <- RUNDEF[[run]]; smp <- shared_samples(comps)
  ts(sprintf("[%s] compartments %s | shared eligible samples %d (%d AML, %d healthy)", run, paste(comps, collapse = "+"), length(smp), SM[sample %in% smp & status == "AML", .N], SM[sample %in% smp & status == "healthy", .N]))
  chk(paste(run, "has >= 10 shared samples"), TRUE, length(smp) >= 10)
  if (length(smp) < 10) next
  res <- run_dialogue(run, comps, smp, full = TRUE); R <- res$R
  saveRDS(R, file.path(OUTD, run, paste0("ws_dialogue_", run, ".rds")))
  ts("    result components: ", paste(names(R), collapse = ", "))
  G <- sig_table(R, run); P <- pref_table(R, run)
  if (!is.null(G)) GENES[[run]] <- G; if (!is.null(P)) PREF[[run]] <- P
  ts(sprintf("    signatures: %s", if (is.null(G)) "none parsed" else paste(G[, .N, by = .(mcp, compartment, direction)][, sprintf("%s %s %s %d", mcp, compartment, direction, N)], collapse = "; ")))
  # sample-level scores: mean of the per-cell MCP scores per sample and compartment (R$scores: list per cell type, cells x MCPs)
  if (!is.null(R$scores)) for (g4 in names(R$scores)) { sc <- R$scores[[g4]]; if (is.null(dim(sc))) next
    d <- as.data.table(sc, keep.rownames = "cell"); d <- melt(d, id.vars = "cell", variable.name = "mcp", value.name = "score")
    d[, sample := cells[match(cell, cell), sample]]; d <- d[, .(score = mean(score), cells = .N), by = .(mcp, sample)]
    d <- merge(d, SM[, .(sample, status, age_group)], by = "sample"); d[, `:=`(run = run, compartment = g4)]; SCORES[[paste(run, g4)]] <- d }
  # stability: lite runs on 80% of the samples, Jaccard of the top-50 genes per MCP x compartment x direction against the full run
  if (!is.null(G)) {
    top <- G[rank <= TOPN, .(genes = list(gene)), by = .(mcp, compartment, direction)]
    for (i in seq_len(N_RES)) { ss <- sample(smp, ceiling(0.8 * length(smp)))
      Gi <- tryCatch(sig_table(run_dialogue(paste0(run, "_res", i), comps, ss, full = FALSE, seed = SEED + i)$R, run), error = function(e) NULL)
      if (is.null(Gi)) next
      ti <- Gi[rank <= TOPN, .(genes = list(gene)), by = .(mcp, compartment, direction)]
      j <- merge(top, ti, by = c("mcp", "compartment", "direction"), suffixes = c("", "_i"))
      STAB[[paste(run, i)]] <- j[, .(run = run, resample = i, mcp, compartment, direction, jaccard = mapply(function(a, b) length(intersect(a, b)) / length(union(a, b)), genes, genes_i))] }
    # permutation: sample labels permuted independently within each compartment (a bijection per compartment)
    for (i in seq_len(N_PERM)) { rl <- lapply(setNames(comps, comps), function(g4) setNames(sample(smp), smp))
      Pi <- tryCatch(pref_table(run_dialogue(paste0(run, "_perm", i), comps, smp, full = FALSE, seed = SEED + 100 + i, relabel = rl)$R, run), error = function(e) NULL)
      if (!is.null(Pi)) { Pi[, perm := i]; PERM[[paste(run, i)]] <- Pi } }
  }
}
GENES <- rbindlist(GENES); PREF <- rbindlist(PREF, fill = TRUE); SCORES <- rbindlist(SCORES); STAB <- rbindlist(STAB); PERM <- rbindlist(PERM, fill = TRUE)
for (nm in c("GENES", "PREF", "SCORES", "STAB", "PERM")) { x <- get(nm); if (nrow(x)) fwrite(x, file.path(TABWS, paste0("ws_s7_", c(GENES = "mcp_genes", PREF = "pref", SCORES = "sample_scores", STAB = "stability", PERM = "permutation")[[nm]], ".csv"))) }
## -- 4. reported overlap with the S5 T imprint, and the reading rule where the tables allow ------------------------
if (nrow(GENES)) {
  timp <- IMP[compartment == "T", gene]
  ov <- GENES[compartment == "T", .(n_genes = .N, in_T_imprint = sum(gene %in% timp)), by = .(run, mcp, direction)]
  ts("    T-compartment MCP genes inside the S5 T imprint: ", paste(ov[, sprintf("%s %s %s %d/%d", run, mcp, direction, in_T_imprint, n_genes)], collapse = "; "))
  fwrite(ov, file.path(TABWS, "ws_s7_overlap_T_imprint.csv"))
}
if (nrow(STAB)) { st <- STAB[, .(median_jaccard = median(jaccard), n = .N), by = .(run, mcp)]; ts("    stability (median Jaccard of top-50): ", paste(st[, sprintf("%s %s %.2f", run, mcp, median_jaccard)], collapse = "; ")); READ$stability <- st }
if (nrow(PREF)) ts("    pref columns: ", paste(names(PREF), collapse = ", "))
C <- rbindlist(CHK); fwrite(C, file.path(TABWS, "ws_s7_checks.csv"))
ts(sprintf("[done] %d checks, %d PASS, %d FAIL -- the reading rule (section 12) is evaluated from ws_s7_pref/permutation/stability once their layout is confirmed", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
