#!/usr/bin/env Rscript
# INPUT  : scripts/ws/config_ws.R  (section 10, frozen and committed in 4c54f6f before this script existed)
#          results/tables/01_preprocess/00_curated_manifest.csv               (eligibility, as T3c)
#          /LARGE1/.../02_seurat_objects/04_annotated/<dataset>/<sample>.rds   (RNA counts, hierarchy_bin, high_error)
#          results/tables/ws/ws_doublet_scores_percell.csv.gz                  (scDblFinder score, union flag)
#          /LARGE1/gr10634/gaozy/reference/cytosig/signature.centroid         (sha256-checked by ws_assert_b0_frozen)
#          progeny::getModel (bundled, offline); CellChatDB.human (fake-ligand exclusion list)
# OUTPUT : results/tables/ws/ws_b0_inputs<SFX>.csv    one row per (sample, arm): cell counts, covariates, ligands, footprints
#          results/tables/ws/ws_b0_results<SFX>.csv   one row per (analysis, ligand, arm): n, t, p_FL, fake P95, pass
#          results/tables/ws/ws_b0_fakes<SFX>.csv.gz  fake-ligand t per (analysis, ligand, arm)
#          results/tables/ws/ws_b0_tier<SFX>.csv      the verdict (gate dataset only)
#          results/tables/ws/ws_b0_checks<SFX>.csv
#          scripts/ws/GATES.md                        "B0 result" block between markers (gate run, not smoke)
#          /LARGE1/.../aml_niche_net/ws_b0_cache/<dataset>/<sample>.rds   per-sample pseudobulk sums (kept)
# what it does: GATE_B0, handoff v2.0 §8.3. Per AML patient it sums counts over the leukaemic progenitor
#   bins (HSC_MPP + LMPP_GMP) and over T_NK. Outcome = CytoSig footprint of the progenitor pseudobulk
#   (decoupleR ulm on gene-centred log2 CPM, ligand and receptor genes removed). Exposure = log1p of the
#   T_NK counts per 10k for the ligand. Covariates = prog_frac, ambient_z, dbl_mean, paediatric. The ligand
#   t is tested by a one-sided Freedman-Lane permutation and compared with the same t for 200 fake ligands
#   from the same sender. Every number is read from config_ws.R section 10.
#   --smoke: 4 samples; builds caches and footprints and runs the synthetic known-answer checks, but fits
#   nothing on real outcomes. --dataset=GSE289435|GSE227903: direction report, never a gate.
# This script deletes nothing and does not modify the 04_annotated objects.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
.a <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
SMOKE   <- "--smoke" %in% .a
DATASET <- .get("dataset", "GSE185381")

suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
ws_assert_frozen(); ws_assert_b0_frozen()
GATE <- DATASET == WS_B0_GATE_DATASET
if (!GATE && !DATASET %in% WS_B0_REPORT_DATASETS) stop("dataset not named in config_ws.R section 10: ", DATASET)
SFX    <- paste0(if (GATE) "" else paste0("_", DATASET), if (SMOKE) "_smoke" else "")
AN     <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated", DATASET)
CACHE  <- file.path(LARGE1_DIR, "ws_b0_cache", DATASET); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
TABWS  <- file.path(TAB_DIR, "ws")
GATES  <- here::here("scripts", "ws", "GATES.md")
N_PERM <- if (SMOKE) 500L else WS_B0_N_PERM
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-60s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}
ts(sprintf("GATE_B0 | dataset=%s (%s) | smoke=%s | perms=%d | fakes=%d", DATASET,
           if (GATE) "GATE" else "direction report", SMOKE, N_PERM, WS_FAKE_N))

## -- 1. inputs ----------------------------------------------------------------------------------
ts("[1] inputs")
man <- fread(file.path(TAB_DIR, "01_preprocess", "00_curated_manifest.csv"))
EL  <- man[dataset == DATASET & disease == "AML" & timepoint == "Diagnosis" & tissue_r == "BM" &
           sorting_r %in% c("viability_only", "none"), .(sample, blast_pct_clinical)]
EL[, paediatric := as.numeric(DATASET == "GSE185381" & !grepl("^AML[0-9]{4}$", sample))]   # rule: config section 10
if (GATE) { chk("eligible samples (T3c rule)", WS_W11_ELIGIBLE_SAMPLES, nrow(EL))
            chk("paediatric samples (v2.0 5.2)", 22L, as.integer(sum(EL$paediatric))) }
if (SMOKE) EL <- EL[seq_len(min(4L, .N))]
hdr <- strsplit(readLines(WS_B0_CYTOSIG_FILE, n = 1L), "\t")[[1]]          # the header line has no gene column
CS  <- fread(WS_B0_CYTOSIG_FILE, skip = 1L, header = FALSE, sep = "\t"); setnames(CS, c("gene", hdr))
chk("CytoSig genes x cytokines", "4881x43", sprintf("%dx%d", nrow(CS), length(hdr)))
chk("CytoSig carries the three signatures", 3L, sum(WS_B0_PAIRS$cytosig %in% hdr))
PG <- melt(as.data.table(progeny::getModel("Human", top = WS_B0_PROGENY_TOP), keep.rownames = "gene"),
           id.vars = "gene", variable.name = "source", value.name = "mor")[mor != 0]
PG[, source := as.character(source)]
suppressPackageStartupMessages(library(CellChat)); data("CellChatDB.human", package = "CellChat")
cc   <- as.data.table(CellChatDB.human$interaction); cplx <- CellChatDB.human$complex
expand <- function(r) if (r %in% rownames(cplx)) unlist(strsplit(gsub(" ", "", paste(cplx[r, ], collapse = ",")), ",")) else r
LR_GENES  <- setdiff(unique(c(unlist(lapply(cc$ligand, expand)), unlist(lapply(cc$receptor, expand)))), "")
RECEPTORS <- setNames(strsplit(WS_B0_PAIRS$receptors, ";"), WS_B0_PAIRS$ligand)
DBL <- fread(file.path(TABWS, WS_DBL_SCORES_TABLE))[dataset == DATASET,
        .(cell, sample, sc_score, dbl_union = (sc_class %in% c(TRUE, "TRUE")) | (df_class %in% c(TRUE, "TRUE")))]
chk("doublet scores cover the eligible samples", TRUE, all(EL$sample %chin% unique(DBL$sample)))
DBL[, sample := NULL]; setkey(DBL, cell)

## -- 2. per-sample pseudobulk sums, cached -------------------------------------------------------
ts("[2] per-sample pseudobulk (cache: ", CACHE, ")")
GROUPS <- list(recv = WS_B0_RECEIVER_BINS, tnk = WS_B0_SENDER_BIN, mono = WS_B0_CLONE_INTERNAL[["sender"]], b = "B_Plasma")
pb_one <- function(sm) {
  cf <- file.path(CACHE, paste0(sm, ".rds")); if (file.exists(cf)) return(readRDS(cf))
  f <- file.path(AN, paste0(sm, ".rds")); if (!file.exists(f)) return(NULL)
  o <- readRDS(f); md <- o@meta.data
  b <- as.character(md$hierarchy_bin); b[b %in% ""] <- NA
  if (WS_DROP_HIGH_ERROR) b[!is.na(md$high_error) & md$high_error == 1] <- NA
  cts <- LayerData(o, assay = "RNA", layer = "counts"); tot <- Matrix::colSums(cts)
  d <- DBL[J(colnames(cts))]; dsc <- d$sc_score; dsc[is.na(dsc)] <- 0; dun <- d$dbl_union %in% TRUE
  out <- list(sample = sm, genes = rownames(cts), n_cells = ncol(cts), dbl_mean = mean(dsc), n_dbl_union = sum(dun))
  for (arm in c("main", "sens")) {
    bb <- if (arm == "main") b else replace(b, dun, NA)       # sens arm: union-flagged doublets removed
    A <- list(n_binned = sum(!is.na(bb)))
    for (g in names(GROUPS)) {
      j <- which(bb %in% GROUPS[[g]]); A[[paste0("n_", g)]] <- length(j)
      if (!length(j)) next
      A[[paste0("sum_", g)]] <- as.numeric(Matrix::rowSums(cts[, j, drop = FALSE]))
      if (g %in% c("tnk", "mono")) A[[paste0("det_", g)]] <- as.numeric(Matrix::rowMeans(cts[, j, drop = FALSE] > 0))
    }
    out[[arm]] <- A
  }
  hb <- intersect(c("HBB", "HBA1", "HBA2"), rownames(cts)); gr <- intersect(c("MPO", "ELANE", "AZU1", "PRTN3"), rownames(cts))
  j1 <- which(!is.na(b) & b != "Erythroid"); j2 <- which(b %in% c("T_NK", "B_Plasma"))
  out$amb_hb   <- if (length(j1)) sum(cts[hb, j1]) / sum(tot[j1]) else NA_real_
  out$amb_gran <- if (length(j2)) sum(cts[gr, j2]) / sum(tot[j2]) else NA_real_
  saveRDS(out, cf)                                          # cache FIRST, before anything downstream can fail
  ts(sprintf("    %-10s recv=%d tnk=%d mono=%d b=%d", sm, out$main$n_recv, out$main$n_tnk, out$main$n_mono, out$main$n_b))
  rm(o, cts); gc(FALSE); out
}
PB <- lapply(EL$sample, pb_one); names(PB) <- EL$sample; PB <- Filter(Negate(is.null), PB)
chk("samples cached", nrow(EL), length(PB))
G <- Reduce(union, lapply(PB, `[[`, "genes"))
chk("one gene axis across samples", TRUE, all(vapply(PB, function(x) identical(x$genes, PB[[1]]$genes), TRUE)))
vec <- function(x, arm, field) {                            # a cached vector placed on the union gene axis G
  v <- x[[arm]][[field]]; if (is.null(v)) return(rep(NA_real_, length(G)))
  z <- numeric(length(G)); z[match(x$genes, G)] <- v; z
}

## -- 3. inputs table, footprints, exposures ------------------------------------------------------
ts("[3] inputs, footprints, exposures")
IN <- rbindlist(lapply(PB, function(x) rbindlist(lapply(c("main", "sens"), function(arm) data.table(
  sample = x$sample, arm = arm, n_binned = x[[arm]]$n_binned, n_recv = x[[arm]]$n_recv, n_tnk = x[[arm]]$n_tnk,
  n_mono = x[[arm]]$n_mono, n_b = x[[arm]]$n_b, dbl_mean = x$dbl_mean, amb_hb = x$amb_hb, amb_gran = x$amb_gran)))))
IN <- merge(IN, EL, by = "sample")
IN[, prog_frac := n_recv / n_binned]
elig_main <- IN[arm == "main" & n_recv >= WS_B0_MIN_RECEIVER & n_tnk >= WS_B0_MIN_SENDER, sample]
IN[, eligible := sample %in% elig_main & n_recv >= WS_B0_MIN_RECEIVER & n_tnk >= WS_B0_MIN_SENDER]
zs <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
IN[eligible == TRUE, ambient_z := rowMeans(cbind(zs(amb_hb), zs(amb_gran)), na.rm = TRUE), by = arm]
ts(sprintf("    eligible: main %d, sens %d, of %d cached", IN[arm == "main", sum(eligible)], IN[arm == "sens", sum(eligible)], length(PB)))

pb_mat <- function(arm, field, samples) {                   # genes x samples: filtered, log2(CPM+1), gene-centred
  M <- vapply(samples, function(s) vec(PB[[s]], arm, field), numeric(length(G))); rownames(M) <- G
  cpm <- t(t(M) / colSums(M)) * 1e6
  L <- log2(cpm[rowMeans(cpm > 1) >= 0.5, , drop = FALSE] + 1)   # WS_B0_GENE_FILTER
  L - rowMeans(L)
}
fp_cytosig <- function(M, cytokine, drop) {
  net <- data.frame(source = cytokine, target = CS$gene, mor = CS[[cytokine]])
  net <- net[net$mor != 0 & !net$target %in% drop & net$target %in% rownames(M), ]
  r <- decoupleR::run_ulm(mat = M, network = net, minsize = 5L)
  setNames(r$score, r$condition)[colnames(M)]
}
fp_progeny <- function(M, pathway, drop) {
  net <- as.data.frame(PG[!gene %in% drop & gene %in% rownames(M), .(source, target = gene, mor)])
  r <- decoupleR::run_mlm(mat = M, network = net, minsize = 5L); r <- r[r$source == pathway, ]
  setNames(r$score, r$condition)[colnames(M)]
}
lig <- function(arm, field, gene, samples, frac = FALSE) vapply(samples, function(s) {
  v <- vec(PB[[s]], arm, field); cp <- v[match(gene, G)] / sum(v) * 1e4
  if (frac) cp <- cp * PB[[s]][[arm]][[sub("sum_", "n_", field)]] / PB[[s]][[arm]]$n_binned   # x sender fraction
  log1p(cp) }, 1)
fake_pool <- function(arm, field, samples, ligand, cytokine) {
  S <- vapply(samples, function(s) { v <- vec(PB[[s]], arm, field); v / sum(v) * 1e4 }, numeric(length(G)))
  D <- vapply(samples, function(s) vec(PB[[s]], arm, sub("sum_", "det_", field)), numeric(length(G)))
  cand <- data.table(gene = G, det = rowMeans(D), IQR = apply(S, 1, IQR), nz = rowSums(S > 0)); ref <- cand[gene == ligand]
  excl <- unique(c(LR_GENES, CS$gene[order(-CS[[cytokine]])][seq_len(WS_B0_FAKE_TOPSIG)], ligand, unlist(RECEPTORS)))
  pool <- function(tol) cand[!gene %in% excl & nz >= 3 & abs(det - ref$det) <= tol * ref$det & abs(IQR - ref$IQR) <= tol * ref$IQR, gene]
  fk <- pool(WS_FAKE_TOL); tol <- WS_FAKE_TOL
  if (length(fk) < WS_FAKE_MIN_POOL) { fk <- pool(WS_FAKE_TOL_RELAXED); tol <- WS_FAKE_TOL_RELAXED }
  set.seed(WS_FAKE_SEED); if (length(fk) > WS_FAKE_N) fk <- sample(fk, WS_FAKE_N)
  A <- log1p(S[match(fk, G), , drop = FALSE]); colnames(A) <- samples
  list(genes = fk, tol = tol, A = A)
}

## -- model: ligand t with a Freedman-Lane permutation ---------------------------------------------
perm_idx <- function(n, B) {                                # fixed permutations; the global RNG state is restored
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(SEED); P <- vapply(seq_len(B), function(i) sample.int(n), integer(n))
  if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv); P
}
tstat <- function(Yres, ra, df) {                           # t of the exposure, by Frisch-Waugh, one column per outcome
  sa <- sum(ra^2); beta <- colSums(ra * Yres) / sa
  beta / sqrt(colSums((Yres - outer(ra, beta))^2) / df / sa)
}
fl <- function(y, a, X, B, Afake = NULL) {
  ok <- is.finite(y) & is.finite(a) & stats::complete.cases(X)
  y <- y[ok]; a <- a[ok]; X <- X[ok, , drop = FALSE]
  v <- apply(X, 2, var); X <- X[, c(1L, which(v[-1] > 0) + 1L), drop = FALSE]   # keep the intercept, drop constants
  q <- qr(X); n <- length(y); df <- n - q$rank - 1L
  er <- qr.resid(q, y); ra <- qr.resid(q, a)               # reduced-model residuals and the residualised exposure
  t_obs <- tstat(matrix(er), ra, df)
  P <- perm_idx(n, B); t_perm <- tstat(qr.resid(q, matrix(er[P], n)), ra, df)   # y* = fitted + permuted residuals
  t_fake <- if (is.null(Afake)) numeric(0) else apply(Afake[, ok, drop = FALSE], 1, function(af) tstat(matrix(er), qr.resid(q, af), df))
  list(n = n, df = df, t = t_obs, p = (1 + sum(t_perm >= t_obs)) / (B + 1), t_fake = t_fake, covs = colnames(X)[-1])
}
Xmat <- function(d) cbind(`(Intercept)` = 1, as.matrix(d[, WS_B0_COVARIATES, with = FALSE]))

## -- known-answer checks on synthetic data (every mode) -------------------------------------------
ts("[k] known-answer checks")
local({
  set.seed(SEED); n <- 32L; Xs <- cbind(1, matrix(rnorm(n * 4), n)); as_ <- rnorm(n) + 0.5 * Xs[, 2]
  ys <- as.numeric(0.8 * as_ + Xs[, 2:5] %*% c(0.5, -0.3, 0.2, 0) + rnorm(n))
  f1 <- fl(ys, as_, Xs, 2000L)
  chk("FL ligand t equals lm t", TRUE, isTRUE(all.equal(f1$t, unname(summary(lm(ys ~ Xs[, -1] + as_))$coefficients["as_", "t value"]), tolerance = 1e-8)))
  chk("FL finds a planted effect (p < 0.01)", TRUE, f1$p < 0.01)
  pn <- vapply(1:300, function(i) fl(as.numeric(Xs[, 2:5] %*% c(0.5, -0.3, 0.2, 0) + rnorm(n)), as_, Xs, 500L)$p, 1)
  chk("FL null rate of p < 0.05 within [0.02, 0.09]", TRUE, mean(pn < 0.05) >= 0.02 && mean(pn < 0.05) <= 0.09)
})

RES <- list(); FK <- list()
run <- function(analysis, ligand, arm, y, a, X, fk, note = "") {
  r <- fl(y, a, X, N_PERM, fk$A[, names(a), drop = FALSE])
  p95 <- if (length(r$t_fake) >= 20L) unname(quantile(r$t_fake, WS_B0_FAKE_PCTL)) else NA_real_
  pass <- is.finite(r$t) && r$p < WS_B0_P_ONESIDED && is.finite(p95) && r$t > p95
  RES[[length(RES) + 1L]] <<- data.table(analysis, ligand, arm, n = r$n, df = r$df, covariates = paste(r$covs, collapse = "+"),
    t = round(r$t, 3), p_FL = signif(r$p, 4), n_perm = N_PERM, fake_n = length(r$t_fake), fake_tol = fk$tol,
    fake_p95 = round(p95, 3), pass = pass, note = note)
  FK[[length(FK) + 1L]] <<- data.table(analysis, ligand, arm, gene = fk$genes, t = round(r$t_fake, 4))
  ts(sprintf("    %-18s %-5s %-4s n=%2d | t=%6.2f p=%.4f | fake P95=%6.2f | %s", analysis, ligand, arm, r$n, r$t, r$p, p95, if (pass) "PASS" else "-"))
  fwrite(rbindlist(RES), file.path(TABWS, paste0("ws_b0_results", SFX, ".csv")))   # incremental
}

## -- 4. footprints and models ----------------------------------------------------------------------
ts("[4] footprints", if (SMOKE) " (smoke: no model is fitted on real outcomes)" else " and models")
for (AR in c("main", "sens")) {
  D <- IN[arm == AR & eligible == TRUE]; smp <- D$sample
  if (length(smp) < (if (SMOKE) 3L else 10L)) { ts("    ", AR, ": only ", length(smp), " eligible samples -> skipped"); next }
  Mr <- pb_mat(AR, "sum_recv", smp); X <- Xmat(D)
  for (i in seq_len(nrow(WS_B0_PAIRS))) {
    L <- WS_B0_PAIRS$ligand[i]; ck <- WS_B0_PAIRS$cytosig[i]; drop <- c(L, RECEPTORS[[L]])
    y <- fp_cytosig(Mr, ck, drop); a <- lig(AR, "sum_tnk", L, smp)
    IN[arm == AR & sample %in% smp, c(paste0("fp_", L), paste0("lig_", L)) := list(y[sample], a[sample])]
    if (i == 1L && AR == "main") {                          # decoupleR ulm = lm over all matrix genes (zero weight off-signature)
      w <- numeric(nrow(Mr)); net <- CS[gene %in% rownames(Mr) & !gene %in% drop & get(ck) != 0]
      w[match(net$gene, rownames(Mr))] <- net[[ck]]
      chk("decoupleR ulm score equals the lm t", TRUE, isTRUE(all.equal(unname(y[smp[1]]), summary(lm(Mr[, smp[1]] ~ w))$coefficients[2, 3], tolerance = 1e-6)))
    }
    if (SMOKE) next
    fk <- fake_pool(AR, "sum_tnk", smp, L, ck)
    run("primary", L, AR, y[smp], a[smp], X, fk)
    if (AR != "main") next
    run("progeny", L, AR, fp_progeny(Mr, WS_B0_PAIRS$progeny[i], drop)[smp], a[smp], X, fk)
    run("ligand_x_fraction", L, AR, y[smp], lig(AR, "sum_tnk", L, smp, frac = TRUE), X, fk)
    run("no_IEG", L, AR, fp_cytosig(Mr, ck, c(drop, WS_B0_IEG))[smp], a[smp], X, fk)
    sb <- D[n_b >= WS_B0_MIN_B, sample]
    if (length(sb) >= 10L) {
      yb <- fp_cytosig(pb_mat(AR, "sum_b", sb), ck, drop)
      run("B_tier", L, AR, y[sb], a[sb], cbind(Xmat(D[sample %in% sb]), b_fp = yb[sb]), fk, note = paste(length(sb), "samples with >= 30 B_Plasma"))
    }
  }
  if (!SMOKE && AR == "main") {                             # TNF sent by Mono_DC: clone-internal, reported only
    sm2 <- D[n_mono >= WS_B0_MIN_SENDER, sample]
    if (length(sm2) >= 10L) {
      Mm <- pb_mat(AR, "sum_recv", sm2); y2 <- fp_cytosig(Mm, "TNFA", c("TNF", RECEPTORS[["TNF"]]))
      run("clone_internal_MonoDC", "TNF", AR, y2[sm2], lig(AR, "sum_mono", "TNF", sm2), Xmat(D[sample %in% sm2]),
          fake_pool(AR, "sum_mono", sm2, "TNF", "TNFA"), note = "sender Mono_DC, labelled clone-internal")
    }
  }
}
fwrite(IN, file.path(TABWS, paste0("ws_b0_inputs", SFX, ".csv")))
if (GATE && !SMOKE) {
  bl <- IN[arm == "main" & eligible == TRUE & is.finite(blast_pct_clinical)]
  if (nrow(bl) >= 5L) ts(sprintf("    prog_frac vs clinical blast%%: Spearman %.2f (n = %d)", cor(bl$prog_frac, bl$blast_pct_clinical, method = "spearman"), nrow(bl)))
}

## -- 5. verdict, GATES.md, checks ------------------------------------------------------------------
if (!SMOKE && length(RES)) {
  R <- rbindlist(RES); fwrite(rbindlist(FK), file.path(TABWS, paste0("ws_b0_fakes", SFX, ".csv.gz")), compress = "gzip")
  tier_of <- function(x) { p <- setNames(x$pass, x$ligand)
    if (isTRUE(p["IFNG"]) && sum(p, na.rm = TRUE) >= 2) "pass" else
    if (!isTRUE(p["IFNG"]) && isTRUE(p["TNF"]) && isTRUE(p["TGFB1"])) "conditional" else "fail" }
  tm <- tier_of(R[analysis == "primary" & arm == "main"]); tsn <- tier_of(R[analysis == "primary" & arm == "sens"])
  tier <- if (tm %in% c("pass", "conditional") && tsn != tm) "one_arm" else tm
  btl <- R[analysis == "B_tier", .(ligand, b_label = fifelse(p_FL < WS_B0_P_ONESIDED, "leukaemia-specific", "marrow-wide"))]
  pm  <- R[analysis == "primary" & arm == "main"]
  TIER <- data.table(date = format(Sys.Date()), dataset = DATASET, gate = GATE, n_perm = N_PERM,
    IFNG = pm[ligand == "IFNG", pass], TNF = pm[ligand == "TNF", pass], TGFB1 = pm[ligand == "TGFB1", pass],
    tier_main = tm, tier_sens = tsn, tier = if (GATE) tier else "report only",
    meaning = if (GATE) WS_B0_TIER[[tier]] else "direction report, never a gate")
  fwrite(TIER, file.path(TABWS, paste0("ws_b0_tier", SFX, ".csv"))); print(TIER)
  if (GATE) {
    block <- c("<!-- B0 result start -->", sprintf("## GATE_B0 footprint positive control -- RESULT %s", format(Sys.Date())), "",
      "Script: scripts/ws/21_gate_b0_footprint.R | config: scripts/ws/config_ws.R section 10 (frozen 2026-10-16, commit 4c54f6f)", "",
      "| analysis | ligand | arm | n | t | p (Freedman-Lane) | fake P95 | pass |", "|---|---|---|---|---|---|---|---|",
      R[, sprintf("| %s | %s | %s | %d | %.2f | %.4f | %.2f | %s |", analysis, ligand, arm, n, t, p_FL, fake_p95, pass)], "",
      sprintf("Primary, main arm -- IFNG: %s | TNF: %s | TGFB1: %s | sens-arm tier: %s", TIER$IFNG, TIER$TNF, TIER$TGFB1, tsn),
      if (nrow(btl)) paste("B-cell tier (reported):", paste(btl$ligand, btl$b_label, sep = " ", collapse = "; ")) else "B-cell tier: not run",
      sprintf("**Verdict: %s** -- %s", toupper(tier), TIER$meaning), "<!-- B0 result end -->")
    old <- readLines(GATES); s0 <- grep("<!-- B0 result start -->", old); e0 <- grep("<!-- B0 result end -->", old)
    writeLines(if (length(s0) && length(e0)) c(old[seq_len(s0 - 1)], block, old[-seq_len(e0)]) else c(old, "", block), GATES)
    ts("    GATES.md B0 result block written")
  }
}
C <- rbindlist(CHK); fwrite(C, file.path(TABWS, paste0("ws_b0_checks", SFX, ".csv")))
ts(sprintf("[5] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", if (SMOKE) "smoke run: caches, footprints and known-answer checks only" else if (GATE) paste("tier =", tier) else paste("direction report on", DATASET))
