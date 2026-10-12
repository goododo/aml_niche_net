#!/usr/bin/env Rscript
# INPUT  : scripts/ws/config_ws.R                                   (frozen thresholds; guard must pass)
#          /LARGE1/.../02_seurat_objects/04_annotated/GSE185381/<sample>.rds   (42 eligible samples)
#          results/tables/ws/ws_doublet_scores_percell.csv.gz             (scDblFinder + DoubletFinder, per cell)
#          /LARGE1/gr10634/gaozy/reference/h.all.v2024.1.Hs.symbols.gmt   (HALLMARK, pinned; see .PROVENANCE)
#          CellChatDB.human (CellChat package)                            (fake-ligand exclusion list)
#          results/tables/01_preprocess/00_curated_manifest.csv          (eligibility, is_multiplexed)
# OUTPUT : results/tables/ws/ws_t3c_results.csv          one row per (combo, arm): z, p_perm, fake P95, recovered, diagnostic
#          results/tables/ws/ws_t3c_pairs.csv            one row per ligand pair: recovered under each arm
#          results/tables/ws/ws_t3c_tier.csv             the GATE_W1.1 verdict per 5.6 item 8
#          results/tables/ws/ws_t3c_gene_effects.csv.gz  b[g] per (combo, arm, gene), observed only
#          results/tables/ws/ws_t3c_checks.csv
#          scripts/ws/GATES.md                           W1.1 block written/replaced between markers
#          /LARGE1/.../aml_niche_net/ws_t3c_cache/<sample>.rds   stage-1 cache per patient (kept)
# what it does: T3c, GATE_W1.1 -- the positive control for the two-stage test. Asks whether three
#   known cytokine -> response pairs (IFNG primary; TNF, TGFB1 auxiliary) are RECOVERED when the
#   sender is a different compartment from the receiver (5.4 step 9). Every threshold is read from
#   config_ws.R and nothing here chooses a number after seeing a result.
#   Stage 1 (per patient, per receiver bin, 5.4 step 4): closed-form linear model. Y = log-normalised
#   expression of the receiver's own response genes (5.4 step 0); Z = [1, CDR, bmm_fine, Phase]
#   (+ scDblFinder score in the MAIN arm); receivers = bin cells NOT expressing the ligand (step 3);
#   beta[g] and se[g] for "receptor detected" after residualising on Z. Z columns with zero variance
#   are dropped per fit (WS_Z_DROP_ZERO_VARIANCE) and rank is taken from the QR, so a bmm_fine level
#   present in one cell cannot make Z'Z singular.
#   Stage 2 (across patients, 5.4 step 5): per gene, WLS beta[g,p] ~ A[L,p] + A_other[L,p] +
#   log(n_cells) + multiplexed + receiver_share, weights 1/se^2 -- per gene, exactly, via the
#   Frisch-Waugh partialling: for each gene the projection off the non-ligand covariates is
#   precomputed once, so a permutation of A costs one vectorised pass over all genes rather than
#   8,000 regressions. Gene-set statistic = rank-sum z of b[g] for genes in the HALLMARK set vs the
#   rest (5.6 item 4). Significance = patient-level permutation of A[L,.], one-sided. The real z must
#   also exceed the 95th percentile of z under WS_FAKE_N fake ligands drawn from the SAME sender
#   compartment with matched detection and inter-patient IQR (5.6 item 4).
#   Fake-ligand approximation (disclosed): fakes reuse stage-1 fits made on ALL receiver cells
#   rather than re-excluding each fake's positive cells -- the approximation 5.4 step 4 allows; a
#   cross-compartment ligand is seldom expressed in the receiver bin, so the two cell sets differ by
#   a handful of cells.
#   Two doublet arms per 5.3 item 1: MAIN keeps every cell and puts the scDblFinder score in Z;
#   SENS removes cells flagged by either caller. GATE_W1.1 requires agreement.
#   Diagnostic (5.6 item 9, reported never scored): per patient, mean log-expression of the gene set
#   in receiver cells vs A[L,p], Spearman.
#   The expensive, reusable part (stage-1 fits per patient) is cached to disk BEFORE anything that
#   can fail downstream -- three earlier jobs lost long reads to a check placed above the write.
# This script deletes nothing and does not modify the 04_annotated objects.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
.a <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
SMOKE <- "--smoke" %in% .a
ARMS  <- strsplit(.get("arms", if (SMOKE) "main" else "main,sens"), ",")[[1]]
SFX   <- if (SMOKE) "_smoke" else ""

suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
ws_assert_frozen()
set.seed(SEED)

AN     <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated/GSE185381")
GMT    <- file.path(REF_DIR, "h.all.v2024.1.Hs.symbols.gmt")
CACHE  <- file.path(LARGE1_DIR, paste0("ws_t3c_cache", SFX)); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
TABWS  <- file.path(TAB_DIR, "ws"); dir.create(TABWS, showWarnings = FALSE, recursive = TRUE)
GATES  <- here::here("scripts", "ws", "GATES.md")
N_PERM <- if (SMOKE) 50L else 1000L
N_FAKE <- if (SMOKE) 20L else WS_FAKE_N
GENE_CAP <- if (SMOKE) 600L else Inf
COMBOS <- ws_w11_runnable(); if (SMOKE) COMBOS <- COMBOS[COMBOS$row == 1L, ]
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-60s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}
ts(sprintf("T3c | smoke=%s | arms=%s | perms=%d | fakes=%d | combos=%s", SMOKE, paste(ARMS, collapse="+"), N_PERM, N_FAKE, paste(COMBOS$row, collapse=",")))

## -- 1. eligibility, gene sets, exclusion list ------------------------------------------------
ts("[1] inputs")
man  <- fread(file.path(TAB_DIR, "01_preprocess", "00_curated_manifest.csv"))
EL   <- man[dataset == "GSE185381" & disease == "AML" & timepoint == "Diagnosis" & tissue_r == "BM" &
            sorting_r %in% c("viability_only","none"), .(sample, multiplexed = is_multiplexed %in% c(TRUE, "TRUE"))]
chk("eligible samples", WS_W11_ELIGIBLE_SAMPLES, nrow(EL))
gmt <- strsplit(readLines(GMT), "\t"); names(gmt) <- vapply(gmt, `[`, "", 1L); gmt <- lapply(gmt, function(x) x[-(1:2)])
GS  <- lapply(WS_W11_GENESETS, function(n) gmt[[n]])
chk("all three HALLMARK sets found", 3L, sum(lengths(GS) > 0))
suppressPackageStartupMessages(library(CellChat)); data("CellChatDB.human", package = "CellChat")
cc  <- as.data.table(CellChatDB.human$interaction)
cplx <- CellChatDB.human$complex
expand <- function(r) if (r %in% rownames(cplx)) unlist(strsplit(gsub(" ", "", paste(cplx[r, ], collapse = ",")), ",")) else r
LR_GENES <- unique(c(unlist(lapply(cc$ligand, expand)), unlist(lapply(cc$receptor, expand)))); LR_GENES <- LR_GENES[LR_GENES != ""]
RECEPTORS_ALL <- unique(WS_W11_COMBOS$receptor)
DBL <- fread(file.path(TABWS, WS_DBL_SCORES_TABLE))[dataset == "GSE185381", .(cell, sc_score, dbl_union = (sc_class %in% c(TRUE,"TRUE")) | (df_class %in% c(TRUE,"TRUE")))]
setkey(DBL, cell)
chk("doublet scores cover the eligible samples", TRUE, uniqueN(fread(file.path(TABWS, WS_DBL_SCORES_TABLE))[dataset=="GSE185381", sample]) >= nrow(EL))

## -- 2. stage 1 per patient, cached ----------------------------------------------------------
ts("[2] stage 1 per patient (cache: ", CACHE, ")")
fit_stage1 <- function(Y, Z, Rpos) {
  # Y cells x genes (dense), Z cells x k (with intercept), Rpos logical. Returns beta, se per gene.
  v <- apply(Z, 2, var); keep <- c(TRUE, v[-1] > 0)               # drop zero-variance columns, keep intercept
  Z <- Z[, keep, drop = FALSE]
  q <- qr(Z); k <- q$rank; Q <- qr.Q(q)[, seq_len(k), drop = FALSE]
  Yres <- Y - Q %*% crossprod(Q, Y)
  r <- as.numeric(Rpos) - as.numeric(Q %*% crossprod(Q, as.numeric(Rpos)))
  rr <- sum(r * r); n <- nrow(Y)
  if (rr <= 1e-10 || n - k - 1 <= 1) return(list(beta = rep(NA_real_, ncol(Y)), se = rep(NA_real_, ncol(Y)), k = k))
  beta <- as.numeric(crossprod(r, Yres)) / rr
  sig2 <- (colSums(Yres^2) - beta^2 * rr) / (n - k - 1)
  list(beta = beta, se = sqrt(pmax(sig2, 0) / rr), k = k)
}
STAGE1 <- list()
for (i in seq_len(nrow(EL))) {
  sm <- EL$sample[i]; cf <- file.path(CACHE, paste0(sm, ".rds"))
  if (file.exists(cf)) { STAGE1[[sm]] <- readRDS(cf); next }
  f <- file.path(AN, paste0(sm, ".rds")); if (!file.exists(f)) next
  o <- readRDS(f); md <- o@meta.data
  b <- md$hierarchy_bin; if (WS_DROP_HIGH_ERROR) b[!is.na(md$high_error) & md$high_error == 1] <- NA
  cts <- LayerData(o, assay = "RNA", layer = "counts"); tot <- Matrix::colSums(cts)
  CDR <- as.numeric(Matrix::colSums(cts > 0) / nrow(cts))
  o <- NormalizeData(o, verbose = FALSE)
  o <- tryCatch(CellCycleScoring(o, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE),
                error = function(e) { o$Phase <- "NA"; o })
  Phase <- as.character(o$Phase); Phase[is.na(Phase)] <- "NA"
  bmm <- as.character(md$bmm_fine); bmm[is.na(bmm)] <- "NA"
  sc  <- DBL[J(colnames(o))]; sc_score <- sc$sc_score; sc_score[is.na(sc_score)] <- 0; dbl_union <- sc$dbl_union %in% TRUE
  # sender-side: mean CPM and detection fraction for ALL genes per compartment with >= 30 cells
  SEND <- list()
  for (bb in unique(b[!is.na(b)])) { j <- which(b == bb); if (length(j) < 30L) next
    m <- cts[, j, drop = FALSE]
    SEND[[bb]] <- data.table(gene = rownames(cts), n = length(j),
      mean_cpm = as.numeric(Matrix::rowSums(m) / sum(tot[j]) * 1e4), det = as.numeric(Matrix::rowSums(m > 0) / length(j))) }
  n_binned <- sum(!is.na(b))
  # receiver side per bin
  RECV <- list()
  for (rb in WS_W11_RECEIVERS) {
    j <- which(b == rb); if (length(j) < WS_MIN_CELLS_BIN) next
    oth <- setdiff(which(!is.na(b)), j)
    det <- Matrix::rowMeans(cts[, j, drop = FALSE] > 0)
    mu_in <- Matrix::rowMeans(cts[, j, drop = FALSE]); mu_out <- if (length(oth)) Matrix::rowMeans(cts[, oth, drop = FALSE]) else 0
    rg <- rownames(cts)[det >= WS_RESPONSE_GENE_MIN_DETECT & mu_out <= WS_RESPONSE_GENE_MAX_FOLD_OTHER * mu_in]
    rg <- setdiff(rg, RECEPTORS_ALL); if (is.finite(GENE_CAP)) rg <- head(rg, GENE_CAP)
    Y <- t(as.matrix(log1p(t(t(cts[rg, j, drop = FALSE]) / tot[j] * 1e4))))   # cells x genes, LogNormalize
    Zb <- model.matrix(~ CDR + bmm + Phase, data = data.frame(CDR = CDR[j], bmm = bmm[j], Phase = Phase[j]))
    fits <- list()
    for (ci in which(COMBOS$receiver == rb)) {
      L <- COMBOS$ligand[ci]; Rg <- COMBOS$receptor[ci]; key <- as.character(COMBOS$row[ci])
      if (!Rg %in% rownames(cts)) next
      Rpos <- as.numeric(cts[Rg, j]) > 0
      Lneg <- if (L %in% rownames(cts)) as.numeric(cts[L, j]) == 0 else rep(TRUE, length(j))
      for (arm in c("main","sens")) {
        cellsel <- if (arm == "main") Lneg else (Lneg & !dbl_union[j])
        Z <- if (arm == "main") cbind(Zb, sc_score = sc_score[j]) else Zb
        ok <- sum(Rpos[cellsel]) >= WS_ARM_PRIMARY && sum(!Rpos[cellsel]) >= WS_ARM_PRIMARY
        f1 <- if (ok) fit_stage1(Y[cellsel, , drop = FALSE], Z[cellsel, , drop = FALSE], Rpos[cellsel]) else NULL
        # all-receiver-cell fit for the fake-ligand approximation
        cellall <- if (arm == "main") rep(TRUE, length(j)) else !dbl_union[j]
        okA <- sum(Rpos[cellall]) >= WS_ARM_PRIMARY && sum(!Rpos[cellall]) >= WS_ARM_PRIMARY
        fA <- if (okA) fit_stage1(Y[cellall, , drop = FALSE], Z[cellall, , drop = FALSE], Rpos[cellall]) else NULL
        gs <- intersect(GS[[L]], rg)
        fits[[paste(key, arm)]] <- list(row = COMBOS$row[ci], arm = arm, genes = rg, n_recv = sum(cellsel),
          beta = if (!is.null(f1)) f1$beta else NULL, se = if (!is.null(f1)) f1$se else NULL,
          beta_all = if (!is.null(fA)) fA$beta else NULL, se_all = if (!is.null(fA)) fA$se else NULL,
          k = if (!is.null(f1)) f1$k else NA_integer_,
          set_score = if (length(gs) && sum(cellsel)) mean(Y[cellsel, gs, drop = FALSE]) else NA_real_)
      }
    }
    RECV[[rb]] <- list(n_bin = length(j), fits = fits)
  }
  STAGE1[[sm]] <- list(sample = sm, SEND = SEND, RECV = RECV, n_binned = n_binned)
  saveRDS(STAGE1[[sm]], cf)                       # cache FIRST
  ts(sprintf("    %-10s bins: %s | sender compartments: %d", sm,
     paste(sprintf("%s=%d", names(RECV), vapply(RECV, function(x) x$n_bin, 1L)), collapse = " "), length(SEND)))
  rm(list = intersect(c("o","cts","Y"), ls())); gc(FALSE)   # Y is absent when no receiver bin passed
}
chk("patients with stage-1 cache", nrow(EL), length(STAGE1))

## -- 3. stage 2 per combo x arm --------------------------------------------------------------
ts("[3] stage 2")
rank_z <- function(b, inset) {               # Wilcoxon rank-sum z, in-set vs out-set, positive = in-set larger
  ok <- is.finite(b); b <- b[ok]; inset <- inset[ok]
  n1 <- sum(inset); n2 <- sum(!inset); if (n1 < 5 || n2 < 5) return(NA_real_)
  R1 <- sum(rank(b)[inset]); U <- R1 - n1 * (n1 + 1) / 2
  (U - n1 * n2 / 2) / sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
}
sender_rank <- fread(file.path(TABWS, "ws_t3c_preflight_ligands.csv"))[, .(patients = .N, IQR = IQR(mean_expr)), by = .(gene, bin)][patients >= WS_MIN_PATIENTS][order(gene, -IQR)]
RES <- list(); GE <- list()
for (ci in seq_len(nrow(COMBOS))) for (arm in ARMS) {
  L <- COMBOS$ligand[ci]; S <- COMBOS$sender[ci]; Rb <- COMBOS$receiver[ci]; Rg <- COMBOS$receptor[ci]; key <- paste(COMBOS$row[ci], arm)
  other <- sender_rank[gene == L & bin != S & bin != Rb][1, bin]
  # patients: have the fit AND >= 30 sender cells (sender usable)
  pats <- names(STAGE1)[vapply(STAGE1, function(x) !is.null(x$RECV[[Rb]]) && !is.null(x$RECV[[Rb]]$fits[[key]]$beta) && !is.null(x$SEND[[S]]), TRUE)]
  if (length(pats) < WS_W11_MIN_PATIENTS) { ts(sprintf("    row %d %s: only %d patients -> not run", COMBOS$row[ci], arm, length(pats))); next }
  genes <- Reduce(union, lapply(pats, function(p) STAGE1[[p]]$RECV[[Rb]]$fits[[key]]$genes))
  getm <- function(field) { M <- matrix(NA_real_, length(pats), length(genes), dimnames = list(pats, genes))
    for (p in pats) { ff <- STAGE1[[p]]$RECV[[Rb]]$fits[[key]]; if (!is.null(ff[[field]])) M[p, ff$genes] <- ff[[field]] }; M }
  BETA <- getm("beta"); SE <- getm("se"); BETA_ALL <- getm("beta_all"); SE_ALL <- getm("se_all")
  A      <- vapply(pats, function(p) STAGE1[[p]]$SEND[[S]][gene == L, mean_cpm], 1)
  A_oth  <- vapply(pats, function(p) { s <- STAGE1[[p]]$SEND[[other]]; if (is.null(s)) 0 else s[gene == L, mean_cpm] }, 1)
  logn   <- log(vapply(pats, function(p) STAGE1[[p]]$RECV[[Rb]]$fits[[key]]$n_recv, 1))
  mux    <- as.numeric(EL$multiplexed[match(pats, EL$sample)])
  share  <- vapply(pats, function(p) STAGE1[[p]]$RECV[[Rb]]$n_bin / STAGE1[[p]]$n_binned, 1)
  # exact per-gene WLS via Frisch-Waugh: projection off O precomputed per gene
  build <- function(BETA, SE, O) {
    W <- 1 / SE^2; W[!is.finite(W)] <- 0; W[is.na(BETA)] <- 0
    keep <- colSums(W > 0) >= WS_W11_MIN_PATIENTS
    W <- W[, keep, drop = FALSE]; B <- BETA[, keep, drop = FALSE]; B[is.na(B)] <- 0
    n <- nrow(B); G <- ncol(B)
    Btil <- matrix(0, n, G); S_j <- vector("list", n); for (j in seq_len(n)) S_j[[j]] <- matrix(0, n, G)
    good <- rep(TRUE, G)
    for (g in seq_len(G)) {
      w <- W[, g]; OW <- O * w; M <- crossprod(O, OW)
      if (rcond(M) < 1e-10) { good[g] <- FALSE; next }
      H <- O %*% solve(M, t(OW)); P <- diag(n) - H
      Btil[, g] <- P %*% B[, g]; for (j in seq_len(n)) S_j[[j]][, g] <- P[, j]
    }
    list(W = W[, good, drop = FALSE], Btil = Btil[, good, drop = FALSE], S_j = lapply(S_j, function(m) m[, good, drop = FALSE]),
         genes = colnames(B)[good], n = n)
  }
  bcoef <- function(st, Avec) {
    Atil <- Reduce(`+`, Map(function(Sj, a) Sj * a, st$S_j, Avec))
    colSums(Atil * st$W * st$Btil) / colSums(Atil * st$W * Atil)
  }
  O  <- cbind(1, A_oth, logn, mux, share)
  st <- build(BETA, SE, O)
  inset <- st$genes %in% GS[[L]]
  b_obs <- bcoef(st, A); z_obs <- rank_z(b_obs, inset)
  z_perm <- vapply(seq_len(N_PERM), function(i) rank_z(bcoef(st, sample(A)), inset), 1)
  p_perm <- (1 + sum(z_perm >= z_obs, na.rm = TRUE)) / (1 + sum(is.finite(z_perm)))
  # fake ligands from the SAME sender compartment, matched on detection and inter-patient IQR
  sendtab <- rbindlist(lapply(pats, function(p) STAGE1[[p]]$SEND[[S]][, .(sample = p, gene, mean_cpm, det)]))
  cand <- sendtab[, .(det = mean(det), IQR = IQR(mean_cpm), nz = sum(mean_cpm > 0)), by = gene]
  ref  <- cand[gene == L]
  excl <- unique(c(LR_GENES, GS[[L]], L, RECEPTORS_ALL))
  pool <- function(tol) cand[!gene %in% excl & nz >= 3 & abs(det - ref$det) <= tol * ref$det & abs(IQR - ref$IQR) <= tol * ref$IQR, gene]
  fk <- pool(WS_FAKE_TOL); tol_used <- WS_FAKE_TOL
  if (length(fk) < WS_FAKE_MIN_POOL) { fk <- pool(WS_FAKE_TOL_RELAXED); tol_used <- WS_FAKE_TOL_RELAXED }
  set.seed(WS_FAKE_SEED); fk <- if (length(fk) > N_FAKE) sample(fk, N_FAKE) else fk
  stA <- build(BETA_ALL, SE_ALL, O)      # approximation: all-receiver-cell fits, O with the REAL A_other
  z_fake <- vapply(fk, function(fg) {
    Af <- vapply(pats, function(p) STAGE1[[p]]$SEND[[S]][gene == fg, mean_cpm], 1)
    Aof <- vapply(pats, function(p) { s <- STAGE1[[p]]$SEND[[other]]; if (is.null(s)) 0 else s[gene == fg, mean_cpm] }, 1)
    stF <- build(BETA_ALL, SE_ALL, cbind(1, Aof, logn, mux, share))
    rank_z(bcoef(stF, Af), stF$genes %in% GS[[L]])
  }, 1)
  fake_p95 <- if (sum(is.finite(z_fake)) >= 5) quantile(z_fake, WS_W11_FAKE_PCTL, na.rm = TRUE) else NA_real_
  recovered <- is.finite(z_obs) && p_perm < WS_W11_P_ONESIDED && is.finite(fake_p95) && z_obs > fake_p95
  # diagnostic
  setsc <- vapply(pats, function(p) STAGE1[[p]]$RECV[[Rb]]$fits[[key]]$set_score, 1)
  dg <- suppressWarnings(cor.test(setsc, A, method = "spearman"))
  RES[[key]] <- data.table(row = COMBOS$row[ci], ligand = L, sender = S, receiver = Rb, receptor = Rg, arm = arm,
    n_patients = length(pats), n_genes = length(st$genes), n_set_genes = sum(inset), A_other_from = other,
    z_obs = round(z_obs, 3), p_perm = round(p_perm, 4), n_perm = N_PERM,
    fake_n = length(fk), fake_tol = tol_used, fake_p95 = round(fake_p95, 3), fake_median = round(median(z_fake, na.rm = TRUE), 3),
    recovered = recovered, diag_spearman = round(unname(dg$estimate), 3), diag_p = signif(dg$p.value, 3))
  GE[[key]] <- data.table(row = COMBOS$row[ci], arm = arm, gene = st$genes, b = round(b_obs, 5), in_set = inset)
  ts(sprintf("    row %d %-5s %-4s %-8s->%-8s %-8s | n=%2d genes=%5d | z=%6.2f p=%.4f | fake P95=%6.2f (n=%d, tol=%.2f) | %s | diag rho=%.2f",
     COMBOS$row[ci], L, arm, S, Rb, Rg, length(pats), length(st$genes), z_obs, p_perm, fake_p95, length(fk), tol_used,
     if (recovered) "RECOVERED" else "not recovered", unname(dg$estimate)))
  fwrite(rbindlist(RES), file.path(TABWS, paste0("ws_t3c_results", SFX, ".csv")))   # incremental
}
RESULTS <- rbindlist(RES); fwrite(rbindlist(GE), file.path(TABWS, paste0("ws_t3c_gene_effects", SFX, ".csv.gz")), compress = "gzip")

## -- 4. pairs, tier, GATES.md ---------------------------------------------------------------
ts("[4] pairs and tier")
PAIRS <- RESULTS[, .(recovered = any(recovered)), by = .(ligand, arm)]
PW <- dcast(PAIRS, ligand ~ arm, value.var = "recovered")
fwrite(PW, file.path(TABWS, paste0("ws_t3c_pairs", SFX, ".csv"))); print(PW)
rec_main <- if ("main" %in% names(PW)) setNames(PW$main, PW$ligand) else NULL
rec_sens <- if ("sens" %in% names(PW)) setNames(PW$sens, PW$ligand) else NULL
both <- !is.null(rec_main) && !is.null(rec_sens)
agree <- if (both) all(rec_main[names(rec_sens)] == rec_sens, na.rm = TRUE) else NA
rec <- if (!is.null(rec_main)) rec_main else rec_sens
tier <- if (!both) "incomplete (one arm only)" else if (!agree) "one_arm" else if (isTRUE(rec["IFNG"]) && sum(rec, na.rm = TRUE) >= 2) "pass" else
        if (!isTRUE(rec["IFNG"]) && isTRUE(rec["TNF"]) && isTRUE(rec["TGFB1"])) "conditional" else "fail"
TIER <- data.table(date = format(Sys.Date()), smoke = SMOKE, arms = paste(ARMS, collapse = "+"), n_perm = N_PERM, n_fake = N_FAKE,
  IFNG = rec["IFNG"], TNF = rec["TNF"], TGFB1 = rec["TGFB1"], arms_agree = agree, tier = tier,
  meaning = if (tier %in% names(WS_W11_TIER)) WS_W11_TIER[[tier]] else tier)
fwrite(TIER, file.path(TABWS, paste0("ws_t3c_tier", SFX, ".csv"))); print(TIER)
if (!SMOKE) {
  block <- c("<!-- W1.1 start -->",
    sprintf("## GATE_W1.1 positive control -- %s", format(Sys.Date())),
    "", "Script: scripts/ws/20_t3c_positive_control.R | config: scripts/ws/config_ws.R (frozen 2026-10-12)", "",
    "| row | ligand | sender -> receiver | receptor | arm | n | z | p_perm | fake P95 | recovered | diag rho |", "|---|---|---|---|---|---|---|---|---|---|---|",
    RESULTS[, sprintf("| %d | %s | %s -> %s | %s | %s | %d | %.2f | %.4f | %.2f | %s | %.2f |", row, ligand, sender, receiver, receptor, arm, n_patients, z_obs, p_perm, fake_p95, recovered, diag_spearman)],
    "", sprintf("Pairs recovered -- IFNG: %s | TNF: %s | TGFB1: %s | arms agree: %s", rec["IFNG"], rec["TNF"], rec["TGFB1"], agree),
    sprintf("**Verdict: %s** -- %s", toupper(tier), TIER$meaning), "<!-- W1.1 end -->")
  old <- if (file.exists(GATES)) readLines(GATES) else c("# GATES.md -- one entry per gate: date, numbers, pass or fail", "")
  s0 <- grep("<!-- W1.1 start -->", old); e0 <- grep("<!-- W1.1 end -->", old)
  new <- if (length(s0) && length(e0)) c(old[seq_len(s0 - 1)], block, old[-(seq_len(e0))]) else c(old, block)
  writeLines(new, GATES); ts("    GATES.md W1.1 block written")
}
C <- rbindlist(CHK); fwrite(C, file.path(TABWS, paste0("ws_t3c_checks", SFX, ".csv")))
ts(sprintf("[5] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] tier = ", tier)
