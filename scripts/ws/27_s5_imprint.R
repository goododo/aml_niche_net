#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../aml_niche_net/ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds   (S4: raw sums per compartment and subtype)
#          results/tables/ws/ws_s4_pseudobulk_meta_GSE185381.csv, ws_s4_subtype_composition_GSE185381.csv
#          results/tables/ws/decontx/GSE185381/ws_s1_ambient_blacklist.csv          (S1: genes in >= 5 channels' top 100)
#          /LARGE1/gr10634/gaozy/reference/gencode_GRCh38_gene_order.txt           (chrY genes; XIST)
#          scripts/ws/PREREGISTRATION_S5_imprint.md (tag prereg-s5-v1) and config_ws.R sections 10-12
# OUTPUT : results/tables/ws/ws_s5_de_<compartment>.csv       every tested gene: logFC, logCPM, PValue, FDR, flags
#          results/tables/ws/ws_s5_imprint_genes.csv           the imprint (FDR < 0.10, |log2FC| > 0.5), up and down, per compartment
#          results/tables/ws/ws_s5_scores.csv                  sample x compartment: imprint score (4.1), status, covariates
#          results/tables/ws/ws_s5_known_answers.csv           planted-DE recovery, sex contrast, label-permutation counts
#          results/tables/ws/ws_s5_sensitivity.csv             limma-voom, IEG-removed, multiplexed covariate, subtype-level arms
#          results/tables/ws/ws_s5_checks.csv
# what it does: PREREGISTRATION_S5_imprint.md sections 3, 4.1 and 8 (items 1-3), nothing else. Per compartment (T, NK,
#   B, Mono) it contrasts the eligible adult AML patients with the eligible healthy donors of GSE185381 by edgeR
#   quasi-likelihood with the pre-registered design, removes the blacklists BEFORE filtering, calls the imprint
#   genes, and scores every eligible sample. The known answers (planted DE, sex contrast, label permutation) run
#   on the same code path; the sensitivity arms re-run the same contrast under the declared variations and are
#   reported, never used to pick genes. The paediatric cohort is not touched (prereg section 10).
#   usage: Rscript scripts/ws/27_s5_imprint.R [--smoke]     (--smoke: 20 permutations, 50 planted genes)
# This script deletes nothing and does not modify its inputs.

suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma) })
.a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% .a
suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
set.seed(SEED)
DATASET <- WS_S4_DEFINE; TABWS <- file.path(TAB_DIR, "ws")
PB   <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", DATASET, paste0("ws_pb_", DATASET, ".rds")))
META <- fread(file.path(TABWS, paste0("ws_s4_pseudobulk_meta_", DATASET, ".csv")))
COMP <- fread(file.path(TABWS, paste0("ws_s4_subtype_composition_", DATASET, ".csv")))
BL   <- fread(file.path(TABWS, "decontx", DATASET, "ws_s1_ambient_blacklist.csv"))
GO   <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))   # inferCNV gene order (config_malignancy.R)
FDR_CUT <- 0.10; LFC_CUT <- 0.5; MIN_IMPRINT <- 10L; N_PERM <- if (SMOKE) 20L else 200L; N_PLANT <- if (SMOKE) 50L else 200L
SUBTYPE_SHARE <- list(T = "^CD8", Mono = "^CD16 Mono$", NK = "^NK CD56high$", B = "^Immature B$")   # prereg 3.2
SEX_GENES <- c("XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY")
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-64s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- 1. blacklists (prereg 3.3) ---------------------------------------------------------------------------
G <- PB$genes
BLACK <- unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene],
                  "XIST", GO[chr == "chrY", gene], grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))
chk("ambient blacklist has 79 genes", 79L, BL[n_units >= 5, .N])
chk("chrY genes found in the gene-order file", TRUE, GO[chr == "chrY", .N] > 20)
ts(sprintf("[1] blacklist: %d genes (%d present on the gene axis)", length(BLACK), sum(BLACK %in% G)))

## -- 2. samples and covariates per compartment -----------------------------------------------------------------
META[, status := fifelse(role == "control", "healthy", fifelse(role == "AML_Dx" & paediatric == FALSE, "AML", NA_character_))]
META[, sex := fifelse(!is.na(sex_meta) & sex_meta != "", sex_meta, sex_inferred)]
share_of <- function(sm, g4) { x <- COMP[sample == sm & group4 == g4]; if (!nrow(x)) return(NA_real_); sum(x[grepl(SUBTYPE_SHARE[[g4]], fine), n]) / sum(x$n) }
design_table <- function(g4) {
  d <- META[group4 == g4 & eligible == TRUE & !is.na(status)]
  d[, subtype_share := vapply(sample, share_of, 1, g4 = g4)]
  d[, multiplexed := as.integer(n_channels >= 1 & grepl(";", channels) | grepl("lib__", library_id))]
  d <- d[order(status, sample), .(sample, status = factor(status, c("healthy", "AML")), hb = hb_share_nonery, gran = gran_share_lymph,
                                   sex = factor(sex), subtype_share, multiplexed, dbl_mean, n_cells = n)]
  bad <- d[is.na(hb) | is.na(gran) | is.na(sex) | sex %in% "unclear" | is.na(subtype_share), sample]   # model.matrix would drop them silently
  if (length(bad)) { ts("    dropped for a missing covariate: ", paste(bad, collapse = " ")); d <- d[!sample %in% bad]; d[, sex := droplevels(sex)] }
  d
}

## -- 3. the DE engine, one code path for everything ------------------------------------------------------------
fml <- function(terms) as.formula(paste("~", paste(terms, collapse = " + ")))      # explicit terms; no update()/dot formulas
fit_de <- function(counts, d, terms, coef, engine = "edgeR", drop_black = TRUE, return_y = FALSE) {
  keep_g <- if (drop_black) !(rownames(counts) %in% BLACK) else rep(TRUE, nrow(counts))
  y <- DGEList(counts[keep_g, d$sample, drop = FALSE], group = d$status)
  keep <- filterByExpr(y, group = d$status) & (rowMeans(cpm(y) > 1) >= 0.5)                 # prereg 3.1
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM")
  X <- model.matrix(fml(terms), data = as.data.frame(d)); stopifnot(nrow(X) == nrow(d), qr(X)$rank == ncol(X))
  cn <- grep(coef, colnames(X), value = TRUE)[1]
  if (engine == "edgeR") {
    y <- estimateDisp(y, X, robust = TRUE); f <- glmQLFit(y, X, robust = TRUE); t <- glmQLFTest(f, coef = cn)
    r <- as.data.table(topTags(t, n = Inf, sort.by = "none")$table, keep.rownames = "gene")[, .(gene, logFC, logCPM, PValue, FDR)]
  } else {
    v <- voom(y, X); f <- eBayes(lmFit(v, X)); r <- as.data.table(topTable(f, coef = cn, n = Inf, sort.by = "none"), keep.rownames = "gene")
    r <- r[, .(gene, logFC, logCPM = AveExpr, PValue = P.Value, FDR = adj.P.Val)]
  }
  r[, hit := FDR < FDR_CUT & abs(logFC) > LFC_CUT]
  if (return_y) list(table = r, y = y) else r
}
MAIN <- c("status", "hb", "gran", "sex", "subtype_share")

## -- 4. per compartment ----------------------------------------------------------------------------------------
DE <- list(); IMP <- list(); SC <- list(); KA <- list(); SENS <- list()
for (g4 in WS_S4_GROUPS) {
  d <- design_table(g4); cnt <- PB$raw[[g4]]
  ts(sprintf("[%s] %d AML vs %d healthy", g4, d[status == "AML", .N], d[status == "healthy", .N]))
  chk(paste(g4, "has >= 12 eligible AML and >= 9 healthy"), TRUE, d[status == "AML", .N] >= 12 && d[status == "healthy", .N] >= 9)
  terms <- if (d[, uniqueN(sex)] < 2) { ts("    one sex only: sex dropped from the design"); setdiff(MAIN, "sex") } else MAIN
  rr <- fit_de(cnt, d, terms, "statusAML", return_y = TRUE); r <- rr$table; y <- rr$y
  r[, `:=`(compartment = g4, ieg = gene %in% WS_B0_IEG, direction = fifelse(logFC > 0, "up", "down"))]
  DE[[g4]] <- r; fwrite(r[order(PValue)], file.path(TABWS, paste0("ws_s5_de_", g4, ".csv")))
  up <- r[hit == TRUE & logFC > 0, gene]; dn <- r[hit == TRUE & logFC < 0, gene]
  has_imprint <- length(up) + length(dn) >= MIN_IMPRINT
  IMP[[g4]] <- data.table(compartment = g4, gene = c(up, dn), direction = rep(c("up", "down"), c(length(up), length(dn))),
                          logFC = r[match(c(up, dn), gene), logFC], FDR = r[match(c(up, dn), gene), FDR], ieg = c(up, dn) %in% WS_B0_IEG,
                          n_tested = nrow(r), imprint_called = has_imprint)
  ts(sprintf("    tested %d genes | up %d, down %d | imprint %s", nrow(r), length(up), length(dn), if (has_imprint) "CALLED" else "NOT called (< 10 genes)"))
  # -- 4.1 scores: z of logCPM across the eligible samples, mean z(up) - mean z(down)
  L <- cpm(y, log = TRUE, prior.count = 1); Z <- t(scale(t(L)))
  sc <- if (has_imprint) colMeans(Z[intersect(up, rownames(Z)), , drop = FALSE]) - (if (length(dn)) colMeans(Z[intersect(dn, rownames(Z)), , drop = FALSE]) else 0) else rep(NA_real_, ncol(Z))
  SC[[g4]] <- cbind(d, data.table(compartment = g4, score = as.numeric(sc)[match(d$sample, colnames(Z))], n_up = length(up), n_down = length(dn)))
  # -- 8.1 planted DE on this compartment's real data (AML samples x 2 or x 0.5 for random kept genes)
  for (lfc in c(1, 1.5)) {
    pool <- setdiff(r$gene, c(up, dn)); pg <- sample(pool, min(N_PLANT, length(pool))); half <- seq_len(length(pg) %/% 2)
    cp <- cnt; aml <- d[status == "AML", sample]
    cp[pg[half], aml] <- round(cp[pg[half], aml] * 2^lfc); cp[pg[-half], aml] <- round(cp[pg[-half], aml] * 2^-lfc)
    rp <- fit_de(cp, d, terms, "statusAML"); rec <- rp[gene %in% pg]
    rec[, ok := FDR < FDR_CUT & abs(logFC) > LFC_CUT & sign(logFC) == ifelse(gene %in% pg[half], 1, -1)]
    KA[[length(KA) + 1]] <- data.table(compartment = g4, test = "planted_DE", planted_log2FC = lfc, n_planted = length(pg), recovered = mean(rec$ok),
                                       detail = sprintf("up %.2f, down %.2f", mean(rec[gene %in% pg[half], ok]), mean(rec[!gene %in% pg[half], ok])))
    ts(sprintf("    planted |log2FC| %.1f: recovered %.2f", lfc, mean(rec$ok)))
  }
  # -- 8.2 sex contrast on the unblacklisted gene set
  if (d[, uniqueN(sex)] == 2) {
    rs <- fit_de(cnt, d, c("sex", "status", "hb", "gran", "subtype_share"), "^sex", drop_black = FALSE)[order(PValue)]
    top10 <- head(rs$gene, 10); n_sex <- sum(SEX_GENES %in% top10)
    KA[[length(KA) + 1]] <- data.table(compartment = g4, test = "sex_contrast_top10", planted_log2FC = NA_real_, n_planted = NA_integer_, recovered = n_sex,
                                       detail = paste(intersect(top10, SEX_GENES), collapse = ","))
    chk(paste(g4, "sex contrast puts >= 3 sex genes in the top 10"), TRUE, n_sex >= 3)
  }
  # -- 8.3 label permutation: how many "imprint genes" does a shuffled status produce
  # a random relabelling that nearly inverts the true labels re-creates the real contrast, so this null is
  # bimodal when a real effect exists; that is why it is read at the 95th percentile over 200 draws, not the max
  np <- vapply(seq_len(N_PERM), function(i) { dp <- copy(d); dp[, status := sample(status)]; sum(fit_de(cnt, dp, terms, "statusAML")$hit) }, 1)
  KA[[length(KA) + 1]] <- data.table(compartment = g4, test = "label_permutation", planted_log2FC = NA_real_, n_planted = N_PERM,
                                     recovered = as.numeric(quantile(np, 0.95)), detail = sprintf("real %d; perm median %.0f, P95 %.0f, max %.0f", length(up) + length(dn), median(np), quantile(np, .95), max(np)))
  chk(paste(g4, "real imprint count exceeds the permutation P95"), TRUE, unname((length(up) + length(dn)) > quantile(np, 0.95)))   # unname: all.equal() rejects the "95%" name
  # -- sensitivity arms (reported only)
  arm <- function(name, rr) { ov <- length(intersect(rr[hit == TRUE, gene], c(up, dn)))
    SENS[[length(SENS) + 1]] <<- data.table(compartment = g4, arm = name, n_hits = sum(rr$hit), overlap_with_imprint = ov, imprint_size = length(up) + length(dn)) }
  arm("limma_voom", fit_de(cnt, d, terms, "statusAML", engine = "limma"))
  arm("IEG_removed", { ri <- copy(r); ri[gene %in% WS_B0_IEG, hit := FALSE]; ri })
  if (d[, uniqueN(multiplexed)] == 2) arm("plus_multiplexed", fit_de(cnt, d, c(terms, "multiplexed"), "statusAML"))
  if (d[!is.na(dbl_mean), .N] == nrow(d)) arm("plus_doublet_mean", fit_de(cnt, d, c(terms, "dbl_mean"), "statusAML"))
  for (fn in intersect(names(PB$fine), c("CD4 Naive", "CD8 Effector Memory 2", "CD14 Mono"))) if (g4 == (if (fn == "CD14 Mono") "Mono" else "T")) {
    ncell <- PB$cells[kept == TRUE & fine == fn, .N, by = sample]; ok <- ncell[N >= WS_S4_MIN_CELLS, sample]
    ds <- d[sample %in% ok]; if (ds[status == "AML", .N] >= 8 && ds[status == "healthy", .N] >= 5) {
      rf <- fit_de(PB$fine[[fn]], ds, setdiff(terms, "subtype_share"), "statusAML")
      gi <- rf[gene %in% c(up, dn)]                                                     # plain vectors: inside r[...] 'gene' would resolve to r's own column
      conc <- mean(sign(gi$logFC) == sign(r$logFC[match(gi$gene, r$gene)]))
      SENS[[length(SENS) + 1]] <- data.table(compartment = g4, arm = paste0("subtype:", fn), n_hits = sum(rf$hit), overlap_with_imprint = round(conc, 3), imprint_size = length(up) + length(dn)) } }
}
DEALL <- rbindlist(DE); IMP <- rbindlist(IMP); SC <- rbindlist(SC); KA <- rbindlist(KA); SENS <- rbindlist(SENS)
chk("no blacklisted gene among the imprint genes", 0L, IMP[gene %in% BLACK, .N])
chk("scores exist for every eligible sample of a compartment with an imprint", 0L, SC[compartment %in% IMP[imprint_called == TRUE, unique(compartment)] & is.na(score), .N])
fwrite(IMP, file.path(TABWS, "ws_s5_imprint_genes.csv")); fwrite(SC, file.path(TABWS, "ws_s5_scores.csv"))
fwrite(KA, file.path(TABWS, "ws_s5_known_answers.csv")); fwrite(SENS, file.path(TABWS, "ws_s5_sensitivity.csv"))
print(IMP[, .(n_genes = .N, up = sum(direction == "up"), down = sum(direction == "down"), called = imprint_called[1]), by = compartment])
print(SENS)
C <- rbindlist(CHK); fwrite(C, file.path(TABWS, "ws_s5_checks.csv"))
ts(sprintf("[done] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
