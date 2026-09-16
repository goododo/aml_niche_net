#!/usr/bin/env Rscript
# 03_permute.R ----
# PREREG 8.3. Decides whether Stage B's hit counts mean anything. Two permutation nulls and one
# planted-effect positive control, per analysed bin. Without this, a hit count is a number with no
# scale and a zero is indistinguishable from a pipeline that cannot detect anything.
#
# INPUT  : DIR_PSEUDOBULK/pb_sample_manifest.csv, PB_RDS_DIR/<ds>/<sample>__pseudobulk.rds,
#          DIR_PSEUDOBULK/depth/*.csv           (all from 01_aggregate.R)
#          DIR_PSEUDOBULK/de_summary.csv        (which bins passed the 8.2 gate, from 02)
# OUTPUT : DIR_PSEUDOBULK/perm_null.csv         hit count per permutation, per bin, per null
#          DIR_PSEUDOBULK/perm_summary.csv      observed vs null, min attainable p, permutable units
#          DIR_PSEUDOBULK/power.csv             realised power per bin + PREREG 8.3 classification
# Usage  : Rscript scripts/12_pseudobulk_de/03_permute.R [--n_perm=1000]
#
#   PREREG 8.3  N_PERM = 1000. TWO nulls, both reported:
#               (1) PRIMARY -- permute `arm` among samples WITHIN dataset. This matches the
#                   exchangeability the blocked model assumes (samples exchangeable CONDITIONAL on
#                   block), which is why it is primary.
#               (2) CONSERVATIVE -- permute WHOLE libraries within dataset, restricted to arm-pure
#                   libraries; mixed libraries held fixed. Needed because arm is a SAMPLE attribute,
#                   not a library one: 5 of GSE185381's 22 libraries hold both an AML and a healthy
#                   sample, so a pure library-level shuffle is undefined for them. Meanwhile all four
#                   of Control1/2/3/4 sit in lib__P02, so a naive sample shuffle inside that lane is
#                   not exchangeable either. Neither null alone is right; both are reported.
#               Report the full hit-count DISTRIBUTION, not its mean. 0-vs-0 carries no information.
#
#               POSITIVE CONTROL, numerically pre-registered: 200 genes from the MIDDLE EXPRESSION
#               TERTILE, logFC = 1.0 spiked into the AML arm, full pipeline rerun. Realised power =
#               fraction of those 200 recovered under the FULL hit definition (q < 0.05 AND
#               |logFC| >= 1). power < 0.50 -> bin is POWER-LIMITED and its null is uninformative,
#               not a null. power >= 0.80 -> the strong statement is available. In between, the null
#               is reported with the power attached and no strong statement.
#
#   PREREG 8.4  Per-test seeds via .perm_seed(label), NOT one set.seed(SEED) for the process. A
#               single stream makes a test's p-value depend on how many tests ran before it; that
#               bug is written up at 05_ccc/04_stemness_purity_sweep.R:95-110. Copied rather than
#               sourced because that file is a runnable analysis script, not a library.
#
# TWO IMPLEMENTATION CHOICES NOT FIXED BY THE PREREG, DECLARED HERE BEFORE THE RUN.
#   (a) The gene set and the voom mean-variance trend are computed ONCE on the observed design and
#       held fixed across permutations; only lmFit/eBayes are recomputed per permutation. The
#       mean-variance trend is a property of the counts, which permutation does not touch, and
#       refitting it 1000x per bin would cost ~6 h per bin for no change in what is being tested.
#   (b) duplicateCorrelation is estimated once on the observed design and its consensus reused
#       across permutations, for the same reason. Both choices make the null SLIGHTLY OPTIMISTIC
#       (the observed fit's weights are tuned to the observed labels), which biases toward finding
#       the observed count unremarkable -- i.e. against declaring a hit. That direction is the safe
#       one and is stated so it cannot be presented as neutral.

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(here); library(limma); library(edgeR); library(statmod)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "config_malignancy.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--n_perm", type = "integer", default = 1000L),
  make_option("--n_planted", type = "integer", default = 200L),
  make_option("--planted_lfc", type = "double", default = 1.0)
)))
N_PERM <- opt$n_perm; Q_CUT <- 0.05; LFC_CUT <- 1.0
MIN_CELLS <- CCC_MIN_CELLS_PER_OCCUPIED_BIN

# PREREG 8.4 -- verbatim from 05_ccc/04_stemness_purity_sweep.R:111-116. The polynomial rolling hash
# is collision-free on this label space; the obvious SEED + sum(utf8ToInt(s) * seq_along(...))
# template collides, and a collision hands two different tests the same permutations.
.perm_seed <- function(label) {
  M <- 2147483647; h <- 0
  for (ch in utf8ToInt(label)) h <- (h * 131 + ch) %% M
  as.integer((h + SEED) %% M)
}

## -------------------------------------------------------------- rebuild ----
MAN <- fread(file.path(DIR_PSEUDOBULK, "pb_sample_manifest.csv"))[in_analysis == TRUE]
DEP <- rbindlist(lapply(list.files(file.path(DIR_PSEUDOBULK, "depth"), full.names = TRUE), fread))
PB  <- lapply(seq_len(nrow(MAN)), function(i)
  readRDS(file.path(PB_RDS_DIR, MAN$dataset[i], paste0(MAN$sample[i], "__pseudobulk.rds"))))
names(PB) <- paste(MAN$dataset, MAN$sample, sep = "|")
DISC <- MAN[split == "Discovery"]
base <- Reduce(intersect, lapply(unique(DISC$dataset),
                                 function(d) PB[[paste(DISC[dataset == d][1]$dataset, DISC[dataset == d][1]$sample, sep = "|")]]$genes))
GO   <- fread(INFERCNV_GENE_ORDER, header = FALSE, col.names = c("gene", "chr", "start", "end"))
UNIV <- setdiff(base, union(GO[chr %in% c("chrY", "Y"), unique(gene)], "XIST"))

SUMM <- fread(file.path(DIR_PSEUDOBULK, "de_summary.csv"))[gate == "PASS"]
message("[bins] ", nrow(SUMM), " passed the 8.2 gate: ", paste(SUMM$hierarchy_bin, collapse = ", "))

bin_matrix <- function(bin, roster) {
  keep <- DEP[hierarchy_bin == bin & n_cells >= MIN_CELLS,
              .(dataset, sample, n_cells, median_umi)][roster, on = c("dataset", "sample"), nomatch = 0]
  m <- vapply(seq_len(nrow(keep)), function(i)
    PB[[paste(keep$dataset[i], keep$sample[i], sep = "|")]]$counts[UNIV, bin], numeric(length(UNIV)))
  dimnames(m) <- list(UNIV, paste(keep$dataset, keep$sample, sep = "|"))
  keep[, arm := relevel(factor(arm), ref = "healthy")]
  keep[, `:=`(dataset = factor(dataset), library_id = factor(library_id))]
  list(counts = m, meta = keep)
}

## --------------------------------------------------------- permutations ----
## PRIMARY null: shuffle arm among samples, within dataset. Preserves each dataset's arm counts.
perm_sample <- function(md) {
  a <- as.character(md$arm)
  for (d in unique(md$dataset)) { i <- which(md$dataset == d); a[i] <- sample(a[i]) }
  factor(a, levels = levels(md$arm))
}

## CONSERVATIVE null: shuffle arm among ARM-PURE libraries, within dataset; mixed libraries keep
## their real per-sample labels. A library's samples all move together.
lib_purity <- function(md) md[, .(n = .N, n_arm = uniqueN(arm), arm1 = as.character(arm)[1]),
                              by = .(dataset, library_id)][, pure := n_arm == 1L][]
perm_library <- function(md, LP) {
  a <- as.character(md$arm)
  for (d in unique(md$dataset)) {
    P <- LP[dataset == d & pure == TRUE]
    if (nrow(P) < 2L) next
    newlab <- sample(P$arm1)
    for (k in seq_len(nrow(P))) a[md$dataset == d & md$library_id == P$library_id[k]] <- newlab[k]
  }
  factor(a, levels = levels(md$arm))
}

## ------------------------------------------------------------- the fits ----
## Observed fit; returns the frozen pieces the permutations reuse (see declaration (a)/(b) above).
observed_fit <- function(B, tier) {
  md <- B$meta
  rhs <- if (tier == "primary") "~ dataset + arm" else "~ arm"
  design <- model.matrix(as.formula(rhs), data = md)
  y <- DGEList(B$counts); keep <- filterByExpr(y, design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  if (tier == "primary") {
    v <- voom(y, design); c1 <- duplicateCorrelation(v, design, block = md$library_id)
    v <- voom(y, design, block = md$library_id, correlation = c1$consensus)
    c2 <- duplicateCorrelation(v, design, block = md$library_id); cons <- c2$consensus
    fit <- eBayes(lmFit(v, design, block = md$library_id, correlation = cons))
  } else { v <- voom(y, design); cons <- NA_real_; fit <- eBayes(lmFit(v, design)) }
  tt <- as.data.table(topTable(fit, coef = "armAML", number = Inf), keep.rownames = "gene")
  list(v = v, y = y, cons = cons, tier = tier, rhs = rhs,
       n_hits = sum(tt$adj.P.Val < Q_CUT & abs(tt$logFC) >= LFC_CUT), tt = tt)
}

## One refit under a supplied arm vector, reusing the frozen voom object.
refit_hits <- function(O, md, arm_new, E = NULL) {
  md2 <- copy(md)[, arm := arm_new]
  design <- model.matrix(as.formula(O$rhs), data = md2)
  if (!"armAML" %in% colnames(design)) return(NA_integer_)   # degenerate draw: one arm emptied
  v <- O$v
  if (!is.null(E)) v$E <- E                                   # planted-effect variant
  fit <- if (O$tier == "primary") lmFit(v, design, block = md$library_id, correlation = O$cons)
         else lmFit(v, design)
  fit <- eBayes(fit)
  tt <- as.data.table(topTable(fit, coef = "armAML", number = Inf), keep.rownames = "gene")
  list(n = sum(tt$adj.P.Val < Q_CUT & abs(tt$logFC) >= LFC_CUT), tt = tt)
}

## ------------------------------------------------------------------ run ----
null_rows <- list(); summ_rows <- list(); pow_rows <- list()
for (i in seq_len(nrow(SUMM))) {
  bin <- SUMM$hierarchy_bin[i]; tier <- SUMM$tier[i]
  roster <- if (tier == "primary") DISC else DISC[dataset == "Chen2023"]
  B <- bin_matrix(bin, roster); md <- B$meta
  O <- observed_fit(B, tier)
  LP <- lib_purity(md)

  # PREREG 8.3 requires the permutable-unit count and the minimum attainable p PER BIN, not assumed.
  n_units_samp <- md[, .N, by = dataset][, prod(choose(N, 1))]   # placeholder, replaced below
  n_units_samp <- prod(vapply(split(md, md$dataset), function(x)
    choose(nrow(x), sum(x$arm == "healthy")), numeric(1)))
  pure_units <- LP[pure == TRUE, .N, by = dataset]
  n_units_lib <- prod(vapply(split(LP[pure == TRUE], LP[pure == TRUE]$dataset), function(x)
    choose(nrow(x), sum(x$arm1 == "healthy")), numeric(1)))

  for (nullname in c("sample_within_dataset", "library_arm_pure")) {
    set.seed(.perm_seed(paste0("perm|", nullname, "|", tier, "|", bin)))   # PREREG 8.4
    hits <- vapply(seq_len(N_PERM), function(k) {
      a <- if (nullname == "sample_within_dataset") perm_sample(md) else perm_library(md, LP)
      r <- refit_hits(O, md, a); if (is.list(r)) r$n else NA_integer_
    }, numeric(1))
    null_rows[[paste(bin, nullname)]] <- data.table(tier, hierarchy_bin = bin, null = nullname,
                                                    perm = seq_len(N_PERM), n_hits = hits)
    p_emp <- (1 + sum(hits >= O$n_hits, na.rm = TRUE)) / (1 + sum(!is.na(hits)))
    summ_rows[[paste(bin, nullname)]] <- data.table(
      tier, hierarchy_bin = bin, null = nullname, observed_hits = O$n_hits,
      null_median = median(hits, na.rm = TRUE), null_mean = round(mean(hits, na.rm = TRUE), 1),
      null_p95 = quantile(hits, 0.95, na.rm = TRUE), null_max = max(hits, na.rm = TRUE),
      p_empirical = p_emp, n_valid_perm = sum(!is.na(hits)),
      permutable_units = if (nullname == "sample_within_dataset") n_units_samp else n_units_lib,
      min_attainable_p = 1 / (1 + min(N_PERM, if (nullname == "sample_within_dataset") n_units_samp else n_units_lib)))
    message(sprintf("[%s/%s] %s: observed %d | null median %g mean %.1f p95 %g max %d | p_emp %.4f | units %g",
                    tier, bin, nullname, O$n_hits, median(hits, na.rm = TRUE), mean(hits, na.rm = TRUE),
                    quantile(hits, 0.95, na.rm = TRUE), max(hits, na.rm = TRUE), p_emp,
                    if (nullname == "sample_within_dataset") n_units_samp else n_units_lib))
  }

  ## PREREG 8.3 positive control: middle expression tertile, logFC = 1.0 into the AML arm.
  set.seed(.perm_seed(paste0("planted|", bin)))
  mu <- rowMeans(O$v$E)
  ter <- cut(rank(mu, ties.method = "first"), 3, labels = FALSE)
  pool <- rownames(O$v$E)[ter == 2L]
  planted <- sample(pool, min(opt$n_planted, length(pool)))
  E2 <- O$v$E
  E2[planted, md$arm == "AML"] <- E2[planted, md$arm == "AML"] + opt$planted_lfc   # voom E is log2-CPM
  r <- refit_hits(O, md, md$arm, E = E2)
  rec <- r$tt[gene %in% planted & adj.P.Val < Q_CUT & abs(logFC) >= LFC_CUT, .N]
  power <- rec / length(planted)
  cls <- if (power < 0.50) "POWER-LIMITED: null is uninformative, not a null"
         else if (power >= 0.80) "WELL-POWERED: strong statement available"
         else "INTERMEDIATE: report null with power attached, no strong statement"
  pow_rows[[bin]] <- data.table(tier, hierarchy_bin = bin, n_planted = length(planted),
                                planted_lfc = opt$planted_lfc, n_recovered = rec,
                                realised_power = round(power, 3), classification = cls,
                                observed_hits = O$n_hits)
  message(sprintf("[%s/%s] PLANTED lfc=%.1f n=%d -> recovered %d, power %.3f -- %s",
                  tier, bin, opt$planted_lfc, length(planted), rec, power, cls))
}

fwrite(rbindlist(null_rows), file.path(DIR_PSEUDOBULK, "perm_null.csv"))
PS <- rbindlist(summ_rows); fwrite(PS, file.path(DIR_PSEUDOBULK, "perm_summary.csv"))
PW <- rbindlist(pow_rows);  fwrite(PW, file.path(DIR_PSEUDOBULK, "power.csv"))
message("\n================ PREREG 8.3 ================"); print(PS); print(PW)
