#!/usr/bin/env Rscript
# 09_cellstate_score.R ----
# Per-cell scores for the functional axes the node features do not currently carry.
#
# WHY. ccc_node_features.csv describes each (sample, bin) with frac_malignant, mean_stemness and
# n_cells. That says how much of a node is malignant and how big it is; it says nothing about what
# the node is DOING. A T_NK node of 121,067 cells is represented by a count. The Stromal node
# resolves 17 subtypes and carries no functional state at all. The panel here fills that in:
#
#   apoptosis            apoptotic_priming (BH3-only) + antiapoptotic (BCL2/MCL1/BCL2L1/...)
#                        The niche protects blasts by shifting BCL2 -> MCL1 dependence, which is
#                        the mechanism of venetoclax resistance. Absent from the pipeline entirely.
#   metabolism           oxphos / glycolysis / fatty_acid_oxidation
#                        LSCs are OXPHOS-dependent, and adipocytes fuel blasts by FAO. This cohort
#                        has 6,707 adipocytes -- its LARGEST stromal subtype -- so the
#                        adipocyte -> blast FAO axis is directly testable here, not borrowed.
#   tnk_function         cytotoxicity / exhaustion / naive_memory
#   stromal_state        sasp / senescence
#   interferon           interferon_I / interferon_II
#                        Dual use: real niche biology AND a known batch confounder, so it can serve
#                        as a covariate in 08_scoring rather than hiding inside the effect.
#   antigen_presentation mhc_ii / mhc_i   (HLA-DR loss is a documented AML immune-escape route and
#                        is the direct mechanistic link between a blast node and the T_NK node)
#   niche_retention      retention_ligand / retention_receptor (CXCL12-CXCR4, VCAM1-VLA4, KITLG-KIT)
#
# SCORING. UCell, not AddModuleScore -- see CELLSTATE_* in config_hierarchy.R for the measurement
# that motivates it. UCell ranks genes within each cell, so a score depends on that cell alone and
# is comparable across samples without joint normalisation.
#
# THE POSITIVE CONTROL IS THE POINT. A gene panel with a stale symbol, or a bin vocabulary that has
# drifted, still produces a perfectly ordinary-looking number for every cell. So each signature
# with an unambiguous home compartment (CELLSTATE_EXPECTED_BIN) must PEAK there, and the script
# fails if it does not. Signatures with no such expectation (metabolism, interferon) are scored but
# cannot be checked this way, and are reported as such rather than silently trusted.
#
# OUTPUT: HIER_PROJ_DIR/<dataset>/<sample>__cellstate_percell.csv
#           cell, sample, dataset, hierarchy_bin, malignant, <signature scores>, <single genes>
#
#   Rscript scripts/03_hierarchy/09_cellstate_score.R [--dataset X] [--limit N] [--force]

suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(UCell); library(Matrix) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = ""),
  make_option("--limit",   type = "integer",   default = 0L),
  make_option("--force",   action = "store_true", default = FALSE)
)))

SIG <- fread(CELLSTATE_TABLE_TSV)
stopifnot(all(c("signature", "gene", "family") %in% names(SIG)), nrow(SIG) > 0)
sigs <- split(SIG$gene, SIG$signature)
message(sprintf("[0] %d signatures over %d genes in %d families",
                length(sigs), uniqueN(SIG$gene), uniqueN(SIG$family)))

R <- qc_rds_roster(on_extra = "ignore")
if (nzchar(opt$dataset)) R <- R[dataset %in% strsplit(opt$dataset, ",")[[1]]]
if (opt$limit > 0) R <- head(R, opt$limit)
stopifnot(nrow(R) > 0)

dst_of <- function(ds, sid) file.path(CELLSTATE_DIR, ds, paste0(sid, "__cellstate_percell.csv"))
proj_of <- function(ds, sid) file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
cons_of <- function(ds, sid) file.path(DIR_MALIGNANCY, ds, paste0(sid, "__consensus_percell.csv"))

cov_report <- NULL
score_one <- function(rds, ds, sid, verbose) {
  seu <- readRDS(rds)
  seu <- NormalizeData(seu, verbose = FALSE)
  dat <- SeuratObject::LayerData(seu, assay = "RNA", layer = "data")
  g_have <- rownames(dat)

  # GENE COVERAGE IS REPORTED, NOT ASSUMED. Datasets differ in annotation; a signature reduced to
  # two surviving genes is not the signature, and the score would still be a number.
  cov <- data.table(signature = names(sigs),
                    n_total = vapply(sigs, length, 1L),
                    n_present = vapply(sigs, function(g) sum(g %in% g_have), 1L))
  cov[, frac := round(n_present / n_total, 3)][, `:=`(dataset = ds, sample = sid)]
  cov_report <<- rbind(cov_report, cov)
  if (verbose) { message("  [coverage] worst signatures in this sample:"); print(head(cov[order(frac)], 4)) }

  use <- sigs[vapply(sigs, function(g) sum(g %in% g_have) >= 3L, logical(1))]
  U <- UCell::ScoreSignatures_UCell(dat, features = use, maxRank = CELLSTATE_MAXRANK,
                                    ncores = 1, name = "")
  out <- data.table(cell = colnames(dat))
  for (nm in colnames(U)) out[, (nm) := as.numeric(U[, nm])]
  for (nm in setdiff(names(sigs), colnames(U))) out[, (nm) := NA_real_]

  # single genes: mean log-normalised expression, so the MCL1-vs-BCL2 balance is directly readable
  for (g in CELLSTATE_SINGLE_GENES) {
    col <- paste0("expr_", gsub("-", "_", g))
    out[, (col) := if (g %in% g_have) as.numeric(dat[g, ]) else NA_real_]
  }

  pf <- proj_of(ds, sid)
  bin <- if (file.exists(pf)) fread(pf, select = c("cell", "hierarchy_bin")) else NULL
  cf <- cons_of(ds, sid)
  mal <- if (file.exists(cf)) fread(cf, select = c("cell", "malignant")) else NULL
  if (!is.null(bin)) out <- merge(out, bin, by = "cell", all.x = TRUE)
  if (!is.null(mal)) out <- merge(out, mal, by = "cell", all.x = TRUE) else out[, malignant := NA_integer_]
  out[, `:=`(dataset = ds, sample = sid)]
  rm(seu, dat, U); gc(verbose = FALSE)
  out
}

for (i in seq_len(nrow(R))) {
  ds <- R$dataset[i]; sid <- R$sample[i]; dst <- dst_of(ds, sid)
  if (file.exists(dst) && !opt$force) { message(sprintf("[%d/%d] %s::%s done", i, nrow(R), ds, sid)); next }
  message(sprintf("[%d/%d] %s::%s", i, nrow(R), ds, sid))
  x <- tryCatch(score_one(R$rds[i], ds, sid, verbose = (i == 1)),
                error = function(e) { message("  [FAIL] ", conditionMessage(e)); NULL })
  if (is.null(x)) next
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  fwrite_safe(x, dst)
}

## ---------------------------------------------------------------------------------------------
## validation
## ---------------------------------------------------------------------------------------------
done <- R[file.exists(dst_of(dataset, sample))]
if (!nrow(done)) { message("[!] nothing scored"); quit(save = "no", status = 1) }
A <- rbindlist(lapply(seq_len(nrow(done)), function(i)
  fread(dst_of(done$dataset[i], done$sample[i]))), fill = TRUE)

cat(sprintf("\n================ scored %d cells across %d samples ================\n",
            nrow(A), uniqueN(A$sample)))

cat("\n================ gene coverage per signature ================\n")
# cov_report only covers samples SCORED THIS RUN -- a resumed run skips everything and would print
# an empty table, which reads like "no coverage problems" rather than "not measured".
if (!is.null(cov_report)) {
  CV <- cov_report[, .(n_total = n_total[1], min_present = min(n_present), min_frac = min(frac)), by = signature]
  setorder(CV, min_frac); print(CV)
  if (CV[min_frac < 0.5, .N]) cat("  [!] signatures above are reduced by >50% in at least one dataset\n")
} else {
  cat("  (not measured: every sample was resumed from disk. Re-run with --force to remeasure.)\n")
}
# This check does NOT depend on what ran this time. A signature that lost too many genes to score
# is written as all-NA for that sample, so an all-NA column is the symptom that survives to disk.
cat("\n  per-signature NA rate across every scored cell (a high rate = the panel is not measurable):\n")
NAR <- data.table(signature = intersect(names(sigs), names(A)))
NAR[, na_frac := vapply(signature, function(s) mean(is.na(A[[s]])), 1.0)]
NAR[, samples_all_NA := vapply(signature, function(s)
  A[, .(allna = all(is.na(get(s)))), by = sample][allna == TRUE, .N], 1L)]
setorder(NAR, -na_frac); print(NAR[na_frac > 0 | samples_all_NA > 0])
if (NAR[na_frac == 0 & samples_all_NA == 0, .N] == nrow(NAR)) cat("  every signature scored in every cell\n")

cat("\n================ POSITIVE CONTROLS (ordered contrasts against settled biology) ================\n")
FAIL <- 0L
sc_cols <- intersect(names(sigs), names(A))
BM <- A[!is.na(hierarchy_bin) & nzchar(hierarchy_bin)]
bin_mean <- function(sg, b) BM[hierarchy_bin == b, mean(get(sg), na.rm = TRUE)]
for (i in seq_len(nrow(CELLSTATE_CONTRASTS))) {
  cc <- CELLSTATE_CONTRASTS[i]
  if (!cc$signature %in% sc_cols) next
  hi <- bin_mean(cc$signature, cc$high_bin); lo <- bin_mean(cc$signature, cc$low_bin)
  if (is.na(hi) || is.na(lo)) { cat(sprintf("  [SKIP] %-18s %s vs %s -- a bin is absent here\n",
                                            cc$signature, cc$high_bin, cc$low_bin)); next }
  okc <- hi > lo
  if (!okc) FAIL <- FAIL + 1L
  cat(sprintf("  [%s] %-18s %s (%.4f) > %s (%.4f)   %s\n",
              if (okc) "PASS" else "FAIL", cc$signature, cc$high_bin, hi, cc$low_bin, lo, cc$why))
}
unchecked <- setdiff(sc_cols, unique(CELLSTATE_CONTRASTS$signature))
cat(sprintf("\n  NOT checked by any contrast (no settled compartment expectation): %s\n",
            paste(unchecked, collapse = ", ")))
cat("  These are scored but UNVALIDATED. A green suite above says nothing about them.\n")

cat("\n================ mean score by bin (all signatures) ================\n")
M <- A[!is.na(hierarchy_bin) & nzchar(hierarchy_bin),
       lapply(.SD, function(v) round(mean(v, na.rm = TRUE), 3)), by = hierarchy_bin, .SDcols = sc_cols]
print(M)

cat("\n================ the venetoclax axis: MCL1 vs BCL2 by bin ================\n")
if (all(c("expr_MCL1", "expr_BCL2") %in% names(A))) {
  V <- A[!is.na(hierarchy_bin) & nzchar(hierarchy_bin),
         .(n = .N, MCL1 = round(mean(expr_MCL1, na.rm = TRUE), 3),
           BCL2 = round(mean(expr_BCL2, na.rm = TRUE), 3),
           BCL2L1 = round(mean(expr_BCL2L1, na.rm = TRUE), 3)), by = hierarchy_bin]
  V[, MCL1_minus_BCL2 := round(MCL1 - BCL2, 3)]
  print(V[order(-MCL1_minus_BCL2)])
}
if ("malignant" %in% names(A) && A[!is.na(malignant), .N] > 0) {
  cat("\n  ...and malignant vs normal:\n")
  print(A[!is.na(malignant), .(n = .N,
         priming = round(mean(apoptotic_priming, na.rm = TRUE), 4),
         antiapop = round(mean(antiapoptotic, na.rm = TRUE), 4),
         oxphos = round(mean(oxphos, na.rm = TRUE), 4),
         fao = round(mean(fatty_acid_oxidation, na.rm = TRUE), 4)), by = malignant])
}

fwrite(M, file.path(HIER_TAB_DIR, "cellstate_by_bin.csv"))
if (!is.null(cov_report)) fwrite(cov_report, file.path(HIER_TAB_DIR, "cellstate_gene_coverage.csv"))
cat(sprintf("\n[done] %s\n", file.path(HIER_TAB_DIR, "cellstate_by_bin.csv")))
if (FAIL > 0) { cat(sprintf("\n[!] %d positive control(s) FAILED -- do not use these scores yet.\n", FAIL))
                quit(save = "no", status = 1) }
cat("positive controls PASS\n")
