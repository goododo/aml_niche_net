#!/usr/bin/env Rscript
# 10_progeny_score.R ----
# PROGENy pathway activity per cell (14 pathways), via decoupleR's multivariate linear model.
#
# WHY THIS IS THE HIGHEST-VALUE ADDITION FOR A CCC PAPER, and it is not "one more feature":
#
#   CellChat infers "A signals to B" from the CO-EXPRESSION of a ligand in A and its receptor in B.
#   That is an expression pattern, not evidence that anything was transmitted. PROGENy scores the
#   DOWNSTREAM TRANSCRIPTIONAL FOOTPRINT of a pathway in the receiver -- the consequence, not the
#   machinery. An inferred edge whose receiver shows the matching pathway footprint is a different
#   class of claim from one that does not.
#
# So this is used two ways: as a node feature, and as an independent check on the CCC edges that
# carry the paper's central argument.
#
# WHY PROGENy AND NOT A PATHWAY GENE SET. Membership sets (KEGG/Hallmark) score the pathway's own
# components, which are frequently NOT transcriptionally regulated by the pathway. PROGENy weights
# are fitted from perturbation experiments on genes that actually MOVE when the pathway is
# manipulated, which is what makes it a footprint rather than an inventory.
#
# THE MODEL IS BUNDLED with the progeny package (5,782 genes x 14 pathways) -- no OmniPath call, so
# this runs on a compute node with no network and gives the same answer every time.
#
# VALIDATION IS CROSS-METHOD, not self-consistency. PROGENy's JAK-STAT footprint and the
# interferon_I / interferon_II UCell panels from 09 are built from different genes by different
# methods; interferon signals through JAK-STAT, so they must correlate. If they do not, one of the
# two is broken and neither should be used. That is a far stronger test than checking PROGENy
# against itself.
#
# OUTPUT: HIER_PROJ_DIR/<dataset>/<sample>__progeny_percell.csv
#           cell, sample, dataset, hierarchy_bin, malignant, <14 pathway activities>
#
#   Rscript scripts/03_hierarchy/10_progeny_score.R [--dataset X] [--limit N] [--force]

suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(progeny); library(decoupleR); library(Matrix)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = ""),
  make_option("--limit",   type = "integer",   default = 0L),
  make_option("--force",   action = "store_true", default = FALSE)
)))

## -- the model, as a decoupleR net ----
M   <- progeny::getModel("Human", top = PROGENY_TOP)
net <- as.data.table(M, keep.rownames = "target")
net <- melt(net, id.vars = "target", variable.name = "source", value.name = "weight")[weight != 0]
net[, source := as.character(source)]
message(sprintf("[0] PROGENy: %d pathways, %d gene-weight pairs (top=%d per pathway)",
                uniqueN(net$source), nrow(net), PROGENY_TOP))
PATHWAYS <- sort(unique(net$source))

R <- qc_rds_roster(on_extra = "ignore")
if (nzchar(opt$dataset)) R <- R[dataset %in% strsplit(opt$dataset, ",")[[1]]]
if (opt$limit > 0) R <- head(R, opt$limit)
stopifnot(nrow(R) > 0)

dst_of  <- function(ds, sid) file.path(CELLSTATE_DIR, ds, paste0(sid, "__progeny_percell.csv"))
proj_of <- function(ds, sid) file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
cons_of <- function(ds, sid) file.path(DIR_MALIGNANCY, ds, paste0(sid, "__consensus_percell.csv"))

score_one <- function(rds, ds, sid) {
  seu <- readRDS(rds); seu <- NormalizeData(seu, verbose = FALSE)
  dat <- SeuratObject::LayerData(seu, assay = "RNA", layer = "data")
  # decoupleR needs the intersection non-trivial; a near-empty overlap means a symbol convention
  # mismatch and would still return numbers for every cell.
  ov <- intersect(rownames(dat), unique(net$target))
  if (length(ov) < PROGENY_MIN_OVERLAP)
    stop(sprintf("only %d genes overlap the PROGENy model (need >= %d) -- gene symbol mismatch?",
                 length(ov), PROGENY_MIN_OVERLAP))
  # SUBSET FIRST, AND STAY SPARSE. Handing decoupleR the full as.matrix() object densified 33,694
  # genes x N cells when only ~5,400 rows are in the model at all; measured 26s vs >5min per sample
  # for an identical result (decoupleR intersects against the net internally either way, so the
  # dropped rows could never have contributed). At 214 samples that is the difference between
  # ~2.4 h and not finishing.
  acts <- decoupleR::run_mlm(mat = dat[ov, , drop = FALSE], net = net, .source = "source",
                             .target = "target", .mor = "weight", minsize = PROGENY_MINSIZE)
  A <- dcast(as.data.table(acts)[statistic == "mlm"], condition ~ source, value.var = "score")
  setnames(A, "condition", "cell")
  for (p in setdiff(PATHWAYS, names(A))) A[, (p) := NA_real_]

  pf <- proj_of(ds, sid); cf <- cons_of(ds, sid)
  if (file.exists(pf)) A <- merge(A, fread(pf, select = c("cell", "hierarchy_bin")), by = "cell", all.x = TRUE)
  if (file.exists(cf)) A <- merge(A, fread(cf, select = c("cell", "malignant")), by = "cell", all.x = TRUE)
  else A[, malignant := NA_integer_]
  A[, `:=`(dataset = ds, sample = sid)]
  rm(seu, dat, acts); gc(verbose = FALSE)
  A[, c("cell", "dataset", "sample", "hierarchy_bin", "malignant", PATHWAYS), with = FALSE]
}

t0 <- Sys.time()
n_fail <- 0L; failed_samples <- character(0)
for (i in seq_len(nrow(R))) {
  ds <- R$dataset[i]; sid <- R$sample[i]; dst <- dst_of(ds, sid)
  # FRESHNESS, not existence. 33 of these files (GSE147989 4, GSE185991 29; 78,315 cells) were
  # written BEFORE their consensus existed and carry malignant = NA permanently; the validation
  # below filters !is.na(malignant) and so silently reports on a subset of the cohort.
  .ins <- c(R$rds[i], proj_of(ds, sid), cons_of(ds, sid))
  if (!is_stale(dst, .ins, force = opt$force)) {
    message(sprintf("[%d/%d] %s::%s current", i, nrow(R), ds, sid)); next
  }
  if (file.exists(dst))
    message(sprintf("[%d/%d] %s::%s RECOMPUTE -- %s", i, nrow(R), ds, sid, stale_reason(dst, .ins, force = opt$force)))
  message(sprintf("[%d/%d] %s::%s", i, nrow(R), ds, sid))
  x <- tryCatch(score_one(R$rds[i], ds, sid),
                error = function(e) { message("  [FAIL] ", conditionMessage(e)); NULL })
  # see 09_cellstate_score.R: a failed sample keeps its stale file, so validating afterwards mixes
  # fresh and stale outputs and reports PASS on a rebuild that did not happen
  if (is.null(x)) { n_fail <- n_fail + 1L; failed_samples <- c(failed_samples, paste0(ds, "::", sid)); next }
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  fwrite_safe(x, dst)
  if (i == 1) message(sprintf("  [timing] first sample took %.1f min",
                              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ---------------------------------------------------------------------------------------------
## validation
## ---------------------------------------------------------------------------------------------
if (n_fail > 0) {
  cat(sprintf("\n[!] %d of %d samples FAILED and kept their previous output: %s\n",
              n_fail, nrow(R), paste(utils::head(failed_samples, 10), collapse = ", ")))
  cat("    Refusing to validate: the tables below would mix fresh and stale files, which is how a\n")
  cat("    broken rebuild reports PASS. Failures here are usually transient -- re-run and retry.\n")
  quit(save = "no", status = 1)
}
done <- R[file.exists(dst_of(dataset, sample))]
if (!nrow(done)) { message("[!] nothing scored"); quit(save = "no", status = 1) }
P <- rbindlist(lapply(seq_len(nrow(done)), function(i)
  fread(dst_of(done$dataset[i], done$sample[i]))), fill = TRUE)
cat(sprintf("\n================ %d cells across %d samples ================\n", nrow(P), uniqueN(P$sample)))

FAIL <- 0L
cat("\n================ CROSS-METHOD CONTROL: PROGENy JAK-STAT vs UCell interferon ================\n")
cs_files <- file.path(CELLSTATE_DIR, done$dataset, paste0(done$sample, "__cellstate_percell.csv"))
cs_files <- cs_files[file.exists(cs_files)]
if (!length(cs_files)) {
  cat("  [SKIP] 09_cellstate_score.R has not run on these samples -- the check cannot be made.\n")
  cat("  Treat PROGENy as UNVALIDATED until it can.\n")
} else {
  CS <- rbindlist(lapply(cs_files, fread, select = c("cell", "sample", "interferon_I", "interferon_II")), fill = TRUE)
  J <- merge(P[, .(cell, sample, `JAK-STAT`, TNFa, NFkB)], CS, by = c("cell", "sample"))
  if (nrow(J) < 1000) {
    cat(sprintf("  [SKIP] only %d cells joined -- not enough to test\n", nrow(J)))
  } else {
    r1 <- cor(J$`JAK-STAT`, J$interferon_I,  use = "complete.obs", method = "spearman")
    r2 <- cor(J$`JAK-STAT`, J$interferon_II, use = "complete.obs", method = "spearman")
    # a NEGATIVE control from the same data: JAK-STAT should track interferon more than an
    # unrelated pathway does, otherwise the correlation is just shared library-size structure
    r0 <- cor(P$`JAK-STAT`, P$Androgen, use = "complete.obs", method = "spearman")
    cat(sprintf("  cells joined: %d\n", nrow(J)))
    cat(sprintf("  spearman(JAK-STAT, interferon_I)  = %+.3f\n", r1))
    cat(sprintf("  spearman(JAK-STAT, interferon_II) = %+.3f\n", r2))
    cat(sprintf("  spearman(JAK-STAT, Androgen)      = %+.3f   <- unrelated-pathway baseline\n", r0))
    okc <- r1 > PROGENY_MIN_IFN_COR && r2 > PROGENY_MIN_IFN_COR && max(r1, r2) > abs(r0)
    if (!okc) FAIL <- FAIL + 1L
    cat(sprintf("  [%s] interferon signals through JAK-STAT: both correlations > %.2f and above the baseline\n",
                if (okc) "PASS" else "FAIL", PROGENY_MIN_IFN_COR))
  }
}

cat("\n================ pathway activity by bin ================\n")
B <- P[!is.na(hierarchy_bin) & nzchar(hierarchy_bin),
       c(lapply(.SD, function(v) round(mean(v, na.rm = TRUE), 3)), .(n = .N)),
       by = hierarchy_bin, .SDcols = PATHWAYS]
print(B)

cat("\n================ malignant vs normal ================\n")
if (P[!is.na(malignant), .N] > 0)
  print(P[!is.na(malignant), c(lapply(.SD, function(v) round(mean(v, na.rm = TRUE), 3)), .(n = .N)),
          by = malignant, .SDcols = PATHWAYS])

fwrite(B, file.path(HIER_TAB_DIR, "progeny_by_bin.csv"))
cat(sprintf("\n[done] %s\n", file.path(HIER_TAB_DIR, "progeny_by_bin.csv")))
if (FAIL > 0) { cat("\n[!] the cross-method control FAILED -- PROGENy and the interferon panel disagree.\n")
                cat("    One of the two is wrong; do not use either until it is resolved.\n")
                quit(save = "no", status = 1) }
cat("cross-method control PASS\n")
