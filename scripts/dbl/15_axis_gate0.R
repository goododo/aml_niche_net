#!/usr/bin/env Rscript
# 15_axis_gate0.R ----
# input : DBL_LANE_DIR/<lane>.rds + <lane>__labels.csv.gz for the two smoke lanes (class-A singlets only; no real doublet is touched)
# output: DBL_TAB_DIR/dbl_axis_gate0_smoke.csv  per lane x partner x thinning: n, Spearman rho(parent axis, doublet axis), median shift
#         DBL_TAB_DIR/dbl_axis_gate0_singlets_smoke.csv  axis spread over leukemic singlets by state (does the ruler separate states?)
# what it does: Gate 0 of the redesigned doublet axis (post-hoc probe, 2026-11-25, not pre-registered). Axis = log ratio of raw
#   UMI in a stem program vs a myeloid program, both restricted to genes leukemic cells express >= 5x more than T cells (pooled CPM LEUK > 5).
#   Adding a T (or NK/B) cell then adds ~0 to both sums, and a ratio ignores how much of the leukemic half was captured, so the
#   axis should survive mixing. Synthetic doublets = raw-count sums of a class-A partner and a class-A leukemic singlet (as stage
#   04), then binomially thinned (keep 1 / 0.5 / 0.3) to mimic capture loss. Bar (same as stage 07): rho >= 0.7.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
STEM <- c("CD34","SPINK2","PROM1","HOPX","CRHBP","MEIS1","HLF","AVP","MLLT3","NPR3","BAALC","MMRN1","CD200","ANGPT1","EGFL7","GATA2","FAM30A","MYCN")
MYEL <- c("MPO","ELANE","AZU1","PRTN3","CTSG","LYZ","S100A8","S100A9","CST3","FCN1","CD14","MS4A6A","CEBPE","CSF3R","FCER1G","CLEC12A","MNDA","SERPINA1")
N_PAIR <- 1000L; KEEP <- c(1, 0.5, 0.3); set.seed(20261125)
axis_of <- function(M, gs, gm) log((Matrix::colSums(M[gs, , drop = FALSE]) + 0.5) / (Matrix::colSums(M[gm, , drop = FALSE]) + 0.5))
thin <- function(M, p) { if (p == 1) return(M); M@x <- as.numeric(rbinom(length(M@x), as.integer(M@x), p)); drop0(M) }  # binomial thinning of each UMI count

RES <- list(); SING <- list()
for (L in DBL_SMOKE_LANES) {
  x <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds"))); M <- x$counts
  lab <- fread(file.path(DBL_LANE_DIR, paste0(L, "__labels.csv.gz")))[match(colnames(M), cell)][class == "A"]
  leuk <- lab[compartment %in% DBL_LEUK_STATES, cell]
  cpm <- function(cc) { s <- M[, cc, drop = FALSE]; Matrix::rowSums(s) / sum(s) * 1e6 }
  cT <- cpm(lab[compartment == "T", cell]); cL <- cpm(leuk)
  # relative filter >= 5x (cT < 1 and then >= 10x left < 3 genes in the pooled lane: its T cells carry ambient myeloid RNA)
  ok <- names(cT)[cL > 5 & cL >= 5 * cT]; gs <- intersect(STEM, ok); gm <- intersect(MYEL, ok)
  message(sprintf("[%s] T-silent stem genes %d: %s | myeloid genes %d: %s", L, length(gs), paste(gs, collapse=","), length(gm), paste(gm, collapse=",")))
  stopifnot(length(gs) >= 3, length(gm) >= 3)                                  # sanity: both programs must exist after the T filter
  ax <- axis_of(M[, leuk], gs, gm)
  SING[[L]] <- lab[cell %in% leuk, .(n = .N, axis_median = median(ax[cell]), axis_iqr = IQR(ax[cell])), by = compartment][, lane := L]
  st_rank <- match(lab[match(leuk, cell), compartment], DBL_LEUK_STATES)
  SING[[L]][, rho_axis_vs_state_order := cor(ax, st_rank, method = "spearman")]   # stem -> myeloid order should give rho < 0
  for (pt in intersect(c("T", "NK", "B"), unique(lab$compartment))) {
    pa <- lab[compartment == pt, cell]; a <- sample(pa, N_PAIR, replace = length(pa) < N_PAIR); b <- sample(leuk, N_PAIR, replace = TRUE)
    D <- M[, a] + M[, b]                                                          # synthetic doublet = raw-count sum (stage 04)
    for (p in KEEP) { ad <- axis_of(thin(D, p), gs, gm); ap <- ax[b]
      RES[[length(RES) + 1]] <- data.table(lane = L, partner = pt, keep = p, n = N_PAIR, n_partner_cells = length(pa),
        rho = round(cor(ap, ad, method = "spearman"), 3), median_shift = round(median(ad - ap), 3), iqr_parent = round(IQR(ap), 3)) } }
  rm(x, M); gc()
}
RES <- rbindlist(RES); SING <- rbindlist(SING)
fwrite(RES, file.path(DBL_TAB_DIR, "dbl_axis_gate0_smoke.csv")); fwrite(SING, file.path(DBL_TAB_DIR, "dbl_axis_gate0_singlets_smoke.csv"))
print(RES); print(SING); message("verdict (T partner, keep 0.5): ", paste(RES[partner == "T" & keep == 0.5, sprintf("%s rho %.2f %s", lane, rho, ifelse(rho >= 0.7, "PASS", "FAIL"))], collapse = " | "))
