# 13_tnk_composition_diagnostic.R ----
# INPUT  : HIER_PROJ_DIR/<ds>/<sample>__bmm_percell.csv        (bmm_broad, mapping_error, high_error)
#          DIR_MALIGNANCY/<ds>/<sample>__consensus_percell.csv (per-cell malignant call)
#          HIER_PROJ_DIR/<ds>/<sample>__stemness_percell.csv   (LSC17)
#          CYTOTRACE_DIR/<ds>/<sample>__cytotrace2_percell.csv (CytoTRACE2_Score, absolute scale)
#          DIR_CCC/ccc_sample_manifest.csv                     (healthy label, CCC eligibility)
# OUTPUT : DIR_TABLES/03_hierarchy/tnk_composition_diagnostic.csv
# WHAT   : Splits the T_NK node by its BoneMarrowMap sub-label and asks which sub-label carries the
#          AML-vs-healthy stemness signal, and whether that sub-label is well projected.
#
# WHY. 05_ccc/04's per-bin sweep left T_NK as the only bin with anything (beta=+0.021, p=0.081), and
# purification STRENGTHENED it (+0.043 at q=50), so it is not inferCNV under-detection. But T_NK is also
# the bin already known to hold misassigned blasts (malignant "T_NK" CytoTRACE2 potency 0.135 vs normal
# 0.096), and misassignment is a PROJECTION error, not a CNV error -- a blast with low CNV burden that
# landed in T_NK is invisible to a burden-based purification. This script separates the two.
#
# WHY NOT SPLIT THE GRAPH INSTEAD. The admission screen (2026-08-21) rejected all 22 candidate node
# splits: at fine granularity, node PRESENCE tracks the disease label (Pre-B +57pp toward healthy,
# Early GMP -50pp toward AML), so a finer vocabulary would let FGW separate on node-set rather than
# topology. This diagnostic answers the same question without touching CCC_NODES.
#
# READING IT:
#   signal in CD4/CD8/NK, mapping_error comparable to the rest -> real T/NK biology
#   signal in Early Lymphoid (bin map flags it AMBIGUOUS: CLP-like, upstream of BOTH T and B),
#     or concentrated in high_error cells                      -> misassignment, not a finding
#
# Usage  : Rscript scripts/03_hierarchy/13_tnk_composition_diagnostic.R [--n_perm 10000]
suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--n_perm",   type = "integer", default = 10000L),
  make_option("--min_cells", type = "integer", default = 20L)
)))
set.seed(SEED)
out_csv <- file.path(DIR_HIERARCHY, "tnk_composition_diagnostic.csv")
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

man <- fread(file.path(DIR_CCC, "ccc_sample_manifest.csv"))[ccc_eligible == TRUE]
man[, healthy := Timepoint == "Healthy"]
message("[1] eligible samples: ", nrow(man), " (healthy ", sum(man$healthy), " / AML ", sum(!man$healthy), ")")

## -- Step 1. per-cell join, T_NK only ----
pull_one <- function(ds, smp) {
  fb <- file.path(HIER_PROJ_DIR,  ds, paste0(smp, "__bmm_percell.csv"))
  fc <- file.path(DIR_MALIGNANCY, ds, paste0(smp, "__consensus_percell.csv"))
  fs <- file.path(HIER_PROJ_DIR,  ds, paste0(smp, "__stemness_percell.csv"))
  fy <- file.path(CYTOTRACE_DIR,  ds, paste0(smp, "__cytotrace2_percell.csv"))
  if (!all(file.exists(fb, fc, fs))) return(NULL)
  b <- fread(fb, select = c("cell", "bmm_broad", "mapping_error", "high_error", "hierarchy_bin"))
  b <- b[hierarchy_bin == "T_NK"]                       # the node under investigation
  if (!nrow(b)) return(NULL)
  d <- merge(b, fread(fc, select = c("cell", "malignant")), by = "cell")
  d <- merge(d, fread(fs, select = c("cell", "LSC17")),     by = "cell")
  # CytoTRACE2 is OPTIONAL: it is the independent replication arm, and a sample without it must lose
  # only that column, not its rows. An inner join here would silently shrink the cohort.
  if (file.exists(fy)) {
    d <- merge(d, fread(fy, select = c("cell", "CytoTRACE2_Score")), by = "cell", all.x = TRUE)
  } else d[, CytoTRACE2_Score := NA_real_]
  d[, `:=`(dataset = ds, sample = smp)][]
}
cells <- rbindlist(lapply(seq_len(nrow(man)), function(i) pull_one(man$dataset[i], man$sample[i])), fill = TRUE)
stopifnot(nrow(cells) > 0)
message("[2] T_NK cells: ", nrow(cells), " over ", uniqueN(cells[, .(dataset, sample)]), " samples")
message("    with CytoTRACE2: ", sum(!is.na(cells$CytoTRACE2_Score)),
        " (", round(100 * mean(!is.na(cells$CytoTRACE2_Score))), "%)")

H <- unique(man[, .(dataset, sample, healthy)])
cells <- merge(cells, H, by = c("dataset", "sample"))

# THE FINDING IS ABOUT NON-MALIGNANT CELLS. Restrict to them, exactly as mean_stemness_normal does.
norm <- cells[malignant == 0]
message("[2] non-malignant T_NK cells: ", nrow(norm), " of ", nrow(cells),
        " (", round(100 * nrow(norm) / nrow(cells)), "%)")

## -- Step 2. composition and projection quality per sub-label ----
message("\n[3] T_NK composition, and how well each sub-label projects")
comp <- norm[, .(cells = .N,
                 pct_of_TNK = NA_real_,
                 mapping_error = mean(mapping_error, na.rm = TRUE),
                 pct_high_error = 100 * mean(high_error == 1 | high_error == TRUE, na.rm = TRUE),
                 LSC17 = mean(LSC17, na.rm = TRUE),
                 cyto = mean(CytoTRACE2_Score, na.rm = TRUE)), by = bmm_broad]
comp[, pct_of_TNK := 100 * cells / sum(cells)]
setorder(comp, -cells)
print(comp[, .(bmm_broad, cells, pct_of_TNK = round(pct_of_TNK, 1),
               mapping_error = round(mapping_error, 2), pct_high_error = round(pct_high_error, 1),
               LSC17 = round(LSC17, 4), cyto = round(cyto, 4))])

## -- Step 3. which sub-label carries the AML-vs-healthy signal ----
# Same model as 05_ccc/04: dataset fixed effects + blast_proxy, label permuted WITHIN dataset.
BLAST_BINS <- c("HSC_MPP", "LMPP_GMP", "Mono_DC")
bp_src <- rbindlist(lapply(seq_len(nrow(man)), function(i) {
  f <- file.path(HIER_PROJ_DIR, man$dataset[i], paste0(man$sample[i], "__bmm_percell.csv"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, select = "hierarchy_bin")
  data.table(dataset = man$dataset[i], sample = man$sample[i],
             blast_proxy = sum(d$hierarchy_bin %in% BLAST_BINS) / nrow(d))
}))
ds_both <- H[, .(nh = sum(healthy), na = sum(!healthy)), by = dataset][nh > 0 & na > 0]$dataset
message("\n[4] datasets with both labels: ", paste(ds_both, collapse = ", "))

# PER-TEST PERMUTATION STREAMS, KEYED ON A STABLE LABEL. Without this a single set.seed(SEED) at the
# top of the script serves every fit_perm() call in the process, so a test's p-value depends on how
# many tests ran before it: reversing the order here moved 35 of 36 per-bin p-values while beta stayed
# bit-identical. No verdict flipped (shifts are MC-noise scale, max |dp| 0.045 at n_perm=2000), but
# none of the published p-values were reproducible without re-running the identical --levels list in
# the identical order.
#
# WHY THIS HASH AND NOT A SIMPLER ONE. The obvious base-R template,
#   SEED + sum(utf8ToInt(s) * seq_along(utf8ToInt(s)))
# COLLIDES on this script's real label space -- "bin|q=17|T_NK" and "bin|q=90|T_NK" both map to
# 499398, 21 collisions over the integer-q labels and 5037 if fractional --levels are passed. A
# collision hands two DIFFERENT tests the SAME permutations, which is worse than the bug being fixed
# and would invalidate any BH correction later applied across the per-bin family. The polynomial
# rolling hash below is collision-free on every label space these scripts generate (verified), and
# modulo 2^31-1 keeps it inside set.seed's 32-bit signed range without overflow: h*131+ch stays under
# 2.8e11, exact in a double.
.perm_seed <- function(label) {
  M <- 2147483647
  h <- 0
  for (ch in utf8ToInt(label)) h <- (h * 131 + ch) %% M
  as.integer((h + SEED) %% M)
}

fit_perm <- function(D, n_perm, label) {
  set.seed(.perm_seed(label))
  if (uniqueN(D$dataset) < 2 || uniqueN(D$healthy) < 2 || nrow(D) < 10) return(c(beta = NA, p = NA))
  X <- model.matrix(~ factor(dataset) + blast_proxy, data = D)
  Q <- qr.Q(qr(X)); rs <- function(v) v - Q %*% (t(Q) %*% v)
  ry <- rs(D$y); ra <- rs(as.numeric(!D$healthy))
  den <- sum(ra * ra); if (den <= 0) return(c(beta = NA, p = NA))
  beta <- sum(ry * ra) / den
  idx <- split(seq_len(nrow(D)), D$dataset); ex <- 0L
  for (k in seq_len(n_perm)) {
    a <- as.numeric(!D$healthy); for (g in idx) a[g] <- sample(a[g])
    rp <- rs(a); dd <- sum(rp * rp); if (dd <= 0) next
    ex <- ex + (abs(sum(ry * rp) / dd) >= abs(beta))
  }
  c(beta = beta, p = (1 + ex) / (n_perm + 1))
}

test_one <- function(sub, measure) {
  S <- norm[bmm_broad == sub, .(y = mean(get(measure), na.rm = TRUE), n = .N),
            by = .(dataset, sample, healthy)]
  S <- S[n >= opt$min_cells & is.finite(y)]
  S <- merge(S, bp_src, by = c("dataset", "sample"))
  D <- S[dataset %in% ds_both]
  ok <- D[, .(k = uniqueN(healthy)), by = dataset][k == 2]$dataset   # a sub-label can lose a label
  D <- D[dataset %in% ok]
  fp <- fit_perm(D, opt$n_perm, paste(sub, measure, sep = "|"))
  data.table(bmm_broad = sub, measure = measure, n_samples = nrow(D),
             n_healthy = sum(D$healthy), n_aml = sum(!D$healthy),   # was sum(!D$healthy == TRUE & D$healthy),
             # which R parses as !(D$healthy == TRUE) & D$healthy -- identically FALSE, so it read 0 in all 15 rows

             mean_healthy = D[healthy == TRUE, mean(y)], mean_aml = D[healthy == FALSE, mean(y)],
             beta = fp["beta"], p = fp["p"])
}

subs <- comp[cells >= 100]$bmm_broad     # a sub-label seen <100x cohort-wide cannot support a test
message("[4] testing sub-labels: ", paste(subs, collapse = ", "))
res <- rbindlist(c(lapply(subs, test_one, measure = "LSC17"),
                   lapply(subs, test_one, measure = "CytoTRACE2_Score"),
                   lapply(subs, test_one, measure = "mapping_error")), fill = TRUE)

message("\n[5] AML-vs-healthy within-dataset effect, per T_NK sub-label")
for (m in c("LSC17", "CytoTRACE2_Score", "mapping_error")) {
  message("  -- ", m, " --")
  r <- res[measure == m][order(p)]
  for (i in seq_len(nrow(r)))
    message(sprintf("     %-18s n=%-3d  healthy=%+.4f  AML=%+.4f  beta=%+.5f  p=%.4f",
                    r$bmm_broad[i], r$n_samples[i], r$mean_healthy[i], r$mean_aml[i], r$beta[i], r$p[i]))
}

## -- Step 4. verdict ----
L <- res[measure == "LSC17"][order(p)]
top <- L$bmm_broad[1]
AMBIG <- c("Early Lymphoid")                       # flagged AMBIGUOUS in bmm_bin_map.tsv
me_top <- comp[bmm_broad == top, mapping_error]
me_rest <- comp[bmm_broad != top, weighted.mean(mapping_error, cells)]
message("\n[6] VERDICT")
message("    strongest LSC17 sub-label: ", top, "  (beta=", signif(L$beta[1], 3), ", p=", signif(L$p[1], 4), ")")
message("    its mapping_error ", round(me_top, 2), " vs ", round(me_rest, 2), " for the rest of T_NK")
if (top %in% AMBIG) {
  message("    -> the carrier is an AMBIGUOUS label. This is a projection artefact, not T/NK biology.")
} else if (is.finite(me_top) && is.finite(me_rest) && me_top > 1.5 * me_rest) {
  message("    -> the carrier projects markedly WORSE than the rest of T_NK: treat as misassignment.")
} else if (L$p[1] < 0.05) {
  message("    -> the signal sits in a well-defined, well-projected lymphoid subset. Real T/NK biology.")
} else {
  message("    -> no sub-label reaches significance. The T_NK effect does not localise; report it as a")
  message("       node-level observation only, and do not name a cell type.")
}
fwrite_safe(res, out_csv)
fwrite_safe(comp, file.path(DIR_HIERARCHY, "tnk_composition_summary.csv"))
message("\n[done] wrote ", out_csv)
