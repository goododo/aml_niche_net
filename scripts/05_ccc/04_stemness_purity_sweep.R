# 04_stemness_purity_sweep.R ----
# INPUT  : LARGE1/05_cnv_snv/infercnv_burden/<ds>/<sample>_infercnv_burden.csv  (per-cell CNV burden)
#          DIR_MALIGNANCY/<ds>/<sample>__consensus_percell.csv                 (per-cell malignant call)
#          CCC_BMM_DIR/<ds>/<sample>__stemness_percell.csv                     (per-cell LSC17 + bin)
#          DIR_CCC/ccc_sample_manifest.csv                                     (healthy label, eligibility)
# OUTPUT : DIR_CCC/stemness_purity_sweep.csv  (one row per purification level x model)
# WHAT   : Asks whether "AML has higher stemness among NON-malignant cells" survives progressive
#          purification of the non-malignant pool, or decays away as an inferCNV under-detection artifact.
#
# WHY THIS EXISTS. 08_scoring/07 found the H2 feature signal is carried by mean_stemness_normal
# (within-dataset p=0.0070) and not by mean_stemness_malignant (p=0.139). "Normal" there means
# "consensus did not call this cell malignant" -- which is not the same as "this cell is not a blast".
# In the 8 samples where an ALLELE-based caller (numbat) also ran, 24.1% of the cells inferCNV called
# normal were called malignant by numbat, and those missed cells carry roughly double the LSC17
# (0.0697 vs 0.0350). 130 of the 138 CCC-eligible samples -- including ALL 37 healthy ones -- have
# inferCNV only. Healthy marrow has few true blasts to miss, so the bias is ONE-DIRECTIONAL: it
# inflates the AML arm alone, which is exactly the shape of the finding.
#
# THE SWEEP. At level q, drop the top q% of the currently-normal cells BY CNV BURDEN, per sample, then
# recompute the per-sample mean LSC17 over what remains. Contaminating blasts sit at the high-burden end
# (that is what the caller is thresholding on), so rising q buys purity at the cost of sample size.
#   effect stable / growing across q  -> the signal is not the missed blasts
#   effect decaying monotonically     -> the signal IS the missed blasts
#
# THE PROCEDURE IS SYMMETRIC. Healthy samples are purified identically. If burden and stemness are
# coupled for ordinary biological reasons, that coupling is removed from BOTH arms, so the contrast
# stays fair. An asymmetric version would manufacture the very decay it claims to measure.
#
# Usage  : Rscript scripts/05_ccc/04_stemness_purity_sweep.R [--levels 0,10,20,30,40,50] [--n_perm 10000]
suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--levels", type = "character", default = "0,10,20,30,40,50"),
  make_option("--n_perm", type = "integer",   default = 10000L),
  make_option("--force",  action = "store_true", default = FALSE)
)))
LEVELS  <- as.numeric(strsplit(opt$levels, ",")[[1]])
BURDEN  <- file.path(LARGE1_DIR, "05_cnv_snv", "infercnv_burden")
out_csv <- file.path(DIR_CCC, "stemness_purity_sweep.csv")
set.seed(SEED)

man <- fread(file.path(DIR_CCC, "ccc_sample_manifest.csv"))[ccc_eligible == TRUE]
man[, healthy := Timepoint == "Healthy"]
message("[1] eligible samples: ", nrow(man), "  (healthy ", sum(man$healthy), " / AML ", sum(!man$healthy), ")")

## -- Step 1. per-cell join: burden + consensus call + LSC17 ----
# Every cell must carry all three or it is dropped, and the drop is COUNTED. A silent inner join here
# would let a sample with no burden file quietly leave the cohort at one level and not another, which
# would show up as a dose-response that is really a changing denominator.
pull_one <- function(ds, smp) {
  fb <- file.path(BURDEN, ds, paste0(smp, "_infercnv_burden.csv"))
  fc <- file.path(DIR_MALIGNANCY, ds, paste0(smp, "__consensus_percell.csv"))
  fs <- file.path(CCC_BMM_DIR, ds, paste0(smp, "__stemness_percell.csv"))
  if (!all(file.exists(fb, fc, fs))) return(NULL)
  b <- fread(fb, select = c("cell", "infercnv_burden", "group"))
  cn <- fread(fc, select = c("cell", "malignant"))
  st <- fread(fs, select = c("cell", "LSC17", "hierarchy_bin"))
  d <- merge(merge(b, cn, by = "cell"), st, by = "cell")
  d <- d[!is.na(malignant) & !is.na(LSC17) & !is.na(infercnv_burden)]
  if (!nrow(d)) return(NULL)
  d[, `:=`(dataset = ds, sample = smp)]
  d[]
}
cells <- rbindlist(lapply(seq_len(nrow(man)), function(i) pull_one(man$dataset[i], man$sample[i])), fill = TRUE)
got <- unique(cells[, .(dataset, sample)])
message("[2] cells joined: ", nrow(cells), " over ", nrow(got), " / ", nrow(man), " samples")
lost <- fsetdiff(man[, .(dataset, sample)], got)
if (nrow(lost)) message("    samples with no complete per-cell join (dropped): ", nrow(lost), " -> ",
                        paste(head(lost[, paste(dataset, sample)], 5), collapse = ", "))
stopifnot(nrow(cells) > 0, nrow(got) > 10)

# Reference cells belong to the inferCNV run, not to the sample. Leaving them in would mix a second
# donor's marrow into the "normal" pool of every sample that used an external reference.
n_ref <- sum(cells$group != "observation")
cells <- cells[group == "observation"]
message("[2] dropped ", n_ref, " reference cells; ", nrow(cells), " observation cells remain")

## -- Step 2. the sweep ----
H <- unique(man[, .(dataset, sample, healthy)])
cells <- merge(cells, H, by = c("dataset", "sample"))

# Within-dataset test, matching every other H2 result: only datasets holding BOTH labels contribute,
# so platform cannot produce the effect. This is the number 08_scoring/07 reports as p_strat.
ds_both <- H[, .(nh = sum(healthy), na = sum(!healthy)), by = dataset][nh > 0 & na > 0]$dataset
message("[3] datasets with both labels: ", paste(ds_both, collapse = ", "))

# blast_proxy from the FULL cell set (not the purified one) -- it is a property of the sample, and
# recomputing it per level would make the covariate move with the treatment.
BLAST_BINS <- c("HSC_MPP", "LMPP_GMP", "Mono_DC")
BP <- cells[, .(blast_proxy = sum(hierarchy_bin %in% BLAST_BINS) / .N), by = .(dataset, sample)]

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
  # is_aml effect on per-sample mean LSC17, dataset absorbed as a fixed effect (FWL), permuting the
  # label WITHIN dataset so the null keeps each dataset's composition intact.
  # blast_proxy is in the model because every H2 result controls for it. Without it this would be a
  # DIFFERENT estimand from the one being defended, and a decay (or its absence) would not transfer.
  X <- model.matrix(~ factor(dataset) + blast_proxy, data = D)
  Q <- qr.Q(qr(X)); res <- function(v) v - Q %*% (t(Q) %*% v)
  ry <- res(D$stem); ra <- res(as.numeric(!D$healthy))
  den <- sum(ra * ra); if (den <= 0) return(c(beta = NA_real_, p = NA_real_))
  beta <- sum(ry * ra) / den
  idx <- split(seq_len(nrow(D)), D$dataset)
  ex <- 0L
  for (k in seq_len(n_perm)) {
    a <- as.numeric(!D$healthy)
    for (g in idx) a[g] <- sample(a[g])
    rp <- res(a); dd <- sum(rp * rp); if (dd <= 0) next
    ex <- ex + (abs(sum(ry * rp) / dd) >= abs(beta))
  }
  c(beta = beta, p = (1 + ex) / (n_perm + 1))
}

res <- rbindlist(lapply(LEVELS, function(q) {
  # purify: within each sample, drop the top q% of NORMAL cells by burden
  keep <- cells[malignant == 0]
  keep[, cut := if (q <= 0) Inf else quantile(infercnv_burden, 1 - q / 100, na.rm = TRUE), by = .(dataset, sample)]
  kept <- keep[infercnv_burden <= cut]
  S <- kept[, .(stem = mean(LSC17), n_cells = .N), by = .(dataset, sample, healthy)]
  S <- S[n_cells >= 20]                                  # a mean over <20 cells is noise, not a measurement
  S <- merge(S, BP, by = c("dataset", "sample"))
  D <- S[dataset %in% ds_both]
  fp <- fit_perm(D, opt$n_perm, sprintf("pooled|q=%s", q))
  data.table(level_pct = q,
             n_samples = nrow(S), n_within = nrow(D),
             n_cells_kept = nrow(kept), frac_cells_kept = nrow(kept) / nrow(keep),
             mean_healthy = S[healthy == TRUE,  mean(stem)],
             mean_aml     = S[healthy == FALSE, mean(stem)],
             diff_raw     = S[healthy == FALSE, mean(stem)] - S[healthy == TRUE, mean(stem)],
             beta_strat   = fp["beta"], p_strat = fp["p"])
}))

message("\n[4] PURITY SWEEP (within-dataset model, ", opt$n_perm, " permutations)")
print(res[, .(level_pct, n_within, frac_cells_kept = round(frac_cells_kept, 3),
              mean_healthy = round(mean_healthy, 5), mean_aml = round(mean_aml, 5),
              diff_raw = round(diff_raw, 5), beta_strat = round(beta_strat, 5),
              p_strat = round(p_strat, 5))])

## -- Step 3. read the shape, not a single row ----
b <- res$beta_strat
message("\n[5] VERDICT")
message("    beta at q=0   : ", signif(b[1], 4), "   p=", signif(res$p_strat[1], 4))
message("    beta at q=", max(LEVELS), "  : ", signif(b[length(b)], 4), "   p=", signif(res$p_strat[length(b)], 4))
retained <- b[length(b)] / b[1]
# READ BETA, NOT P. Purification throws away cells, so p necessarily worsens with q whatever the truth
# is; conditioning the verdict on the terminal p (as a first draft of this script did) would report
# "partial contamination" for a trajectory that is dead flat. The dose-response question is entirely
# about whether the EFFECT SIZE decays.
if (!is.finite(retained)) {
  message("    -> beta undefined at one end; inspect the table above.")
} else if (retained >= 0.7) {
  message("    -> SURVIVES purification: ", round(100 * retained), "% of the effect retained.")
  message("       Removing the highest-burden ", max(LEVELS), "% of the 'normal' pool does not shrink the")
  message("       AML-vs-healthy gap, so the gap is not the blasts inferCNV missed.")
  if (res$p_strat[length(b)] >= 0.05)
    message("       (Terminal p=", signif(res$p_strat[length(b)], 3), " reflects the discarded sample size, not decay.)")
} else if (retained <= 0.3) {
  message("    -> DECAYS to ", round(100 * retained), "% of the effect. Consistent with the signal being")
  message("       undetected blasts in the 'normal' pool. It cannot be reported as a microenvironment result.")
} else {
  message("    -> PARTIAL: ", round(100 * retained), "% retained. Some of the effect is contamination and some")
  message("       is not; report the purified effect size, not the q=0 one.")
}
# Power drops as q rises, so a rising p ALONE is not evidence of decay -- beta is the quantity to read.
message("    n_within at q=0 -> q=", max(LEVELS), ": ", res$n_within[1], " -> ", res$n_within[length(b)],
        " | cells retained: ", round(100 * res$frac_cells_kept[length(b)]), "%")

fwrite_safe(res, out_csv)

## -- Step 4. the same sweep PER BIN ----
# The pooled mean above is the wrong estimand for what 08_scoring/07 actually tests. That script feeds
# mean_stemness_normal as a per-NODE feature, so it sees the 7-bin PROFILE; the pooled mean collapses it.
# On the 3 both-label datasets the per-bin differences point BOTH ways (LMPP_GMP/T_NK/Mono_DC/Mk up,
# HSC_MPP/Erythroid/B_Plasma down) and 66% of the total movement cancels on pooling -- which is why the
# pooled test sits at p~0.14 while the profile test sits at p=0.007. Sweeping per bin is therefore the
# apples-to-apples contamination check; the pooled one above answers a strictly weaker question.
res_bin <- rbindlist(lapply(LEVELS, function(q) {
  keep <- cells[malignant == 0]
  keep[, cut := if (q <= 0) Inf else quantile(infercnv_burden, 1 - q / 100, na.rm = TRUE), by = .(dataset, sample)]
  kept <- keep[infercnv_burden <= cut]
  SB <- kept[, .(stem = mean(LSC17), n_cells = .N), by = .(dataset, sample, healthy, hierarchy_bin)]
  SB <- SB[n_cells >= 20]
  SB <- merge(SB, BP, by = c("dataset", "sample"))
  rbindlist(lapply(CCC_NODES, function(b) {
    D <- SB[hierarchy_bin == b & dataset %in% ds_both]
    if (uniqueN(D$healthy) < 2 || nrow(D) < 10) return(NULL)
    ok <- D[, .(n = uniqueN(healthy)), by = dataset][n == 2]$dataset   # a bin can lose a label at high q
    D <- D[dataset %in% ok]
    if (uniqueN(D$dataset) < 2 || uniqueN(D$healthy) < 2) return(NULL)
    fp <- fit_perm(D, opt$n_perm, sprintf("bin|q=%s|%s", q, b))
    data.table(level_pct = q, hierarchy_bin = b, n = nrow(D),
               beta_strat = fp["beta"], p_strat = fp["p"])
  }))
}))

message("\n[6] PER-BIN SWEEP -- beta by purification level (the estimand 08_scoring/07 actually uses)")
W <- dcast(res_bin, hierarchy_bin ~ level_pct, value.var = "beta_strat")
print(W)
message("\n[6] retention of beta from q=0 to q=", max(LEVELS), " per bin:")
for (b in W$hierarchy_bin) {
  v <- as.numeric(W[hierarchy_bin == b, -1, with = FALSE]); v <- v[is.finite(v)]
  if (length(v) < 2 || v[1] == 0) next
  message(sprintf("    %-15s %+.5f -> %+.5f   (%d%% retained)  p(q=0)=%.4f",
                  b, v[1], v[length(v)], round(100 * v[length(v)] / v[1]),
                  res_bin[hierarchy_bin == b & level_pct == 0, p_strat][1]))
}
out_bin <- file.path(DIR_CCC, "stemness_purity_sweep_bybin.csv")
fwrite_safe(res_bin, out_bin)
message("\n[done] wrote ", out_csv, "\n[done] wrote ", out_bin)
