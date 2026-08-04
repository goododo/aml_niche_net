# 01_parse_infercnv_regions.R ----
# Turn per-sample inferCNV HMM output into arm-level CNV event calls per malignant subclone.
# This is the technical work the HANDOFF flagged as unverified: "these two sources have NOT been
# validated for the specific query 'extract 17p status'".
#
# INPUT  : INFERCNV_ROOT/<ds>/<sample>/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat
#          INFERCNV_ROOT/<ds>/<sample>/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.genes_used.dat
#          INFERCNV_ROOT/<ds>/<sample>/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings
#          DIR_MALIGNANCY/<ds>/<sample>__consensus_percell.csv     (main-line labels, READ-ONLY)
#          LCC_PANEL_DIR/chrom_arms_grch38.tsv
# OUTPUT : LCC_TAB_DIR/01_percnv/<ds>__<sample>__subclone.csv      (per-subclone summary)
#          LCC_TAB_DIR/01_percnv/<ds>__<sample>__armevents.csv     (per-subclone x arm, non-zero only)
#          LCC_TAB_DIR/01_subclone_cnv.csv                         (rollup, all samples)
#          LCC_TAB_DIR/01_arm_events_long.csv                      (rollup, all samples)
#          LCC_TAB_DIR/01_parse_qc.csv                             (per-sample parse QC)
# Usage  : Rscript LCC_proj/scripts/01_parse_infercnv_regions.R [--dataset=Petti2019] [--force]
#
# WHY gene-fraction and not bp-fraction for arm coverage: inferCNV only observes the genome where
# expressed genes are. A bp-based arm coverage would systematically under-call gene-poor arms and
# is not comparable across samples with different gene universes (6.3k-8k genes used here).
#
# WHY reference_normal subclones are kept and not filtered: they are inferCNV's own internal
# normal cells, so any arm event called in them is a per-sample FALSE-POSITIVE readout. This is the
# built-in negative control for the whole proxy and 02 uses it to flag unreliable samples.
#
# NOTE [main-line known issue 2]: malignancy labels are from 07-14 while the QC objects were
# refreshed 07-15, so inferCNV ran on a LARGER cell set than the current consensus covers.
# Verified here on the cohort: cons_covered_frac == 1.000 for every sample -- every consensus-
# labelled cell IS present in the inferCNV run (consensus is a strict subset). The consequence is
# therefore NOT missing labels but a reduced effective n per subclone: in 25/130 samples only
# 13-78% of inferCNV cells carry a label (icnv_labelled_frac). Downstream must gate subclones on
# n_cells_labelled, NOT n_cells, and re-test excluding low-icnv_labelled_frac samples.

suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("LCC_proj", "scripts", "config_lcc.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = NULL,
              help = "restrict to one dataset dir under INFERCNV_ROOT (default: all)"),
  make_option("--force", action = "store_true", default = FALSE,
              help = "recompute samples whose per-sample output already exists")
)))

PERCNV_DIR <- file.path(LCC_TAB_DIR, "01_percnv")
dir.create(PERCNV_DIR, recursive = TRUE, showWarnings = FALSE)

## -- Step 1. chromosome arm scaffold ----
message("[1] loading GRCh38 arm boundaries")
ARMS <- load_chrom_arms()
setkey(ARMS, chrom)

# gene -> arm ("chr17" + midpoint 7.67Mb -> "17p"). Genes whose midpoint falls inside the
# centromere gap are assigned to the nearer side; that band holds essentially no expressed genes.
assign_arm <- function(chr, start, end) {
  mid <- (as.numeric(start) + as.numeric(end)) / 2
  a   <- ARMS[.(chr)]
  ifelse(is.na(a$cen_start), NA_character_,
         paste0(sub("^chr", "", chr),
                ifelse(mid < a$cen_start, "p",
                       ifelse(mid > a$cen_end, "q",
                              ifelse(mid < (a$cen_start + a$cen_end) / 2, "p", "q")))))
}

## -- Step 2. enumerate samples ----
ds_dirs <- list.dirs(INFERCNV_ROOT, recursive = FALSE)
if (!is.null(opt$dataset)) ds_dirs <- ds_dirs[basename(ds_dirs) == opt$dataset]
if (!length(ds_dirs)) stop("[01] no dataset dirs matched under ", INFERCNV_ROOT)

jobs <- rbindlist(lapply(ds_dirs, function(d) {
  data.table(dataset = basename(d), sample = basename(list.dirs(d, recursive = FALSE)), dir = list.dirs(d, recursive = FALSE))
}))
message("[2] ", nrow(jobs), " sample dirs found")

## -- Step 3. per-sample parse ----
parse_one <- function(dataset, sample, dir) {
  f_genes <- file.path(dir, LCC_ICNV_GENES_FILE)
  f_used  <- file.path(dir, LCC_ICNV_USED_FILE)
  f_group <- file.path(dir, LCC_ICNV_GROUP_FILE)
  f_cons  <- file.path(DIR_MALIGNANCY, dataset, paste0(sample, "__consensus_percell.csv"))
  miss    <- c(f_genes, f_used, f_group)[!file.exists(c(f_genes, f_used, f_group))]
  if (length(miss)) return(list(qc = data.table(dataset, sample, status = paste0("missing:", basename(miss[1])))))

  ## 3a. gene universe of THIS sample -> arm denominators.
  # genes_used.dat has a 3-name header ("chr start stop") over 4 data columns (gene is unnamed).
  used <- fread(f_used, header = FALSE, skip = 1, col.names = c("gene", "chr", "start", "end"))
  used[, arm := assign_arm(chr, start, end)]
  used <- used[!is.na(arm) & !arm %in% LCC_ARMS_EXCLUDE]
  arm_denom <- used[, .(n_genes_arm = .N), by = arm]

  ## 3b. per-subclone gene-level HMM states (already sparse: only predicted CNV regions).
  gp <- fread(f_genes, select = c("cell_group_name", "state", "gene"))
  gp[, state := as.integer(state)]
  gp <- gp[state != LCC_STATE_NEUTRAL]
  gp[used, arm := i.arm, on = "gene"]
  gp <- gp[!is.na(arm)]
  gp[, direction := fifelse(state %in% LCC_STATE_LOSS, "loss",
                     fifelse(state %in% LCC_STATE_GAIN, "gain", NA_character_))]
  gp <- gp[!is.na(direction)]

  ## 3c. subclone sizes + malignant fraction from the main-line consensus labels.
  grp <- fread(f_group, col.names = c("cell_group_name", "cell"))
  n_cells_grouped <- nrow(grp)
  # Two DIFFERENT quantities, both needed (see header note on known issue 2):
  #   icnv_labelled_frac = share of inferCNV cells that carry a consensus label  -> drives effective n
  #   cons_covered_frac  = share of consensus cells present in the inferCNV run  -> must be 1.0,
  #                        anything less means the two label sets genuinely disagree.
  icnv_labelled_frac <- NA_real_; cons_covered_frac <- NA_real_; cons_n <- NA_integer_
  if (file.exists(f_cons)) {
    cons <- fread(f_cons, select = c("cell", "malignant"))
    cons_n <- nrow(cons)
    cons_covered_frac <- mean(cons$cell %in% grp$cell)
    grp[cons, malignant := i.malignant, on = "cell"]
    icnv_labelled_frac <- mean(!is.na(grp$malignant))
  } else {
    grp[, malignant := NA_integer_]
  }
  sub <- grp[, .(n_cells = .N,
                 n_cells_labelled = sum(!is.na(malignant)),
                 n_malignant = sum(malignant == 1L, na.rm = TRUE)), by = cell_group_name]
  sub[, frac_malignant := fifelse(n_cells_labelled > 0, n_malignant / n_cells_labelled, NA_real_)]
  sub[, is_reference := grepl("^reference_normal", cell_group_name)]
  # subclone share of the sample's malignant pool (drives LCC_CLONE_FRAC_MIN downstream)
  tot_mal <- sum(sub[is_reference == FALSE]$n_malignant)
  sub[, clone_frac_of_malignant := if (tot_mal > 0) n_malignant / tot_mal else NA_real_]

  ## 3d. arm coverage per subclone x arm x direction.
  arm_hits <- gp[, .(n_genes_nonneutral = .N), by = .(cell_group_name, arm, direction)]
  arm_hits[arm_denom, n_genes_arm := i.n_genes_arm, on = "arm"]
  arm_hits[, frac_arm := n_genes_nonneutral / n_genes_arm]
  arm_hits[, is_event := frac_arm >= LCC_ARM_FRAC_MIN]
  arm_hits[, `:=`(dataset = dataset, sample = sample)]

  ## 3e. per-subclone rollup. An arm counts once even if both directions clear the floor
  ## (a whole-arm loss plus a focal gain is still one aberrant arm).
  ev <- arm_hits[is_event == TRUE]
  by_sub <- ev[, .(n_arm_events      = uniqueN(arm),
                   n_arm_events_loss = uniqueN(arm[direction == "loss"]),
                   n_arm_events_gain = uniqueN(arm[direction == "gain"]),
                   arms_hit = paste(sort(unique(paste0(arm, ":", substr(direction, 1, 1)))), collapse = ";")),
               by = cell_group_name]
  sub[by_sub, `:=`(n_arm_events = i.n_arm_events, n_arm_events_loss = i.n_arm_events_loss,
                   n_arm_events_gain = i.n_arm_events_gain, arms_hit = i.arms_hit),
      on = "cell_group_name"]
  sub[is.na(n_arm_events), `:=`(n_arm_events = 0L, n_arm_events_loss = 0L, n_arm_events_gain = 0L, arms_hit = "")]

  ## 3f. TP53 locus readout: the state of the TP53 gene itself in each subclone.
  tp53_in_universe <- LCC_TP53_GENE %in% used$gene
  tp53 <- gp[gene == LCC_TP53_GENE, .(tp53_state = state[1], tp53_direction = direction[1]), by = cell_group_name]
  sub[tp53, `:=`(tp53_state = i.tp53_state, tp53_direction = i.tp53_direction), on = "cell_group_name"]
  # TP53 absent from a subclone's predicted-CNV rows means neutral there, not missing.
  sub[is.na(tp53_state) & tp53_in_universe, `:=`(tp53_state = LCC_STATE_NEUTRAL, tp53_direction = "neutral")]
  sub[, tp53_loss := !is.na(tp53_state) & tp53_state %in% LCC_STATE_LOSS]
  # whole-17p support for the TP53 call (a single-gene state is noisy on its own)
  a17p <- arm_hits[arm == LCC_TP53_ARM & direction == "loss", .(arm17p_frac_loss = frac_arm), by = cell_group_name]
  sub[a17p, arm17p_frac_loss := i.arm17p_frac_loss, on = "cell_group_name"]
  sub[is.na(arm17p_frac_loss), arm17p_frac_loss := 0]

  sub[, `:=`(dataset = dataset, sample = sample)]
  setnames(sub, "cell_group_name", "subclone")
  setnames(arm_hits, "cell_group_name", "subclone")

  qc <- data.table(dataset, sample, status = "ok",
                   n_genes_used = nrow(used), tp53_in_universe = tp53_in_universe,
                   n_subclone_obs = sum(!sub$is_reference), n_subclone_ref = sum(sub$is_reference),
                   n_cells_grouped = n_cells_grouped, n_cells_consensus = cons_n,
                   icnv_labelled_frac = icnv_labelled_frac, cons_covered_frac = cons_covered_frac,
                   frac_genes_nonneutral = nrow(gp) / max(1, nrow(used) * uniqueN(gp$cell_group_name)))
  list(sub = sub, arms = arm_hits, qc = qc)
}

message("[3] parsing (resume-safe; per-sample outputs under ", PERCNV_DIR, ")")
qc_all <- vector("list", nrow(jobs))
for (i in seq_len(nrow(jobs))) {
  ds <- jobs$dataset[i]; sm <- jobs$sample[i]
  f_sub <- file.path(PERCNV_DIR, sprintf("%s__%s__subclone.csv", ds, sm))
  f_arm <- file.path(PERCNV_DIR, sprintf("%s__%s__armevents.csv", ds, sm))
  f_qc  <- file.path(PERCNV_DIR, sprintf("%s__%s__qc.csv", ds, sm))
  if (!opt$force && file.exists(f_sub) && file.exists(f_qc)) {
    qc_all[[i]] <- fread(f_qc); next                                        # resume-skip
  }
  res <- parse_one(ds, sm, jobs$dir[i])
  if (is.null(res$sub)) { message("    [skip] ", ds, "/", sm, " -- ", res$qc$status); qc_all[[i]] <- res$qc; next }
  fwrite_safe(res$sub, f_sub)
  fwrite_safe(res$arms[n_genes_nonneutral > 0], f_arm)
  fwrite_safe(res$qc, f_qc)
  qc_all[[i]] <- res$qc
  message(sprintf("    [write] %-14s %-28s subclones=%3d(ref %3d)  labelled=%.2f  covered=%.2f  TP53gene=%s",
                  ds, sm, res$qc$n_subclone_obs, res$qc$n_subclone_ref,
                  res$qc$icnv_labelled_frac, res$qc$cons_covered_frac, res$qc$tp53_in_universe))
}

## -- Step 4. rollup ----
message("[4] rolling up")
roll <- function(pat) rbindlist(lapply(list.files(PERCNV_DIR, pattern = pat, full.names = TRUE), fread),
                                fill = TRUE, use.names = TRUE)
fwrite_safe(roll("__subclone\\.csv$"),   file.path(LCC_TAB_DIR, "01_subclone_cnv.csv"))
fwrite_safe(roll("__armevents\\.csv$"),  file.path(LCC_TAB_DIR, "01_arm_events_long.csv"))
fwrite_safe(rbindlist(qc_all, fill = TRUE), file.path(LCC_TAB_DIR, "01_parse_qc.csv"))
message("[done] 01_parse_infercnv_regions")
