# 04_aggregate_scan.R ----
# Roll the 220 per-sample scan outputs up into the cohort tables the analysis actually runs on,
# and answer the three gate questions Batch 1 was built to answer:
#   G1 NECTIN4/PVRL4 detection  -> is the PI's Nectin-4 question answerable at all?
#   G2 macrophage cell counts   -> should Macrophage become its own CCC node, or stay a feature?
#   G3 fibrosis/ECM detection   -> which panel genes are measurable in which compartment?
#
# INPUT  : LCC_TAB_DIR/03_detect/*.csv, 03_myeloid/*.csv, 03_sample_manifest.csv
#          LCC_PERCELL_DIR/<ds>/<sample>__lcc_percell.csv.gz
# OUTPUT : LCC_TAB_DIR/04_detection_by_sample.csv    gene x sample x stratum (the full grid)
#          LCC_TAB_DIR/04_detection_cohort.csv       gene x stratum cohort summary
#          LCC_TAB_DIR/04_nectin_gate.csv            G1: Nectin family, per dataset
#          LCC_TAB_DIR/04_myeloid_gate.csv           G2: macrophage counts + the node verdict
#          LCC_TAB_DIR/04_pathway_sample_bin.csv     sample x bin x pathway mean UCell (analysis input)
#          LCC_TAB_DIR/04_pathway_sample.csv         sample-level, all cells + malignant-only
# Usage  : Rscript LCC_proj/scripts/04_aggregate_scan.R
#
# SYMBOL DRIFT IS REAL AND HANDLED: the cohort's 13 datasets were annotated against different
# GENCODE releases. Chen2023 carries PVRL4, the GRCh38-2020-A reference carries NECTIN4, and a
# sample only ever has ONE of them. G1 therefore reports detection under the UNION of both symbols
# per sample, plus which symbol each sample actually uses -- collapsing them would read as dropout.

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("LCC_proj", "scripts", "config_lcc.R"))

## -- Step 1. detection tables ----
message("[1] rolling up per-sample detection")
det <- rbindlist(lapply(list.files(file.path(LCC_TAB_DIR, "03_detect"), full.names = TRUE), fread),
                 fill = TRUE, use.names = TRUE)
man <- load_sample_meta(include_stroma_ref = TRUE)   # timepoint-CORRECTED; see config_lcc.R
det[man, `:=`(timepoint = i.timepoint, uid_patient = i.uid_patient), on = c("dataset", "sample")]
det[, gene_present := !is.na(n_nonzero)]
fwrite_safe(det, file.path(LCC_TAB_DIR, "04_detection_by_sample.csv"))

panel <- load_gene_panel()
cat_of <- rbind(panel[, .(gene, category)], panel[nzchar(alias), .(gene = alias, category)])
coh <- det[gene_present == TRUE & n_cells >= 20,
           .(n_samples = .N,
             median_pct_nonzero = as.numeric(median(pct_nonzero)),
             p90_pct_nonzero    = as.numeric(quantile(pct_nonzero, 0.9)),
             frac_samples_gt1pct = mean(pct_nonzero > 1),
             mean_lognorm       = mean(mean_lognorm)),
           by = .(gene, stratum, level)]
coh[cat_of, category := i.category, on = "gene"]
fwrite_safe(coh[order(level, stratum, -median_pct_nonzero)], file.path(LCC_TAB_DIR, "04_detection_cohort.csv"))

## -- Step 2. G1 -- Nectin gate, resolving the PVRL4/NECTIN4 symbol drift ----
message("[2] G1: Nectin family detection")
nec <- det[level == "sample" & gene %in% LCC_NECTIN_GENES]
# per sample, which symbol for the target gene is in the reference at all
tgt <- nec[gene %in% c("NECTIN4", "PVRL4")]
g1 <- tgt[, .(symbol_in_ref = paste(sort(gene[gene_present]), collapse = "/"),
              n_cells       = max(n_cells),
              n_nonzero     = sum(n_nonzero, na.rm = TRUE),
              pct_nonzero   = sum(n_nonzero, na.rm = TRUE) / max(n_cells) * 100),
          by = .(dataset, sample, timepoint)]
g1[symbol_in_ref == "", symbol_in_ref := "ABSENT_FROM_REFERENCE"]
setorder(g1, -pct_nonzero)
fwrite_safe(g1, file.path(LCC_TAB_DIR, "04_nectin_gate.csv"))

fam <- nec[gene_present == TRUE, .(n_samples = .N, median_pct = as.numeric(median(pct_nonzero)),
                                   max_pct = max(pct_nonzero)), by = gene][order(-median_pct)]
message("    Nectin family cohort-wide (median % non-zero cells per sample):")
print(fam)
message("    NECTIN4/PVRL4: ", sum(g1$pct_nonzero > 1), "/", nrow(g1),
        " samples exceed 1% detection; max = ", sprintf("%.2f%%", max(g1$pct_nonzero)))

## -- Step 3. G2 -- macrophage node gate ----
message("[3] G2: macrophage counts")
my <- rbindlist(lapply(list.files(file.path(LCC_TAB_DIR, "03_myeloid"), full.names = TRUE), fread), fill = TRUE)
my[man, `:=`(timepoint = i.timepoint, uid_patient = i.uid_patient), on = c("dataset", "sample")]
my[, mac_node_usable := n_macrophage_like >= LCC_MAC_MIN_CELLS]
fwrite_safe(my[order(-n_macrophage_like)], file.path(LCC_TAB_DIR, "04_myeloid_gate.csv"))
my[man, in_selection := i.in_selection, on = c("dataset", "sample")]
n_dg_ok <- my[timepoint == "Diagnosis" & mac_node_usable == TRUE, .N]
n_dg_ok_sel <- my[timepoint == "Diagnosis" & mac_node_usable == TRUE & in_selection == TRUE, .N]
message(sprintf("    Diagnosis samples with >= %d macrophage-like cells: %d cohort-wide, %d in the selected datasets (gate needs >= %d)",
                LCC_MAC_MIN_CELLS, n_dg_ok, n_dg_ok_sel, LCC_MAC_GATE_MIN_SAMPLES))
message(sprintf("    VERDICT: %s",
        if (n_dg_ok >= LCC_MAC_GATE_MIN_SAMPLES) "re-run CellChat with Macrophage as its own node"
        else "keep the locked 7-bin vocabulary; macrophage stays a node FEATURE"))

## -- Step 4. pathway scores aggregated to the analysis units ----
message("[4] aggregating UCell scores to sample x bin")
uc_cols <- NULL
agg <- rbindlist(lapply(seq_len(nrow(man)), function(i) {
  f <- file.path(LCC_PERCELL_DIR, man$dataset[i], paste0(man$sample[i], "__lcc_percell.csv.gz"))
  if (!file.exists(f)) return(NULL)
  p <- fread(f)
  if (is.null(uc_cols)) uc_cols <<- grep("_UCell$", names(p), value = TRUE)
  bybin <- p[, c(list(n_cells = .N), lapply(.SD, mean)), by = hierarchy_bin, .SDcols = uc_cols]
  bybin[, `:=`(dataset = man$dataset[i], sample = man$sample[i], unit = "bin")]
  setnames(bybin, "hierarchy_bin", "group")
  bybin
}), fill = TRUE)
agg[man, `:=`(timepoint = i.timepoint, uid_patient = i.uid_patient, study = i.study),
    on = c("dataset", "sample")]
fwrite_safe(agg, file.path(LCC_TAB_DIR, "04_pathway_sample_bin.csv"))

message("[5] aggregating UCell scores to sample level (all / malignant / nonmalignant)")
aggs <- rbindlist(lapply(seq_len(nrow(man)), function(i) {
  f <- file.path(LCC_PERCELL_DIR, man$dataset[i], paste0(man$sample[i], "__lcc_percell.csv.gz"))
  if (!file.exists(f)) return(NULL)
  p <- fread(f)
  uc <- grep("_UCell$", names(p), value = TRUE)
  parts <- list(all = p, malignant = p[malignant == 1L], nonmalignant = p[malignant == 0L])
  rbindlist(lapply(names(parts), function(k) {
    d <- parts[[k]]
    if (!nrow(d)) return(NULL)
    cbind(data.table(dataset = man$dataset[i], sample = man$sample[i], unit = "sample",
                     group = k, n_cells = nrow(d)),
          d[, lapply(.SD, mean), .SDcols = uc])
  }))
}), fill = TRUE)
aggs[man, `:=`(timepoint = i.timepoint, uid_patient = i.uid_patient, study = i.study),
     on = c("dataset", "sample")]
fwrite_safe(aggs, file.path(LCC_TAB_DIR, "04_pathway_sample.csv"))
message("[done] 04_aggregate_scan")
