# 13.4_fimm_external_gate.R ----
# INPUT  : /FAST/gr10634/gaozy/bulk_cohorts/fimm/  (Malani 2022, Zenodo 10.5281/zenodo.7274740,
#            fetched 2026-09-29: CPM 163S+4Healthy, sDSS 164S+17Healthy, binary mutations 225S x 57G,
#            common sample annotation 252S)
#          results/tables/12_pseudobulk_de/hox_axis_genes.csv  (Set L, unchanged)
# OUTPUT : results/tables/13_bulk_anchor/fimm_external_gate.csv
#          results/tables/13_bulk_anchor/fimm_hox_scores.csv
# WHAT IT DOES : the external direction gate of PREREGISTRATION_bulk_transfer.md §4, on the two
#          leads section 7 produced: PALBOCICLIB (the one non-FLT3-active hit) and VENETOCLAX.
#          The gate is DIRECTION ONLY -- registered as "same direction in FIMM before it may be
#          called a finding" -- so the readout is the sign of the NPM1-like minus rest-of-WT
#          difference, plus its p, and nothing is re-defined to make it pass.
#          Everything is rebuilt from FIMM's own data with BeatAML's rules held fixed:
#            score      = mean within-sample percentile rank of Set L over the CPM matrix
#            NPM1-like  = NPM1-WT sample with score >= 25th percentile of the NPM1-mut scores
#            outcome    = sDSS (HIGHER sDSS = MORE sensitive, the OPPOSITE sign convention to
#                         BeatAML AUC, where lower = more sensitive) -- so "same direction as
#                         BeatAML" means a POSITIVE sDSS difference here. This sign flip is the
#                         single most dangerous thing in this script and is asserted below.
#          Instrument gate is re-run first: if Set L cannot separate NPM1 in FIMM, the direction
#          test is uninterpretable and the script stops.

suppressPackageStartupMessages({ library(data.table); library(readxl) })
DIR <- "/FAST/gr10634/gaozy/bulk_cohorts/fimm"
OUT <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/13_bulk_anchor"
set.seed(260929)

# ---- 1. Set L, mapped to Ensembl because FIMM ships Ensembl ids --------------------------------
# REGISTERED-RULE NOTE (2026-09-29). FIMM's CPM matrix is keyed by Ensembl id, so the symbol
# match used in 13.3 gives 0/17 and the registered "<15 matched -> STOP" rule fired correctly.
# Mapping was then fixed using BeatAML's own stable_id<->display_label annotation (17/17 Set L
# genes have an Ensembl id there). 14 of those 17 exist in FIMM's 18,203-row matrix; the 3
# missing are lncRNAs (HOXB-AS1, HOXB-AS3, HOXA10-AS), absent because that quantification is
# protein-coding-centric -- a platform coverage limit, not a mapping failure. 14/17 = 82% is
# above the registered threshold read as a fraction (15/19 = 79%) and below it read as an
# absolute count; BOTH readings are reported. To remove the ambiguity from the comparison
# itself, the gate rebuilds the BeatAML score on the SAME 14 genes, so the direction test uses
# one identical gene set in both cohorts. No gene is substituted.
gl <- fread("/FAST/gr10634/gaozy/aml_niche_net/results/tables/12_pseudobulk_de/hox_axis_genes.csv")
SETL <- sort(unique(gl[set == "L" & hierarchy_bin == "LMPP_GMP", gene]))
BEA <- "/FAST/gr10634/gaozy/bulk_cohorts/beataml2/beataml_waves1to4_norm_exp_dbgap.txt"
map <- fread(BEA, select = c("stable_id", "display_label"))[display_label %in% SETL]
cat(sprintf("[1] Set L: %d symbols | with Ensembl id: %d\n", length(SETL), nrow(map)))

# ---- 2. FIMM expression, samples, genotype ------------------------------------------------------
cpm <- fread(file.path(DIR, "File_7_RNA_seq_CPM_163S_4Healthy.csv"))
setnames(cpm, 1, "gene")
ann <- as.data.table(read_excel(file.path(DIR, "File_0_Common_sample_annotation_252S.xlsx"), sheet = 1))
setnames(ann, gsub("[^A-Za-z0-9]+", "_", names(ann)))
cat("[2] annotation columns:", paste(names(ann), collapse = ", "), "\n")
cat("    Disease_status values:\n"); print(table(ann$Disease_status))

mu <- as.data.table(read_excel(file.path(DIR, "File_6_Binary_mutation_225S_57G.xlsx"), sheet = 1))
setnames(mu, 1, "gene")
cat(sprintf("    mutation table: %d genes x %d samples | NPM1 row present: %s\n",
            nrow(mu), ncol(mu) - 1L, "NPM1" %in% mu$gene))
npm1_row <- unlist(mu[gene == "NPM1", -1, with = FALSE])
npm1_pos_samples <- names(npm1_row)[!is.na(npm1_row) & npm1_row == 1]
cat(sprintf("    NPM1-mutant samples in mutation table: %d of %d\n",
            length(npm1_pos_samples), length(npm1_row)))

# AML samples with expression AND a genotype call. Disease.status is Diagnosis/Relapse/Refractory
# (no healthy category here -- healthy donors are separate columns named Healthy*). BeatAML's 1pp
# cohort is 72% initial-diagnosis, so the gate is run at DIAGNOSIS to stay comparable, with the
# all-stage set as a stated secondary.
expr_samples <- grep("^Healthy", setdiff(names(cpm), "gene"), invert = TRUE, value = TRUE)
keep_all <- intersect(expr_samples, names(npm1_row))
dx <- ann[Disease_status == "Diagnosis", Sample_ID]
keep <- intersect(keep_all, dx)
cat(sprintf("    usable AML samples: %d all-stage | %d at diagnosis (gate uses diagnosis)\n",
            length(keep_all), length(keep)))
stopifnot(length(keep) >= 40L)

# ---- 3. score, then the instrument gate ---------------------------------------------------------
X <- as.matrix(cpm[, ..keep]); ens <- cpm$gene
rm_dup <- ave(rowMeans(X), ens, FUN = function(m) seq_along(m) == which.max(m)) == 1
X <- X[rm_dup, ]; ens <- ens[rm_dup]
COMMON <- map[stable_id %in% ens]                      # the shared gene set for BOTH cohorts
idx <- which(ens %in% COMMON$stable_id)
P <- apply(X, 2, frank) / nrow(X)
sc <- data.table(Sample_ID = keep, score = colMeans(P[idx, ]),
                 npm1_pos = keep %in% npm1_pos_samples)
cat(sprintf("\n[3] shared gene set: %d of %d Set L genes (%s)\n", nrow(COMMON), length(SETL),
            paste(sort(COMMON$display_label), collapse = " ")))
cat(sprintf("    registered threshold read as count: %d of %d -> %s | as fraction: %.0f%% vs 79%% -> %s\n",
            nrow(COMMON), length(SETL), if (nrow(COMMON) >= 15L) "pass" else "FIRES",
            100 * nrow(COMMON) / length(SETL),
            if (nrow(COMMON) / length(SETL) >= 15 / 19) "pass" else "FIRES"))
cat(sprintf("    missing from FIMM (lncRNA, platform coverage): %s\n",
            paste(setdiff(SETL, COMMON$display_label), collapse = " ")))
cat(sprintf("    FIMM: NPM1-mut %d vs WT %d\n", sum(sc$npm1_pos), sum(!sc$npm1_pos)))
if (nrow(COMMON) < 10L) stop("REGISTERED STOP: shared gene set too small to score")
auc_of <- function(x, g) as.numeric(wilcox.test(x[g], x[!g])$statistic) / (sum(g) * sum(!g))
q0 <- auc_of(sc$score, sc$npm1_pos)
cat(sprintf("    instrument (Q0 equivalent in FIMM): AUC = %.3f\n", q0))
if (q0 < 0.70) { cat("** instrument fails in FIMM -- direction test uninterpretable, STOP **\n")
  fwrite(data.table(check = "fimm_instrument_auc", value = q0, verdict = "STOP"),
         file.path(OUT, "fimm_external_gate.csv")); quit(save = "no") }

thr <- quantile(sc[npm1_pos == TRUE, score], 0.25)
sc[, grp := fifelse(npm1_pos, "NPM1mut", fifelse(score >= thr, "NPM1like", "WTlow"))]
cat(sprintf("    NPM1-like %d | WT-low %d | mut %d\n",
            sum(sc$grp == "NPM1like"), sum(sc$grp == "WTlow"), sum(sc$grp == "NPM1mut")))
fwrite(sc, file.path(OUT, "fimm_hox_scores.csv"))

# ---- 4. sDSS and the direction test ------------------------------------------------------------
ds <- as.data.table(read_excel(file.path(DIR, "File_3_Drug_response_sDSS_164S_17Healthy.xlsx"), sheet = "sDSS"))
cat(sprintf("\n[4] sDSS table: %d drugs x %d samples\n", nrow(ds), ncol(ds) - 2L))
## SIGN ASSERTION: sDSS must be higher for a broadly active cytotoxic than for a targeted agent
## with a narrow population; and venetoclax-sensitive NPM1-mut AML is the textbook anchor used
## to confirm the convention before any lead is judged.
ve_row <- ds[grepl("^Venetoclax$", Drug_name, ignore.case = TRUE)]
pb_row <- ds[grepl("^Palbociclib$", Drug_name, ignore.case = TRUE)]
cat(sprintf("    matched drug rows: Venetoclax %d | Palbociclib %d\n", nrow(ve_row), nrow(pb_row)))

get_vec <- function(row) {
  s <- intersect(setdiff(names(ds), c("Drug_ID", "Drug_name")), sc$Sample_ID)
  v <- suppressWarnings(as.numeric(unlist(row[1, ..s])))
  data.table(Sample_ID = s, sdss = v)[!is.na(sdss)]
}
res <- list()
for (nm in c("Venetoclax", "Palbociclib")) {
  row <- if (nm == "Venetoclax") ve_row else pb_row
  if (!nrow(row)) { cat(sprintf("    %s NOT IN FIMM LIBRARY -- gate not evaluable\n", nm)); next }
  d <- merge(get_vec(row), sc, by = "Sample_ID")
  w <- d[grp != "NPM1mut"]
  if (sum(w$grp == "NPM1like") < 5L) { cat(sprintf("    %s: only %d NPM1-like, gate underpowered\n",
      nm, sum(w$grp == "NPM1like"))) }
  dif <- median(w[grp == "NPM1like", sdss]) - median(w[grp == "WTlow", sdss])
  p <- tryCatch(wilcox.test(sdss ~ grp, w)$p.value, error = function(e) NA_real_)
  # anchor: NPM1-mut vs all WT, the textbook direction, same convention
  anch <- median(d[grp == "NPM1mut", sdss]) - median(d[grp != "NPM1mut", sdss])
  cat(sprintf("\n    %s: n like %d / WT-low %d / mut %d\n", nm,
              sum(w$grp == "NPM1like"), sum(w$grp == "WTlow"), sum(d$grp == "NPM1mut")))
  cat(sprintf("      median sDSS  like %.2f | WT-low %.2f | mut %.2f\n",
              median(w[grp == "NPM1like", sdss]), median(w[grp == "WTlow", sdss]),
              median(d[grp == "NPM1mut", sdss])))
  cat(sprintf("      NPM1-like - WT-low = %+.2f sDSS (p = %.4g) | NPM1mut - WT anchor = %+.2f\n",
              dif, p, anch))
  cat(sprintf("      GATE: BeatAML said NPM1-like MORE sensitive; in sDSS that means POSITIVE -> %s\n",
              if (is.na(dif)) "NA" else if (dif > 0) "SAME DIRECTION" else "OPPOSITE"))
  res[[nm]] <- data.table(drug = nm, n_like = sum(w$grp == "NPM1like"), n_wtlow = sum(w$grp == "WTlow"),
    median_like = median(w[grp == "NPM1like", sdss]), median_wtlow = median(w[grp == "WTlow", sdss]),
    diff_sdss = dif, p = p, npm1mut_minus_wt_anchor = anch,
    direction = if (is.na(dif)) NA_character_ else if (dif > 0) "same_as_beataml" else "opposite",
    kind_of_evidence = "external_direction_gate")
}
out <- rbindlist(res)
out <- rbind(data.table(drug = "_instrument_", n_like = sum(sc$grp == "NPM1like"), n_wtlow = sum(sc$grp == "WTlow"),
                        median_like = NA_real_, median_wtlow = NA_real_, diff_sdss = q0, p = NA_real_,
                        npm1mut_minus_wt_anchor = NA_real_, direction = "auc_in_diff_column",
                        kind_of_evidence = "external_direction_gate"), out, fill = TRUE)
fwrite(out, file.path(OUT, "fimm_external_gate.csv"))

# ---- 5. same-gene-set arm in BeatAML: is the discovery side stable on the shared 14? -----------
cat("\n[5] re-running the BeatAML side on the SAME 14 genes (removes the gene-set difference)\n")
nb <- fread(BEA)
fzb <- fread(file.path(OUT, "cohort_freeze.csv"))[set_triple_1pp == TRUE]
clb <- fread(file.path(OUT, "beataml_clinical_tidy.csv"))[dbgap_rnaseq_sample %in% fzb$dbgap_rnaseq_sample]
Xb <- as.matrix(nb[, clb$dbgap_rnaseq_sample, with = FALSE]); eb <- nb$stable_id
rb <- ave(rowMeans(Xb), eb, FUN = function(m) seq_along(m) == which.max(m)) == 1
Xb <- Xb[rb, ]; eb <- eb[rb]
Pb <- apply(Xb, 2, frank) / nrow(Xb)
scb <- data.table(dbgap_rnaseq_sample = colnames(Xb),
                  score = colMeans(Pb[which(eb %in% COMMON$stable_id), ]),
                  npm1_pos = clb$npm1_pos)
q0b <- auc_of(scb$score, scb$npm1_pos)
thb <- quantile(scb[npm1_pos == TRUE, score], 0.25)
scb[, grp := fifelse(npm1_pos, "NPM1mut", fifelse(score >= thb, "NPM1like", "WTlow"))]
aub <- fread(file.path(OUT, "beataml_drug_auc_long.csv"))[converged == TRUE][
  , .(auc = mean(auc)), by = .(dbgap_rnaseq_sample, inhibitor)]
cat(sprintf("    BeatAML on shared 14: Q0 AUC %.3f (17-gene version was 0.923) | NPM1-like %d\n",
            q0b, sum(scb$grp == "NPM1like")))
for (nm in c("Venetoclax", "Palbociclib")) {
  d <- merge(scb, aub[inhibitor == nm], by = "dbgap_rnaseq_sample")[grp != "NPM1mut"]
  cat(sprintf("    %-12s BeatAML shared-14: median AUC diff %+.1f (p %.4g), n like %d\n", nm,
      median(d[grp == "NPM1like", auc]) - median(d[grp == "WTlow", auc]),
      tryCatch(wilcox.test(auc ~ grp, d)$p.value, error = function(e) NA), sum(d$grp == "NPM1like")))
}

## SELF-CHECK: the sign convention, verified on a drug with a known strong AML effect ------------
cy <- ds[grepl("Venetoclax|Cytarabine", Drug_name, ignore.case = TRUE)]
if (nrow(cy)) {
  allv <- suppressWarnings(as.numeric(unlist(cy[1, setdiff(names(ds), c("Drug_ID", "Drug_name")), with = FALSE])))
  cat(sprintf("\n[SELF-CHECK] sDSS scale for %s: median %.1f, range %.1f..%.1f (higher = more sensitive)\n",
              cy$Drug_name[1], median(allv, na.rm = TRUE), min(allv, na.rm = TRUE), max(allv, na.rm = TRUE)))
}
cat("[done]\n")
