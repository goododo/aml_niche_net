# =====================================================================
# f04_refnorm_overview.R  --  autologous-reference (SingleR) stage overview
# ---------------------------------------------------------------------
# Source: ref_norm_summary.csv (from c05_ref_norm_identify.R).
# Panels: routing composition / val_lineage_pos precision audit /
# reference availability (n_ref, frac_ref) with MIN_REF_CELLS line.
# Column names auto-resolved (case-insensitive) and printed for review.
# =====================================================================
.here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                  error = function(e) ".")
if (length(.here) == 0 || .here == "") .here <- "."
source(file.path(.here, "f00_fig_config.R"))

REFNORM_MIN_REF <- 30L     # MIN_REF_CELLS in c00_refnorm_config.R
VAL_LOW_CUT     <- 0.60    # flag suspiciously low precision

stopifnot(file.exists(REFNORM_SUMMARY_CSV))
rn <- fread(REFNORM_SUMMARY_CSV)
message("[f04] columns in ref_norm summary:")
print(names(rn))

pick <- function(cands, required = TRUE) {
  for (cc in cands) { h <- names(rn)[tolower(names(rn)) == tolower(cc)]; if (length(h)) return(h[1]) }
  if (required) stop(sprintf("missing column; tried: %s", paste(cands, collapse = ", ")))
  NA_character_
}
C_ds   <- pick(c("dataset", "Dataset"))
C_dec  <- pick(c("decision"))
C_meth <- pick(c("method"), required = FALSE)
C_nref <- pick(c("n_ref", "n_ref_cells"), required = FALSE)
C_frac <- pick(c("frac_ref", "ref_frac"), required = FALSE)
C_val  <- pick(c("val_lineage_pos", "val_pos", "lineage_pos"), required = FALSE)
C_ntb  <- pick(c("n_singler_TB", "n_TB", "n_singler_tb"), required = FALSE)

setnames(rn, C_ds, "Dataset")
rn[, Dataset := factor(Dataset, levels = sort(unique(Dataset)))]

## -- route exactly like c40: decision + method ----
dec  <- tolower(as.character(rn[[C_dec]]))
meth <- if (!is.na(C_meth)) tolower(as.character(rn[[C_meth]])) else rep("", nrow(rn))
rn[, route := fifelse(dec == "autologous_ok", "autologous",
              fifelse(meth == "sorted_guard", "skip", "external"))]
rn[, route := factor(route, levels = c("autologous", "external", "skip"))]

panels <- list()

## -- A. routing composition per dataset ----
rt <- rn[, .N, by = .(Dataset, route)]
panels$A <- ggplot(rt, aes(Dataset, N, fill = route)) +
  geom_col() +
  scale_fill_manual(values = ROUTE_COLS, name = "reference route") +
  labs(title = "A. Reference routing per dataset", x = NULL, y = "n samples") +
  theme_qc()

## -- B. precision audit: val_lineage_pos ----
if (!is.na(C_val)) {
  rn[, .val := as.numeric(get(C_val))]
  vd <- rn[!is.na(.val)]
  panels$B <- ggplot(vd, aes(Dataset, .val)) +
    geom_boxplot(outlier.shape = NA, fill = "#e9eef2") +
    geom_jitter(aes(color = .val < VAL_LOW_CUT), width = 0.15, size = 0.7, alpha = 0.7) +
    geom_hline(yintercept = VAL_LOW_CUT, linetype = 2, color = "#c0392b") +
    scale_color_manual(values = c(`FALSE` = "#34495e", `TRUE` = "#c0392b"),
                       name = NULL, labels = c("ok", sprintf("< %.2f", VAL_LOW_CUT))) +
    labs(title = "B. SingleR precision audit (val_lineage_pos)",
         x = NULL, y = "lineage-positive fraction") +
    theme_qc()
}

## -- C. reference availability ----
if (!is.na(C_nref) && !is.na(C_frac)) {
  rn[, `:=`(.nref = as.numeric(get(C_nref)), .frac = as.numeric(get(C_frac)))]
  panels$C <- ggplot(rn[is.finite(.nref)], aes(.nref, .frac, color = route)) +
    geom_vline(xintercept = REFNORM_MIN_REF, linetype = 2, color = "#c0392b") +
    geom_point(size = 0.9, alpha = 0.7) +
    scale_x_log10(labels = scales::comma) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = ROUTE_COLS, name = "route") +
    annotate("text", x = REFNORM_MIN_REF, y = Inf, vjust = 1.3, hjust = -0.1,
             label = sprintf("MIN_REF_CELLS = %d", REFNORM_MIN_REF), size = 3, color = "#c0392b") +
    labs(title = "C. Autologous reference availability",
         x = "n reference cells (log10)", y = "fraction of sample") +
    theme_qc() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

## -- D. optional: blast-gate drop (n_singler_TB vs n_ref) ----
if (!is.na(C_ntb) && !is.na(C_nref)) {
  rn[, .ntb := as.numeric(get(C_ntb))]
  panels$D <- ggplot(rn[is.finite(.ntb) & is.finite(.nref)], aes(.ntb, .nref, color = route)) +
    geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey60") +
    geom_point(size = 0.9, alpha = 0.7) +
    scale_x_log10(labels = scales::comma) + scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = ROUTE_COLS, name = "route") +
    labs(title = "D. Blast-gate effect (T/B labeled vs kept as reference)",
         x = "n SingleR T/B (log10)", y = "n reference cells (log10)") +
    theme_qc() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combo <- Reduce(`/`, panels) +
    plot_annotation(title = "AML niche_net  --  autologous reference (SingleR) overview")
  save_fig(combo, "f04_refnorm_overview", w = 11, h = 4 * length(panels), subdir = "04_malignancy/refnorm")
} else {
  for (nm in names(panels)) save_fig(panels[[nm]], paste0("f04_refnorm_overview_", nm), w = 11, h = 4, subdir = "04_malignancy/refnorm")
}
message("[f04] done.")
