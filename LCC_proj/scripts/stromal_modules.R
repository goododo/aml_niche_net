# stromal_modules.R ----
# Single source of truth for the de-novo stromal cell-identity modules and the subtype rule.
#
# WHY THIS FILE EXISTS: 10_stromal_denovo.R calibrates and applies these modules; 16_dsp_handoff.R
# has to reproduce exactly the same per-cell calls to build the SpatialDecon stromal profiles. Two
# copies of the definitions would silently diverge the first time either is edited, and the resulting
# reference matrix would no longer describe the cells the paper's numbers came from. So both source
# this file and neither carries its own copy.
#
# The definitions themselves are unchanged from the version 10 was validated with; this file is a
# move, not an edit.

MOD <- list(
  # fibroblast / MSC core. DCN+LUM+COL1A1/2/3 is the canonical marrow-fibroblast signature; the
  # niche factors (CXCL12/LEPR/KITLG) mark the CAR-cell / peri-sinusoidal MSC specifically.
  fibro_msc = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "COL6A3", "PDGFRB", "THY1",
                "SFRP2", "SFRP4", "FAP", "GREM1", "POSTN", "MMP2", "CXCL12", "LEPR", "KITLG", "NGFR"),
  endothelial = c("PECAM1", "CDH5", "EMCN", "VWF", "EGFL7", "CLDN5", "RAMP2", "AQP1"),
  pericyte    = c("ACTA2", "TAGLN", "RGS5", "MYH11", "NOTCH3", "PDGFRB"),
  # SPP1 and ALPL are deliberately NOT here. SPP1 (osteopontin) is expressed by macrophages and by
  # much of the MSC pool, ALPL broadly by MSC; including them made osteolineage swallow the MSC.
  osteolineage= c("BGLAP", "RUNX2", "SP7", "IBSP"),
  # Adipocyte. CAVEAT THAT MUST TRAVEL WITH ANY RESULT FROM THIS MODULE: mature marrow adipocytes are
  # large, buoyant and fragile, so droplet scRNA-seq loses them systematically -- they float off in
  # processing and exceed the droplet size limit. A near-zero count here is therefore NOT evidence
  # that the marrow lacked adipocytes. What droplet data CAN see is the adipogenic-primed
  # LepR+/adipo-CAR MSC, which is scored separately below.
  # STRICT, and the strictness was learned the hard way: a first version also carried CFD, LPL,
  # FABP4 and PPARG, and it relabelled 63% of the stroma-enriched compartment as adipocyte. Those
  # four are not adipocyte-restricted in marrow -- CFD (adipsin) is a general fibroblast/MSC gene,
  # LPL and FABP4 are also macrophage, PPARG is broad. Only genes whose marrow expression is
  # essentially confined to the adipocyte lineage are allowed to make the CALL.
  adipocyte   = c("ADIPOQ", "PLIN1", "PLIN4", "CIDEC", "LEP", "AQP7"),
  # the loose set, used ONLY to grade adipogenic priming inside the MSC pool, never to call a subtype
  adipo_prime = c("CFD", "LPL", "FABP4", "PPARG", "CEBPA", "APOE", "ADIRF"),
  chondrocyte = c("ACAN", "SOX9", "COL2A1", "COL9A3", "COMP"),
  # the exclusion module: every marrow haematopoietic cell carries several of these. A true stromal
  # cell carries none. This is what separates a real MSC from a blast sitting in ambient collagen.
  haematopoietic = c("PTPRC", "LAPTM5", "SRGN", "CORO1A", "CD52", "ARHGDIB", "HCST", "LCP1")
)
SCREEN_GENES <- c("DCN", "LUM", "COL1A1", "CXCL12")

# Count, per cell, how many genes of each module are detected (count > 0). Detection rather than an
# expression cutoff: a binary readout is far more robust to the depth differences across 13 studies.
# `cnt` is a genes x cells sparse matrix.
stromal_module_counts <- function(cnt) {
  gn <- rownames(cnt)
  nd <- lapply(MOD, function(gs) {
    idx <- match(intersect(gs, gn), gn)
    if (!length(idx)) return(rep(0L, ncol(cnt)))
    as.integer(Matrix::colSums(cnt[idx, , drop = FALSE] > 0))
  })
  d <- data.table::as.data.table(nd)
  d[, n_genes := as.integer(Matrix::colSums(cnt > 0))]
  d[]
}

# The subtype rule, verbatim from 10. Order matters: adipocyte and osteolineage are tested BEFORE the
# MSC default, because both arise FROM the MSC and so still carry the fibro_msc genes that would
# otherwise capture them.
stromal_subtype_of <- function(d) {
  data.table::fifelse(d$adipocyte >= 2, "adipocyte",
   data.table::fifelse(d$osteolineage >= 2, "osteolineage",
    data.table::fifelse(d$chondrocyte >= 2, "chondrocyte",
     data.table::fifelse(d$endothelial >= 2 & d$fibro_msc < 4, "endothelial",
      data.table::fifelse(d$pericyte >= 3 & d$fibro_msc < 4, "pericyte", "MSC_fibroblast")))))
}

# The operating point is DERIVED, not hard-coded: read it back from the calibration table 10 wrote,
# by the same rule 10 used (loosest grid point whose CD34-sorted false-positive rate stays under
# 1 in 10,000 cells). Anything downstream therefore tracks a re-calibration automatically.
stromal_operating_point <- function(calibration_csv) {
  cal <- data.table::fread(calibration_csv)
  cal[, acceptable := fpr_sorted < 1e-4]
  data.table::setorder(cal, -acceptable, -rate_enrich)
  op <- cal[acceptable == TRUE][1]
  if (!nrow(op) || is.na(op$min_fibro)) op <- cal[order(fpr_sorted)][1]
  list(min_fibro = op$min_fibro, max_haem = op$max_haem, fpr_sorted = op$fpr_sorted)
}
