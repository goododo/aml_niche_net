#!/usr/bin/env Rscript
# =============================================================================
# 00_node_status.R   --   emit the (sample x node) three-state table
#
# Decides, for every sample and every hierarchy bin, whether a missing node is
# missing BIOLOGICALLY or TECHNICALLY. See config_nodestatus.R for the rules and
# their evidence.
#
# This runs BEFORE the projection, on protocol metadata alone, and deliberately
# so: the states must not be inferable from the cell counts they will be applied
# to, or the mask becomes a function of the thing it is supposed to protect.
# Counts are joined in afterwards only to split `present` from `absent_*`.
#
# Input  : 01_preprocess/00_curated_manifest.csv   (protocol fields, resolved)
#          03_hierarchy/perbin/<ds>/<sample>_perbin.csv   (optional; if the
#          projection has run, counts refine present vs absent)
# Output : 03_hierarchy/00_node_status.csv         one row per sample x node
#          03_hierarchy/00_node_status_summary.csv per dataset x node
#
# Usage: Rscript scripts/03_hierarchy/00_node_status.R
# =============================================================================

suppressPackageStartupMessages({library(data.table); library(here)})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))   # HIER_TAB_DIR (per_bin_malignant.csv)
source(here::here("scripts", "config", "config_nodestatus.R"))

MAN <- file.path(DIR_PREPROCESS, "00_curated_manifest.csv")
if (!file.exists(MAN)) stop("run 01_preprocess/00_curated_manifest.R first: ", MAN)
M <- fread(MAN)
message("[1] ", nrow(M), " samples from the curated manifest")

# ---------------------------------------------------------------------------
# [2] Apply the mask rules
# ---------------------------------------------------------------------------
grid <- CJ(idx = seq_len(nrow(M)), node = NODE_ALL)
grid[, `:=`(dataset = M$dataset[idx], sample = M$sample[idx])]
grid[, `:=`(status = "present", masked_by = NA_character_)]

fieldval <- function(f) {
  # the manifest carries the RESOLVED columns under _r suffixes where a config
  # fallback exists; fall back to the raw curated column otherwise
  for (cand in c(paste0(f, "_r"), f)) if (cand %in% names(M)) return(as.character(M[[cand]]))
  rep(NA_character_, nrow(M))
}

for (r in seq_len(nrow(NODE_MASK_RULES))) {
  rule <- NODE_MASK_RULES[r]
  v    <- fieldval(rule$field)
  hit  <- if (rule$regex) grepl(rule$value, v, ignore.case = TRUE)
          else            !is.na(v) & tolower(v) == tolower(rule$value)
  if (!any(hit)) { message("    [rule ", rule$rule_id, "] matched 0 samples"); next }
  nodes <- trimws(strsplit(rule$mask, ",")[[1]])
  grid[idx %in% which(hit) & node %in% nodes,
       `:=`(status = "NA_technical",
            masked_by = fifelse(is.na(masked_by), rule$rule_id,
                                paste(masked_by, rule$rule_id, sep = "+")))]
  message("    [rule ", rule$rule_id, "] ", sum(hit), " samples x ", length(nodes), " nodes masked")
}

# STROMAL: default NA_technical, kept only where the protocol positively supports
# it. Applied after the mask rules because it is a different KIND of rule -- the
# others mask on evidence of removal, this one masks on absence of evidence of
# retention. See the note in config_nodestatus.R.
.srt <- fieldval("sorting"); .cp <- fieldval("cell_prep")
stromal_ok <- (!is.na(.srt) & .srt %in% STROMAL_EVIDENCE_SORTING) |
              (!is.na(.cp)  & grepl(STROMAL_EVIDENCE_CELL_PREP, .cp, ignore.case = TRUE))
grid[node == "Stromal" & !(idx %in% which(stromal_ok)),
     `:=`(status = "NA_technical",
          masked_by = fifelse(is.na(masked_by), "STROMAL_NO_EVIDENCE",
                              paste(masked_by, "STROMAL_NO_EVIDENCE", sep = "+")))]
message("    [rule STROMAL_NO_EVIDENCE] ", sum(!stromal_ok), " samples keep no Stromal node; ",
        sum(stromal_ok), " retain it (", paste(unique(M$dataset[stromal_ok]), collapse = ", "), ")")

# A sample with no protocol information at all cannot be shown to have captured
# ANY node, so nothing can be called biologically absent for it. Loud, because
# masking every node is the most destructive outcome available here.
no_info <- which(is.na(fieldval("sorting")) & is.na(fieldval("cell_prep")))
if (length(no_info)) {
  grid[idx %in% no_info, `:=`(status = "NA_technical", masked_by = "NO_PROTOCOL_INFO")]
  warning(length(no_info), " sample(s) have neither sorting nor cell_prep and are fully masked: ",
          paste(unique(M$sample[no_info]), collapse = ", "), call. = FALSE)
}

grid[, abundance_unreliable :=
       fieldval("sorting")[idx] %in% ABUNDANCE_UNRELIABLE_SORTING]

# ---------------------------------------------------------------------------
# [3] Split present from absent_biological using observed counts, if available
# ---------------------------------------------------------------------------
# NOTE the asymmetry: counts can only ever DEMOTE `present` to
# `absent_biological`. They can never overturn `NA_technical` -- a node the
# protocol removed does not become present because a handful of cells leaked
# through the gate. Cells in a masked node are a purity readout, not evidence.
# WHERE THE COUNTS ACTUALLY COME FROM. This used to read per-sample
# PERBIN_DIR/<ds>/<sample>_perbin.csv. Nothing in the current pipeline writes
# there -- that path belongs to the archived 以前06_hierarchy/d30 script, and
# 03_hierarchy/02_per_bin_malignant.R writes ONE consolidated long table instead.
# So the directory was permanently empty, this block never ran, and the
# three-state table silently degraded to two states: every unmasked node stayed
# `present` and nothing could ever be called biologically absent. It did not
# error, which is exactly why it survived.
pb_one <- file.path(HIER_TAB_DIR, "per_bin_malignant.csv")
pb_files <- list.files(PERBIN_DIR, "_perbin\\.csv$", recursive = TRUE, full.names = TRUE)
read_counts <- function() {
  if (file.exists(pb_one)) {                      # current pipeline
    x <- fread(pb_one)
    if (all(c("dataset", "sample", "hierarchy_bin", "n_cells") %in% names(x))) {
      # FRESHNESS. Counts demote nodes from `present` to `absent_biological`, so a
      # stale table does not merely add noise -- it asserts that specific nodes hold
      # no cells, on the strength of a run against a cohort that no longer exists.
      # The July table covers 130 samples of which 62 are gone and 146 of the
      # current 214 are missing; using it would mark those 146 samples' nodes absent.
      x[, .k := paste(dataset, sample)]
      have <- sum(paste(M$dataset, M$sample) %in% x$.k)
      if (have < nrow(M)) {
        warning(sprintf(paste0("per_bin_malignant.csv covers %d of %d manifest samples (dated %s) ",
                "-- IGNORED. Counts that cover only part of the cohort would mark every ",
                "uncovered node `absent_biological`. Re-run 02_per_bin_malignant.R first."),
                have, nrow(M), format(file.mtime(pb_one), "%Y-%m-%d")), call. = FALSE)
        return(NULL)
      }
      message("[3] counts from ", pb_one, " (", nrow(x), " sample x bin rows)")
      return(x[, .(dataset, sample = as.character(sample),
                   node = hierarchy_bin, n_cells = as.integer(n_cells))])
    }
    warning("per_bin_malignant.csv lacks the expected columns; counts skipped", call. = FALSE)
  }
  if (!length(pb_files)) return(NULL)
  message("[3] counts from ", length(pb_files), " legacy per-sample perbin files")
  rbindlist(lapply(pb_files, function(f) {        # legacy layout, kept as a fallback
    x <- tryCatch(fread(f), error = function(e) NULL); if (is.null(x)) return(NULL)
    nm <- names(x)
    dsc <- intersect(c("dataset", "Dataset"), nm)[1]
    smc <- intersect(c("sample", "Sample"), nm)[1]
    ndc <- intersect(c("hierarchy_bin", "node", "bin"), nm)[1]
    ncc <- intersect(c("n_cells", "n", "count"), nm)[1]
    if (any(is.na(c(dsc, smc, ndc, ncc)))) return(NULL)
    x[, .(dataset = get(dsc), sample = as.character(get(smc)),
          node = get(ndc), n_cells = as.integer(get(ncc)))]
  }), fill = TRUE)
}
cnt <- read_counts()
if (!is.null(cnt) && nrow(cnt)) {
  {
    grid <- merge(grid, cnt, by = c("dataset", "sample", "node"), all.x = TRUE)
    grid[status == "present" & (is.na(n_cells) | n_cells < NODE_PRESENT_MIN_CELLS),
         status := "absent_biological"]
    leak <- grid[status == "NA_technical" & !is.na(n_cells) & n_cells >= NODE_PRESENT_MIN_CELLS]
    message("    counts joined for ", uniqueN(cnt[, paste(dataset, sample)]), " sample(s)")
    if (nrow(leak))
      message("    ", nrow(leak), " masked node(s) still hold >= ", NODE_PRESENT_MIN_CELLS,
              " cells -- gate purity readout, status unchanged")
  }
} else {
  grid[, n_cells := NA_integer_]
  message("[3] no per-bin counts yet (", pb_one, " absent) -- every unmasked node stays ",
          "'present'. This table is TWO-state until 03_hierarchy/02_per_bin_malignant.R has run; ",
          "re-run then to split present vs absent_biological.")
}

# ---------------------------------------------------------------------------
# [4] Write + report
# ---------------------------------------------------------------------------
out <- grid[, .(dataset, sample, node, status, masked_by, abundance_unreliable, n_cells)]
setorder(out, dataset, sample, node)
fwrite_safe(out, file.path(DIR_HIERARCHY, "00_node_status.csv"))

message("\n[4] status counts")
print(out[, .N, by = status][order(-N)])

message("\n    NA_technical nodes per dataset:")
print(dcast(out[status == "NA_technical", .N, by = .(dataset, node)],
            dataset ~ node, value.var = "N", fill = 0))

message("\n    usable CCC nodes per sample (7-node graph, Stromal excluded):")
ccc <- out[node != "Stromal"]
per <- ccc[, .(n_usable = sum(status != "NA_technical")), by = .(dataset, sample)]
print(per[, .N, by = .(dataset, n_usable)][order(dataset, n_usable)])

summ <- out[, .(n_samples = .N,
                present = sum(status == "present"),
                absent_biological = sum(status == "absent_biological"),
                NA_technical = sum(status == "NA_technical")), by = .(dataset, node)]
fwrite_safe(summ, file.path(DIR_HIERARCHY, "00_node_status_summary.csv"))

message("\n    samples with abundance_unreliable (marginal is a pipetting decision):")
print(unique(out[abundance_unreliable == TRUE, .(dataset, sample)])[, .N, by = dataset])

message("\n[done] ", file.path(DIR_HIERARCHY, "00_node_status.csv"))
