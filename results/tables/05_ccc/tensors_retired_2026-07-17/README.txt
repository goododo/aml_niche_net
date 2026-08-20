Retired CellChat tensors, written 2026-07-17 against a superseded 220-row CCC manifest.

Why they were moved out of tensors/:
  06_distance/01_intensity_to_distance.R globs CCC_TENSOR_DIR with no manifest
  filter, so every file here would have entered edge_distance.csv and then FGW.
  05_ccc/03_node_features.R derives has_graph from the same glob.

Of the 70 files:
  55 name samples that are not in ccc_sample_manifest.csv at all (pre-merge
     sublibrary names, e.g. Chen2023/AML103_Niche_Immune, and 39 GSE185381
     per-count-matrix names).
  15 name samples that are in the manifest but ccc_eligible = FALSE
     (exclude_reason = below_min_occupied_bins).

None of the 70 overlaps the 138 eligible samples, so no eligible sample was
served a stale graph. Kept rather than deleted so the count is auditable.
