#!/usr/bin/env python
# 04_node_feature_distance.py ----
# INPUT  : <root>/07_fgw/fgw_nodes_long.csv                      (138 samples x 7 bins, global-z features)
#          <root>/06_distance/scaccordion_distance__prop7.csv    (canonical sample set + order)
#          <root>/05_ccc/ccc_node_presence.csv                   (present-bin counts, sparsity readout)
# OUTPUT : <root>/06_distance/scaccordion_distance__node_feat.csv     (sample x sample, primary)
#          <root>/06_distance/scaccordion_distance__node_feat150.csv  (sample x sample, sensitivity)
#          <root>/06_distance/node_feature_distance_qc.csv            (one row per arm)
# WHAT IT DOES : adds rung 1b to the pre-registered ablation ladder -- a plain Euclidean distance
#          between corresponding node FEATURE vectors, no transport and no communication. Writes in
#          03_scaccordion_distance.py's output contract so 08_scoring/17 picks it up by filename.
#
# WHY THIS EXISTS. PREREGISTRATION_scaccordion.md AMENDMENT 1 (2026-09-24), which is post hoc and
# says so. The ladder as registered goes from 7 bin cell COUNTS (rung 1) straight to EDGE slots
# (rung 2), so it has no rung for node features richer than counting cells -- which is exactly the
# comparison this project's summary claim rests on ("in the limit, counting cells, are sufficient").
# The prediction fixed before this ran is that node_feat does NOT beat prop7.
#
# WHY NOT EDIT 03. Rungs 5-6 import the vendored scACCorDiON; this rung needs none of it, and
# rerunning that chain to add a Euclidean distance would risk the ladder for no gain.
#
# Usage : python scripts/06_distance/04_node_feature_distance.py [--root ...]
import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

CCC_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
# excluded from every arm, for two different reasons, both stated in AMENDMENT 1:
#   n_cells        IS rung 1 (prop7). Including it makes the comparison circular.
#   frac_malignant is forced to 0 in healthy samples by FGW_ZERO_HEALTHY_MAL.
ALWAYS_EXCLUDE = ["n_cells", "frac_malignant"]
META = ["dataset", "sample", "hierarchy_bin", "timepoint", "healthy",
        "n_cells_raw", "present", "mass", "sparse_flag"]
N_PANEL_EXPECTED = 47

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
args = ap.parse_args()
D_FGW = os.path.join(args.root, "07_fgw")
D_DST = os.path.join(args.root, "06_distance")
D_CCC = os.path.join(args.root, "05_ccc")

# ---- 1. load, and take the sample set from the ladder rather than from this table --------------
nd = pd.read_csv(os.path.join(D_FGW, "fgw_nodes_long.csv"))
prop7 = pd.read_csv(os.path.join(D_DST, "scaccordion_distance__prop7.csv"), index_col=0)
SAMPLES = list(prop7.index)

## -- SELF-CHECK 1: the two sample sets must be identical, not merely overlapping ----------------
## 17_scaccordion_gate1.py takes its sample list from ARMS[0] and reindexes every other arm onto
## it. A silent mismatch would be reindexed to NaN and scored as if it were a distance.
here = set(nd["sample"].unique())
there = set(SAMPLES)
print("[1] ladder samples %d | node table samples %d | symmetric difference %d"
      % (len(there), len(here), len(here ^ there)))
if here != there:
    raise SystemExit("sample sets differ: %s" % sorted(here ^ there)[:10])

# ---- 2. the two feature sets ------------------------------------------------------------------
num = [c for c in nd.columns if c not in META]
panel = [c for c in num
         if c.split("_")[0] in ("st", "pg", "cs", "mt", "pt", "mp")
         and not (c.endswith("_normal") or c.endswith("_malignant"))]
full = [c for c in num if c not in ALWAYS_EXCLUDE]

print("[2] all-cells panel features: %d  (expected %d)" % (len(panel), N_PANEL_EXPECTED))
if len(panel) != N_PANEL_EXPECTED:
    raise SystemExit("panel selector drifted; AMENDMENT 1 fixes this set at %d" % N_PANEL_EXPECTED)
print("    sensitivity set (all numeric minus %s): %d" % (ALWAYS_EXCLUDE, len(full)))
for bad in ALWAYS_EXCLUDE:
    assert bad not in panel and bad not in full, "%s leaked into a feature set" % bad

ARMS = {"node_feat": panel, "node_feat150": full}

# ---- 3. build one 7*k vector per sample, in a fixed bin order ---------------------------------
def stack(feats):
    """samples x (7*k), bins in CCC_NODES order. Absent bins are already mean-imputed upstream."""
    w = nd.pivot_table(index="sample", columns="hierarchy_bin", values=feats, dropna=False)
    # pivot_table gives a (feature, bin) column MultiIndex; fix both orders explicitly
    w = w.reindex(index=SAMPLES, columns=pd.MultiIndex.from_product([feats, CCC_NODES]))
    return w

## -- SELF-CHECK 2: every (sample, bin) cell must exist. 138 x 7 = 966 rows, no holes ------------
cnt = nd.groupby("sample")["hierarchy_bin"].nunique()
print("\n[SELF-CHECK 2] node table completeness")
print("    rows %d (expect %d) | samples with all 7 bins: %d of %d"
      % (len(nd), len(SAMPLES) * 7, int((cnt == 7).sum()), len(SAMPLES)))
if len(nd) != len(SAMPLES) * 7 or not (cnt == 7).all():
    raise SystemExit("node table is not a complete samples x 7 grid")

# ---- 4. distances -----------------------------------------------------------------------------
pres = (pd.read_csv(os.path.join(D_CCC, "ccc_node_presence.csv"))
        .query("present == True").groupby("sample").size().reindex(SAMPLES).fillna(0).to_numpy())
iu = np.triu_indices(len(SAMPLES), 1)
d_prop7 = prop7.reindex(index=SAMPLES, columns=SAMPLES).to_numpy(float)
tv_f = os.path.join(D_DST, "scaccordion_distance__tv.csv")
d_tv = (pd.read_csv(tv_f, index_col=0).reindex(index=SAMPLES, columns=SAMPLES).to_numpy(float)
        if os.path.exists(tv_f) else None)
d_pres = np.abs(pres[:, None] - pres[None, :])          # difference in number of present bins

qc = []
for arm, feats in ARMS.items():
    X = stack(feats).to_numpy(float)
    if not np.isfinite(X).all():
        raise SystemExit("%s: %d non-finite cells in the feature matrix" % (arm, (~np.isfinite(X)).sum()))
    D = squareform(pdist(X, "euclidean"))
    pd.DataFrame(D, index=SAMPLES, columns=SAMPLES).to_csv(
        os.path.join(D_DST, "scaccordion_distance__%s.csv" % arm))

    r_prop7 = float(stats.spearmanr(D[iu], d_prop7[iu]).statistic)
    r_spars = float(stats.spearmanr(D[iu], d_pres[iu]).statistic)
    r_tv = float(stats.spearmanr(D[iu], d_tv[iu]).statistic) if d_tv is not None else np.nan
    qc.append(dict(arm=arm, n_features=len(feats), n_dims=X.shape[1], n_samples=len(SAMPLES),
                   spearman_vs_prop7=r_prop7, spearman_vs_tv=r_tv,
                   spearman_vs_present_bin_diff=r_spars))
    print("\n[4] %s: %d features x 7 bins = %d dims" % (arm, len(feats), X.shape[1]))
    print("    rho vs prop7 (composition)      %+.4f" % r_prop7)
    print("    rho vs tv (edge mass)           %+.4f" % r_tv)
    print("    rho vs present-bin-count diff   %+.4f   <- AMENDMENT 1 mandatory sparsity readout"
          % r_spars)

pd.DataFrame(qc).to_csv(os.path.join(D_DST, "node_feature_distance_qc.csv"), index=False)
print("\n[done] wrote %d arms + node_feature_distance_qc.csv" % len(ARMS))
print("       next: add the arm(s) to ARMS in 08_scoring/17_scaccordion_gate1.py and rerun")
