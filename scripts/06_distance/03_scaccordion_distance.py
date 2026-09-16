#!/usr/bin/env python
# 03_scaccordion_distance.py ----
# INPUT  : <root>/05_ccc/tensors/<ds>/<sample>__ccc_cellchat.csv  (LR-level: prob, pval, bins)
#          <root>/05_ccc/ccc_node_presence.csv                    (per-sample per-bin n_cells)
#          <root>/06_distance/edge_distance.csv                   (production weight_probsum; check only)
# OUTPUT : <root>/06_distance/scaccordion_pmat.csv               (line-graph slot x sample)
#          <root>/06_distance/scaccordion_cost__<arm>.csv         (slot x slot, OT arms only)
#          <root>/06_distance/scaccordion_distance__<arm>.csv     (sample x sample, every arm)
#          <root>/06_distance/scaccordion_qc.csv                  (one row per arm)
# WHAT IT DOES : builds the six sample-to-sample distances of the pre-registered ablation ladder
#          (PREREGISTRATION_scaccordion.md), from the bare cell-composition control at one end to
#          the full published scACCorDiON construction at the other, so that a later benchmark can
#          attribute any separation to a specific ingredient rather than to "optimal transport".
#
# WHY THIS EXISTS. FINDINGS section A closes AML-vs-healthy topology across eleven configurations
# of this project's own FGW pipeline. The open question that closes was NOT "is there signal" but
# "is the null a property of the cohort or of our chosen statistic". A published method that
# recovered disease labels on seven cohorts answers that -- but only if it is run faithfully and
# only if what it does is measured rather than assumed.
#
# WHY A LADDER AND NOT ONE NUMBER. Pre-flight measurement on 2026-09-16 found that the ground
# cost scACCorDiON computes is nearly constant off-diagonal, so its distance reproduces plain
# total variation on the same mass vectors at Spearman 0.997 here -- AND at 0.9995 on the authors'
# own published Peng PDAC cohort, which has ten cell types rather than seven. The degeneracy is a
# property of the construction (a near-complete line graph has near-uniform stationary
# distribution), not of this cohort's resolution. So "scACCorDiON separates the groups" would not
# by itself license any claim about optimal transport, and the rungs below are what make the
# attribution possible. SELF-CHECK 4 re-measures this on whatever data it is actually given.
#
# Usage : PYTHONPATH=/FAST/gr10634/gaozy/external/pylibs \
#         python scripts/06_distance/03_scaccordion_distance.py [--root ...] [--jobs 1]
import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config"))
from scaccordion_geometry import (                                          # noqa: E402
    tables_from_tensors, pmat, variance_filter, stg, htd_cost, emd_matrix,
    scaccordion_on_path, PSEUDO, VAR_FILTER_Q, HTD_BETA, TELEPORT, SCACC_REPO)

CCC_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
PVAL_THRESH = 0.05        # 06_distance/01's locked DIST_PVAL_THRESH; see scaccordion_geometry
COMP_PSEUDO = 0.5         # cell-count pseudocount for the CLR in the prop7 rung
COST_ABS_MAX = 1e6        # any |C| above this means getCTD's uninitialised-memory path fired

ARMS = ["prop7", "tabular", "tv", "corrot", "dwot_shipped", "dwot_fixed"]

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--pval", type=float, default=PVAL_THRESH)
ap.add_argument("--filter_q", type=float, default=VAR_FILTER_Q)
args = ap.parse_args()
D_CCC = os.path.join(args.root, "05_ccc")
D_DST = os.path.join(args.root, "06_distance")

sha = "unknown"
sha_f = os.path.join(SCACC_REPO, ".VENDORED_SHA")
if os.path.exists(sha_f):
    sha = open(sha_f).read().strip()
print("[0] scACCorDiON vendored at %s\n    commit %s" % (SCACC_REPO, sha))

# ---- 1. adapter -------------------------------------------------------------------------------
paths = sorted(glob.glob(os.path.join(D_CCC, "tensors", "*", "*__ccc_cellchat.csv")))
tbls = tables_from_tensors(paths, pval_thresh=args.pval)
empty = [k for k, v in tbls.items() if len(v) == 0]
if empty:
    print("[1] %d samples have no significant LR pair and are dropped: %s" % (len(empty), empty))
    tbls = {k: v for k, v in tbls.items() if len(v)}
SAMPLES = sorted(tbls)
print("[1] %d samples -> scACCorDiON input tables (pval < %.3g)" % (len(SAMPLES), args.pval))

P_FULL = pmat(tbls, nodes=CCC_NODES)
P = variance_filter(P_FULL, q=args.filter_q)[SAMPLES]
dropped = [s for s in P_FULL.index if s not in set(P.index)]
print("[1] line-graph slots: %d built -> %d kept after Accordion's filter=%.2f"
      % (P_FULL.shape[0], P.shape[0], args.filter_q))
print("    dropped by the variance filter (%d): %s" % (len(dropped), ", ".join(dropped)))
P.to_csv(os.path.join(D_DST, "scaccordion_pmat.csv"))

## -- SELF-CHECK 1: their aggregation must reproduce our production edge weights ----------------
## If this fails, the distance below is not being computed on the pipeline's own edges and none
## of it is comparable to the eleven configurations already in FINDINGS section A.
ed = pd.read_csv(os.path.join(D_DST, "edge_distance.csv"))
prod = (ed.pivot_table(index=["sender_bin", "receiver_bin"], columns="sample",
                       values="weight_probsum", fill_value=0.0))
prod.index = ["%s$%s" % (a, b) for a, b in prod.index]
common_s = [s for s in SAMPLES if s in prod.columns]
common_e = [e for e in P_FULL.index if e in prod.index]
delta = np.abs(P_FULL.loc[common_e, common_s].to_numpy() - prod.loc[common_e, common_s].to_numpy())
print("\n[SELF-CHECK 1] scACCorDiON's groupby-sum vs production weight_probsum")
print("               compared %d slots x %d samples ; max |difference| = %.3e"
      % (len(common_e), len(common_s), delta.max()))
if delta.max() > 1e-8:
    raise SystemExit("the adapter does not reproduce production edge weights; stop")

# ---- 2. the ladder ----------------------------------------------------------------------------
A = P.to_numpy(float)
A = A / A.sum(0, keepdims=True)                       # per-sample distribution over slots
TV = squareform(pdist(A.T, "cityblock")) / 2.0        # the rung every OT arm is measured against

D, C_OF = {}, {}

# rung 1 -- composition only. No communication at all: if this separates the groups as well as
# anything above it, the ladder has no methodological finding and the signal is blast fraction.
np_ = pd.read_csv(os.path.join(D_CCC, "ccc_node_presence.csv"))
comp = (np_.pivot_table(index="sample", columns="hierarchy_bin", values="n_cells", fill_value=0)
        .reindex(index=SAMPLES, columns=CCC_NODES).fillna(0.0).to_numpy(float))
comp = comp + COMP_PSEUDO
clr = np.log(comp) - np.log(comp).mean(1, keepdims=True)
D["prop7"] = squareform(pdist(clr, "euclidean"))

# rungs 2-3 -- the same edge vector, no transport. Their paper's "Tabular" competitor, and the
# total-variation distance that a linear OT with a constant ground cost is proportional to.
D["tabular"] = squareform(pdist(A.T, "euclidean"))
D["tv"] = TV

# rung 4 -- transport, but with a cost that has no graph geometry in it (their CORR-OT).
C_OF["corrot"] = squareform(pdist(P.to_numpy(float), "correlation"))
D["corrot"] = emd_matrix(P, C_OF["corrot"])

# rung 5 -- THEIR CODE, RUN AS SHIPPED. Imported, not reimplemented.
scaccordion_on_path()
from scaccordion.tools.Accordion import Accordion                           # noqa: E402
acc = Accordion({k: tbls[k] for k in SAMPLES}, weight="lr_means",
                filter=args.filter_q, filter_mode="edge", pseudo=PSEUDO)
acc.compute_cost(mode="HTD", beta=HTD_BETA)
acc.compute_wassestein(cost="HTD_%s" % HTD_BETA, algorithm="emd")
D["dwot_shipped"] = acc.wdist["HTD_%s" % HTD_BETA].reindex(index=SAMPLES, columns=SAMPLES).to_numpy(float)
C_OF["dwot_shipped"] = np.asarray(acc.Cs["HTD_%s" % HTD_BETA], float)
import networkx as nx                                                       # noqa: E402
shipped_order = list(acc.expgraph.nodes())
n_misordered = sum(1 for a, b in zip(shipped_order, list(acc.p.index)) if a != b)

# rung 6 -- the construction the PAPER describes. Three documented deviations from rung 5, plus
# the node-order alignment their code loses. See scaccordion_geometry.py's header.
Pt, idx = stg(P, arm="paper", teleport=TELEPORT)
assert idx == list(P.index), "paper-arm cost is not aligned to its own marginals"
C_OF["dwot_fixed"] = htd_cost(Pt, beta=HTD_BETA)
D["dwot_fixed"] = emd_matrix(P, C_OF["dwot_fixed"])

## -- SELF-CHECK 2: the published cost function has a known uninitialised-memory path ----------
## `-np.log10(Aht, where=Aht != 0)` is called without `out=`, so entries where the symmetrised
## matrix is exactly zero are whatever was in memory. Do not patch it; detect it.
print("\n[SELF-CHECK 2] cost matrix integrity")
bad = False
for arm, C in C_OF.items():
    fin = np.isfinite(C).all()
    sym = float(np.abs(C - C.T).max())
    dia = float(np.abs(np.diag(C)).max())
    neg = int((C < -1e-12).sum())
    big = int((np.abs(C) > COST_ABS_MAX).sum())
    print("    %-13s finite=%-5s max|C-C'|=%.2e max|diag|=%.2e negatives=%d |C|>%.0e: %d"
          % (arm, fin, sym, dia, neg, COST_ABS_MAX, big))
    bad |= (not fin) or big > 0
if bad:
    raise SystemExit("a cost matrix is non-finite or absurdly large; getCTD's `where=` path fired")

## -- SELF-CHECK 3: the node-ordering defect, measured rather than assumed ----------------------
print("\n[SELF-CHECK 3] node ordering inside the shipped implementation")
print("    Accordion.compute_cost builds C in networkx insertion order; compute_wassestein")
print("    passes marginals in p.index order. Positions that differ: %d of %d"
      % (n_misordered, len(shipped_order)))
print("    (the dwot_fixed arm is aligned by construction and asserts it above)")

## -- SELF-CHECK 4: THE DECISIVE ONE. Is the ground geometry doing any work at all? -------------
## A linear OT whose cost is constant off the diagonal is exactly proportional to total variation.
## Measure the spread of each cost and how much of the sample ordering survives beyond TV.
iu = np.triu_indices(len(SAMPLES), 1)
print("\n[SELF-CHECK 4] does the ground cost contribute beyond total variation?")
print("    %-13s %10s %9s %11s %11s" % ("arm", "offdiagCV", "max/min", "r vs TV", "rho vs TV"))
qc = []
for arm in ARMS:
    C = C_OF.get(arm)
    if C is None:
        cv = mx = np.nan
    else:
        off = C[~np.eye(C.shape[0], dtype=bool)]
        off = off[off > 0]
        cv, mx = float(off.std() / off.mean()), float(off.max() / off.min())
    r = float(stats.pearsonr(D[arm][iu], TV[iu])[0])
    rho = float(stats.spearmanr(D[arm][iu], TV[iu]).statistic)
    print("    %-13s %10.4f %9.3f %11.5f %11.5f" % (arm, cv, mx, r, rho))
    qc.append(dict(arm=arm, n_samples=len(SAMPLES), n_slots=P.shape[0], offdiag_cv=cv,
                   cost_max_over_min=mx, pearson_vs_tv=r, spearman_vs_tv=rho,
                   vendored_commit=sha))
print("    PRE-REGISTERED READING: cost_max_over_min < 2.0 OR spearman_vs_tv > 0.98 means the")
print("    geometry is inert and that arm may not be credited with an optimal-transport result.")

# ---- 3. write ---------------------------------------------------------------------------------
for arm in ARMS:
    pd.DataFrame(D[arm], index=SAMPLES, columns=SAMPLES).to_csv(
        os.path.join(D_DST, "scaccordion_distance__%s.csv" % arm))
for arm, C in C_OF.items():
    pd.DataFrame(C, index=list(P.index), columns=list(P.index)).to_csv(
        os.path.join(D_DST, "scaccordion_cost__%s.csv" % arm))
qcdf = pd.DataFrame(qc)
qcdf["n_misordered_shipped"] = n_misordered
qcdf["filter_q"] = args.filter_q
qcdf["pval_thresh"] = args.pval
qcdf.to_csv(os.path.join(D_DST, "scaccordion_qc.csv"), index=False)
print("\n[done] %d distance matrices, %d cost matrices, %d samples x %d slots"
      % (len(ARMS), len(C_OF), len(SAMPLES), P.shape[0]))
