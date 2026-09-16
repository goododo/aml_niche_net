#!/usr/bin/env python
# 02_graph_geometry.py ----
# INPUT  : <root>/06_distance/edge_distance.csv   (weight_probsum, 49 directed edges per sample)
# OUTPUT : <root>/06_distance/graph_geometry.csv  (per sample x directed pair: reach, C_geo)
#          <root>/06_distance/graph_geometry_mass.csv (per sample x node: signal mass)
# WHAT IT DOES : turns each sample's local edge weights into a GLOBAL graph geometry, by running
#          a damped random walk on the directed graph, and puts the communication signal into the
#          node mass instead of the cost. It does not touch 01_intensity_to_distance.R; the
#          production C is unchanged and both live side by side.
#
# WHY. 08_scoring/13 measured that the production C carries 2.4x more between-sample variation in
# WHERE its values sit than in WHAT those values are: within-sample rank leaves every sample with
# nearly the same 49 numbers, differing mainly in arrangement. Gromov-Wasserstein re-optimises its
# coupling, which is precisely the operation that discards arrangement, so the pipeline encoded its
# signal in the one dimension its statistic cannot read. Four weight-to-cost formulas failed for
# this shared reason, and 13's D1 measured the consequence: the coupling collapses to within 0.15
# of the product coupling's 1/7 diagonal mass.
#
# THE TWO ROLES WERE SWAPPED. In an OT comparison of graphs the structure should supply the COST
# and the signal should supply the MASS (Nagai et al., Bioinformatics 2025, btaf288, do exactly
# this with a hitting-time cost on a line graph). This pipeline had it the other way round: the LR
# signal was the cost and the cell counts were the mass, and cell composition is the quantity
# already known here to be driven by sample preparation rather than biology.
#
# NO WITHIN-SAMPLE NORMALISATION, AND THAT IS THE POINT. Row-normalising the weight matrix into a
# transition matrix already removes the absolute sequencing-depth scale, which is what the rank
# transform was introduced to achieve (see 01's D2 note). The reachability matrix is a row-
# stochastic probability, so it is comparable across samples with no further rescaling, and the
# value spectrum survives.
#
# Usage : python scripts/06_distance/02_graph_geometry.py [--alpha 0.85]
import os, sys, argparse
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config"))
from graph_geometry import walk_geometry, WALK_ALPHA, WALK_TELEPORT, EPS_MASS

CCC_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
N = len(CCC_NODES)
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
LOG_FLOOR = 1e-12        # only to keep -log finite; the count of entries that need it is reported

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--alpha", type=float, default=WALK_ALPHA,
                help="damping: probability the walk continues at each step")
ap.add_argument("--teleport", type=float, default=WALK_TELEPORT,
                help="at each step, mix this much uniform jump into the transition matrix")
args = ap.parse_args()
D_DST = os.path.join(args.root, "06_distance")

edges = pd.read_csv(os.path.join(D_DST, "edge_distance.csv"))
mi = pd.MultiIndex.from_product([CCC_NODES, CCC_NODES])
W = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                       values="weight_probsum", fill_value=0.0)
     .reindex(columns=mi).fillna(0.0))
KEYS = list(W.index)
print("[0] %d samples x %d directed edges from edge_distance.csv" % (len(KEYS), W.shape[1]))
print("[0] damping alpha = %.2f, teleport = %.3f" % (args.alpha, args.teleport))


# The construction itself lives in scripts/config/graph_geometry.py, imported above, because
# 08_scoring/10 must rebuild it from planted weights with the identical formula. See that file
# for why it is not inlined here.


rows, mrows = [], []
n_dangle_tot = 0
n_floor = 0
rowsum_dev = 0.0
for k in KEYS:
    Wm = W.loc[k].to_numpy(float).reshape(N, N)
    dang = int((Wm.sum(1) <= 0).sum())
    R, Cg, m = walk_geometry(Wm, args.alpha, args.teleport)
    n_dangle_tot += dang
    rowsum_dev = max(rowsum_dev, float(np.abs(R.sum(1) - 1.0).max()))
    n_floor += int((R < LOG_FLOOR).sum())
    stre = Wm.sum(1) + Wm.sum(0)      # reported alongside the mass for inspectability
    for i, a in enumerate(CCC_NODES):
        mrows.append(dict(dataset=k[0], sample=k[1], hierarchy_bin=a,
                          signal_strength=float(stre[i]), signal_mass=float(m[i])))
        for j, b in enumerate(CCC_NODES):
            rows.append(dict(dataset=k[0], sample=k[1], sender_bin=a, receiver_bin=b,
                             reach=float(R[i, j]), C_geo=float(Cg[i, j])))
G = pd.DataFrame(rows)
M = pd.DataFrame(mrows)
G.to_csv(os.path.join(D_DST, "graph_geometry.csv"), index=False)
M.to_csv(os.path.join(D_DST, "graph_geometry_mass.csv"), index=False)

## -- SELF-CHECK 1: R must be row-stochastic, or the walk is not a probability ----
print("\n[SELF-CHECK 1] max |rowsum(R) - 1| over all samples = %.3e" % rowsum_dev)
if rowsum_dev > 1e-10:
    raise SystemExit("the reachability matrix is not row-stochastic; the construction is wrong")
print("               dangling rows completed to uniform : %d of %d node-rows (%.2f%%)"
      % (n_dangle_tot, len(KEYS) * N, 100.0 * n_dangle_tot / (len(KEYS) * N)))
print("               entries needing the log floor      : %d of %d (%.3f%%)"
      % (n_floor, len(KEYS) * N * N, 100.0 * n_floor / (len(KEYS) * N * N)))

## -- SELF-CHECK 2: THE DECISIVE ONE. Does the new cost put its information back into the
## -- value spectrum, where GW can read it, instead of into the arrangement, where it cannot?
Cg_all = G.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                       values="C_geo").reindex(columns=mi).to_numpy(float)
Cr_all = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                            values="C", fill_value=1.0).reindex(columns=mi).fillna(1.0)
          .loc[KEYS].to_numpy(float))
print("\n[SELF-CHECK 2] where does the between-sample information live?")
print("               %-14s %14s %16s %8s" % ("cost matrix", "SD sorted", "SD at edge", "ratio"))
for lab, A in [("production C", Cr_all), ("C_geo (new)", Cg_all)]:
    ss = np.sort(A, axis=1).std(0, ddof=1).mean()
    se = A.std(0, ddof=1).mean()
    print("               %-14s %14.4f %16.4f %8.2f" % (lab, ss, se, se / ss if ss > 0 else np.inf))
print("               ratio near 1 means the value spectrum itself varies between samples,")
print("               which is the part a re-optimised GW coupling can still read.")

## -- SELF-CHECK 3: directionality must survive. rho near 1 would mean the walk symmetrised
## -- the graph, which would throw away the one thing 49 edges have over 28.
iu = np.triu_indices(N, 1)
rr, rc = [], []
for t, k in enumerate(KEYS):
    Cg = Cg_all[t].reshape(N, N)
    Cr = Cr_all[t].reshape(N, N)
    for store, Am in ((rr, Cg), (rc, Cr)):
        a, b = Am[iu], Am.T[iu]
        if np.ptp(a) > 0 and np.ptp(b) > 0:
            store.append(stats.spearmanr(a, b).statistic)
# A HIGH rho here means the two directions are REDUNDANT, which is what a shared abundance
# component looks like, not what good directional information looks like. So report it, but
# judge the geometry on whether it tracks the RAW directional signal instead.
W_all = W.loc[KEYS].to_numpy(float)
lr, lw = [], []
for t in range(len(KEYS)):
    Wm = W_all[t].reshape(N, N); Cg = Cg_all[t].reshape(N, N)
    for a, b in zip(*iu):
        if Wm[a, b] > 0 and Wm[b, a] > 0:
            lr.append(Cg[b, a] - Cg[a, b])          # -log R_ab + log R_ba, i.e. log(R_ab/R_ba)
            lw.append(np.log(Wm[a, b] / Wm[b, a]))
r_sig = stats.pearsonr(np.array(lr), np.array(lw))[0]
print("\n[SELF-CHECK 3] directionality")
print("               redundancy between the two directions (high = the graph is effectively")
print("               undirected, and a shared abundance component is inflating both):")
print("                 production C  median rho = %+.3f" % np.median(rc))
print("                 C_geo         median rho = %+.3f" % np.median(rr))
print("               the test that matters, does C_geo still track the RAW directional signal:")
print("                 corr( log(R_ab/R_ba), log(W_ab/W_ba) ) = %+.3f  over %d pairs"
      % (r_sig, len(lr)))
if r_sig < 0.5:
    raise SystemExit("the geometry has lost the raw directional signal; do not build on it")

print("\n[done] wrote graph_geometry.csv (%d rows) and graph_geometry_mass.csv (%d rows)"
      % (len(G), len(M)))
