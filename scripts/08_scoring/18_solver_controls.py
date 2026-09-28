#!/usr/bin/env python
# 18_solver_controls.py ----
# INPUT  : 06_distance/edge_distance.csv, 07_fgw/{fgw_input_index,fgw_nodes_long}.csv
# OUTPUT : 08_scoring/solver_controls_{A,B,C}.csv
# WHAT IT DOES : three controls that separate "the GW objective cannot distinguish node
#          arrangements" from "the POT solver returned its initialisation".
#
# WHY. The C2 claim was originally narrated as "the coupling collapses to the product, where the
# moment identity applies". That narration was WRONG twice over, caught 2026-09-24: (i) the
# converged cost is far below the cost at the product coupling, so the optimum is not the product
# and the moment identity (which holds only AT that point) cannot carry the claim; (ii) diagonal
# mass near 1/7 says "does not favour identity", not "is near p (x) q". The defensible mechanism is
# simpler: GW minimises over ALL couplings, so it is relabelling-invariant BY CONSTRUCTION, and on
# a labelled space it always finds a label-ignoring matching that makes two graphs look at least as
# similar as the label-respecting one. Control C quantifies exactly that. Controls A and B answer
# the solver-artifact worry (entropic epsilon is moot: every call in this project is the exact
# conditional gradient, ot.gromov.fused_gromov_wasserstein, loss_fun="square_loss", no epsilon).
#
#   A  one-hot node features, alpha in {->0 (exact EMD), 0.1, 0.5}: the coupling must reach the
#      attainable diagonal cap sum_i min(p_i, q_i). Tests that the harness moves when the
#      objective has a gradient. (q is the uniform barycentre marginal, so the cap is < 1 by
#      design of the production pipeline, not by solver failure.)
#   B  alpha=1, seven initialisations per sample (product default, max-diagonal, 5 random
#      transport-polytope vertices): if final costs agree, the OPTIMUM is init-independent.
#   C  the load-bearing number: cost at the UNOPTIMISED label-respecting coupling vs the
#      converged optimum, per sample. "Re-matching makes the graphs look X% more similar than
#      matching by label" is the absorption mechanism, quantified.
#
# Usage : python scripts/08_scoring/18_solver_controls.py
import os, sys
import numpy as np
import pandas as pd
import ot

ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config"))
from distance_variants import weights_to_C

FGW_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
D_OUT = os.path.join(ROOT, "08_scoring")
SEED = 491638
rng = np.random.default_rng(SEED)
print("[0] POT %s | every solve: ot.gromov.fused_gromov_wasserstein, square_loss, exact CG" % ot.__version__)

# ---- load exactly as 13_gw_blindness.py D1 does (same files, same rank arm, same barycentre) ----
edges = pd.read_csv(os.path.join(ROOT, "06_distance", "edge_distance.csv"))
index = pd.read_csv(os.path.join(ROOT, "07_fgw", "fgw_input_index.csv"))
nodes_all = pd.read_csv(os.path.join(ROOT, "07_fgw", "fgw_nodes_long.csv"))
NODE_P = {}
for (ds, smp), g in nodes_all.groupby(["dataset", "sample"]):
    pm = np.nan_to_num(g.set_index("hierarchy_bin").reindex(FGW_NODES)["mass"].to_numpy(float), nan=1e-6)
    NODE_P[(ds, smp)] = pm / pm.sum()
W = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                       values="weight_probsum", fill_value=0.0)
     .reindex(columns=pd.MultiIndex.from_product([FGW_NODES, FGW_NODES])).fillna(0.0))
NL = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                        values="n_lr_sig", fill_value=0.0)
      .reindex(columns=pd.MultiIndex.from_product([FGW_NODES, FGW_NODES])).fillna(0.0))
keys = [(r.dataset, r["sample"]) for _, r in index.iterrows()
        if (r.dataset, r["sample"]) in W.index and (r.dataset, r["sample"]) in NODE_P]
CMAT = {k: weights_to_C(W.loc[k].to_numpy(float), "rank", NL.loc[k].to_numpy(float)).reshape(7, 7)
        for k in keys}
sparse = nodes_all.groupby(["dataset", "sample"])["sparse_flag"].first().fillna(False).to_dict()
is_heal = index.set_index(["dataset", "sample"])["timepoint"].eq("Healthy").to_dict()
heal = [k for k in keys if is_heal.get(k) and not sparse.get(k, False)]
q = np.ones(7) / 7.0
Cs = [CMAT[k] for k in heal]
out = ot.gromov.fgw_barycenters(7, [np.zeros((7, 1)) for _ in heal], Cs, [NODE_P[k] for k in heal],
                                lambdas=[1.0 / len(heal)] * len(heal), alpha=1.0,
                                loss_fun="square_loss", symmetric=False, max_iter=1000, p=q,
                                init_C=np.mean(np.stack(Cs), 0), init_X=np.zeros((7, 1)),
                                random_state=SEED, log=True)
Cb = np.asarray(out[1])
print("[0] %d samples | healthy barycentre from %d donors -- identical construction to gw_blindness D1"
      % (len(keys), len(heal)))

I7 = np.eye(7)
def gw_cost(C1, C2, T):
    L = (C1[:, None, :, None] - C2[None, :, None, :]) ** 2
    return float(np.einsum("ijkl,ij,kl->", L, T, T))

# ---- A: one-hot snap test --------------------------------------------------------------------
M1 = ot.dist(I7, I7); M1 = M1 / (M1.max() + 1e-9)
ra = []
for k in keys:
    p = NODE_P[k]
    row = dict(dataset=k[0], sample=k[1], cap=float(np.minimum(p, q).sum()),
               diag_emd=float(np.trace(ot.emd(p, q, M1))))
    for al in (0.1, 0.5):
        T = np.asarray(ot.gromov.fused_gromov_wasserstein(M1, CMAT[k], Cb, p, q,
                       loss_fun="square_loss", alpha=al, symmetric=False))
        row["diag_a%g" % al] = float(np.trace(T))
    ra.append(row)
A = pd.DataFrame(ra); A.to_csv(os.path.join(D_OUT, "solver_controls_A.csv"), index=False)
print("\n[A] one-hot features: median attainable cap %.3f" % A.cap.median())
for c in ["diag_emd", "diag_a0.1", "diag_a0.5"]:
    print("    %-10s median diag mass %.3f | fraction of cap %.3f" % (c, A[c].median(), (A[c] / A.cap).median()))
# NOTE: at alpha->0 with the discrete one-hot metric, EMD(p, q, 1-I) equals the total-variation
# distance TV(p, q). With identity pinned, the feature end of this distance IS a composition
# distance -- checked numerically here:
tv_err = max(abs(A.loc[i, "diag_emd"] - (1.0 - 0.5 * np.abs(NODE_P[k2] - q).sum()) * 0 - A.loc[i, "cap"])
             for i, k2 in enumerate(keys))
tv_chk = max(abs((1.0 - A.loc[i, "diag_emd"]) - 0.5 * np.abs(NODE_P[k2] - q).sum())
             for i, k2 in enumerate(keys))
print("    identity check: (1 - diag mass at EMD) == TV(p, uniform) to %.1e" % tv_chk)

# ---- B: multi-init ---------------------------------------------------------------------------
M0 = np.zeros((7, 7))
rb = []
for k in keys:
    p = NODE_P[k]
    inits = {"product": None, "maxdiag": ot.emd(p, q, 1.0 - I7)}
    for j in range(5):
        inits["rand%d" % j] = ot.emd(p, q, rng.random((7, 7)))
    costs = {}
    diags = {}
    for nm, g0 in inits.items():
        T = np.asarray(ot.gromov.fused_gromov_wasserstein(M0, CMAT[k], Cb, p, q,
                       loss_fun="square_loss", alpha=1.0, symmetric=False, G0=g0))
        costs[nm], diags[nm] = gw_cost(CMAT[k], Cb, T), float(np.trace(T))
    c = np.array(list(costs.values()))
    rb.append(dict(dataset=k[0], sample=k[1], cost_prod=costs["product"],
                   cost_best=float(c.min()), rel_spread=float((c.max() - c.min()) / max(abs(np.median(c)), 1e-30)),
                   rel_gain_over_prod=float((costs["product"] - c.min()) / max(abs(costs["product"]), 1e-30)),
                   diag_from_prod=diags["product"], diag_from_maxdiag=diags["maxdiag"]))
B = pd.DataFrame(rb); B.to_csv(os.path.join(D_OUT, "solver_controls_B.csv"), index=False)
print("\n[B] alpha=1, 7 inits x %d samples:" % len(B))
print("    final-cost relative spread: median %.2e | max %.2e | inits beating product by >1e-4: %d/%d"
      % (B.rel_spread.median(), B.rel_spread.max(), (B.rel_gain_over_prod > 1e-4).sum(), len(B)))
print("    converged diag mass: from product %.3f | from max-diag %.3f (solver LEAVES the diagonal)"
      % (B.diag_from_prod.median(), B.diag_from_maxdiag.median()))

# ---- C: the load-bearing number --------------------------------------------------------------
rc = []
for k in keys:
    p = NODE_P[k]
    T_lab = ot.emd(p, q, 1.0 - I7)
    T_opt = np.asarray(ot.gromov.fused_gromov_wasserstein(M0, CMAT[k], Cb, p, q,
                       loss_fun="square_loss", alpha=1.0, symmetric=False))
    c_lab, c_opt, c_prd = gw_cost(CMAT[k], Cb, T_lab), gw_cost(CMAT[k], Cb, T_opt), gw_cost(CMAT[k], Cb, np.outer(p, q))
    rc.append(dict(dataset=k[0], sample=k[1], cost_label=c_lab, cost_opt=c_opt, cost_product=c_prd,
                   excess_label=(c_lab - c_opt) / c_opt, excess_product=(c_prd - c_opt) / c_opt))
C = pd.DataFrame(rc); C.to_csv(os.path.join(D_OUT, "solver_controls_C.csv"), index=False)
print("\n[C] label-respecting coupling vs converged optimum (alpha=1, rank arm):")
print("    median excess cost of matching BY LABEL: %+.1f%%  (IQR %+.1f%%..%+.1f%%)"
      % (100 * C.excess_label.median(), 100 * C.excess_label.quantile(.25), 100 * C.excess_label.quantile(.75)))
print("    label coupling strictly worse than optimum: %d/%d | worse than even the product: %d/%d"
      % ((C.cost_label > C.cost_opt * (1 + 1e-9)).sum(), len(C), (C.cost_label > C.cost_product).sum(), len(C)))
print("    product coupling excess: %+.1f%% median (the optimum is NOT the product either)"
      % (100 * C.excess_product.median()))
print("\n[done] wrote solver_controls_{A,B,C}.csv to %s" % D_OUT)
