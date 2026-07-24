#!/usr/bin/env python
# recon_fgw_directed_probe.py ----  RECON ONLY (not part of 07_fgw)
# The ONE point still untested before writing 07 (decisions 1-3 locked):
#   [A] pairwise FGW on a DIRECTED (asymmetric) 7x7 C with symmetric=False  -> for HDS/ATS
#   [B] DIRECTED FGW *barycenter* with symmetric=False, and does Cbar stay ASYMMETRIC (direction kept)?
#   [C] eps-mass trick (absent node ~ 0 mass, fixed 7-node vocab) runs without NaN/Inf?
# If all pass, 07's directed 7x7 + eps-mass + balanced-barycenter design is fully de-risked.
# OUTPUT: recon_fgw_directed_probe.txt   (upload this file back)
# Usage:  python scripts/07_fgw/recon_fgw_directed_probe.py
import sys, traceback
import numpy as np
import ot

OUT = "recon_fgw_directed_probe.txt"; L = []
def rep(*a):
    s = " ".join(str(x) for x in a); L.append(s); print(s)

rep("[0] POT", ot.__version__, "| numpy", np.__version__)
rng = np.random.default_rng(491638)
n = 7

def toy_directed(eps_node=None):
    F = rng.random((n, 3))
    A = rng.random((n, n)); np.fill_diagonal(A, 0.0)   # ASYMMETRIC on purpose (directed): A != A.T
    C = A
    p = np.ones(n) / n
    if eps_node is not None:                            # simulate an absent node -> near-0 mass
        p = p.copy(); p[eps_node] = 1e-6; p = p / p.sum()
    return F, C, p

## -- [A] pairwise directed FGW ----
rep("\n[A] pairwise FGW on DIRECTED C (symmetric=False)")
try:
    F1, C1, p1 = toy_directed(); F2, C2, p2 = toy_directed()
    rep("[A] C1 asymmetric?", not np.allclose(C1, C1.T))
    M = ot.dist(F1, F2); M = M / (M.max() + 1e-9)
    T = np.asarray(ot.gromov.fused_gromov_wasserstein(
        M, C1, C2, p1, p2, loss_fun="square_loss", alpha=0.5, symmetric=False))
    d = float(ot.gromov.fused_gromov_wasserstein2(
        M, C1, C2, p1, p2, loss_fun="square_loss", alpha=0.5, symmetric=False))
    rep("[A] plan sum:", round(float(T.sum()), 4), "| FGW2 dist:", round(d, 6),
        "| NaN in plan:", bool(np.isnan(T).any()))
except Exception as e:
    rep("[A] FAILED:", repr(e)); rep(traceback.format_exc())

## -- [B] directed FGW barycenter ----
rep("\n[B] DIRECTED FGW barycenter (symmetric=False), 4 directed graphs")
try:
    Fs, Cs, ps = [], [], []
    for _ in range(4):
        F, C, p = toy_directed(); Fs.append(F); Cs.append(C); ps.append(p)
    res = ot.gromov.fgw_barycenters(n, Fs, Cs, ps, lambdas=[0.25]*4, alpha=0.5,
                                    loss_fun="square_loss", symmetric=False,
                                    random_state=491638)
    # return order = (X_features, C_structure) per prior probe (len 2)
    Xb = np.asarray(res[0]); Cb = np.asarray(res[1])
    rep("[B] barycenter feature shape:", Xb.shape, "| structure shape:", Cb.shape)
    rep("[B] barycenter structure ASYMMETRIC (direction kept)?", not np.allclose(Cb, Cb.T),
        "| NaN:", bool(np.isnan(Cb).any()))
    # for contrast: what does symmetric=True give?
    resS = ot.gromov.fgw_barycenters(n, Fs, Cs, ps, lambdas=[0.25]*4, alpha=0.5,
                                     loss_fun="square_loss", symmetric=True, random_state=491638)
    CbS = np.asarray(resS[1])
    rep("[B] (symmetric=True barycenter asymmetric?)", not np.allclose(CbS, CbS.T),
        "-- expect False; confirms symmetric flag actually matters")
except Exception as e:
    rep("[B] FAILED:", repr(e)); rep(traceback.format_exc())

## -- [C] eps-mass robustness (absent node) ----
rep("\n[C] eps-mass robustness (fixed 7-node vocab; some graphs miss node 4)")
try:
    Fs, Cs, ps = [], [], []
    for k in range(4):
        F, C, p = toy_directed(eps_node=(4 if k < 2 else None))  # 2 of 4 graphs lack node 4
        Fs.append(F); Cs.append(C); ps.append(p)
    res = ot.gromov.fgw_barycenters(n, Fs, Cs, ps, lambdas=[0.25]*4, alpha=0.5,
                                    loss_fun="square_loss", symmetric=False, random_state=491638)
    Cb = np.asarray(res[1])
    rep("[C] barycenter w/ eps-mass nodes | NaN:", bool(np.isnan(Cb).any()),
        "| Inf:", bool(np.isinf(Cb).any()))
    F1, C1, p1 = toy_directed(eps_node=4); F2, C2, p2 = toy_directed()
    M = ot.dist(F1, F2); M = M / (M.max() + 1e-9)
    d = float(ot.gromov.fused_gromov_wasserstein2(M, C1, C2, p1, p2, alpha=0.5, symmetric=False))
    rep("[C] pairwise FGW2 with eps node:", round(d, 6), "| NaN:", bool(np.isnan(d)))
except Exception as e:
    rep("[C] FAILED:", repr(e)); rep(traceback.format_exc())

open(OUT, "w").write("\n".join(L) + "\n")
print("[done] wrote", OUT, "-- please upload it")
