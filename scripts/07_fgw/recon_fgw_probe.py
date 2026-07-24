#!/usr/bin/env python
# recon_fgw_probe.py ----  RECON ONLY (not part of the 07_fgw pipeline)
# Verify the installed Python OT (POT) library BEFORE writing 07_fgw, per project discipline:
#   [0] is `ot` (POT) installed? version?
#   [1] real signatures of fused_gromov_wasserstein / *2 / gromov_wasserstein / fgw_barycenters / ...
#       (these vary a LOT across POT versions), + which UNBALANCED FGW variant exists (needed for
#       missing nodes, e.g. a sample lacking Megakaryocyte).
#   [2] does a tiny end-to-end FGW actually compute on toy 7-node graphs?
#   [3] does fgw_barycenters actually compute on a few toy graphs, and what does it return?
# OUTPUT: recon_fgw_probe.txt   (upload this file back)
# Usage:  python scripts/07_fgw/recon_fgw_probe.py
#   (if it reports ot missing: pip install POT   -- pypi is allowed)
import sys, traceback, inspect

OUT = "recon_fgw_probe.txt"
lines = []
def rep(*a):
    s = " ".join(str(x) for x in a)
    lines.append(s); print(s)

## -- Section 0. import + version ----
try:
    import numpy as np
    import ot
    rep("[0] POT (ot) version:", ot.__version__)
    rep("[0] numpy version:", np.__version__)
except Exception as e:
    rep("[0] FATAL: cannot import ot/numpy ->", repr(e))
    rep("    -> install with: pip install POT")
    open(OUT, "w").write("\n".join(lines) + "\n"); sys.exit(1)

## -- Section 1. function availability + signatures ----
def sig(path):
    obj = ot
    try:
        for part in path.split("."):
            obj = getattr(obj, part)
        return str(inspect.signature(obj))
    except Exception as e:
        return "NOT FOUND (%r)" % (e,)

for fn in ["gromov.fused_gromov_wasserstein",
           "gromov.fused_gromov_wasserstein2",
           "gromov.gromov_wasserstein",
           "gromov.gromov_wasserstein2",
           "gromov.fgw_barycenters",
           "gromov.gromov_barycenters",
           "gromov.entropic_fused_gromov_wasserstein",
           "gromov.entropic_gromov_wasserstein",
           # unbalanced candidates (names differ by version):
           "gromov.fused_unbalanced_gromov_wasserstein",
           "gromov.fused_unbalanced_gromov_wasserstein2",
           "gromov.unbalanced_gromov_wasserstein",
           "gromov.entropic_fused_unbalanced_gromov_wasserstein",
           "dist"]:
    rep("[1] ot.%s :: %s" % (fn, sig(fn)))

# scan the gromov submodule for any relevant functions actually present in THIS version
try:
    unb = sorted(n for n in dir(ot.gromov) if "unbalanc" in n.lower())
    bar = sorted(n for n in dir(ot.gromov) if "barycenter" in n.lower())
    rep("[1] ot.gromov 'unbalanc' functions:", unb)
    rep("[1] ot.gromov 'barycenter' functions:", bar)
except Exception as e:
    rep("[1] dir(ot.gromov) failed:", repr(e))

## -- Section 2. tiny end-to-end FGW on toy 7-node graphs ----
rep("\n[2] toy FGW (7 nodes, 3 features)")
rng = np.random.default_rng(491638)
n = 7
def toy():
    F = rng.random((n, 3))                                  # node features (frac_mal, stemness, n_cells-ish)
    A = rng.random((n, n)); C = (A + A.T) / 2; np.fill_diagonal(C, 0.0)  # symmetric structure matrix
    p = np.ones(n) / n                                       # uniform node mass
    return F, C, p
try:
    F1, C1, p1 = toy(); F2, C2, p2 = toy()
    M = ot.dist(F1, F2); M = M / (M.max() + 1e-9)           # feature cost n1 x n2, normalized
    T = ot.gromov.fused_gromov_wasserstein(M, C1, C2, p1, p2, loss_fun="square_loss", alpha=0.5)
    T = np.asarray(T)
    rep("[2] FGW plan shape:", T.shape, "| total mass:", round(float(T.sum()), 4))
    d = ot.gromov.fused_gromov_wasserstein2(M, C1, C2, p1, p2, loss_fun="square_loss", alpha=0.5)
    rep("[2] FGW2 distance (scalar):", round(float(d), 6))
except Exception as e:
    rep("[2] FGW toy FAILED:", repr(e)); rep(traceback.format_exc())

## -- Section 3. tiny FGW barycenter on 3 toy graphs ----
rep("\n[3] toy FGW barycenter (3 graphs -> N=7 nodes)")
try:
    Fs, Cs, ps = [], [], []
    for _ in range(3):
        F, C, p = toy(); Fs.append(F); Cs.append(C); ps.append(p)
    lambdas = [1/3, 1/3, 1/3]
    res = ot.gromov.fgw_barycenters(7, Fs, Cs, ps, lambdas, alpha=0.5)
    rep("[3] fgw_barycenters return type:", type(res).__name__,
        "| len:", (len(res) if hasattr(res, "__len__") else "NA"))
    try:
        X = np.asarray(res[0]); Cb = np.asarray(res[1])
        rep("[3] barycenter feature shape:", X.shape, "| structure shape:", Cb.shape)
    except Exception as e2:
        rep("[3] unpacking barycenter result:", repr(e2))
except Exception as e:
    rep("[3] FGW barycenter toy FAILED:", repr(e)); rep(traceback.format_exc())

open(OUT, "w").write("\n".join(lines) + "\n")
print("[done] wrote", OUT, "-- please upload it")
