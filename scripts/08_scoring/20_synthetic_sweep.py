#!/usr/bin/env python
# 20_synthetic_sweep.py ----
# INPUT  : none (synthetic; seeds fixed below)
# OUTPUT : 08_scoring/synthetic_sweep_inertness.csv   (C3 across node count x density)
#          08_scoring/synthetic_sweep_planting.csv    (C2 across node count x fraction x amplitude)
# WHAT IT DOES : re-derives the two STRUCTURAL claims of the benchmark on synthetic graphs, so
#          that "observed on AML marrow" becomes "holds across regimes" (paper skeleton sec 8
#          path A; approved 2026-09-24).
#
# DESIGN, FIXED BEFORE THE FIRST RUN. Two sweeps, both on the same generator:
#   W_ij(sample) = exp(mu_ij + sigma_w * eps),  mu_ij ~ N(0, sigma_b) shared across samples on a
#   fixed support mask of the chosen density. sigma_b = 1.0, sigma_w = 0.35 -- chosen to put the
#   K=7 cell near the production regime and NOT tuned afterwards (SELF-CHECK 1 verifies the
#   anchor, it does not adjust it).
#
#   SWEEP 1 (C3, inertness of the hitting-time construction):
#     K in {5, 7, 10, 15, 22, 30} node types x density in {0.35, 0.65, 0.95}, 30 samples each.
#     Build the paper-arm shared-topology graph + HTD cost on the pooled slots (the same
#     scaccordion_geometry code path as 06_distance/03), then measure off-diagonal CV, max/min,
#     and Spearman(EMD distance, total variation) across the 435 sample pairs.
#     PREDICTION (from the structural argument): all three move TOWARD inertness as K and
#     density grow; the registered inertness rule (max/min < 2 or rho > 0.98) fires ever earlier.
#
#   SWEEP 2 (C2, what GW can and cannot see; uniform marginals so composition contributes zero):
#     K in {5, 7, 10, 15, 22}; 20 control vs 20 planted samples; per-sample cost = within-sample
#     rank transform of W (the production "rank" arm's operation).
#     Two modes (v3):
#       perm      group 2 gets ONE fixed node permutation pi (C -> pi C pi^T) -- the exact
#                 transformation GW quotients out. Readout: the DIRECT invariance check,
#                 median |GW(Ci, pi Cj pi^T) - GW(Ci, Cj)| over 50 pairs (prediction ~1e-9),
#                 plus per-edge recovery at the moved positions (prediction: high).
#       additive  log-space shift a in {0.5, 2.0} on a fraction f in {0.10, 0.25, 0.50} of edges,
#                 swept over BETWEEN-SAMPLE HETEROGENEITY h in {0.0, 0.5, 1.0}: each sample gets
#                 its own mu_i = mu + N(0, h). Readouts: per-edge recovery (BH q<0.05) and the GW
#                 omnibus permutation p (999 shuffles, mean cross - mean within).
#     PREDICTION (v3): per-edge recovery is robust to h (n=20/group absorbs it); GW omnibus
#     detection COLLAPSES as h grows -- the sensitivity gap (per-edge detects, GW does not) opens
#     at high h, which is the regime real cohorts occupy (FINDINGS B.6: between-patient variation
#     is 3.9-5.6x measurement noise).
#
#   DESIGN ITERATIONS, ALL LOGGED (logs kept in 00_project/logs/):
#     v1 "swap" mode permuted VALUES between arbitrary edge positions. Own SELF-CHECK failed,
#        correctly: a value-swap is not a node relabelling and GW is not blind to it. The v1
#        design conflated "arrangement over labelled slots" with "structure up to relabelling".
#     v2 fixed the mode (node permutation) but kept two flaws its run exposed: (a) a group-label
#        omnibus is a noisy, indirect way to show an exact invariance -- the direct pairwise
#        check replaces it; (b) the "spectrum shift" readout is void under the rank transform,
#        which makes every sample's value multiset identical BY CONSTRUCTION -- and the fact
#        that GW still separated additive groups at p=0.001 with identical spectra is itself
#        informative: GW is not spectrum-blind, it is relabelling-invariant, and the real-data
#        blindness (6/6 per-edge vs omnibus p=0.55) is a REGIME phenomenon. Hence v3's h axis.
#
# Usage : PYTHONPATH=/FAST/gr10634/gaozy/external/pylibs \
#         python scripts/08_scoring/20_synthetic_sweep.py
import os
import sys
import itertools
import numpy as np
import pandas as pd
import ot
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config"))
from scaccordion_geometry import stg, htd_cost, emd_matrix, variance_filter, HTD_BETA, TELEPORT, VAR_FILTER_Q

D_OUT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables/08_scoring"
SEED = 491638
SIGMA_B, SIGMA_W = 1.0, 0.35
N_SAMP_INERT, N_GRP, N_PERM = 30, 20, 999

def gen_W(rng, K, density, n, mask=None, mu=None):
    """n samples of a K x K weight matrix on a shared support mask (no self-loops)."""
    if mask is None:
        off = ~np.eye(K, dtype=bool)
        mask = off & (rng.random((K, K)) < density)
        while mask.sum() < 2:                          # degenerate draw at tiny K x density
            mask = off & (rng.random((K, K)) < density)
    if mu is None:
        mu = rng.normal(0.0, SIGMA_B, (K, K))
    W = np.exp(mu[None] + SIGMA_W * rng.normal(size=(n, K, K))) * mask[None]
    return W, mask, mu

# ================================ SWEEP 1: C3 inertness ========================================
rng = np.random.default_rng(SEED)
rows1 = []
print("=== SWEEP 1: hitting-time inertness across (K, density) ===")
print("%4s %8s %6s %10s %9s %9s %8s" % ("K", "density", "slots", "offdiagCV", "max/min", "rho_TV", "INERT?"))
for K, dens in itertools.product([5, 7, 10, 15, 22, 30], [0.35, 0.65, 0.95]):
    W, mask, mu = gen_W(rng, K, dens, N_SAMP_INERT)
    nodes = ["T%02d" % i for i in range(K)]
    slots = ["%s$%s" % (a, b) for a in nodes for b in nodes]
    P_full = pd.DataFrame(W.reshape(N_SAMP_INERT, K * K).T, index=slots,
                          columns=["s%02d" % i for i in range(N_SAMP_INERT)])
    P = variance_filter(P_full, q=VAR_FILTER_Q)
    A = P.to_numpy(float); A = A / A.sum(0, keepdims=True)
    from scipy.spatial.distance import pdist, squareform
    TVm = squareform(pdist(A.T, "cityblock")) / 2.0
    Pt, idx = stg(P, arm="paper", teleport=TELEPORT)
    C = htd_cost(Pt, beta=HTD_BETA)
    D = emd_matrix(P, C)
    off = C[~np.eye(C.shape[0], dtype=bool)]; off = off[off > 0]
    cv, mx = float(off.std() / off.mean()), float(off.max() / off.min())
    iu = np.triu_indices(N_SAMP_INERT, 1)
    rho = float(stats.spearmanr(D[iu], TVm[iu]).statistic)
    inert = (mx < 2.0) or (rho > 0.98)
    print("%4d %8.2f %6d %10.4f %9.3f %9.5f %8s" % (K, dens, P.shape[0], cv, mx, rho, "YES" if inert else "no"))
    rows1.append(dict(K=K, density=dens, n_slots=P.shape[0], offdiag_cv=cv,
                      cost_max_over_min=mx, spearman_vs_tv=rho, inert_rule_fires=inert))
pd.DataFrame(rows1).to_csv(os.path.join(D_OUT, "synthetic_sweep_inertness.csv"), index=False)

## -- SELF-CHECK 1: the K=7 anchor must land near the production regime, or the generator is
## -- telling us about itself rather than about the construction.
a7 = [r for r in rows1 if r["K"] == 7 and r["density"] == 0.65][0]
print("\n[SELF-CHECK 1] K=7, density=0.65 anchor: rho_TV=%.4f (production 7-bin: 0.9968; 22-bin: 0.99972)"
      % a7["spearman_vs_tv"])

# ================================ SWEEP 2: C2 planting =========================================
def rank_C(W):
    """Within-sample rank transform of the off-diagonal weights -> [0,1] cost. Production 'rank'."""
    K = W.shape[0]
    C = np.zeros((K, K))
    off = ~np.eye(K, dtype=bool)
    v = W[off]
    r = stats.rankdata(v) / len(v)
    C[off] = r
    return C

def gw_d(C1, C2, q):
    return float(ot.gromov.gromov_wasserstein2(C1, C2, q, q, loss_fun="square_loss", symmetric=False))

rng = np.random.default_rng(SEED + 1)
rows2 = []
print("\n=== SWEEP 2 (v3): invariance check + heterogeneity-swept planting (uniform marginals) ===")
print("%4s %-9s %4s %5s %5s | %9s %10s" % ("K", "mode", "h", "frac", "amp", "edge_rec", "gw_readout"))
for K in [5, 7, 10, 15, 22]:
    q = np.ones(K) / K
    # ---- perm mode: the DIRECT invariance check ----
    Wc, mask, mu = gen_W(rng, K, 0.65, N_GRP)
    pi = rng.permutation(K)
    while (pi == np.arange(K)).all():
        pi = rng.permutation(K)
    Cc = np.stack([rank_C(w) for w in Wc])
    Wp = gen_W(rng, K, 0.65, N_GRP, mask=mask, mu=mu)[0]
    Cp0 = np.stack([rank_C(w) for w in Wp])            # unpermuted twins
    Cp = Cp0[:, pi][:, :, pi]                          # relabelled group 2
    dev = []
    for t in range(50):
        i, j = rng.integers(N_GRP), rng.integers(N_GRP)
        dev.append(abs(gw_d(Cc[i], Cp[j], q) - gw_d(Cc[i], Cp0[j], q)))
    dev = float(np.median(dev))
    moved = np.argwhere(mask | mask[np.ix_(pi, pi)])
    pv = np.array([stats.ttest_ind(Cc[:, i, j], Cp[:, i, j], equal_var=False).pvalue
                   for i, j in moved])
    pv = pv[np.isfinite(pv)]
    o = np.argsort(pv); qv = np.empty(len(pv))
    qv[o] = np.minimum.accumulate((pv[o] * len(pv) / (np.arange(len(pv)) + 1))[::-1])[::-1]
    rec = float((qv < 0.05).mean())
    print("%4d %-9s %4s %5s %5s | %9.2f %10.2e  <- median |GW(C,piC'pi)-GW(C,C')|"
          % (K, "perm", "-", "-", "-", rec, dev))
    rows2.append(dict(K=K, mode="perm", h=np.nan, frac=np.nan, amp=np.nan,
                      edge_recovery=rec, gw_readout=dev, readout_kind="invariance_dev"))
    # ---- additive mode across heterogeneity ----
    for h in (0.0, 0.5, 1.0):
        for frac in (0.10, 0.25, 0.50):
            for amp in (0.5, 2.0):
                muA = mu[None] + rng.normal(0.0, h, (N_GRP, K, K))
                muB = mu[None] + rng.normal(0.0, h, (N_GRP, K, K))
                Wa = np.exp(muA + SIGMA_W * rng.normal(size=(N_GRP, K, K))) * mask[None]
                Wb = np.exp(muB + SIGMA_W * rng.normal(size=(N_GRP, K, K))) * mask[None]
                pos = np.argwhere(mask)
                n_pl = max(2, int(round(frac * len(pos))))
                pick = pos[rng.choice(len(pos), n_pl, replace=False)]
                Wb[:, pick[:, 0], pick[:, 1]] *= np.exp(amp)
                Ca = np.stack([rank_C(w) for w in Wa])
                Cb = np.stack([rank_C(w) for w in Wb])
                pv = np.array([stats.ttest_ind(Ca[:, i, j], Cb[:, i, j], equal_var=False).pvalue
                               for i, j in pick])
                o = np.argsort(pv); qv = np.empty(n_pl)
                qv[o] = np.minimum.accumulate((pv[o] * n_pl / (np.arange(n_pl) + 1))[::-1])[::-1]
                rec = float((qv < 0.05).mean())
                allC = np.concatenate([Ca, Cb]); n = 2 * N_GRP
                Dg = np.zeros((n, n))
                for i in range(n):
                    for j in range(i + 1, n):
                        Dg[i, j] = Dg[j, i] = gw_d(allC[i], allC[j], q)
                lab = np.array([0] * N_GRP + [1] * N_GRP)
                def stat(l):
                    cross = Dg[np.ix_(l == 0, l == 1)].mean()
                    within = (Dg[np.ix_(l == 0, l == 0)].sum() + Dg[np.ix_(l == 1, l == 1)].sum()) \
                             / (N_GRP * (N_GRP - 1))
                    return cross - within
                obs = stat(lab)
                null = np.array([stat(rng.permutation(lab)) for _ in range(N_PERM)])
                gw_p = float((1 + (null >= obs).sum()) / (1 + N_PERM))
                print("%4d %-9s %4.1f %5.2f %5.1f | %9.2f %10.3f"
                      % (K, "additive", h, frac, amp, rec, gw_p))
                rows2.append(dict(K=K, mode="additive", h=h, frac=frac, amp=amp,
                                  edge_recovery=rec, gw_readout=gw_p, readout_kind="omnibus_p"))
pd.DataFrame(rows2).to_csv(os.path.join(D_OUT, "synthetic_sweep_planting.csv"), index=False)

## -- SELF-CHECK 2: the design's own falsifiers, printed not adjusted ---------------------------
R2 = pd.DataFrame(rows2)
pm = R2[R2["mode"] == "perm"]; ad = R2[R2["mode"] == "additive"]
print("\n[SELF-CHECK 2] the v3 predictions, scored on the table above:")
print("  invariance: median |GW(C,piC'pi^T)-GW(C,C')| <= 1e-9 at every K : %s (max %.2e)"
      % (bool((pm.gw_readout <= 1e-9).all()), pm.gw_readout.max()))
print("  per-edge recovery >= 0.5 in every perm cell : %s (min %.2f)"
      % (bool((pm.edge_recovery >= 0.5).all()), pm.edge_recovery.min()))
for h in (0.0, 0.5, 1.0):
    g = ad[ad.h == h]
    gap = g[(g.edge_recovery >= 0.8) & (g.gw_readout > 0.05)]
    print("  h=%.1f: GW detects %2d/%2d cells | sensitivity gap (edge>=0.8, GW p>0.05): %d"
          % (h, int((g.gw_readout < 0.05).sum()), len(g), len(gap)))
print("\n[done] wrote synthetic_sweep_{inertness,planting}.csv")
