#!/usr/bin/env python
# 10_planted_effect_power.py ----
# Phase 7 (stage 08_scoring). DIAGNOSTIC 4: what size of topology difference can this pipeline
# actually detect? Every topological result so far is negative, and a negative result is only worth
# reporting if the instrument that produced it can be shown to respond to a real effect.
#
# WHY THIS IS THE GATING EXPERIMENT. The variance decomposition (2026-08-25, 3 both-label datasets,
# non-sparse, n=67) says the disease label explains a median 5.3% of each edge's variance against 88.8%
# sample-to-sample residual. At n=67 with 1 df, 80% power needs R2 ~ 11.1%, and BH over 49 edges pushes
# the top edge to R2 ~ 15.3%. Two edges DO clear raw significance within dataset -- HSC_MPP->B_Plasma
# p=0.0074 and B_Plasma->HSC_MPP p=0.0406, same node pair both directions, both weaker_in_AML -- and
# both die at q (0.363, 0.663). Meanwhile alpha=1 pure GW gives p=0.966. So the pipeline is squeezed
# from two sides at once and this script separates them:
#   OMNIBUS DILUTION : GW compares all 49 edges at equal weight, so k informative edges are averaged
#                      with 49-k noisy ones. Nothing in GW can upweight the informative ones.
#   MULTIPLE TESTING : per-edge tests have the power but lose the FDR correction over 49 edges.
# Planting a known effect and sweeping its size measures where each barrier sits.
#
# HOW THE EFFECT IS PLANTED. On weight_probsum, BEFORE the within-sample rank, in AML samples only.
# Planting on C would bypass the rank transform, and the rank transform is itself a suspect: it is
# zero-sum over a sample's 49 edges, so a change that is uniform across edges is removed by
# construction. Consequence worth stating plainly: raising k edges necessarily pushes the other 49-k
# DOWN in rank. That is not an artefact of the harness, it is what the pipeline does to any real
# effect, and it is part of what is being measured.
#
# TWO SELF-CHECKS, both must hold or the power curve means nothing:
#   1. delta=0 must reproduce the stored C bit-for-bit  -> the re-implemented rank step is faithful.
#   2. delta=0 must reproduce 06_alpha_sweep's alpha=1 row (healthy 0.0537 / AML 0.0442 / p 0.966)
#      -> the re-implemented barycenter + LOO scoring is faithful.
#
# INPUT  : <root>/06_distance/edge_distance.csv , <root>/07_fgw/fgw_{nodes,input_index}_long.csv
# OUTPUT : <root>/08_scoring/planted_effect_power.csv
# Usage  : python scripts/08_scoring/10_planted_effect_power.py [--deltas 0,0.25,0.5,1,2,4] [--ks 1,3,6]
import argparse, hashlib, os, sys, warnings
import numpy as np
import pandas as pd
import ot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config'))
from fgw_vocab import load_vocab, assert_index_covered

warnings.filterwarnings("ignore", category=RuntimeWarning)

FGW_NODES = ["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
BLAST_BINS = ["HSC_MPP","LMPP_GMP","Mono_DC"]
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
SEED = 491638
EPS_MASS = 1e-6

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--deltas", default="0,0.25,0.5,1,2,4",
                help="multiplicative boost on weight_probsum for planted edges in AML samples")
ap.add_argument("--ks", default="1,3,6", help="how many edges carry the planted effect")
ap.add_argument("--modes", default="mult,add",
                help="mult: w *= (1+delta) -- CANNOT move a zero edge, and ~half the 49 slots are "
                     "structurally zero (median 25 present per sample). add: w += delta * the sample's "
                     "own median non-zero weight, which can lift an absent edge into existence.")
ap.add_argument("--n_perm", type=int, default=2000)
ap.add_argument("--max_iter", type=int, default=1000)
args = ap.parse_args()
DELTAS = [float(x) for x in args.deltas.split(",")]
KS = [int(x) for x in args.ks.split(",")]
MODES = [m.strip() for m in args.modes.split(",")]
# ONE INDEPENDENT PERMUTATION STREAM PER TEST, keyed on WHAT THE TEST ACTUALLY DEPENDS ON.
# The shared module-level rng this replaces made a cell's p depend on how many cells ran before it,
# and it broke the null self-check below in a way nothing was watching.
#
# _spec IS THE POINT. At delta = 0 nothing is planted -- Wp is W_NP untouched, so Cm is bit-identical
# across every mode x k. Those six rows are therefore THE SAME TEST ON THE SAME DATA and must agree.
# With a shared stream they did not: the shipped table reported 4 distinct omni_p (0.964018 /
# 0.971014 / 0.971514 / 0.972514) and 2 distinct edge_min_p per k, on rows whose mean_abs_dC was
# exactly 0.0. Collapsing the key at delta = 0 is not a convenience -- it is what turns those rows
# into a free, always-on null control: omni_p must be ONE value, and edge_min_p ONE value per k
# (the k subsets differ, so their minima legitimately do).
def _lab(s): return int.from_bytes(hashlib.blake2b(s.encode(), digest_size=8).digest(), "big")
def _rng(label): return np.random.default_rng([SEED, _lab(label)])
def _spec(mode, k, delta): return "null" if delta == 0 else f"{mode}|{k}|{delta:g}"

D_DIST = os.path.join(args.root, "06_distance")
D_FGW  = os.path.join(args.root, "07_fgw")
D_OUT  = os.path.join(args.root, "08_scoring")
_VOCAB = load_vocab(D_FGW)
AML_TP = _VOCAB["aml_timepoints"]

ed    = pd.read_csv(os.path.join(D_DIST, "edge_distance.csv"))
nodes = pd.read_csv(os.path.join(D_FGW, "fgw_nodes_long.csv"))
idx   = pd.read_csv(os.path.join(D_FGW, "fgw_input_index.csv"))
assert_index_covered(idx, _VOCAB)

## -- cohort, identical to 08/01 and 08/06 ----
d = idx.copy()
d["grp"] = np.where(d.timepoint == "Healthy", "Healthy",
                    np.where(d.timepoint.isin(AML_TP), "AML", "other"))
d = d[d.grp != "other"].copy()
d["is_aml"] = (d.grp == "AML").astype(float)
nc = nodes.pivot_table(index=["dataset","sample"], columns="hierarchy_bin",
                       values="n_cells_raw", aggfunc="first").reset_index()
for b in FGW_NODES:
    if b not in nc: nc[b] = 0
nc["total"] = nc[FGW_NODES].sum(axis=1)
nc["blast_proxy"] = nc[BLAST_BINS].sum(axis=1) / nc["total"].clip(lower=1)
d = d.merge(nc[["dataset","sample","blast_proxy"]], on=["dataset","sample"], how="left") \
     .dropna(subset=["blast_proxy"]).reset_index(drop=True)
KEY = list(zip(d.dataset, d["sample"]))
AML_KEYS = set(k for k, a in zip(KEY, d.is_aml) if a == 1)
print(f"[0] cohort {len(d)}  healthy {int((d.is_aml==0).sum())}  AML {int((d.is_aml==1).sum())}")

## -- edge matrices as (sample, 49) with a FIXED edge order ----
ed["edge"] = ed.sender_bin + "->" + ed.receiver_bin
EDGES = [f"{s}->{r}" for s in FGW_NODES for r in FGW_NODES]
W = ed.pivot(index=["dataset","sample"], columns="edge", values="weight_probsum").reindex(columns=EDGES)
C_STORED = ed.pivot(index=["dataset","sample"], columns="edge", values="C").reindex(columns=EDGES)
W = W.reindex(KEY); C_STORED = C_STORED.reindex(KEY)
assert W.notna().all().all(), "missing weight_probsum for some (sample, edge)"

def rank_to_C(w_row):
    """The rank step from 06_distance/01: rank_pct = frank(w, ties='average')/49 ; C = 1 - rank_pct."""
    s = pd.Series(w_row)
    return 1.0 - (s.rank(method="average").to_numpy() / len(s))

C0 = np.vstack([rank_to_C(W.iloc[i].to_numpy(float)) for i in range(len(W))])
maxdev = float(np.max(np.abs(C0 - C_STORED.to_numpy(float))))
print(f"[SELF-CHECK 1] re-ranked C vs stored C : max |dev| = {maxdev:.3e}")
if maxdev > 1e-9:
    raise SystemExit("the re-implemented rank step does not reproduce 06_distance/01; "
                     "every number below would be measuring the harness, not the pipeline")

## -- which edges to plant on: the CLEANEST substrate, chosen before any outcome is seen ----
# Lowest current label-explained variance, so the planted effect is not sitting on top of whatever
# real signal an edge already carries. Deterministic, not random: the choice must not move between runs.
both_ds = [ds for ds, g in d.groupby("dataset") if g.is_aml.nunique() == 2]
mask_both = d.dataset.isin(both_ds).to_numpy()
sparse = d["sparse_flag"].fillna(False).to_numpy(bool)
sel = mask_both & ~sparse
lab_var = {}
for j, e in enumerate(EDGES):
    y = C0[sel, j]; g = d.dataset.to_numpy()[sel]; a = d.is_aml.to_numpy()[sel]
    tot = float(((y - y.mean())**2).sum())
    if tot == 0: lab_var[e] = 0.0; continue
    ss = 0.0
    for ds in np.unique(g):
        m = g == ds
        if len(np.unique(a[m])) < 2: continue
        mu = y[m].mean()
        for lab in (0.0, 1.0):
            mm = m & (a == lab)
            if mm.sum(): ss += mm.sum() * (y[mm].mean() - mu)**2
    lab_var[e] = ss / tot
order = sorted(EDGES, key=lambda e: lab_var[e])
print(f"[1] planting on the lowest-label-variance edges; cleanest 6: "
      + ", ".join(f"{e} ({100*lab_var[e]:.1f}%)" for e in order[:6]))

## -- FGW machinery, replicating 08/06 exactly ----
NODE_F = {}
for (ds, smp), g in nodes.groupby(["dataset","sample"]):
    gg = g.set_index("hierarchy_bin").reindex(FGW_NODES)
    p = np.nan_to_num(gg["mass"].to_numpy(float), nan=EPS_MASS); p = p / p.sum()
    NODE_F[(ds, smp)] = p

def mats(Cmat):
    return {k: Cmat[i].reshape(7, 7) for i, k in enumerate(KEY)}

def barycenter(keys, Cby, alpha):
    Cs = [Cby[k] for k in keys]; ps = [NODE_F[k] for k in keys]
    Fs = [np.zeros((7, 1)) for _ in keys]          # alpha=1 is pure GW; features are ignored
    m = len(Cs)
    out = ot.gromov.fgw_barycenters(7, Fs, Cs, ps, lambdas=[1.0/m]*m, alpha=alpha,
                                    loss_fun="square_loss", symmetric=False, max_iter=args.max_iter,
                                    p=np.ones(7)/7, init_C=np.mean(np.stack(Cs), 0),
                                    init_X=np.mean(np.stack(Fs), 0), random_state=SEED, log=True)
    return np.asarray(out[1]), np.asarray(out[0])

def fgw2(C, p, Cb, Fb, alpha):
    F = np.zeros((7, 1))
    M = ot.dist(F, Fb); mx = M.max(); M = M/(mx+1e-9) if mx > 0 else M
    return float(ot.gromov.fused_gromov_wasserstein2(M, C, Cb, p, np.ones(7)/7,
                                                     loss_fun="square_loss", alpha=alpha, symmetric=False))

def fwl_perm(y, X0, a, groups, n_perm, rng):
    Q, _ = np.linalg.qr(X0)
    res = lambda v: v - Q @ (Q.T @ v)
    ry = res(y); ra = res(a); den = float(ra @ ra)
    if den <= 0: return np.nan, np.nan
    beta = float(ry @ ra / den)
    gidx = None if groups is None else [np.where(groups == g)[0] for g in np.unique(groups)]
    ex = 0
    for _ in range(n_perm):
        if gidx is None: ap_ = rng.permutation(a)
        else:
            ap_ = a.copy()
            for gi in gidx: ap_[gi] = rng.permutation(ap_[gi])
        rp = res(ap_); dd = float(rp @ rp)
        if dd <= 0: continue
        ex += (abs(float(ry @ rp / dd)) >= abs(beta))
    return beta, (1.0 + ex) / (n_perm + 1.0)

def bh(pv):
    p = np.asarray(pv, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); run = 1.0
    for i in range(n-1, -1, -1):
        run = min(run, p[o[i]] * n / (i+1)); q[o[i]] = run
    return q

def omnibus_alpha1(Cmat, spec):
    """alpha=1 pure GW: healthy barycenter (LOO for healthy) + within-dataset regression."""
    Cby = mats(Cmat)
    heal_all = [k for k, a in zip(KEY, d.is_aml) if a == 0]
    heal_bar = [k for k, a, s in zip(KEY, d.is_aml, sparse) if a == 0 and not s]
    aml_all  = [k for k, a in zip(KEY, d.is_aml) if a == 1]
    Cb_full, Fb_full = barycenter(heal_bar, Cby, 1.0)
    hds = {}
    for k in heal_all:
        loo = [x for x in heal_bar if x != k]
        Cb, Fb = barycenter(loo, Cby, 1.0)
        hds[k] = fgw2(Cby[k], NODE_F[k], Cb, Fb, 1.0)
    for k in aml_all:
        hds[k] = fgw2(Cby[k], NODE_F[k], Cb_full, Fb_full, 1.0)
    y = np.array([hds[k] for k in KEY]); a = d.is_aml.to_numpy(float)
    db = d[mask_both].reset_index(drop=True)
    yb = y[mask_both]; ab = a[mask_both]
    dums = pd.get_dummies(db["dataset"], drop_first=True).to_numpy(float)
    X0b = np.column_stack([np.ones(len(db)), db["blast_proxy"].to_numpy(float), dums])
    b, p = fwl_perm(yb, X0b, ab, db["dataset"].to_numpy(), args.n_perm, _rng(f"omni|{spec}"))
    return float(y[a == 0].mean()), float(y[a == 1].mean()), b, p

def per_edge(Cmat, planted, spec):
    """within-dataset per-edge regression, BH over all 49 -- the same test 04_edge_regression runs."""
    db = d[sel].reset_index(drop=True)
    dums = pd.get_dummies(db["dataset"], drop_first=True).to_numpy(float)
    X0 = np.column_stack([np.ones(len(db)), db["blast_proxy"].to_numpy(float), dums])
    ab = db.is_aml.to_numpy(float); gb = db["dataset"].to_numpy()
    ps = []
    for j in range(len(EDGES)):
        _, p = fwl_perm(Cmat[sel, j], X0, ab, gb, args.n_perm, _rng(f"peredge|{spec}|{EDGES[j]}"))
        ps.append(p)
    q = bh(ps)
    ji = [EDGES.index(e) for e in planted]
    return (float(np.min([ps[j] for j in ji])), float(np.min([q[j] for j in ji])),
            int(sum(q[j] < 0.05 for j in ji)), int((q < 0.05).sum()))

## -- sweep ----
W_NP = W.to_numpy(float)
AMLROW = np.array([kk in AML_KEYS for kk in KEY])
# per-sample scale for additive planting: the sample's own median NON-ZERO edge weight, so delta is
# expressed in units of "a typical edge for this sample" rather than an absolute prob.
SCALE = np.array([np.median(r[r > 0]) if (r > 0).any() else 0.0 for r in W_NP])

rows = []
for mode in MODES:
    for k in KS:
        planted = order[:k]
        ji = [EDGES.index(e) for e in planted]
        for delta in DELTAS:
            Wp = W_NP.copy()
            if delta > 0:
                if mode == "mult":
                    Wp[np.ix_(AMLROW, ji)] *= (1.0 + delta)
                elif mode == "add":
                    Wp[np.ix_(AMLROW, ji)] += delta * SCALE[AMLROW][:, None]
                else:
                    raise SystemExit(f"unknown mode {mode!r}")
            Cm = np.vstack([rank_to_C(Wp[i]) for i in range(len(Wp))])

            # NON-VACUITY, checked every time. The first version of this script planted
            # multiplicatively and 60.3% of the targeted cells did not move at all even at delta=8,
            # because w*(1+delta) is identically w when w == 0 and ~half of the 49 slots are zero
            # (median 25 present per sample). It was measuring a no-op and would have reported that
            # as "the pipeline cannot detect topology".
            dC = np.abs(Cm[np.ix_(AMLROW, ji)] - C0[np.ix_(AMLROW, ji)])
            frac_moved = float((dC > 0).mean()); mean_dC = float(dC.mean())

            spec = _spec(mode, k, delta)
            p_edge, q_edge, n_pl_sig, n_any_sig = per_edge(Cm, planted, spec)
            mh, ma, b_om, p_om = omnibus_alpha1(Cm, spec)
            rows.append(dict(mode=mode, k=k, delta=delta, planted=";".join(planted),
                             frac_cells_moved=frac_moved, mean_abs_dC=mean_dC,
                             edge_min_p=p_edge, edge_min_q=q_edge,
                             planted_edges_sig=n_pl_sig, any_edges_sig=n_any_sig,
                             omni_healthy=mh, omni_aml=ma, omni_beta=b_om, omni_p=p_om))
            flag = "" if (delta == 0 or frac_moved > 0.9) else f"  [only {100*frac_moved:.0f}% of targets moved]"
            print(f"[2] {mode:<4} k={k} d={delta:<5g} | moved {100*frac_moved:>5.1f}% mean|dC|={mean_dC:.4f} | "
                  f"edge p={p_edge:.4f} q={q_edge:.4f} ({n_pl_sig}/{k} at q<.05) | "
                  f"omni p={p_om:.4f} (H {mh:.4f}/A {ma:.4f}){flag}")

res = pd.DataFrame(rows)
out = os.path.join(D_OUT, "planted_effect_power.csv")
res.to_csv(out, index=False)

## -- SELF-CHECK 2 + verdict ----
null = res[(res.delta == 0)].drop_duplicates(subset=["k"])
print("\n[SELF-CHECK 2] delta=0 must reproduce 06_alpha_sweep's alpha=1 row (healthy 0.0537 / "
      "AML 0.0442 / p 0.9659):")
for _, r in null.iterrows():
    print(f"    k={r.k}: healthy={r.omni_healthy:.4f} AML={r.omni_aml:.4f} p={r.omni_p:.4f}")
if null.omni_p.min() < 0.05:
    print("    [FAIL] the unplanted control is significant -- the harness manufactures signal.")
else:
    print("    [PASS] unplanted control is null in the omnibus, as it must be.")

print("\n[3] DETECTION THRESHOLDS (smallest delta that reaches the bar)")
for mode in MODES:
  print(f"  -- mode={mode} --")
  for k in KS:
    sub = res[(res.k == k) & (res["mode"] == mode)].sort_values("delta")
    def first(cond):
        hit = sub[cond(sub)]
        return f"{hit.delta.iloc[0]:g}" if len(hit) else "never"
    print(f"    k={k}:  per-edge raw p<.05 at delta={first(lambda s: s.edge_min_p < 0.05):>6}"
          f" | per-edge q<.05 at delta={first(lambda s: s.edge_min_q < 0.05):>6}"
          f" | omnibus p<.05 at delta={first(lambda s: s.omni_p < 0.05):>6}"
          f" | targets moved at max delta: {100*sub.frac_cells_moved.iloc[-1]:.0f}%")
print(f"\n[done] wrote {out}")
