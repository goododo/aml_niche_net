#!/usr/bin/env python
# 02_fgw_align.py ----
# Phase 6 (stage 07_fgw), Python side. Reads the 3 tidy CSVs from 01, builds per-sample DIRECTED FGW
# inputs (C 7x7, F 7x3 z-scored, p 7), computes healthy & AML consensus barycenters (symmetric=False,
# balanced, non-sparse samples), and pairwise patient scores:
#     HDS = FGW2(patient, B_healthy)      ATS = FGW2(patient, B_AML)
# for ALL samples (sparse ones flagged, not dropped).
#
# LOCKED (config_fgw.R, mirrored here as constants -- keep in sync; reindex-by-name makes row order safe):
#   directed 7x7 C (symmetric=False), fixed 7-node vocab, eps-mass absent nodes, mass mode from 01,
#   alpha main 0.5. Self-loops (diagonal) kept. Barycenter init = aligned element-wise mean (nodes are
#   already label-aligned across samples -> keeps barycenter node i ~ bin i, interpretable + stable;
#   FGW's OT still allows soft feature/structure re-matching = the cross-heterogeneity convergence claim).
#
# Verification prints (the checklist from the directed-barycenter probe, now on REAL data):
#   (1) barycenter converged, no NaN;
#   (2) barycenter DIRECTIONALITY survived -- report asymmetry + a specific primitive->immune edge
#       (LMPP_GMP->B_Plasma vs reverse) in both barycenters;
#   (3) HDS/ATS sane -- healthy mean HDS should be low; HDS by timepoint previews H3.
#
# INPUT  : <input_dir>/fgw_{nodes,edges,input_index}_long.csv   (from 01; default DIR_FGW)
# OUTPUT : <input_dir>/barycenters.npz        (C_healthy,F_healthy,C_aml,F_aml,p,nodes,features,alpha)
#          <input_dir>/patient_scores.csv     (dataset,sample,timepoint,sparse_flag,HDS,ATS,
#                                             HDS_norm,HDS_par,HDS_perp,cos_theta,alpha,mass_mode)
# Usage  : python scripts/07_fgw/02_fgw_align.py [--input_dir DIR] [--alpha 0.5] [--test N]
import argparse, os, sys, json
import numpy as np
import pandas as pd
import ot

# The timepoint vocabulary is LOADED, not literal. A hard-coded set here silently deleted 34 of
# 214 samples (64% of the treated arm) after CANONICAL_TIMEPOINTS changed on 2026-08-04.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config'))
from fgw_vocab import load_vocab, load_features, assert_index_covered


# ---- LOCKED constants (must match config_fgw.R) ----
FGW_NODES    = ["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
DEFAULT_INPUT_DIR = "/FAST/gr10634/gaozy/aml_niche_net/results/tables/07_fgw"
AML_TIMEPOINTS = None  # set from fgw_vocab.json below -- see scripts/config/fgw_vocab.py
SEED = 491638
EXCLUDE_SPARSE_FROM_BARY = True

ap = argparse.ArgumentParser()
ap.add_argument("--input_dir", default=DEFAULT_INPUT_DIR)
ap.add_argument("--alpha", type=float, default=0.5)
ap.add_argument("--test", type=int, default=0, help="if >0, use at most N samples per group (smoke test)")
ap.add_argument("--max_iter", type=int, default=1000)
args = ap.parse_args()

nodes = pd.read_csv(os.path.join(args.input_dir, "fgw_nodes_long.csv"))
edges = pd.read_csv(os.path.join(args.input_dir, "fgw_edges_long.csv"))
index = pd.read_csv(os.path.join(args.input_dir, "fgw_input_index.csv"))

_VOCAB = load_vocab(args.input_dir)
# The FEATURE set is LOADED, not literal -- same reason as the timepoints above. Eight files
# carried the hard-coded triple while 01 wrote the real set into fgw_vocab.json; changing
# FGW_FEATURES in config_fgw.R silently left every one of them scoring the old model.
FGW_FEATURES = load_features(args.input_dir)
AML_TIMEPOINTS = _VOCAB["aml_timepoints"]
# Fail loudly on a label the vocabulary does not cover. Without this an unrecognised
# timepoint is not an error -- it is a sample that quietly stops being AML and stops
# being healthy, and therefore stops existing for every test below.
assert_index_covered(index, _VOCAB)

print(f"[0] loaded {len(index)} samples | alpha={args.alpha} | test={args.test}")

## -- build per-sample C (7x7 directed), F (7x3), p (7); reindex by NAME (order-safe) ----
def build_one(ds, smp):
    e = edges[(edges["dataset"] == ds) & (edges["sample"] == smp)]
    C = (e.pivot(index="sender_bin", columns="receiver_bin", values="C")
           .reindex(index=FGW_NODES, columns=FGW_NODES))
    C = C.to_numpy(dtype=float)
    if np.isnan(C).any(): C = np.nan_to_num(C, nan=1.0)          # defensive: any unfilled pair -> max dist

    nd = nodes[(nodes["dataset"] == ds) & (nodes["sample"] == smp)].set_index("hierarchy_bin").reindex(FGW_NODES)
    F = nd[FGW_FEATURES].to_numpy(dtype=float)
    F = np.nan_to_num(F, nan=0.0)
    p = nd["mass"].to_numpy(dtype=float)
    p = np.nan_to_num(p, nan=1e-6); p = p / p.sum()
    return C, F, p

# group membership
idx = index.copy()
idx["is_healthy"] = idx["timepoint"].astype(str) == "Healthy"
idx["is_aml"]     = idx["timepoint"].astype(str).isin(AML_TIMEPOINTS)
idx["bary_ok"]    = (~idx["sparse_flag"].fillna(False)) if EXCLUDE_SPARSE_FROM_BARY else True

def group_inputs(mask):
    sub = idx[mask]
    if args.test > 0: sub = sub.head(args.test)
    Cs, Fs, ps = [], [], []
    for _, r in sub.iterrows():
        C, F, p = build_one(r.dataset, r["sample"]); Cs.append(C); Fs.append(F); ps.append(p)
    return Cs, Fs, ps, sub

## -- barycenter (balanced, directed, aligned-mean init) ----
def barycenter(Cs, Fs, ps, name):
    n = len(FGW_NODES); m = len(Cs)
    lambdas = [1.0/m]*m
    initC = np.mean(np.stack(Cs), axis=0)          # aligned mean init (nodes already label-aligned)
    initX = np.mean(np.stack(Fs), axis=0)
    pbar  = np.ones(n)/n
    out = ot.gromov.fgw_barycenters(n, Fs, Cs, ps, lambdas=lambdas, alpha=args.alpha,
                                    loss_fun="square_loss", symmetric=False, max_iter=args.max_iter,
                                    p=pbar, init_C=initC, init_X=initX, random_state=SEED, log=True)
    X, C = np.asarray(out[0]), np.asarray(out[1])
    print(f"[1] barycenter {name}: built from {m} graphs | NaN(C)={bool(np.isnan(C).any())} "
          f"NaN(F)={bool(np.isnan(X).any())} | asym=max|C-C.T|={float(np.abs(C-C.T).max()):.4f}")
    return C, X, pbar

print("[1] building barycenters (non-sparse samples)")
Ch, Fh, ph, subh = group_inputs(idx.is_healthy & idx.bary_ok)
Ca, Fa, pa, suba = group_inputs(idx.is_aml     & idx.bary_ok)
print(f"[1] healthy barycenter n={len(Ch)} | AML barycenter n={len(Ca)}")
C_healthy, F_healthy, p_bar = barycenter(Ch, Fh, ph, "healthy")
C_aml,     F_aml,     _     = barycenter(Ca, Fa, pa, "aml")

## -- (2) directionality check on a primitive->immune edge ----
i, j = FGW_NODES.index("LMPP_GMP"), FGW_NODES.index("B_Plasma")
print(f"[2] directionality LMPP_GMP<->B_Plasma  (lower C = stronger):")
print(f"    B_healthy: fwd C[LMPP_GMP,B_Plasma]={C_healthy[i,j]:.3f}  rev={C_healthy[j,i]:.3f}")
print(f"    B_AML:     fwd C[LMPP_GMP,B_Plasma]={C_aml[i,j]:.3f}  rev={C_aml[j,i]:.3f}")

## -- pairwise HDS / ATS for ALL samples ----
def fgw_dist(C, F, p, Cb, Fb, pb):
    M = ot.dist(F, Fb); M = M / (M.max() + 1e-9)
    return float(ot.gromov.fused_gromov_wasserstein2(M, C, Cb, p, pb, loss_fun="square_loss",
                                                     alpha=args.alpha, symmetric=False))

def fgw_coupling(C, F, p, Cb, Fb, pb):
    """FGW distance AND the optimal coupling T that produced it."""
    M = ot.dist(F, Fb); M = M / (M.max() + 1e-9)
    T, log = ot.gromov.fused_gromov_wasserstein(M, C, Cb, p, pb, loss_fun="square_loss",
                                                alpha=args.alpha, symmetric=False, log=True)
    return T, float(log.get("fgw_dist", np.nan))


def pullback(C, T, pb):
    """Push a sample's cost matrix onto the barycenter's node vocabulary through the coupling.

    C_tilde[j,l] = sum_{i,k} T[i,j] T[k,l] C[i,k] / (pb[j] pb[l])

    This is what makes the difference vector below live in a FIXED coordinate system: without it,
    d_p for two samples would be expressed in two different node orderings and could not be
    projected onto a common axis.
    """
    den = np.outer(pb, pb)
    den[den <= 0] = 1e-12
    return (T.T @ C @ T) / den


# M6 (BLUEPRINT_v1.1_PATCH.md section 7). HDS as a single scalar answers only "how far", and the
# patch is explicit about the cost: "没有方向分量，B 和 C 无法与 A 区分" -- a sample that has moved
# FARTHER along the healthy->AML axis and one that has moved OFF that axis get the same number.
#
# DEVIATION, deliberate: the patch writes vec_triangle (upper triangle). Our C is DIRECTED and 55%
# of node pairs disagree between the two directions, so an upper-triangle vectorisation would throw
# away exactly the asymmetry M6 exists to measure. We vectorise the full 7x7 = 49 entries.
#
# Also: HDS_norm below is ||d_p||, which the patch calls "equivalent to the original GW distance".
# It is NOT identical -- the FGW distance also carries the feature term and is a different
# functional -- so both are reported and neither is presented as the other.
U = (C_aml - C_healthy).reshape(-1)
u_norm = float(np.linalg.norm(U))
if u_norm <= 0:
    raise SystemExit("the two barycenters have identical structure; there is no healthy->AML axis "
                     "to project onto, and the direction components below would be meaningless.")
U = U / u_norm
print(f"[3] healthy->AML structural axis: ||C_aml - C_healthy|| = {u_norm:.4f} over 49 directed edges")

print("[3] scoring all samples (HDS to B_healthy, ATS to B_AML, + direction components)")
rows = []
scan = idx.head(args.test*2) if args.test > 0 else idx
for _, r in scan.iterrows():
    C, F, p = build_one(r.dataset, r["sample"])
    T, hds = fgw_coupling(C, F, p, C_healthy, F_healthy, p_bar)
    ats = fgw_dist(C, F, p, C_aml, F_aml, p_bar)
    d = (pullback(C, T, p_bar) - C_healthy).reshape(-1)
    mag = float(np.linalg.norm(d))
    par = float(d @ U)
    perp = float(np.sqrt(max(0.0, mag**2 - par**2)))
    cos = par / mag if mag > 0 else np.nan
    rows.append(dict(dataset=r.dataset, sample=r["sample"], timepoint=r.timepoint,
                     sparse_flag=bool(r.sparse_flag) if pd.notna(r.sparse_flag) else None,
                     HDS=hds, ATS=ats,
                     HDS_norm=mag, HDS_par=par, HDS_perp=perp, cos_theta=cos,
                     alpha=args.alpha, mass_mode=r.get("mass_mode","")))
scores = pd.DataFrame(rows)

# NON-VACUITY: if the pullback collapsed, every d_p is the same vector and the decomposition is a
# constant dressed as four readouts.
if scores["HDS_norm"].std() < 1e-9:
    raise SystemExit("HDS_norm is constant across samples -- the pullback produced the same matrix "
                     "for every sample. The coupling is degenerate; do not report the components.")
print(f"[3] direction components: cos_theta ranges {scores.cos_theta.min():+.3f} to "
      f"{scores.cos_theta.max():+.3f}, HDS_perp / HDS_norm median "
      f"{(scores.HDS_perp / scores.HDS_norm).median():.3f}")

## -- (3) sanity: mean HDS by timepoint (healthy should be low; disease -> higher = H3 preview) ----
print("[3] mean HDS / ATS by timepoint:")
print(scores.groupby("timepoint")[["HDS","ATS"]].mean().round(4).to_string())

## -- save ----
np.savez(os.path.join(args.input_dir, "barycenters.npz"),
         C_healthy=C_healthy, F_healthy=F_healthy, C_aml=C_aml, F_aml=F_aml,
         p=p_bar, nodes=np.array(FGW_NODES), features=np.array(FGW_FEATURES), alpha=args.alpha)
scores.to_csv(os.path.join(args.input_dir, "patient_scores.csv"), index=False)
print(f"[done] wrote barycenters.npz + patient_scores.csv ({len(scores)} samples) to {args.input_dir}")