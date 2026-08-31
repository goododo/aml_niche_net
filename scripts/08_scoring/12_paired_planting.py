#!/usr/bin/env python
# 12_paired_planting.py ----
# Phase 7 (stage 08_scoring). The planted-effect harness, run INSIDE the paired samples.
#
# WHY A PAIRED VERSION. 10_planted_effect_power.py plants into AML samples and tests the
# cross-sectional AML-vs-healthy contrast, where donor, batch and protocol are all noise.
# Here the effect is planted into the SECOND timepoint of a patient who supplied the first,
# so all three are controlled by construction. The same statistic as GATE 2 in
# 11_paired_gate.py is then applied, which makes the two directly comparable:
#
#   GATE 2 measured the pipeline against a REAL effect that is near-certain by biology
#   (diagnosis to post-treatment). It passed at alpha 0.5 in all four arms and failed at
#   alpha 1 in all four. This script answers the follow-up question that leaves open:
#   HOW LARGE would an effect have to be for alpha = 1 to see it?
#
# The healthy barycentre is never planted into, so it is computed once per (arm, alpha) and
# reused across every delta. That is what makes a 4-arm sweep cheap here.
#
# INPUT  : <root>/06_distance/edge_distance.csv
#          <root>/07_fgw/fgw_{nodes,input_index}_long.csv , fgw_vocab.json
#          <root>/01_preprocess/00_curated_manifest.csv , 03_qc_report__ALL.csv
# OUTPUT : <root>/08_scoring/paired_planting.csv
# Usage  : python scripts/08_scoring/12_paired_planting.py [--arms rank,mask,const,logs]
import argparse, os, sys, warnings
import numpy as np
import pandas as pd
import ot
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config'))
from fgw_vocab import load_vocab, load_features, assert_index_covered
from distance_variants import weights_to_C, ARMS as DV_ARMS

from contextlib import contextmanager
@contextmanager
def pot_quiet():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered")
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered")
        yield

FGW_NODES = ["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
SEED = 491638
REL_TP = {"Relapse", "Relapse2"}
TRT_TP = {"On_treatment", "Post_induction", "Post_consolidation",
          "Post_treatment_unspecified", "Refractory"}

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--arms", default="rank,mask,const,logs")
ap.add_argument("--alphas", default="0.5,1")
ap.add_argument("--deltas", default="0,0.25,0.5,1,2,4")
ap.add_argument("--ks", default="1,3,6")
ap.add_argument("--max_iter", type=int, default=1000)
args = ap.parse_args()
ARMS   = [a.strip() for a in args.arms.split(",")]
ALPHAS = [float(x) for x in args.alphas.split(",")]
DELTAS = [float(x) for x in args.deltas.split(",")]
KS     = [int(x) for x in args.ks.split(",")]

D_FGW = os.path.join(args.root, "07_fgw"); D_OUT = os.path.join(args.root, "08_scoring")
D_DST = os.path.join(args.root, "06_distance"); D_PRE = os.path.join(args.root, "01_preprocess")

edges = pd.read_csv(os.path.join(D_DST, "edge_distance.csv"))
nodes = pd.read_csv(os.path.join(D_FGW, "fgw_nodes_long.csv"))
index = pd.read_csv(os.path.join(D_FGW, "fgw_input_index.csv"))
VOCAB = load_vocab(D_FGW); FEATS = load_features(D_FGW)
assert_index_covered(index, VOCAB)

man = (pd.read_csv(os.path.join(D_PRE, "00_curated_manifest.csv"))
       [["dataset","sample","patient_id"]].drop_duplicates(["dataset","sample"]))
qc = (pd.read_csv(os.path.join(D_PRE, "03_qc_report__ALL.csv")).rename(columns={"Sample":"sample"})
      [["dataset","sample","med_ncount_final"]].drop_duplicates(["dataset","sample"]))
idx = index.merge(man, on=["dataset","sample"]).merge(qc, on=["dataset","sample"])
idx["logdepth"] = np.log10(idx.med_ncount_final.clip(lower=1))
KEY = [(r.dataset, r["sample"]) for _, r in idx.iterrows()]

GRID = pd.MultiIndex.from_product([FGW_NODES, FGW_NODES], names=["sender_bin","receiver_bin"])
EDGE_NAMES = ["%s->%s" % (a, b) for a in FGW_NODES for b in FGW_NODES]
W, NLR = {}, {}
for k, g in edges.groupby(["dataset","sample"]):
    sub = g.set_index(["sender_bin","receiver_bin"]).reindex(GRID).reset_index()
    W[k]   = sub["weight_probsum"].to_numpy(float)
    NLR[k] = sub["n_lr_sig"].to_numpy(float)
ND = {}
for k, g in nodes.groupby(["dataset","sample"]):
    nd = g.set_index("hierarchy_bin").reindex(FGW_NODES)
    F = np.nan_to_num(nd[FEATS].to_numpy(float), nan=0.0)
    p = np.nan_to_num(nd["mass"].to_numpy(float), nan=1e-6); p = p / p.sum()
    ND[k] = (F, p)

## -- SELF-CHECK: the rank arm must reproduce the stored C ----
dev = max(float(np.abs(weights_to_C(W[k], "rank", NLR[k])
                       - edges[(edges.dataset==k[0]) & (edges["sample"]==k[1])]
                         .set_index(["sender_bin","receiver_bin"]).reindex(GRID)["C"].to_numpy(float)).max())
          for k in KEY)
print("[0] SELF-CHECK rank arm vs stored C: max |diff| = %.3e" % dev)
if dev > 1e-12:
    raise SystemExit("the shared transform no longer reproduces 06_distance/01; aborting.")

## -- the pairs, same rule as 11_paired_gate.py ----
pairs = []
for (ds, pt), g in idx.groupby(["dataset","patient_id"]):
    tps = dict(zip(g.timepoint, zip(g.dataset, g["sample"])))
    if "Diagnosis" not in tps: continue
    for kind, want in [("Dx_to_Relapse", REL_TP), ("Dx_to_Treatment", TRT_TP)]:
        hit = sorted([t for t in g.timepoint if t in want])
        if hit:
            pairs.append(dict(kind=kind, a=tps["Diagnosis"], b=tps[hit[0]]))
P = pd.DataFrame(pairs)
print("[1] pairs: %s" % P.kind.value_counts().to_dict())

## -- which edges to plant on ----
# The SAME rule and the SAME pinned transform as 10_planted_effect_power.py: rank edges by how
# little of their cost variance the AML/healthy label explains WITHIN dataset, on the production
# rank C, over the both-arm non-sparse samples. Pinning to the rank transform means every arm
# plants on the same edges, so the arms differ only in the transform.
AML_SET = set(VOCAB["aml_timepoints"])
_g = idx.copy()
_g["grp"] = np.where(_g.timepoint == "Healthy", "Healthy",
                     np.where(_g.timepoint.isin(AML_SET), "AML", "other"))
_g["is_aml"] = (_g.grp == "AML").astype(float)
_both = [ds for ds, gg in _g.groupby("dataset") if gg.is_aml.nunique() == 2]
_sel = (_g.dataset.isin(_both) & ~_g.sparse_flag.fillna(False)).to_numpy()
C0_rank = np.vstack([weights_to_C(W[k], "rank", NLR[k]) for k in KEY])
lab_var = {}
for j, e in enumerate(EDGE_NAMES):
    y = C0_rank[_sel, j]; gd = _g.dataset.to_numpy()[_sel]; a = _g.is_aml.to_numpy()[_sel]
    tot = float(((y - y.mean())**2).sum())
    if tot == 0: lab_var[e] = 0.0; continue
    ss = 0.0
    for ds in np.unique(gd):
        m = gd == ds
        if len(np.unique(a[m])) < 2: continue
        mu = y[m].mean()
        for l in (0.0, 1.0):
            mm = m & (a == l)
            if mm.sum(): ss += mm.sum() * (y[mm].mean() - mu)**2
    lab_var[e] = ss / tot
ORDER = sorted(EDGE_NAMES, key=lambda e: lab_var[e])
print("[1] planting on the lowest-label-variance edges (pinned to the rank transform); cleanest 6:")
print("    " + ", ".join("%s (%.1f%%)" % (e, 100*lab_var[e]) for e in ORDER[:6]))
EXPECTED6 = ["B_Plasma->Erythroid","HSC_MPP->HSC_MPP","Megakaryocyte->HSC_MPP",
             "T_NK->Erythroid","B_Plasma->Mono_DC","Megakaryocyte->T_NK"]
if ORDER[:6] != EXPECTED6:
    raise SystemExit("planted-edge selection has drifted from 10_planted_effect_power.py:\n"
                     "  got      %s\n  expected %s" % (ORDER[:6], EXPECTED6))
print("    [PASS] identical to the cross-sectional harness, so the two are comparable")

def fgw2(Ca, Fa, pa, Cb, Fb, pb, alpha):
    M = ot.dist(Fa, Fb); mx = M.max(); M = M/(mx+1e-9) if mx > 0 else M
    with pot_quiet():
        return float(ot.gromov.fused_gromov_wasserstein2(M, Ca, Cb, pa, pb,
                     loss_fun="square_loss", alpha=alpha, symmetric=False))

def C_of(k, arm, w=None):
    ww = W[k] if w is None else w
    return weights_to_C(ww, arm, NLR[k]).reshape(len(FGW_NODES), len(FGW_NODES))

dep = dict(zip(KEY, idx.logdepth))
SCALE = {k: (np.median(W[k][W[k] > 0]) if (W[k] > 0).any() else 0.0) for k in KEY}
rows = []
for arm in ARMS:
    for alpha in ALPHAS:
        # the barycentre is built from healthy samples, which are never planted into,
        # so it is computed once here and reused for every delta below.
        heal = [k for k, t, s in zip(KEY, idx.timepoint, idx.sparse_flag.fillna(False))
                if t == "Healthy" and not s]
        Cs = [C_of(k, arm) for k in heal]; Fs = [ND[k][0] for k in heal]; ps = [ND[k][1] for k in heal]
        m = len(Cs); n = len(FGW_NODES)
        with pot_quiet():
            out = ot.gromov.fgw_barycenters(n, Fs, Cs, ps, lambdas=[1.0/m]*m, alpha=alpha,
                  loss_fun="square_loss", symmetric=False, max_iter=args.max_iter,
                  p=np.ones(n)/n, init_C=np.mean(np.stack(Cs),0), init_X=np.mean(np.stack(Fs),0),
                  random_state=SEED, log=True)
        Cb, Fb, pb = np.asarray(out[1]), np.asarray(out[0]), np.ones(n)/n
        hds_a = {k: fgw2(C_of(k, arm), *ND[k], Cb, Fb, pb, alpha) for k in set(P.a)}
        for k_edges in KS:
            ji = [EDGE_NAMES.index(e) for e in ORDER[:k_edges]]
            for delta in DELTAS:
                dh, dd, moved = [], [], []
                for _, r in P.iterrows():
                    kb = r.b
                    wb = W[kb].copy()
                    if delta > 0: wb[ji] += delta * SCALE[kb]
                    Cm = C_of(kb, arm, wb)
                    moved.append(float((np.abs(Cm.ravel()[ji] - C_of(kb, arm).ravel()[ji]) > 0).mean()))
                    dh.append(fgw2(Cm, *ND[kb], Cb, Fb, pb, alpha) - hds_a[r.a])
                    dd.append(dep[kb] - dep[r.a])
                for kind in P.kind.unique():
                    sel = (P.kind == kind).to_numpy()
                    y = np.array(dh)[sel]; x = np.array(dd)[sel]
                    X = np.column_stack([np.ones(len(y)), x])
                    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                    res = y - X @ beta; dof = len(y) - 2
                    se = np.sqrt((res @ res / dof) * np.linalg.inv(X.T @ X)[0,0])
                    t = beta[0]/se
                    rows.append(dict(arm=arm, alpha=alpha, kind=kind, k=k_edges, delta=delta,
                                     n=int(sel.sum()), frac_targets_moved=float(np.mean(np.array(moved)[sel])),
                                     intercept=float(beta[0]), p_two_sided=float(2*stats.t.sf(abs(t), dof)),
                                     wilcoxon_p=float(stats.wilcoxon(y).pvalue) if np.any(y != 0) else 1.0))
            print("[2] %-5s a=%-3g k=%d done" % (arm, alpha, k_edges), flush=True)

R = pd.DataFrame(rows)
os.makedirs(D_OUT, exist_ok=True)
R.to_csv(os.path.join(D_OUT, "paired_planting.csv"), index=False)
# GATE 2 (Dx_to_Treatment) is registered one-sided: HDS(treatment) < HDS(diagnosis), i.e.
# intercept < 0. GATE 3 (Dx_to_Relapse) has no pre-specified sign, so it stays two-sided.
# Reporting the two-sided p for both understated GATE 2 by exactly a factor of two.
R["p_registered"] = np.where(
    R.kind == "Dx_to_Treatment",
    np.where(R.intercept < 0, R.p_two_sided / 2, 1 - R.p_two_sided / 2),
    R.p_two_sided)
R.to_csv(os.path.join(D_OUT, "paired_planting.csv"), index=False)
print("\n[3] smallest delta reaching the REGISTERED test at p < 0.05")
print("    (delta=0 means the real effect is already detected without any planting)")
for (arm, alpha, kind, k_edges), g in R.groupby(["arm","alpha","kind","k"]):
    hit = g[g.p_registered < 0.05].sort_values("delta")
    print("    %-5s a=%-3g %-16s k=%d : %s" % (arm, alpha, kind, k_edges,
          ("delta=%g" % hit.delta.iloc[0]) if len(hit) else "never"))
print("\n[done] wrote paired_planting.csv (%d rows)" % len(R))
