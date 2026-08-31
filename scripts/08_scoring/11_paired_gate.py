#!/usr/bin/env python
# 11_paired_gate.py ----
# Phase 7 (stage 08_scoring). The paired-sample gate on the FGW transport score.
# Everything below is fixed in scripts/08_scoring/PREREGISTRATION_paired_gate.md,
# written and committed BEFORE this script produced a number.
#
# WHY. Every topological result so far is negative, and 10_planted_effect_power showed the
# alpha=1 score has no power against the class of difference this pipeline makes. Two
# explanations remain and have never been separated:
#   (1) the graph REPRESENTATION carries no usable signal, or
#   (2) the edge-to-cost TRANSFORM destroys it, because 46.9% of every cost matrix is a
#       placeholder whose value is set by how sparse the sample is (C = 1 - (k+1)/98).
# Paired samples are the cleanest material for (1). Four distance arms address (2). Crossed.
#
# THREE GATES, in order:
#   GATE 1  identity     is a patient's own second graph closer than other patients' graphs?
#                        Direct FGW2(graph_i, graph_j). NO healthy barycentre, no label, no
#                        covariate. Candidate pool = same dataset, OTHER patients only --
#                        without that restriction the test is won by batch.
#   GATE 2  positive     diagnosis -> post-treatment, 11 pairs. The marrow goes from
#           control      blast-packed to regenerating, so a difference is near-certain.
#                        Sign pre-specified: HDS_treatment < HDS_diagnosis.
#   GATE 3  the question diagnosis -> relapse, 11 pairs. NOT a gate. Reported, never used
#                        to stop anything.
# Pairing does NOT remove depth (measured within-pair ratios 0.23x to 3.54x), so GATE 2/3
# regress the within-pair HDS difference on the within-pair log10-depth difference and test
# the INTERCEPT.
#
# INPUT  : <root>/06_distance/edge_distance.csv
#          <root>/07_fgw/fgw_{nodes,input_index}_long.csv , fgw_vocab.json
#          <root>/01_preprocess/00_curated_manifest.csv    (patient_id -> the pairing)
#          <root>/01_preprocess/03_qc_report__ALL.csv      (med_ncount_final -> depth)
#          <root>/07_fgw/patient_scores.csv                (self-check 2 only)
# OUTPUT : <root>/08_scoring/paired_gate.csv        one row per (arm, alpha, gate, test)
#          <root>/08_scoring/paired_gate_ranks.csv  per patient per arm per alpha
# Usage  : python scripts/08_scoring/11_paired_gate.py [--arms rank,mask,const,logs]
import argparse, os, sys, warnings
import numpy as np
import pandas as pd
import ot
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config'))
from fgw_vocab import load_vocab, load_features, assert_index_covered
from distance_variants import weights_to_C, LOG_EPS

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
REL_TP  = {"Relapse", "Relapse2"}
TRT_TP  = {"On_treatment", "Post_induction", "Post_consolidation",
           "Post_treatment_unspecified", "Refractory"}

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--arms", default="rank,mask,const,logs")
ap.add_argument("--alphas_gate1", default="0,0.5,1")
ap.add_argument("--alphas_gate23", default="0.5,1")
ap.add_argument("--max_iter", type=int, default=1000)
ap.add_argument("--skip_selfcheck", action="store_true",
                help="DIAGNOSTIC ONLY. The registration says a failed self-check voids the run.")
args = ap.parse_args()
ARMS = [a.strip() for a in args.arms.split(",")]
D_FGW = os.path.join(args.root, "07_fgw"); D_OUT = os.path.join(args.root, "08_scoring")
D_DST = os.path.join(args.root, "06_distance"); D_PRE = os.path.join(args.root, "01_preprocess")

## -- load ----
edges = pd.read_csv(os.path.join(D_DST, "edge_distance.csv"))
nodes = pd.read_csv(os.path.join(D_FGW, "fgw_nodes_long.csv"))
index = pd.read_csv(os.path.join(D_FGW, "fgw_input_index.csv"))
VOCAB = load_vocab(D_FGW); FEATS = load_features(D_FGW)
assert_index_covered(index, VOCAB)
AML_TP = set(VOCAB["aml_timepoints"])

man = (pd.read_csv(os.path.join(D_PRE, "00_curated_manifest.csv"))
       [["dataset", "sample", "patient_id"]].drop_duplicates(["dataset", "sample"]))
qc = (pd.read_csv(os.path.join(D_PRE, "03_qc_report__ALL.csv")).rename(columns={"Sample": "sample"})
      [["dataset", "sample", "med_ncount_final"]].drop_duplicates(["dataset", "sample"]))
idx = index.merge(man, on=["dataset", "sample"], how="left").merge(qc, on=["dataset", "sample"], how="left")
if idx.patient_id.isna().any():
    raise SystemExit("%d graph sample(s) carry no patient_id; pairing is undefined for them."
                     % int(idx.patient_id.isna().sum()))
if idx.med_ncount_final.isna().any():
    raise SystemExit("%d graph sample(s) carry no depth; GATE 2/3 cannot run."
                     % int(idx.med_ncount_final.isna().sum()))
idx["logdepth"] = np.log10(idx.med_ncount_final.clip(lower=1))
KEY = [(r.dataset, r["sample"]) for _, r in idx.iterrows()]
print("[0] %d graphs | %d patients | arms %s" % (len(idx), idx.patient_id.nunique(), ARMS))

## -- the four cost transforms. Same w in, four C out. ----
GRID = pd.MultiIndex.from_product([FGW_NODES, FGW_NODES], names=["sender_bin", "receiver_bin"])

def build_C(sub, arm):
    """sub: the 49 rows of one sample, already on the fixed grid order. Returns 7x7.
    The transforms themselves live in scripts/config/distance_variants.py so that this
    script, 06_alpha_sweep and 10_planted_effect_power cannot drift apart."""
    c = weights_to_C(sub["weight_probsum"].to_numpy(float), arm,
                     sub["n_lr_sig"].to_numpy(float))
    return c.reshape(len(FGW_NODES), len(FGW_NODES))

## -- per-sample inputs ----
E = {k: g.set_index(["sender_bin", "receiver_bin"]).reindex(GRID).reset_index()
     for k, g in edges.groupby(["dataset", "sample"])}
missing = [k for k in KEY if k not in E]
if missing:
    raise SystemExit("%d graph sample(s) have no rows in edge_distance.csv: %s" % (len(missing), missing[:5]))

ND = {}
for k, g in nodes.groupby(["dataset", "sample"]):
    nd = g.set_index("hierarchy_bin").reindex(FGW_NODES)
    F = np.nan_to_num(nd[FEATS].to_numpy(float), nan=0.0)
    p = np.nan_to_num(nd["mass"].to_numpy(float), nan=1e-6); p = p / p.sum()
    ND[k] = (F, p)

CBANK = {arm: {k: build_C(E[k], arm) for k in KEY} for arm in ARMS}

## -- SELF-CHECK 1: the rank arm must reproduce the stored C bit for bit ----
if "rank" in ARMS and not args.skip_selfcheck:
    dev = 0.0
    for k in KEY:
        stored = E[k]["C"].to_numpy(float).reshape(7, 7)
        dev = max(dev, float(np.abs(CBANK["rank"][k] - stored).max()))
    print("[1] SELF-CHECK 1  rank arm vs stored C: max |diff| = %.3e" % dev)
    if dev > 1e-12:
        raise SystemExit("SELF-CHECK 1 FAILED. The re-implemented transform is not the production one.")

## -- FGW ----
def fgw2(Ca, Fa, pa, Cb, Fb, pb, alpha):
    M = ot.dist(Fa, Fb); mx = M.max()
    M = M / (mx + 1e-9) if mx > 0 else M
    with pot_quiet():
        return float(ot.gromov.fused_gromov_wasserstein2(
            M, Ca, Cb, pa, pb, loss_fun="square_loss", alpha=alpha, symmetric=False))

def barycenter(keys, arm, alpha):
    Cs = [CBANK[arm][k] for k in keys]; Fs = [ND[k][0] for k in keys]; ps = [ND[k][1] for k in keys]
    m = len(Cs); n = len(FGW_NODES)
    with pot_quiet():
        out = ot.gromov.fgw_barycenters(
            n, Fs, Cs, ps, lambdas=[1.0/m]*m, alpha=alpha, loss_fun="square_loss",
            symmetric=False, max_iter=args.max_iter, p=np.ones(n)/n,
            init_C=np.mean(np.stack(Cs), 0), init_X=np.mean(np.stack(Fs), 0),
            random_state=SEED, log=True)
    return np.asarray(out[1]), np.asarray(out[0]), np.ones(n)/n

## -- SELF-CHECK 2: the rank arm at alpha 0.5 must reproduce patient_scores.csv ----
def hds_all(arm, alpha):
    heal = [k for k, t in zip(KEY, idx.timepoint) if t == "Healthy"]
    sparse = dict(zip(KEY, idx.sparse_flag.fillna(False)))
    bar_keys = [k for k in heal if not sparse[k]]
    Cb, Fb, pb = barycenter(bar_keys, arm, alpha)
    return {k: fgw2(CBANK[arm][k], ND[k][0], ND[k][1], Cb, Fb, pb, alpha) for k in KEY}

if "rank" in ARMS and not args.skip_selfcheck:
    ref = pd.read_csv(os.path.join(D_FGW, "patient_scores.csv"))
    ref = dict(zip(zip(ref.dataset, ref["sample"]), ref.HDS))
    got = hds_all("rank", 0.5)
    dev = max(abs(got[k] - ref[k]) for k in KEY if k in ref)
    print("[1] SELF-CHECK 2  rank arm HDS vs patient_scores.csv: max |diff| = %.3e" % dev)
    if dev > 1e-6:
        raise SystemExit("SELF-CHECK 2 FAILED. The re-implemented barycentre/scoring is not the production one.")

## -- the pairs ----
pairs = []
for (ds, pt), g in idx.groupby(["dataset", "patient_id"]):
    tps = dict(zip(g.timepoint, zip(g.dataset, g["sample"])))
    if "Diagnosis" not in tps: continue
    for kind, want in [("Dx_to_Relapse", REL_TP), ("Dx_to_Treatment", TRT_TP)]:
        hit = [t for t in g.timepoint if t in want]
        if hit:
            pairs.append(dict(kind=kind, dataset=ds, patient=pt,
                              a=tps["Diagnosis"], b=tps[sorted(hit)[0]], tp_b=sorted(hit)[0]))
P = pd.DataFrame(pairs)
print("[2] pairs: %s" % P.kind.value_counts().to_dict())

## -- GATE 1: identity, same dataset, OTHER patients only ----
rank_rows, res = [], []
for arm in ARMS:
    for alpha in [float(x) for x in args.alphas_gate1.split(",")]:
        for kind, grp in P.groupby("kind"):
            pcts, top1 = [], 0
            for _, r in grp.iterrows():
                ka, kb = r.a, r.b
                pool = [k for k, ds, pid, tp in zip(KEY, idx.dataset, idx.patient_id, idx.timepoint)
                        if ds == r.dataset and pid != r.patient and tp in AML_TP]
                if len(pool) < 3:
                    print("      [skip] %s/%s only %d same-dataset distractors" % (r.dataset, r.patient, len(pool)))
                    continue
                d_true = fgw2(CBANK[arm][ka], *ND[ka], CBANK[arm][kb], *ND[kb], alpha)
                d_pool = np.array([fgw2(CBANK[arm][ka], *ND[ka], CBANK[arm][g], *ND[g], alpha) for g in pool])
                pct = float((d_pool < d_true).sum()) / len(d_pool)
                pcts.append(pct); top1 += int(pct == 0.0)
                rank_rows.append(dict(arm=arm, alpha=alpha, kind=kind, dataset=r.dataset,
                                      patient=r.patient, n_pool=len(pool), d_true=d_true,
                                      d_pool_median=float(np.median(d_pool)), pct=pct))
            if len(pcts) < 5: continue
            w = stats.wilcoxon(np.array(pcts) - 0.5, alternative="less")
            res.append(dict(gate="GATE1_identity", arm=arm, alpha=alpha, kind=kind, n=len(pcts),
                            stat="median percentile of the true partner", value=float(np.median(pcts)),
                            p=float(w.pvalue), extra="top1=%d/%d" % (top1, len(pcts))))
            print("[3] GATE1 %-5s a=%-3g %-16s n=%2d median pct %.3f  p=%.4f  top1 %d/%d"
                  % (arm, alpha, kind, len(pcts), np.median(pcts), w.pvalue, top1, len(pcts)))

## -- GATE 2 and 3: paired HDS difference, depth as covariate, test the intercept ----
dep = dict(zip(KEY, idx.logdepth))
for arm in ARMS:
    for alpha in [float(x) for x in args.alphas_gate23.split(",")]:
        H = hds_all(arm, alpha)
        for kind, grp in P.groupby("kind"):
            gate = "GATE2_positive_control" if kind == "Dx_to_Treatment" else "GATE3_question"
            side = "less" if kind == "Dx_to_Treatment" else "two-sided"
            dh = np.array([H[r.b] - H[r.a] for _, r in grp.iterrows()])
            dd = np.array([dep[r.b] - dep[r.a] for _, r in grp.iterrows()])
            X = np.column_stack([np.ones(len(dh)), dd])
            beta, *_ = np.linalg.lstsq(X, dh, rcond=None)
            resid = dh - X @ beta; dof = len(dh) - X.shape[1]
            se = np.sqrt((resid @ resid / dof) * np.linalg.inv(X.T @ X)[0, 0])
            t = beta[0] / se
            p = stats.t.cdf(t, dof) if side == "less" else 2 * stats.t.sf(abs(t), dof)
            wp = stats.wilcoxon(dh, alternative=side).pvalue
            res.append(dict(gate=gate, arm=arm, alpha=alpha, kind=kind, n=len(dh),
                            stat="intercept of dHDS ~ dLogDepth", value=float(beta[0]), p=float(p),
                            extra="slope=%+.4f wilcoxon_nocovar_p=%.4f median_dHDS=%+.4f"
                                  % (beta[1], wp, float(np.median(dh)))))
            print("[4] %s %-5s a=%-3g n=%2d intercept %+.5f  p=%.4f  (raw Wilcoxon p=%.4f)"
                  % (gate.split("_")[0], arm, alpha, len(dh), beta[0], p, wp))

os.makedirs(D_OUT, exist_ok=True)
pd.DataFrame(res).to_csv(os.path.join(D_OUT, "paired_gate.csv"), index=False)
pd.DataFrame(rank_rows).to_csv(os.path.join(D_OUT, "paired_gate_ranks.csv"), index=False)
print("\n[done] wrote paired_gate.csv (%d rows) and paired_gate_ranks.csv (%d rows)"
      % (len(res), len(rank_rows)))
