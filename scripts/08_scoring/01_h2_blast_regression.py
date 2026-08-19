#!/usr/bin/env python
# 01_h2_blast_regression.py ----
# Phase 7 (stage 08_scoring). H2 core test vs the #1 reviewer attack (R5): "HDS just proxies blast
# burden / composition." Fit HDS ~ blast_proxy + is_aml and test whether is_aml stays significant
# AFTER controlling for composition. If yes -> network topology carries disease information BEYOND cell
# composition -> H2 supported.
#
# Two correctness safeguards:
#   (1) LOO-clean HDS: B_healthy was BUILT from healthy samples, so scoring those same samples is
#       circular (they sit near a barycenter they defined). For each healthy sample, HDS is recomputed
#       against a barycenter that EXCLUDES it (leave-one-out). AML samples are never in B_healthy -> their
#       existing HDS is already clean.
#   (2) blast_proxy = (HSC_MPP + LMPP_GMP + Mono_DC cells) / total cells -- a COMPOSITION proxy that does
#       NOT depend on inferCNV (avoids the malignancy-call pollution). frac_malignant used only as a
#       sensitivity covariate.
# Significance of the is_aml partial effect via LABEL PERMUTATION (dependency-free), not a t-table.
#
# INPUT  : <root>/07_fgw/fgw_{nodes,edges,input_index}_long.csv , <root>/07_fgw/patient_scores.csv (AML HDS)
#          <root>/05_ccc/ccc_node_features.csv (n_malignant/n_evaluable for the frac_malignant sensitivity)
# OUTPUT : <root>/08_scoring/h2_blast_regression.csv (per-sample HDS_clean, blast_proxy, mal_frac, is_aml)
# Usage  : python scripts/08_scoring/01_h2_blast_regression.py [--root .../results/tables] [--alpha 0.5] [--n_perm 2000]
import argparse, os
import sys
import numpy as np
import pandas as pd
import ot

# The timepoint vocabulary is LOADED, not literal. A hard-coded set here silently deleted 34 of
# 214 samples (64% of the treated arm) after CANONICAL_TIMEPOINTS changed on 2026-08-04.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config'))
from fgw_vocab import load_vocab, assert_index_covered


FGW_NODES = ["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
FGW_FEATURES = ["frac_malignant","mean_stemness","n_cells"]
BLAST_BINS = ["HSC_MPP","LMPP_GMP","Mono_DC"]
AML_TP = None  # set from fgw_vocab.json below -- see scripts/config/fgw_vocab.py
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"; SEED = 491638

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--alpha", type=float, default=0.5)
ap.add_argument("--n_perm", type=int, default=2000)
args = ap.parse_args()
rng = np.random.default_rng(SEED)
D_FGW = os.path.join(args.root, "07_fgw"); D_CCC = os.path.join(args.root, "05_ccc"); D_OUT = os.path.join(args.root, "08_scoring")
_VOCAB = load_vocab(D_FGW)
AML_TP = _VOCAB["aml_timepoints"]

os.makedirs(D_OUT, exist_ok=True)

edges = pd.read_csv(os.path.join(D_FGW, "fgw_edges_long.csv"))
nodes = pd.read_csv(os.path.join(D_FGW, "fgw_nodes_long.csv"))
idx   = pd.read_csv(os.path.join(D_FGW, "fgw_input_index.csv"))
# Fail loudly on a label the vocabulary does not cover. Without this an unrecognised
# timepoint is not an error -- it is a sample that quietly stops being AML and stops
# being healthy, and therefore stops existing for every test below.
assert_index_covered(idx, _VOCAB)

sc    = pd.read_csv(os.path.join(D_FGW, "patient_scores.csv"))
nfraw = pd.read_csv(os.path.join(D_CCC, "ccc_node_features.csv"))

def key(df, ds, smp): return df[(df["dataset"]==ds) & (df["sample"]==smp)]
def build_CFp(ds, smp):
    e = key(edges, ds, smp)
    C = e.pivot(index="sender_bin", columns="receiver_bin", values="C").reindex(FGW_NODES, columns=FGW_NODES).to_numpy(float)
    C = np.nan_to_num(C, nan=1.0)
    nd = key(nodes, ds, smp).set_index("hierarchy_bin").reindex(FGW_NODES)
    F = np.nan_to_num(nd[FGW_FEATURES].to_numpy(float), nan=0.0)
    p = np.nan_to_num(nd["mass"].to_numpy(float), nan=1e-6); p = p/p.sum()
    return C, F, p
def barycenter(keys):
    Cs=[];Fs=[];ps=[]
    for ds,smp in keys:
        C,F,p=build_CFp(ds,smp); Cs.append(C);Fs.append(F);ps.append(p)
    m=len(Cs); n=len(FGW_NODES)
    out=ot.gromov.fgw_barycenters(n,Fs,Cs,ps,lambdas=[1.0/m]*m,alpha=args.alpha,loss_fun="square_loss",
        symmetric=False,max_iter=1000,p=np.ones(n)/n,init_C=np.mean(np.stack(Cs),0),init_X=np.mean(np.stack(Fs),0),
        random_state=SEED,log=True)
    return np.asarray(out[1]),np.asarray(out[0]),np.ones(n)/n
def fgw2(C,F,p,Cb,Fb,pb):
    M=ot.dist(F,Fb); M=M/(M.max()+1e-9)
    return float(ot.gromov.fused_gromov_wasserstein2(M,C,Cb,p,pb,loss_fun="square_loss",alpha=args.alpha,symmetric=False))

## -- assemble per-sample covariates + labels ----
df = idx.copy()
df["is_aml"]  = df["timepoint"].astype(str).isin(AML_TP).astype(int)
df["is_heal"] = (df["timepoint"].astype(str)=="Healthy").astype(int)
df = df[(df.is_aml==1) | (df.is_heal==1)].reset_index(drop=True)      # drop Unknown etc.

# blast_proxy from RAW n_cells (composition, inferCNV-free)
nc = nodes.pivot_table(index=["dataset","sample"], columns="hierarchy_bin", values="n_cells_raw", aggfunc="first").reset_index()
for b in FGW_NODES:
    if b not in nc: nc[b]=0
nc["total"] = nc[FGW_NODES].sum(axis=1)
nc["blast_proxy"] = nc[BLAST_BINS].sum(axis=1) / nc["total"].clip(lower=1)
df = df.merge(nc[["dataset","sample","blast_proxy"]], on=["dataset","sample"], how="left")

# malignant fraction (sensitivity covariate) from ccc_node_features
mf = nfraw.groupby(["dataset","sample"]).apply(
    lambda g: pd.Series({"mal_frac": g["n_malignant"].sum()/max(g["n_evaluable"].sum(),1)})).reset_index()
df = df.merge(mf, on=["dataset","sample"], how="left")

## -- HDS: LOO-clean for healthy, existing for AML ----
sc_key = sc.set_index(["dataset","sample"])["HDS"]
heal_all = [(r.dataset, r["sample"]) for _,r in df[df.is_heal==1].iterrows()]
heal_nonsparse = [k for k in heal_all if not bool(idx.set_index(["dataset","sample"]).loc[k,"sparse_flag"])]
print(f"[1] healthy total={len(heal_all)} nonsparse-in-barycenter={len(heal_nonsparse)} | LOO-recomputing healthy HDS")

hds_clean = {}
for k in heal_all:
    loo = [x for x in heal_nonsparse if x != k]
    Cb,Fb,pb = barycenter(loo)
    C,F,p = build_CFp(*k)
    hds_clean[k] = fgw2(C,F,p, Cb,Fb,pb)
for _,r in df[df.is_aml==1].iterrows():
    k=(r.dataset,r["sample"]); hds_clean[k]=float(sc_key.loc[k])     # AML already clean

df["HDS_clean"] = [hds_clean[(r.dataset,r["sample"])] for _,r in df.iterrows()]
df.to_csv(os.path.join(D_OUT,"h2_blast_regression.csv"), index=False)

## -- OLS via lstsq + permutation p for the is_aml partial effect ----
def ols(y, X):
    beta,_,_,_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X@beta; ss_res=float(((y-yhat)**2).sum()); ss_tot=float(((y-y.mean())**2).sum())
    return beta, (1-ss_res/ss_tot if ss_tot>0 else np.nan)

def run_model(cov_name):
    d = df.dropna(subset=["HDS_clean", cov_name]).copy()
    y = d["HDS_clean"].to_numpy(float)
    cov = d[cov_name].to_numpy(float); aml = d["is_aml"].to_numpy(float)
    n = len(y); ones = np.ones(n)
    X_full = np.column_stack([ones, cov, aml])
    X_red  = np.column_stack([ones, cov])
    beta_full, r2_full = ols(y, X_full)
    _,        r2_red  = ols(y, X_red)
    obs_aml_coef = beta_full[2]
    # permutation null: shuffle is_aml, refit, collect the is_aml coefficient
    perm = np.empty(args.n_perm)
    for i in range(args.n_perm):
        a = rng.permutation(aml)
        b,_ = ols(y, np.column_stack([ones, cov, a])); perm[i]=b[2]
    p_greater = (1 + int((perm >= obs_aml_coef).sum())) / (args.n_perm + 1)
    print(f"\n[2] model: HDS_clean ~ {cov_name} + is_aml   (n={n})")
    print(f"    beta[{cov_name}]={beta_full[1]:+.4f}  beta[is_aml]={obs_aml_coef:+.4f}")
    print(f"    R2(cov+is_aml)={r2_full:.3f}  vs  R2(cov only)={r2_red:.3f}   (increment from topology-label={r2_full-r2_red:+.3f})")
    print(f"    is_aml partial effect > 0 (topology beats composition), permutation p={p_greater:.4f}  [n_perm={args.n_perm}]")
    print(f"    group mean HDS_clean: healthy={y[aml==0].mean():.4f}  AML={y[aml==1].mean():.4f}")

print(f"[1] samples: healthy={int((df.is_heal==1).sum())}  AML={int((df.is_aml==1).sum())}")
run_model("blast_proxy")     # primary (inferCNV-free composition)
run_model("mal_frac")        # sensitivity (inferCNV-based malignant fraction)
print(f"\n[done] wrote {os.path.join(D_OUT,'h2_blast_regression.csv')}")
