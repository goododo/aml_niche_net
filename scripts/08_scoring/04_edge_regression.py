#!/usr/bin/env python
# 04_edge_regression.py ----
# Phase 7 (stage 08_scoring). Per-edge inference using SAMPLES AS REPLICATES -- the properly powered
# way to ask "which edges differ between AML and healthy".
#
# WHY THIS, AFTER 02/03 WERE NULL: 02 (per-edge) and 03 (per-block) both collapsed 113 samples into TWO
# barycenters and built the null from random 20/93 splits. That discards all within-group replication and
# the small healthy group (n=20) makes the null very wide -- effects came out <1 SD of null. Here each
# sample is a replicate (n=113), which is the same design that DID work for HDS in 08/01 (p=0.0005).
# This is NOT post-hoc selection: all 49 edges are tested, BH-corrected, with the pre-specified labels.
#
# MODELS (both reported; they trade power against platform strictness):
#   A "global"     : C_edge ~ 1 + blast_proxy + is_aml           over all healthy+AML samples.
#                    Permutation: shuffle is_aml globally. Higher power; does NOT control dataset.
#   B "stratified" : C_edge ~ 1 + blast_proxy + dataset_dummies + is_aml, restricted to datasets that
#                    contain BOTH healthy and AML; permutation shuffles is_aml WITHIN dataset. The
#                    is_aml effect is then identified only from WITHIN-dataset contrasts -> platform is
#                    controlled by design (the strict confirmation of A).
#
# SIGN: C is rank-distance, LOWER = stronger. beta_is_aml < 0 => that edge is STRONGER in AML.
# Inference by permutation via Frisch-Waugh-Lovell residualization (exact, and fast enough that 10k
# permutations run in seconds -- no SLURM needed).
# CAVEAT retained from 03: C is a within-sample rank -> the 49 edges are zero-sum, so edges are
# compositionally coupled; results describe RE-PRIORITIZATION of the hierarchy, not absolute gain/loss.
#
# INPUT  : <root>/07_fgw/fgw_{nodes,edges,input_index}_long.csv
# OUTPUT : <root>/08_scoring/edge_regression.csv
# Usage  : python scripts/08_scoring/04_edge_regression.py [--n_perm 10000] [--include_sparse]
import argparse, os
import sys
import numpy as np
import pandas as pd

# The timepoint vocabulary is LOADED, not literal. A hard-coded set here silently deleted 34 of
# 214 samples (64% of the treated arm) after CANONICAL_TIMEPOINTS changed on 2026-08-04.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config'))
from fgw_vocab import load_vocab, assert_index_covered


FGW_NODES=["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
BLAST_BINS=["HSC_MPP","LMPP_GMP","Mono_DC"]
AML_TP = None  # set from fgw_vocab.json below -- see scripts/config/fgw_vocab.py
DEFAULT_ROOT="/FAST/gr10634/gaozy/aml_niche_net/results/tables"; SEED=491638

ap=argparse.ArgumentParser()
ap.add_argument("--root",default=DEFAULT_ROOT)
ap.add_argument("--n_perm",type=int,default=10000)
ap.add_argument("--include_sparse",action="store_true",
                help="include sparse-flagged graphs (default: exclude, matching the barycenter cohort)")
args=ap.parse_args(); rng=np.random.default_rng(SEED)
D_FGW=os.path.join(args.root,"07_fgw"); D_OUT=os.path.join(args.root,"08_scoring"); os.makedirs(D_OUT,exist_ok=True)
_VOCAB = load_vocab(D_FGW)
AML_TP = _VOCAB["aml_timepoints"]


edges=pd.read_csv(os.path.join(D_FGW,"fgw_edges_long.csv"))
nodes=pd.read_csv(os.path.join(D_FGW,"fgw_nodes_long.csv"))
idx=pd.read_csv(os.path.join(D_FGW,"fgw_input_index.csv"))

# Fail loudly on a label the vocabulary does not cover. Without this an unrecognised
# timepoint is not an error -- it is a sample that quietly stops being AML and stops
# being healthy, and therefore stops existing for every test below.
assert_index_covered(idx, _VOCAB)

## -- sample table: label + blast_proxy ----
s=idx.copy()
s["grp"]=np.where(s.timepoint=="Healthy","Healthy",np.where(s.timepoint.isin(AML_TP),"AML","other"))
s=s[s.grp!="other"].copy()
if not args.include_sparse: s=s[~s["sparse_flag"].fillna(False)]
s["is_aml"]=(s.grp=="AML").astype(float)

nc=nodes.pivot_table(index=["dataset","sample"],columns="hierarchy_bin",values="n_cells_raw",aggfunc="first").reset_index()
for b in FGW_NODES:
    if b not in nc: nc[b]=0
nc["total"]=nc[FGW_NODES].sum(axis=1)
nc["blast_proxy"]=nc[BLAST_BINS].sum(axis=1)/nc["total"].clip(lower=1)
s=s.merge(nc[["dataset","sample","blast_proxy"]],on=["dataset","sample"],how="left")
s=s.dropna(subset=["blast_proxy"]).reset_index(drop=True)
print(f"[1] samples: healthy={int((s.is_aml==0).sum())} AML={int((s.is_aml==1).sum())} "
      f"(sparse {'included' if args.include_sparse else 'excluded'})")

## -- edge matrix: samples x 49 edges ----
edges["edge"]=edges.sender_bin+"->"+edges.receiver_bin
W=edges.pivot_table(index=["dataset","sample"],columns="edge",values="C",aggfunc="first").reset_index()
mrg=s.merge(W,on=["dataset","sample"],how="left")
edge_cols=[c for c in W.columns if c not in ("dataset","sample")]
Y=mrg[edge_cols].to_numpy(float)              # n_samples x 49
Y=np.nan_to_num(Y,nan=1.0)
print(f"[1] edge matrix: {Y.shape[0]} samples x {Y.shape[1]} edges")

def fwl_perm(Y, X0, a, groups, n_perm, rng):
    """Frisch-Waugh-Lovell: beta_a from Y ~ [X0, a] equals beta from resid(Y|X0) ~ resid(a|X0).
    groups=None -> global shuffle of a; else shuffle a WITHIN each group (stratified)."""
    Q,_=np.linalg.qr(X0)                       # orthonormal basis of the nuisance space
    def resid(v): return v - Q@(Q.T@v)
    rY=Y - Q@(Q.T@Y)                           # n x m residualized responses
    ra=resid(a)
    denom=float(ra@ra)
    if denom<=0: raise ValueError("is_aml has no residual variance after covariates")
    beta_obs=(rY.T@ra)/denom                   # m
    exceed=np.zeros(Y.shape[1])
    if groups is None:
        for _ in range(n_perm):
            ap_=rng.permutation(a)
            rp=resid(ap_); d=float(rp@rp)
            if d<=0: continue
            b=(rY.T@rp)/d
            exceed+=(np.abs(b)>=np.abs(beta_obs))
    else:
        gidx=[np.where(groups==g)[0] for g in np.unique(groups)]
        for _ in range(n_perm):
            ap_=a.copy()
            for gi in gidx: ap_[gi]=rng.permutation(ap_[gi])
            rp=resid(ap_); d=float(rp@rp)
            if d<=0: continue
            b=(rY.T@rp)/d
            exceed+=(np.abs(b)>=np.abs(beta_obs))
    p=(1.0+exceed)/(n_perm+1.0)
    return beta_obs,p

def bh(p):
    p=np.asarray(p); m=len(p); o=np.argsort(p); q=np.empty(m)
    q[o]=np.minimum.accumulate((p[o]*m/np.arange(1,m+1))[::-1])[::-1]
    return np.clip(q,0,1)

## -- Model A: global ----
n=len(mrg)
X0=np.column_stack([np.ones(n), mrg["blast_proxy"].to_numpy(float)])
a=mrg["is_aml"].to_numpy(float)
betaA,pA=fwl_perm(Y,X0,a,None,args.n_perm,rng)
qA=bh(pA)
print(f"[2] model A (global, n={n}, blast-controlled): edges q<0.05 = {int((qA<0.05).sum())}/49")

## -- Model B: within-dataset (only datasets containing BOTH labels) ----
both=[d_ for d_,g in mrg.groupby("dataset") if g.is_aml.nunique()==2]
mb=mrg[mrg.dataset.isin(both)].reset_index(drop=True)
print(f"[3] model B datasets with both labels: {both} -> n={len(mb)} "
      f"(healthy={int((mb.is_aml==0).sum())}, AML={int((mb.is_aml==1).sum())})")
Yb=np.nan_to_num(mb[edge_cols].to_numpy(float),nan=1.0)
dums=pd.get_dummies(mb["dataset"],drop_first=True).to_numpy(float)
X0b=np.column_stack([np.ones(len(mb)), mb["blast_proxy"].to_numpy(float), dums])
ab=mb["is_aml"].to_numpy(float)
betaB,pB=fwl_perm(Yb,X0b,ab,mb["dataset"].to_numpy(),args.n_perm,rng)
qB=bh(pB)
print(f"[3] model B (within-dataset, platform-controlled): edges q<0.05 = {int((qB<0.05).sum())}/49")

## -- assemble ----
res=pd.DataFrame(dict(edge=edge_cols,
                      beta_global=betaA,p_global=pA,q_global=qA,
                      beta_strat=betaB,p_strat=pB,q_strat=qB))
res[["sender_bin","receiver_bin"]]=res["edge"].str.split("->",expand=True)
res["direction_global"]=np.where(res.beta_global<0,"stronger_in_AML","weaker_in_AML")
res["concordant_sign"]=np.sign(res.beta_global)==np.sign(res.beta_strat)
res["robust"]=(res.q_global<0.05)&(res.q_strat<0.05)&res.concordant_sign
res=res.sort_values("q_global").reset_index(drop=True)
res.to_csv(os.path.join(D_OUT,"edge_regression.csv"),index=False)

cols=["sender_bin","receiver_bin","beta_global","q_global","beta_strat","q_strat","concordant_sign"]
print("\n[4] top 12 edges by global q:")
print(res[cols].head(12).round(5).to_string(index=False))
print(f"\n[4] ROBUST edges (q<0.05 in BOTH models, same sign): {int(res.robust.sum())}")
rb=res[res.robust].sort_values("beta_global")
print(rb[cols].round(5).to_string(index=False) if len(rb) else "    (none)")
print(f"\n[done] wrote {os.path.join(D_OUT,'edge_regression.csv')}")
