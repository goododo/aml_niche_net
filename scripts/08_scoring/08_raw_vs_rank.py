#!/usr/bin/env python
# 08_raw_vs_rank.py ----
# Phase 7 (stage 08_scoring). DIAGNOSTIC 2, following the alpha sweep.
#
# THE QUESTION: alpha=1 (pure topology on the rank-distance C) showed no AML/healthy separation, and the
# per-edge regression on C found 0/49. But C is a WITHIN-SAMPLE rank-percentile, chosen (D2) to remove
# platform depth offsets. Rank is a lossy, zero-sum transform: it keeps only the ORDERING of a sample's 49
# edges and discards all absolute magnitude. If the disease signal lives in absolute communication
# strength, our own preprocessing destroyed it -- and "topology carries no signal" would be a statement
# about our distance definition, not about the biology.
#
# THE TEST: rerun the same inference on RAW edge strength, log1p(weight_probsum) (the CellChat
# aggregateNet weight computed in 06 before ranking), and run the IDENTICAL analysis on the rank C values
# so the two are directly comparable. Platform is controlled with dataset dummies + within-dataset
# permutation, exactly as in 04/05 (raw weights are platform-sensitive, so this control is essential).
#
# TWO LEVELS, for each of {raw log-weight, rank C}:
#   (i)  PER-EDGE : y_edge ~ 1 + blast_proxy + dataset + is_aml, FWL permutation, BH-FDR over 49 edges.
#   (ii) OMNIBUS  : is there ANY multivariate association across the 49-edge profile? Statistic =
#                   sum of squared FWL betas over column-standardized edges, permutation null. This is the
#                   aggregate analogue of HDS but computed directly on edges, with no FGW and no node
#                   features -- so it cannot be contaminated by the frac_malignant encoding.
#
# READING THE RESULT:
#   raw omnibus significant, rank omnibus null -> the rank transform destroyed real signal; the distance
#       definition must be revised (e.g. keep magnitude, control platform statistically instead).
#   both null -> edge-level communication genuinely does not separate AML from healthy in this cohort;
#       the negative topology conclusion stands and is not an artifact of ranking.
#   raw significant only in the GLOBAL model but not within-dataset -> raw magnitude is platform-driven,
#       which is exactly why ranking was adopted; no rescue.
#
# INPUT  : <root>/06_distance/edge_distance.csv (weight_probsum + C) , <root>/07_fgw/fgw_{nodes,input_index}_long.csv
# OUTPUT : <root>/08_scoring/raw_vs_rank.csv (per-edge results) + printed omnibus table
# Usage  : python scripts/08_scoring/08_raw_vs_rank.py [--n_perm 10000]
import argparse, os
import numpy as np
import pandas as pd

FGW_NODES=["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
BLAST_BINS=["HSC_MPP","LMPP_GMP","Mono_DC"]
AML_TP={"Diagnosis","MRD","Post_treatment","Relapse","Relapse2"}
DEFAULT_ROOT="/FAST/gr10634/gaozy/aml_niche_net/results/tables"; SEED=491638

ap=argparse.ArgumentParser()
ap.add_argument("--root",default=DEFAULT_ROOT)
ap.add_argument("--n_perm",type=int,default=10000)
ap.add_argument("--include_sparse",action="store_true")
args=ap.parse_args(); rng=np.random.default_rng(SEED)
D_DIST=os.path.join(args.root,"06_distance"); D_FGW=os.path.join(args.root,"07_fgw")
D_OUT=os.path.join(args.root,"08_scoring"); os.makedirs(D_OUT,exist_ok=True)

ed=pd.read_csv(os.path.join(D_DIST,"edge_distance.csv"))
nodes=pd.read_csv(os.path.join(D_FGW,"fgw_nodes_long.csv"))
idx=pd.read_csv(os.path.join(D_FGW,"fgw_input_index.csv"))

## -- cohort + covariates ----
s=idx.copy()
s["grp"]=np.where(s.timepoint=="Healthy","Healthy",np.where(s.timepoint.isin(AML_TP),"AML","other"))
s=s[s.grp!="other"].copy()
if not args.include_sparse: s=s[~s["sparse_flag"].fillna(False)]
s["is_aml"]=(s.grp=="AML").astype(float)
nc=nodes.pivot_table(index=["dataset","sample"],columns="hierarchy_bin",values="n_cells_raw",aggfunc="first").reset_index()
for b in FGW_NODES:
    if b not in nc: nc[b]=0
nc["total"]=nc[FGW_NODES].sum(axis=1); nc["blast_proxy"]=nc[BLAST_BINS].sum(axis=1)/nc["total"].clip(lower=1)
s=s.merge(nc[["dataset","sample","blast_proxy"]],on=["dataset","sample"],how="left").dropna(subset=["blast_proxy"]).reset_index(drop=True)

ed["edge"]=ed.sender_bin+"->"+ed.receiver_bin
W_raw=ed.pivot_table(index=["dataset","sample"],columns="edge",values="weight_probsum",aggfunc="first").reset_index()
W_rank=ed.pivot_table(index=["dataset","sample"],columns="edge",values="C",aggfunc="first").reset_index()
edge_cols=[c for c in W_raw.columns if c not in ("dataset","sample")]

m_raw=s.merge(W_raw,on=["dataset","sample"],how="left")
m_rank=s.merge(W_rank,on=["dataset","sample"],how="left")
Y_raw=np.log1p(np.nan_to_num(m_raw[edge_cols].to_numpy(float),nan=0.0))   # log1p: weights are >=0, heavy-tailed
Y_rank=np.nan_to_num(m_rank[edge_cols].to_numpy(float),nan=1.0)
print(f"[1] samples: healthy={int((s.is_aml==0).sum())} AML={int((s.is_aml==1).sum())} | edges={len(edge_cols)}")
print(f"[1] raw weight range: {np.expm1(Y_raw).min():.4f} .. {np.expm1(Y_raw).max():.4f}")

def zcol(Y):
    sd=Y.std(axis=0); sd[sd==0]=1.0
    return (Y-Y.mean(axis=0))/sd

def fwl_setup(X0):
    Q,_=np.linalg.qr(X0)
    return Q

def analyze(Y, X0, a, groups, n_perm, rng, label):
    """returns per-edge (beta,p) and the omnibus (stat,p)."""
    Q=fwl_setup(X0)
    def resid_mat(M): return M-Q@(Q.T@M)
    def resid_vec(v): return v-Q@(Q.T@v)
    rY=resid_mat(zcol(Y)); ra=resid_vec(a); den=float(ra@ra)
    beta=(rY.T@ra)/den
    omni=float(np.sum(beta**2))
    exceed_edge=np.zeros(Y.shape[1]); exceed_omni=0
    gidx=None if groups is None else [np.where(groups==g)[0] for g in np.unique(groups)]
    for _ in range(n_perm):
        if gidx is None: ap_=rng.permutation(a)
        else:
            ap_=a.copy()
            for gi in gidx: ap_[gi]=rng.permutation(ap_[gi])
        rp=resid_vec(ap_); dd=float(rp@rp)
        if dd<=0: continue
        b=(rY.T@rp)/dd
        exceed_edge+=(np.abs(b)>=np.abs(beta))
        exceed_omni+=(float(np.sum(b**2))>=omni)
    p_edge=(1.0+exceed_edge)/(n_perm+1.0)
    p_omni=(1.0+exceed_omni)/(n_perm+1.0)
    return beta,p_edge,omni,p_omni

def bh(p):
    p=np.asarray(p); M=len(p); o=np.argsort(p); q=np.empty(M)
    q[o]=np.minimum.accumulate((p[o]*M/np.arange(1,M+1))[::-1])[::-1]
    return np.clip(q,0,1)

## -- model matrices ----
n=len(s); a=s["is_aml"].to_numpy(float)
X0_glob=np.column_stack([np.ones(n),s["blast_proxy"].to_numpy(float)])
both=[ds for ds,g in s.groupby("dataset") if g.is_aml.nunique()==2]
mask=s.dataset.isin(both).to_numpy()
sb=s[mask].reset_index(drop=True); ab=sb["is_aml"].to_numpy(float)
dums=pd.get_dummies(sb["dataset"],drop_first=True).to_numpy(float)
X0_strat=np.column_stack([np.ones(len(sb)),sb["blast_proxy"].to_numpy(float),dums])
print(f"[2] within-dataset model uses {both} -> n={len(sb)}")

results={}; omni_rows=[]
for label,Y in (("raw_logweight",Y_raw),("rank_C",Y_rank)):
    bG,pG,oG,poG=analyze(Y,X0_glob,a,None,args.n_perm,rng,label)
    bS,pS,oS,poS=analyze(Y[mask],X0_strat,ab,sb["dataset"].to_numpy(),args.n_perm,rng,label)
    results[label]=dict(bG=bG,pG=pG,qG=bh(pG),bS=bS,pS=pS,qS=bh(pS))
    omni_rows.append(dict(matrix=label,model="global",stat=oG,p=poG,n=n))
    omni_rows.append(dict(matrix=label,model="within_dataset",stat=oS,p=poS,n=len(sb)))
    print(f"[3] {label:14s} | per-edge q<0.05: global={int((bh(pG)<0.05).sum())}/49 "
          f"within-ds={int((bh(pS)<0.05).sum())}/49 | OMNIBUS p: global={poG:.5f} within-ds={poS:.5f}",flush=True)

rows=[]
for j,e in enumerate(edge_cols):
    r=dict(edge=e)
    for label in results:
        R=results[label]
        r[f"{label}_beta_global"]=R["bG"][j]; r[f"{label}_q_global"]=R["qG"][j]
        r[f"{label}_beta_strat"]=R["bS"][j];  r[f"{label}_q_strat"]=R["qS"][j]
    rows.append(r)
res=pd.DataFrame(rows).sort_values("raw_logweight_q_strat").reset_index(drop=True)
res.to_csv(os.path.join(D_OUT,"raw_vs_rank.csv"),index=False)

om=pd.DataFrame(omni_rows)
print("\n[4] OMNIBUS multivariate test (49-edge profile; no FGW, no node features):")
print(om.round(5).to_string(index=False))
print("\n[5] top 8 edges by raw-weight within-dataset q:")
cols=["edge","raw_logweight_beta_strat","raw_logweight_q_strat","rank_C_beta_strat","rank_C_q_strat"]
print(res[cols].head(8).round(5).to_string(index=False))

pr=float(om[(om.matrix=="raw_logweight")&(om.model=="within_dataset")].p.iloc[0])
pk=float(om[(om.matrix=="rank_C")&(om.model=="within_dataset")].p.iloc[0])
prg=float(om[(om.matrix=="raw_logweight")&(om.model=="global")].p.iloc[0])
print("\n[6] VERDICT")
print(f"    omnibus within-dataset: raw log-weight p={pr:.5f} | rank C p={pk:.5f}")
if pr<0.05 and pk>=0.05:
    print("    -> RANK DESTROYED SIGNAL: absolute edge strength separates AML/healthy even within dataset,")
    print("       while the rank-transformed version does not. The distance definition (D2) must be revised.")
elif pr>=0.05 and pk>=0.05:
    if prg<0.05:
        print("    -> raw weights separate globally but NOT within dataset: absolute magnitude is platform-driven,")
        print("       which is precisely why ranking was adopted. No rescue; the negative result stands.")
    else:
        print("    -> NEITHER raw nor ranked edge profiles separate AML from healthy. Edge-level communication")
        print("       genuinely carries no detectable disease signal here; the negative topology conclusion is")
        print("       NOT an artifact of the rank transform.")
else:
    print("    -> both carry signal; compare effect sizes and revisit the distance definition.")
print(f"\n[done] wrote {os.path.join(D_OUT,'raw_vs_rank.csv')}")
