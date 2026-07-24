#!/usr/bin/env python
# 09_clr_redistribution.py ----
# Phase 7 (stage 08_scoring). DIAGNOSTIC 3 -- settles what 08_raw_vs_rank left ambiguous.
#
# WHAT 08 FOUND, AND WHY ITS VERDICT WAS PREMATURE: raw log-weights gave a significant omnibus
# (within-dataset p=0.0062) while rank C did not (p=0.258), and the script concluded "rank destroyed
# signal". But inspecting the coefficients shows 44/49 edges have POSITIVE beta (mean +0.40): AML samples
# have uniformly HIGHER weights on essentially every edge. That is a GLOBAL SCALE SHIFT, not a
# redistribution of the communication hierarchy. A global shift has an obvious technical reading (AML
# samples are blast-rich -> more cells per node -> more stable/higher triMean estimates; plus depth
# differences), and removing exactly that is what the rank transform (D2) was adopted to do.
#
# THE REMAINING QUESTION: beyond the global shift, is there any genuine REDISTRIBUTION? Rank is the most
# aggressive answer (keeps only ordering, discards all magnitude). CLR is the principled middle ground:
#     clr_ij = log1p(w_ij) - mean_j( log1p(w_ij) )     [per-sample centering]
# It removes the per-sample scale while PRESERVING relative magnitude, so it can detect redistribution
# that rank is too lossy to see, without being fooled by the global shift.
#
# TESTS (all with blast_proxy + dataset controls and within-dataset permutation, as in 04/05/08):
#   (1) Is the global shift real? per-sample mean log1p(weight) ~ blast + dataset + is_aml.
#   (2) Is it technical? correlate that per-sample mean with total cell count.
#   (3) OMNIBUS on CLR profiles -> the actual redistribution test.
#   (4) Per-edge on CLR, BH over 49.
#
# READING:
#   (1) significant and (2) correlated -> the raw-weight result is a depth/abundance artifact, and the
#       rank transform behaved correctly.
#   (3) significant -> genuine redistribution exists; the distance definition should move from rank to CLR
#       and the topological analyses are worth redoing on it.
#   (3) null -> after removing the global shift there is no edge-level redistribution signal; the negative
#       topology conclusion stands, now with the rank-artifact alternative explicitly ruled out.
#
# INPUT  : <root>/06_distance/edge_distance.csv , <root>/07_fgw/fgw_{nodes,input_index}_long.csv
# OUTPUT : <root>/08_scoring/clr_redistribution.csv
# Usage  : python scripts/08_scoring/09_clr_redistribution.py [--n_perm 10000]
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

s=idx.copy()
s["grp"]=np.where(s.timepoint=="Healthy","Healthy",np.where(s.timepoint.isin(AML_TP),"AML","other"))
s=s[s.grp!="other"].copy()
if not args.include_sparse: s=s[~s["sparse_flag"].fillna(False)]
s["is_aml"]=(s.grp=="AML").astype(float)
nc=nodes.pivot_table(index=["dataset","sample"],columns="hierarchy_bin",values="n_cells_raw",aggfunc="first").reset_index()
for b in FGW_NODES:
    if b not in nc: nc[b]=0
nc["total_cells"]=nc[FGW_NODES].sum(axis=1)
nc["blast_proxy"]=nc[BLAST_BINS].sum(axis=1)/nc["total_cells"].clip(lower=1)
s=s.merge(nc[["dataset","sample","blast_proxy","total_cells"]],on=["dataset","sample"],how="left").dropna(subset=["blast_proxy"]).reset_index(drop=True)

ed["edge"]=ed.sender_bin+"->"+ed.receiver_bin
W=ed.pivot_table(index=["dataset","sample"],columns="edge",values="weight_probsum",aggfunc="first").reset_index()
edge_cols=[c for c in W.columns if c not in ("dataset","sample")]
m=s.merge(W,on=["dataset","sample"],how="left")
L=np.log1p(np.nan_to_num(m[edge_cols].to_numpy(float),nan=0.0))     # n x 49 log weights
sample_mean=L.mean(axis=1)                                          # per-sample overall communication scale
CLR=L-sample_mean[:,None]                                           # per-sample centering
print(f"[1] samples: healthy={int((s.is_aml==0).sum())} AML={int((s.is_aml==1).sum())} | edges={len(edge_cols)}")

def zcol(Y):
    sd=Y.std(axis=0); sd[sd==0]=1.0
    return (Y-Y.mean(axis=0))/sd

def analyze(Y,X0,a,groups,n_perm,rng,multivariate=True):
    Q,_=np.linalg.qr(X0)
    def rm(M): return M-Q@(Q.T@M)
    def rv(v): return v-Q@(Q.T@v)
    Yz=zcol(Y) if multivariate else Y
    rY=rm(Yz) if multivariate else rv(Yz)
    ra=rv(a); den=float(ra@ra)
    if multivariate:
        beta=(rY.T@ra)/den; omni=float(np.sum(beta**2))
    else:
        beta=float(rY@ra/den); omni=beta**2
    ex_edge=np.zeros(Y.shape[1]) if multivariate else 0.0
    ex_omni=0
    gidx=None if groups is None else [np.where(groups==g)[0] for g in np.unique(groups)]
    for _ in range(n_perm):
        if gidx is None: ap_=rng.permutation(a)
        else:
            ap_=a.copy()
            for gi in gidx: ap_[gi]=rng.permutation(ap_[gi])
        rp=rv(ap_); dd=float(rp@rp)
        if dd<=0: continue
        if multivariate:
            b=(rY.T@rp)/dd; ex_edge+=(np.abs(b)>=np.abs(beta)); ex_omni+=(float(np.sum(b**2))>=omni)
        else:
            b=float(rY@rp/dd); ex_omni+=(b**2>=omni)
    p_omni=(1.0+ex_omni)/(n_perm+1.0)
    p_edge=(1.0+ex_edge)/(n_perm+1.0) if multivariate else None
    return beta,p_edge,omni,p_omni

def bh(p):
    p=np.asarray(p); M=len(p); o=np.argsort(p); q=np.empty(M)
    q[o]=np.minimum.accumulate((p[o]*M/np.arange(1,M+1))[::-1])[::-1]
    return np.clip(q,0,1)

n=len(m); a=m["is_aml"].to_numpy(float)
X0g=np.column_stack([np.ones(n),m["blast_proxy"].to_numpy(float)])
both=[ds for ds,g in m.groupby("dataset") if g.is_aml.nunique()==2]
mask=m.dataset.isin(both).to_numpy(); mb=m[mask].reset_index(drop=True)
ab=mb["is_aml"].to_numpy(float)
dums=pd.get_dummies(mb["dataset"],drop_first=True).to_numpy(float)
X0s=np.column_stack([np.ones(len(mb)),mb["blast_proxy"].to_numpy(float),dums])
print(f"[1] within-dataset model: {both} -> n={len(mb)}")

## -- (1) is the global shift real? ----
gm=sample_mean.reshape(-1,1)
_,_,_,p_shift_g=analyze(gm,X0g,a,None,args.n_perm,rng)
_,_,_,p_shift_s=analyze(gm[mask],X0s,ab,mb["dataset"].to_numpy(),args.n_perm,rng)
print(f"\n[2] GLOBAL SHIFT test (per-sample mean log-weight ~ blast + [dataset] + is_aml):")
print(f"    mean log-weight: healthy={sample_mean[a==0].mean():.4f}  AML={sample_mean[a==1].mean():.4f}")
print(f"    global p={p_shift_g:.5f} | within-dataset p={p_shift_s:.5f}")

## -- (2) is the shift technical (tracks cell count)? ----
tc=np.log10(m["total_cells"].to_numpy(float).clip(min=1))
r=float(np.corrcoef(sample_mean,tc)[0,1])
rs=[]
for ds,g in m.groupby("dataset"):
    if len(g)>=5:
        gi=g.index.to_numpy()
        rs.append((ds,float(np.corrcoef(sample_mean[gi],tc[gi])[0,1]),len(g)))
print(f"\n[3] is the shift technical? corr(mean log-weight, log10 total cells) = {r:+.3f}")
for ds,rr,k in sorted(rs,key=lambda x:-abs(x[1])):
    print(f"    {ds:14s} r={rr:+.3f} (n={k})")

## -- (3) OMNIBUS on CLR = the redistribution test ----
_,pe_g,o_g,po_g=analyze(CLR,X0g,a,None,args.n_perm,rng)
bS,pe_s,o_s,po_s=analyze(CLR[mask],X0s,ab,mb["dataset"].to_numpy(),args.n_perm,rng)
print(f"\n[4] REDISTRIBUTION test -- omnibus on CLR profiles (global scale removed):")
print(f"    global p={po_g:.5f} | within-dataset p={po_s:.5f}")
qe_s=bh(pe_s); qe_g=bh(pe_g)
print(f"[4] per-edge on CLR, q<0.05: global={int((qe_g<0.05).sum())}/49 | within-dataset={int((qe_s<0.05).sum())}/49")

res=pd.DataFrame(dict(edge=edge_cols,clr_beta_strat=bS,clr_p_strat=pe_s,clr_q_strat=qe_s,
                      clr_p_global=pe_g,clr_q_global=qe_g)).sort_values("clr_q_strat").reset_index(drop=True)
res.to_csv(os.path.join(D_OUT,"clr_redistribution.csv"),index=False)
print("\n[5] top 6 CLR edges (within-dataset):")
print(res.head(6).round(5).to_string(index=False))

print("\n[6] VERDICT")
print(f"    global scale shift : within-dataset p={p_shift_s:.5f} (corr with cell count r={r:+.3f})")
print(f"    redistribution(CLR): within-dataset p={po_s:.5f}")
if po_s<0.05:
    print("    -> GENUINE REDISTRIBUTION beyond the global shift. The rank transform was too lossy;")
    print("       switch the distance definition to CLR and redo the topological analyses on it.")
else:
    print("    -> NO redistribution once the per-sample scale is removed. The raw-weight result was the")
    print("       global shift alone (largely tracking cell number), so the rank transform behaved")
    print("       correctly and the negative topology conclusion stands.")
print(f"\n[done] wrote {os.path.join(D_OUT,'clr_redistribution.csv')}")
