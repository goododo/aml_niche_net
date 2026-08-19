#!/usr/bin/env python
# 07_feature_decomposition.py ----
# Phase 7 (stage 08_scoring). DIAGNOSTIC 1, following the alpha sweep.
#
# CONTEXT: 06_alpha_sweep showed a perfectly monotone dose-response -- separation scales with the FEATURE
# weight and vanishes at alpha=1 (pure topology: healthy 0.048 vs AML 0.052, p=0.271). So whatever signal
# HDS carries lives in the node-FEATURE term, not the edge structure. That raises a circularity worry:
# one of the three features is frac_malignant, which we FORCE to 0 for healthy samples by design
# (FGW_ZERO_HEALTHY_MAL). If that single feature drives everything, the significance is an artifact of our
# own encoding and cannot be reported as a finding at all.
#
# THIS SCRIPT re-runs the whole HDS pipeline (LOO healthy barycenter + scoring + the two regression models)
# separately for each FEATURE SUBSET:
#   all3          : frac_malignant, mean_stemness, n_cells      (the current definition)
#   no_frac_mal   : mean_stemness, n_cells                      <- KEY: signal here is NOT circular
#   only_frac_mal : frac_malignant                              <- KEY: isolates the encoded label
#   only_stemness : mean_stemness                               (transcriptional, CNV-independent)
#   only_ncells   : n_cells                                     (pure composition; already covered by blast_proxy)
#
# READING THE RESULT:
#   only_frac_mal strong AND no_frac_mal null  -> CIRCULAR. The result is an artifact of our encoding.
#   no_frac_mal still significant              -> a real (if non-topological) feature signal survives;
#                                                 then check WHICH: stemness (biology) vs n_cells (composition).
#   only_ncells carries it                     -> trivially compositional, and blast_proxy already covers it.
# Run at alpha=0 (pure feature, cleanest attribution) and alpha=0.5 (the default setting).
#
# INPUT  : <root>/07_fgw/fgw_{nodes,edges,input_index}_long.csv
# OUTPUT : <root>/08_scoring/feature_decomposition.csv
# Usage  : python scripts/08_scoring/07_feature_decomposition.py [--n_perm 10000]
import argparse, os, warnings
import sys
import numpy as np
import pandas as pd
import ot

# The timepoint vocabulary is LOADED, not literal. A hard-coded set here silently deleted 34 of
# 214 samples (64% of the treated arm) after CANONICAL_TIMEPOINTS changed on 2026-08-04.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config'))
from fgw_vocab import load_vocab, assert_index_covered

warnings.filterwarnings("ignore", category=RuntimeWarning)   # alpha=0 divide-by-zero in POT's log term only

FGW_NODES=["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
ALL_FEATURES=["frac_malignant","mean_stemness","n_cells"]
BLAST_BINS=["HSC_MPP","LMPP_GMP","Mono_DC"]
AML_TP = None  # set from fgw_vocab.json below -- see scripts/config/fgw_vocab.py
DEFAULT_ROOT="/FAST/gr10634/gaozy/aml_niche_net/results/tables"; SEED=491638

SUBSETS={
 "all3":          ["frac_malignant","mean_stemness","n_cells"],
 "no_frac_mal":   ["mean_stemness","n_cells"],
 "only_frac_mal": ["frac_malignant"],
 "only_stemness": ["mean_stemness"],
 "only_ncells":   ["n_cells"],
}

ap=argparse.ArgumentParser()
ap.add_argument("--root",default=DEFAULT_ROOT)
ap.add_argument("--n_perm",type=int,default=10000)
ap.add_argument("--alphas",default="0,0.5")
ap.add_argument("--max_iter",type=int,default=1000)
args=ap.parse_args(); rng=np.random.default_rng(SEED)
ALPHAS=[float(x) for x in args.alphas.split(",")]
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

_cache={}
def build_one(ds,smp,feats):
    k=(ds,smp,tuple(feats))
    if k in _cache: return _cache[k]
    e=edges[(edges["dataset"]==ds)&(edges["sample"]==smp)]
    C=e.pivot(index="sender_bin",columns="receiver_bin",values="C").reindex(FGW_NODES,columns=FGW_NODES).to_numpy(float)
    C=np.nan_to_num(C,nan=1.0)
    nd=nodes[(nodes["dataset"]==ds)&(nodes["sample"]==smp)].set_index("hierarchy_bin").reindex(FGW_NODES)
    F=np.nan_to_num(nd[feats].to_numpy(float),nan=0.0)
    if F.ndim==1: F=F.reshape(-1,1)
    p=np.nan_to_num(nd["mass"].to_numpy(float),nan=1e-6); p=p/p.sum()
    _cache[k]=(C,F,p); return _cache[k]

def barycenter(keys,alpha,feats):
    Cs=[];Fs=[];ps=[]
    for k in keys:
        C,F,p=build_one(*k,feats); Cs.append(C);Fs.append(F);ps.append(p)
    m=len(Cs); n=7
    out=ot.gromov.fgw_barycenters(n,Fs,Cs,ps,lambdas=[1.0/m]*m,alpha=alpha,loss_fun="square_loss",
        symmetric=False,max_iter=args.max_iter,p=np.ones(n)/n,init_C=np.mean(np.stack(Cs),0),
        init_X=np.mean(np.stack(Fs),0),random_state=SEED,log=True)
    return np.asarray(out[1]),np.asarray(out[0]),np.ones(n)/n

def fgw2(C,F,p,Cb,Fb,pb,alpha):
    M=ot.dist(F,Fb); mx=M.max(); M=M/(mx+1e-9) if mx>0 else M
    return float(ot.gromov.fused_gromov_wasserstein2(M,C,Cb,p,pb,loss_fun="square_loss",alpha=alpha,symmetric=False))

def fwl_perm(y,X0,a,groups,n_perm,rng):
    Q,_=np.linalg.qr(X0)
    def resid(v): return v-Q@(Q.T@v)
    ry=resid(y); ra=resid(a); den=float(ra@ra)
    if den<=0: return np.nan,np.nan
    beta=float(ry@ra/den)
    gidx=None if groups is None else [np.where(groups==g)[0] for g in np.unique(groups)]
    exceed=0
    for _ in range(n_perm):
        if gidx is None: ap_=rng.permutation(a)
        else:
            ap_=a.copy()
            for gi in gidx: ap_[gi]=rng.permutation(ap_[gi])
        rp=resid(ap_); dd=float(rp@rp)
        if dd<=0: continue
        exceed+=(abs(float(ry@rp/dd))>=abs(beta))
    return beta,(1.0+exceed)/(n_perm+1.0)

## -- cohort ----
d=idx.copy()
d["grp"]=np.where(d.timepoint=="Healthy","Healthy",np.where(d.timepoint.isin(AML_TP),"AML","other"))
d=d[d.grp!="other"].copy(); d["is_aml"]=(d.grp=="AML").astype(float)
nc=nodes.pivot_table(index=["dataset","sample"],columns="hierarchy_bin",values="n_cells_raw",aggfunc="first").reset_index()
for b in FGW_NODES:
    if b not in nc: nc[b]=0
nc["total"]=nc[FGW_NODES].sum(axis=1); nc["blast_proxy"]=nc[BLAST_BINS].sum(axis=1)/nc["total"].clip(lower=1)
d=d.merge(nc[["dataset","sample","blast_proxy"]],on=["dataset","sample"],how="left").dropna(subset=["blast_proxy"]).reset_index(drop=True)
heal_all=[(r.dataset,r["sample"]) for _,r in d[d.is_aml==0].iterrows()]
heal_bary=[(r.dataset,r["sample"]) for _,r in d[(d.is_aml==0)&(~d["sparse_flag"].fillna(False))].iterrows()]
aml_all=[(r.dataset,r["sample"]) for _,r in d[d.is_aml==1].iterrows()]
print(f"[1] healthy={len(heal_all)} (barycenter {len(heal_bary)}) | AML={len(aml_all)} | alphas={ALPHAS}")

rows=[]
for alpha in ALPHAS:
    for name,feats in SUBSETS.items():
        Cb_full,Fb_full,pb=barycenter(heal_bary,alpha,feats)
        hds={}
        for k in heal_all:
            loo=[x for x in heal_bary if x!=k]
            Cb,Fb,_=barycenter(loo,alpha,feats)
            C,F,p=build_one(*k,feats); hds[k]=fgw2(C,F,p,Cb,Fb,pb,alpha)
        for k in aml_all:
            C,F,p=build_one(*k,feats); hds[k]=fgw2(C,F,p,Cb_full,Fb_full,pb,alpha)
        dd=d.copy(); dd["HDS"]=[hds[(r.dataset,r["sample"])] for _,r in dd.iterrows()]
        y=dd["HDS"].to_numpy(float); a=dd["is_aml"].to_numpy(float)
        X0=np.column_stack([np.ones(len(dd)),dd["blast_proxy"].to_numpy(float)])
        bA,pA=fwl_perm(y,X0,a,None,args.n_perm,rng)
        both=[ds for ds,g in dd.groupby("dataset") if g.is_aml.nunique()==2]
        db=dd[dd.dataset.isin(both)].reset_index(drop=True)
        yb=db["HDS"].to_numpy(float); ab=db["is_aml"].to_numpy(float)
        dums=pd.get_dummies(db["dataset"],drop_first=True).to_numpy(float)
        X0b=np.column_stack([np.ones(len(db)),db["blast_proxy"].to_numpy(float),dums])
        bB,pB=fwl_perm(yb,X0b,ab,db["dataset"].to_numpy(),args.n_perm,rng)
        rows.append(dict(alpha=alpha,feature_set=name,features="+".join(feats),
                         mean_healthy=float(y[a==0].mean()),mean_aml=float(y[a==1].mean()),
                         beta_global=bA,p_global=pA,beta_strat=bB,p_strat=pB))
        print(f"[2] alpha={alpha:<4g} {name:14s} | healthy={y[a==0].mean():.4f} AML={y[a==1].mean():.4f} "
              f"| global b={bA:+.5f} p={pA:.5f} | within-ds b={bB:+.5f} p={pB:.5f}",flush=True)

res=pd.DataFrame(rows)
res.to_csv(os.path.join(D_OUT,"feature_decomposition.csv"),index=False)
print("\n[3] FEATURE DECOMPOSITION")
print(res[["alpha","feature_set","mean_healthy","mean_aml","beta_global","p_global","beta_strat","p_strat"]].round(5).to_string(index=False))

print("\n[4] VERDICT (alpha=0, within-dataset model)")
z=res[(res.alpha==0.0)].set_index("feature_set")
def g(k,c):
    return float(z.loc[k,c]) if k in z.index else float("nan")
p_no=g("no_frac_mal","p_strat"); p_only=g("only_frac_mal","p_strat")
print(f"    only_frac_mal : p={p_only:.5f}   (the feature we force to 0 for healthy)")
print(f"    no_frac_mal   : p={p_no:.5f}   (stemness + n_cells, no encoded label)")
if p_only<0.05 and p_no>=0.05:
    print("    -> CIRCULAR: the separation comes from the feature we encoded by design.")
    print("       It cannot be reported as a finding.")
elif p_no<0.05:
    print("    -> A non-circular feature signal survives. Check only_stemness vs only_ncells below to see")
    print("       whether it is biological (stemness) or merely compositional (n_cells, already covered by blast_proxy).")
    print(f"       only_stemness p={g('only_stemness','p_strat'):.5f} | only_ncells p={g('only_ncells','p_strat'):.5f}")
else:
    print("    -> neither subset is significant on its own; the signal needs the full feature combination.")
print(f"\n[done] wrote {os.path.join(D_OUT,'feature_decomposition.csv')}")
