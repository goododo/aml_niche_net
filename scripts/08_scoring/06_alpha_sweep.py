#!/usr/bin/env python
# 06_alpha_sweep.py ----
# Phase 7 (stage 08_scoring). THE DECISIVE TEST for whether H2 is about TOPOLOGY or about NODE FEATURES.
#
# THE PROBLEM THIS RESOLVES: HDS is an FGW distance, and at the default alpha=0.5 half the weight sits on
# the node-FEATURE term. One of those features is frac_malignant, which we deliberately force to 0 for
# healthy samples (FGW_ZERO_HEALTHY_MAL). So part of the healthy-vs-AML separation could be near-circular
# ("healthy samples differ because we encoded them as having no malignant cells"), which is exactly the
# R5 attack surface. Until this is settled we cannot claim the signal is topological.
#
# POT convention: FGW loss = (1-alpha)*<T,M_features> + alpha*GW(C_structure).
#   alpha = 1  -> PURE GW: node features are IGNORED ENTIRELY. frac_malignant cannot contribute at all,
#                 so the circularity concern is structurally impossible. If H2 holds here, "topological
#                 deviation" is established.
#   alpha = 0  -> PURE feature OT: structure ignored. If the signal lives only here, the claim must be
#                 restated as compositional/feature-based, not topological.
# Sweep = the blueprint's pre-specified {0, 0.25, 0.5, 0.75, 1}; nothing is chosen after seeing results.
#
# At every alpha the FULL 08/01 + 08/05 pipeline is re-run: leave-one-out healthy barycenter (removing the
# circularity of scoring healthy samples against a barycenter they helped build), HDS for all samples, then
#   model A "global"     : HDS ~ 1 + blast_proxy + is_aml           (global permutation of is_aml)
#   model B "stratified" : HDS ~ 1 + blast_proxy + dataset + is_aml (permutation WITHIN dataset;
#                          platform controlled by design, restricted to datasets holding both labels)
# Inference by FWL residualization + permutation, identical to 04/05 so numbers are directly comparable.
#
# INPUT  : <root>/07_fgw/fgw_{nodes,edges,input_index}_long.csv
# OUTPUT : <root>/08_scoring/alpha_sweep.csv  (per alpha: beta/p for both models + group means)
# Usage  : python scripts/08_scoring/06_alpha_sweep.py [--n_perm 10000] [--alphas 0,0.25,0.5,0.75,1]
import argparse, hashlib, os, warnings
from contextlib import contextmanager
import sys
import numpy as np
import pandas as pd
import ot

# The timepoint vocabulary is LOADED, not literal. A hard-coded set here silently deleted 34 of
# 214 samples (64% of the treated arm) after CANONICAL_TIMEPOINTS changed on 2026-08-04.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config'))
from fgw_vocab import load_vocab, load_features, assert_index_covered


FGW_NODES=["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
BLAST_BINS=["HSC_MPP","LMPP_GMP","Mono_DC"]
AML_TP = None  # set from fgw_vocab.json below -- see scripts/config/fgw_vocab.py
DEFAULT_ROOT="/FAST/gr10634/gaozy/aml_niche_net/results/tables"; SEED=491638

ap=argparse.ArgumentParser()
ap.add_argument("--root",default=DEFAULT_ROOT)
ap.add_argument("--n_perm",type=int,default=10000)
ap.add_argument("--alphas",default="0,0.25,0.5,0.75,1")
ap.add_argument("--max_iter",type=int,default=1000)
# node MASS mode. "ncells" = p_i proportional to cell count (canonical, feeds alpha_sweep.csv).
# "uniform" = equal mass on every present node -- removes the cell-count channel that even alpha=1
# retains through the masses, so alpha=1 + uniform is the fully cell-count-free topology test.
ap.add_argument("--mass_mode",choices=["ncells","uniform"],default="ncells")
ap.add_argument("--out",default="alpha_sweep.csv")
args=ap.parse_args()
# ONE INDEPENDENT PERMUTATION STREAM PER TEST, keyed on the test's stable identity (alpha, mass mode,
# model) and never on loop position. With a single shared module-level rng a given alpha's p depended
# on how many alphas ran before it: --alphas 0.5 and --alphas 0,0.25,0.5 gave p_strat 0.03598 vs
# 0.03998 with beta bit-identical. MASS_MODE is in the key so that the --mass_mode uniform run is a
# genuinely independent control rather than sharing the ncells run's Monte-Carlo error.
# Do NOT "fix" this by reseeding one rng inside the loop: that gives every test the SAME permutations,
# which is worse. And do not key on the loop INDEX -- inserting an alpha ahead of another would then
# silently re-roll every p below it.
def _lab(s): return int.from_bytes(hashlib.blake2b(s.encode(),digest_size=8).digest(),"big")
def _rng(alpha,mass_mode,model): return np.random.default_rng([SEED,_lab(f"{alpha:g}|{mass_mode}|{model}")])
MASS_MODE=args.mass_mode; EPS_MASS=1e-6
ALPHAS=[float(x) for x in args.alphas.split(",")]
D_FGW=os.path.join(args.root,"07_fgw"); D_OUT=os.path.join(args.root,"08_scoring"); os.makedirs(D_OUT,exist_ok=True)
_VOCAB = load_vocab(D_FGW)
# The FEATURE set is LOADED, not literal -- same reason as the timepoints above. Eight files
# carried the hard-coded triple while 01 wrote the real set into fgw_vocab.json; changing
# FGW_FEATURES in config_fgw.R silently left every one of them scoring the old model.
FGW_FEATURES = load_features(D_FGW)
AML_TP = _VOCAB["aml_timepoints"]


edges=pd.read_csv(os.path.join(D_FGW,"fgw_edges_long.csv"))
nodes=pd.read_csv(os.path.join(D_FGW,"fgw_nodes_long.csv"))
idx=pd.read_csv(os.path.join(D_FGW,"fgw_input_index.csv"))

# Fail loudly on a label the vocabulary does not cover. Without this an unrecognised
# timepoint is not an error -- it is a sample that quietly stops being AML and stops
# being healthy, and therefore stops existing for every test below.
assert_index_covered(idx, _VOCAB)

## -- per-sample inputs (cached; alpha-independent) ----
_cache={}
def build_one(ds,smp):
    k=(ds,smp)
    if k in _cache: return _cache[k]
    e=edges[(edges["dataset"]==ds)&(edges["sample"]==smp)]
    C=e.pivot(index="sender_bin",columns="receiver_bin",values="C").reindex(FGW_NODES,columns=FGW_NODES).to_numpy(float)
    C=np.nan_to_num(C,nan=1.0)
    nd=nodes[(nodes["dataset"]==ds)&(nodes["sample"]==smp)].set_index("hierarchy_bin").reindex(FGW_NODES)
    F=np.nan_to_num(nd[FGW_FEATURES].to_numpy(float),nan=0.0)
    if MASS_MODE=="uniform":                       # equal mass on present nodes (cell count removed)
        present=nd["present"].fillna(False).to_numpy(bool)
        p=np.where(present,1.0,EPS_MASS)
    else:                                          # canonical: p_i proportional to cell count
        p=np.nan_to_num(nd["mass"].to_numpy(float),nan=EPS_MASS)
    p=p/p.sum()
    _cache[k]=(C,F,p); return _cache[k]

# Silence ONLY the expected POT alpha=0 log-term 0/0, and only at the two call sites that raise it.
# Not process-wide: a blanket warnings.filterwarnings would also hide overflow or an invalid value in
# our own arithmetic, and the chain gates on unexpected warnings in this log.
@contextmanager
def pot_quiet():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered")
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered")
        yield

def barycenter(keys,alpha):
    Cs=[];Fs=[];ps=[]
    for k in keys:
        C,F,p=build_one(*k); Cs.append(C);Fs.append(F);ps.append(p)
    m=len(Cs); n=7
    with pot_quiet():
        out=ot.gromov.fgw_barycenters(n,Fs,Cs,ps,lambdas=[1.0/m]*m,alpha=alpha,loss_fun="square_loss",
            symmetric=False,max_iter=args.max_iter,p=np.ones(n)/n,init_C=np.mean(np.stack(Cs),0),
            init_X=np.mean(np.stack(Fs),0),random_state=SEED,log=True)
    return np.asarray(out[1]),np.asarray(out[0]),np.ones(n)/n

def fgw2(C,F,p,Cb,Fb,pb,alpha):
    M=ot.dist(F,Fb); mx=M.max()
    M=M/(mx+1e-9) if mx>0 else M
    with pot_quiet():
        return float(ot.gromov.fused_gromov_wasserstein2(M,C,Cb,p,pb,loss_fun="square_loss",
                                                         alpha=alpha,symmetric=False))

## -- cohort definition (identical to 08/01) ----
d=idx.copy()
d["grp"]=np.where(d.timepoint=="Healthy","Healthy",np.where(d.timepoint.isin(AML_TP),"AML","other"))
d=d[d.grp!="other"].copy()
d["is_aml"]=(d.grp=="AML").astype(float)
nc=nodes.pivot_table(index=["dataset","sample"],columns="hierarchy_bin",values="n_cells_raw",aggfunc="first").reset_index()
for b in FGW_NODES:
    if b not in nc: nc[b]=0
nc["total"]=nc[FGW_NODES].sum(axis=1)
nc["blast_proxy"]=nc[BLAST_BINS].sum(axis=1)/nc["total"].clip(lower=1)
d=d.merge(nc[["dataset","sample","blast_proxy"]],on=["dataset","sample"],how="left").dropna(subset=["blast_proxy"]).reset_index(drop=True)

heal_all=[(r.dataset,r["sample"]) for _,r in d[d.is_aml==0].iterrows()]
heal_bary=[(r.dataset,r["sample"]) for _,r in d[(d.is_aml==0)&(~d["sparse_flag"].fillna(False))].iterrows()]
aml_all =[(r.dataset,r["sample"]) for _,r in d[d.is_aml==1].iterrows()]
print(f"[1] healthy={len(heal_all)} (in barycenter: {len(heal_bary)}) | AML={len(aml_all)} | alphas={ALPHAS}")

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
        exceed+= (abs(float(ry@rp/dd))>=abs(beta))
    return beta,(1.0+exceed)/(n_perm+1.0)

rows=[]
for alpha in ALPHAS:
    tag=f"alpha={alpha:g}"
    try:
        Cb_full,Fb_full,pb=barycenter(heal_bary,alpha)
        hds={}
        for k in heal_all:                                   # leave-one-out for healthy
            loo=[x for x in heal_bary if x!=k]
            Cb,Fb,_=barycenter(loo,alpha)
            C,F,p=build_one(*k); hds[k]=fgw2(C,F,p,Cb,Fb,pb,alpha)
        for k in aml_all:                                    # AML never in the healthy barycenter
            C,F,p=build_one(*k); hds[k]=fgw2(C,F,p,Cb_full,Fb_full,pb,alpha)
    except Exception as ex:
        print(f"[!] {tag} FAILED: {ex!r}")
        rows.append(dict(alpha=alpha,status=f"failed: {ex}")); continue

    dd=d.copy()
    dd["HDS"]=[hds[(r.dataset,r["sample"])] for _,r in dd.iterrows()]
    y=dd["HDS"].to_numpy(float); a=dd["is_aml"].to_numpy(float)
    X0=np.column_stack([np.ones(len(dd)),dd["blast_proxy"].to_numpy(float)])
    bA,pA=fwl_perm(y,X0,a,None,args.n_perm,_rng(alpha,MASS_MODE,"global"))

    both=[ds for ds,g in dd.groupby("dataset") if g.is_aml.nunique()==2]
    db=dd[dd.dataset.isin(both)].reset_index(drop=True)
    yb=db["HDS"].to_numpy(float); ab=db["is_aml"].to_numpy(float)
    dums=pd.get_dummies(db["dataset"],drop_first=True).to_numpy(float)
    X0b=np.column_stack([np.ones(len(db)),db["blast_proxy"].to_numpy(float),dums])
    bB,pB=fwl_perm(yb,X0b,ab,db["dataset"].to_numpy(),args.n_perm,_rng(alpha,MASS_MODE,"strat"))

    rows.append(dict(alpha=alpha,status="ok",n_global=len(dd),n_strat=len(db),
                     mean_healthy=float(y[a==0].mean()),mean_aml=float(y[a==1].mean()),
                     beta_global=bA,p_global=pA,beta_strat=bB,p_strat=pB))
    print(f"[2] {tag:11s} | healthy={y[a==0].mean():.4f} AML={y[a==1].mean():.4f} "
          f"| global beta={bA:+.5f} p={pA:.5f} | within-dataset beta={bB:+.5f} p={pB:.5f}",flush=True)

res=pd.DataFrame(rows)
res.to_csv(os.path.join(D_OUT,args.out),index=False)
print("\n[3] ALPHA SWEEP (alpha=1 -> pure topology, features ignored; alpha=0 -> pure features):")
show=[c for c in ["alpha","mean_healthy","mean_aml","beta_global","p_global","beta_strat","p_strat"] if c in res.columns]
print(res[show].round(5).to_string(index=False))

print("\n[4] VERDICT")
ok=res[res.get("status","")=="ok"]
r1=ok[ok.alpha==1.0]
if len(r1):
    b1=float(r1.beta_strat.iloc[0]); p1=float(r1.p_strat.iloc[0])
    if p1<0.05 and b1>0:
        print(f"    alpha=1 (PURE GW, node features ignored): within-dataset beta={b1:+.5f} p={p1:.5f}")
        print("    -> H2 holds with ZERO contribution from node features. frac_malignant cannot be driving")
        print("       it (structurally excluded at alpha=1). The deviation is TOPOLOGICAL. Claim stands.")
    else:
        print(f"    alpha=1 (PURE GW): within-dataset beta={b1:+.5f} p={p1:.5f} -- NOT significant.")
        print("    -> the separation depends on the node-feature term; H2 must be restated as")
        print("       feature/composition-driven rather than topological.")
r0=ok[ok.alpha==0.0]
if len(r0):
    print(f"    alpha=0 (pure features) for contrast: within-dataset beta={float(r0.beta_strat.iloc[0]):+.5f} "
          f"p={float(r0.p_strat.iloc[0]):.5f}")
print(f"\n[done] wrote {os.path.join(D_OUT,args.out)}  (mass_mode={MASS_MODE})")
