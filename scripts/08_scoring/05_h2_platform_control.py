#!/usr/bin/env python
# 05_h2_platform_control.py ----
# Phase 7 (stage 08_scoring). CRITICAL VERIFICATION of H2 under the same scrutiny that killed H1.
#
# 04 showed that per-EDGE AML-vs-healthy effects vanish (and flip sign) once dataset/platform is
# controlled: 0/49 edges at q<0.05 in both the global and within-dataset models, with HSC_MPP->T_NK
# going from beta=-0.095 (global) to +0.004 (within-dataset). Those edge effects were platform-driven.
#
# H2 (08/01) was fitted as HDS_clean ~ blast_proxy + is_aml WITHOUT a dataset term, so the same doubt
# applies to it and MUST be resolved: if HDS is also platform-driven, the project's main claim falls.
# Here HDS_clean is run through the IDENTICAL two-model design used for edges in 04:
#   A "global"     : HDS ~ 1 + blast_proxy + is_aml, global permutation of is_aml (reproduces 08/01)
#   B "stratified" : restricted to datasets containing BOTH labels; HDS ~ 1 + blast_proxy +
#                    dataset_dummies + is_aml, permutation shuffles is_aml WITHIN dataset, so the
#                    is_aml effect is identified only from within-dataset contrasts (platform controlled
#                    by design).
#
# INTERPRETATION KEY:
#   B significant  -> aggregate topological deviation survives platform control while no individual edge
#                     does => the disease signal is DISTRIBUTED across the network, not localizable to a
#                     single edge (the H2 "topology, not single pairs" claim, now with a real control).
#   B null         -> H2 cannot be separated from platform in this cohort; the main claim must be
#                     restated as platform-confounded and downgraded. (Report honestly either way.)
#   HDS_clean is the leave-one-out corrected score from 08/01 (healthy samples scored against a healthy
#   barycenter that excludes them), so the circularity is already removed.
#
# INPUT  : <root>/08_scoring/h2_blast_regression.csv (HDS_clean, blast_proxy, is_aml, dataset, sample)
#          <root>/07_fgw/fgw_input_index.csv         (sparse_flag, if not already present)
# OUTPUT : <root>/08_scoring/h2_platform_control.csv
# Usage  : python scripts/08_scoring/05_h2_platform_control.py [--n_perm 10000] [--exclude_sparse]
import argparse, os
import numpy as np
import pandas as pd

DEFAULT_ROOT="/FAST/gr10634/gaozy/aml_niche_net/results/tables"; SEED=491638
ap=argparse.ArgumentParser()
ap.add_argument("--root",default=DEFAULT_ROOT)
ap.add_argument("--n_perm",type=int,default=10000)
ap.add_argument("--exclude_sparse",action="store_true",
                help="drop sparse-flagged graphs (default: keep, matching 08/01)")
args=ap.parse_args(); rng=np.random.default_rng(SEED)
D_OUT=os.path.join(args.root,"08_scoring"); D_FGW=os.path.join(args.root,"07_fgw")

d=pd.read_csv(os.path.join(D_OUT,"h2_blast_regression.csv"))
if "sparse_flag" not in d.columns:
    idx=pd.read_csv(os.path.join(D_FGW,"fgw_input_index.csv"))
    d=d.merge(idx[["dataset","sample","sparse_flag"]],on=["dataset","sample"],how="left")
if args.exclude_sparse: d=d[~d["sparse_flag"].fillna(False)]
d=d.dropna(subset=["HDS_clean","blast_proxy","is_aml"]).reset_index(drop=True)
print(f"[1] samples: healthy={int((d.is_aml==0).sum())} AML={int((d.is_aml==1).sum())} "
      f"(sparse {'excluded' if args.exclude_sparse else 'included'})")

def fwl_perm(y, X0, a, groups, n_perm, rng):
    """beta_a from y ~ [X0,a] == beta from resid(y|X0) ~ resid(a|X0). groups=None -> global shuffle."""
    Q,_=np.linalg.qr(X0)
    def resid(v): return v - Q@(Q.T@v)
    ry=resid(y); ra=resid(a); den=float(ra@ra)
    if den<=0: raise ValueError("is_aml has no residual variance after covariates")
    beta=float(ry@ra/den)
    if groups is None:
        gidx=None
    else:
        gidx=[np.where(groups==g)[0] for g in np.unique(groups)]
    exceed=0
    for _ in range(n_perm):
        if gidx is None:
            ap_=rng.permutation(a)
        else:
            ap_=a.copy()
            for gi in gidx: ap_[gi]=rng.permutation(ap_[gi])
        rp=resid(ap_); dd=float(rp@rp)
        if dd<=0: continue
        b=float(ry@rp/dd)
        exceed+= (abs(b)>=abs(beta))
    return beta,(1.0+exceed)/(n_perm+1.0)

rows=[]
## -- Model A: global (reproduces 08/01, now with a two-sided permutation p) ----
y=d["HDS_clean"].to_numpy(float); blast=d["blast_proxy"].to_numpy(float); a=d["is_aml"].to_numpy(float)
X0=np.column_stack([np.ones(len(d)),blast])
bA,pA=fwl_perm(y,X0,a,None,args.n_perm,rng)
print(f"\n[2] model A (global, n={len(d)}, blast-controlled):")
print(f"    beta[is_aml]={bA:+.5f}  permutation p={pA:.5f}")
print(f"    mean HDS_clean: healthy={y[a==0].mean():.4f}  AML={y[a==1].mean():.4f}")
rows.append(dict(model="A_global",n=len(d),beta_is_aml=bA,perm_p=pA,
                 mean_healthy=float(y[a==0].mean()),mean_aml=float(y[a==1].mean())))

## -- Model B: within-dataset (platform controlled by design) ----
both=[ds for ds,g in d.groupby("dataset") if g.is_aml.nunique()==2]
db=d[d.dataset.isin(both)].reset_index(drop=True)
print(f"\n[3] model B datasets with both labels: {both}")
print(f"    n={len(db)} (healthy={int((db.is_aml==0).sum())}, AML={int((db.is_aml==1).sum())})")
yb=db["HDS_clean"].to_numpy(float); ab=db["is_aml"].to_numpy(float)
dums=pd.get_dummies(db["dataset"],drop_first=True).to_numpy(float)
X0b=np.column_stack([np.ones(len(db)),db["blast_proxy"].to_numpy(float),dums])
bB,pB=fwl_perm(yb,X0b,ab,db["dataset"].to_numpy(),args.n_perm,rng)
print(f"    beta[is_aml]={bB:+.5f}  permutation p={pB:.5f}  (within-dataset contrasts only)")
rows.append(dict(model="B_within_dataset",n=len(db),beta_is_aml=bB,perm_p=pB,
                 mean_healthy=float(yb[ab==0].mean()),mean_aml=float(yb[ab==1].mean())))

## -- per-dataset breakdown (descriptive, for the figure) ----
print("\n[4] per-dataset HDS_clean (datasets with both labels):")
for ds in both:
    g=db[db.dataset==ds]
    print(f"    {ds:14s} healthy n={int((g.is_aml==0).sum()):2d} mean={g[g.is_aml==0].HDS_clean.mean():.4f} | "
          f"AML n={int((g.is_aml==1).sum()):2d} mean={g[g.is_aml==1].HDS_clean.mean():.4f}")
    rows.append(dict(model=f"desc_{ds}",n=len(g),beta_is_aml=np.nan,perm_p=np.nan,
                     mean_healthy=float(g[g.is_aml==0].HDS_clean.mean()),
                     mean_aml=float(g[g.is_aml==1].HDS_clean.mean())))

pd.DataFrame(rows).to_csv(os.path.join(D_OUT,"h2_platform_control.csv"),index=False)
print("\n[5] VERDICT")
print(f"    aggregate HDS: global p={pA:.5f} | within-dataset p={pB:.5f}")
print(f"    individual edges (from 04): 0/49 at q<0.05 in BOTH models")
if pB<0.05:
    print("    -> H2 SURVIVES platform control while no single edge does: the signal is DISTRIBUTED")
    print("       across the network rather than localizable to any one edge.")
else:
    print("    -> H2 does NOT separate from platform in the within-dataset test; the main claim must be")
    print("       restated with this limitation (note model B has much lower n, so this is not proof of absence).")
print(f"\n[done] wrote {os.path.join(D_OUT,'h2_platform_control.csv')}")
