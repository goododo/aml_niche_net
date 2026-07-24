#!/usr/bin/env python
# 03_block_permutation.py ----
# Phase 7 (stage 08_scoring). BLOCK-level permutation test, the principled follow-up to 02.
#
# WHY: 02 tested all 49 directed edges individually -> 0/49 significant, with a p-distribution
# indistinguishable from uniform. That is a RESOLUTION/POWER problem, not necessarily absence of signal:
# the effect-size pattern (HSC_MPP-sender edges up, immune-immune edges down) is BLOCK-shaped, and the
# per-edge null is inflated by cross-dataset heterogeneity + a small healthy group (n=20) + 49 tests.
# Aggregating pre-specified blocks cuts the tests 49 -> 5 and pools signal across edges.
#
# PRE-SPECIFICATION (critical, this is NOT post-hoc fishing): the blocks are derived from constants that
# were locked in config_ccc.R BEFORE any of these results existed --
#   CCC_SENDER_LSC  = {HSC_MPP, LMPP_GMP}          (the SCS sender set)
#   CCC_RECV_IMMUNE = {T_NK, B_Plasma, Mono_DC}    (the SCS receiver set)
# The primary block (primitive->immune) IS the blueprint's SCS axis. Report as a pre-specified block test.
#
# STATISTIC: mean dC over the edges of a block, where dC = C_AML_barycenter - C_healthy_barycenter.
#   dC < 0 => that block is STRONGER in AML (C is rank-distance: lower = stronger).
# NULL: identical global label shuffle as 02 (healthy count fixed), rebuilding both barycenters.
#
# INTERPRETIVE CAVEAT (stated in the output): C is a WITHIN-SAMPLE rank-percentile, so the 49 edges are
# zero-sum -- blocks are compositionally COUPLED (one block rising forces others to fall). The permutation
# null carries the same constraint, so the test is calibrated; but results must be read as a
# RE-PRIORITIZATION of the communication hierarchy, never as absolute gain/loss of signalling.
#
# INPUT  : <root>/07_fgw/fgw_{nodes,edges,input_index}_long.csv
# OUTPUT : <root>/08_scoring/block_permutation.csv  (per block: n_edges, mean_dC, direction, perm_p, q_bh)
#          <root>/08_scoring/block_null.npz         (null distributions, for plotting)
# Usage  : python scripts/08_scoring/03_block_permutation.py [--n_perm 10000] [--max_iter 300]
import argparse, os
import numpy as np
import pandas as pd
import ot

FGW_NODES=["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
FGW_FEATURES=["frac_malignant","mean_stemness","n_cells"]
AML_TP={"Diagnosis","MRD","Post_treatment","Relapse","Relapse2"}
DEFAULT_ROOT="/FAST/gr10634/gaozy/aml_niche_net/results/tables"; SEED=491638

# --- pre-specified sets (mirror config_ccc.R; locked before results existed) ---
PRIMITIVE=["HSC_MPP","LMPP_GMP"]          # CCC_SENDER_LSC
IMMUNE   =["T_NK","B_Plasma","Mono_DC"]   # CCC_RECV_IMMUNE
OTHER    =["Erythroid","Megakaryocyte"]

ap=argparse.ArgumentParser()
ap.add_argument("--root",default=DEFAULT_ROOT); ap.add_argument("--alpha",type=float,default=0.5)
ap.add_argument("--n_perm",type=int,default=10000); ap.add_argument("--max_iter",type=int,default=300)
args=ap.parse_args(); rng=np.random.default_rng(SEED)
D_FGW=os.path.join(args.root,"07_fgw"); D_OUT=os.path.join(args.root,"08_scoring"); os.makedirs(D_OUT,exist_ok=True)

edges=pd.read_csv(os.path.join(D_FGW,"fgw_edges_long.csv"))
nodes=pd.read_csv(os.path.join(D_FGW,"fgw_nodes_long.csv"))
idx=pd.read_csv(os.path.join(D_FGW,"fgw_input_index.csv"))

## -- block masks over the 7x7 directed grid ----
ii={n:i for i,n in enumerate(FGW_NODES)}
def mask(senders,receivers):
    m=np.zeros((7,7),dtype=bool)
    for s in senders:
        for r in receivers: m[ii[s],ii[r]]=True
    return m
BLOCKS={
 "primitive_to_immune": mask(PRIMITIVE,IMMUNE),   # PRIMARY: the SCS axis (LSC -> microenvironment)
 "immune_to_primitive": mask(IMMUNE,PRIMITIVE),
 "immune_to_immune"   : mask(IMMUNE,IMMUNE),      # incl. self-loops
 "primitive_autocrine": mask(PRIMITIVE,PRIMITIVE),# incl. self-loops
}
used=np.zeros((7,7),dtype=bool)
for m in BLOCKS.values(): used|=m
BLOCKS["other_lineages"]=~used                    # erythroid / megakaryocyte involvement
assert sum(m.sum() for m in BLOCKS.values())==49, "blocks must partition the 49 directed edges"
print("[0] blocks (edges):", {k:int(v.sum()) for k,v in BLOCKS.items()})

## -- per-sample inputs (cached) ----
Ccache={};Fcache={};pcache={}
def build_one(ds,smp):
    k=(ds,smp)
    if k in Ccache: return Ccache[k],Fcache[k],pcache[k]
    e=edges[(edges["dataset"]==ds)&(edges["sample"]==smp)]
    C=e.pivot(index="sender_bin",columns="receiver_bin",values="C").reindex(FGW_NODES,columns=FGW_NODES).to_numpy(float)
    C=np.nan_to_num(C,nan=1.0)
    nd=nodes[(nodes["dataset"]==ds)&(nodes["sample"]==smp)].set_index("hierarchy_bin").reindex(FGW_NODES)
    F=np.nan_to_num(nd[FGW_FEATURES].to_numpy(float),nan=0.0)
    p=np.nan_to_num(nd["mass"].to_numpy(float),nan=1e-6); p=p/p.sum()
    Ccache[k],Fcache[k],pcache[k]=C,F,p; return C,F,p

def bary_C(keys):
    Cs=[];Fs=[];ps=[]
    for k in keys:
        C,F,p=build_one(*k); Cs.append(C);Fs.append(F);ps.append(p)
    m=len(Cs); n=7
    out=ot.gromov.fgw_barycenters(n,Fs,Cs,ps,lambdas=[1.0/m]*m,alpha=args.alpha,loss_fun="square_loss",
        symmetric=False,max_iter=args.max_iter,p=np.ones(n)/n,init_C=np.mean(np.stack(Cs),0),
        init_X=np.mean(np.stack(Fs),0),random_state=SEED,log=True)
    return np.asarray(out[1])

d=idx.copy()
d["grp"]=np.where(d.timepoint=="Healthy","Healthy",np.where(d.timepoint.isin(AML_TP),"AML","other"))
d=d[d.grp!="other"].copy(); d["nonsparse"]=~d["sparse_flag"].fillna(False)
dn=d[d.nonsparse]
pool=[(r.dataset,r["sample"]) for _,r in dn.iterrows()]
is_h=np.array([r.grp=="Healthy" for _,r in dn.iterrows()])
n_h=int(is_h.sum())
heal_keys=[pool[i] for i in range(len(pool)) if is_h[i]]
aml_keys =[pool[i] for i in range(len(pool)) if not is_h[i]]
print(f"[1] pool: healthy={len(heal_keys)} AML={len(aml_keys)} | n_perm={args.n_perm}")

## -- real block statistics ----
dC_real=bary_C(aml_keys)-bary_C(heal_keys)
names=list(BLOCKS.keys())
real_stat={b: float(dC_real[BLOCKS[b]].mean()) for b in names}
print("[1] real block mean dC (negative => block STRONGER in AML):")
for b in names: print(f"    {b:22s} n_edges={int(BLOCKS[b].sum()):2d}  mean_dC={real_stat[b]:+.5f}")

## -- permutation null on block means ----
print(f"[2] permuting (rebuilds 2 barycenters/perm)...")
null=np.zeros((args.n_perm,len(names)))
step=max(1,args.n_perm//20)
for it in range(args.n_perm):
    pi=rng.permutation(len(pool))
    ph=[pool[j] for j in pi[:n_h]]; pa=[pool[j] for j in pi[n_h:]]
    dC=bary_C(pa)-bary_C(ph)
    for bi,b in enumerate(names): null[it,bi]=dC[BLOCKS[b]].mean()
    if (it+1)%step==0: print(f"    ... {it+1}/{args.n_perm}",flush=True)

## -- two-sided p + BH-FDR over the 5 blocks ----
rows=[]
for bi,b in enumerate(names):
    obs=real_stat[b]
    p=(1.0+float((np.abs(null[:,bi])>=abs(obs)).sum()))/(args.n_perm+1.0)
    rows.append(dict(block=b,n_edges=int(BLOCKS[b].sum()),mean_dC=obs,
                     direction=("stronger_in_AML" if obs<0 else "weaker_in_AML"),
                     null_sd=float(null[:,bi].std()), perm_p=p))
res=pd.DataFrame(rows).sort_values("perm_p").reset_index(drop=True)
m=len(res)
res["q_bh"]=(res["perm_p"]*m/np.arange(1,m+1)).clip(upper=1.0)
res["q_bh"]=res["q_bh"][::-1].cummin()[::-1]
res["significant_fdr"]=res["q_bh"]<0.05
res.to_csv(os.path.join(D_OUT,"block_permutation.csv"),index=False)
np.savez(os.path.join(D_OUT,"block_null.npz"),null=null,names=np.array(names),
         real=np.array([real_stat[b] for b in names]))

print("\n[3] BLOCK RESULTS (BH-FDR over 5 blocks):")
print(res.round(5).to_string(index=False))
print(f"\n[3] FDR-significant blocks: {int(res.significant_fdr.sum())}/5")
print("[3] CAVEAT: C is a within-sample rank -> the 49 edges are zero-sum -> blocks are COUPLED.")
print("    Read results as RE-PRIORITIZATION of the communication hierarchy, not absolute gain/loss.")
print(f"\n[done] wrote {os.path.join(D_OUT,'block_permutation.csv')} + block_null.npz")
