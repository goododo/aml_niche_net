#!/usr/bin/env python
# 02_permutation_emergent.py ----
# Phase 7 (stage 08_scoring). Emergent edge identification + significance (scheme X):
#   real dC = C_AML_barycenter - C_healthy_barycenter   (49 directed edges; C=rank-distance, LOWER=stronger)
#   dC<0 => edge STRONGER in AML;  dC>0 => WEAKER in AML.
#
# (1) PERMUTATION NULL (global): shuffle healthy/AML labels among all non-sparse healthy+AML samples
#     (keep the healthy count fixed), rebuild BOTH barycenters, recompute dC^perm. Per-edge two-sided
#     permutation p = P(|dC^perm| >= |dC_real|). Full-cohort barycenter -> full statistical power.
# (2) CROSS-DATASET REPRODUCIBILITY (scheme b = H1 conserved-axis test): for each edge, every AML
#     dataset votes its direction = sign( mean_over_its_AML(C[edge]) - C_healthy_real[edge] ), using the
#     COMMON healthy barycenter as reference. An edge is CONSERVED if >=2 datasets agree with sign(dC_real).
#     This replaces stratified shuffling (which would have shrunk the barycenter to 3 datasets); it is
#     stricter (requires multi-cohort agreement) AND is the H1 evidence, without losing barycenter power.
#
# EMERGENT = permutation-significant (p<0.05) AND cross-dataset-reproducible (>=2 concordant). BH-FDR q also reported.
#
# INPUT  : <root>/07_fgw/fgw_{nodes,edges,input_index}_long.csv
# OUTPUT : <root>/08_scoring/emergent_edges.csv (per edge: dC_real, perm_p, q, n_concordant, is_emergent, direction)
# Usage  : python scripts/08_scoring/02_permutation_emergent.py [--root ...] [--alpha 0.5] [--n_perm 500] [--max_iter 500]
#          FIRST run with --n_perm 20 to time it (each perm rebuilds a ~93-graph barycenter).
import argparse, os
import sys
import numpy as np
import pandas as pd
import ot

# The timepoint vocabulary is LOADED, not literal. A hard-coded set here silently deleted 34 of
# 214 samples (64% of the treated arm) after CANONICAL_TIMEPOINTS changed on 2026-08-04.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config'))
from fgw_vocab import load_vocab, load_features, assert_index_covered


FGW_NODES=["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
AML_TP = None  # set from fgw_vocab.json below -- see scripts/config/fgw_vocab.py
DEFAULT_ROOT="/FAST/gr10634/gaozy/aml_niche_net/results/tables"; SEED=491638
MIN_DS_AML=3   # a dataset needs >= this many non-sparse AML samples to cast a reproducibility vote

ap=argparse.ArgumentParser()
ap.add_argument("--root",default=DEFAULT_ROOT); ap.add_argument("--alpha",type=float,default=0.5)
ap.add_argument("--n_perm",type=int,default=500); ap.add_argument("--max_iter",type=int,default=500)
args=ap.parse_args(); rng=np.random.default_rng(SEED)
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

# per-sample C matrix cache (7x7, node-ordered)
Ccache={}; Fcache={}; pcache={}
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
    m=len(Cs); n=len(FGW_NODES)
    out=ot.gromov.fgw_barycenters(n,Fs,Cs,ps,lambdas=[1.0/m]*m,alpha=args.alpha,loss_fun="square_loss",
        symmetric=False,max_iter=args.max_iter,p=np.ones(n)/n,init_C=np.mean(np.stack(Cs),0),
        init_X=np.mean(np.stack(Fs),0),random_state=SEED,log=True)
    return np.asarray(out[1])   # structure C only

## -- label pools (barycenter uses NON-SPARSE only, matching 02) ----
d=idx.copy()
d["grp"]=np.where(d.timepoint=="Healthy","Healthy",np.where(d.timepoint.isin(AML_TP),"AML","other"))
d=d[d.grp!="other"].copy()
d["nonsparse"]=~d["sparse_flag"].fillna(False)
pool=[(r.dataset,r["sample"]) for _,r in d[d.nonsparse].iterrows()]
labels=np.array([1 if r.grp=="Healthy" else 0 for _,r in d[d.nonsparse].iterrows()])  # 1=healthy
n_h=int(labels.sum())
heal_keys=[pool[i] for i in range(len(pool)) if labels[i]==1]
aml_keys =[pool[i] for i in range(len(pool)) if labels[i]==0]
print(f"[1] barycenter pool (non-sparse): healthy={len(heal_keys)} AML={len(aml_keys)} | shuffling {len(pool)} labels, n_healthy fixed={n_h}")

## -- real dC ----
C_h=bary_C(heal_keys); C_a=bary_C(aml_keys)
dC_real=C_a - C_h                 # 7x7
print(f"[1] real dC computed | max|dC|={np.abs(dC_real).max():.4f} | edges with dC<0 (stronger in AML): {int((dC_real<0).sum())}/49")

## -- permutation null ----
print(f"[2] permutation null: n_perm={args.n_perm}, max_iter={args.max_iter} (rebuilds 2 barycenters/perm)")
exceed=np.zeros((7,7)); pool_arr=np.array(pool,dtype=object)
for it in range(args.n_perm):
    pi=rng.permutation(len(pool))
    ph=[pool[j] for j in pi[:n_h]]; pa=[pool[j] for j in pi[n_h:]]
    dC=bary_C(pa)-bary_C(ph)          # pseudo-AML - pseudo-healthy, mirroring real (C_a - C_h)
    exceed += (np.abs(dC) >= np.abs(dC_real)).astype(float)
    if (it+1)%25==0: print(f"    ... {it+1}/{args.n_perm}")
perm_p=(1.0+exceed)/(args.n_perm+1.0)

## -- cross-dataset reproducibility (scheme b: common healthy reference) ----
# per (dataset, edge): mean C over that dataset's non-sparse AML samples
aml_ds=d[(d.grp=="AML")&(d.nonsparse)].copy()
ds_edge_meanC={}  # dataset -> 7x7 mean C
for ds,grp in aml_ds.groupby("dataset"):
    keys=[(ds,s) for s in grp["sample"]]
    if len(keys)<MIN_DS_AML: continue
    ds_edge_meanC[ds]=np.mean(np.stack([build_one(*k)[0] for k in keys]),axis=0)
print(f"[3] reproducibility voting datasets (>= {MIN_DS_AML} non-sparse AML): {sorted(ds_edge_meanC.keys())}")

sign_real=np.sign(dC_real)          # sign of (C_aml - C_healthy)
concordant=np.zeros((7,7),dtype=int)
for ds,Cm in ds_edge_meanC.items():
    sign_ds=np.sign(Cm - C_h)       # dataset AML mean C - common healthy barycenter C
    concordant += (sign_ds==sign_real).astype(int)

## -- assemble + BH-FDR ----
rows=[]
for i,si in enumerate(FGW_NODES):
    for j,rj in enumerate(FGW_NODES):
        rows.append(dict(sender_bin=si,receiver_bin=rj,dC_real=dC_real[i,j],
                         direction=("stronger_in_AML" if dC_real[i,j]<0 else "weaker_in_AML"),
                         perm_p=perm_p[i,j],n_concordant=int(concordant[i,j])))
res=pd.DataFrame(rows)
# BH-FDR
order=res["perm_p"].rank(method="first").astype(int); m=len(res)
res=res.sort_values("perm_p").reset_index(drop=True)
res["q_bh"]=(res["perm_p"]*m/(np.arange(1,m+1))).clip(upper=1.0)
res["q_bh"]=res["q_bh"][::-1].cummin()[::-1]
# THE VERDICT USES q_bh, NOT perm_p. 49 edges are tested at once; the BH correction is computed two
# lines up and used to be discarded here, so "emergent" meant "uncorrected p < 0.05" -- a label that
# 2.45 edges earn by chance. That is not hypothetical: the count on disk went 5 -> 1 -> 5 across three
# runs of this script while the corrected answer was 0 every time. Keep the uncorrected column so the
# nominal count stays visible, but do not let it name anything.
res["is_nominal"]=(res["perm_p"]<0.05)&(res["n_concordant"]>=2)
res["is_emergent"]=(res["q_bh"]<0.05)&(res["n_concordant"]>=2)
res.to_csv(os.path.join(D_OUT,"emergent_edges.csv"),index=False)

n_sig=int((res.perm_p<0.05).sum()); n_nom=int(res.is_nominal.sum()); n_emg=int(res.is_emergent.sum())
print(f"\n[4] edges of {len(res)}: nominally significant (uncorrected p<0.05) = {n_sig}"
      f"  [expected by chance at this cutoff: {0.05*len(res):.1f}]")
print(f"[4]   ...of those, cross-dataset-concordant(>=2) = {n_nom}   <- NOMINAL, not a result")
print(f"[4] EMERGENT (BH q<0.05 AND >=2 datasets concordant) = {n_emg}")
cols=["sender_bin","receiver_bin","dC_real","direction","perm_p","q_bh","n_concordant"]
em=res[res.is_emergent].sort_values("dC_real")
print(em[cols].round(4).to_string(index=False) if len(em) else "    (none survive correction)")
if n_emg==0 and n_nom>0:
    print(f"\n[4] The {n_nom} nominal edge(s), shown ONLY so the gap between the two counts is visible.")
    print("    They do not survive multiple-testing correction and must not be reported as findings:")
    print(res[res.is_nominal].sort_values("perm_p")[cols].round(4).to_string(index=False))
    print(f"    smallest q_bh over all {len(res)} edges = {res.q_bh.min():.4f}")
print(f"\n[done] wrote {os.path.join(D_OUT,'emergent_edges.csv')}")
