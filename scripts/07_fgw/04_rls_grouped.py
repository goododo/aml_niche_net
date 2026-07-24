#!/usr/bin/env python
# 04_rls_grouped.py ----
# Phase 6 (stage 07_fgw). GROUPED-barycenter RLS = blueprint's original form (A), run ONLY as a
# double-confirmation of the H3 negative result already seen with the paired within-patient form (B, 03).
# GSE227903 is single-platform (10x) so the cross-patient barycenter does NOT reintroduce platform
# confound here (the one place A is defensible). Expected: still non-significant + near-tautological
# ("a relapse graph sits closer to the relapse barycenter"). If A ALSO fails -> H3 negative is robust to
# the RLS definition (reported as such). NO further investment beyond this one run.
#
#   B_Dx      = FGW barycenter of GSE227903 Diagnosis graphs
#   B_Relapse = FGW barycenter of GSE227903 Relapse graphs
#   RLS_grp(sample) = FGW2(G, B_Dx) - FGW2(G, B_Relapse)   ( >0 => closer to relapse consensus )
# Tested on the 8 paired patients' Relapse graphs (sign of RLS_grp; expect >0 if relapse converges to
# a relapse-like consensus). Small n, single dataset -> descriptive/confirmatory only.
#
# INPUT  : <root>/07_fgw/fgw_{nodes,edges,input_index}_long.csv
# OUTPUT : <root>/07_fgw/rls_grouped.csv  (per relapse sample: d_to_Bdx, d_to_Brel, RLS_grp)
# Usage  : python scripts/07_fgw/04_rls_grouped.py [--root .../results/tables] [--alpha 0.5]
import argparse, os
import numpy as np
import pandas as pd
import ot

FGW_NODES = ["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
FGW_FEATURES = ["frac_malignant","mean_stemness","n_cells"]
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
DS = "GSE227903"; SEED = 491638

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--alpha", type=float, default=0.5)
args = ap.parse_args()
D_FGW = os.path.join(args.root, "07_fgw")

edges = pd.read_csv(os.path.join(D_FGW, "fgw_edges_long.csv"))
nodes = pd.read_csv(os.path.join(D_FGW, "fgw_nodes_long.csv"))
idx   = pd.read_csv(os.path.join(D_FGW, "fgw_input_index.csv"))

def build_CFp(smp):
    e = edges[(edges["dataset"]==DS) & (edges["sample"]==smp)]
    C = e.pivot(index="sender_bin", columns="receiver_bin", values="C").reindex(FGW_NODES, columns=FGW_NODES).to_numpy(float)
    C = np.nan_to_num(C, nan=1.0)
    nd = nodes[(nodes["dataset"]==DS) & (nodes["sample"]==smp)].set_index("hierarchy_bin").reindex(FGW_NODES)
    F = np.nan_to_num(nd[FGW_FEATURES].to_numpy(float), nan=0.0)
    p = np.nan_to_num(nd["mass"].to_numpy(float), nan=1e-6); p = p/p.sum()
    return C, F, p

def barycenter(smps):
    Cs=[]; Fs=[]; ps=[]
    for s in smps:
        C,F,p = build_CFp(s); Cs.append(C); Fs.append(F); ps.append(p)
    m=len(Cs); n=len(FGW_NODES)
    out = ot.gromov.fgw_barycenters(n, Fs, Cs, ps, lambdas=[1.0/m]*m, alpha=args.alpha,
                                    loss_fun="square_loss", symmetric=False, max_iter=1000,
                                    p=np.ones(n)/n, init_C=np.mean(np.stack(Cs),0), init_X=np.mean(np.stack(Fs),0),
                                    random_state=SEED, log=True)
    return np.asarray(out[1]), np.asarray(out[0]), np.ones(n)/n   # C, F, p

def fgw2(C, F, p, Cb, Fb, pb):
    M = ot.dist(F, Fb); M = M/(M.max()+1e-9)
    return float(ot.gromov.fused_gromov_wasserstein2(M, C, Cb, p, pb, loss_fun="square_loss", alpha=args.alpha, symmetric=False))

g = idx[idx["dataset"]==DS].copy()
g["patient"] = g["sample"].str.extract(r'^(\d+)_')
g["tp"] = g["sample"].str.extract(r'_(Dg|MRD|R2|R)$')
dg_samples  = g[g.tp=="Dg"]["sample"].tolist()
rel_samples = g[g.tp=="R"]["sample"].tolist()
print(f"[1] B_Dx from {len(dg_samples)} Dg graphs | B_Relapse from {len(rel_samples)} Relapse graphs")

C_dx, F_dx, p_dx    = barycenter(dg_samples)
C_rel, F_rel, p_rel = barycenter(rel_samples)
print(f"[1] barycenters built | asym B_Dx={np.abs(C_dx-C_dx.T).max():.3f} B_Rel={np.abs(C_rel-C_rel.T).max():.3f}")

# paired patients (Dg+R) -> score their Relapse graph
present = g.groupby("patient")["tp"].apply(set)
pairs = sorted(p for p,t in present.items() if "Dg" in t and "R" in t)
rows=[]
for pt in pairs:
    s = f"{pt}_R"
    C,F,p = build_CFp(s)
    d_dx  = fgw2(C,F,p, C_dx,  F_dx,  p_dx)
    d_rel = fgw2(C,F,p, C_rel, F_rel, p_rel)
    rows.append(dict(patient=pt, sample=s, d_to_Bdx=d_dx, d_to_Brel=d_rel, RLS_grp=d_dx-d_rel))
res = pd.DataFrame(rows)
res.to_csv(os.path.join(D_FGW, "rls_grouped.csv"), index=False)
print("[2] per relapse sample (RLS_grp>0 => closer to Relapse consensus):")
print(res.round(4).to_string(index=False))

x = res.RLS_grp.to_numpy()
n_up = int((x>0).sum()); n=len(x)
p_w = "NA(no scipy)"
try:
    from scipy.stats import wilcoxon
    p_w = round(float(wilcoxon(x, alternative="greater").pvalue), 4)
except Exception: pass
print(f"\n[3] RLS_grp > 0 (relapse closer to relapse consensus): {n_up}/{n} up | median={np.median(x):.4f} | wilcoxon_p(greater)={p_w}")
print("[3] interpretation: even a POSITIVE result here is partly tautological (relapse graph vs relapse")
print("    barycenter it helped build). Compare to 03 (paired form B). If both non-sig -> H3 robustly null.")
print(f"[done] wrote {os.path.join(D_FGW,'rls_grouped.csv')}")
