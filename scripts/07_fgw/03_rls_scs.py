#!/usr/bin/env python
# 03_rls_scs.py ----
# Phase 6 (stage 07_fgw). WITHIN-PATIENT paired analysis on GSE227903 Dg+Relapse pairs (n=8). Paired
# design => same patient, same platform => platform confound cancels (the cleanest H3 test; blueprint M1
# keeps relapse as an EXPLORATORY endpoint, so n=8 is reported honestly).
#
# Per patient (Diagnosis vs Relapse):
#   RLS   = FGW2(G_Dx, G_Relapse)            -- magnitude the communication topology moved (directed, alpha main)
#   dHDS  = HDS_Relapse - HDS_Dx             -- did relapse move AWAY from healthy?  (H3: > 0)
#   dATS  = ATS_Relapse - ATS_Dx             -- did relapse move TOWARD the AML barycenter? (expect < 0)
#   SCS   = sum_{i in primitive} frac_mal_i(RAW) * sum_{j in immune} weight_probsum[i,j]
#           dSCS = SCS_Relapse - SCS_Dx      -- LSC->immune malignant-weighted communication change (H: > 0)
#
# [locked] SCS weight = RAW frac_malignant [decision (i), conservative: inferCNV under-calls quiescent
#   relapse LSCs -> UNDER-estimates relapse SCS -> biases dSCS toward null, never false-positive].
#   SCS uses RAW edge strength (weight_probsum), NOT rank-distance C (C is inverted/relative).
#   RLS/dHDS/dATS use the SAME z-scored F + rank C + mass as 02 (consistency with the barycenters).
# 4978 has a sparse-flagged Relapse graph (8 edges) -> flagged; tests reported WITH and WITHOUT it.
#
# INPUT  : <root>/06_distance/edge_distance.csv (C + weight_probsum) , <root>/07_fgw/fgw_nodes_long.csv
#          (z F + mass) , <root>/05_ccc/ccc_node_features.csv (RAW frac_malignant) ,
#          <root>/07_fgw/fgw_input_index.csv (sparse_flag) , <root>/07_fgw/patient_scores.csv (HDS/ATS)
# OUTPUT : <root>/07_fgw/paired_rls_scs.csv  (per patient: RLS, HDS/ATS + deltas, SCS + delta, sparse)
# Usage  : python scripts/07_fgw/03_rls_scs.py [--root .../results/tables] [--alpha 0.5]
import argparse, os
import numpy as np
import pandas as pd
import ot

FGW_NODES = ["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
FGW_FEATURES = ["frac_malignant","mean_stemness","n_cells"]
PRIMITIVE = ["HSC_MPP","LMPP_GMP"]
IMMUNE    = ["T_NK","B_Plasma","Mono_DC"]
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--alpha", type=float, default=0.5)
args = ap.parse_args()
D_FGW = os.path.join(args.root, "07_fgw"); D_DIST = os.path.join(args.root, "06_distance"); D_CCC = os.path.join(args.root, "05_ccc")

ed    = pd.read_csv(os.path.join(D_DIST, "edge_distance.csv"))        # C + weight_probsum
nodes = pd.read_csv(os.path.join(D_FGW, "fgw_nodes_long.csv"))        # z F + mass
nfraw = pd.read_csv(os.path.join(D_CCC, "ccc_node_features.csv"))     # RAW frac_malignant
idx   = pd.read_csv(os.path.join(D_FGW, "fgw_input_index.csv"))
sc    = pd.read_csv(os.path.join(D_FGW, "patient_scores.csv"))

DS = "GSE227903"
def build_CFp(smp):
    e = ed[(ed["dataset"]==DS) & (ed["sample"]==smp)]
    C = e.pivot(index="sender_bin", columns="receiver_bin", values="C").reindex(FGW_NODES, columns=FGW_NODES).to_numpy(float)
    C = np.nan_to_num(C, nan=1.0)
    nd = nodes[(nodes["dataset"]==DS) & (nodes["sample"]==smp)].set_index("hierarchy_bin").reindex(FGW_NODES)
    F = np.nan_to_num(nd[FGW_FEATURES].to_numpy(float), nan=0.0)
    p = np.nan_to_num(nd["mass"].to_numpy(float), nan=1e-6); p = p/p.sum()
    return C, F, p

def scs(smp):
    e = ed[(ed["dataset"]==DS) & (ed["sample"]==smp)]
    W = e.pivot(index="sender_bin", columns="receiver_bin", values="weight_probsum").reindex(FGW_NODES, columns=FGW_NODES).to_numpy(float)
    W = np.nan_to_num(W, nan=0.0)
    fr = nfraw[(nfraw["dataset"]==DS) & (nfraw["sample"]==smp)].set_index("hierarchy_bin").reindex(FGW_NODES)["frac_malignant"].to_numpy(float)
    fr = np.nan_to_num(fr, nan=0.0)
    total = 0.0
    for i in PRIMITIVE:
        wi = fr[FGW_NODES.index(i)]
        for j in IMMUNE:
            total += wi * W[FGW_NODES.index(i), FGW_NODES.index(j)]
    return float(total)

def fgw2(smpA, smpB):
    Ca,Fa,pa = build_CFp(smpA); Cb,Fb,pb = build_CFp(smpB)
    M = ot.dist(Fa, Fb); M = M/(M.max()+1e-9)
    return float(ot.gromov.fused_gromov_wasserstein2(M, Ca, Cb, pa, pb, loss_fun="square_loss", alpha=args.alpha, symmetric=False))

def hdsats(smp):
    r = sc[(sc["dataset"]==DS) & (sc["sample"]==smp)]
    return (float(r["HDS"].iloc[0]), float(r["ATS"].iloc[0])) if len(r) else (np.nan, np.nan)

sparse_set = set(idx[(idx["dataset"]==DS) & (idx["sparse_flag"]==True)]["sample"])

# GSE227903 patients with both Dg and R
g = sc[sc["dataset"]==DS].copy()
g["patient"] = g["sample"].str.extract(r'^(\d+)_')
present = g.groupby("patient")["sample"].apply(lambda s: set(x.split("_",1)[1] for x in s))
pairs = [p for p,tps in present.items() if "Dg" in tps and "R" in tps]
print(f"[1] GSE227903 Dg+Relapse paired patients: {len(pairs)} -> {sorted(pairs)}")

rows = []
for pt in sorted(pairs):
    dg, rr = f"{pt}_Dg", f"{pt}_R"
    hdg,adg = hdsats(dg); hrr,arr = hdsats(rr)
    rows.append(dict(patient=pt,
        RLS=fgw2(dg, rr),
        HDS_Dg=hdg, HDS_R=hrr, dHDS=hrr-hdg,
        ATS_Dg=adg, ATS_R=arr, dATS=arr-adg,
        SCS_Dg=scs(dg), SCS_R=scs(rr), dSCS=scs(rr)-scs(dg),
        sparse_involved=(dg in sparse_set or rr in sparse_set)))
res = pd.DataFrame(rows)
res.to_csv(os.path.join(D_FGW, "paired_rls_scs.csv"), index=False)
print("[2] per-patient:")
print(res.round(4).to_string(index=False))

# -- paired tests (signed-rank if scipy, else sign test), WITH and WITHOUT the sparse patient --
def test_delta(x, alt):
    x = np.asarray([v for v in x if np.isfinite(v)])
    n_up = int((x>0).sum()); n = len(x); med = float(np.median(x))
    p_w = None
    try:
        from scipy.stats import wilcoxon
        if n >= 1 and np.any(x != 0):
            p_w = float(wilcoxon(x, alternative=alt, zero_method="wilcox").pvalue)
    except Exception:
        p_w = None
    return dict(n=n, n_up=n_up, median=round(med,4), wilcoxon_p=(round(p_w,4) if p_w is not None else "NA(no scipy)"))

def report(label, df):
    print(f"\n[3] paired tests ({label}, n={len(df)}):")
    print(f"    dHDS (H3 relapse moves AWAY from healthy; expect >0): {test_delta(df.dHDS, 'greater')}")
    print(f"    dATS (relapse moves TOWARD AML barycenter; expect <0): {test_delta(df.dATS, 'less')}")
    print(f"    dSCS (LSC->immune malignant-wt comm up; expect >0)    : {test_delta(df.dSCS, 'greater')}")
    print(f"    RLS  (Dg->Relapse topological move, descriptive): median={df.RLS.median():.4f} range=[{df.RLS.min():.4f},{df.RLS.max():.4f}]")

report("all 8 patients", res)
report("excluding sparse-relapse patient 4978", res[~res.sparse_involved])
print(f"\n[done] wrote {os.path.join(D_FGW,'paired_rls_scs.csv')}")
