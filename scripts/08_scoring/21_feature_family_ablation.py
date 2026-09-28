#!/usr/bin/env python
# 21_feature_family_ablation.py ----
# INPUT  : 07_fgw/fgw_nodes_long.csv                     (138 samples x 7 bins, global-z features)
#          06_distance/scaccordion_distance__prop7.csv    (canonical sample set + order, baseline arm)
#          05_ccc/ccc_sample_manifest.csv                 (uid_patient, Timepoint, dataset)
#          07_fgw/fgw_vocab.json                          (locked AML timepoint vocabulary)
# OUTPUT : 08_scoring/feature_family_ablation.csv          (one row per arm x contrast)
#          08_scoring/feature_family_null.csv              (size-matched random null distributions)
# WHAT IT DOES : asks WHICH of the 47 per-population features carry the GATE 1 patient-matching
#          signal -- by family and one by one -- with a size-matched random null so that "more
#          features" cannot masquerade as "better features".
#
# WHY THIS EXISTS. 04_node_feature_distance.py established that the 47-feature distance beats
# cell-type composition on GATE 1 (q = 0.004 vs 0.026), falsifying its own registered prediction.
# It does NOT say which features do the work, and the families are biologically very different
# things: stemness, pathway activity, apoptotic priming, metabolism, pseudotime, meta-programs.
# The answer changes what the next experiment is, so it is worth measuring rather than guessing.
#
# THIS IS A SCREEN, NOT A CONFIRMATORY TEST, AND IT IS LABELLED SO IN EVERY OUTPUT ROW.
# It is run after seeing that the full 47 wins, on the same 11 patients, over 60 arms. No p-value
# here may be reported as a test of a pre-specified hypothesis. What it can legitimately do is
# rank the families and single features for follow-up, and rule out the ones that clearly do
# nothing. Anything it nominates has to be re-tested on data these 11 patients are not in.
#
# THE ONE CONFOUND THAT WOULD MAKE THE WHOLE SCREEN MEANINGLESS, and how it is handled.
# The families differ in size by 14x (pg_ 14 features, pt_ 1). Euclidean distance on a 7*k vector
# changes character with k, so "pg_ beats pt_" could mean nothing but "14 features beat 1".
# Every arm is therefore scored TWICE: its own Wilcoxon p, and its rank inside a null of
# N_RAND random subsets OF THE SAME SIZE drawn from the same 47. Only the second number is
# comparable across arms of different size, and it is the one the summary reads.
#
# Usage : python scripts/08_scoring/21_feature_family_ablation.py [--n_rand 500]
import argparse
import json
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
CCC_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
META = ["dataset", "sample", "hierarchy_bin", "timepoint", "healthy",
        "n_cells_raw", "present", "mass", "sparse_flag"]
FAMILIES = ["st", "pg", "cs", "mt", "pt", "mp"]
REL_TP = {"Relapse", "Relapse2"}
TRT_TP = {"On_treatment", "Post_induction", "Post_consolidation",
          "Post_treatment_unspecified", "Refractory"}
MIN_POOL, MIN_PATIENTS = 3, 5
SEED = 491638

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=ROOT)
ap.add_argument("--n_rand", type=int, default=500)
args = ap.parse_args()
D_DST = os.path.join(args.root, "06_distance")
D_SCO = os.path.join(args.root, "08_scoring")
rng = np.random.default_rng(SEED)

# ---- 1. features, exactly the 04_node_feature_distance.py panel selector -----------------------
nd = pd.read_csv(os.path.join(args.root, "07_fgw", "fgw_nodes_long.csv"))
prop7 = pd.read_csv(os.path.join(D_DST, "scaccordion_distance__prop7.csv"), index_col=0)
SAMPLES = list(prop7.index)
num = [c for c in nd.columns if c not in META]
PANEL = [c for c in num
         if c.split("_")[0] in FAMILIES
         and not (c.endswith("_normal") or c.endswith("_malignant"))]
assert len(PANEL) == 47, "panel selector drifted: %d" % len(PANEL)
FAM = {f: [c for c in PANEL if c.split("_")[0] == f] for f in FAMILIES}
print("[1] 47 panel features over %d families: %s"
      % (len(FAM), ", ".join("%s=%d" % (f, len(v)) for f, v in FAM.items())))

WIDE = nd.pivot_table(index="sample", columns="hierarchy_bin", values=PANEL, dropna=False)
WIDE = WIDE.reindex(index=SAMPLES, columns=pd.MultiIndex.from_product([PANEL, CCC_NODES]))
assert np.isfinite(WIDE.to_numpy(float)).all(), "non-finite cells in the feature matrix"
COL = {f: [i for i, (feat, _) in enumerate(WIDE.columns) if feat == f] for f in PANEL}
XALL = WIDE.to_numpy(float)

# ---- 2. the GATE 1 rule, copied from 17_scaccordion_gate1.py ----------------------------------
man = pd.read_csv(os.path.join(args.root, "05_ccc", "ccc_sample_manifest.csv")).set_index("sample")
idx = man.reindex(SAMPLES)[["dataset", "uid_patient", "Timepoint"]]
idx.columns = ["dataset", "patient", "timepoint"]
AML_TP = set(json.load(open(os.path.join(args.root, "07_fgw", "fgw_vocab.json")))["aml_timepoints"])
pairs = []
for (ds, pt), g in idx.groupby(["dataset", "patient"]):
    tps = {t: s for s, t in zip(g.index, g.timepoint)}
    if "Diagnosis" not in tps:
        continue
    for kind, want in [("Dx_to_Relapse", REL_TP), ("Dx_to_Treatment", TRT_TP)]:
        hit = [t for t in g.timepoint if t in want]
        if hit:
            pairs.append(dict(kind=kind, dataset=ds, patient=pt,
                              a=tps["Diagnosis"], b=tps[sorted(hit)[0]]))
P = pd.DataFrame(pairs)
POS = {s: i for i, s in enumerate(SAMPLES)}
POOLS = {(r.kind, r.patient): np.array(
            [POS[s] for s in SAMPLES
             if idx.at[s, "dataset"] == r.dataset and idx.at[s, "patient"] != r.patient
             and idx.at[s, "timepoint"] in AML_TP], dtype=int)
         for _, r in P.iterrows()}
print("[2] pairs: %s | %d patients" % (P.kind.value_counts().to_dict(), P.patient.nunique()))

def gate1(D, grp):
    """Percentile rank of the true partner among same-dataset other-patient AML graphs."""
    pcts = []
    for _, r in grp.iterrows():
        pool = POOLS[(r.kind, r.patient)]
        if len(pool) < MIN_POOL:
            continue
        ia, ib = POS[r.a], POS[r.b]
        d_pool = D[ia, pool[pool != ib]]
        if len(d_pool) < MIN_POOL:
            continue
        pcts.append(round(float((d_pool < D[ia, ib]).sum()) / len(d_pool), 12))
    return np.array(pcts)

def score(cols, grp):
    """Euclidean distance on the selected columns, then the GATE 1 statistic."""
    D = squareform(pdist(XALL[:, cols], "euclidean"))
    p = gate1(D, grp)
    if len(p) < MIN_PATIENTS:
        return np.nan, np.nan, 0
    return (float(np.median(p)),
            float(stats.wilcoxon(p - 0.5, alternative="less").pvalue), len(p))

def cols_of(feats):
    out = []
    for f in feats:
        out.extend(COL[f])
    return out

# ---- 3. the arms --------------------------------------------------------------------------------
ARMS = {"full47": PANEL}
for f in FAMILIES:
    ARMS["only_" + f] = FAM[f]
    ARMS["drop_" + f] = [c for c in PANEL if c not in FAM[f]]
for c in PANEL:
    ARMS["single:" + c] = [c]
print("[3] %d arms: 1 reference + 6 leave-one-family-out + 6 single-family + %d single-feature"
      % (len(ARMS), len(PANEL)))

# ---- 4. size-matched null, the mandatory control -----------------------------------------------
SIZES = sorted({len(v) for v in ARMS.values()})
print("[4] size-matched nulls for k in %s, %d draws each" % (SIZES, args.n_rand))
NULL = {}
nullrows = []
for kind, grp in P.groupby("kind"):
    for k in SIZES:
        if k == len(PANEL):                       # nothing to draw against the full set
            continue
        meds = np.empty(args.n_rand)
        for t in range(args.n_rand):
            pick = list(rng.choice(PANEL, k, replace=False))
            meds[t] = score(cols_of(pick), grp)[0]
        NULL[(kind, k)] = meds
        nullrows.append(pd.DataFrame(dict(contrast=kind, k=k, median_pct=meds)))
    print("    %s done" % kind)
pd.concat(nullrows).to_csv(os.path.join(D_SCO, "feature_family_null.csv"), index=False)

# ---- 5. score every arm ------------------------------------------------------------------------
d_prop7 = prop7.reindex(index=SAMPLES, columns=SAMPLES).to_numpy(float)
rows = []
for kind, grp in P.groupby("kind"):
    base_med, base_p, base_n = (float(np.median(gate1(d_prop7, grp))),
                                float(stats.wilcoxon(gate1(d_prop7, grp) - 0.5,
                                                     alternative="less").pvalue),
                                len(gate1(d_prop7, grp)))
    rows.append(dict(contrast=kind, arm="prop7_baseline", block="baseline", k=7,
                     median_pct=base_med, wilcoxon_p=base_p, n_patients=base_n,
                     null_pct_of_arm=np.nan, beats_prop7=np.nan, kind_of_evidence="baseline"))
    for arm, feats in ARMS.items():
        med, p, n = score(cols_of(feats), grp)
        k = len(feats)
        nl = NULL.get((kind, k))
        # where does this arm sit inside its own size-matched null? low = better than chance for k
        npc = float((nl <= med).mean()) if nl is not None and med == med else np.nan
        block = ("reference" if arm == "full47" else
                 "leave_one_family_out" if arm.startswith("drop_") else
                 "single_family" if arm.startswith("only_") else "single_feature")
        rows.append(dict(contrast=kind, arm=arm, block=block, k=k, median_pct=med,
                         wilcoxon_p=p, n_patients=n, null_pct_of_arm=npc,
                         beats_prop7=bool(med < base_med) if med == med else np.nan,
                         kind_of_evidence="exploratory_screen"))
R = pd.DataFrame(rows)
# BH within contrast, over the scored arms only -- reported, and dependent by construction
for kind in R.contrast.unique():
    m = (R.contrast == kind) & R.wilcoxon_p.notna() & (R.block != "baseline")
    pv = R.loc[m, "wilcoxon_p"].to_numpy()
    o = np.argsort(pv); adj = np.empty(len(pv))
    adj[o] = np.minimum.accumulate((pv[o] * len(pv) / (np.arange(len(pv)) + 1))[::-1])[::-1]
    R.loc[m, "p_bh_within_contrast"] = np.minimum(adj, 1.0)
R.to_csv(os.path.join(D_SCO, "feature_family_ablation.csv"), index=False)

# ---- 6. print the screen ------------------------------------------------------------------------
for kind in ["Dx_to_Relapse", "Dx_to_Treatment"]:
    g = R[R.contrast == kind]
    b = g[g.block == "baseline"].iloc[0]
    print("\n" + "=" * 96)
    print("%s   (n = %d patients; composition baseline: median %.3f, p = %.4f)"
          % (kind, int(b.n_patients), b.median_pct, b.wilcoxon_p))
    print("=" * 96)
    for block, label in [("reference", "REFERENCE"),
                         ("single_family", "SINGLE FAMILY -- is it sufficient?"),
                         ("leave_one_family_out", "LEAVE ONE OUT -- is it necessary?")]:
        print("\n  %s" % label)
        print("  %-14s %3s %10s %9s %9s %14s" % ("arm", "k", "median", "wilcox p", "BH", "vs size-null"))
        for _, r in g[g.block == block].sort_values("median_pct").iterrows():
            print("  %-14s %3d %10.3f %9.4f %9.4f %13s"
                  % (r.arm, r.k, r.median_pct, r.wilcoxon_p, r.p_bh_within_contrast,
                     "%.3f" % r.null_pct_of_arm if r.null_pct_of_arm == r.null_pct_of_arm else "n/a"))
    sf = g[g.block == "single_feature"].sort_values("median_pct")
    print("\n  SINGLE FEATURE -- top 10 of 47 (screen only; %d arms tested on %d patients)"
          % (len(sf), int(b.n_patients)))
    print("  %-34s %10s %9s %9s %14s" % ("feature", "median", "wilcox p", "BH", "vs size-null"))
    for _, r in sf.head(10).iterrows():
        print("  %-34s %10.3f %9.4f %9.4f %13s"
              % (r.arm.split(":", 1)[1], r.median_pct, r.wilcoxon_p, r.p_bh_within_contrast,
                 "%.3f" % r.null_pct_of_arm if r.null_pct_of_arm == r.null_pct_of_arm else "n/a"))

print("\n[done] wrote feature_family_ablation.csv (%d rows) and feature_family_null.csv" % len(R))
print("       REMINDER: every non-baseline row is kind_of_evidence = exploratory_screen.")
