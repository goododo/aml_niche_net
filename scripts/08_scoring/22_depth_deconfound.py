#!/usr/bin/env python
# 22_depth_deconfound.py ----
# INPUT  : 07_fgw/fgw_nodes_long.csv                     (138 x 7, global-z per-population features)
#          06_distance/scaccordion_distance__{prop7,tv,tabular}.csv   (baseline + edge arms)
#          05_ccc/ccc_sample_manifest.csv , 07_fgw/fgw_vocab.json     (GATE 1 roster + vocabulary)
#          00_project/metadata/ALL__samples.tsv                        (med_nFeature = library depth)
# OUTPUT : 08_scoring/depth_deconfound.csv       (one row per arm x contrast x deconfounding mode)
# WHAT IT DOES : asks whether the per-population feature distance still beats cell-type composition
#          on GATE 1 once library sequencing depth has been removed from it.
#
# WHY. 04_node_feature_distance.py found the 47-feature distance beating composition on GATE 1
# (median percentile 0.043 vs 0.174; q 0.004 vs 0.026). A follow-up probe then found that this
# distance is the MOST depth-dependent arm in the whole ladder: Spearman(distance, |delta median
# nFeature|) = +0.207 within dataset, against +0.042 for composition and -0.09 to -0.12 for every
# edge-based arm. Gene-set scores computed from sparse counts depend on detection depth, so this
# has an explanation that is not biology. If the advantage does not survive removing depth, the
# node-feature line closes -- the same standard this project applied to the OT distances.
#
# THREE MODES, ALL FIXED HERE BEFORE THE RUN. They fail in different ways on purpose:
#
#   raw        no deconfounding. The published comparison, carried for reference.
#   residual   each of the 47 features is regressed on log10(med_nFeature) ACROSS SAMPLES, WITHIN
#              (dataset, hierarchy_bin), and replaced by its residual. Within-stratum because depth
#              differs by dataset and because a bin's score scale differs by bin; a pooled
#              regression would remove dataset and bin structure along with depth.
#   matched    no feature is altered. Instead the GATE 1 candidate pool for each patient is
#              restricted to distractors whose depth is CLOSER to the diagnosis sample than the
#              true partner's is. A true partner that still ranks first cannot be doing it on
#              depth, because every remaining distractor is a better depth match than it is.
#              This is the assumption-free version and it is the one that decides.
#
#   Two arms are carried through all three modes: node_feat (47) and prop7 (composition). Both must
#   be deconfounded the same way or the comparison is not like-for-like -- composition is depth
#   dependent too (+0.042), just far less.
#
# THE PLACEBO, AND HOW TO READ IT. Residualising on ANY covariate removes one degree of freedom
# from every feature, so some degradation is generic rather than informative. The control replaces
# depth with a random reshuffle of the same values and re-runs the residualisation. Three readings,
# fixed here so the outcome cannot be re-interpreted afterwards:
#   shuffled degrades as much as real  -> the degradation is generic; `residual` says nothing about
#                                         depth and only `matched_pool` is readable.
#   real degrades MORE than shuffled   -> depth carried matching signal beyond a random direction.
#                                         That is the confound being removed, and it is the healthy
#                                         outcome for this control, NOT a failure.
#   real degrades LESS than shuffled   -> depth is protective; would need explaining before use.
# In every case the verdict on the arm comes from whether it still beats composition AFTER removal,
# not from the size of the drop.
#
# Usage : python scripts/08_scoring/22_depth_deconfound.py [--n_placebo 200]
import argparse
import json
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
META_TSV = "/FAST/gr10634/gaozy/aml_niche_net/00_project/metadata/ALL__samples.tsv"
CCC_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
META = ["dataset", "sample", "hierarchy_bin", "timepoint", "healthy",
        "n_cells_raw", "present", "mass", "sparse_flag"]
FAMILIES = ["st", "pg", "cs", "mt", "pt", "mp"]
REL_TP = {"Relapse", "Relapse2"}
TRT_TP = {"On_treatment", "Post_induction", "Post_consolidation",
          "Post_treatment_unspecified", "Refractory"}
MIN_POOL, MIN_PATIENTS, SEED = 3, 5, 491638

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=ROOT)
ap.add_argument("--n_placebo", type=int, default=200)
args = ap.parse_args()
D_DST = os.path.join(args.root, "06_distance")
D_SCO = os.path.join(args.root, "08_scoring")
rng = np.random.default_rng(SEED)

# ---- roster, features, depth -------------------------------------------------------------------
nd = pd.read_csv(os.path.join(args.root, "07_fgw", "fgw_nodes_long.csv"))
prop7 = pd.read_csv(os.path.join(D_DST, "scaccordion_distance__prop7.csv"), index_col=0)
SAMPLES = list(prop7.index)
PANEL = [c for c in nd.columns if c not in META
         and c.split("_")[0] in FAMILIES
         and not (c.endswith("_normal") or c.endswith("_malignant"))]
assert len(PANEL) == 47, len(PANEL)

meta = pd.read_csv(META_TSV, sep="\t").set_index("sample")
depth = pd.to_numeric(meta.reindex(SAMPLES).med_nFeature, errors="coerce")
print("[0] %d samples | depth (med_nFeature) present for %d | median %.0f, range %.0f-%.0f"
      % (len(SAMPLES), int(depth.notna().sum()), depth.median(), depth.min(), depth.max()))
LOGD = np.log10(depth.to_numpy(float))
if not np.isfinite(LOGD).all():
    raise SystemExit("depth missing for %d samples; GATE 1 cannot be depth-matched" % int((~np.isfinite(LOGD)).sum()))

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

# ---- the two GATE 1 variants -------------------------------------------------------------------
def gate1(D, grp, depth_matched=False):
    pcts = []
    for _, r in grp.iterrows():
        pool = POOLS[(r.kind, r.patient)]
        if len(pool) < MIN_POOL:
            continue
        ia, ib = POS[r.a], POS[r.b]
        cand = pool[pool != ib]
        if depth_matched:
            # keep only distractors at least as close to the anchor in depth as the true partner
            gap_true = abs(LOGD[ia] - LOGD[ib])
            cand = cand[np.abs(LOGD[ia] - LOGD[cand]) <= gap_true]
        if len(cand) < MIN_POOL:
            continue
        pcts.append(round(float((D[ia, cand] < D[ia, ib]).sum()) / len(cand), 12))
    return np.array(pcts)

def summarise(D, grp, depth_matched=False):
    p = gate1(D, grp, depth_matched)
    if len(p) < MIN_PATIENTS:
        return dict(median_pct=np.nan, wilcoxon_p=np.nan, n_patients=len(p))
    return dict(median_pct=float(np.median(p)),
                wilcoxon_p=float(stats.wilcoxon(p - 0.5, alternative="less").pvalue),
                n_patients=len(p))

# ---- feature matrix builders -------------------------------------------------------------------
def wide(frame):
    w = frame.pivot_table(index="sample", columns="hierarchy_bin", values=PANEL, dropna=False)
    return w.reindex(index=SAMPLES, columns=pd.MultiIndex.from_product([PANEL, CCC_NODES]))

def residualise(frame, logd_by_sample):
    """Regress every feature on log10 depth within (dataset, bin); keep residuals. OLS, 1 covariate."""
    out = frame.copy()
    out["_ld"] = out["sample"].map(logd_by_sample)
    for (_, _), g in out.groupby(["dataset", "hierarchy_bin"], sort=False):
        x = g["_ld"].to_numpy(float)
        if len(g) < 4 or not np.isfinite(x).all() or x.std() == 0:
            continue                                   # too few samples to fit; leave as is
        xc = x - x.mean()
        denom = float((xc ** 2).sum())
        for c in PANEL:
            y = g[c].to_numpy(float)
            beta = float((xc * (y - y.mean())).sum()) / denom
            out.loc[g.index, c] = y - beta * xc
    return out.drop(columns="_ld")

def eucl(frame):
    X = wide(frame).to_numpy(float)
    assert np.isfinite(X).all()
    return squareform(pdist(X, "euclidean"))

D_RAW = eucl(nd)
LD_MAP = dict(zip(SAMPLES, LOGD))
D_RES = eucl(residualise(nd, LD_MAP))
D_P7 = prop7.reindex(index=SAMPLES, columns=SAMPLES).to_numpy(float)

# composition's own depth residualisation: prop7 is a CLR of 7 counts, so residualise the distance's
# inputs the same way -- regress the 7 CLR coordinates on log depth within dataset.
pres = pd.read_csv(os.path.join(args.root, "05_ccc", "ccc_node_presence.csv"))
comp = (pres.pivot_table(index="sample", columns="hierarchy_bin", values="n_cells", fill_value=0)
        .reindex(index=SAMPLES, columns=CCC_NODES).fillna(0.0).to_numpy(float)) + 0.5
clr = np.log(comp) - np.log(comp).mean(1, keepdims=True)
dsv = idx.dataset.to_numpy()
clr_r = clr.copy()
for d in np.unique(dsv):
    m = dsv == d
    if m.sum() < 4:
        continue
    xc = LOGD[m] - LOGD[m].mean()
    dn = float((xc ** 2).sum())
    if dn == 0:
        continue
    for j in range(clr.shape[1]):
        y = clr[m, j]
        clr_r[m, j] = y - (float((xc * (y - y.mean())).sum()) / dn) * xc
D_P7R = squareform(pdist(clr_r, "euclidean"))

# ---- depth dependence of each distance, for the record -----------------------------------------
iu = np.triu_indices(len(SAMPLES), 1)
same_ds = (dsv[:, None] == dsv[None, :])[iu]
dnf = np.abs(depth.to_numpy(float)[:, None] - depth.to_numpy(float)[None, :])[iu]
print("\n[1] depth dependence of each distance, within-dataset pairs (n=%d)" % int(same_ds.sum()))
DEP = {}
for nm, D in [("node_feat raw", D_RAW), ("node_feat residual", D_RES),
              ("prop7 raw", D_P7), ("prop7 residual", D_P7R)]:
    r = float(stats.spearmanr(D[iu][same_ds], dnf[same_ds]).statistic)
    DEP[nm] = r
    print("    %-22s Spearman(distance, |d med_nFeature|) = %+.3f" % (nm, r))

# ---- the table ---------------------------------------------------------------------------------
rows = []
ARMS = [("node_feat", D_RAW, D_RES), ("prop7", D_P7, D_P7R)]
for kind, grp in P.groupby("kind"):
    for arm, Draw, Dres in ARMS:
        for mode, D, dm in [("raw", Draw, False), ("residual", Dres, False),
                            ("matched_pool", Draw, True),
                            ("residual+matched_pool", Dres, True)]:
            s = summarise(D, grp, dm)
            rows.append(dict(contrast=kind, arm=arm, deconf_mode=mode,
                             depth_rho=DEP.get("%s %s" % (arm, "residual" if "residual" in mode else "raw")),
                             **s))
R = pd.DataFrame(rows)

# ---- the placebo: shuffled depth must NOT reproduce the residual result ------------------------
print("\n[2] placebo: residualising on SHUFFLED depth, %d draws" % args.n_placebo)
plac = {}
for kind, grp in P.groupby("kind"):
    meds = np.empty(args.n_placebo)
    for t in range(args.n_placebo):
        sh = dict(zip(SAMPLES, rng.permutation(LOGD)))
        meds[t] = summarise(eucl(residualise(nd, sh)), grp)["median_pct"]
    plac[kind] = meds
    obs = float(R[(R.contrast == kind) & (R.arm == "node_feat") & (R.deconf_mode == "residual")].median_pct.iloc[0])
    pc = float((meds <= obs).mean())
    raw = float(R[(R.contrast == kind) & (R.arm == "node_feat") & (R.deconf_mode == "raw")].median_pct.iloc[0])
    print("    %-16s real-depth residual median %.3f | shuffled-depth median %.3f (IQR %.3f-%.3f)"
          % (kind, obs, np.median(meds), np.percentile(meds, 25), np.percentile(meds, 75)))
    verdict = ("real degrades MORE than shuffled -> depth carried signal (healthy)" if pc > 0.5
               else "shuffled degrades as much or more -> degradation is generic")
    print("    %-16s raw %.3f -> real %.3f | shuffled <= real in %.1f%% of draws | %s"
          % ("", raw, obs, 100 * pc, verdict))
    R.loc[(R.contrast == kind) & (R.arm == "node_feat") & (R.deconf_mode == "residual"),
          "placebo_frac_shuffled_as_good"] = pc
R.to_csv(os.path.join(D_SCO, "depth_deconfound.csv"), index=False)

# ---- read it out --------------------------------------------------------------------------------
for kind in ["Dx_to_Relapse", "Dx_to_Treatment"]:
    g = R[R.contrast == kind]
    print("\n" + "=" * 92)
    print("%s" % kind)
    print("=" * 92)
    print("  %-22s %-12s %5s %10s %10s" % ("mode", "arm", "n", "median", "wilcox p"))
    for mode in ["raw", "residual", "matched_pool", "residual+matched_pool"]:
        for arm in ["node_feat", "prop7"]:
            r = g[(g.arm == arm) & (g.deconf_mode == mode)]
            if not len(r):
                continue
            r = r.iloc[0]
            print("  %-22s %-12s %5d %10.3f %10.4f"
                  % (mode, arm, int(r.n_patients), r.median_pct, r.wilcoxon_p))
        nf = g[(g.arm == "node_feat") & (g.deconf_mode == mode)].iloc[0]
        p7 = g[(g.arm == "prop7") & (g.deconf_mode == mode)].iloc[0]
        if nf.median_pct == nf.median_pct and p7.median_pct == p7.median_pct:
            print("  %-22s %s" % ("", "-> node_feat %s composition (%.3f vs %.3f)"
                  % ("BEATS" if nf.median_pct < p7.median_pct else
                     "TIES" if nf.median_pct == p7.median_pct else "LOSES TO",
                     nf.median_pct, p7.median_pct)))
print("\n[done] wrote depth_deconfound.csv (%d rows)" % len(R))
