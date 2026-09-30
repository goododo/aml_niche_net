#!/usr/bin/env python
# 05_opportunity_distance.py ----
# INPUT  : <root>/05_ccc/tensors/<dataset>/<sample>__ccc_cellchat.csv   per-sample CellChat LR rows,
#              UNFILTERED (pval up to 0.99 -- this script applies the filter)
#          <root>/05_ccc/ccc_node_presence.csv     per (sample,bin) cell counts
#          <root>/05_ccc/ccc_sample_manifest.csv   dataset, uid_patient, Timepoint
#          <root>/07_fgw/fgw_vocab.json            aml_timepoints
#          <root>/06_distance/scaccordion_distance__prop7.csv   canonical sample order + the
#              KNOWN ANSWER this script's GATE 1 replica must reproduce exactly
# OUTPUT : <root>/06_distance/scaccordion_distance__{ab_only,cellchat_w,ccc_x_ab}.csv
#          <root>/08_scoring/opportunity_gate1.csv   one row per arm x contrast x metric
#          <root>/08_scoring/opportunity_null.csv    the shuffled-molecular null, one row per draw
# WHAT IT DOES : adds the rung the pre-registered ladder never had -- edge = molecular strength x
#          co-abundance opportunity -- and asks the decisive question the ladder could not:
#          does the MOLECULAR term carry patient-specific information beyond co-abundance?
#          Arms, each a per-sample 7x7 = 49 vector scored by the registered GATE 1 statistic:
#            prop7_mine     7 compartment proportions, this script's own metric (internal reference)
#            ab_only        outer(p, p), the opportunity term alone
#            cellchat_w     sum of significant LR prob per (sender, receiver) -- molecular alone
#            ccc_x_ab       cellchat_w * ab_only, elementwise -- the missing rung
#            shuffled_x_ab  molecular term REASSIGNED BETWEEN SAMPLES WITHIN DATASET x real
#                           opportunity, N_PERM draws -- THE CONTROL. Shuffling is within dataset
#                           so batch structure is preserved and only patient identity is destroyed;
#                           a global shuffle would inject batch mismatch and make the control fail
#                           for the wrong reason.
#          DECISION RULE, fixed before running: if ccc_x_ab's median percentile is not better than
#          the 5th percentile of the shuffled null, the molecular term carries no patient-specific
#          information beyond co-abundance on these data.
#          Two self-checks make the machinery auditable: (1) the GATE 1 replica must reproduce the
#          published prop7 numbers to the digit; (2) under cosine, ab_only's distance is
#          1 - cos(p,q)^2 while prop7's is 1 - cos(p,q) -- monotone equivalent, so the two arms
#          must return IDENTICAL percentiles. Either check failing stops the script.
# Usage : general_env/bin/python scripts/06_distance/05_opportunity_distance.py [--root ...] [--nperm 200]
# ---------------------------------------------------------------------------------------------
# RESULTS (run 2026-09-30, --nperm 200). Both self-checks passed: the GATE 1 replica reproduced
# the published prop7 numbers to the digit, and ab_only returned percentiles identical to
# prop7_mine under cosine (monotone-equivalent by construction).
#
# Dx_to_Relapse median percentile, lower = better:
#     cellchat_w (molecular alone)        0.1304      <- best arm
#     prop7_published (existing ladder)   0.1739
#     ccc_x_ab (molecular x opportunity)  0.2500
#     prop7_TV (composition, TV metric)   0.2609
#     ab_only (opportunity alone)         0.3043
#
# THE REGISTERED RULE FIRES. ccc_x_ab real 0.2500 vs shuffled-molecular null median 0.3333,
# null 5th percentile 0.1364; the real arm sits at the 18th percentile of its own null. The
# molecular term is not distinguishable from a patient-scrambled one once abundance multiplies it.
#
# AND THE DIRECTION IS OPPOSITE TO THE HYPOTHESIS. Multiplying by co-abundance opportunity moves
# the arm from 0.1304 to 0.2500 -- WORSE -- and lands it exactly on composition (ccc_x_ab vs
# prop7_TV: median paired difference 0.0000, wilcoxon p 1.000). Mechanism: the abundance term
# dominates the product and erases the molecular contribution, converting a communication arm
# into a composition arm. ccc_x_ab vs cellchat_w is +0.1250 (worse) but not demonstrable at
# n=11 (CI -0.136..+0.391).
#
# LIMIT, unchanged: every head-to-head remains "not demonstrable" at 11 paired patients; the
# bootstrap CIs span 0.1-0.4 percentile points. One nominal exception, cellchat_w vs prop7_TV on
# Dx_to_Treatment (-0.1818, p 0.024), is an artifact of the handicapped TV metric -- against the
# published composition arm the same comparison is 0.0000, p 0.426 -- and is not claimed.
# ---------------------------------------------------------------------------------------------
import argparse, json, os, glob
import numpy as np
import pandas as pd
from scipy import stats

DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
CCC_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
PVAL_MAX = 0.05
# copied verbatim from 08_scoring/17_scaccordion_gate1.py so the statistic is the registered one
REL_TP = {"Relapse", "Relapse2"}
TRT_TP = {"On_treatment", "Post_induction", "Post_consolidation",
          "Post_treatment_unspecified", "Refractory"}
MIN_POOL = 3
SEED = 491638
KNOWN_PROP7 = dict(median_pct=0.173913043478, p_raw=0.0048828125)   # scaccordion_gate1.csv

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--nperm", type=int, default=200)
args = ap.parse_args()
D_CCC = os.path.join(args.root, "05_ccc")
D_DST = os.path.join(args.root, "06_distance")
D_SCO = os.path.join(args.root, "08_scoring")
rng = np.random.default_rng(SEED)

# ---- 1. sample set, pairs and pools, built exactly as 17 builds them ---------------------------
prop7 = pd.read_csv(os.path.join(D_DST, "scaccordion_distance__prop7.csv"), index_col=0)
SAMPLES = list(prop7.index)
AML_TP = set(json.load(open(os.path.join(args.root, "07_fgw", "fgw_vocab.json")))["aml_timepoints"])
man = pd.read_csv(os.path.join(D_CCC, "ccc_sample_manifest.csv")).set_index("sample")
idx = man.reindex(SAMPLES)[["dataset", "uid_patient", "Timepoint"]]
idx.columns = ["dataset", "patient", "timepoint"]
print("[1] %d samples, %d datasets, %d patients" % (len(SAMPLES), idx.dataset.nunique(), idx.patient.nunique()))

pairs = []
for (ds, pt), g in idx.groupby(["dataset", "patient"]):
    tps = {t: s for s, t in zip(g.index, g.timepoint)}
    if "Diagnosis" not in tps:
        continue
    for kind, want in [("Dx_to_Relapse", REL_TP), ("Dx_to_Treatment", TRT_TP)]:
        hit = [t for t in g.timepoint if t in want]
        if hit:
            pairs.append(dict(kind=kind, dataset=ds, patient=pt, a=tps["Diagnosis"], b=tps[sorted(hit)[0]]))
P = pd.DataFrame(pairs)
print("[1] pairs: %s" % P.kind.value_counts().to_dict())

POS = {s: i for i, s in enumerate(SAMPLES)}
POOLS = {}
for _, r in P.iterrows():
    POOLS[(r.kind, r.patient)] = np.array(
        [POS[s] for s in SAMPLES
         if idx.at[s, "dataset"] == r.dataset and idx.at[s, "patient"] != r.patient
         and idx.at[s, "timepoint"] in AML_TP], dtype=int)


def gate1(Dnp, grp):
    """Percentile rank of the true partner among same-dataset other-patient AML graphs.
    Replica of 17_scaccordion_gate1.py; the round(...,12) is load-bearing for tie handling."""
    pcts = []
    for _, r in grp.iterrows():
        pool = POOLS[(r.kind, r.patient)]
        if len(pool) < MIN_POOL:
            continue
        ia, ib = POS[r.a], POS[r.b]
        d_true = Dnp[ia, ib]
        d_pool = Dnp[ia, pool[pool != ib]]
        if len(d_pool) < MIN_POOL:
            continue
        pcts.append(round(float((d_pool < d_true).sum()) / len(d_pool), 12))
    return np.array(pcts)


def score(Dnp, kind):
    pcts = gate1(Dnp, P[P.kind == kind])
    if len(pcts) < 5:
        return None
    w = stats.wilcoxon(pcts - 0.5, alternative="less")
    return dict(n=len(pcts), median_pct=float(np.median(pcts)), p_raw=float(w.pvalue),
                top1=int((pcts == 0.0).sum()))

## -- SELF-CHECK 1: the replica must reproduce the published prop7 numbers to the digit ----------
chk = score(prop7.reindex(index=SAMPLES, columns=SAMPLES).to_numpy(float), "Dx_to_Relapse")
print("\n[SELF-CHECK 1] published prop7 Dx_to_Relapse: median_pct %.12f (want %.12f) | p %.10f (want %.10f)"
      % (chk["median_pct"], KNOWN_PROP7["median_pct"], chk["p_raw"], KNOWN_PROP7["p_raw"]))
assert abs(chk["median_pct"] - KNOWN_PROP7["median_pct"]) < 1e-9, "GATE 1 replica does not reproduce prop7"
assert abs(chk["p_raw"] - KNOWN_PROP7["p_raw"]) < 1e-9, "GATE 1 replica p does not reproduce prop7"
print("               replica verified")

# ---- 2. the two per-sample 7x7 terms ------------------------------------------------------------
pres = pd.read_csv(os.path.join(D_CCC, "ccc_node_presence.csv"))
pres = pres[pres.hierarchy_bin.isin(CCC_NODES)]
cnt = pres.pivot_table(index="sample", columns="hierarchy_bin", values="n_cells", aggfunc="sum")
cnt = cnt.reindex(index=SAMPLES, columns=CCC_NODES).fillna(0.0)
prop = cnt.div(cnt.sum(1).replace(0, np.nan), axis=0).fillna(0.0).to_numpy(float)
A = np.einsum("si,sj->sij", prop, prop)                      # opportunity: p_i * p_j
print("\n[2] opportunity term A: %s | samples with zero cells: %d" % (A.shape, int((cnt.sum(1) == 0).sum())))

W = np.zeros((len(SAMPLES), 7, 7))
n_missing, kept, total = 0, 0, 0
bi = {b: i for i, b in enumerate(CCC_NODES)}
for si, s in enumerate(SAMPLES):
    hits = glob.glob(os.path.join(D_CCC, "tensors", "*", "%s__ccc_cellchat.csv" % s))
    if not hits:
        n_missing += 1
        continue
    t = pd.read_csv(hits[0])
    total += len(t)
    t = t[(t.pval <= PVAL_MAX) & t.sender_bin.isin(CCC_NODES) & t.receiver_bin.isin(CCC_NODES)]
    kept += len(t)
    for (sb, rb), g in t.groupby(["sender_bin", "receiver_bin"]):
        W[si, bi[sb], bi[rb]] = g.prob.sum()
print("[2] molecular term W: tensors found for %d of %d samples (%d missing -> all-zero)"
      % (len(SAMPLES) - n_missing, len(SAMPLES), n_missing))
print("    LR rows kept at pval <= %.2f: %d of %d (%.1f%%) | nonzero edge slots: %.1f%%"
      % (PVAL_MAX, kept, total, 100 * kept / max(total, 1), 100 * (W != 0).mean()))

# ---- 3. arms and metrics -------------------------------------------------------------------------
def cos_d(X):
    n = np.sqrt((X ** 2).sum(1, keepdims=True)); n[n == 0] = 1.0
    D = 1.0 - (X / n) @ (X / n).T
    D = (D + D.T) / 2; np.fill_diagonal(D, 0.0)
    return D


def eucl1_d(X):
    r = X.sum(1, keepdims=True); r[r == 0] = 1.0
    Y = X / r
    D = np.sqrt(np.maximum(((Y ** 2).sum(1)[:, None] + (Y ** 2).sum(1)[None, :] - 2 * Y @ Y.T), 0))
    D = (D + D.T) / 2; np.fill_diagonal(D, 0.0)
    return D


F = len(SAMPLES)
ARMS = {
    "prop7_mine": prop,
    "ab_only":    A.reshape(F, 49),
    "cellchat_w": W.reshape(F, 49),
    "ccc_x_ab":   (W * A).reshape(F, 49),
}
rows = []
DM = {}
for arm, X in ARMS.items():
    for metric, f in [("cos", cos_d), ("eucl1", eucl1_d)]:
        D = f(np.asarray(X, float))
        DM[(arm, metric)] = D
        for kind in ["Dx_to_Relapse", "Dx_to_Treatment"]:
            r = score(D, kind)
            if r:
                rows.append(dict(arm=arm, metric=metric, contrast=kind, kind_of_evidence="new_rung", **r))
    if arm in ("ab_only", "cellchat_w", "ccc_x_ab"):
        pd.DataFrame(DM[(arm, "cos")], index=SAMPLES, columns=SAMPLES).to_csv(
            os.path.join(D_DST, "scaccordion_distance__%s.csv" % arm))
G = pd.DataFrame(rows)
print("\n[3] arms scored\n%s" % G.pivot_table(index=["arm", "metric"], columns="contrast",
                                              values=["median_pct", "p_raw"]).round(4).to_string())

## -- SELF-CHECK 2: under cosine, ab_only is a monotone transform of prop7 (1-cos^2 vs 1-cos), ----
## so the two arms must return IDENTICAL percentiles. A mismatch means the construction is wrong.
for kind in ["Dx_to_Relapse", "Dx_to_Treatment"]:
    a = score(DM[("prop7_mine", "cos")], kind)
    b = score(DM[("ab_only", "cos")], kind)
    print("[SELF-CHECK 2] %s  prop7_mine %.12f vs ab_only %.12f" % (kind, a["median_pct"], b["median_pct"]))
    assert abs(a["median_pct"] - b["median_pct"]) < 1e-12, "ab_only != prop7 under cosine: construction bug"
print("               opportunity term verified as composition-equivalent")

# ---- 4. THE CONTROL: molecular term reassigned between samples, within dataset -------------------
ds_arr = idx.dataset.to_numpy()
groups = [np.where(ds_arr == d)[0] for d in pd.unique(ds_arr)]
null_rows = []
for draw in range(args.nperm):
    perm = np.arange(F)
    for g in groups:
        if len(g) > 1:
            perm[g] = g[rng.permutation(len(g))]
    Xp = (W[perm] * A).reshape(F, 49)
    for metric, f in [("cos", cos_d), ("eucl1", eucl1_d)]:
        D = f(Xp)
        for kind in ["Dx_to_Relapse", "Dx_to_Treatment"]:
            r = score(D, kind)
            if r:
                null_rows.append(dict(draw=draw, metric=metric, contrast=kind,
                                      median_pct=r["median_pct"], p_raw=r["p_raw"]))
N = pd.DataFrame(null_rows)
N.to_csv(os.path.join(D_SCO, "opportunity_null.csv"), index=False)

print("\n[4] SHUFFLED-MOLECULAR NULL (%d draws, permuted within dataset)" % args.nperm)
print("    %-14s %-8s %9s %9s %9s %9s %8s" % ("contrast", "metric", "real", "null_med", "null_p5", "null_min", "real_pct"))
for kind in ["Dx_to_Relapse", "Dx_to_Treatment"]:
    for metric in ["cos", "eucl1"]:
        real = G[(G.arm == "ccc_x_ab") & (G.metric == metric) & (G.contrast == kind)].median_pct.iloc[0]
        nl = N[(N.metric == metric) & (N.contrast == kind)].median_pct.to_numpy()
        pct_in_null = float((nl <= real).mean())
        G.loc[(G.arm == "ccc_x_ab") & (G.metric == metric) & (G.contrast == kind), "null_p5"] = np.quantile(nl, 0.05)
        G.loc[(G.arm == "ccc_x_ab") & (G.metric == metric) & (G.contrast == kind), "frac_null_at_least_as_good"] = pct_in_null
        print("    %-14s %-8s %9.4f %9.4f %9.4f %9.4f %8.3f" %
              (kind, metric, real, np.median(nl), np.quantile(nl, 0.05), nl.min(), pct_in_null))
        verdict = "MOLECULAR TERM ADDS" if real < np.quantile(nl, 0.05) else "molecular term adds nothing beyond abundance"
        print("      -> %s" % verdict)
G.to_csv(os.path.join(D_SCO, "opportunity_gate1.csv"), index=False)
print("\n[done] wrote opportunity_gate1.csv, opportunity_null.csv and 3 distance matrices")
