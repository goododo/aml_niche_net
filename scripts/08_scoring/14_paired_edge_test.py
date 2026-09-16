#!/usr/bin/env python
# 14_paired_edge_test.py ----
# INPUT  : <root>/06_distance/edge_distance.csv , <root>/07_fgw/fgw_input_index.csv ,
#          <root>/01_preprocess/00_curated_manifest.csv , <root>/01_preprocess/03_qc_report__ALL.csv ,
#          <root>/01_preprocess/02_study_split.csv , <root>/08_scoring/paired_roster.csv
# OUTPUT : <root>/08_scoring/paired_edge_test.csv   (one row per directed edge per contrast)
# WHAT IT DOES : asks, edge by edge, whether any of the 49 directed edges moves in a CONSISTENT
#          direction between a patient's two timepoints. This is the paired counterpart of
#          04_edge_regression.py, and it is not an omnibus.
#
# WHY THIS AND NOT ANOTHER METRIC. D4c measured that the shared timepoint effect is only 3-5% of
# the variance of C while patient identity is 53-74%. D5 then measured that the within-patient
# change is real, 3.88x the split-half noise floor for Dx->Treatment. Real change with no shared
# direction is the average over 49 edges; individual edges can still move consistently and be
# washed out by that average. This script looks for them. It needs no optimal transport, so it is
# unaffected by whatever D1 finds about the FGW statistic.
#
# THIS IS EXPLORATORY AND CANNOT BE CONFIRMATORY. 02_study_split.csv puts every longitudinal
# dataset in the Validation arm, so there is no held-out paired cohort left to confirm anything
# found here. Nothing in this file is pre-registered. Findings are hypotheses for an external
# cohort, and the leave-one-dataset-out block below is a consistency check, NOT a validation.
#
# Usage : python scripts/08_scoring/14_paired_edge_test.py [--n_perm 10000]
import os, sys, argparse
import numpy as np
import pandas as pd
from scipy import stats

FGW_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
EDGE_NAMES = ["%s->%s" % (a, b) for a in FGW_NODES for b in FGW_NODES]
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
REL_TP = {"Relapse", "Relapse2"}
TRT_TP = {"On_treatment", "Post_induction", "Post_consolidation",
          "Post_treatment_unspecified", "Refractory"}
SEED = 491638

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--n_perm", type=int, default=10000)
args = ap.parse_args()

D_DST = os.path.join(args.root, "06_distance")
D_FGW = os.path.join(args.root, "07_fgw")
D_PRE = os.path.join(args.root, "01_preprocess")
D_OUT = os.path.join(args.root, "08_scoring")

edges = pd.read_csv(os.path.join(D_DST, "edge_distance.csv"))
index = pd.read_csv(os.path.join(D_FGW, "fgw_input_index.csv"))
man = (pd.read_csv(os.path.join(D_PRE, "00_curated_manifest.csv"))
       [["dataset", "sample", "patient_id"]].drop_duplicates(["dataset", "sample"]))
qc = (pd.read_csv(os.path.join(D_PRE, "03_qc_report__ALL.csv")).rename(columns={"Sample": "sample"})
      [["dataset", "sample", "med_ncount_final"]].drop_duplicates(["dataset", "sample"]))
spl = pd.read_csv(os.path.join(D_PRE, "02_study_split.csv"))[["dataset", "split"]]
idx = index.merge(man, on=["dataset", "sample"]).merge(qc, on=["dataset", "sample"])
idx["logdepth"] = np.log10(idx.med_ncount_final.clip(lower=1))

## -- pairs, same rule as 11/12/13 ----
pairs = []
for (ds, pt), g in idx.groupby(["dataset", "patient_id"]):
    tps = dict(zip(g.timepoint, zip(g.dataset, g["sample"])))
    if "Diagnosis" not in tps:
        continue
    for kind, want in [("Dx_to_Relapse", REL_TP), ("Dx_to_Treatment", TRT_TP)]:
        hit = sorted([t for t in g.timepoint if t in want])
        if hit:
            pairs.append(dict(kind=kind, patient=pt, dataset=ds, a=tps["Diagnosis"], b=tps[hit[0]]))
P = pd.DataFrame(pairs)

## -- ROSTER, printed before any result ----
KEYS = sorted({k for _, r in P.iterrows() for k in (r.a, r.b)})
arm = spl.set_index("dataset")["split"].to_dict()
print("=" * 78)
print("SAMPLE SET: PAIRED SAMPLES ONLY")
print("=" * 78)
print("%d pairs, %d patients, %d samples, %d datasets."
      % (len(P), P.patient.nunique(), len(KEYS), P.dataset.nunique()))
for ds, g in P.groupby("dataset"):
    print("  %-12s %-11s %2d pairs  (%s)"
          % (ds, arm.get(ds, "?"), len(g), ", ".join("%s %d" % kv for kv in g.kind.value_counts().items())))
n_val = sum(len(g) for ds, g in P.groupby("dataset") if arm.get(ds) == "Validation")
print("\n  %d of %d pairs sit in the Validation arm. There is no held-out paired cohort left,"
      % (n_val, len(P)))
print("  so everything below is EXPLORATORY and cannot be called validated.\n")

## -- C, per sample, straight from production ----
CST = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                         values="C", fill_value=1.0)
       .reindex(columns=pd.MultiIndex.from_product([FGW_NODES, FGW_NODES])).fillna(1.0))
occ = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                         values="weight_probsum", fill_value=0.0)
       .reindex(columns=pd.MultiIndex.from_product([FGW_NODES, FGW_NODES])).fillna(0.0))
dep = idx.set_index(["dataset", "sample"])["logdepth"].to_dict()
miss = [k for k in KEYS if k not in CST.index]
if miss:
    raise SystemExit("paired samples missing from edge_distance.csv: %s" % miss)


def contrast(sub):
    """dC and d(logdepth) per pair, aligned rows."""
    dC = np.vstack([CST.loc[r.b].to_numpy(float) - CST.loc[r.a].to_numpy(float)
                    for _, r in sub.iterrows()])
    dD = np.array([dep[r.b] - dep[r.a] for _, r in sub.iterrows()])
    return dC, dD


def edge_stats(y, x):
    """Paired test on one edge. Two statistics, both reported, neither privileged:
       raw   = Wilcoxon signed-rank on dC (no covariate)
       adj   = intercept of dC ~ d(logdepth), the same form GATE 2 used, because within-pair
               sequencing depth ratios in this cohort run 0.23x to 3.54x."""
    if np.allclose(y, 0):
        return np.nan, 1.0, np.nan, 1.0
    w_p = float(stats.wilcoxon(y).pvalue)
    X = np.column_stack([np.ones(len(y)), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    dof = len(y) - 2
    if dof <= 0:
        return float(np.median(y)), w_p, float(beta[0]), 1.0
    se = np.sqrt((res @ res / dof) * np.linalg.inv(X.T @ X)[0, 0])
    a_p = float(2 * stats.t.sf(abs(beta[0] / se), dof)) if se > 0 else 1.0
    return float(np.median(y)), w_p, float(beta[0]), a_p


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); run = 1.0
    for i in range(n - 1, -1, -1):
        run = min(run, p[o[i]] * n / (i + 1)); q[o[i]] = run
    return q


rows = []
for kind in sorted(P.kind.unique()):
    sub = P[P.kind == kind].reset_index(drop=True)
    dC, dD = contrast(sub)
    occ_pair = np.vstack([((occ.loc[r.a].to_numpy(float) > 0) | (occ.loc[r.b].to_numpy(float) > 0))
                          for _, r in sub.iterrows()]).mean(0)
    st = [edge_stats(dC[:, j], dD) for j in range(49)]
    q_raw = bh([s[1] for s in st]); q_adj = bh([s[3] for s in st])
    for j, e in enumerate(EDGE_NAMES):
        y = dC[:, j]
        same = max((y > 0).mean(), (y < 0).mean()) if np.any(y != 0) else np.nan
        rows.append(dict(kind=kind, edge=e, n_pairs=len(sub), occupancy=occ_pair[j],
                         median_dC=st[j][0], frac_same_direction=same,
                         wilcoxon_p=st[j][1], wilcoxon_q=q_raw[j],
                         depth_adj_intercept=st[j][2], depth_adj_p=st[j][3], depth_adj_q=q_adj[j]))
R = pd.DataFrame(rows)
os.makedirs(D_OUT, exist_ok=True)
R.to_csv(os.path.join(D_OUT, "paired_edge_test.csv"), index=False)

## -- SELF-CHECK: sign-flip permutation. Flipping which timepoint is "second" destroys any real
## -- direction, so the significant-edge count must fall to the nominal rate. If it does not,
## -- the test is manufacturing direction and nothing below can be reported.
# Flipping which timepoint is "second" destroys any real direction, so the number of edges
# reaching p<.05 must fall to the nominal rate. The point of doing this at the FAMILY level is
# that the 49 edges are correlated (they share a within-sample ranking), so the expected count
# under the null is not simply 0.05*49 and has to be measured.
#
# Under a sign flip the ranks of |dC| do not change, only which ranks are counted as positive.
# So the ranks are computed ONCE per edge and every permutation is a matrix product. The same
# hand-rolled statistic scores the observed data and the null, so the comparison is internally
# consistent whatever scipy would do with ties.
def signed_rank_p_matrix(Y, S):
    """Y (n_pairs, 49) of dC; S (n_perm, n_pairs) of +-1. -> (n_perm, 49) two-sided p."""
    n, m = Y.shape
    A = np.abs(Y)
    Rk = np.vstack([stats.rankdata(A[:, j]) for j in range(m)]).T      # (n, 49), fixed under flips
    Rk = np.where(Y == 0, 0.0, Rk)                                     # zeros carry no rank
    pos = (S[:, :, None] * Y[None, :, :]) > 0                          # (n_perm, n, 49)
    Wp = np.einsum("pnj,nj->pj", pos.astype(float), Rk)                # (n_perm, 49)
    tot = Rk.sum(0)                                                    # per edge
    mu = tot / 2.0
    # exact-enough two-sided tail by the normal approximation to the signed-rank null; used
    # identically for observed and permuted values, so it cancels in the comparison
    nz = (Rk > 0).sum(0)
    sd = np.sqrt(np.array([ (Rk[:, j][Rk[:, j] > 0] ** 2).sum() / 4.0 for j in range(m) ]))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(sd > 0, (Wp - mu) / sd, 0.0)
    p = 2 * stats.norm.sf(np.abs(z))
    p[:, nz < 3] = 1.0
    return np.clip(p, 0.0, 1.0)

rng = np.random.default_rng(SEED)
print("=" * 78)
print("SELF-CHECK  sign-flip permutation (%d draws)" % args.n_perm)
print("=" * 78)
for kind in sorted(P.kind.unique()):
    sub = P[P.kind == kind].reset_index(drop=True)
    dC, dD = contrast(sub)
    n = len(sub)
    obs = int((signed_rank_p_matrix(dC, np.ones((1, n))) < 0.05).sum())
    S = rng.choice([-1.0, 1.0], size=(args.n_perm, n))
    null = (signed_rank_p_matrix(dC, S) < 0.05).sum(1)
    print("  %-17s observed %2d edges at p<.05 | null mean %.2f, 95th pct %.0f | p = %.4f"
          % (kind, obs, null.mean(), np.percentile(null, 95),
             (1.0 + (null >= obs).sum()) / (args.n_perm + 1.0)))

## -- what actually survives ----
print("\n" + "=" * 78)
print("RESULT  edges moving consistently between timepoints")
print("=" * 78)
for kind in sorted(P.kind.unique()):
    g = R[R.kind == kind]
    print("\n%s  (n=%d pairs)" % (kind, g.n_pairs.iloc[0]))
    print("  edges at raw p<.05 : wilcoxon %d/49, depth-adjusted %d/49"
          % (int((g.wilcoxon_p < .05).sum()), int((g.depth_adj_p < .05).sum())))
    print("  edges at BH q<.05  : wilcoxon %d/49, depth-adjusted %d/49"
          % (int((g.wilcoxon_q < .05).sum()), int((g.depth_adj_q < .05).sum())))
    top = g.sort_values("wilcoxon_p").head(5)
    print("  %-28s %6s %8s %9s %8s %8s" % ("edge", "occ", "med dC", "same dir", "wilc p", "wilc q"))
    for _, r in top.iterrows():
        print("  %-28s %5.0f%% %+8.3f %8.0f%% %8.4f %8.3f"
              % (r.edge, 100 * r.occupancy, r.median_dC, 100 * r.frac_same_direction,
                 r.wilcoxon_p, r.wilcoxon_q))

## -- consistency across datasets. NOT a validation: every dataset here is in the Validation arm. ----
print("\n" + "=" * 78)
print("CONSISTENCY  GSE227903 (the 16-pair majority) vs the rest")
print("=" * 78)
for kind in sorted(P.kind.unique()):
    sub = P[P.kind == kind].reset_index(drop=True)
    big = sub[sub.dataset == "GSE227903"].reset_index(drop=True)
    rest = sub[sub.dataset != "GSE227903"].reset_index(drop=True)
    if len(big) < 3 or len(rest) < 3:
        print("  %-17s not splittable (%d vs %d pairs)" % (kind, len(big), len(rest)))
        continue
    dCb, _ = contrast(big); dCr, _ = contrast(rest)
    mb = np.median(dCb, 0); mr = np.median(dCr, 0)
    rho = stats.spearmanr(mb, mr).statistic
    print("  %-17s %d vs %d pairs | Spearman of the 49 per-edge median dC = %+.3f"
          % (kind, len(big), len(rest), rho))
    print("                    %d/49 edges agree in sign" % int(np.sum(np.sign(mb) == np.sign(mr))))

print("\n[done] wrote %s (%d rows)" % (os.path.join(D_OUT, "paired_edge_test.csv"), len(R)))
