#!/usr/bin/env python
# 15_scaccordion_benchmark.py ----
# INPUT  : <root>/06_distance/scaccordion_distance__<arm>.csv   (06_distance/03 output, 6 arms)
#          <root>/06_distance/scaccordion_qc.csv                (the inertness readout)
#          <root>/05_ccc/tensors/<ds>/<sample>__ccc_cellchat.csv (re-read, for planting)
# OUTPUT : <root>/08_scoring/scaccordion_null_ari.csv           STEP 1
#          <root>/08_scoring/scaccordion_positive_control.csv   STEP 2
#          <root>/08_scoring/scaccordion_ari__dataset.csv       STEP 3
#          <root>/08_scoring/scaccordion_ari__disease.csv       STEP 4
#          <root>/08_scoring/scaccordion_support_confound.csv   STEP 5
# WHAT IT DOES : scores the six-rung ablation ladder from 06_distance/03 against the labels, in a
#          LOCKED ORDER registered in PREREGISTRATION_scaccordion.md: the simulated null first,
#          then the planted-effect positive control, then a label known to exist (dataset), and
#          only then the disease label the project actually cares about.
#
# WHY THE ORDER IS LOCKED. Section 11's rule, applied to somebody else's method: a representation
# that cannot see an effect deliberately planted in it may not be used to report that an effect is
# absent. The steps are separate output files and each refuses to run until the previous one is on
# disk, so the sequence is enforced by the filesystem rather than by intention.
#
# WHY ARI IS REPORTED TWICE. scACCorDiON's benchmark scores k-medoids at k = the number of true
# classes. On a cohort with ten datasets, a clustering that recovers only the study of origin
# already reaches a substantial ARI against disease status, so a bare ARI is not interpretable.
# Every ARI here is therefore printed next to (a) the simulated null for that cell and k, and
# (b) the batch ceiling: the best ARI reachable by ANY clustering that knows only `dataset`.
#
# Usage : PYTHONPATH=/FAST/gr10634/gaozy/external/pylibs \
#         python scripts/08_scoring/15_scaccordion_benchmark.py --step {1,2,3,4,5|all}
import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config"))
from scaccordion_geometry import (                                          # noqa: E402
    tables_from_tensors, pmat, variance_filter, stg, htd_cost, emd_matrix,
    plant, scaccordion_on_path, HTD_BETA, TELEPORT, VAR_FILTER_Q)

CCC_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
ARMS = ["prop7", "tabular", "tv", "corrot", "dwot_shipped", "dwot_fixed"]
KS = list(range(2, 8))                     # scACCorDiON's own sweep range
SEED = 491638                              # config_paths.sh SEED, as everywhere else here
N_NULL = 5000

# The SAME six edges 10_planted_effect_power.py plants on, recovered from its own output so the
# two positive controls target identical biology. Four of the six do not survive scACCorDiON's
# default variance filter, which is precisely why STEP 2 runs at filter_q = 0 as well as 0.2.
PINNED6 = [("B_Plasma", "Erythroid"), ("HSC_MPP", "HSC_MPP"), ("Megakaryocyte", "HSC_MPP"),
           ("T_NK", "Erythroid"), ("B_Plasma", "Mono_DC"), ("Megakaryocyte", "T_NK")]
DELTAS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]   # 10_planted_effect_power.py's registered grid

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--step", default="all")
ap.add_argument("--n_null", type=int, default=N_NULL)
args = ap.parse_args()
D_DST = os.path.join(args.root, "06_distance")
D_SCO = os.path.join(args.root, "08_scoring")
D_CCC = os.path.join(args.root, "05_ccc")
STEPS = ARMS and (["1", "2", "3", "4", "5"] if args.step == "all" else [args.step])
rng = np.random.default_rng(SEED)


def kmedoid_labels(Dm, k):
    """Their clusterer: kmedoids.KMedoids(method='fasterpam'), as utils.performance_eval uses."""
    scaccordion_on_path()
    import kmedoids
    return kmedoids.KMedoids(n_clusters=k, method="fasterpam",
                             random_state=SEED).fit(np.asarray(Dm, float)).labels_


def load_D(arm):
    d = pd.read_csv(os.path.join(D_DST, "scaccordion_distance__%s.csv" % arm), index_col=0)
    return d


def roster(samples):
    """Sample -> labels, built from the tensors themselves so it cannot drift from the input."""
    rows = []
    for f in sorted(glob.glob(os.path.join(D_CCC, "tensors", "*", "*__ccc_cellchat.csv"))):
        h = pd.read_csv(f, nrows=1)
        rows.append(dict(sample=str(h["sample"].iloc[0]), dataset=str(h["dataset"].iloc[0]),
                         timepoint=str(h["timepoint"].iloc[0])))
    R = pd.DataFrame(rows).drop_duplicates("sample").set_index("sample").reindex(samples)
    R["disease"] = np.where(R.timepoint == "Healthy", "healthy", "AML")
    R["tp4"] = R.timepoint.replace({"Post_induction": "Treatment", "On_treatment": "Treatment",
                                    "Post_treatment_unspecified": "Treatment",
                                    "Post_consolidation": "Treatment", "Refractory": "Treatment",
                                    "Relapse2": "Relapse"})
    return R


def require(path, step):
    if not os.path.exists(path):
        raise SystemExit("STEP %s refuses to run: %s is not on disk. The order is locked; see "
                         "PREREGISTRATION_scaccordion.md." % (step, os.path.basename(path)))


D0 = load_D(ARMS[0])
SAMPLES = list(D0.index)
R = roster(SAMPLES)
print("[0] %d samples, %d datasets, arms: %s" % (len(SAMPLES), R.dataset.nunique(), ", ".join(ARMS)))

F_NULL = os.path.join(D_SCO, "scaccordion_null_ari.csv")
F_POS = os.path.join(D_SCO, "scaccordion_positive_control.csv")
F_DS = os.path.join(D_SCO, "scaccordion_ari__dataset.csv")
F_DIS = os.path.join(D_SCO, "scaccordion_ari__disease.csv")
F_SUP = os.path.join(D_SCO, "scaccordion_support_confound.csv")

# ================================================================================================
# STEP 1 -- what does chance look like? Label structure only; no distance is touched.
# ================================================================================================
if "1" in STEPS:
    print("\n=== STEP 1: simulated null ARI (%d draws) ===" % args.n_null)
    out = []
    for lname in ["disease", "tp4", "dataset"]:
        y = R[lname].to_numpy()
        a2, am = np.empty(args.n_null), np.empty(args.n_null)
        for t in range(args.n_null):
            vs = [adjusted_rand_score(y, rng.integers(0, k, len(y))) for k in KS]
            a2[t], am[t] = vs[0], max(vs)
        out.append(dict(label=lname, n=len(y), n_classes=len(set(y)),
                        k2_mean=a2.mean(), k2_sd=a2.std(ddof=1), k2_p95=np.percentile(a2, 95),
                        k2_p99=np.percentile(a2, 99), maxk_p95=np.percentile(am, 95),
                        maxk_p99=np.percentile(am, 99)))
        print("    %-8s n=%3d classes=%d  k=2 p95 %.4f  max-over-k p95 %.4f"
              % (lname, len(y), len(set(y)), out[-1]["k2_p95"], out[-1]["maxk_p95"]))
    pd.DataFrame(out).to_csv(F_NULL, index=False)

# ================================================================================================
# STEP 2 -- POSITIVE CONTROL. Plant a known effect at the ligand-receptor level and rebuild the
# entire published construction from the perturbed input. An arm that cannot recover this may not
# be used in STEP 4 to report an absence.
# ================================================================================================
if "2" in STEPS:
    require(F_NULL, "2")
    print("\n=== STEP 2: planted-effect positive control ===")
    null = pd.read_csv(F_NULL).set_index("label")
    paths = sorted(glob.glob(os.path.join(D_CCC, "tensors", "*", "*__ccc_cellchat.csv")))
    tbls = {k: v for k, v in tables_from_tensors(paths).items() if len(v)}
    tbls = {k: tbls[k] for k in SAMPLES}
    half = set(rng.choice(SAMPLES, size=len(SAMPLES) // 2, replace=False))
    ytrue = np.array([1 if s in half else 0 for s in SAMPLES])
    print("    planting on %d pinned edges, %d of %d samples, deltas %s"
          % (len(PINNED6), len(half), len(SAMPLES), DELTAS))
    rows = []
    for fq in (VAR_FILTER_Q, 0.0):
        for delta in DELTAS:
            pt = tbls if delta == 0 else plant(tbls, PINNED6, delta, half)
            Pp = variance_filter(pmat(pt, nodes=CCC_NODES), q=fq)[SAMPLES] if fq > 0 else \
                pmat(pt, nodes=CCC_NODES)[SAMPLES]
            Pp = Pp.loc[Pp.sum(axis=1) != 0, :]
            kept = [e for e in PINNED6 if "%s$%s" % e in set(Pp.index)]
            Ptr, idx = stg(Pp, arm="paper", teleport=TELEPORT)
            C = htd_cost(Ptr, beta=HTD_BETA)
            A = Pp.to_numpy(float); A = A / A.sum(0, keepdims=True)
            bank = {"dwot_fixed": emd_matrix(Pp, C),
                    "tv": squareform(pdist(A.T, "cityblock")) / 2.0,
                    "tabular": squareform(pdist(A.T, "euclidean"))}
            for arm, Dm in bank.items():
                ar = {k: adjusted_rand_score(ytrue, kmedoid_labels(Dm, k)) for k in KS}
                rows.append(dict(filter_q=fq, delta=delta, arm=arm, n_slots=Pp.shape[0],
                                 pinned_surviving=len(kept), ari_k2=ar[2],
                                 ari_maxk=max(ar.values()), best_k=max(ar, key=ar.get)))
            b = [r for r in rows if r["delta"] == delta and r["filter_q"] == fq]
            print("    filter=%.2f delta=%-5g slots=%2d pinned_kept=%d/6 | "
                  % (fq, delta, b[0]["n_slots"], b[0]["pinned_surviving"])
                  + "  ".join("%s k2=%.3f maxk=%.3f" % (r["arm"], r["ari_k2"], r["ari_maxk"])
                              for r in b))
    P = pd.DataFrame(rows)
    P["null_k2_p95"] = float(null.loc["disease", "k2_p95"])
    P["null_maxk_p95"] = float(null.loc["disease", "maxk_p95"])
    P.to_csv(F_POS, index=False)
    ok = P[(P.arm == "dwot_fixed") & (P.delta == max(DELTAS)) & (P.filter_q == 0.0)]
    print("    REGISTERED PASS RULE: dwot_fixed at delta=%g, filter_q=0 must reach ARI >= 0.30."
          % max(DELTAS))
    print("    observed ari_maxk = %.4f -> %s"
          % (ok.ari_maxk.iloc[0], "PASS" if ok.ari_maxk.iloc[0] >= 0.30 else "FAIL"))

# ================================================================================================
# STEP 3 -- a label that is KNOWN to be present. If no arm recovers `dataset`, the distances are
# at the noise floor and STEP 4 cannot distinguish "no effect" from "no measurement".
# ================================================================================================
if "3" in STEPS:
    require(F_POS, "3")
    print("\n=== STEP 3: ARI against `dataset` (a label known to exist) ===")
    y = R["dataset"].to_numpy()
    rows = []
    for arm in ARMS:
        Dm = load_D(arm).to_numpy(float)
        ar = {k: adjusted_rand_score(y, kmedoid_labels(Dm, k)) for k in KS}
        k10 = adjusted_rand_score(y, kmedoid_labels(Dm, R.dataset.nunique()))
        rows.append(dict(arm=arm, label="dataset", ari_k2=ar[2], ari_maxk=max(ar.values()),
                         best_k=max(ar, key=ar.get), ari_at_n_classes=k10))
        print("    %-13s k=2 %.4f   max over k=2..7 %.4f (k=%d)   at k=%d %.4f"
              % (arm, ar[2], max(ar.values()), max(ar, key=ar.get), R.dataset.nunique(), k10))
    pd.DataFrame(rows).to_csv(F_DS, index=False)

# ================================================================================================
# STEP 4 -- THE REGISTERED PRIMARY OUTCOME. Runs last, on purpose.
# ================================================================================================
if "4" in STEPS:
    require(F_DS, "4")
    print("\n=== STEP 4: ARI against disease status (PRIMARY) ===")
    null = pd.read_csv(F_NULL).set_index("label")
    cells = {"full_cohort": np.ones(len(R), bool),
             "GSE185381": (R.dataset == "GSE185381").to_numpy()}
    rows = []
    for cname, mask in cells.items():
        y = R.loc[mask, "disease"].to_numpy()
        if len(set(y)) < 2:
            continue
        for arm in ARMS:
            Dm = load_D(arm).to_numpy(float)[np.ix_(mask, mask)]
            ar = {k: adjusted_rand_score(y, kmedoid_labels(Dm, k)) for k in KS}
            rows.append(dict(cell=cname, n=int(mask.sum()), arm=arm, label="disease",
                             ari_k2=ar[2], ari_maxk=max(ar.values()), best_k=max(ar, key=ar.get),
                             null_k2_p95=float(null.loc["disease", "k2_p95"]),
                             null_maxk_p95=float(null.loc["disease", "maxk_p95"])))
            print("    %-12s %-13s k=2 %.4f   max over k %.4f   (null p95 %.4f / %.4f)"
                  % (cname, arm, ar[2], max(ar.values()),
                     rows[-1]["null_k2_p95"], rows[-1]["null_maxk_p95"]))
    pd.DataFrame(rows).to_csv(F_DIS, index=False)

# ================================================================================================
# STEP 5 -- the companion readout the registration makes mandatory. 86% of this distance's rank
# information is which edges were callable at all, so any separation must be shown NOT to be a
# restatement of per-sample detection sparsity.
# ================================================================================================
if "5" in STEPS:
    require(F_DIS, "5")
    print("\n=== STEP 5: support-sparsity companion ===")
    Pm = pd.read_csv(os.path.join(D_DST, "scaccordion_pmat.csv"), index_col=0)[SAMPLES]
    support = (Pm.to_numpy() > 0).sum(0)
    ysup = (support > np.median(support)).astype(int)
    rows = []
    for arm in ARMS:
        Dm = load_D(arm).to_numpy(float)
        ar_s = max(adjusted_rand_score(ysup, kmedoid_labels(Dm, k)) for k in KS)
        iu = np.triu_indices(len(SAMPLES), 1)
        B = (Pm.to_numpy() > 0).astype(float)
        jac = squareform(pdist(B.T, "jaccard"))
        rows.append(dict(arm=arm, ari_support_split=ar_s,
                         spearman_vs_support_jaccard=float(
                             stats.spearmanr(Dm[iu], jac[iu]).statistic)))
        print("    %-13s ARI vs support median-split %.4f   rho vs support Jaccard %.4f"
              % (arm, ar_s, rows[-1]["spearman_vs_support_jaccard"]))
    d = pd.DataFrame(rows)
    d["support_min"], d["support_median"], d["support_max"] = support.min(), np.median(support), support.max()
    d.to_csv(F_SUP, index=False)
    print("    REGISTERED READING: an arm whose ARI against disease does not exceed its ARI")
    print("    against the support split has not shown a communication result.")

print("\n[done] steps run: %s" % ", ".join(STEPS))
