#!/usr/bin/env python
# 13_gw_blindness.py ----
# INPUT  : <root>/06_distance/edge_distance.csv , <root>/07_fgw/fgw_input_index.csv ,
#          <root>/01_preprocess/00_curated_manifest.csv
# OUTPUT : <root>/08_scoring/gw_blindness_D4.csv  (one row per directed edge per arm)
#          <root>/08_scoring/gw_blindness_D5.csv  (one row per paired sample)
#          <root>/08_scoring/gw_blindness_D1.csv  (one row per sample: GW coupling geometry)
# WHAT IT DOES : the blindness diagnostics. D4 asks whether the 7-bin cost matrix C carries
#          any usable between-sample signal AT ALL, inside the paired samples only. If it does
#          not, no metric on top of C can work and D1/D2/D3 are moot.
#
# WHY THIS RUNS FIRST. 10_planted_effect_power.py showed that changing 98% of the targeted
# entries of C by a mean of 0.64 leaves the alpha=1 FGW score unmoved, in all four registered
# arms. Two explanations survive: the FGW statistic throws the information away, or C never
# had discriminative structure to begin with. D4 tests the second, and it needs no optimal
# transport at all, only arithmetic. It is therefore the cheapest way to possibly end the
# whole line of inquiry.
#
# SAMPLE SET. Paired samples only, per the 2026-08-31 instruction. The roster is printed in
# full before any number is reported. Pairs are built with the SAME rule as 11_paired_gate.py
# and 12_paired_planting.py; a mismatch there would make the three incomparable, so the rule
# is re-derived from the shared vocabulary rather than re-listed.
#
# Usage : python scripts/08_scoring/13_gw_blindness.py [--arms rank,mask,const,logs]
import os, sys, argparse
from scipy import stats
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config"))
from distance_variants import weights_to_C, ARMS as DV_ARMS

FGW_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
EDGE_NAMES = ["%s->%s" % (a, b) for a in FGW_NODES for b in FGW_NODES]
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
REL_TP = {"Relapse", "Relapse2"}
TRT_TP = {"On_treatment", "Post_induction", "Post_consolidation",
          "Post_treatment_unspecified", "Refractory"}

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=DEFAULT_ROOT)
ap.add_argument("--arms", default="rank,mask,const,logs")
args = ap.parse_args()
ARMS = [a.strip() for a in args.arms.split(",")]
for a in ARMS:
    if a not in DV_ARMS:
        raise SystemExit("unknown arm %r; expected from %s" % (a, DV_ARMS))

D_DST = os.path.join(args.root, "06_distance")
D_FGW = os.path.join(args.root, "07_fgw")
D_PRE = os.path.join(args.root, "01_preprocess")
D_OUT = os.path.join(args.root, "08_scoring")

edges = pd.read_csv(os.path.join(D_DST, "edge_distance.csv"))
index = pd.read_csv(os.path.join(D_FGW, "fgw_input_index.csv"))
man = (pd.read_csv(os.path.join(D_PRE, "00_curated_manifest.csv"))
       [["dataset", "sample", "patient_id"]].drop_duplicates(["dataset", "sample"]))
idx = index.merge(man, on=["dataset", "sample"])

## -- the pairs, same rule as 11_paired_gate.py / 12_paired_planting.py ----
pairs = []
for (ds, pt), g in idx.groupby(["dataset", "patient_id"]):
    tps = dict(zip(g.timepoint, zip(g.dataset, g["sample"])))
    if "Diagnosis" not in tps:
        continue
    for kind, want in [("Dx_to_Relapse", REL_TP), ("Dx_to_Treatment", TRT_TP)]:
        hit = sorted([t for t in g.timepoint if t in want])
        if hit:
            pairs.append(dict(kind=kind, patient=pt, dataset=ds,
                              a=tps["Diagnosis"], b=tps[hit[0]], tp_b=hit[0]))
P = pd.DataFrame(pairs)
if P.empty:
    raise SystemExit("no pairs found; the timepoint vocabulary has drifted")

## -- ROSTER. Printed in full, before any result. ----
PAIRED_KEYS = sorted({k for _, r in P.iterrows() for k in (r.a, r.b)})
print("=" * 78)
print("SAMPLE SET: PAIRED SAMPLES ONLY")
print("=" * 78)
print("%d pairs (%s) drawn from %d distinct patients in %d datasets,"
      % (len(P), ", ".join("%s %d" % (k, v) for k, v in P.kind.value_counts().items()),
         P.patient.nunique(), P.dataset.nunique()))
print("covering %d distinct samples. The other %d samples in the index are NOT used here."
      % (len(PAIRED_KEYS), len(idx) - len(PAIRED_KEYS)))
# The pairing rule lives HERE and is written out, so 05_ccc/06_split_half_reliability.R reads the
# roster instead of re-implementing the rule. Two copies of a rule is how the timepoint vocabulary
# and the family map both silently diverged in this project.
os.makedirs(D_OUT, exist_ok=True)
pd.DataFrame([{"dataset": k[0], "sample": k[1]} for k in PAIRED_KEYS]).to_csv(
    os.path.join(D_OUT, "paired_roster.csv"), index=False)
print("roster written to %s" % os.path.join(D_OUT, "paired_roster.csv"))
print("\n%-14s %-12s %-18s %-16s %s" % ("dataset", "patient", "kind", "diagnosis", "second timepoint"))
for _, r in P.sort_values(["dataset", "patient", "kind"]).iterrows():
    print("%-14s %-12s %-18s %-16s %s (%s)" % (r.dataset, r.patient, r.kind, r.a[1], r.b[1], r.tp_b))
print()

## -- C for the paired samples, per arm. SELF-CHECK: the rank arm must reproduce stored C. ----
W = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                       values="weight_probsum", fill_value=0.0)
     .reindex(columns=pd.MultiIndex.from_product([FGW_NODES, FGW_NODES])).fillna(0.0))
NLR = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                         values="n_lr_sig", fill_value=0.0)
       .reindex(columns=pd.MultiIndex.from_product([FGW_NODES, FGW_NODES])).fillna(0.0))
CST = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                         values="C", fill_value=1.0)
       .reindex(columns=pd.MultiIndex.from_product([FGW_NODES, FGW_NODES])).fillna(1.0))
missing = [k for k in PAIRED_KEYS if k not in W.index]
if missing:
    raise SystemExit("paired samples absent from edge_distance.csv: %s" % missing)
W, NLR, CST = W.loc[PAIRED_KEYS], NLR.loc[PAIRED_KEYS], CST.loc[PAIRED_KEYS]

C = {}
for arm in ARMS:
    C[arm] = np.vstack([weights_to_C(W.iloc[i].to_numpy(float), arm, NLR.iloc[i].to_numpy(float))
                        for i in range(len(W))])
dev = np.abs(C["rank"] - CST.to_numpy(float)).max() if "rank" in C else np.nan
print("[SELF-CHECK] rank arm vs stored C over the paired samples: max |dev| = %.3e" % dev)
if "rank" in C and dev > 1e-9:
    raise SystemExit("the rank transform no longer reproduces 06_distance/01; stop and fix that first")

## -- D4a. How many edges even exist in these samples? ----
det = (W.to_numpy(float) > 0)
occ = det.mean(0)
print("\n" + "=" * 78)
print("D4a  EDGE OCCUPANCY over the %d paired samples" % len(W))
print("=" * 78)
print("  detected in every sample      : %d / 49" % int((occ == 1).sum()))
print("  detected in >= 90%% of samples : %d / 49" % int((occ >= 0.9).sum()))
print("  detected in <= 10%% of samples : %d / 49" % int((occ <= 0.1).sum()))
print("  never detected in any sample  : %d / 49" % int((occ == 0).sum()))
print("  median occupancy              : %.1f%%" % (100 * np.median(occ)))

## -- D4b. Within-pair change vs between-patient spread, per edge, per arm. ----
# The paired design can only detect a change that is large relative to how much C differs
# between patients to begin with. This ratio is that comparison, and it needs no model.
ia = {k: i for i, k in enumerate(PAIRED_KEYS)}
rows = []
for arm in ARMS:
    Ca = C[arm]
    for kind in sorted(P.kind.unique()):
        sub = P[P.kind == kind]
        d_within = np.abs(np.vstack([Ca[ia[r.b]] - Ca[ia[r.a]] for _, r in sub.iterrows()]))
        dx = np.vstack([Ca[ia[r.a]] for _, r in sub.iterrows()])          # diagnosis samples
        sd_between = dx.std(0, ddof=1)
        for j, e in enumerate(EDGE_NAMES):
            rows.append(dict(arm=arm, kind=kind, edge=e, occupancy=occ[j],
                             median_within_pair_abs_dC=float(np.median(d_within[:, j])),
                             mean_within_pair_abs_dC=float(d_within[:, j].mean()),
                             sd_between_patients_at_dx=float(sd_between[j]),
                             snr=float(np.median(d_within[:, j]) / sd_between[j])
                                 if sd_between[j] > 0 else np.nan))
R = pd.DataFrame(rows)
os.makedirs(D_OUT, exist_ok=True)
R.to_csv(os.path.join(D_OUT, "gw_blindness_D4.csv"), index=False)

print("\n" + "=" * 78)
print("D4b  IS THERE ANYTHING FOR A PAIRED TEST TO SEE?")
print("=" * 78)
print("snr = median within-pair |dC|  /  between-patient SD of C at diagnosis.")
print("snr near 0 means the within-patient change is invisible against patient-to-patient")
print("spread, and no metric built on C can recover it.\n")
print("%-7s %-17s %10s %10s %10s %14s" % ("arm", "kind", "med snr", "p90 snr", "edges>1", "dead edges"))
for arm in ARMS:
    for kind in sorted(P.kind.unique()):
        g = R[(R.arm == arm) & (R.kind == kind)]
        dead = int((g.sd_between_patients_at_dx < 1e-12).sum())
        s = g.snr.dropna()
        print("%-7s %-17s %10.3f %10.3f %10d %14d"
              % (arm, kind, s.median(), s.quantile(0.9), int((s > 1).sum()), dead))

## -- D4c. Where does the variance of C actually live? ----
# Patient identity vs timepoint. If patient swamps timepoint, a paired test is the right
# design but the effect it chases is small; if BOTH are near zero, C is simply flat.
print("\n" + "=" * 78)
print("D4c  VARIANCE OF C: patient vs timepoint")
print("=" * 78)
print("%-7s %-17s %12s %12s %12s" % ("arm", "kind", "var(C)", "frac patient", "frac timepoint"))
for arm in ARMS:
    Ca = C[arm]
    for kind in sorted(P.kind.unique()):
        sub = P[P.kind == kind].reset_index(drop=True)
        A = np.vstack([Ca[ia[r.a]] for _, r in sub.iterrows()])
        B = np.vstack([Ca[ia[r.b]] for _, r in sub.iterrows()])
        X = np.stack([A, B])                       # (2 timepoints, n pairs, 49 edges)
        gm = X.mean((0, 1))
        tot = ((X - gm) ** 2).mean()
        pat = ((X.mean(0) - gm) ** 2).mean()        # patient main effect
        tim = ((X.mean(1) - gm) ** 2).mean()        # timepoint main effect
        print("%-7s %-17s %12.5f %11.1f%% %13.1f%%"
              % (arm, kind, tot, 100 * pat / tot if tot > 0 else 0,
                 100 * tim / tot if tot > 0 else 0))

print("\n[done] wrote %s (%d rows)" % (os.path.join(D_OUT, "gw_blindness_D4.csv"), len(R)))


## =====================================================================================
## D5. SPLIT-HALF RELIABILITY OF C
## =====================================================================================
# D4b measured a median within-pair |dC| as large as the between-patient SD, with only a
# 3-5% shared timepoint effect. Large-but-directionless is what real idiosyncratic biology
# looks like AND what a noisy estimate looks like. D5 separates them: two halves of ONE
# sample differ by measurement error alone, so their |dC| IS the noise floor. If the
# within-patient |dC| of D4b does not clear that floor, D4b measured noise.
#
# Inputs are written by 05_ccc/06_split_half.sbatch (which calls 02_run_cellchat.py's R
# runner with --cell_subset / --out_suffix). Absent -> D5 is skipped with a message.
D_SPLIT = {h: os.path.join(args.root, "05_ccc", "tensors__split_%s" % h) for h in ("A", "B")}
D_PROD = os.path.join(args.root, "05_ccc", "tensors")

# The p-value threshold is READ FROM THE R CONFIG, not restated here. A second copy of this
# number is exactly how the timepoint vocabulary and the family map silently diverged.
_cfg = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config", "config_distance.R")
import re
_m = re.search(r"^DIST_PVAL_THRESH\s*<-\s*([0-9.]+)", open(_cfg).read(), re.M)
if not _m:
    raise SystemExit("cannot read DIST_PVAL_THRESH from %s" % _cfg)
PVAL_THRESH = float(_m.group(1))
MIN_CELLS = 10   # CCC_MIN_CELLS_PER_NODE; verified against config_ccc.R below
_mc = re.search(r"^CCC_MIN_CELLS_PER_NODE\s*<-\s*(\d+)L?",
                open(os.path.join(os.path.dirname(_cfg), "config_ccc.R")).read(), re.M)
if _mc and int(_mc.group(1)) != MIN_CELLS:
    raise SystemExit("CCC_MIN_CELLS_PER_NODE is %s in config, %d here" % (_mc.group(1), MIN_CELLS))

_PAIR_IDX = [(a, b) for a in FGW_NODES for b in FGW_NODES]


def tensor_to_wnl(path):
    """One LR tensor -> (weight_probsum, n_lr_sig, bins_present) on the fixed 49-edge grid.

    This replicates 06_distance/01_intensity_to_distance.R. It is SELF-CHECKED below against
    the production edge_distance.csv, so the replication cannot drift unnoticed."""
    t = pd.read_csv(path)
    sig = t[t.pval < PVAL_THRESH]
    w = np.zeros(49)
    nl = np.zeros(49)
    if len(sig):
        agg = sig.groupby(["sender_bin", "receiver_bin"])["prob"].agg(["sum", "size"])
        for j, key in enumerate(_PAIR_IDX):
            if key in agg.index:
                w[j] = agg.loc[key, "sum"]
                nl[j] = agg.loc[key, "size"]
    # a bin counts as present if CellChat ever reported it at >= MIN_CELLS on either side
    present = set()
    for col, ncol in (("sender_bin", "n_sender"), ("receiver_bin", "n_receiver")):
        if ncol in t.columns:
            ok = t.loc[t[ncol] >= MIN_CELLS, col]
            present |= set(ok.unique())
    return w, nl, present


missing = {h: [k for k in PAIRED_KEYS
               if not os.path.exists(os.path.join(D_SPLIT[h], k[0], "%s__ccc_cellchat.csv" % k[1]))]
           for h in ("A", "B")}
if missing["A"] or missing["B"]:
    print("\n" + "=" * 78)
    print("D5  SKIPPED: %d of 37 half-A and %d of 37 half-B tensors are not on disk yet."
          % (len(missing["A"]), len(missing["B"])))
    print("    Run: sbatch --array=1-74%8 scripts/05_ccc/06_split_half.sbatch")
    print("=" * 78)
else:
    ## -- SELF-CHECK: the python replication must reproduce the production C exactly ----
    devs = []
    for k in PAIRED_KEYS:
        w, nl, _ = tensor_to_wnl(os.path.join(D_PROD, k[0], "%s__ccc_cellchat.csv" % k[1]))
        devs.append(np.abs(weights_to_C(w, "rank", nl) - CST.loc[k].to_numpy(float)).max())
    dv = float(np.max(devs))
    print("\n" + "=" * 78)
    print("D5  SPLIT-HALF RELIABILITY OF C")
    print("=" * 78)
    print("[SELF-CHECK] python re-aggregation of the PRODUCTION tensors vs stored edge_distance.csv:")
    print("             max |dev| over %d samples = %.3e" % (len(PAIRED_KEYS), dv))
    if dv > 1e-9:
        raise SystemExit("the python aggregation does not reproduce 06_distance/01; D5 would be "
                         "measuring the re-implementation, not the data")
    print("             [PASS] the replication is exact, so half-tensors are on the same footing.\n")

    ## -- reliability, restricted to edges both halves could actually see ----
    d5 = []
    TESTABLE = {}
    HALF_C = {}
    for k in PAIRED_KEYS:
        wa, nla, pa = tensor_to_wnl(os.path.join(D_SPLIT["A"], k[0], "%s__ccc_cellchat.csv" % k[1]))
        wb, nlb, pb = tensor_to_wnl(os.path.join(D_SPLIT["B"], k[0], "%s__ccc_cellchat.csv" % k[1]))
        Ca = weights_to_C(wa, "rank", nla)
        Cb = weights_to_C(wb, "rank", nlb)
        both = pa & pb
        testable = np.array([(a in both) and (b in both) for a, b in _PAIR_IDX])
        rho = (stats.spearmanr(wa[testable], wb[testable]).statistic
               if testable.sum() >= 3 and np.ptp(wa[testable]) > 0 and np.ptp(wb[testable]) > 0
               else np.nan)
        TESTABLE[k] = testable
        HALF_C[k] = (Ca, Cb)
        d5.append(dict(dataset=k[0], sample=k[1],
                       n_bins_A=len(pa), n_bins_B=len(pb), n_bins_both=len(both),
                       n_testable_edges=int(testable.sum()),
                       spearman_w=rho,
                       median_abs_dC_halves=float(np.median(np.abs(Ca - Cb)[testable]))
                           if testable.any() else np.nan,
                       detect_agree=float((( wa > 0) == (wb > 0))[testable].mean())
                           if testable.any() else np.nan))
    D5 = pd.DataFrame(d5)
    D5.to_csv(os.path.join(D_OUT, "gw_blindness_D5.csv"), index=False)

    print("Restricted to edges whose sender AND receiver clear %d cells in BOTH halves," % MIN_CELLS)
    print("so node dropout is reported separately rather than counted as disagreement.\n")
    print("  samples                        : %d" % len(D5))
    print("  bins surviving in both halves  : median %.0f of 7 (production graphs use 7)"
          % D5.n_bins_both.median())
    print("  testable edges                 : median %.0f of 49" % D5.n_testable_edges.median())
    print("  Spearman(w_A, w_B)             : median %.3f   [IQR %.3f, %.3f]"
          % (D5.spearman_w.median(), D5.spearman_w.quantile(.25), D5.spearman_w.quantile(.75)))
    print("  edge detected/absent agreement : median %.1f%%" % (100 * D5.detect_agree.median()))

    ## -- THE DECISIVE COMPARISON ----
    # noise floor  = |dC| between two halves of the SAME sample
    # claimed signal = |dC| between two timepoints of the SAME patient (D4b)
    noise = D5.set_index(["dataset", "sample"]).median_abs_dC_halves
    print("\n" + "-" * 78)
    print("NOISE FLOOR vs THE WITHIN-PATIENT CHANGE D4b REPORTED  (rank arm)")
    print("-" * 78)
    print("%-17s %14s %14s %10s %10s %8s"
          % ("kind", "within-patient", "noise floor", "ratio", "n pairs", "n edges"))
    Ca_rank = C["rank"]
    for kind in sorted(P.kind.unique()):
        sub = P[P.kind == kind]
        sig_v, noi_v, nedge = [], [], []
        for _, r in sub.iterrows():
            # SAME EDGES ON BOTH SIDES. A pair's testable set is the intersection of the two
            # samples' testable sets; the noise floor for that pair is recomputed on it too.
            m = TESTABLE[r.a] & TESTABLE[r.b]
            if m.sum() < 3:
                continue
            nedge.append(int(m.sum()))
            sig_v.append(np.median(np.abs(Ca_rank[ia[r.b]] - Ca_rank[ia[r.a]])[m]))
            for kk in (r.a, r.b):
                noi_v.append(np.median(np.abs(HALF_C[kk][0] - HALF_C[kk][1])[m]))
        s, n_ = float(np.median(sig_v)), float(np.median(noi_v))
        print("%-17s %14.4f %14.4f %10.2f %10d %8.0f" % (kind, s, n_, s / n_ if n_ > 0 else np.inf,
                                                         len(sig_v), np.median(nedge)))
    print("\nratio <= 1 means the within-patient change does not clear the measurement noise,")
    print("and D4b's snr would then be a statement about CellChat's variance, not biology.")
    print("\n[done] wrote %s" % os.path.join(D_OUT, "gw_blindness_D5.csv"))


## =====================================================================================
## D1. WHAT DOES THE GW COUPLING ACTUALLY LOOK LIKE?   (+ D3, the Frobenius baseline)
## =====================================================================================
# The planting harness showed the alpha=1 FGW score does not move when 98% of the targeted
# entries of C change by a mean of 0.64. Two mechanisms could do that, and they need opposite
# fixes, so they have to be told apart by MEASUREMENT rather than argued about:
#
#   mode A, REARRANGEMENT : the optimal coupling T drifts toward a permutation, so a change
#           planted on edge (i,j) is absorbed by re-matching which node plays which role.
#           Signature: T far from diag(p) AND far from the product coupling.
#   mode B, DEGENERACY    : T collapses toward the product coupling p (x) q, which is the
#           uninformative solution. Then GW only sees the DISTRIBUTION of values in C and is
#           blind to where in the matrix they sit. Signature: T close to p (x) q.
#
# D3 then asks the practical question the two modes share: does the simplest possible
# alternative, a plain mass-weighted squared difference between the two cost matrices with NO
# transport at all, separate AML from healthy where GW does not?
#
# SAMPLE SET. This section CANNOT be paired-only. The score is a distance to the healthy
# barycentre, and every paired sample is AML, so the healthy donors have to be loaded to build
# the reference. The paired samples are reported as their own block; the healthy block is the
# reference, not a result.
import ot

D_PRE_ = os.path.join(args.root, "01_preprocess")
nodes_all = pd.read_csv(os.path.join(D_FGW, "fgw_nodes_long.csv"))
EPS_MASS = 1e-6
ALL_KEYS = [(r.dataset, r["sample"]) for _, r in index.iterrows()]

NODE_P = {}
for (ds, smp), g in nodes_all.groupby(["dataset", "sample"]):
    gg = g.set_index("hierarchy_bin").reindex(FGW_NODES)
    pm = np.nan_to_num(gg["mass"].to_numpy(float), nan=EPS_MASS)
    NODE_P[(ds, smp)] = pm / pm.sum()

W_ALL = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                           values="weight_probsum", fill_value=0.0)
         .reindex(columns=pd.MultiIndex.from_product([FGW_NODES, FGW_NODES])).fillna(0.0))
NL_ALL = (edges.pivot_table(index=["dataset", "sample"], columns=["sender_bin", "receiver_bin"],
                            values="n_lr_sig", fill_value=0.0)
          .reindex(columns=pd.MultiIndex.from_product([FGW_NODES, FGW_NODES])).fillna(0.0))
ALL_KEYS = [k for k in ALL_KEYS if k in W_ALL.index and k in NODE_P]
CMAT = {k: weights_to_C(W_ALL.loc[k].to_numpy(float), "rank",
                        NL_ALL.loc[k].to_numpy(float)).reshape(7, 7) for k in ALL_KEYS}

sparse = nodes_all.groupby(["dataset", "sample"])["sparse_flag"].first().fillna(False).to_dict()
is_heal = index.set_index(["dataset", "sample"])["timepoint"].eq("Healthy").to_dict()
heal_bar = [k for k in ALL_KEYS if is_heal.get(k) and not sparse.get(k, False)]
q = np.ones(7) / 7.0

print("\n" + "=" * 78)
print("D1  GW COUPLING GEOMETRY  (+ D3 Frobenius baseline)")
print("=" * 78)
print("healthy barycentre built from %d healthy donors (sparse-flagged excluded);" % len(heal_bar))
print("%d samples scored against it, of which %d are the paired samples.\n"
      % (len(ALL_KEYS), len(PAIRED_KEYS)))

Cs = [CMAT[k] for k in heal_bar]
out = ot.gromov.fgw_barycenters(7, [np.zeros((7, 1)) for _ in heal_bar], Cs,
                                [NODE_P[k] for k in heal_bar],
                                lambdas=[1.0 / len(heal_bar)] * len(heal_bar), alpha=1.0,
                                loss_fun="square_loss", symmetric=False, max_iter=1000,
                                p=q, init_C=np.mean(np.stack(Cs), 0),
                                init_X=np.zeros((7, 1)), random_state=491638, log=True)
Cb = np.asarray(out[1])

d1 = []
for k in ALL_KEYS:
    p_ = NODE_P[k]
    M = np.zeros((7, 7))                      # alpha=1: the feature cost never enters
    T = ot.gromov.fused_gromov_wasserstein(M, CMAT[k], Cb, p_, q, loss_fun="square_loss",
                                           alpha=1.0, symmetric=False)
    T = np.asarray(T)
    prod = np.outer(p_, q)                    # the uninformative coupling
    diagT = np.diag(p_ * 1.0)                 # the identity coupling, scaled to the marginals
    # normalised so the two distances are on one scale: 0 = identical, 1 = as far as prod is
    # from diag, which is the only natural yardstick available inside a 7x7 simplex
    ref = np.linalg.norm(prod - diagT)
    frob = float(((CMAT[k] - Cb) ** 2 * np.outer(p_, p_)).sum())   # D3, no transport at all
    d1.append(dict(dataset=k[0], sample=k[1],
                   healthy=bool(is_heal.get(k)), paired=k in set(PAIRED_KEYS),
                   diag_mass=float(np.trace(T)),
                   dist_to_identity=float(np.linalg.norm(T - diagT) / ref) if ref > 0 else np.nan,
                   dist_to_product=float(np.linalg.norm(T - prod) / ref) if ref > 0 else np.nan,
                   frobenius_score=frob))
D1 = pd.DataFrame(d1)
D1.to_csv(os.path.join(D_OUT, "gw_blindness_D1.csv"), index=False)

print("A coupling that re-matches nodes sits far from BOTH references (mode A).")
print("A coupling that has collapsed to the uninformative solution sits at ~0 from the")
print("product coupling (mode B). Distances are scaled so 1.0 = the gap between the two.\n")
print("%-22s %6s %11s %14s %14s" % ("block", "n", "diag mass", "-> identity", "-> product"))
for lab, m in [("healthy (reference)", D1.healthy),
               ("paired samples", D1.paired),
               ("all AML", ~D1.healthy)]:
    g = D1[m]
    print("%-22s %6d %11.3f %14.3f %14.3f"
          % (lab, len(g), g.diag_mass.median(), g.dist_to_identity.median(),
             g.dist_to_product.median()))

## -- D3. does the transport-free baseline separate AML from healthy? ----
h2 = pd.read_csv(os.path.join(D_OUT, "h2_blast_regression.csv"))
h2 = h2[h2.mass_mode == "ncells"][["dataset", "sample", "blast_proxy", "is_aml"]]
M3 = D1.merge(h2, on=["dataset", "sample"], how="inner")
print("\n" + "-" * 78)
print("D3  FROBENIUS BASELINE: sum_ik (C_sample - C_bary)^2 * p_i * p_k, no transport")
print("-" * 78)
print("  n = %d samples (%d AML, %d healthy)"
      % (len(M3), int((M3.is_aml == 1).sum()), int((M3.is_aml == 0).sum())))
print("  mean healthy %.4f | mean AML %.4f"
      % (M3.loc[M3.is_aml == 0, "frobenius_score"].mean(),
         M3.loc[M3.is_aml == 1, "frobenius_score"].mean()))
rng1 = np.random.default_rng(491638)
y = M3.frobenius_score.to_numpy(float); a = M3.is_aml.to_numpy(float)
dums = pd.get_dummies(M3["dataset"], drop_first=True).to_numpy(float)
X0 = np.column_stack([np.ones(len(M3)), M3.blast_proxy.to_numpy(float), dums])
Qm, _ = np.linalg.qr(X0)
rs = lambda v: v - Qm @ (Qm.T @ v)
ry, ra = rs(y), rs(a)
beta = float(ry @ ra / (ra @ ra))
grp = M3["dataset"].to_numpy()
gidx = [np.where(grp == gg)[0] for gg in np.unique(grp)]
ex = 0
NP = 20000
for _ in range(NP):
    ap_ = a.copy()
    for gi in gidx:
        ap_[gi] = rng1.permutation(ap_[gi])
    rp = rs(ap_)
    dd = float(rp @ rp)
    if dd > 0 and abs(float(ry @ rp / dd)) >= abs(beta):
        ex += 1
print("  within-dataset beta = %+.5f, permutation p = %.4f  (%d draws)"
      % (beta, (1.0 + ex) / (NP + 1.0), NP))
print("\n  compare: alpha=1 FGW on the same cohort gave beta -0.0015, p 0.966 (rank arm).")
print("\n[done] wrote %s" % os.path.join(D_OUT, "gw_blindness_D1.csv"))
