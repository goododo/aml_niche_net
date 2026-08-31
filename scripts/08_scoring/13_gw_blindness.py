#!/usr/bin/env python
# 13_gw_blindness.py ----
# INPUT  : <root>/06_distance/edge_distance.csv , <root>/07_fgw/fgw_input_index.csv ,
#          <root>/01_preprocess/00_curated_manifest.csv
# OUTPUT : <root>/08_scoring/gw_blindness_D4.csv  (one row per directed edge per arm)
# WHAT IT DOES : D4 of the blindness diagnostics. Asks whether the 7-bin cost matrix C carries
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
