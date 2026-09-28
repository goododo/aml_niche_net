#!/usr/bin/env python
# 23_node_feature_noise_floor.py ----
# INPUT  : 07_fgw__half_A/fgw_nodes_long.csv , 07_fgw__half_B/fgw_nodes_long.csv
#            (built by 05_ccc/03 --split_half_dir --half=<A|B> then 07_fgw/01 --node_features)
#          05_ccc/ccc_sample_manifest.csv                  (patient identity, timepoint)
# OUTPUT : 08_scoring/node_feature_noise_floor.csv          (one row per arm)
# WHAT IT DOES : measures the split-half MEASUREMENT NOISE FLOOR of the per-population node-feature
#          distance, and of cell-type composition, on the same 37 paired samples and with the same
#          statistics 16_scaccordion_noise_floor.py used for the five edge-based arms -- so the
#          numbers are directly comparable to the published 3.9-5.6x.
#
# WHY. 04_node_feature_distance.py found the 47-feature distance beating composition on GATE 1
# (median percentile 0.043 vs 0.174, q 0.004 vs 0.026) and 22_depth_deconfound.py showed the
# advantage survives removing library depth in all four modes. Neither says the distance is
# reproducible. 16 covered ONLY the edge arms: its output file has five rows, none of them prop7 or
# node_feat. So the arm this project's positive claim rests on is the one arm with no floor. That is
# the gap this closes, and it can close the line: a distance whose within-sample retest is as large
# as its between-sample spread is measuring noise, whatever its p-value.
#
# THE ONE DEPARTURE, AND WHY IT IS NECESSARY. FGW_FEATURE_SCALE is "global_z", and 07_fgw/01 applies
# it WITHIN whatever table it is given -- so half A was standardised across half A, and half B
# across half B. Two independently standardised spaces cannot be compared point to point, and
# d(A_i, B_i) would carry that scale difference as if it were noise. Every feature column is
# therefore re-standardised here on the POOLED A+B values before any distance is computed. This
# fixes a shared-scale problem 01 was never asked to solve; it does not re-implement any feature.
# The size of the correction is printed (SELF-CHECK 2) rather than assumed negligible: random halves
# of the same cells should give nearly identical means and SDs, and if they do not, that is itself
# the finding.
#
# Usage : python scripts/08_scoring/23_node_feature_noise_floor.py
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

A = pd.read_csv(os.path.join(ROOT, "07_fgw__half_A", "fgw_nodes_long.csv"))
B = pd.read_csv(os.path.join(ROOT, "07_fgw__half_B", "fgw_nodes_long.csv"))
PANEL = [c for c in A.columns if c not in META
         and c.split("_")[0] in FAMILIES
         and not (c.endswith("_normal") or c.endswith("_malignant"))]
assert len(PANEL) == 47, len(PANEL)
SAMPLES = sorted(set(A["sample"]) & set(B["sample"]))
print("[0] %d samples present in both halves | %d panel features" % (len(SAMPLES), len(PANEL)))

## -- SELF-CHECK 1: both halves must be a complete samples x 7 grid ----
for nm, F in [("A", A), ("B", B)]:
    n = F[F["sample"].isin(SAMPLES)].groupby("sample")["hierarchy_bin"].nunique()
    assert (n == 7).all(), "half %s is not a complete grid: %s" % (nm, n[n != 7].to_dict())
print("[0] both halves are complete %d x 7 grids" % len(SAMPLES))

def wide(F, feats):
    w = F[F["sample"].isin(SAMPLES)].pivot_table(index="sample", columns="hierarchy_bin",
                                                 values=feats, dropna=False)
    return w.reindex(index=SAMPLES, columns=pd.MultiIndex.from_product([feats, CCC_NODES]))

XA, XB = wide(A, PANEL), wide(B, PANEL)
assert list(XA.columns) == list(XB.columns)

## -- SELF-CHECK 2: how much does the shared-scale correction actually change? ----
ma, sa = XA.to_numpy(float).mean(0), XA.to_numpy(float).std(0)
mb, sb = XB.to_numpy(float).mean(0), XB.to_numpy(float).std(0)
print("[SELF-CHECK 2] independent standardisation of the two halves:")
print("    per-column |mean_A - mean_B| : median %.4f, max %.4f" % (np.median(np.abs(ma - mb)), np.abs(ma - mb).max()))
# Guard the ratio: a feature that is constant inside a bin has SD 0 there, which is itself worth
# counting -- a degenerate column contributes no information to the 329-dimensional distance.
deg = (sa == 0) | (sb == 0)
ok = ~deg
print("    degenerate columns (SD 0 in a half): %d of %d" % (int(deg.sum()), len(sa)))
if ok.any():
    rt = sa[ok] / sb[ok]
    print("    per-column SD ratio A/B (non-degenerate): median %.4f, range %.4f-%.4f"
          % (np.median(rt), rt.min(), rt.max()))
P = np.vstack([XA.to_numpy(float), XB.to_numpy(float)])
mu, sd = P.mean(0), np.where(P.std(0) > 0, P.std(0), 1.0)
ZA, ZB = (XA.to_numpy(float) - mu) / sd, (XB.to_numpy(float) - mu) / sd
print("    -> both halves re-standardised on the pooled A+B values")

# composition from the halves' own cell counts; CLR is invariant to the common 1/2 factor
def comp_of(F):
    c = (F[F["sample"].isin(SAMPLES)]
         .pivot_table(index="sample", columns="hierarchy_bin", values="n_cells_raw", fill_value=0)
         .reindex(index=SAMPLES, columns=CCC_NODES).fillna(0.0).to_numpy(float)) + 0.5
    return np.log(c) - np.log(c).mean(1, keepdims=True)
CA, CB = comp_of(A), comp_of(B)

man = pd.read_csv(os.path.join(ROOT, "05_ccc", "ccc_sample_manifest.csv")).set_index("sample")
info = man.reindex(SAMPLES)[["uid_patient", "Timepoint"]]
pat = info.uid_patient.to_numpy()
same_pat = (pat[:, None] == pat[None, :])
n = len(SAMPLES)
iu = np.triu_indices(n, 1)
# within-patient = two DIFFERENT samples of one patient; between = different patients
off_same = same_pat[iu]

rows = []
for arm, Ma, Mb in [("node_feat", ZA, ZB), ("prop7", CA, CB)]:
    # one space, so a retest distance and a between-sample distance are in the same units
    both = np.vstack([Ma, Mb])
    D = squareform(pdist(both, "euclidean"))
    retest = np.array([D[i, n + i] for i in range(n)])          # same sample, half A vs half B
    DA = D[:n, :n]
    between = DA[iu][~off_same]
    within_pat = DA[iu][off_same]
    med_r, med_b = float(np.median(retest)), float(np.median(between))
    med_w = float(np.median(within_pat)) if len(within_pat) else np.nan
    # would the answer be the same from the other half? rank agreement of the two distance matrices
    DB = D[n:, n:]
    rho_ab = float(stats.spearmanr(DA[iu], DB[iu]).statistic)
    w = stats.wilcoxon(retest - np.median(between), alternative="less")
    rows.append(dict(arm=arm, n_samples=n, retest_half_median=med_r,
                     between_half_median=med_b, within_patient_half_median=med_w,
                     snr_between=med_b / med_r if med_r > 0 else np.inf,
                     snr_within_patient=med_w / med_r if med_r > 0 else np.inf,
                     spearman_A_vs_B_matrix=rho_ab, wilcoxon_p=float(w.pvalue),
                     n_within_patient_pairs=int(off_same.sum())))
R = pd.DataFrame(rows)
R.to_csv(os.path.join(ROOT, "08_scoring", "node_feature_noise_floor.csv"), index=False)

print("\n" + "=" * 92)
print("%-11s %6s %10s %10s %10s %9s %9s %11s" % ("arm", "n", "retest", "between", "within-pt",
                                                 "SNR_btw", "SNR_wp", "rho(A,B mat)"))
print("=" * 92)
for _, r in R.iterrows():
    print("%-11s %6d %10.3f %10.3f %10.3f %9.2f %9.2f %11.3f"
          % (r.arm, r.n_samples, r.retest_half_median, r.between_half_median,
             r.within_patient_half_median, r.snr_between, r.snr_within_patient,
             r.spearman_A_vs_B_matrix))
print("\nreference -- 16_scaccordion_noise_floor.py, the five EDGE arms on these same 37 samples:")
f = os.path.join(ROOT, "08_scoring", "scaccordion_noise_floor.csv")
if os.path.exists(f):
    E = pd.read_csv(f)
    for _, r in E.iterrows():
        print("  %-13s SNR_between %.2f | SNR_within_patient %.2f"
              % (r.arm, r.snr_between, r.snr_within_patient))
print("\nREADING: SNR_between near 1 means two halves of ONE sample differ as much as two different")
print("         samples do -- the distance would then be measuring noise. The GATE 1 comparison")
print("         between node_feat and prop7 is only interpretable if BOTH clear that bar.")
print("\n[done] wrote node_feature_noise_floor.csv")
