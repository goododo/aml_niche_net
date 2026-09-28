#!/usr/bin/env python
# 19_gate1_equivalence.py ----
# INPUT  : 08_scoring/scaccordion_gate1_ranks.csv   (per patient per arm percentiles, from 17)
# OUTPUT : 08_scoring/gate1_equivalence.csv
# WHAT IT DOES : puts a confidence bound on HOW MUCH better than composition each communication
#          arm could plausibly be, given only 11 paired patients.
#
# WHY. The GATE 1 comparison "no arm improves on prop7, paired p = 0.23-0.56" is, at n = 11,
# absence of evidence and must not be written as evidence of absence. The honest statement has the
# form "we can rule out median improvements larger than X percentile points", and X comes from a
# bootstrap interval on the median paired difference. The roster composition (9 of 11 patients
# from one dataset) is printed next to it, not hidden under it.
#
# Usage : python scripts/08_scoring/19_gate1_equivalence.py
import os
import numpy as np
import pandas as pd

import argparse
_ap = argparse.ArgumentParser()
_ap.add_argument("--suffix", type=str, default="")
_args = _ap.parse_args()
SFX = ("__" + _args.suffix) if _args.suffix else ""
ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables/08_scoring"
N_BOOT = 10000
SEED = 491638
rng = np.random.default_rng(SEED)

R = pd.read_csv(os.path.join(ROOT, "scaccordion_gate1_ranks%s.csv" % SFX))
out = []
for contrast, g in R.groupby("contrast"):
    wide = g.pivot_table(index=["dataset", "patient"], columns="arm", values="pct")
    wide = wide.dropna(subset=["prop7"])
    ds_counts = wide.reset_index().dataset.value_counts()
    print("=" * 78)
    print("%s : n = %d patients | dataset composition: %s"
          % (contrast, len(wide), ", ".join("%s %d" % kv for kv in ds_counts.items())))
    print("(improvement = prop7 percentile - arm percentile; positive = the arm beats composition)")
    print("%-14s %4s %10s %26s %12s" % ("arm", "n", "obs median", "bootstrap 95% CI", "rule out >"))
    for arm in [c for c in wide.columns if c != "prop7"]:
        d = (wide["prop7"] - wide[arm]).dropna().to_numpy()
        if len(d) < 5:
            continue
        obs = float(np.median(d))
        bs = np.array([np.median(rng.choice(d, len(d), replace=True)) for _ in range(N_BOOT)])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        # the equivalence statement: improvements larger than `hi` are excluded at 95%
        print("%-14s %4d %+10.3f %12s [%+.3f, %+.3f] %11.1f pts"
              % (arm, len(d), obs, "", lo, hi, 100 * hi))
        out.append(dict(contrast=contrast, arm=arm, n=len(d), obs_median_improvement=obs,
                        ci95_lo=float(lo), ci95_hi=float(hi),
                        rule_out_improvement_pct_points=float(100 * hi),
                        n_patients_arm_better=int((d > 0).sum()),
                        top_dataset=str(ds_counts.index[0]), top_dataset_n=int(ds_counts.iloc[0])))
E = pd.DataFrame(out)
E.to_csv(os.path.join(ROOT, "gate1_equivalence%s.csv" % SFX), index=False)
print("\n[done] wrote gate1_equivalence.csv (%d rows)" % len(E))
