#!/usr/bin/env python3
"""cooption_redistribution.py ----
Probe for CO-OPTION vs EMERGENCE in the top CCC pathways (MIF, Galectin-9/LGALS9, CD99).

Motivation
----------
The H1 test (`emergent_edges`, dC = C_AML - C_healthy) is a DIFFERENCE test: it only
flags a link when its magnitude changes between AML and healthy marrow. It is therefore
blind to a *co-opted* dependency -- a pathway that is present in BOTH conditions at a
similar level but is "used differently" in AML (redirected to different senders/receivers).
A conserved+co-opted axis produces dC ~= 0, i.e. exactly the null we observed.

This script tests the co-option hypothesis on the three dominant ligands using a
composition measure that is normalized WITHIN each pathway WITHIN each sample, so it is
scale-free (unaffected by cell-count / edge-count differences between AML and healthy):

  For each sample s and pathway L (ligand), over its significant edges (pval < 0.05):
    total(s,L)              = sum(prob)
    primitive_sender_share  = sum(prob | sender in PRIMITIVE) / total          # who sends it
    immune_receiver_share   = sum(prob | receiver in IMMUNE)  / total          # who receives it
    prim_to_immune_share    = sum(prob | sender in PRIMITIVE & receiver in IMMUNE) / total  # SCS axis

We then ask two things per pathway:
  (1) CONSERVED?  prevalence (fraction of samples with >=1 significant edge) in AML vs healthy.
  (2) SHIFTED?    is the composition share distribution different AML vs healthy (Mann-Whitney U)?
                  A conserved-but-shifted pathway = high prevalence in both + a significant
                  share shift toward the primitive/LSC compartment in AML = candidate co-option.
  (3) CONVERGE?   does the LSC/primitive compartment both SEND and RECEIVE the pathway in AML?
                  primitive_sender_share, primitive_receiver_share, and especially
                  primitive_autocrine_share (LSC->LSC) test whether a conserved ligand collapses
                  into an autocrine LSC loop in AML -- the receiver-side evidence for the
                  "convergence onto the same LSC" hypothesis (M10 / blueprint v1.2).

Because healthy and AML samples come from different studies, a pooled test is confounded
by study. We therefore also report the test restricted to datasets that contain BOTH
healthy and AML samples (the only within-study contrast available).

INPUT  : results/tables/05_ccc/tensors/<dataset>/<sample>__ccc_cellchat.csv
         results/tables/07_fgw/fgw_input_index.csv          (healthy flag per sample)
OUTPUT : results/tables/05_ccc/cooption/pathway_summary.csv       (headline: conserved? shifted?)
         results/tables/05_ccc/cooption/sender_composition.csv    (pathway x group x sender_bin mean share)
         results/tables/05_ccc/cooption/prevalence_top_ligands.csv(conserved-presence table, top ligands)
         results/tables/05_ccc/cooption/persample_shares.csv      (raw per-sample shares, for transparency)
Usage  : python3 scripts/10_figures/cooption_redistribution.py
"""
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = pathlib.Path("/FAST/gr10634/gaozy/aml_niche_net")
TENSOR_DIR = ROOT / "results" / "tables" / "05_ccc" / "tensors"
INDEX = ROOT / "results" / "tables" / "07_fgw" / "fgw_input_index.csv"
OUTDIR = ROOT / "results" / "tables" / "05_ccc" / "cooption"

PVAL_CUT = 0.05
TARGET_LIGANDS = ["MIF", "LGALS9", "CD99"]            # Galectin-9 = LGALS9
LIGAND_LABEL = {"MIF": "MIF", "LGALS9": "Galectin-9 (LGALS9)", "CD99": "CD99"}

ALL_BINS = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid",
            "Megakaryocyte", "T_NK", "B_Plasma"]
PRIMITIVE = {"HSC_MPP", "LMPP_GMP"}                    # LSC-like / stem-progenitor compartment
IMMUNE = {"T_NK", "B_Plasma", "Mono_DC"}


def load_index() -> pd.DataFrame:
    idx = pd.read_csv(INDEX, dtype={"dataset": str, "sample": str})
    idx["healthy"] = idx["healthy"].astype(str).str.upper().eq("TRUE")
    idx["group"] = np.where(idx["healthy"], "healthy", "AML")
    return idx[["dataset", "sample", "group"]]


def read_all_tensors() -> pd.DataFrame:
    files = sorted(TENSOR_DIR.glob("*/*__ccc_cellchat.csv"))
    if not files:
        print(f"ERROR: no tensor files under {TENSOR_DIR}", file=sys.stderr)
        raise SystemExit(1)
    frames = []
    for f in files:
        df = pd.read_csv(f)
        frames.append(df)
    ccc = pd.concat(frames, ignore_index=True)
    # Purely-numeric sample IDs (e.g. Petti2019 "508084") are inferred as int64 within a
    # single file; force str so the merge with the (string) index never silently drops them.
    ccc["dataset"] = ccc["dataset"].astype(str)
    ccc["sample"] = ccc["sample"].astype(str)
    ccc = ccc[ccc["pval"] < PVAL_CUT].copy()
    ccc["prob"] = pd.to_numeric(ccc["prob"], errors="coerce")
    ccc = ccc.dropna(subset=["prob"])
    ccc = ccc[ccc["prob"] > 0]
    return ccc


def per_sample_shares(ccc: pd.DataFrame, idx: pd.DataFrame) -> pd.DataFrame:
    """One row per (sample, ligand-present) with composition shares."""
    rows = []
    # denominator of samples per group (for prevalence): every sample in the index
    for (dataset, sample, ligand), g in ccc.groupby(["dataset", "sample", "ligand"]):
        if ligand not in TARGET_LIGANDS:
            continue
        total = g["prob"].sum()
        if total <= 0:
            continue
        is_prim_send = g["sender_bin"].isin(PRIMITIVE)
        is_prim_recv = g["receiver_bin"].isin(PRIMITIVE)
        is_imm_recv = g["receiver_bin"].isin(IMMUNE)
        rows.append({
            "dataset": dataset,
            "sample": sample,
            "ligand": ligand,
            "n_edges": len(g),
            "total_prob": total,
            "primitive_sender_share": g.loc[is_prim_send, "prob"].sum() / total,       # who SENDS it
            "primitive_receiver_share": g.loc[is_prim_recv, "prob"].sum() / total,      # who RECEIVES it
            "primitive_autocrine_share": g.loc[is_prim_send & is_prim_recv, "prob"].sum() / total,  # LSC->LSC loop
            "immune_receiver_share": g.loc[is_imm_recv, "prob"].sum() / total,
            "prim_to_immune_share": g.loc[is_prim_send & is_imm_recv, "prob"].sum() / total,
        })
    ps = pd.DataFrame(rows)
    ps = ps.merge(idx, on=["dataset", "sample"], how="left")
    unmatched = ps["group"].isna().sum()
    if unmatched:
        print(f"WARN: {unmatched} sample-rows not found in index (dropped)", file=sys.stderr)
        ps = ps.dropna(subset=["group"])
    return ps


def mw(a: np.ndarray, b: np.ndarray):
    """Two-sided Mann-Whitney U; returns (p, median_a, median_b) or (nan,...) if too few."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 3 or len(b) < 3:
        return np.nan, (np.median(a) if len(a) else np.nan), (np.median(b) if len(b) else np.nan)
    try:
        _, p = mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        p = np.nan
    return p, float(np.median(a)), float(np.median(b))


def summarize(ps: pd.DataFrame, idx: pd.DataFrame) -> pd.DataFrame:
    n_group = idx.groupby("group")["sample"].nunique().to_dict()      # prevalence denominators
    both_ds = (idx.groupby("dataset")["group"].nunique() == 2)
    both_datasets = sorted(both_ds[both_ds].index)                    # within-study contrast pool

    metrics = ["primitive_sender_share", "primitive_receiver_share", "primitive_autocrine_share",
               "prim_to_immune_share", "immune_receiver_share"]
    out = []
    for lig in TARGET_LIGANDS:
        sub = ps[ps["ligand"] == lig]
        aml = sub[sub["group"] == "AML"]
        hlt = sub[sub["group"] == "healthy"]
        rec = {
            "ligand": lig,
            "label": LIGAND_LABEL[lig],
            "n_AML_present": len(aml),
            "n_healthy_present": len(hlt),
            "prevalence_AML": len(aml) / n_group.get("AML", np.nan),
            "prevalence_healthy": len(hlt) / n_group.get("healthy", np.nan),
        }
        for m in metrics:
            p, med_a, med_h = mw(aml[m].values, hlt[m].values)
            rec[f"{m}__median_AML"] = med_a
            rec[f"{m}__median_healthy"] = med_h
            rec[f"{m}__p_pooled"] = p
            rec[f"{m}__shift"] = ("higher_in_AML" if med_a > med_h else
                                  "higher_in_healthy" if med_a < med_h else "equal")
        # within-study contrast (datasets with both groups)
        subw = sub[sub["dataset"].isin(both_datasets)]
        amlw = subw[subw["group"] == "AML"]
        hltw = subw[subw["group"] == "healthy"]
        rec["n_AML_bothDS"] = len(amlw)
        rec["n_healthy_bothDS"] = len(hltw)
        for m in metrics:
            p, med_a, med_h = mw(amlw[m].values, hltw[m].values)
            rec[f"{m}__p_bothDS"] = p
            rec[f"{m}__median_AML_bothDS"] = med_a
            rec[f"{m}__median_healthy_bothDS"] = med_h
        out.append(rec)
    summ = pd.DataFrame(out)
    summ.attrs["both_datasets"] = both_datasets
    summ.attrs["n_group"] = n_group
    return summ


def sender_composition(ccc: pd.DataFrame, idx: pd.DataFrame) -> pd.DataFrame:
    """Mean within-pathway sender-bin share per (ligand, group, sender_bin) -- for a stacked bar."""
    sub = ccc[ccc["ligand"].isin(TARGET_LIGANDS)].copy()
    # per-sample within-pathway share by sender bin
    tot = sub.groupby(["dataset", "sample", "ligand"])["prob"].transform("sum")
    sub["share"] = sub["prob"] / tot
    persamp = (sub.groupby(["dataset", "sample", "ligand", "sender_bin"])["share"]
               .sum().reset_index())
    # Reindex to the full (sample x ligand) x 7-bin grid so a sender_bin that sends NONE of
    # the pathway in a sample contributes share 0 -- otherwise the per-bin mean is taken over
    # only the samples where that bin fires, and the bins no longer sum to 1 (cannot be stacked).
    keys = persamp[["dataset", "sample", "ligand"]].drop_duplicates()
    grid = keys.merge(pd.DataFrame({"sender_bin": ALL_BINS}), how="cross")
    persamp = grid.merge(persamp, on=["dataset", "sample", "ligand", "sender_bin"], how="left")
    persamp["share"] = persamp["share"].fillna(0.0)
    persamp = persamp.merge(idx, on=["dataset", "sample"], how="left").dropna(subset=["group"])
    comp = (persamp.groupby(["ligand", "group", "sender_bin"])["share"]
            .mean().reset_index().rename(columns={"share": "mean_sender_share"}))
    # keep a stable bin order (means now sum to 1 within each ligand x group)
    comp["sender_bin"] = pd.Categorical(comp["sender_bin"], categories=ALL_BINS, ordered=True)
    comp = comp.sort_values(["ligand", "group", "sender_bin"]).reset_index(drop=True)
    return comp


def prevalence_top(ccc: pd.DataFrame, idx: pd.DataFrame, top_n: int = 15) -> pd.DataFrame:
    """Prevalence (fraction of samples with >=1 sig edge) per ligand, AML vs healthy.

    Directly evidences the 'conserved not emergent' framing: dominant AML ligands that are
    also highly prevalent in healthy marrow are conserved, so H1's difference test cannot
    see them.
    """
    n_group = idx.groupby("group")["sample"].nunique().to_dict()
    pres = (ccc.groupby(["dataset", "sample", "ligand"]).size().reset_index(name="n"))
    pres = pres.merge(idx, on=["dataset", "sample"], how="left").dropna(subset=["group"])
    prev = (pres.groupby(["ligand", "group"])["sample"].nunique().reset_index(name="n_present"))
    prev = prev.pivot(index="ligand", columns="group", values="n_present").fillna(0)
    for grp in ("AML", "healthy"):
        if grp not in prev.columns:
            prev[grp] = 0
    prev["prevalence_AML"] = prev["AML"] / n_group.get("AML", np.nan)
    prev["prevalence_healthy"] = prev["healthy"] / n_group.get("healthy", np.nan)
    prev = prev.sort_values("AML", ascending=False).head(top_n).reset_index()
    prev = prev.rename(columns={"AML": "n_present_AML", "healthy": "n_present_healthy"})
    return prev[["ligand", "n_present_AML", "prevalence_AML",
                 "n_present_healthy", "prevalence_healthy"]]


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    idx = load_index()
    ccc = read_all_tensors()
    print(f"[cooption] {ccc['sample'].nunique()} samples, {len(ccc):,} significant edges (pval<{PVAL_CUT})")

    ps = per_sample_shares(ccc, idx)
    ps.to_csv(OUTDIR / "persample_shares.csv", index=False)

    summ = summarize(ps, idx)
    summ.to_csv(OUTDIR / "pathway_summary.csv", index=False)

    comp = sender_composition(ccc, idx)
    comp.to_csv(OUTDIR / "sender_composition.csv", index=False)

    prev = prevalence_top(ccc, idx)
    prev.to_csv(OUTDIR / "prevalence_top_ligands.csv", index=False)

    # ---- console report ----
    ng = summ.attrs["n_group"]; bds = summ.attrs["both_datasets"]
    print(f"[cooption] group sizes: {ng}")
    print(f"[cooption] datasets with BOTH groups (within-study contrast): {bds}")
    print("\n=== CONSERVED? (prevalence AML vs healthy) ===")
    for _, r in summ.iterrows():
        print(f"  {r['label']:<22} AML {r['prevalence_AML']:.0%}  healthy {r['prevalence_healthy']:.0%}"
              f"   (n present: {r['n_AML_present']} / {r['n_healthy_present']})")
    print("\n=== SHIFTED? (composition share: median AML vs healthy, Mann-Whitney) ===")
    for _, r in summ.iterrows():
        print(f"  {r['label']}")
        for m in ("primitive_sender_share", "primitive_receiver_share",
                  "primitive_autocrine_share", "prim_to_immune_share", "immune_receiver_share"):
            print(f"    {m:<26} AML {r[m+'__median_AML']:.3f}  healthy {r[m+'__median_healthy']:.3f}"
                  f"  p_pooled={r[m+'__p_pooled']:.3g}  p_bothDS={r[m+'__p_bothDS']:.3g}"
                  f"  [{r[m+'__shift']}]")
    print("\n=== CONVERGENCE onto the LSC (autocrine LSC->LSC share): the M10 headline ===")
    for _, r in summ.iterrows():
        m = "primitive_autocrine_share"
        print(f"  {r['label']:<22} AML {r[m+'__median_AML']:.3f}  healthy {r[m+'__median_healthy']:.3f}"
              f"  p_pooled={r[m+'__p_pooled']:.3g}")
    print(f"\n[cooption] wrote 4 tables to {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
