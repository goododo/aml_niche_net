#!/usr/bin/env python3
"""agg_ccc_lr_summary.py ----
Aggregate the 148 per-sample CellChat tensors into small cohort-level summary
tables for the lab-meeting CCC figure (what the communication graphs actually contain).

INPUT  : results/tables/05_ccc/tensors/<dataset>/<sample>__ccc_cellchat.csv
         results/tables/05_ccc/ccc_edge_distance.csv     (bin-bin weight per sample)
         results/tables/07_fgw/fgw_input_index.csv       (healthy flag per sample)
OUTPUT : results/tables/05_ccc/ccc_bin_strength_cohort.csv   (7x7 who-talks-to-whom, AML)
         results/tables/05_ccc/ccc_lr_overall.csv           (top ligand-receptor pairs, AML)
         results/tables/05_ccc/ccc_lr_scs_axis.csv          (top LSC/primitive -> immune L-R, AML)
Scope  : AML samples only (healthy == FALSE); significant interactions only (pval < 0.05).
Usage  : python3 scripts/10_figures/agg_ccc_lr_summary.py
"""
import csv
import glob
import pathlib
import statistics as st
from collections import defaultdict

ROOT = pathlib.Path("/FAST/gr10634/gaozy/aml_niche_net")
TBL = ROOT / "results" / "tables"
CCC = TBL / "05_ccc"

NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
PRIMITIVE = {"HSC_MPP", "LMPP_GMP"}          # CCC_SENDER_LSC
IMMUNE = {"T_NK", "B_Plasma", "Mono_DC"}     # CCC_RECV_IMMUNE
PVAL_MAX = 0.05
MIN_SAMPLES = 12                             # keep L-R rows present in >= this many AML samples


def load_healthy():
    """(dataset, sample) -> True/False from fgw_input_index; default AML if absent."""
    h = {}
    with open(TBL / "07_fgw" / "fgw_input_index.csv") as fh:
        for d in csv.DictReader(fh):
            h[(d["dataset"], d["sample"])] = (d["healthy"].strip().upper() == "TRUE")
    return h


def main():
    healthy = load_healthy()

    def is_aml(ds, sm):
        return not healthy.get((ds, sm), False)   # unknown -> treat as AML

    # ---- Panel A: 7x7 who-talks-to-whom (from bin-level edge table) --------------
    binstr = defaultdict(lambda: {"w": [], "nlr": [], "samp": set()})
    with open(CCC / "ccc_edge_distance.csv") as fh:
        for d in csv.DictReader(fh):
            if not is_aml(d["dataset"], d["sample"]):
                continue
            key = (d["sender_bin"], d["receiver_bin"])
            binstr[key]["w"].append(float(d["weight_probsum"]))
            binstr[key]["nlr"].append(float(d["n_lr_sig"]))
            binstr[key]["samp"].add((d["dataset"], d["sample"]))
    aml_edge_samples = len({s for v in binstr.values() for s in v["samp"]})

    with open(CCC / "ccc_bin_strength_cohort.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sender_bin", "receiver_bin", "mean_weight", "mean_n_lr_sig", "n_samples"])
        for s in NODES:
            for r in NODES:
                v = binstr.get((s, r))
                if v:
                    w.writerow([s, r, f"{st.mean(v['w']):.6f}", f"{st.mean(v['nlr']):.4f}", len(v["samp"])])
                else:
                    w.writerow([s, r, "0", "0", "0"])

    # ---- Panels from L-R tensors -------------------------------------------------
    files = sorted(glob.glob(str(CCC / "tensors" / "*" / "*.csv")))
    lr_full = defaultdict(lambda: {"samp": set(), "prob": []})   # (s,r,lig,rec)
    lr_over = defaultdict(lambda: {"samp": set(), "prob": []})   # (lig,rec)
    n_aml_tensor = set()
    for fp in files:
        with open(fp) as fh:
            for d in csv.DictReader(fh):
                if float(d["pval"]) >= PVAL_MAX:
                    continue
                if not is_aml(d["dataset"], d["sample"]):
                    continue
                sm = (d["dataset"], d["sample"])
                n_aml_tensor.add(sm)
                s, r, lg, rc = d["sender_bin"], d["receiver_bin"], d["ligand"], d["receptor"]
                p = float(d["prob"])
                lr_full[(s, r, lg, rc)]["samp"].add(sm); lr_full[(s, r, lg, rc)]["prob"].append(p)
                lr_over[(lg, rc)]["samp"].add(sm); lr_over[(lg, rc)]["prob"].append(p)
    n_aml = len(n_aml_tensor)

    # overall top L-R pairs (across all bin pairs)
    with open(CCC / "ccc_lr_overall.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ligand", "receptor", "n_samples", "frac_samples", "mean_prob"])
        rows = sorted(lr_over.items(), key=lambda kv: len(kv[1]["samp"]), reverse=True)
        for (lg, rc), v in rows:
            if len(v["samp"]) < MIN_SAMPLES:
                continue
            w.writerow([lg, rc, len(v["samp"]), f"{len(v['samp'])/n_aml:.4f}", f"{st.mean(v['prob']):.4f}"])

    # SCS axis: primitive sender -> immune receiver, per specific L-R
    with open(CCC / "ccc_lr_scs_axis.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sender_bin", "receiver_bin", "ligand", "receptor", "n_samples", "frac_samples", "mean_prob"])
        scs = {k: v for k, v in lr_full.items() if k[0] in PRIMITIVE and k[1] in IMMUNE}
        rows = sorted(scs.items(), key=lambda kv: len(kv[1]["samp"]), reverse=True)
        for (s, r, lg, rc), v in rows:
            if len(v["samp"]) < MIN_SAMPLES:
                continue
            w.writerow([s, r, lg, rc, len(v["samp"]), f"{len(v['samp'])/n_aml:.4f}", f"{st.mean(v['prob']):.4f}"])

    print(f"[ccc-agg] AML samples: edge-table {aml_edge_samples}, tensors {n_aml}")
    print(f"[ccc-agg] wrote ccc_bin_strength_cohort.csv, ccc_lr_overall.csv, ccc_lr_scs_axis.csv -> {CCC}")


if __name__ == "__main__":
    main()
