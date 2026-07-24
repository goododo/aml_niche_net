#!/usr/bin/env python
# e10_run_cnmf.py ----
# Phase 3 step 1 (Python): per-sample cNMF over K=4..9, keeping ALL K programs for the aggregation
# pool (no single-K selection -- cross-sample SRRS is the real filter). Reads the e05 mtx triple,
# lets cNMF pick HVGs (--numhvgenes), factorizes with N_ITER restarts, and for each K writes the
# consensus gene-spectra (program x gene loadings) + the top-N gene set per program.
#
# cNMF (Kotliar 2019) v1.7.1 pipeline: prepare -> factorize -> combine -> consensus(k) -> load spectra.
#
# INPUT  : CNMF_INPUT_DIR/<ds>/<sample>/{counts.mtx, barcodes.tsv, genes.tsv, meta.json}
# OUTPUT : CNMF_PROG_DIR/<ds>/<sample>__programs.csv   (long: dataset,sample,K,program,rank,gene,loading)
#          + CNMF_PROG_DIR/<ds>/<sample>__programs_meta.json
# cNMF working dir (big, scratch): CNMF_RUN_DIR/<ds>/<sample>/
# Per-sample, serial, resume-skip.
# Usage:
#   python e10_run_cnmf.py --dataset GSE227903
#   python e10_run_cnmf.py --dataset GSE239721 --sample PT10A
#
# NOTE: constants mirror d00_config.R (single source of truth is R; kept in sync here by hand):
#   K_RANGE=4..9, N_ITER=50, NUMHVG=2000, TOP_N=50, SEED=20260605.

import os, sys, json, argparse, logging
import numpy as np
import pandas as pd
import scipy.io as sio
import scanpy as sc
from cnmf import cNMF

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("e10")

# -- paths / constants (keep in sync with d00_config.R) ----
LARGE1 = "/LARGE1/gr10634/gaozy/aml_niche_net"
FAST   = "/FAST/gr10634/gaozy/aml_niche_net"
CNMF_INPUT_DIR = os.path.join(LARGE1, "06_cnmf", "input")
CNMF_RUN_DIR   = os.path.join(LARGE1, "06_cnmf", "runs")
CNMF_PROG_DIR  = os.path.join(FAST, "results/tables/05_cnmf/programs")
K_RANGE  = list(range(4, 10))   # 4..9
N_ITER   = 50
NUMHVG   = 2000
TOP_N    = 50
SEED     = 20260605
DENSITY_THRESHOLD = 2.00        # consensus default (keep all replicates; no aggressive outlier filter)


def load_mtx(sample_dir):
    """Read e05 mtx triple (genes x cells) into an AnnData (cells x genes), raw counts."""
    m = sio.mmread(os.path.join(sample_dir, "counts.mtx")).tocsr()   # genes x cells
    genes = [l.strip() for l in open(os.path.join(sample_dir, "genes.tsv"))]
    cells = [l.strip() for l in open(os.path.join(sample_dir, "barcodes.tsv"))]
    adata = sc.AnnData(X=m.T.tocsr())          # cells x genes
    adata.var_names = genes
    adata.obs_names = cells
    adata.var_names_make_unique()
    return adata


def run_one(ds, sample):
    in_dir  = os.path.join(CNMF_INPUT_DIR, ds, sample)
    meta_in = os.path.join(in_dir, "meta.json")
    if not os.path.exists(meta_in):
        log.info("    [no meta] %s/%s -- run e05 first", ds, sample); return
    meta = json.load(open(meta_in))
    if meta.get("input_mode") == "skip_too_few":
        log.info("    [skip] %s  (e05 marked skip_too_few)", sample); return

    prog_csv = os.path.join(CNMF_PROG_DIR, ds, f"{sample}__programs.csv")
    if os.path.exists(prog_csv):
        log.info("    [skip done] %s", sample); return

    adata = load_mtx(in_dir)
    n_cells, n_genes = adata.shape
    # numhvgenes can't exceed available genes
    numhvg = min(NUMHVG, n_genes)
    log.info("    [%s] cells=%d genes=%d mode=%s -> cNMF K=%s iter=%d hvg=%d",
             sample, n_cells, n_genes, meta.get("input_mode"), K_RANGE, N_ITER, numhvg)

    run_dir = os.path.join(CNMF_RUN_DIR, ds)
    os.makedirs(run_dir, exist_ok=True)
    name = sample

    # cNMF wants an .h5ad path as counts input
    h5ad_path = os.path.join(in_dir, "counts.h5ad")
    adata.write(h5ad_path)

    cnmf_obj = cNMF(output_dir=run_dir, name=name)
    cnmf_obj.prepare(counts_fn=h5ad_path, components=K_RANGE, n_iter=N_ITER,
                     num_highvar_genes=numhvg, seed=SEED)
    cnmf_obj.factorize(worker_i=0, total_workers=1)
    cnmf_obj.combine()
    # k_selection_plot is optional/diagnostic; skip to stay headless

    rows = []
    cnmf_sample_dir = os.path.join(run_dir, name)
    for k in K_RANGE:
        import glob
        existing = glob.glob(os.path.join(cnmf_sample_dir, f"{name}.gene_spectra_score.k_{k}.dt_*.txt"))
        if not existing:
            try:
                cnmf_obj.consensus(k=k, density_threshold=DENSITY_THRESHOLD, show_clustering=False)
            except Exception as e:
                log.info("      [K=%d consensus failed] %s: %s", k, sample, str(e)[:120]); continue

        # Read the gene-spectra score file cNMF wrote, instead of load_results() which in 1.7.1
        # chokes on Ensembl-style gene names (int('AL627309.1')). File = programs x genes z-scores.
        g = glob.glob(os.path.join(cnmf_sample_dir, f"{name}.gene_spectra_score.k_{k}.dt_*.txt"))
        score_fn = g[0] if g else None
        if score_fn is None:
            log.info("      [K=%d no score file] %s", k, sample); continue

        sp = pd.read_csv(score_fn, sep="\t", index_col=0)          # rows=program(1..k), cols=gene
        for prog in sp.index:
            ranked = sp.loc[prog].sort_values(ascending=False)
            topg = ranked.head(TOP_N)
            for rank, (gene, loading) in enumerate(topg.items(), start=1):
                rows.append((ds, sample, int(k), int(prog), rank, str(gene), float(loading)))

    if not rows:
        log.info("    [no programs] %s -- all K failed", sample); return

    out = pd.DataFrame(rows, columns=["dataset", "sample", "K", "program", "rank", "gene", "loading"])
    os.makedirs(os.path.dirname(prog_csv), exist_ok=True)
    out.to_csv(prog_csv, index=False)
    json.dump({"dataset": ds, "sample": sample, "input_mode": meta.get("input_mode"),
               "n_cells": int(n_cells), "n_genes": int(n_genes), "numhvg": int(numhvg),
               "K_range": K_RANGE, "n_iter": N_ITER, "top_n": TOP_N, "seed": SEED,
               "n_programs": int(out[["K", "program"]].drop_duplicates().shape[0])},
              open(prog_csv.replace("__programs.csv", "__programs_meta.json"), "w"), indent=2)
    log.info("    [write] %s  programs=%d -> %s",
             sample, out[["K", "program"]].drop_duplicates().shape[0], prog_csv)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample", default=None)
    args = ap.parse_args()

    ds_dir = os.path.join(CNMF_INPUT_DIR, args.dataset)
    if not os.path.isdir(ds_dir):
        sys.exit(f"[e10] no e05 input dir for {args.dataset} (run e05 first)")
    samples = [args.sample] if args.sample else sorted(
        d for d in os.listdir(ds_dir) if os.path.isdir(os.path.join(ds_dir, d)))
    log.info("[1] %d sample(s) for %s", len(samples), args.dataset)
    for s in samples:
        run_one(args.dataset, s)
    log.info("[done] e10 for %s -> %s", args.dataset, os.path.join(CNMF_PROG_DIR, args.dataset))


if __name__ == "__main__":
    main()