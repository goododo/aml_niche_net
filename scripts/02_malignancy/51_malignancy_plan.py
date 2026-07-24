#!/usr/bin/env python3
# 51_malignancy_plan.py ----
# STATUS SCANNER (not a data-availability oracle): for every AML sample, report which malignancy
# arms are actually computed, which evidence TIER it can reach, and whether a >=2-type consensus
# is possible NOW. Does NOT query GEO/SRA -- the "is FASTQ usable for Numbat" decision is a
# per-dataset human call recorded in dataset_roles.tsv (column 'fastq').
#
# Emits: <prefix>_malignancy_plan.tsv, <prefix>_pilot_roster.txt, _generated/run_consensus.sh
#
# 3a rewire (no logic change): renamed; the generated command chain now targets the new script
#   numbers (41/43/42/50); DEF_MAL points at 02_malignancy; the generated runner is written to
#   _generated/run_consensus.sh (gitignored). Tier strings are unchanged (match 50; the tier
#   vocabulary is unified with 50 in batch 3b).
#
# Evidence TYPES (match 50): expression {inferCNV, copykat, scevan, author}, allele {Numbat},
# SNV {VarTrix}. Tier A if allele/snv; B if >=2 expression methods; else C.
#
# Usage:
#   python 51_malignancy_plan.py --roles dataset_roles.tsv \
#     --scripts_dir /FAST/.../scripts/02_malignancy \
#     --infercnv_root /LARGE1/gr10634/gaozy/aml_niche_net/05_cnv_snv/infercnv \
#     [--cnv_root /LARGE1/.../03_cnv_snv] [--qc_root ...] [--mal_tab ...] \
#     [--with_copykat_fallback] [--min_tools 1]

import argparse, csv, os, glob, sys

DEF_QC   = "/LARGE1/gr10634/gaozy/aml_niche_net/02_seurat_objects/01_per_sample_qc"
DEF_CNV  = "/LARGE1/gr10634/gaozy/aml_niche_net/03_cnv_snv"          # numbat / vartrix / author / copykat
DEF_ICNV = "/LARGE1/gr10634/gaozy/aml_niche_net/05_cnv_snv/infercnv" # inferCNV lives here (separate root)
DEF_MAL  = "/FAST/gr10634/gaozy/aml_niche_net/results/tables/02_malignancy"
SKIP_ROLES = {"healthy","reference","sorted","comparator","receiver_only","skip"}

def load_roles(path):
    roles = {}
    if not path or not os.path.exists(path):
        print(f"[warn] roles '{path}' missing; defaults ...", file=sys.stderr); return roles
    with open(path) as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    fields = ["dataset","role","author_cnv_col","has_snv_vcf","fastq"]
    for r in csv.DictReader(lines, fieldnames=fields, delimiter="\t"):
        ds = (r.get("dataset") or "").strip()
        if not ds: continue
        roles[ds] = {
            "role": (r.get("role") or "aml_discovery").strip(),
            "author_cnv_col": (r.get("author_cnv_col") or "").strip(),
            "has_snv_vcf": (r.get("has_snv_vcf") or "").strip().lower() in ("1","y","yes","true"),
            "fastq": (r.get("fastq") or "no").strip().lower() in ("1","y","yes","true"),
        }
    return roles
    
def find_infercnv_obj(ic_dir):
    if not os.path.isdir(ic_dir): return None
    direct = os.path.join(ic_dir, "run.final.infercnv_obj")
    if os.path.exists(direct): return direct
    for pat in ("**/run.final.infercnv_obj", "**/*.infercnv_obj"):
        hits = glob.glob(os.path.join(ic_dir, pat), recursive=True)
        if hits: return hits[0]
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roles", default="dataset_roles.tsv")
    ap.add_argument("--scripts_dir", required=True)
    ap.add_argument("--qc_root",  default=DEF_QC)
    ap.add_argument("--cnv_root", default=DEF_CNV)
    ap.add_argument("--infercnv_root", default=DEF_ICNV)
    ap.add_argument("--infercnv_burden_root",
                    default="/LARGE1/gr10634/gaozy/aml_niche_net/05_cnv_snv/infercnv_burden")
    ap.add_argument("--mal_tab",  default=DEF_MAL)
    ap.add_argument("--out_prefix", default="51")
    ap.add_argument("--min_tools", type=int, default=1)
    ap.add_argument("--with_copykat_fallback", action="store_true",
                    help="for inferCNV-only samples, actively run c44 CopyKAT to reach Tier B")
    args = ap.parse_args()
    roles = load_roles(args.roles)
    SD = args.scripts_dir.rstrip("/")

    rows, ready, skipped = [], [], []
    datasets = sorted(d for d in os.listdir(args.qc_root)
                      if os.path.isdir(os.path.join(args.qc_root, d)))
    for ds in datasets:
        cfg = roles.get(ds, {"role":"aml_discovery","author_cnv_col":"","has_snv_vcf":False,"fastq":False})
        if cfg["role"] in SKIP_ROLES:
            skipped.append(f"{ds}:{cfg['role']}"); continue
        qc_dir = os.path.join(args.qc_root, ds)
        for s in sorted(f[:-4] for f in os.listdir(qc_dir) if f.lower().endswith(".rds")):
            qc = os.path.join(qc_dir, f"{s}.rds")
            nb_pc = os.path.join(args.cnv_root,"numbat",ds,s,"numbat",f"{s}__numbat_percell.csv")
            ic_dir = os.path.join(args.infercnv_root, ds, s)
            ic_obj = find_infercnv_obj(ic_dir)
            ic_pc  = os.path.join(ic_dir, f"{s}__infercnv_percell.csv")
            ic_burden = os.path.join(args.infercnv_burden_root, ds, f"{s}_infercnv_burden.csv")
            vt_pc  = os.path.join(args.cnv_root,"vartrix",ds,s,f"{s}__vartrix_percell.csv")
            au_pc  = os.path.join(args.cnv_root,"author", ds,s,f"{s}__author_percell.csv")
            ck_pc  = os.path.join(args.cnv_root,"copykat",ds,s,f"{s}__copykat_percell.csv")
            sc_pc  = os.path.join(args.cnv_root,"scevan", ds,s,f"{s}__scevan_percell.csv")

            allele = os.path.exists(nb_pc)
            snv    = cfg["has_snv_vcf"] and os.path.exists(vt_pc)
            expr = []
            if ic_obj or os.path.exists(ic_pc) or os.path.exists(ic_burden): expr.append("infercnv")
            if cfg["author_cnv_col"]:           expr.append("author")
            if os.path.exists(ck_pc):           expr.append("copykat")
            if os.path.exists(sc_pc):           expr.append("scevan")

            # numbat status for the plan (no scraping; uses your 'fastq' determination)
            nb_status = "done" if allele else ("pending" if cfg["fastq"] else "na")

            # optional CopyKAT fallback to lift inferCNV-only -> Tier B
            plan_ck = (args.with_copykat_fallback and not (allele or snv)
                       and "author" not in expr and "copykat" not in expr and len(expr) < 2)
            if plan_ck: expr.append("copykat")

            types_now = (1 if expr else 0) + (1 if allele else 0) + (1 if snv else 0)
            # 51 predicts the evidence CONFIGURATION reachable now. The final A_concordant vs
            # B_multi_partial split is data-dependent (cross-type conflict) and decided by 50;
            # here we only distinguish single-type (C) from multi-type (A or B, TBD at 50).
            tier = "C_single" if types_now < 2 else "AB_multi"
            arms = (["numbat"] if allele else []) + (["vartrix"] if snv else []) + expr

            # build command chain
            pre, flags = [], [f"--qc_rds {qc}"]
            if allele: flags.append(f"--numbat {nb_pc}")
            if snv:    flags.append(f"--vartrix {vt_pc}")
            if "infercnv" in expr:
                if not os.path.exists(ic_pc):
                    if os.path.exists(ic_burden):
                        pre.append(f"Rscript {SD}/41_infercnv_to_percell.R --sample {s} "
                                   f"--burden_csv {ic_burden} --outdir {ic_dir}")
                    elif ic_obj:
                        pre.append(f"Rscript {SD}/41_infercnv_to_percell.R --sample {s} "
                                   f"--infercnv_dir {ic_dir} --outdir {ic_dir}")
                flags.append(f"--infercnv {ic_pc}")
            if "author" in expr:
                if not os.path.exists(au_pc):
                    pre.append(f"Rscript {SD}/43_author_to_percell.R --sample {s} --qc_rds {qc} "
                               f"--outdir {os.path.dirname(au_pc)} --col {cfg['author_cnv_col']}")
                flags.append(f"--author {au_pc}")
            if plan_ck:
                pre.append(f"Rscript {SD}/42_exprcnv_run.R --sample {s} --method copykat "
                           f"--qc_rds {qc} --outdir {os.path.dirname(ck_pc)}")
                flags.append(f"--copykat {ck_pc}")
            elif "copykat" in expr:
                flags.append(f"--copykat {ck_pc}")
            if "scevan" in expr:
                flags.append(f"--scevan {sc_pc}")
            c50 = (f"Rscript {SD}/50_consensus_malignancy.R --sample {s} --dataset {ds} "
                   f"{' '.join(flags)} --outdir {os.path.join(args.mal_tab, ds)} "
                   f"--min_tools {args.min_tools}")

            ready_now = "YES" if types_now >= 2 else ("pilot1" if types_now == 1 else "no")
            rows.append([ds, s, cfg["role"], cfg["fastq"], nb_status,
                         ("yes" if ("infercnv" in expr) else "no"),
                         (cfg["author_cnv_col"] or "-"),
                         "+".join(arms) if arms else "(none)", types_now, tier, ready_now])
            if types_now >= 2 or (args.min_tools <= 1 and types_now >= 1):
                ready.append((ds, s, tier, pre, c50))

    cols = ["dataset","sample","role","fastq","numbat","infercnv","author_col",
            "arms_available","n_types_now","expected_tier","ready_now"]
    plan = f"{args.out_prefix}_malignancy_plan.tsv"
    with open(plan, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows: f.write("\t".join(str(x) for x in r) + "\n")

    roster = f"{args.out_prefix}_pilot_roster.txt"
    with open(roster, "w") as f:
        f.write("# samples with >=2 evidence types ready NOW (true >=2/3 consensus)\n")
        for ds, s, tier, _, _ in ready:
            if tier != "C_single": f.write(f"{ds}\t{s}\t{tier}\n")

    os.makedirs("_generated", exist_ok=True)
    sh = "_generated/run_consensus.sh"
    with open(sh, "w") as f:
        f.write("#!/usr/bin/env bash\n"
                "# Auto-generated by 51_malignancy_plan.py: 41/43[/42] prereqs then 50.\n"
                "# Run from the PROJECT ROOT so each Rscript's here::here() resolves the config.\n"
                "set -e\n")
        for ds, s, tier, pre, c50 in ready:
            f.write(f"\necho '=== {ds}/{s} (tier {tier}) ==='\n")
            for p in pre: f.write(p + "\n")
            f.write(c50 + "\n")
    os.chmod(sh, 0o755)

    n2 = sum(1 for r in rows if r[-1] == "YES"); n1 = sum(1 for r in rows if r[-1] == "pilot1")
    print(f"[done] {len(rows)} AML samples; ready NOW(>=2 types)={n2}, single-type(pilot)={n1}",
          file=sys.stderr)
    print(f"   plan -> {plan}\n   roster -> {roster}\n   script -> {sh}", file=sys.stderr)
    if skipped: print("   skipped: " + ", ".join(skipped), file=sys.stderr)

if __name__ == "__main__":
    main()
