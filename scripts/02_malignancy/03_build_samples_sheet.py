#!/usr/bin/env python3
# 03_build_samples_sheet.py ----
# Build a one-line-per-sample sheet  <sample_id>\t<srr1,srr2,...>  from an SRA Run Selector
# table (SraRunTable.csv) or a pysradb --detailed tsv.
#
# UNIVERSAL grouping key = GSM (one GEO sample = one library = one biological sample). This
# correctly handles BOTH layouts:
#   - GSE227903: each GSM has ~16 SRR lanes  -> grouped into one sample
#   - GSE239721: each GSM has 1 SRR          -> one sample per run
#
# sample_id naming (pick what matches your per-sample QC RDS basenames):
#   --name_from geo_title : fetch GEO series-matrix titles (e.g. PT9A) and use them  [needs --gse]
#   --name_from donor_tp  : parse "Donor 1216, Diagnosis" -> 1216_Dg   (GSE227903 style)
#   --name_from gsm       : use the GSM accession
#   default = auto: geo_title if --gse resolves, else donor_tp if parseable, else gsm
# ALWAYS pass --rds_dir to cross-check/rename col1 against the real RDS basenames.
#
# 3a-1 rewire: rename only (pure Python helper; no config-source system, no path defaults in code).
#
# Usage:
#   python scripts/02_malignancy/03_build_samples_sheet.py SraRunTable_GSE239721.csv --gse GSE239721 \
#       --out samples_GSE239721.tsv --rds_dir /LARGE1/.../01_per_sample_qc/GSE239721

import argparse, csv, gzip, io, os, re, sys, urllib.request
from collections import defaultdict, OrderedDict

TP_SHORT = {"diagnosis":"Dg","dx":"Dg","newly diagnosed":"Dg","relapse":"R",
            "minimal residual disease":"MRD","mrd":"MRD","residual":"MRD"}

def sniff_delim(path):
    with open(path) as f: head = f.readline()
    return "\t" if head.count("\t") >= head.count(",") else ","

def pick(cols, *cands):
    low = {c.lower().strip(): c for c in cols}
    for cand in cands:
        if cand in low: return low[cand]
    return None

def sanitize(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", (s or "").strip()).strip("_")

def tp_short(s):
    s = (s or "").strip().lower()
    if s in TP_SHORT: return TP_SHORT[s]
    for k, v in TP_SHORT.items():
        if k in s: return v
    return sanitize(s)[:6] or "NA"

def donor_num(s):
    m = re.search(r"(\d+)", s or ""); return m.group(1) if m else (s or "NA").strip()

def parse_donor_tp(*texts):
    for t in texts:
        m = re.search(r"Donor\s+(\w+)\s*,\s*([^;]+)", t or "")
        if m: return f"{donor_num(m.group(1))}_{tp_short(m.group(2))}"
    return None

def find_gsm_col(cols, rows):
    for c in (pick(cols,"library name"), pick(cols,"sample name"), pick(cols,"gsm")):
        if c and any(re.match(r"GSM\d+", str(r.get(c,"")) or "") for r in rows[:10]):
            return c
    return None

def fetch_geo_titles(gse):
    num = re.sub(r"\D", "", gse)
    stub = "GSE" + (num[:-3] if len(num) > 3 else "") + "nnn"
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        raw = urllib.request.urlopen(req, timeout=90).read()
        text = gzip.decompress(raw).decode("utf-8", "replace")
    except Exception as e:
        print(f"[geo] could not fetch series matrix ({e}); GEO titles unavailable",
              file=sys.stderr)
        return {}
    acc = title = None
    relations = []   # each !Sample_relation row -> per-sample values (BioSample / SRA links)
    for line in text.splitlines():
        if line.startswith("!Sample_geo_accession"):
            acc = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            title = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_relation"):
            relations.append([x.strip().strip('"') for x in line.split("\t")[1:]])
    titles = dict(zip(acc, title)) if acc and title else {}
    # link SRX -> GSM via the SRA relation ("SRA: ...term=SRXNNN") so we can group pysradb runs
    srx2gsm = {}
    if acc:
        for i, gsm in enumerate(acc):
            for rel in relations:
                if i < len(rel):
                    m = re.search(r"SRX\d+", rel[i])
                    if m:
                        srx2gsm[m.group(0)] = gsm
    return titles, srx2gsm

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("metadata")
    ap.add_argument("--gse", default=None, help="GEO series id (enables geo_title naming)")
    ap.add_argument("--name_from", choices=["auto","geo_title","donor_tp","gsm"], default="auto")
    ap.add_argument("--out", default=None)
    ap.add_argument("--rds_dir", default=None)
    args = ap.parse_args()

    delim = sniff_delim(args.metadata)
    with open(args.metadata, newline="") as f:
        rows = list(csv.DictReader(f, delimiter=delim))
    if not rows: sys.exit("ERROR: no rows parsed.")
    cols = rows[0].keys()

    c_run = pick(cols, "run", "run_accession")
    if not c_run: sys.exit("ERROR: no run/run_accession column.")
    # Grouping key = one biological sample. Prefer GSM (SRA Run Selector); fall back to
    # experiment_accession/SRX (pysradb --detailed), since 1 SRX = 1 library = 1 GSM.
    c_gsm = find_gsm_col(cols, rows)
    c_srx = pick(cols, "experiment_accession", "experiment")
    group_col = c_gsm or c_srx
    if not group_col:
        sys.exit("ERROR: no GSM column (Library Name/Sample Name) and no experiment_accession "
                 "column found. Provide an SRA Run Selector table or a pysradb --detailed table.")
    by_gsm_mode = bool(c_gsm)
    print(f"[group] grouping runs by {'GSM' if by_gsm_mode else 'experiment_accession (SRX)'} "
          f"column '{group_col}'", file=sys.stderr)
    c_title = pick(cols, "sample_title", "experiment_title", "experiment_alias")
    c_subj  = pick(cols, "subject_id", "subject id")
    c_tp     = pick(cols, "time_point", "time point", "treatment", "disease_state")

    geo, srx2gsm = ({}, {})
    if args.gse and args.name_from in ("auto", "geo_title"):
        geo, srx2gsm = fetch_geo_titles(args.gse)
    if geo:
        print(f"[geo] fetched {len(geo)} GSM->title, {len(srx2gsm)} SRX->GSM mappings", file=sys.stderr)

    # group runs by the chosen key (preserve first-seen order)
    by_key = OrderedDict()
    for r in rows:
        srr = (r.get(c_run) or "").strip()
        key = (r.get(group_col) or "").strip()
        if not srr.startswith("SRR"): continue
        if by_gsm_mode and not key.startswith("GSM"): continue
        if (not by_gsm_mode) and not key.startswith("SRX"): continue
        by_key.setdefault(key, {"srr": [], "row": r})
        by_key[key]["srr"].append(srr)

    def gsm_of(key):
        return key if by_gsm_mode else srx2gsm.get(key, "")

    # choose sample_id per group (geo_title preferred; resolves via SRX->GSM when grouped by SRX)
    def name_for(key, r):
        gsm = gsm_of(key)
        if args.name_from in ("auto","geo_title") and gsm and geo.get(gsm):
            return sanitize(geo[gsm])
        if args.name_from in ("auto","donor_tp"):
            dt = parse_donor_tp(r.get(c_title,""), r.get(c_subj,""))
            if dt: return dt
            if c_subj and c_tp and r.get(c_subj):
                return f"{donor_num(r.get(c_subj))}_{tp_short(r.get(c_tp,''))}"
        t = (r.get(c_title,"") or "").strip()   # last resort: sanitized SRA title, else the key
        return sanitize(t) if t else key

    samples = OrderedDict()   # sample_id -> {srr,key,gsm,title}
    for key, d in by_key.items():
        sid = name_for(key, d["row"])
        if sid in samples and samples[sid]["key"] != key:   # collision -> suffix with key tail
            sid = f"{sid}_{key[-4:]}"
        samples[sid] = {"srr": sorted(set(d["srr"])), "key": key, "gsm": gsm_of(key),
                        "title": geo.get(gsm_of(key), "")}

    out = args.out or (f"samples_{args.gse}.tsv" if args.gse else "samples_built.tsv")
    with open(out, "w") as f:
        f.write("# sample_id\tsrr_list  (one line per biological sample = one GSM)\n")
        f.write("# VERIFY sample_id matches the QC RDS basenames; rename col1 if needed.\n")
        for sid, d in samples.items():
            f.write(f"{sid}\t{','.join(d['srr'])}\n")
    # sidecar map for traceability
    mp = os.path.splitext(out)[0] + ".map.tsv"
    with open(mp, "w") as f:
        f.write("sample_id\tkey\tgsm\tgeo_title\tn_runs\n")
        for sid, d in samples.items():
            f.write(f"{sid}\t{d['key']}\t{d['gsm']}\t{d['title']}\t{len(d['srr'])}\n")

    print(f"[built] {out}: {len(samples)} samples, "
          f"{sum(len(d['srr']) for d in samples.values())} runs (map -> {mp})", file=sys.stderr)
    print(f"{'sample_id':<22}{'n_runs':<8}{'key':<14}{'geo_title'}", file=sys.stderr)
    for sid, d in samples.items():
        print(f"{sid:<22}{len(d['srr']):<8}{d['key']:<14}{d['title']}", file=sys.stderr)

    if args.rds_dir and os.path.isdir(args.rds_dir):
        rds = [os.path.splitext(x)[0] for x in os.listdir(args.rds_dir)
               if x.lower().endswith(".rds")]
        print("\n[rds check] sample_id -> exact RDS match?", file=sys.stderr)
        for sid in samples:
            exact = sid in rds
            near = [b for b in rds if sid in b or b in sid] if not exact else []
            tag = "OK exact" if exact else (f"NEAR {near}" if near else "NO MATCH -> rename")
            print(f"  {sid:<18} {tag}", file=sys.stderr)
        miss = [b for b in rds if b not in samples]
        if miss:
            print(f"\n[rds check] RDS basenames with no sample_id yet: {sorted(miss)}",
                  file=sys.stderr)

if __name__ == "__main__":
    main()