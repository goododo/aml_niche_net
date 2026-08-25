# Snakefile ----
# Step 2 of the five-step Snakemake migration described in CLAUDE.md.
#
#   1. One rule, one sample, hardcoded paths            done (commit 68f8423)
#   2. Add the {sample} wildcard; read the sample list  <-- YOU ARE HERE
#   3. Add the upstream rule so two stages chain
#   4. Move parameters into config.yaml
#   5. Per-rule environment declarations
#
# DEVIATION FROM THE PLAN TEXT, AND WHY. The plan says step 2 reads
# infercnv_tasks.tsv. That file holds 160 of the cohort's 212 samples. The other
# 52 live in five sibling files (_healthy, _retry, _sorted, _sens_eval,
# _matchedref_eval) whose union happens to equal the cohort today, with nothing
# enforcing it. Encoding the 160 here would start the migration by silently
# dropping four datasets (GSE185991, GSE253355, E-MTAB-11536, GSE147989).
#
# So the queue is read from ref_norm_summary.csv -- a GENERATED pipeline output,
# 212 rows, and the same file 45_infercnv_submit.slurm already tells you to
# build the task list from. The hand-maintained files become the CHECK rather
# than the source: the assertion below fails the workflow if the two ever
# disagree, which is the evidence needed before those six files can be deleted.
#
# RUN IT:
#   cd /FAST/gr10634/gaozy/aml_niche_net
#   SM=/FAST/gr10634/gaozy/snakemake_env/bin/snakemake
#   $SM -n
#   $SM --cores 4

import csv
import glob

LARGE1 = "/LARGE1/gr10634/gaozy/aml_niche_net"
ENV    = "/FAST/gr10634/gaozy/general_env"

REFNORM = "results/tables/02_malignancy/ref_norm_summary.csv"
SCRIPT  = "scripts/02_malignancy/41_infercnv_to_percell.R"

# -- the cohort, from the generated summary --
DATASETS, SAMPLES = [], []
with open(REFNORM) as fh:
    for row in csv.DictReader(fh):
        DATASETS.append(row["dataset"])
        SAMPLES.append(row["sample_id"])

# -- the check: does the generated cohort still equal the hand-maintained queue? --
# Not decoration. It is the only thing standing between "Snakemake covers the
# cohort" and "Snakemake covers whatever happened to be in one text file".
_generated = set(zip(DATASETS, SAMPLES))
_queued = set()
for f in glob.glob("infercnv_tasks*.tsv"):
    for line in open(f):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2 and parts[0] and parts[1]:
            _queued.add((parts[0], parts[1]))
if _generated != _queued:
    raise WorkflowError(
        "ref_norm_summary.csv and the infercnv_tasks*.tsv queue disagree.\n"
        f"  only in ref_norm_summary.csv : {sorted(_generated - _queued)}\n"
        f"  only in infercnv_tasks*.tsv  : {sorted(_queued - _generated)}\n"
        "Do not paper over this by editing the tsv. One of them is wrong about "
        "which samples are in the cohort, and that question has to be answered first."
    )

# Sample ids contain underscores (GSM8019220_H36), and the output pattern uses a
# double underscore. Constrain both wildcards to a single path component so the
# regex cannot split a name at the wrong place.
wildcard_constraints:
    dataset = r"[^/]+",
    sample  = r"[^/]+",

# HOW STALENESS IS ACTUALLY DECIDED, measured 2026-08-25 on this workflow.
# Snakemake 9 records a sha256 of every input in .snakemake/metadata, so it does
# NOT work on timestamps once a job has run through it:
#
#   touch a burden CSV (mtime moves, bytes identical) ->   0 jobs
#   touch a config     (mtime moves, bytes identical) ->   0 jobs
#   change one burden CSV's CONTENT                   ->   1 job  (that sample)
#   change config_malignancy.R's CONTENT              -> 212 jobs (all samples)
#
# That is stronger than "update the timestamps and hope they line up". A stray
# touch -- rsync, git checkout, a restored backup -- cannot cause a spurious
# 212-sample rebuild, and a real edit cannot be missed. Before the first run
# there is no metadata, so it falls back to mtime; that is the only regime where
# a bare touch does anything.

# rule all exists only to name the targets. It has no shell command; its input
# list IS the request. Snakemake works backwards from here.
rule all:
    input:
        expand(
            f"{LARGE1}/05_cnv_snv/infercnv/{{dataset}}/{{sample}}/{{sample}}__infercnv_percell.csv",
            zip, dataset=DATASETS, sample=SAMPLES,
        ),

rule infercnv_percell:
    input:
        burden    = f"{LARGE1}/05_cnv_snv/infercnv_burden/{{dataset}}/{{sample}}_infercnv_burden.csv",
        # The script and the configs it sources are inputs on purpose: edit the
        # code or the threshold and every affected output is rebuilt. This is
        # the property a hand-written .ins list keeps losing.
        script    = SCRIPT,
        cfg_paths = "scripts/config/config_paths.R",
        cfg_mal   = "scripts/config/config_malignancy.R",
    output:
        percell = f"{LARGE1}/05_cnv_snv/infercnv/{{dataset}}/{{sample}}/{{sample}}__infercnv_percell.csv",
    params:
        outdir = f"{LARGE1}/05_cnv_snv/infercnv/{{dataset}}/{{sample}}",
    log:
        "logs/snakemake/infercnv_percell.{dataset}.{sample}.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "conda run -p " + ENV + " Rscript {input.script} "
        "  --sample {wildcards.sample} "
        "  --burden_csv {input.burden} "
        "  --outdir {params.outdir} "
        "&> {log}"
