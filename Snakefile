# Snakefile ----
# Step 1 of the five-step Snakemake migration described in CLAUDE.md.
#
#   1. One rule, one sample, hardcoded paths        <-- YOU ARE HERE
#   2. Add {sample} wildcard; read the list from infercnv_tasks.tsv
#   3. Add the upstream rule so two stages chain
#   4. Move parameters into config.yaml
#   5. Per-rule environment declarations
#
# WHY THIS STAGE FIRST. 41_infercnv_to_percell.R is per-sample, takes ~7 s, and
# its input is the burden CSV that 44_infercnv_run_one.R produces. So step 3
# ("chain two stages") is literally adding the 44 rule above this one -- the
# chain is the real pipeline, not a demo. And a 7 s rule can be run twenty times
# while learning; an 8 h inferCNV rule cannot.
#
# WHAT SNAKEMAKE BUYS US, concretely. This repo currently decides "is my output
# stale?" with hand-written .ins lists. On 2026-08-25 one of those lists omitted
# the CCC tensors, so the freshness check printed "[skip] is current" over an
# answer computed from 78 of 138 samples. Snakemake derives the same decision
# from the input: field below, which cannot be silently incomplete in the same
# way -- if a file is not listed, the rule cannot read it either.
#
# RUN IT (SM is the standalone env; the pipeline's own conda env is separate):
#   SM=/FAST/gr10634/gaozy/snakemake_env/bin/snakemake
#   $SM -n            # dry run: say what you WOULD do, touch nothing
#   $SM --cores 1     # actually do it
#   $SM -n            # again: should now say "Nothing to be done"

# -- hardcoded for step 1. Step 2 turns SAMPLE into a {sample} wildcard. --
DATASET = "GSE227903"
SAMPLE  = "1216_Dg"
LARGE1  = "/LARGE1/gr10634/gaozy/aml_niche_net"
ENV     = "/FAST/gr10634/gaozy/general_env"

BURDEN_CSV = f"{LARGE1}/05_cnv_snv/infercnv_burden/{DATASET}/{SAMPLE}_infercnv_burden.csv"
OUT_DIR    = f"{LARGE1}/05_cnv_snv/infercnv/{DATASET}/{SAMPLE}"
SCRIPT     = "scripts/02_malignancy/41_infercnv_to_percell.R"

# Snakemake runs the FIRST rule in the file when you name no target, so this one
# is the default. A rule fires when its output is missing, or when any input is
# NEWER than the output. That is the whole mechanism -- there is nothing else.
rule infercnv_percell:
    input:
        burden = BURDEN_CSV,
        # The script is listed as an input ON PURPOSE. Edit 41 and this rule
        # re-runs; that is the "timestamps always line up" property. Leave it
        # out and Snakemake will happily keep an output produced by code that no
        # longer exists.
        script = SCRIPT,
        # Same reasoning for the config files 41 sources: INFERCNV_SCORE_Q lives
        # in config_malignancy.R, so changing the threshold must invalidate the
        # output.
        cfg_paths = "scripts/config/config_paths.R",
        cfg_mal   = "scripts/config/config_malignancy.R",
    output:
        percell = f"{OUT_DIR}/{SAMPLE}__infercnv_percell.csv",
    log:
        f"logs/snakemake/infercnv_percell.{DATASET}.{SAMPLE}.log",
    shell:
        # conda run -p, matching how 45_infercnv_submit.slurm invokes R today.
        # Step 5 replaces this with a per-rule `conda:` directive.
        "mkdir -p $(dirname {log}) && "
        "conda run -p {ENV} Rscript {input.script} "
        "  --sample {SAMPLE} "
        "  --burden_csv {input.burden} "
        "  --outdir {OUT_DIR} "
        "&> {log}"
