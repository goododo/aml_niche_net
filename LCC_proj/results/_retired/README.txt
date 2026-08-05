RETIRED OUTPUTS -- LCC_proj, 2026-08-04
=======================================

Everything in this directory was moved out of results/ and logs/ because it is either
superseded by a later run or belongs to an analysis line the project rejected.

Nothing was deleted, and no script was deleted -- every item below is regenerable.
Once you have confirmed the current figures and tables are what you want:

    rm -rf /FAST/gr10634/gaozy/aml_niche_net/LCC_proj/results/_retired

--------------------------------------------------------------------------------------
1. inferCNV CNV-proxy raw dumps                                                    25 MB
--------------------------------------------------------------------------------------
    tables/01_percnv/                    (390 files, 13 MB)

WHY: the 17p-loss CNV proxy was tested against the curated genotype in 07 and FAILED
(sensitivity 0/2, specificity 0/1). It defines no group in this project; P1B exists to
say so. This is its per-sample raw dump, read only by 01 itself when it assembles the two
aggregate tables.

NOT RETIRED, though they came off the same line: 01_arm_events_long.csv and
01_subclone_cnv.csv. They were moved here in a first pass and put back. 06_select_datasets.R
reads both (lines 55 and 110) to score datasets on subclonal structure -- a different use
from the refuted TP53 proxy, and retiring them would have broken a live script.

KEPT in results/tables/ deliberately: 01_parse_qc.csv and 02_*_all_datasets.csv. Those
carry the verdict, and P1B plus the clinician README cite them as the evidence that the
proxy failed. Deleting the evidence for a negative result would be the wrong cleanup.

REGENERATE: Rscript LCC_proj/scripts/01_parse_infercnv_regions.R
             (reads the inferCNV outputs, which are still on LARGE1)

--------------------------------------------------------------------------------------
2. p53 transcriptional-axis grouping                                               22 KB
--------------------------------------------------------------------------------------
    tables/05_p53_axis.csv
    tables/05_anchor_validation.csv
    logs/05_p53_axis.log, logs/05_p53_axis.rerun.log

WHY: an early attempt to define TP53 status from a transcriptional signature instead of
from the mutation. Superseded by 07 (curated genotype) + 08 (Numbat 17p LOH) + 09 (final
assignment), which is what every current figure uses. Read by no script downstream of 05
and cited in no document.

REGENERATE: Rscript LCC_proj/scripts/05_p53_axis.R

--------------------------------------------------------------------------------------
3. stale 4-dataset tables                                                          20 KB
--------------------------------------------------------------------------------------
    tables/02_arm_specificity.csv
    tables/02_sensitivity_grid.csv
    tables/02_validation.csv
    tables/02_sample_cnv_proxy.csv

WHY: 02 was re-run with --all_datasets (11 datasets, Jul 30 14:49), which writes to
*_all_datasets.csv so the earlier 4-dataset tables were never overwritten. Both versions
sat side by side, and 07/09 contain a fallback that silently reads the smaller one if the
larger is missing -- a hazard, not a backup.

NOTE: 06_select_datasets.R line 78 read the non-suffixed file as an optional 17p-anchor
column. It was patched to prefer *_all_datasets.csv, so a re-run of 06 now uses the
11-dataset table rather than silently dropping the column.

REGENERATE: Rscript LCC_proj/scripts/02_define_cnv_proxy.R   (without --all_datasets)

--------------------------------------------------------------------------------------
WHAT THIS BREAKS -- one script, deliberately
--------------------------------------------------------------------------------------
05_p53_axis.R will not run: it reads 02_sample_cnv_proxy.csv, and both of its own outputs
are retired here. That is the abandoned line itself, so the break is the point rather than
a side effect. Its header now says so. Every other script was checked and still has all of
its inputs; 06_select_datasets.R was patched to read 02_sample_cnv_proxy_all_datasets.csv.

Checked with:
    grep -ho 'LCC_TAB_DIR, "[^"]*"' LCC_proj/scripts/*.R | sed 's/.*"\(.*\)"/\1/' |
      sort -u | while read f; do [ -e "LCC_proj/results/tables/$f" ] || echo "MISSING $f"; done

--------------------------------------------------------------------------------------
4. CopyKat install log                                                            1.1 MB
--------------------------------------------------------------------------------------
    logs/install_copykat.log

WHY: CopyKat was evaluated as a third CNV caller and never used. inferCNV and Numbat are
what ran.

--------------------------------------------------------------------------------------
5. slurm array logs for the 03 per-cell pass                        440 files, 7.0 MB
--------------------------------------------------------------------------------------
    logs/lcc_scan/lcc_scan.3623864_*.{out,err}

WHY: all 220 array tasks exited clean -- 0 of 220 .err files contain "Error" or
"Execution halted" -- and their outputs (results/tables/03_detect/, 03_myeloid/) are on
disk and in use. The .err files are `set -x` traces of the job wrapper.

--------------------------------------------------------------------------------------
KEPT, though they looked like candidates
--------------------------------------------------------------------------------------
results/tables/03_detect/ (22 MB), 03_myeloid/    active intermediates. Aggregated into
                                                  04_*, but regenerating them means
                                                  re-running a 220-task sbatch array.
results/tables/04_pathway_sample_bin.csv          written by 04 and currently read by
                                                  nothing, but it is the
                                                  compartment-stratified twin of
                                                  04_pathway_sample.csv, which is in use.
                                                  Not refuted -- just not plotted yet.
results/figures/TP53 Kao discussion.pptx          input from the DSP group, not our output.
results/tables/F*, G*, P*, Q*, U*, V*             every current figure's source table.
