#!/bin/bash
# =============================================================================
# submit_00_prebuild_ext_ref.sbatch   (job J0 of the v2 re-run)
#
# Rebuilds the inferCNV external healthy reference cache so it carries the
# per-cell lineage label and hierarchy bin. Without those fields
# .ext_ref_lineage_block() cannot split the reference, and get_external_ref()
# deliberately refuses the stale cache rather than silently falling back to a
# single reference group -- which is the exact behaviour INFERCNV_REF_PER_LINEAGE
# exists to remove.
#
# Deterministic: same INFERCNV_EXT_REF_SEED, so the same cells are selected as
# in the v1 cache; only label/bin are added.
#
# Needs a big-memory node: the source object
# (BoneMarrowMap_Annotated_Dataset_expandedFeatures.rds) is 3.4 GB on disk and
# OOM-kills a login node when loaded.
#
# Usage: sbatch scripts/02_malignancy/submit_00_prebuild_ext_ref.sbatch
# =============================================================================
#SBATCH --job-name=J0_extref
#SBATCH -p gr10634c
#SBATCH --rsc p=1:t=1:c=4:m=512G
#SBATCH --time=4:00:00
#SBATCH --output=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J0_extref_%A.out
#SBATCH --error=/FAST/gr10634/gaozy/aml_niche_net/00_project/logs/J0_extref_%A.err

set -euo pipefail

module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env

cd /FAST/gr10634/gaozy/aml_niche_net
Rscript scripts/02_malignancy/44_infercnv_run_one.R --prebuild_ext_ref

# Report the resulting lineage split so the log is self-verifying.
Rscript -e '
suppressMessages(source(here::here("scripts","02_malignancy","00_infercnv_common.R")))
ext <- get_external_ref()
cat("cache fields:", paste(names(ext), collapse=","), "\n")
cat("cells:", length(ext$cells), "\n")
print(table(ext$bin, useNA="ifany"))
b1 <- .ext_ref_lineage_block(ext)
b2 <- .ext_ref_lineage_block(ext, drop_bins=c("T_NK","B_Plasma"))
cat(sprintf("external route : %d cells / %d groups\n", length(b1$idx), length(unique(b1$groups))))
cat(sprintf("autologous aug : %d cells / %d groups (+1 autologous lymphoid)\n",
            length(b2$idx), length(unique(b2$groups))))
'
