#!/bin/bash
# INPUT : network (run on the LOGIN node; compute nodes have no network); conda env /FAST/gr10634/gaozy/general_env (R 4.4.3)
# OUTPUT: /FAST/gr10634/gaozy/rlib_scdbl_devel/  a SEPARATE R library holding scDblFinder 1.27.12 (GitHub commit 1a0831bd,
#         origin-label fix of issue #145) and the Bioconductor 3.23 packages it needs that Bioc 3.20 cannot provide:
#         assorthead 1.6.4, beachmat 2.28.0, BiocNeighbors 2.6.0, biocmake 1.4.0, Rigraphlib 1.4.0, scrapper 1.6.3.
#         Source tarballs are kept in _src/ with their sha256 in _src/SHA256SUMS; the install log is _install.log.
# what it does: builds the devel scDblFinder for the doublet-pairing smoke comparison (decision of 2026-11-09). It never
#   writes to the user library or to general_env: every install goes to the separate library, which only a job that
#   prepends it to .libPaths() will see. Packages are installed in dependency order; the first failure stops the script.
#   usage: bash scripts/dbl/00_install_scdbl_devel.sh      (the caller wraps it in `timeout 2700` = the 45-minute box)
set -eu -o pipefail
LIB=/FAST/gr10634/gaozy/rlib_scdbl_devel
SRC=$LIB/_src
SHA=1a0831bd096dad402af8c3b3e4d6e472b88ffab8
BIOC=https://bioconductor.org/packages/3.23/bioc/src/contrib
PKGS=(assorthead_1.6.4 beachmat_2.28.0 BiocNeighbors_2.6.0 biocmake_1.4.0 Rigraphlib_1.4.0 scrapper_1.6.3)
mkdir -p "$SRC"
log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "downloading ${#PKGS[@]} Bioconductor 3.23 tarballs and scDblFinder@${SHA:0:8}"
for p in "${PKGS[@]}"; do [[ -s $SRC/$p.tar.gz ]] || curl -fsSL -o "$SRC/$p.tar.gz" "$BIOC/$p.tar.gz"; done
[[ -s $SRC/scDblFinder_$SHA.tar.gz ]] || curl -fsSL -o "$SRC/scDblFinder_$SHA.tar.gz" "https://github.com/plger/scDblFinder/archive/$SHA.tar.gz"
(cd "$SRC" && sha256sum ./*.tar.gz > SHA256SUMS)

module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /FAST/gr10634/gaozy/general_env
export MAKEFLAGS=-j8                                  # parallel C++ compilation inside each package
cd "$SRC"                                             # never run R from inside the repo or next to downloaded code

for p in "${PKGS[@]}"; do
  log "installing $p"
  # .libPaths(c(LIB, ...)): headers of packages installed earlier in this loop (assorthead, beachmat, ...) win over the old ones
  Rscript -e ".libPaths(c('$LIB', .libPaths())); install.packages('$SRC/$p.tar.gz', repos = NULL, type = 'source', lib = '$LIB')"
  Rscript -e ".libPaths(c('$LIB', .libPaths())); stopifnot(requireNamespace('${p%%_*}', lib.loc = '$LIB', quietly = TRUE)); cat('  loaded', '${p%%_*}', as.character(packageVersion('${p%%_*}', lib.loc = '$LIB')), '\n')"
done

log "installing scDblFinder@${SHA:0:8}"
rm -rf "$SRC/scDblFinder-$SHA"; tar -xzf "$SRC/scDblFinder_$SHA.tar.gz" -C "$SRC"
Rscript -e ".libPaths(c('$LIB', .libPaths())); install.packages('$SRC/scDblFinder-$SHA', repos = NULL, type = 'source', lib = '$LIB')"
Rscript -e ".libPaths(c('$LIB', .libPaths())); suppressPackageStartupMessages(library(scDblFinder, lib.loc = '$LIB'));
  cat('  scDblFinder', as.character(packageVersion('scDblFinder')), 'from', find.package('scDblFinder'), '\n');
  cat('  xgboost', as.character(packageVersion('xgboost')), ' scrapper', as.character(packageVersion('scrapper')), '\n')"
log "done"
