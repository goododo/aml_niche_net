#!/usr/bin/env bash
# 02_fetch_srr_table.sh ----
# GEO series -> SRA run table (run_accession + experiment/sample titles), so you can
# map runs to the per-sample QC RDS names by hand. RUN ON A NODE WITH INTERNET.
#
# IMPORTANT: pysradb has NO 'gse-to-srr' subcommand. The reliable path is two steps:
#   (1) gse-to-srp   GSE -> SRP study accession
#   (2) metadata SRP --detailed   -> full run table (run_accession + titles)
#
# Prereq:  pysradb in the active env (you already have it in general_env).
# Usage:   bash scripts/02_malignancy/02_fetch_srr_table.sh GSE227903
#
# 3a-1 rewire (no logic change): ENV_PREFIX / QC_RDS_DIR now come from config_paths.sh (located
#   via the .here anchor) instead of being hard-coded in messages. pysradb logic byte-preserved.

set -uo pipefail

## ---- locate project root via the .here anchor (bash analogue of R here::here) ----
_find_here() {
  local d; d="$(cd "${1:-$PWD}" 2>/dev/null && pwd)" || return 1
  while [[ "$d" != "/" ]]; do
    [[ -e "$d/.here" ]] && { printf '%s\n' "$d"; return 0; }
    d="$(dirname "$d")"
  done; return 1; }
_self="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(_find_here "${_self}")" || true
[[ -n "${PROJECT_ROOT}" ]] || { echo "ERROR: no .here anchor found walking up from ${_self}"; exit 1; }
source "${PROJECT_ROOT}/scripts/config/config_paths.sh"

GSE="${1:?need GSE accession, e.g. GSE227903}"
OUT="srr_table_${GSE}.tsv"

command -v pysradb >/dev/null 2>&1 || {
  echo "ERROR: pysradb not on PATH. Activate the env first:"
  echo "       conda activate ${ENV_PREFIX}"
  exit 1; }

## ---- 1. GSE -> SRP ----
echo "[1] ${GSE} -> SRP"
SRP=$(pysradb gse-to-srp "${GSE}" 2>/dev/null | grep -oE 'SRP[0-9]+' | head -1 || true)
if [[ -z "${SRP}" ]]; then
  echo "ERROR: could not resolve an SRP for ${GSE}. Raw pysradb output below:"
  pysradb gse-to-srp "${GSE}"
  echo
  echo "  Fallback: open NCBI SRA Run Selector for ${GSE}, click 'Metadata',"
  echo "  and save that table as ${OUT} (it has Run + sample columns)."
  exit 1
fi
echo "    SRP=${SRP}"

## ---- 2. SRP -> detailed run table (run_accession + sample/experiment titles) ----
echo "[2] ${SRP} -> detailed run table"
if ! pysradb metadata "${SRP}" --detailed --saveto "${OUT}"; then
  echo "ERROR: 'pysradb metadata ${SRP} --detailed' failed (see message above)."
  echo "  Fallback: SRA Run Selector metadata table for ${GSE}."
  exit 1
fi
[[ -s "${OUT}" ]] || { echo "ERROR: ${OUT} is empty"; exit 1; }

echo "[3] Wrote ${OUT}. Preview (run_accession + titles):"
column -t -s$'\t' "${OUT}" 2>/dev/null | head -40

cat <<NOTE

NEXT STEPS
  1) Find the columns holding the SRR run (run_accession) and the GSM / sample title.
     In --detailed output these are usually: run_accession, experiment_title,
     and a sample_* / *_title column. Use whichever identifies the patient/timepoint.
  2) Match those against the per-sample RDS names under
       ${QC_RDS_DIR}/${GSE}/
  3) Save a clean two-column sheet (tab-separated, '#'-comments allowed):
        <sample_id>\t<srr>
     where <sample_id> EXACTLY matches the RDS basename (no .rds). File: samples_${GSE}.tsv
       - If a sample has several SRR runs, list each SRR on its own line with the SAME
         sample_id; 10's downloader will pull them, but for STARsolo you'll want one
         sample's runs concatenated -- tell me if GSE227903 is multi-run per sample.
  4) Vet read structure BEFORE downloading:  bash scripts/02_malignancy/04_vet_read_structure.sh samples_${GSE}.tsv

ALT (one-shot, also gives FASTQ URLs): once you have SRP=${SRP}, you can pull an ENA
table with run_accession + sample_title + fastq_ftp in a single call:
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRP}&result=read_run&fields=run_accession,sample_title,experiment_title,fastq_ftp,fastq_bytes&format=tsv"
NOTE