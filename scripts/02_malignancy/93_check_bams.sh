#!/usr/bin/env bash
# 93_check_bams.sh ----
# Verify every downloaded submitted BAM under WORK_ROOT/<dataset>/: size vs ENA submitted_bytes,
# structural integrity (samtools quickcheck), and optionally md5 (slow on 16-34GB files).
# REPORT ONLY -- it never deletes; it prints the exact rm commands for anything bad so you decide.
#
# Usage:
#   bash scripts/02_malignancy/93_check_bams.sh GSE289435          # size + quickcheck (fast)
#   bash scripts/02_malignancy/93_check_bams.sh GSE289435 --md5    # + md5 (minutes per file)

set -uo pipefail
HERE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${HERE_DIR%/scripts/*}"
source "${PROJECT_ROOT}/scripts/config/config_paths.sh" >/dev/null 2>&1

DATASET="${1:?usage: 93_check_bams.sh <DATASET> [--md5]}"
DO_MD5=0; [[ "${2:-}" == "--md5" ]] && DO_MD5=1
DIR="${WORK_ROOT}/${DATASET}"
[[ -d "$DIR" ]] || { echo "no work dir: $DIR"; exit 1; }

printf '%-26s %-14s %-16s %-16s %s\n' SAMPLE STATUS ON_DISK EXPECTED NOTE
bad=()
for d in "$DIR"/*/; do
  s=$(basename "$d")
  bam=$(ls "$d"*.bam 2>/dev/null | head -1)
  [[ -n "$bam" ]] || continue
  acc=$(basename "$bam" .bam)

  rep=$(curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=submitted_ftp,submitted_bytes,submitted_md5&format=tsv" 2>/dev/null)
  sub=$(echo "$rep"   | awk -F'\t' 'NR>1{print $2}' | head -1)
  exp=$(echo "$rep"   | awk -F'\t' 'NR>1{print $3}' | tr ';' '\n' | head -1)
  emd5=$(echo "$rep"  | awk -F'\t' 'NR>1{print $4}' | tr ';' '\n' | head -1)
  have=$(stat -c%s "$bam" 2>/dev/null || echo 0)

  status="OK"; note=""
  if [[ -n "$exp" && "$have" != "$exp" ]]; then
    status="PARTIAL"; note="resume with a re-submit (keep the file)"
  elif ! samtools quickcheck -v "$bam" 2>/dev/null; then
    status="CORRUPT"; note="failed quickcheck -> delete + re-download"
    bad+=("$bam")
  elif [[ $DO_MD5 -eq 1 && -n "$emd5" ]]; then
    gmd5=$(md5sum "$bam" | awk '{print $1}')
    if [[ "$gmd5" != "$emd5" ]]; then
      status="CORRUPT"; note="md5 mismatch -> delete + re-download"
      bad+=("$bam")
    else
      note="md5 verified"
    fi
  fi
  printf '%-26s %-14s %-16s %-16s %s\n' "$s" "$status" "$have" "${exp:-?}" "$note"
done

echo
if [[ ${#bad[@]} -gt 0 ]]; then
  echo "=== ${#bad[@]} corrupt file(s). Review, then run these to remove them: ==="
  for f in "${bad[@]}"; do echo "rm -f '$f'"; done
  echo "(a PARTIAL file is fine -- leave it, the pipeline resumes it)"
else
  echo "no corrupt BAMs found (PARTIAL entries, if any, will resume normally)"
fi
