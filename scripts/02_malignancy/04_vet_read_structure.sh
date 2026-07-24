#!/usr/bin/env bash
# 04_vet_read_structure.sh ----
# MANDATORY pre-flight: check read layout BEFORE downloading, so we never waste a full
# download on a barcode-stripped/non-10x deposit (cf. GSE201966). Vets the FIRST lane of
# each sample (all lanes of one library share chemistry). RUN ON A NODE WITH INTERNET.
#
# Sheet format (tab-separated; '#'-comments allowed):  <sample_id>\t<srr1,srr2,...>
# Usage:  bash scripts/02_malignancy/04_vet_read_structure.sh samples_GSE227903.tsv
#
# 3a-1 rewire (no logic change): ENV_PREFIX now comes from config_paths.sh, located via the
#   .here anchor (bash analogue of R here::here) instead of being hard-coded. Vet logic and the
#   `set -uo pipefail` flags (deliberately no -e) are byte-preserved.

set -uo pipefail

## ---- locate project root via the .here anchor (bash analogue of R here::here) ----
# Plain .sh: start the search at THIS script's own on-disk dir (BASH_SOURCE), so resolution is
# independent of CWD. Fail loud if no anchor is found rather than using a wrong path.
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

SHEET="${1:?need sample sheet tsv}"
[[ -s "${SHEET}" ]] || { echo "ERROR: sheet not found: ${SHEET}"; exit 1; }
command -v vdb-dump >/dev/null 2>&1 || source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate "${ENV_PREFIX}" 2>/dev/null || true

printf "%-16s %-6s %-13s %-16s %s\n" "sample_id" "lanes" "first_srr" "READ_LEN" "verdict"
printf '%.0s-' {1..92}; echo

n_ok=0; n_bad=0
while IFS=$'\t' read -r SID SRR_LIST _rest; do
  [[ -z "${SID:-}" || "${SID}" == \#* ]] && continue
  IFS=',' read -ra SRRS <<< "${SRR_LIST}"
  nlane=${#SRRS[@]}
  SRR="$(echo "${SRRS[0]}" | tr -d '[:space:]')"
  rl=$(vdb-dump -R 1 -C READ_LEN "${SRR}" 2>/dev/null | sed -n 's/.*READ_LEN: *//p' | head -1)
  if [[ -z "${rl}" ]]; then
    v="QUERY_FAILED (network? wrong SRR?)"; ((n_bad++))
  elif echo "${rl}" | grep -qE '(^|, )(26|28)(,|$)'; then
    v="OK_10x (26/28bp barcode read present)"; ((n_ok++))
  else
    v="SUSPECT (no 26/28bp barcode read -- inspect before download)"; ((n_bad++))
  fi
  printf "%-16s %-6s %-13s %-16s %s\n" "${SID}" "${nlane}" "${SRR}" "${rl:-NA}" "${v}"
done < "${SHEET}"

echo
echo "Summary: ${n_ok} OK_10x, ${n_bad} needing attention."
echo "Note: 28bp barcode read => v3 (UMI 12); 26bp => v2 (UMI 10)."
echo "      Seq-Well / non-10x (e.g. GSE116256) reads as SUSPECT and needs custom geometry."
