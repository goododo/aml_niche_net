#!/usr/bin/env bash
# 10_prealign_one_sample.sh ----
# Cluster per-sample BAM chain (single env = general_env). MULTI-RUN AWARE: a sample may
# be split into many SRR lanes (e.g. GSE227903); all lanes are downloaded and passed to
# STARsolo together as comma-separated lists. Compute node has internet -> download in-job.
# Steps: ENA-verified download (all lanes) -> STARsolo (sorted BAM w/ CB,UB) ->
#        Numbat pileup_and_phase (cellsnp-lite + Eagle2) -> run_numbat (30) -> cleanup.
#
# Usage:  bash scripts/02_malignancy/10_prealign_one_sample.sh <DATASET> <SAMPLE_ID> <SRR1,SRR2,...>
#   SAMPLE_ID must match the per-sample QC RDS basename so outputs line up.
#
# 3a-1 rewire (no logic change): config now comes from config_paths.sh, located via the .here
#   anchor (was: source ./c00_cluster_config.sh); the run_numbat call targets 30_numbat_run.R via
#   ${PROJECT_ROOT} (was: c30_run_numbat.R sibling of $0). BAM-chain logic byte-preserved.

set -euo pipefail

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
MAL_DIR="${PROJECT_ROOT}/scripts/02_malignancy"
# 10x bamtofastq (restores R1/R2 from a submitted cellranger BAM when SRA normalized fastq lacks
# the barcode read). Override with $BAMTOFASTQ; defaults to the static binary under the FAST parent.
BAMTOFASTQ="${BAMTOFASTQ:-$(dirname "${FAST_DIR}")/bamtofastq_linux}"

DATASET="${1:?need DATASET}"; SAMPLE_ID="${2:?need SAMPLE_ID}"; SRR_LIST="${3:?need comma-separated SRR list}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_PREFIX}"

SDIR="${WORK_ROOT}/${DATASET}/${SAMPLE_ID}"
FQ_DIR="${SDIR}/fastq"; STAR_OUT="${SDIR}/star"; NB_OUT="${SDIR}/numbat"
CNV_OUT="${CNV_ROOT}/numbat/${DATASET}/${SAMPLE_ID}"
SNP_OUT="${CNV_ROOT}/cellsnp/${DATASET}/${SAMPLE_ID}"
TAB_OUT="${MAL_TAB}/${DATASET}"
mkdir -p "${FQ_DIR}" "${STAR_OUT}" "${NB_OUT}" "${CNV_OUT}" "${SNP_OUT}" "${TAB_OUT}"
LOG="${SDIR}/run.log"; exec > >(tee -a "${LOG}") 2>&1
IFS=',' read -ra SRRS <<< "${SRR_LIST}"
echo "==== $(date) | ${DATASET} | ${SAMPLE_ID} | ${#SRRS[@]} lane(s): ${SRR_LIST} ===="

## ---- helpers ----
verify_file(){ local f=$1 eb=${2:-} em=${3:-}; [[ -s $f ]]||return 1
  if [[ -n $eb ]]; then local s; s=$(stat -c%s "$f"); [[ $s == "$eb" ]]||{ echo "      size mismatch $(basename "$f")"; return 1; }; fi
  if [[ -n $em ]]; then local m; m=$(md5sum "$f"|awk '{print $1}'); [[ $m == "$em" ]]||{ echo "      md5 mismatch $(basename "$f")"; return 1; }; fi
  return 0; }

# Download ONE run's mates into FQ_DIR as <ACC>_1.fastq.gz / <ACC>_2.fastq.gz.
# Preference: ENA fastq_ftp -> original 10x BAM (bamtofastq) -> prefetch/normalized SRA.
# NOTE: SRA "Normalized Format" strips 10x barcode reads. When fastq_ftp does not carry a real
# mate pair, the author's submitted cellranger BAM (barcodes in CB/UB tags) is the reliable
# source; bamtofastq restores R1/R2.
#
# TWO TRAPS, both hit by GSE289435/SRR32323361 and both fixed here:
#
#  (1) fastq_ftp NON-EMPTY BUT SINGLE-FILE. The BAM route used to fire only on an EMPTY fastq_ftp.
#      SRR32323361 advertises one 3.4G fastq_ftp file that is 91bp single-end -- the cDNA read with
#      the barcode read discarded. The old gate took the fastq route, downloaded 3.4G, and died at
#      "missing _1/_2" having thrown away the only usable source. A 10x run ALWAYS has >=2 mate
#      files, so a single-file listing is the normalized-format trap, not a valid download.
#
#  (2) ENA submitted_ftp EMPTY WHILE NCBI STILL HOLDS THE ORIGINAL. ENA lists no submitted_ftp for
#      SRR32323361, but NCBI's SRAFiles manifest carries MLL_29532.bam (17.4G, supertype="Original",
#      semantic_name="10X Genomics bam file") on sra-pub-src S3. Mirrors disagree about what they
#      kept, so exhausting ENA is not the same as the data being gone.
fetch_one(){
  local ACC=$1
  # resume: mates already restored (e.g. a prior run got past download/bamtofastq then failed at
  # STARsolo) -> reuse them, don't re-download the (possibly huge) source again
  if [[ -s "${FQ_DIR}/${ACC}_1.fastq.gz" && -s "${FQ_DIR}/${ACC}_2.fastq.gz" ]]; then
    echo "      mates already present for ${ACC} (${FQ_DIR}/${ACC}_1/2.fastq.gz); skip download"
    return
  fi
  local api="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes,submitted_ftp,submitted_format,submitted_bytes,submitted_md5&format=tsv"
  local rep; rep=$(curl -fsSL "${api}" || true)
  mapfile -t FTP < <(echo "${rep}"|awk -F'\t' 'NR>1&&$2!=""{print $2}'|tr ';' '\n'|sed '/^$/d')
  mapfile -t MD5 < <(echo "${rep}"|awk -F'\t' 'NR>1&&$3!=""{print $3}'|tr ';' '\n'|sed '/^$/d')
  mapfile -t BYT < <(echo "${rep}"|awk -F'\t' 'NR>1&&$4!=""{print $4}'|tr ';' '\n'|sed '/^$/d')
  local SUB SUBFMT SUBBYT SUBMD5
  SUB=$(echo "${rep}"|awk -F'\t' 'NR>1{print $5}'|head -1)
  SUBFMT=$(echo "${rep}"|awk -F'\t' 'NR>1{print $6}'|head -1)
  SUBBYT=$(echo "${rep}"|awk -F'\t' 'NR>1{print $7}'|head -1)
  SUBMD5=$(echo "${rep}"|awk -F'\t' 'NR>1{print $8}'|head -1)

  # Resolve the original 10x BAM from EITHER mirror before deciding the route. ENA first (it
  # publishes byte size + md5 in the same call), then NCBI's SRAFiles manifest.
  local BAM_URL="" BAM_BYTES="" BAM_MD5=""
  if [[ -n "${SUB}" && ( "${SUBFMT^^}" == *BAM* || "${SUB}" == *.bam ) ]]; then
    local subbam sub_i
    subbam=$(echo "${SUB}"|tr ';' '\n'|grep -i '\.bam$'|head -1)
    # matching submitted_bytes/md5 entry (fields are ';'-aligned when multiple files)
    sub_i=$(echo "${SUB}"|tr ';' '\n'|grep -in '\.bam$'|head -1|cut -d: -f1)
    BAM_URL="https://${subbam}"
    BAM_BYTES=$(echo "${SUBBYT}"|tr ';' '\n'|sed -n "${sub_i}p")
    BAM_MD5=$(echo "${SUBMD5}"|tr ';' '\n'|sed -n "${sub_i}p")
  else
    # supertype="Original" is the author's submission; the "Primary ETL"/normalized entries in the
    # same manifest are the lossy derivatives we are trying to avoid, so match on it explicitly.
    local nx rec
    nx=$(curl -fsSL "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=${ACC}" || true)
    rec=$(printf '%s' "${nx}" | tr '<' '\n' \
          | grep -i '^SRAFile ' | grep 'supertype="Original"' | grep -iE 'filename="[^"]*\.bam"' | head -1)
    if [[ -n "${rec}" ]]; then
      BAM_URL=$(sed -n 's/.* url="\([^"]*\)".*/\1/p' <<<"${rec}")
      BAM_BYTES=$(sed -n 's/.* size="\([^"]*\)".*/\1/p' <<<"${rec}")
      BAM_MD5=$(sed -n 's/.* md5="\([^"]*\)".*/\1/p' <<<"${rec}")
      [[ -n "${BAM_URL}" ]] && echo "      ENA has no submitted BAM; NCBI SRAFiles does (${BAM_URL##*/})"
    fi
  fi

  # -lt 2, not -eq 0: see trap (1) in the header. One fastq for a 10x run is the normalized trap.
  if [[ ${#FTP[@]} -lt 2 ]]; then
    [[ ${#FTP[@]} -eq 1 ]] && echo "      fastq_ftp lists ONE file for a 10x run (barcode read stripped); ignoring it"
    # (a) original 10x BAM -> bamtofastq (barcodes survive here; normalized SRA fastq does not)
    if [[ -n "${BAM_URL}" ]]; then
      local exp_bytes="${BAM_BYTES}"
      echo "      original 10x BAM (${BAM_URL##*/}, expect ${exp_bytes:-?} bytes); restoring via bamtofastq"
      [[ -x "${BAMTOFASTQ}" ]] || { echo "ERROR: bamtofastq not found/executable at '${BAMTOFASTQ}' (set \$BAMTOFASTQ)"; exit 1; }
      local bam="${SDIR}/${ACC}.bam"; local t=0 stall=0 prev=-1
      # wget fetches ONLY the submitted BAM (prefetch --type all also pulls the normalized SRA + its
      # refseq deps -- slower/wasteful). Resume with -c; accept only when size == submitted_bytes.
      # NEVER delete a partial. Generous retry budget so a flaky ENA link completes hands-off; a
      # truly dead transfer fails PRESERVING the partial, so a plain re-submit resumes where it left.
      while :; do
        local have; have=$(stat -c%s "${bam}" 2>/dev/null || echo 0)
        if [[ -n "${exp_bytes}" && "${have}" == "${exp_bytes}" ]]; then break; fi
        t=$((t+1)); [[ $t -gt 50 ]] && { echo "ERROR: ${ACC} BAM incomplete after 50 tries (${have}/${exp_bytes:-?}); partial KEPT -- re-submit to resume"; exit 1; }
        if [[ "${have}" == "${prev}" ]]; then stall=$((stall+1)); else stall=0; fi
        [[ $stall -ge 15 ]] && { echo "ERROR: ${ACC} BAM stalled 15x at ${have}/${exp_bytes:-?}; partial KEPT -- re-submit to resume"; exit 1; }
        prev="${have}"
        echo "      [dl BAM ${t}] ${BAM_URL##*/}  (have ${have}/${exp_bytes:-?})"
        wget -c -q --tries=10 --timeout=300 --retry-connrefused --waitretry=20 \
             -O "${bam}" "${BAM_URL}" || true
      done
      # INTEGRITY: size match is NOT enough -- a resumed transfer can be byte-complete yet corrupt
      # (seen on GSM8791438: bgzf inflate failure). Verify md5 when the mirror provides it, plus a
      # cheap BAM structural check. A CORRUPT file is deleted and re-fetched (unlike a partial, it
      # is worthless and would only fail again downstream).
      local exp_md5="${BAM_MD5}"
      if [[ -n "${exp_md5}" ]]; then
        echo "      verifying md5 (${exp_md5})"
        local got_md5; got_md5=$(md5sum "${bam}" | awk '{print $1}')
        if [[ "${got_md5}" != "${exp_md5}" ]]; then
          echo "ERROR: ${ACC} BAM md5 mismatch (got ${got_md5}); CORRUPT file removed -- re-submit to download afresh"
          rm -f "${bam}"; exit 1
        fi
      fi
      if ! samtools quickcheck -v "${bam}" 2>/dev/null; then
        echo "ERROR: ${ACC} BAM failed samtools quickcheck; CORRUPT file removed -- re-submit to download afresh"
        rm -f "${bam}"; exit 1
      fi
      echo "      BAM complete + verified (${exp_bytes} bytes); running bamtofastq"
      local btf="${SDIR}/btf_${ACC}"; rm -rf "${btf}"
      "${BAMTOFASTQ}" --nthreads="${THREADS}" "${bam}" "${btf}" \
        || { echo "ERROR: bamtofastq failed for ${ACC}"; exit 1; }
      # concat all R1 chunks -> _1 (barcode), all R2 chunks -> _2 (cDNA); consistent sort keeps pairing
      cat $(find "${btf}" -name '*_R1_*.fastq.gz'|sort) > "${FQ_DIR}/${ACC}_1.fastq.gz"
      cat $(find "${btf}" -name '*_R2_*.fastq.gz'|sort) > "${FQ_DIR}/${ACC}_2.fastq.gz"
      rm -f "${bam}"; rm -rf "${btf}"
      [[ -s "${FQ_DIR}/${ACC}_1.fastq.gz" && -s "${FQ_DIR}/${ACC}_2.fastq.gz" ]] \
        || { echo "ERROR: ${ACC} bamtofastq produced no R1/R2 (is the BAM a 10x cellranger output?)"; exit 1; }
      return
    fi
    # (b) last resort: prefetch normalized SRA (may lack the 10x barcode read -> caller errors out)
    echo "      no usable mate pair and no original BAM on either mirror for ${ACC};"
    echo "      prefetch fallback (WARNING: normalized SRA may drop barcodes)"
    mkdir -p "${SDIR}/sra"
    prefetch --max-size 500G -O "${SDIR}/sra" "${ACC}"
    local sra; sra=$(find "${SDIR}/sra" -name "${ACC}*.sra"|head -1)
    fasterq-dump --split-files --include-technical --threads "${THREADS}" -O "${FQ_DIR}" "${sra}"
    ( pigz -p "${THREADS}" "${FQ_DIR}/${ACC}"*.fastq 2>/dev/null || gzip "${FQ_DIR}/${ACC}"*.fastq )
    rm -rf "${SDIR}/sra"; return
  fi
  local i
  for i in "${!FTP[@]}"; do
    local url=${FTP[$i]} md5=${MD5[$i]:-} byt=${BYT[$i]:-} dest="${FQ_DIR}/$(basename "${FTP[$i]}")"
    local t=0
    until verify_file "${dest}" "${byt}" "${md5}"; do
      t=$((t+1)); [[ $t -gt 3 ]] && { echo "ERROR: ${ACC} $(basename "${dest}") integrity fail x3"; exit 1; }
      [[ $t -ge 2 ]] && rm -f "${dest}"
      echo "      [dl ${t}] $(basename "${dest}")"
      wget -c -q -O "${dest}" "https://${url}"
    done
  done
}

read_len(){ zcat "$1" 2>/dev/null | head -2 | tail -1 | tr -d '\n' | wc -c || true; }

## ---- 1. download all lanes + assign CB/cDNA per lane, build STARsolo comma-lists ----
echo "[1] download + assign mates (${#SRRS[@]} lanes)"
CB_FILES=(); CDNA_FILES=(); CB_LEN=0
for ACC in "${SRRS[@]}"; do
  ACC="$(echo "$ACC" | tr -d '[:space:]')"; [[ -z $ACC ]] && continue
  echo "    lane ${ACC}"
  fetch_one "${ACC}"
  r1="${FQ_DIR}/${ACC}_1.fastq.gz"; r2="${FQ_DIR}/${ACC}_2.fastq.gz"
  [[ -s $r1 && -s $r2 ]] || { echo "ERROR: ${ACC} missing _1/_2 after download"; exit 1; }
  L1=$(read_len "$r1"); L2=$(read_len "$r2")
  if   [[ $L1 -eq 26 || $L1 -eq 28 ]]; then cb=$r1; cdna=$r2; cblen=$L1
  elif [[ $L2 -eq 26 || $L2 -eq 28 ]]; then cb=$r2; cdna=$r1; cblen=$L2
  else echo "ERROR: ${ACC} no 26/28bp barcode read (L1=${L1} L2=${L2}); vet with 04"; exit 1; fi
  echo "      CB=$(basename "$cb") (${cblen}bp)  cDNA=$(basename "$cdna")"
  if [[ $CB_LEN -eq 0 ]]; then CB_LEN=$cblen
  elif [[ $CB_LEN -ne $cblen ]]; then echo "ERROR: mixed barcode lengths across lanes (${CB_LEN} vs ${cblen})"; exit 1; fi
  CB_FILES+=("$cb"); CDNA_FILES+=("$cdna")
done
[[ ${#CDNA_FILES[@]} -gt 0 ]] || { echo "ERROR: no lanes downloaded"; exit 1; }
if [[ $CB_LEN -eq 26 ]]; then CHEM=v2; UMILEN=10; WL=$WHITELIST_V2; else CHEM=v3; UMILEN=12; WL=$WHITELIST_V3; fi
CDNA_CSV=$(IFS=,; echo "${CDNA_FILES[*]}")
CB_CSV=$(IFS=,;   echo "${CB_FILES[*]}")
echo "    chemistry=${CHEM} CB=16 UMI=${UMILEN} lanes=${#CDNA_FILES[@]}"
[[ -s $WL ]] || { echo "ERROR: missing whitelist ${WL}"; exit 1; }

## ---- 2. STARsolo over ALL lanes -> UNSORTED BAM (CB/UB) + filtered matrix, then samtools sort ----
# STARsolo emits CB/UB tags fine on Unsorted output; we do the coordinate sort with samtools (step
# 2a) because STAR's in-RAM sort OOMs on huge 10x BAMs. --outTmpDir still goes on the BIG disk
# (/LARGE1, under SDIR) for STAR's own scratch. STAR creates readFilesCommand FIFOs inside outTmpDir
# and some Lustre mounts reject FIFOs -> self-test mkfifo first; if it fails, fall back to
# pre-decompressing the FASTQs (no readFilesCommand, so no FIFOs needed).
echo "[2] STARsolo (merged lanes; UNSORTED BAM, temp on big disk)"
STAR_TMP="${SDIR}/startmp"; rm -rf "${STAR_TMP}"; mkdir -p "${STAR_TMP}"
USE_FIFO=1
if ! ( mkfifo "${STAR_TMP}/.fifotest" 2>/dev/null ); then USE_FIFO=0; fi
rm -f "${STAR_TMP}/.fifotest"

# BAM sorting strategy, decided AT RUNTIME from what this STAR build supports:
#  (i)  STAR refuses CB/UB on Unsorted output UNLESS --soloAddTagsToUnsorted is available. When it
#       is, we emit Unsorted + tags and coordinate-sort with samtools: disk-backed, BOUNDED RAM,
#       so it never OOMs regardless of read count (this cohort demanded 48G -> 72G -> 143G from
#       STAR's in-RAM sort and kept climbing).
#  (ii) Otherwise fall back to STAR's own sorted output with a large --limitBAMsortRAM. That path
#       needs a FAT node (see 11.slurm); set STAR_SORT_RAM to fit the job's memory.
STAR_SORT_RAM="${STAR_SORT_RAM:-200000000000}"     # only used by the fallback (ii)
if STAR --help 2>&1 | grep -q -- '--soloAddTagsToUnsorted'; then
  SORT_MODE="samtools"
  SORT_ARGS=( --outSAMtype BAM Unsorted --soloAddTagsToUnsorted yes )
  echo "    [sort] STAR supports --soloAddTagsToUnsorted -> Unsorted + samtools sort (bounded RAM)"
else
  SORT_MODE="star"
  SORT_ARGS=( --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM "${STAR_SORT_RAM}" --outBAMsortingBinsN 100 )
  echo "    [sort] STAR lacks --soloAddTagsToUnsorted -> STAR sorted, limitBAMsortRAM=${STAR_SORT_RAM} (needs a fat node)"
fi
STAR_COMMON=( --runMode alignReads --genomeDir "${STAR_INDEX}" --genomeLoad NoSharedMemory
  --soloType CB_UMI_Simple --soloCBwhitelist "${WL}"
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen "${UMILEN}"
  --soloFeatures Gene --soloCellFilter EmptyDrops_CR
  "${SORT_ARGS[@]}" --outSAMattributes NH HI AS nM CB UB
  --runThreadN "${THREADS}" --outTmpDir "${STAR_TMP}/run" --outFileNamePrefix "${STAR_OUT}/" )

if [[ ${USE_FIFO} -eq 1 ]]; then
  STAR "${STAR_COMMON[@]}" --readFilesIn "${CDNA_CSV}" "${CB_CSV}" --readFilesCommand zcat
else
  echo "    [2!] FIFOs unsupported on ${STAR_TMP}; pre-decompressing FASTQs (no readFilesCommand)"
  DEC_CDNA=(); DEC_CB=()
  for f in "${CDNA_FILES[@]}"; do d="${f%.gz}"; [[ -s $d ]] || zcat "$f" > "$d"; DEC_CDNA+=("$d"); done
  for f in "${CB_FILES[@]}";   do d="${f%.gz}"; [[ -s $d ]] || zcat "$f" > "$d"; DEC_CB+=("$d");   done
  DEC_CDNA_CSV=$(IFS=,; echo "${DEC_CDNA[*]}"); DEC_CB_CSV=$(IFS=,; echo "${DEC_CB[*]}")
  STAR "${STAR_COMMON[@]}" --readFilesIn "${DEC_CDNA_CSV}" "${DEC_CB_CSV}"
  rm -f "${DEC_CDNA[@]}" "${DEC_CB[@]}"
fi
rm -rf "${STAR_TMP}"

UNSORTED="${STAR_OUT}/Aligned.out.bam"
BAM="${STAR_OUT}/Aligned.sortedByCoord.out.bam"
FILT="${STAR_OUT}/Solo.out/Gene/filtered"; BARCODES="${FILT}/barcodes.tsv"
if [[ "${SORT_MODE}" == "samtools" ]]; then
  # memory-SAFE coordinate sort: samtools external merge-sort, RAM capped at ~THREADS x 2G with the
  # rest spilled to disk (-T on the big /LARGE1 disk). Unlike STAR's in-RAM sort, this cannot OOM.
  [[ -s $UNSORTED ]] || { echo "ERROR: STARsolo unsorted BAM missing"; exit 1; }
  echo "[2a] samtools sort (external, disk-backed; RAM capped)"
  samtools sort -@ "${THREADS}" -m 2G -T "${SDIR}/samsort" -o "${BAM}" "${UNSORTED}" \
    || { echo "ERROR: samtools sort failed"; exit 1; }
  rm -f "${UNSORTED}"
else
  echo "[2a] STAR already produced the coordinate-sorted BAM"
fi
[[ -s $BAM ]] || { echo "ERROR: sorted BAM missing"; exit 1; }
[[ -f ${BARCODES}.gz ]] && gunzip -k "${BARCODES}.gz"
samtools index -@ "${THREADS}" "${BAM}"
cp -r "${FILT}" "${CNV_OUT}/star_filtered"
cp "${STAR_OUT}/Solo.out/Gene/Summary.csv" "${CNV_OUT}/Solo_Gene_Summary.csv" 2>/dev/null || true

# NOTE (consensus/R2 layer): malignancy is called on STARsolo's own cell set here.
# Intersect these barcodes with the per-sample QC RDS barcodes when merging the 3-way
# consensus, so malignant calls land on the QC-passing cells.

[[ ${CLEAN_FASTQ} == true ]] && { echo "[2b] removing FASTQ"; rm -rf "${FQ_DIR}"; }

## ---- 3. Numbat pileup_and_phase (cellsnp-lite + Eagle2) -> allele counts ----
echo "[3] pileup_and_phase"
PHASE=$(Rscript -e 'cat(system.file("bin/pileup_and_phase.R", package="numbat"))')
[[ -s $PHASE ]] || { echo "ERROR: pileup_and_phase.R not found in ${ENV_PREFIX}"; exit 1; }
Rscript "${PHASE}" --label "${SAMPLE_ID}" --samples "${SAMPLE_ID}" \
  --bams "${BAM}" --barcodes "${BARCODES}" --outdir "${NB_OUT}" \
  --gmap "${NB_GMAP}" --eagle "${NB_EAGLE_BIN}" --snpvcf "${NB_SNP_VCF}" \
  --paneldir "${NB_PANEL_DIR}" --ncores "${THREADS}"
ALLELE=$(find "${NB_OUT}" -name '*allele_counts.tsv.gz'|head -1)
[[ -s $ALLELE ]] || { echo "ERROR: allele_counts not produced (check NB_SNP_VCF contig naming vs BAM)"; exit 1; }
cp "${ALLELE}" "${CNV_OUT}/"
[[ -d "${NB_OUT}/pileup" ]] && cp -r "${NB_OUT}/pileup" "${SNP_OUT}/" 2>/dev/null || true

[[ ${KEEP_BAM} != true ]] && { echo "[3b] removing BAM"; rm -f "${BAM}" "${BAM}.bai"; }

## ---- 4. Numbat CNV/clone inference (30) ----
echo "[4] run_numbat"
Rscript "${MAL_DIR}/30_numbat_run.R" \
  --sample "${SAMPLE_ID}" --matrix "${CNV_OUT}/star_filtered" \
  --allele "${CNV_OUT}/$(basename "${ALLELE}")" --outdir "${CNV_OUT}/numbat" \
  --genome "${NB_GENOME_VERSION}" --ncores "${THREADS}" --ref_rds "${NB_REF_RDS}"
cp "${CNV_OUT}"/numbat/*__numbat_percell.csv "${TAB_OUT}/" 2>/dev/null || true

echo "[done] ${SAMPLE_ID}: outputs -> ${CNV_OUT} ; per-cell table -> ${TAB_OUT}"
rm -rf "${SDIR}"