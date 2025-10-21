#!/usr/bin/env bash
# Safer VCF → consensus FASTA (+optional Bowtie2 & ISSL) for WSL2/low-I/O environments
# - Robust tmp placement (defaults under OUT/.tmp)
# - Controllable extractOfftargets fan-in and retries (to avoid EIO)
# - Defensive ulimit/open-files capping
# - Atomic-ish writes and clearer logging
#
# Usage:
# VCF_To_Consensus_safe.sh \
# --vcf-dir DIR --ref GRCh38.fa --out OUTDIR \
# [--rename TSV] [--threads N] [--build-indexes yes|no] [--haplotypes both|1|2|iupac] \
# [--tmp DIR] [--issl-threads N] [--issl-max-open-files N] [--issl-retries N]
#
# Suggested on WSL2:
# place OUT on Linux FS (e.g., /home/<user>), NOT /mnt/c
# start with: --issl-threads 2 --issl-max-open-files 256


set -euo pipefail
log()  { local ts; ts="$(date +%F' '%T)"; printf '[%s] %s\n' "$ts" "$*" | tee -a "${LOGFILE:-/dev/null}" >&2; }
warn() { local ts; ts="$(date +%F' '%T)"; printf '[%s] WARN: %s\n' "$ts" "$*" | tee -a "${LOGFILE:-/dev/null}" >&2; }
die()  { local ts; ts="$(date +%F' '%T)"; printf '[%s] ERROR: %s\n' "$ts" "$*" | tee -a "${LOGFILE:-/dev/null}" >&2; exit 1; }
run()  { log "$*"; "$@"; }

# --- arg parse (overrides env defaults below) ---
VCFDIR="${VCFDIR:-}"
REF="${REF:-}"
OUT="${OUT:-}"
CHRMAP="${CHRMAP:-}"
THREADS="${THREADS:-8}"
BUILD_INDEX="${BUILD_INDEX:-0}" # 1=yes 0=no
HAPLOTYPES="${HAPLOTYPES:-both}" # both|1|2|iupac (consensus FASTAs always built for H1/H2)
FORCE="${FORCE:-0}"   # set FORCE=1 to rebuild even if outputs exist
PAM="${PAM:-NGG}"   # SpCas9 default; change if needed


TMP_ROOT="${TMP_ROOT:-}" # if empty, will use "$OUT/.tmp"
ISSL_THREADS="${ISSL_THREADS:-2}"
ISSL_MAX_OPEN_FILES="${ISSL_MAX_OPEN_FILES:-}" # default computed from ulimit
ISSL_RETRIES="${ISSL_RETRIES:-3}"

# --- contig utilities ---
ref_contigs_file() {
  # writes sorted unique REF contigs to stdout (requires .fai)
  local fai="${REF}.fai"
  [[ -s "$fai" ]] || run samtools faidx "$REF"
  cut -f1 "$fai" | sort -u
}

vcf_contigs_used() {          # CHROM values actually used (data lines)
  # $1 = vcf.gz
  bcftools view -H "$1" | cut -f1 | sort -u
}

vcf_contigs_header() {        # declared contigs in the header
  # $1 = vcf.gz
  bcftools view -h "$1" | awk -F'[=,>]' '/^##contig=/{print $3}' | sort -u
}

vcf_idxstats() {              # bcftools idxstats counts
  # $1 = vcf.gz
  bcftools idxstats "$1" | awk 'BEGIN{OFS="\t"} {print $1,$3}'
}

tabix_list() {                # sequences in the tabix index
  # $1 = vcf.gz
  tabix -l "$1" | sort -u
}

write_list() {                # util: write lines to file safely
  # $1 = path
  # stdin = content
  local f="$1"
  cat > "$f"
}

compare_lists() {
  # $1 list A, $2 list B
  # prints: INTERSECT/MISSING_A/MISSING_B sections
  local A="$1" B="$2"
  local tmpA tmpB
  tmpA="$(mktemp)"; tmpB="$(mktemp)"
  sort -u "$A" > "$tmpA"; sort -u "$B" > "$tmpB"
  printf "INTERSECT:\n"
  comm -12 "$tmpA" "$tmpB"
  printf "\nMISSING_IN_REF (present_in_VCF_not_in_REF):\n"
  comm -23 "$tmpA" "$tmpB"
  printf "\nMISSING_IN_VCF (present_in_REF_not_in_VCF):\n"
  comm -13 "$tmpA" "$tmpB"
  rm -f "$tmpA" "$tmpB"
}

contig_checks() {
  # $1 = vcf.gz ; $2 = stage tag (raw|renamed|norm)
  local vcf="$1" stage="$2" stem="$REPORT_DIR/${SAMPLE}.${stage}"

  # ensure index (needed for idxstats / tabix -l)
  [[ -s "${vcf}.tbi" || -s "${vcf}.csi" ]] || run tabix -p vcf "$vcf"

  # files we write
  local used="${stem}.contigs.used.txt"
  local hdr="${stem}.contigs.header.txt"
  local tbi="${stem}.tabix_list.txt"
  local idx="${stem}.idxstats.tsv"
  local ref="${REPORT_DIR}/${SAMPLE}.ref.contigs.txt"
  local cmp="${stem}.contig_compare.txt"

  # collect
  vcf_contigs_used   "$vcf" | write_list "$used"
  vcf_contigs_header "$vcf" | write_list "$hdr"
  tabix_list         "$vcf" | write_list "$tbi"
  vcf_idxstats       "$vcf" | write_list "$idx"
  # ref contigs (write once per sample)
  [[ -s "$ref" ]] || ref_contigs_file > "$ref"

  # compare "used" vs reference contigs
  compare_lists "$used" "$ref" > "$cmp"

  # summarise counts into global TSV (variants per contig)
  # idxstats is contig \t n_records
  awk -v s="$SAMPLE" -v st="$stage" 'BEGIN{OFS="\t"} {print s,st,$1,$2}' "$idx" >> "$SUMMARY_TSV"

  # log highlights
  local n_used n_hdr n_tbi n_idx n_overlap n_missing_ref
  n_used=$(wc -l < "$used" 2>/dev/null || echo 0)
  n_hdr=$(wc -l < "$hdr"  2>/dev/null || echo 0)
  n_tbi=$(wc -l < "$tbi"  2>/dev/null || echo 0)
  n_idx=$(wc -l < "$idx"  2>/dev/null || echo 0)
  n_overlap=$(awk '/^INTERSECT:/{p=1;next} /^$/{next} /^MISSING_/{p=0} p{c++} END{print c+0}' "$cmp")
  n_missing_ref=$(awk '/^MISSING_IN_REF/{p=1;next} /^$/{next} /^MISSING_IN_VCF/{p=0} p{c++} END{print c+0}' "$cmp")

  log "VCF[$stage]: used=${n_used}, header=${n_hdr}, tabix=${n_tbi}, idxstats_rows=${n_idx}, overlap_with_REF=${n_overlap}, missing_in_REF=${n_missing_ref}"
  if (( n_overlap == 0 && n_used > 0 )); then
    warn "No contig name overlap with REF at stage=${stage}. Check naming (e.g., chr1 vs 1 vs NC_*). See: $cmp"
  elif (( n_missing_ref > 0 )); then
    warn "Some VCF contigs not in REF at stage=${stage}. See details: $cmp"
  fi
}

sanitize_fa() {
  # $1=in fasta, $2=out fasta
  awk 'BEGIN{IGNORECASE=1}
    /^>/ {sub(/^>/,""); split($0,a,/[ \t]/); print ">" a[1]; next}
    { u=toupper($0); gsub(/[^ACGTN]/,"N",u); print u }' "$1" > "$2"
  samtools faidx "$2"
}

detect_issl_mode() {
  # Sets global ISSL_MODE to "TXT" or "FASTA"
  local help
  help="$("$CREATE_ISSL" -h 2>&1 || true)"
  if grep -Eiq -- '-t[ ,=]|--txt|--offtargets' <<<"$help"; then
    ISSL_MODE="TXT"
  elif grep -Eiq -- '-g[ ,=]|--genome|FASTA' <<<"$help"; then
    ISSL_MODE="FASTA"
  else
    ISSL_MODE="TXT"  # best effort default
  fi
  log "Detected ISSL mode: $ISSL_MODE"
}

detect_issl_iface() {
  local h; h="$("$CREATE_ISSL" -h 2>&1 || true)"
  if grep -qi '\[offtargetSites\.txt\].*\[sequence length\].*\[slice width' <<<"$h"; then
    ISSL_IFACE="POS"   # positional args
  else
    ISSL_IFACE="FLAG"  # -t/-l/-w/-o flags
  fi
  log "Detected ISSL interface: $ISSL_IFACE"
}

usage(){ cat <<EOF
Usage: \
$0 --vcf-dir DIR --ref FASTA --out DIR [--rename TSV] [--threads N] \
[--build-indexes yes|no] [--haplotypes both|1|2|iupac] \
[--tmp DIR] [--issl-threads N] [--issl-max-open-files N] [--issl-retries N]
EOF
exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf-dir) VCFDIR="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --rename|--rename-map|--chrmap) CHRMAP="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --build-indexes)
      case "$2" in yes|1|true) BUILD_INDEX=1;; no|0|false) BUILD_INDEX=0;;
        *) echo "Bad --build-indexes $2"; exit 2;;
      esac
      shift 2;;
    --haplotypes) HAPLOTYPES="$2"; shift 2;;
    --tmp) TMP_ROOT="$2"; shift 2;;
    --issl-threads) ISSL_THREADS="$2"; shift 2;;
    --issl-max-open-files) ISSL_MAX_OPEN_FILES="$2"; shift 2;;
    --issl-retries) ISSL_RETRIES="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done
[[ -n "${VCFDIR:-}" && -n "${REF:-}" && -n "${OUT:-}" ]] || usage

# --- OUT sanity + normalise + report/log dirs (must be after arg parsing) ---
: "${OUT:?--out was empty; pass a real path}"

# Make OUT absolute (best effort)
OUT="$(readlink -m "$OUT" 2>/dev/null || realpath -m "$OUT" 2>/dev/null || printf "%s" "$OUT")"
OUT_PARENT="$(dirname "$OUT")"
[[ -d "$OUT_PARENT" && -w "$OUT_PARENT" ]] || die "Parent dir not writable: $OUT_PARENT"

mkdir -p "$OUT" || die "Could not create OUT: $OUT"

REPORT_DIR="${OUT%/}/reports"
LOG_DIR="${OUT%/}/logs"
mkdir -p "$REPORT_DIR" "$LOG_DIR" || die "Could not create $REPORT_DIR / $LOG_DIR"

SUMMARY_TSV="$REPORT_DIR/contig_summary.tsv"
[[ -s "$SUMMARY_TSV" ]] || printf "sample\tstage\tcontig\tvariants\n" > "$SUMMARY_TSV"

write_versions_report

mapfile -t VCF_FILES < <(find "$VCFDIR" -maxdepth 1 -type f \( -name '*.vcf' -o -name '*.vcf.gz' \) | sort)
(( ${#VCF_FILES[@]} > 0 )) || die "No VCFs found in: $VCFDIR"
log "Found ${#VCF_FILES[@]} VCF file(s) in $VCFDIR"

# ---- config (override via env) ----
CHRMAP="${CHRMAP:-}" # optional TSV "vcfContig<TAB>refContig"
THREADS="${THREADS:-8}"
BUILD_INDEX="${BUILD_INDEX:-0}"
GUIDE_LEN="${GUIDE_LEN:-20}"
SLICE_WIDTH="${SLICE_WIDTH:-8}"


need(){ command -v "$1" >/dev/null 2>&1 || { echo "Missing tool: $1" >&2; exit 2; }; }
for t in samtools bgzip tabix bcftools; do need "$t"; done
[[ "$BUILD_INDEX" != 1 ]] || { need bowtie2-build; need extractOfftargets; }
CREATE_ISSL="$(command -v createIsslIndex || command -v isslCreateIndex || true)"
[[ "$BUILD_INDEX" != 1 || -n "$CREATE_ISSL" ]] || { echo "Missing createIsslIndex/isslCreateIndex" >&2; exit 2; }

detect_issl_mode
detect_issl_iface


mkdir -p "$OUT"
if [[ -z "${TMP_ROOT:-}" ]]; then TMP_ROOT="$OUT/.tmp"; fi
mkdir -p "$TMP_ROOT"


# Prefer tmp under Linux FS to avoid WSL2 /mnt/c I/O flakiness
export TMPDIR="$TMP_ROOT"

# Defensive cap for max open files
_soft_limit=$(ulimit -Sn 2>/dev/null || echo 1024)
log "Soft open-files limit detected: ${_soft_limit}"

if [[ -z "${ISSL_MAX_OPEN_FILES:-}" ]]; then
  ISSL_MAX_OPEN_FILES=$(( _soft_limit - 512 ))
  (( ISSL_MAX_OPEN_FILES < 256 )) && ISSL_MAX_OPEN_FILES=256
fi
ISSL_THREADS="${ISSL_THREADS:-2}"
log "Capping ISSL open files to ${ISSL_MAX_OPEN_FILES}; ISSL threads = ${ISSL_THREADS}"

# --- version helpers ---
version_of() {
  # Prints: "<name>\t<version or status>"
  # Tries common flags; falls back to -h first line.
  local name="$1"
  local path; path="$(command -v "$name" 2>/dev/null || true)"
  if [[ -z "$path" ]]; then
    printf "%s\tNOT FOUND\n" "$name"
    return 0
  fi
  local ver=""
  ver="$("$name" --version 2>&1 | head -n1 || true)"
  [[ -n "$ver" ]] || ver="$("$name" -version 2>&1 | head -n1 || true)"
  [[ -n "$ver" ]] || ver="$("$name" -V 2>&1 | head -n1 || true)"
  [[ -n "$ver" ]] || ver="$("$name" -h 2>&1 | head -n1 || true)"
  printf "%s\t%s\n" "$name" "${ver:-unknown}"
}
write_versions_report() {
  # Writes tool versions to OUT/reports/tool_versions.txt and logs it.
  local out="$REPORT_DIR/tool_versions.txt"
  : > "$out"

  {
    echo "# Tool versions"
    printf "generated_at\t%s\n" "$(date -Is)"
    printf "host\t%s\n" "$(uname -a 2>/dev/null || echo unknown)"
    if command -v lsb_release >/dev/null 2>&1; then
      printf "distro\t%s\n" "$(lsb_release -ds)"
    fi
    echo

    # Core tools
    version_of bash
    version_of awk
    version_of bcftools
    version_of samtools
    version_of bgzip
    version_of tabix

    # Optional indexers (only if present/used)
    if [[ "$BUILD_INDEX" == 1 ]]; then
      version_of bowtie2-build
      version_of extractOfftargets
    fi

    # ISSL creator (if detected)
    if [[ -n "$CREATE_ISSL" ]]; then
      # $CREATE_ISSL may be 'createIsslIndex' or 'isslCreateIndex'
      local n; n="$(basename "$CREATE_ISSL")"
      # Try dedicated flags then help
      local v
      v="$("$CREATE_ISSL" --version 2>&1 | head -n1 || true)"
      [[ -n "$v" ]] || v="$("$CREATE_ISSL" -version 2>&1 | head -n1 || true)"
      [[ -n "$v" ]] || v="$("$CREATE_ISSL" -h 2>&1 | head -n1 || true)"
      printf "%s\t%s\n" "$n" "${v:-unknown}"
    fi

    echo
    echo "# Inputs"
    printf "REF\t%s\n" "$REF"
    if [[ -s "${REF}.fai" || -s "$REF" ]]; then
      # Include a quick checksum for reproducibility if available
      if command -v sha256sum >/dev/null 2>&1; then
        printf "REF_sha256\t%s\n" "$(sha256sum "$REF" | awk '{print $1}')"
      fi
      if [[ -s "${REF}.fai" ]]; then
        printf "REF_n_contigs\t%s\n" "$(cut -f1 "${REF}.fai" | wc -l)"
      fi
    fi
    if [[ -n "${CHRMAP:-}" ]]; then
      printf "CHRMAP\t%s\n" "$CHRMAP"
      if command -v sha256sum >/dev/null 2>&1 && [[ -f "$CHRMAP" ]]; then
        printf "CHRMAP_sha256\t%s\n" "$(sha256sum "$CHRMAP" | awk '{print $1}')"
      fi
    fi
  } | tee -a "$out" | while IFS= read -r line; do log "$line"; done
}

for VCF in "${VCF_FILES[@]}"; do
  BASENAME=$(basename "$VCF")
  SAMPLE=${BASENAME%.vcf.gz}
  SAMPLE=${SAMPLE%.vcf}

  LOGFILE="$LOG_DIR/${SAMPLE}.log"
  log "==== Begin ${SAMPLE} ===="
  if [[ -s "$REPORT_DIR/tool_versions.txt" ]]; then
    while IFS= read -r line; do log "$line"; done < "$REPORT_DIR/tool_versions.txt"
  fi

  # ----- Step 0: Ensure we have a bgzipped, indexed VCF -----
  if [[ "$VCF" =~ \.vcf$ ]]; then
    GZ="${VCF}.gz"
  else
    GZ="$VCF"
  fi
  if [[ $FORCE -eq 1 || ! ( -s "$GZ" && -s "${GZ}.tbi" ) ]]; then
    if [[ "$VCF" =~ \.vcf$ ]]; then
      run bgzip -@ "$THREADS" -c "$VCF" > "$GZ"
    fi
    run tabix -p vcf "$GZ"
  else
    log "Skip bgzip/tabix: $GZ exists"
  fi
  VCF="$GZ"
  contig_checks "$VCF" "raw"


  # ----- Step 1: Optional contig rename to match the FASTA -----
  if [[ -n "${CHRMAP:-}" ]]; then
    REN="$OUT/${SAMPLE}.renamed.vcf.gz"
    if [[ $FORCE -eq 1 || ! ( -s "$REN" && -s "$REN.tbi" ) ]]; then
      run bcftools annotate --rename-chrs "$CHRMAP" "$VCF" -Oz -o "$REN"
      run tabix -p vcf "$REN"
    else
      log "Skip rename: $REN exists"
    fi
    VCF="$REN"
    contig_checks "$VCF" "renamed"

  fi

  # ----- Step 2: Normalize/sort/index -----
  NORM="$OUT/${SAMPLE}.norm.vcf.gz"
  if [[ $FORCE -eq 1 || ! ( -s "$NORM" && -s "$NORM.tbi" ) ]]; then
    run bcftools norm --threads "$THREADS" -f "$REF" -m -both "$VCF" \
      | bcftools sort -T "$TMP_ROOT" -Oz -o "$NORM"
    run tabix -p vcf "$NORM"
  else
    log "Skip norm/sort: $NORM exists"
  fi
  contig_checks "$NORM" "norm"

  # ----- Step 3: Build consensus FASTA(s) -----
  case "$HAPLOTYPES" in
    iupac)
      CONS="$OUT/GRCh38_${SAMPLE}.consensus.fa"
      if [[ $FORCE -eq 1 || ! ( -s "$CONS" && -s "${CONS}.fai" ) ]]; then
        run bcftools consensus -f "$REF" --iupac-codes "$NORM" -o "$CONS"
        run samtools faidx "$CONS"
      else
        log "Skip consensus: $CONS exists"
      fi
      ;;

    1|2)
      CONS="$OUT/GRCh38_${SAMPLE}.H${HAPLOTYPES}.consensus.fa"
      if [[ $FORCE -eq 1 || ! ( -s "$CONS" && -s "${CONS}.fai" ) ]]; then
        run bcftools consensus -f "$REF" -H "$HAPLOTYPES" "$NORM" -o "$CONS"
        run samtools faidx "$CONS"
      else
        log "Skip consensus: $CONS exists"
      fi
      ;;

    both)
      for H in 1 2; do
        CONS="$OUT/GRCh38_${SAMPLE}.H${H}.consensus.fa"
        if [[ $FORCE -eq 1 || ! ( -s "$CONS" && -s "${CONS}.fai" ) ]]; then
          run bcftools consensus -f "$REF" -H "$H" "$NORM" -o "$CONS"
          run samtools faidx "$CONS"
        else
          log "Skip consensus (H$H): $CONS exists"
        fi

        # ----- Optional: per-haplotype indexes -----
        if [[ "$BUILD_INDEX" == 1 ]]; then
          # Bowtie2
          if [[ $FORCE -eq 1 || ! ( -s "$OUT/GRCh38_${SAMPLE}.H${H}.1.bt2" && -s "$OUT/GRCh38_${SAMPLE}.H${H}.rev.2.bt2" ) ]]; then
            run bowtie2-build --threads "$THREADS" "$CONS" "$OUT/GRCh38_${SAMPLE}.H${H}"
          else
            log "Skip Bowtie2 (H$H): index exists"
          fi

          # ISSL via extractOfftargets -> createIsslIndex
          if [[ -n "$CREATE_ISSL" ]]; then
            SAN="${CONS%.fa}.san.fa"
            if [[ $FORCE -eq 1 || ! ( -s "$SAN" && -s "${SAN}.fai" ) ]]; then
              run sanitize_fa "$CONS" "$SAN"
            else
              log "Skip sanitize (H$H): $SAN exists"
            fi

            OFF="$OUT/GRCh38_${SAMPLE}.H${H}_offtargets.txt"
            if [[ $FORCE -eq 1 || ! -s "$OFF" ]]; then
              run extractOfftargets --threads "$ISSL_THREADS" --maxOpenFiles "$ISSL_MAX_OPEN_FILES" "$OFF" "$SAN"
            else
              log "Skip extractOfftargets (H$H): $OFF exists"
            fi

            ISSL_OUT="$OUT/GRCh38_${SAMPLE}${H:+.H${H}}.issl"
            if [[ $FORCE -eq 1 || ! -s "$ISSL_OUT" ]]; then
            ( ulimit -n "${ISSL_MAX_OPEN_FILES}";
                if [[ "$ISSL_IFACE" == "POS" ]]; then
                # Positional-args build: [txt] [guide_len] [slice_width_bits] [out]
                run "$CREATE_ISSL" "$OFF" "$GUIDE_LEN" "$SLICE_WIDTH" "$ISSL_OUT"
                else
                # Flag-args build
                run "$CREATE_ISSL" -t "$OFF" -l "$GUIDE_LEN" -w "$SLICE_WIDTH" -o "$ISSL_OUT"
            fi
  )
else
  log "Skip ISSL${H:+ (H$H)}: $ISSL_OUT exists"
fi
          fi
        fi
      done
      log "Finished ${SAMPLE} (both haplotypes)"
      log "==== End ${SAMPLE} ===="
      continue
      ;;

    *)
      warn "Unknown --haplotypes '${HAPLOTYPES}', defaulting to iupac"
      CONS="$OUT/GRCh38_${SAMPLE}.consensus.fa"
      if [[ $FORCE -eq 1 || ! ( -s "$CONS" && -s "${CONS}.fai" ) ]]; then
        run bcftools consensus -f "$REF" --iupac-codes "$NORM" -o "$CONS"
        run samtools faidx "$CONS"
      else
        log "Skip consensus: $CONS exists"
      fi
      ;;
  esac

  # ----- Step 4: Global indexes (non-both modes) -----
  if [[ "$BUILD_INDEX" == 1 ]]; then
    # Bowtie2
    if [[ $FORCE -eq 1 || ! ( -s "$OUT/GRCh38_${SAMPLE}.1.bt2" && -s "$OUT/GRCh38_${SAMPLE}.rev.2.bt2" ) ]]; then
      run bowtie2-build --threads "$THREADS" "$CONS" "$OUT/GRCh38_${SAMPLE}"
    else
      log "Skip Bowtie2: index exists for $SAMPLE"
    fi

    # ISSL
    if [[ -n "$CREATE_ISSL" ]]; then
      SAN="${CONS%.fa}.san.fa"
      if [[ $FORCE -eq 1 || ! ( -s "$SAN" && -s "${SAN}.fai" ) ]]; then
        run sanitize_fa "$CONS" "$SAN"
      else
        log "Skip sanitize: $SAN exists"
      fi

      OFF="$OUT/GRCh38_${SAMPLE}_offtargets.txt"
      if [[ $FORCE -eq 1 || ! -s "$OFF" ]]; then
        run extractOfftargets --threads "$ISSL_THREADS" --maxOpenFiles "$ISSL_MAX_OPEN_FILES" "$OFF" "$SAN"
      else
        log "Skip extractOfftargets: $OFF exists"
      fi

      ISSL_OUT="$OUT/GRCh38_${SAMPLE}${H:+.H${H}}.issl"
        if [[ $FORCE -eq 1 || ! -s "$ISSL_OUT" ]]; then
        ( ulimit -n "${ISSL_MAX_OPEN_FILES}";
            if [[ "$ISSL_IFACE" == "POS" ]]; then
            # Positional-args build: [txt] [guide_len] [slice_width_bits] [out]
            run "$CREATE_ISSL" "$OFF" "$GUIDE_LEN" "$SLICE_WIDTH" "$ISSL_OUT"
            else
            # Flag-args build
            run "$CREATE_ISSL" -t "$OFF" -l "$GUIDE_LEN" -w "$SLICE_WIDTH" -o "$ISSL_OUT"
            fi
        )
else
  log "Skip ISSL${H:+ (H$H)}: $ISSL_OUT exists"
fi
    fi
  fi

  log "Finished ${SAMPLE}"
  log "==== End ${SAMPLE} ===="
done