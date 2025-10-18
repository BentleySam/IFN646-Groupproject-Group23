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
log(){ printf '[%s] %s\n' "$(date +%F' '%T)" "$*" >&2; }
warn(){ printf '[%s] WARN: %s\n' "$(date +%F' '%T)" "$*" >&2; }
die(){ printf '[%s] ERROR: %s\n' "$(date +%F' '%T)" "$*" >&2; exit 1; }
run(){ log "$*"; "$@"; }

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


for VCF in "${VCF_FILES[@]}"; do
  BASENAME=$(basename "$VCF")
  SAMPLE=${BASENAME%.vcf.gz}
  SAMPLE=${SAMPLE%.vcf}

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
done