#!/usr/bin/env bash
set -euo pipefail

# --- arg parse (overrides env defaults below) ---
VCFDIR="${VCFDIR:-}"
REF="${REF:-}"
OUT="${OUT:-}"
CHRMAP="${CHRMAP:-}"
THREADS="${THREADS:-8}"
BUILD_INDEX="${BUILD_INDEX:-0}"   # 1=yes 0=no
HAPLOTYPES="${HAPLOTYPES:-both}"  # both|1|2|iupac (not fully used yet)

usage(){ cat <<EOF
Usage: $0 --vcf-dir DIR --ref FASTA --out DIR [--rename TSV] [--threads N]
            [--build-indexes yes|no] [--haplotypes both|1|2|iupac]
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
      case "$2" in yes|1|true) BUILD_INDEX=1;; no|0|false) BUILD_INDEX=0;; *) echo "Bad --build-indexes $2"; exit 2;; esac
      shift 2;;
    --haplotypes) HAPLOTYPES="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

[[ -n "${VCFDIR:-}" && -n "${REF:-}" && -n "${OUT:-}" ]] || usage

# ---- config (override via env) ----
CHRMAP="${CHRMAP:-}"            # optional TSV "vcfContig<TAB>refContig"
THREADS="${THREADS:-8}"
BUILD_INDEX="${BUILD_INDEX:-0}" # 1 = also build Bowtie2+ISSL per consensus
GUIDE_LEN="${GUIDE_LEN:-20}"
SLICE_WIDTH="${SLICE_WIDTH:-8}"

need(){ command -v "$1" >/dev/null 2>&1 || { echo "Missing tool: $1" >&2; exit 2; }; }
for t in samtools bgzip tabix bcftools; do need "$t"; done
[[ "$BUILD_INDEX" != 1 ]] || { need bowtie2-build; need extractOfftargets; }

CREATE_ISSL="$(command -v createIsslIndex || command -v isslCreateIndex || true)"
[[ "$BUILD_INDEX" != 1 || -n "$CREATE_ISSL" ]] || { echo "Missing createIsslIndex/isslCreateIndex" >&2; exit 2; }

mkdir -p "$OUT"
samtools faidx "$REF"

log(){ echo "[$(date +%F\ %T)] $*"; }

process_vcf() {
  local vcf="$1" stem gz ren norm
  stem="$(basename "$vcf")"; stem="${stem%.vcf.gz}"; stem="${stem%.vcf}"
  log "VCF: $stem"

  # compress/index
  gz="$OUT/${stem}.vcf.gz"
  if [[ "$vcf" == *.vcf.gz ]]; then cp -f "$vcf" "$gz"; else bgzip -c "$vcf" > "$gz"; fi
  tabix -f -p vcf "$gz"

  # rename contigs if map provided
  if [[ -n "$CHRMAP" ]]; then
    ren="$OUT/${stem}.renamed.vcf.gz"
    bcftools annotate --rename-chrs "$CHRMAP" -Oz -o "$ren" "$gz"
    tabix -f -p vcf "$ren"
  else
    ren="$gz"
  fi

  # normalize/sort
  norm="$OUT/${stem}.norm.vcf.gz"
  bcftools norm -f "$REF" -m -both "$ren" | bcftools sort -Oz -o "$norm"
  tabix -f -p vcf "$norm"

  # samples
  mapfile -t samples < <(bcftools query -l "$norm" || true)
  ((${#samples[@]})) || samples=("SAMPLE")

  for s in "${samples[@]}"; do
    for H in 1 2; do
      local fa="$OUT/$(basename "$REF" .fa)_${stem}.${s}.H${H}.fa"
      if [[ "$s" == "SAMPLE" ]]; then
        bcftools consensus -f "$REF" -H "$H" "$norm" > "$fa"
      else
        bcftools consensus -f "$REF" -s "$s" -H "$H" "$norm" > "$fa"
      fi
      log "FASTA: $(basename "$fa") bytes=$(stat -c%s "$fa")"

      if [[ "$BUILD_INDEX" == 1 ]]; then
        local idx="$OUT/${stem}.${s}.H${H}.idx"
        mkdir -p "$idx"
        if ! [[ -s "$idx/consensus.1.bt2" ]]; then
          log "bowtie2-build -> $idx/consensus"
          bowtie2-build --threads "$THREADS" "$fa" "$idx/consensus"
        fi
        if ! [[ -s "$idx/offtargets.issl" ]]; then
          log "extractOfftargets + ISSL -> $idx/offtargets.issl"
          extractOfftargets "$idx/offtargets.txt" "$fa"
          "$CREATE_ISSL" -t "$idx/offtargets.txt" -l "$GUIDE_LEN" -w "$SLICE_WIDTH" -o "$idx/offtargets.issl"
        fi
      fi
    done
  done
}

shopt -s nullglob
mapfile -t VCFs < <(find "$VCFDIR" -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) | sort)
if [[ ${#VCFs[@]} -eq 0 ]]; then
  echo "No VCFs in \"$VCFDIR\"" >&2
  exit 2
fi

log "REF=$REF  VCFDIR=$VCFDIR  OUT=$OUT  BUILD_INDEX=$BUILD_INDEX"
for v in "${VCFs[@]}"; do process_vcf "$v"; done
log "Done."
