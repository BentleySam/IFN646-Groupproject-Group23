#!/usr/bin/env bash
# Orchestrator for Crackling + VCF consensus pipeline
# Wires together:
#   VCF_To_Consensus.sh
#   Run_Crackling.sh
#   Run_Custom_Crackling.sh
#   Compare_Crackling_Runs.sh
#   Check_Guides_In_Fastas.sh
#
# Keep this thin; all heavy lifting lives in the single-purpose scripts.
set -euo pipefail

# --- ensure conda env 'crackling' is active ---
if command -v conda >/dev/null 2>&1; then
  # Works even in non-interactive shells
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate crackling || {
    echo "ERROR: could not activate 'crackling' env" >&2; exit 2;
  }
else
  # Fallback if 'conda' isn’t on PATH
  if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate crackling || { echo "ERROR: conda activate failed" >&2; exit 2; }
  else
    echo "ERROR: conda not found. Install miniconda and/or add it to PATH." >&2
    exit 2
  fi
fi

# ------------------ defaults ------------------
REF="${REF:-$HOME/genomes/GRCh38.fa}"
GENOME_DIR="${GENOME_DIR:-$HOME/genomes}"
BASE_BOWTIE_PREFIX="${BASE_BOWTIE_PREFIX:-GRCh38}"
BASE_ISSL="${BASE_ISSL:-GRCh38_offtargets.issl}"

VCF_DIR=""
PERSONAL_OUT="${PERSONAL_OUT:-$GENOME_DIR/personalised}"

INPUT_FASTA=""                # target region(s)
CR_OUTDIR="${CR_OUTDIR:-$HOME/cracklingOutput}"
JOBTAG="${JOBTAG:-CracklingRun}"

THREADS="${THREADS:-$(nproc || echo 8)}"
HAPLOTYPES="${HAPLOTYPES:-both}"    # both|1|2|iupac
BUILD_INDEXES="${BUILD_INDEXES:-yes}"  # yes|no

OFFTARGET_MODE="${OFFTARGET_MODE:-avg}"       # avg|and|mit|cfd|or
OFFTARGET_THRESHOLD="${OFFTARGET_THRESHOLD:-65}"

DO_GUIDE_CHECK="${DO_GUIDE_CHECK:-yes}"       # yes|no

# Optional contig rename map (VCF -> FASTA), TSV with 2 columns: old<TAB>new
RENAME_MAP="${RENAME_MAP:-}"

# ------------------ usage ------------------
usage() {
  cat <<EOF
Usage: $0 --vcf-dir DIR --input FASTA [options]

Required:
  --vcf-dir DIR           Directory containing *.vcf or *.vcf.gz (can be nested)
  --input FASTA           Target region FASTA for Crackling (one or many regions)

Optional (paths):
  --ref FASTA             Reference FASTA (default: $REF)
  --genome-dir DIR        Genome assets dir (default: $GENOME_DIR)
  --personal-out DIR      Output root for consensus builds (default: $PERSONAL_OUT)
  --cr-outdir DIR         Crackling outputs root (default: $CR_OUTDIR)

Optional (base genome):
  --base-bowtie-prefix P  Bowtie2 basename for base genome (default: $BASE_BOWTIE_PREFIX)
  --base-issl FILE        ISSL file for base genome (default: $BASE_ISSL)

Optional (VCF to consensus):
  --haplotypes MODE       both|1|2|iupac (default: $HAPLOTYPES)
  --build-indexes yes|no  Build Bowtie2+ISSL per consensus (default: $BUILD_INDEXES)
  --rename-map TSV        Contig rename map (old\\tnew), if VCF != FASTA names

Optional (Crackling):
  --offtarget-mode MODE   avg|and|mit|cfd|or (default: $OFFTARGET_MODE)
  --offtarget-threshold N (default: $OFFTARGET_THRESHOLD)
  --threads N             (default: $THREADS)
  --jobtag NAME           (default: $JOBTAG)

Optional (post-processing):
  --guide-check yes|no    Check guide presence across FASTAs (default: $DO_GUIDE_CHECK)

Notes:
- This script assumes the helper scripts are in PATH:
  VCF_To_Consensus.sh, Run_Crackling.sh, Run_Custom_Crackling.sh,
  Compare_Crackling_Runs.sh, Check_Guides_In_Fastas.sh
EOF
  exit 1
}

# ------------------ arg parse ------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf-dir) VCF_DIR="$2"; shift 2;;
    --input) INPUT_FASTA="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --genome-dir) GENOME_DIR="$2"; shift 2;;
    --personal-out) PERSONAL_OUT="$2"; shift 2;;
    --cr-outdir) CR_OUTDIR="$2"; shift 2;;
    --base-bowtie-prefix) BASE_BOWTIE_PREFIX="$2"; shift 2;;
    --base-issl) BASE_ISSL="$2"; shift 2;;
    --haplotypes) HAPLOTYPES="$2"; shift 2;;
    --build-indexes) BUILD_INDEXES="$2"; shift 2;;
    --rename-map) RENAME_MAP="$2"; shift 2;;
    --offtarget-mode) OFFTARGET_MODE="$2"; shift 2;;
    --offtarget-threshold) OFFTARGET_THRESHOLD="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --jobtag) JOBTAG="$2"; shift 2;;
    --guide-check) DO_GUIDE_CHECK="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done
[[ -n "${VCF_DIR:-}" && -n "${INPUT_FASTA:-}" ]] || usage

# ------------------ helpers ------------------
die(){ echo "ERROR: $*" >&2; exit 2; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing tool: $1"; }
abspath(){
  python3 - "$1" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
}

latest_run_dir() {
  # Picks the most recent run folder under CR_OUTDIR matching "<jobtag>_*"
  local base="$1" tag="$2"
  ls -1dt "${base}/${tag}_"* 2>/dev/null | head -n1 || true
}

timestamp(){ date +%Y%m%d_%H%M%S; }

# ------------------ preflight ------------------
for s in VCF_To_Consensus.sh Run_Crackling.sh Run_Custom_Crackling.sh Compare_Crackling_Runs.sh; do
  need "$s"
done
[[ "${DO_GUIDE_CHECK}" == "no" ]] || need Check_Guides_In_Fastas.sh

[[ -s "$REF" ]] || die "Reference FASTA missing: $REF"
mkdir -p "$CR_OUTDIR" "$PERSONAL_OUT"

REF="$(abspath "$REF")"
VCF_DIR="$(abspath "$VCF_DIR")"
GENOME_DIR="$(abspath "$GENOME_DIR")"
PERSONAL_OUT="$(abspath "$PERSONAL_OUT")"
CR_OUTDIR="$(abspath "$CR_OUTDIR")"
INPUT_FASTA="$(abspath "$INPUT_FASTA")"
[[ -n "${RENAME_MAP:-}" ]] && RENAME_MAP="$(abspath "$RENAME_MAP")"

# ------------------ logging ------------------
RUNSTAMP="$(timestamp)"
MASTER_LOG="${CR_OUTDIR}/${JOBTAG}_${RUNSTAMP}_MASTER.log"
exec > >(tee -a "$MASTER_LOG") 2>&1
set -x
echo "# $(date) :: starting Run_All_Crackling_Pipeline.sh"
echo "# REF=$REF"
echo "# VCF_DIR=$VCF_DIR"
echo "# PERSONAL_OUT=$PERSONAL_OUT"
echo "# CR_OUTDIR=$CR_OUTDIR"
echo "# HAPLOTYPES=$HAPLOTYPES BUILD_INDEXES=$BUILD_INDEXES RENAME_MAP=${RENAME_MAP:-none}"
echo "# OFFTARGET_MODE=$OFFTARGET_MODE OFFTARGET_THRESHOLD=$OFFTARGET_THRESHOLD THREADS=$THREADS"
echo "# BASE_BOWTIE_PREFIX=$BASE_BOWTIE_PREFIX BASE_ISSL=$BASE_ISSL"
echo

# ------------------ 1) Run base genome once ------------------
Run_Crackling.sh \
  --input  "$INPUT_FASTA" \
  --outdir "$CR_OUTDIR" \
  --genome-dir "$GENOME_DIR" \
  --bowtie-prefix "$BASE_BOWTIE_PREFIX" \
  --issl "$BASE_ISSL" \
  --threads "$THREADS" \
  --offtarget-mode "$OFFTARGET_MODE" \
  --offtarget-threshold "$OFFTARGET_THRESHOLD" \
  --jobtag "${JOBTAG}_BASE"

BASE_RUN_DIR="$(latest_run_dir "$CR_OUTDIR" "${JOBTAG}_BASE")"
[[ -n "$BASE_RUN_DIR" ]] || die "Could not find base run dir"
BASE_SUMMARY="${BASE_RUN_DIR}/${JOBTAG}_BASE_final_guides_summary.csv"
BASE_GUIDES="${BASE_RUN_DIR}/${JOBTAG}_BASE-guides.txt"
[[ -s "$BASE_SUMMARY" && -s "$BASE_GUIDES" ]] || die "Base outputs missing"

# ------------------ 2) Build consensus FASTAs from VCFs ------------------
VCF_ARGS=( --vcf-dir "$VCF_DIR" --ref "$REF" --out "$PERSONAL_OUT" --haplotypes "$HAPLOTYPES" --build-indexes "$BUILD_INDEXES" --threads "$THREADS" )
[[ -n "${RENAME_MAP:-}" ]] && VCF_ARGS+=( --rename "$RENAME_MAP" )
VCF_To_Consensus.sh "${VCF_ARGS[@]}"

# Discover per-sample directories created by VCF_To_Consensus.sh
mapfile -t SAMPLE_DIRS < <(find "$PERSONAL_OUT" -maxdepth 1 -type d -not -path "$PERSONAL_OUT" | sort)
[[ ${#SAMPLE_DIRS[@]} -gt 0 ]] || die "No per-sample directories found in $PERSONAL_OUT"

# ------------------ 3) Run Crackling per sample/haplotype ------------------
COMPARE_ROOT="${CR_OUTDIR}/Comparisons/${JOBTAG}_${RUNSTAMP}"
mkdir -p "$COMPARE_ROOT"

for SDIR in "${SAMPLE_DIRS[@]}"; do
  SNAME="$(basename "$SDIR")"

  # Detect bowtie/issl basenames (convention from VCF_To_Consensus.sh)
  # Expect something like:
  #   ${SDIR}/GRCh38_${SNAME}.H1.fa  -> bowtie prefix: ${SDIR}/GRCh38_${SNAME}.H1.consensus
  #   ${SDIR}/${SNAME}_offtargets.issl
  # IUPAC mode might produce GRCh38_${SNAME}.IUPAC.fa
  if [[ "$HAPLOTYPES" == "both" || "$HAPLOTYPES" == "1" || "$HAPLOTYPES" == "2" ]]; then
    HAPS=( H1 H2 )
    [[ "$HAPLOTYPES" == "1" ]] && HAPS=( H1 )
    [[ "$HAPLOTYPES" == "2" ]] && HAPS=( H2 )
    for H in "${HAPS[@]}"; do
      BOWPFX="GRCh38_${SNAME}.${H}.consensus"
      ISSL="${SNAME}_offtargets.issl"

      Run_Custom_Crackling.sh \
        --input  "$INPUT_FASTA" \
        --outdir "$CR_OUTDIR" \
        --personal-dir "$SDIR" \
        --bowtie-prefix "$BOWPFX" \
        --issl "$ISSL" \
        --threads "$THREADS" \
        --offtarget-mode "$OFFTARGET_MODE" \
        --offtarget-threshold "$OFFTARGET_THRESHOLD" \
        --jobtag "${JOBTAG}_${SNAME}_${H}"

      ALT_RUN_DIR="$(latest_run_dir "$CR_OUTDIR" "${JOBTAG}_${SNAME}_${H}")"
      [[ -n "$ALT_RUN_DIR" ]] || die "Missing alt run dir for ${SNAME} ${H}"

      # Compare base vs personalized
      OUTCMP="${COMPARE_ROOT}/${SNAME}_${H}"
      mkdir -p "$OUTCMP"
      Compare_Crackling_Runs.sh \
        --base-run "$BASE_RUN_DIR" \
        --alt-run  "$ALT_RUN_DIR" \
        --outdir   "$OUTCMP"

      # Optional guide presence check across FASTAs (base + personalized consensus)
      if [[ "$DO_GUIDE_CHECK" == "yes" ]]; then
        Check_Guides_In_Fastas.sh \
          --guides "${BASE_SUMMARY}" \
          --fastas "$REF" "${SDIR}/GRCh38_${SNAME}.${H}.fa" \
          --out    "${OUTCMP}/GuidePresence_${SNAME}_${H}.csv"
      fi
    done

  else
    # IUPAC or other single build
    BOWPFX="GRCh38_${SNAME}.IUPAC.consensus"
    ISSL="${SNAME}_offtargets.issl"

    Run_Custom_Crackling.sh \
      --input  "$INPUT_FASTA" \
      --outdir "$CR_OUTDIR" \
      --personal-dir "$SDIR" \
      --bowtie-prefix "$BOWPFX" \
      --issl "$ISSL" \
      --threads "$THREADS" \
      --offtarget-mode "$OFFTARGET_MODE" \
      --offtarget-threshold "$OFFTARGET_THRESHOLD" \
      --jobtag "${JOBTAG}_${SNAME}_IUPAC"

    ALT_RUN_DIR="$(latest_run_dir "$CR_OUTDIR" "${JOBTAG}_${SNAME}_IUPAC")"
    [[ -n "$ALT_RUN_DIR" ]] || die "Missing alt run dir for ${SNAME} IUPAC"

    OUTCMP="${COMPARE_ROOT}/${SNAME}_IUPAC"
    mkdir -p "$OUTCMP"
    Compare_Crackling_Runs.sh \
      --base-run "$BASE_RUN_DIR" \
      --alt-run  "$ALT_RUN_DIR" \
      --outdir   "$OUTCMP"

    if [[ "$DO_GUIDE_CHECK" == "yes" ]]; then
      Check_Guides_In_Fastas.sh \
        --guides "${BASE_SUMMARY}" \
        --fastas "$REF" "${SDIR}/GRCh38_${SNAME}.IUPAC.fa" \
        --out    "${OUTCMP}/GuidePresence_${SNAME}_IUPAC.csv"
    fi
  fi
done

set +x
echo
echo "=== ALL DONE ==="
echo "Master log        : $MASTER_LOG"
echo "Base run          : $BASE_RUN_DIR"
echo "Comparisons root  : $COMPARE_ROOT"
echo "Crackling outdir  : $CR_OUTDIR"
