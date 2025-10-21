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

# ---------- hashing helpers ----------
sha_file() {  # fast and robust: path + size + mtime + sha256 of first 1MB
  local p="$1"
  [[ -f "$p" ]] || { echo "MISSING:$p"; return; }
  local size mtime headhash
  size=$(stat -c%s "$p" 2>/dev/null || echo 0)
  mtime=$(stat -c%Y "$p" 2>/dev/null || echo 0)
  headhash=$(head -c 1048576 "$p" | sha256sum | awk '{print $1}')
  printf "%s|%s|%s|%s\n" "$p" "$size" "$mtime" "$headhash"
}

sha_list_bt2() {  # bowtie2 index family (6 files)
  local prefix="$1" sfx
  for sfx in 1 2 3 4 "rev.1" "rev.2"; do
    sha_file "${prefix}.${sfx}.bt2"
  done
}

sha_text() { printf "%s" "$1" | sha256sum | awk '{print $1}'; }

tool_ver() {
  # best-effort version string, stable enough for fingerprinting
  case "$1" in
    bowtie2) bowtie2 --version | head -n1 ;;
    RNAfold) RNAfold --version 2>/dev/null | head -n1 || echo RNAfold ;;
    Crackling) Crackling --help 2>/dev/null | head -n1 || echo Crackling ;;
    isslScoreOfftargets) command -v isslScoreOfftargets || echo isslScoreOfftargets ;;
    *) command -v "$1" || echo "$1" ;;
  esac
}

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

# ------------------ 1) Run base genome once (hash-based reuse) ------------------
MANIDIR="${CR_OUTDIR}/.manifests"
mkdir -p "$MANIDIR"

# Build a recipe string that is independent of timestamps/paths that vary per run
BASE_RECIPE=$(
  {
    echo "INPUT_FASTA_SHA=$(sha_file "$INPUT_FASTA")"
    echo "BOWTIE_BT2_SHA="
    sha_list_bt2 "${GENOME_DIR}/${BASE_BOWTIE_PREFIX}"
    echo "ISSL_SHA=$(sha_file "${GENOME_DIR}/${BASE_ISSL}")"
    echo "OFFTARGET_MODE=$OFFTARGET_MODE"
    echo "OFFTARGET_THRESHOLD=$OFFTARGET_THRESHOLD"
    echo "THREADS=$THREADS"
    echo "VERSIONS:"
    echo "  $(tool_ver bowtie2)"
    echo "  $(tool_ver RNAfold)"
    echo "  $(tool_ver Crackling)"
    echo "  $(tool_ver isslScoreOfftargets)"
  } | sed 's/[[:space:]]\+$//'
)

BASE_HASH="$(sha_text "$BASE_RECIPE")"
BASE_MANIFEST="${MANIDIR}/${JOBTAG}_BASE_${BASE_HASH}.manifest"

reuse_base=false
if [[ -s "$BASE_MANIFEST" ]]; then
  # shellcheck disable=SC1090
  source "$BASE_MANIFEST" || true
  if [[ -n "${RUN_DIR:-}" && -s "${RUN_DIR}/${JOBTAG}_BASE_final_guides_summary.csv" ]]; then
    reuse_base=true
    BASE_RUN_DIR="$RUN_DIR"
    echo "# Reusing base run via manifest: $BASE_MANIFEST"
  else
    echo "# Stale manifest found; will rebuild base run."
    reuse_base=false
  fi
fi

if ! $reuse_base; then
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
  [[ -n "$BASE_RUN_DIR" ]] || die "Could not find base run dir after run"

  # Write manifest and stable symlink
  {
    echo "HASH=$BASE_HASH"
    echo "RUN_DIR=$BASE_RUN_DIR"
    echo "RECIPE_EPOCH=$(date -u +%s)"
  } > "$BASE_MANIFEST"

  ln -sfn "$BASE_RUN_DIR" "${CR_OUTDIR}/${JOBTAG}_BASE_${BASE_HASH}"
  ln -sfn "$BASE_RUN_DIR" "${CR_OUTDIR}/${JOBTAG}_BASE_latest"
fi

BASE_SUMMARY="${BASE_RUN_DIR}/${JOBTAG}_BASE_final_guides_summary.csv"
BASE_GUIDES="${BASE_RUN_DIR}/${JOBTAG}_BASE-guides.txt"
[[ -s "$BASE_SUMMARY" && -s "$BASE_GUIDES" ]] || die "Base outputs missing (summary/guides)"

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
        ALT_TAG="${JOBTAG}_${SNAME}_${H}"
        ALT_FASTA="${SDIR}/GRCh38_${SNAME}.${H}.fa"
        [[ -s "$ALT_FASTA" ]] || die "Missing consensus FASTA: $ALT_FASTA"

        # --- per-alt hash/manifest reuse ---
        ALT_MANIDIR="${CR_OUTDIR}/.manifests"
        mkdir -p "$ALT_MANIDIR"

        ALT_RECIPE=$(
        {
            echo "INPUT_FASTA_SHA=$(sha_file "$INPUT_FASTA")"
            echo "ALT_FASTA_SHA=$(sha_file "$ALT_FASTA")"
            echo "BOWTIE_BT2_SHA="
            sha_list_bt2 "${SDIR}/${BOWPFX}"
            echo "ISSL_SHA=$(sha_file "${SDIR}/${ISSL}")"
            echo "OFFTARGET_MODE=$OFFTARGET_MODE"
            echo "OFFTARGET_THRESHOLD=$OFFTARGET_THRESHOLD"
            echo "THREADS=$THREADS"
            echo "VERSIONS:"
            echo "  $(tool_ver bowtie2)"
            echo "  $(tool_ver RNAfold)"
            echo "  $(tool_ver Crackling)"
            echo "  $(tool_ver isslScoreOfftargets)"
        } | sed 's/[[:space:]]\+$//'
        )
        ALT_HASH="$(sha_text "$ALT_RECIPE")"
        ALT_MANIFEST="${ALT_MANIDIR}/${ALT_TAG}_${ALT_HASH}.manifest"

        reuse_alt=false
        if [[ -s "$ALT_MANIFEST" ]]; then
        # shellcheck disable=SC1090
        source "$ALT_MANIFEST" || true
        if [[ -n "${RUN_DIR:-}" && -s "${RUN_DIR}/${ALT_TAG}_final_guides_summary.csv" ]]; then
            reuse_alt=true
            ALT_RUN_DIR="$RUN_DIR"
            echo "# Reusing alt run via manifest: $ALT_MANIFEST"
        else
            echo "# Stale alt manifest; will rebuild: $ALT_TAG"
            reuse_alt=false
        fi
        fi

        if ! $reuse_alt; then
        Run_Custom_Crackling.sh \
            --input  "$INPUT_FASTA" \
            --outdir "$CR_OUTDIR" \
            --personal-dir "$SDIR" \
            --bowtie-prefix "$BOWPFX" \
            --issl "$ISSL" \
            --threads "$THREADS" \
            --offtarget-mode "$OFFTARGET_MODE" \
            --offtarget-threshold "$OFFTARGET_THRESHOLD" \
            --jobtag "$ALT_TAG"

        ALT_RUN_DIR="$(latest_run_dir "$CR_OUTDIR" "$ALT_TAG")"
        [[ -n "$ALT_RUN_DIR" ]] || die "Missing alt run dir for $ALT_TAG"

        {
            echo "HASH=$ALT_HASH"
            echo "RUN_DIR=$ALT_RUN_DIR"
            echo "RECIPE_EPOCH=$(date -u +%s)"
        } > "$ALT_MANIFEST"

        ln -sfn "$ALT_RUN_DIR" "${CR_OUTDIR}/${ALT_TAG}_${ALT_HASH}"
        ln -sfn "$ALT_RUN_DIR" "${CR_OUTDIR}/${ALT_TAG}_latest"
        fi

        # --- compare + optional guide presence ---
        OUTCMP="${COMPARE_ROOT}/${SNAME}_${H}"
        mkdir -p "$OUTCMP"
        Compare_Crackling_Runs.sh \
        --base-run "$BASE_RUN_DIR" \
        --alt-run  "$ALT_RUN_DIR" \
        --outdir   "$OUTCMP"

        if [[ "$DO_GUIDE_CHECK" == "yes" ]]; then
        Check_Guides_In_Fastas.sh \
            --guides "${BASE_SUMMARY}" \
            --fastas "$REF" "$ALT_FASTA" \
            --out    "${OUTCMP}/GuidePresence_${SNAME}_${H}.csv"
        fi
    done

  else
    # IUPAC or other single build
    BOWPFX="GRCh38_${SNAME}.IUPAC.consensus"
    ISSL="${SNAME}_offtargets.issl"
    ALT_TAG="${JOBTAG}_${SNAME}_IUPAC"
    ALT_FASTA="${SDIR}/GRCh38_${SNAME}.IUPAC.fa"
    [[ -s "$ALT_FASTA" ]] || die "Missing consensus FASTA: $ALT_FASTA"

    ALT_MANIDIR="${CR_OUTDIR}/.manifests"
    mkdir -p "$ALT_MANIDIR"

    ALT_RECIPE=$(
        {
        echo "INPUT_FASTA_SHA=$(sha_file "$INPUT_FASTA")"
        echo "ALT_FASTA_SHA=$(sha_file "$ALT_FASTA")"
        echo "BOWTIE_BT2_SHA="
        sha_list_bt2 "${SDIR}/${BOWPFX}"
        echo "ISSL_SHA=$(sha_file "${SDIR}/${ISSL}")"
        echo "OFFTARGET_MODE=$OFFTARGET_MODE"
        echo "OFFTARGET_THRESHOLD=$OFFTARGET_THRESHOLD"
        echo "THREADS=$THREADS"
        echo "VERSIONS:"
        echo "  $(tool_ver bowtie2)"
        echo "  $(tool_ver RNAfold)"
        echo "  $(tool_ver Crackling)"
        echo "  $(tool_ver isslScoreOfftargets)"
        } | sed 's/[[:space:]]\+$//'
    )
    ALT_HASH="$(sha_text "$ALT_RECIPE")"
    ALT_MANIFEST="${ALT_MANIDIR}/${ALT_TAG}_${ALT_HASH}.manifest"

    reuse_alt=false
    if [[ -s "$ALT_MANIFEST" ]]; then
        # shellcheck disable=SC1090
        source "$ALT_MANIFEST" || true
        if [[ -n "${RUN_DIR:-}" && -s "${RUN_DIR}/${ALT_TAG}_final_guides_summary.csv" ]]; then
        reuse_alt=true
        ALT_RUN_DIR="$RUN_DIR"
        echo "# Reusing alt run via manifest: $ALT_MANIFEST"
        else
        echo "# Stale alt manifest; will rebuild: $ALT_TAG"
        reuse_alt=false
        fi
    fi

    if ! $reuse_alt; then
        Run_Custom_Crackling.sh \
        --input  "$INPUT_FASTA" \
        --outdir "$CR_OUTDIR" \
        --personal-dir "$SDIR" \
        --bowtie-prefix "$BOWPFX" \
        --issl "$ISSL" \
        --threads "$THREADS" \
        --offtarget-mode "$OFFTARGET_MODE" \
        --offtarget-threshold "$OFFTARGET_THRESHOLD" \
        --jobtag "$ALT_TAG"

        ALT_RUN_DIR="$(latest_run_dir "$CR_OUTDIR" "$ALT_TAG")"
        [[ -n "$ALT_RUN_DIR" ]] || die "Missing alt run dir for $ALT_TAG"

        {
        echo "HASH=$ALT_HASH"
        echo "RUN_DIR=$ALT_RUN_DIR"
        echo "RECIPE_EPOCH=$(date -u +%s)"
        } > "$ALT_MANIFEST"

        ln -sfn "$ALT_RUN_DIR" "${CR_OUTDIR}/${ALT_TAG}_${ALT_HASH}"
        ln -sfn "$ALT_RUN_DIR" "${CR_OUTDIR}/${ALT_TAG}_latest"
    fi

    OUTCMP="${COMPARE_ROOT}/${SNAME}_IUPAC"
    mkdir -p "$OUTCMP"
    Compare_Crackling_Runs.sh \
        --base-run "$BASE_RUN_DIR" \
        --alt-run  "$ALT_RUN_DIR" \
        --outdir   "$OUTCMP"

    if [[ "$DO_GUIDE_CHECK" == "yes" ]]; then
        Check_Guides_In_Fastas.sh \
        --guides "${BASE_SUMMARY}" \
        --fastas "$REF" "$ALT_FASTA" \
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
