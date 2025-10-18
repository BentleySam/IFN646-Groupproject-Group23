#!/usr/bin/env bash
# Reproducible Crackling runner with auto index build (Bowtie2 + ISSL)
set -euo pipefail

# ---------- defaults ----------
INPUT=""
OUTDIR="$HOME/cracklingOutput"
GENOME_DIR="$HOME/genomes"
BOWTIE_PREFIX="GRCh38"               # index basename (e.g. GRCh38.*bt2)
ISSL_FILE="GRCh38_offtargets.issl"   # filename under GENOME_DIR or absolute
REF_FASTA=""                         # if set and indexes missing, build from this FASTA
THREADS="$(nproc || echo 8)"
OFFTARGET_MODE="avg"                 # 'avg' or 'and'
OFFTARGET_THRESHOLD="65"             # 65 lenient, 75 stricter
GUIDE_LEN="20"                       # guide length for ISSL
SLICE_WIDTH="8"                      # ISSL slice width (bits); 8 suits 4 mismatches
JOBTAG="CracklingRun"
CONDA_ENV="crackling"

GFF_PATH="${GFF_PATH:-$GENOME_DIR/GRCh38.gff3}"   # or leave empty to skip
OPTIMISATION="${OPTIMISATION:-high}"              # ultralow|low|medium|high
CONSENSUS_N="${CONSENSUS_N:-2}"
RNAFOLD_THREADS="$THREADS"
BOWTIE_THREADS="$THREADS"
ISSL_THREADS="$THREADS"
PAGE_LEN="${PAGE_LEN:-5000000}"
LOW_E="-30"
HIGH_E="-18"

usage(){ cat <<EOF
Usage: $0 --input FASTA [options]
  --input PATH              Target region FASTA
  --outdir DIR              Output dir (default: $OUTDIR)
  --genome-dir DIR          Dir for genome indices (default: $GENOME_DIR)
  --bowtie-prefix NAME      Bowtie2 basename (default: $BOWTIE_PREFIX)
  --issl FILE               ISSL file (abs or relative to --genome-dir) (default: $ISSL_FILE)
  --ref-fasta PATH          Reference FASTA (used to build indices if missing)
  --threads N               Threads (default: $THREADS)
  --offtarget-mode MODE     'avg' or 'and' (default: $OFFTARGET_MODE)
  --offtarget-threshold N   65 / 75 etc. (default: $OFFTARGET_THRESHOLD)
  --guide-len N             Guide length for ISSL (default: $GUIDE_LEN)
  --slice-width N           ISSL slice width in bits (default: $SLICE_WIDTH)
  --jobtag NAME             Label for outputs (default: $JOBTAG)
  --env NAME                Conda env (default: $CONDA_ENV)
EOF
exit 1; }

# ---------- args ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --genome-dir) GENOME_DIR="$2"; shift 2;;
    --bowtie-prefix) BOWTIE_PREFIX="$2"; shift 2;;
    --issl) ISSL_FILE="$2"; shift 2;;
    --ref-fasta) REF_FASTA="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --offtarget-mode) OFFTARGET_MODE="$2"; shift 2;;
    --offtarget-threshold) OFFTARGET_THRESHOLD="$2"; shift 2;;
    --guide-len) GUIDE_LEN="$2"; shift 2;;
    --slice-width) SLICE_WIDTH="$2"; shift 2;;
    --jobtag) JOBTAG="$2"; shift 2;;
    --env) CONDA_ENV="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done
[[ -n "${INPUT:-}" ]] || usage

# ---------- helpers ----------
die(){ echo "ERROR: $*" >&2; exit 2; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing tool: $1"; }
abspath(){ python3 - "$1" <<'PY'
import os,sys; print(os.path.abspath(sys.argv[1]))
PY
}

# ---------- env/tools ----------
if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV" || die "Failed to activate env '$CONDA_ENV'"
fi
for t in python3 bowtie2 bowtie2-build RNAfold Crackling; do need "$t"; done
# Crackling utils
if ! command -v isslScoreOfftargets >/dev/null 2>&1; then
  # try to add ~/Crackling/bin if present
  [[ -d "$HOME/Crackling/bin" ]] && export PATH="$HOME/Crackling/bin:$PATH"
  need isslScoreOfftargets
fi
need extractOfftargets

CREATE_ISSL="$(command -v createIsslIndex || command -v isslCreateIndex || true)"
[[ -n "$CREATE_ISSL" ]] || die "Missing tool: createIsslIndex (or isslCreateIndex)"

# ---------- provenance ----------
echo "# $(date) :: starting $0"
echo "# Python: $(python3 --version 2>&1)"
echo "# Bowtie2: $(bowtie2 --version | head -n1)"
echo "# RNAfold: $(RNAfold --version 2>/dev/null | head -n1 || echo 'unknown')"
echo "# Crackling bin: $(command -v Crackling)"

# ---------- paths & logs ----------
INPUT="$(abspath "$INPUT")"; [[ -s "$INPUT" ]] || die "Input FASTA missing/empty: $INPUT"
mkdir -p "$OUTDIR"; RUNSTAMP="$(date +%Y%m%d_%H%M%S)"; RUNROOT="$OUTDIR/${JOBTAG}_${RUNSTAMP}"; mkdir -p "$RUNROOT"
MASTER_LOG="$RUNROOT/${JOBTAG}_run.log"; exec > >(tee -a "$MASTER_LOG") 2>&1; set -x
command -v dos2unix >/dev/null 2>&1 && dos2unix -q "$INPUT" || true

# ---------- resolve genome assets ----------
BOWTIE_PREFIX_PATH="$GENOME_DIR/$BOWTIE_PREFIX"
# ISSL path (abs or under genome dir)
if [[ -f "$ISSL_FILE" ]]; then ISSL_PATH="$(abspath "$ISSL_FILE")"
else ISSL_PATH="$(abspath "$GENOME_DIR/$ISSL_FILE")"; fi

have_bt2=true
for sfx in 1 2 3 4 "rev.1" "rev.2"; do
  [[ -s "${BOWTIE_PREFIX_PATH}.${sfx}.bt2" ]] || { have_bt2=false; break; }
done
have_issl=true; [[ -s "$ISSL_PATH" ]] || have_issl=false

# ---------- auto-build if missing ----------
build_from_ref(){
  local ref="$1"; [[ -s "$ref" ]] || die "REF_FASTA missing/empty: $ref"
  command -v samtools >/dev/null 2>&1 && samtools faidx "$ref" || true
  # Bowtie2 index
  bowtie2-build --threads "$THREADS" "$ref" "$BOWTIE_PREFIX_PATH"
  # Off-targets + ISSL
  local offtxt="$GENOME_DIR/${BOWTIE_PREFIX}_offtargets.txt"
  extractOfftargets "$offtxt" "$ref"
  "$CREATE_ISSL" -t "$offtxt" -l "$GUIDE_LEN" -w "$SLICE_WIDTH" -o "$ISSL_PATH"
  echo "# Using createIsslIndex at: $CREATE_ISSL"
}

if ! $have_bt2 || ! $have_issl; then
  [[ -n "$REF_FASTA" ]] || die "Indices missing. Provide --ref-fasta to build."
  mkdir -p "$GENOME_DIR"
  REF_FASTA="$(abspath "$REF_FASTA")"
  build_from_ref "$REF_FASTA"
fi

# Re-check
for sfx in 1 2 3 4 "rev.1" "rev.2"; do
  [[ -s "${BOWTIE_PREFIX_PATH}.${sfx}.bt2" ]] || die "Bowtie2 index still missing: ${BOWTIE_PREFIX_PATH}.${sfx}.bt2"
done
[[ -s "$ISSL_PATH" ]] || die "ISSL file still missing: $ISSL_PATH"

ISSLSCORE_BIN="$(command -v isslScoreOfftargets)"
[[ -n "$ISSLSCORE_BIN" ]] || die "isslScoreOfftargets not found in PATH"
BOWTIE2_BIN="$(command -v bowtie2 || echo bowtie2)"
RNAFOLD_BIN="$(command -v RNAfold || echo RNAfold)"

# try to find the sgRNAScorer2 model automatically; fall back if needed
SG_MODEL="${SG_MODEL:-$(python3 - <<'PY'
import importlib.resources, sys
try:
    import crackling
    with importlib.resources.path('crackling.utils.data', 'model-py3.txt') as p:
        print(p)
except Exception:
    pass
PY
)}"
: "${SG_MODEL:="/mnt/c/Users/Axis3/Github Repos/Crackling/src/crackling/utils/data/model-py3.txt"}"

### uses my windows path as fallback, change to own


echo "# Using genome:"
echo "#   Bowtie2 prefix : $BOWTIE_PREFIX_PATH"
echo "#   ISSL index     : $ISSL_PATH"

[[ -f "$GFF_PATH" ]] || GFF_PATH=""




# ---------- config for this run ----------
CFG="$RUNROOT/${JOBTAG}.ini"
LOGPREFIX="$RUNROOT/${JOBTAG}"
cat > "$CFG" <<EOF
; Auto-generated by runcrackling.sh

[general]
name = ${JOBTAG}
optimisation = ${OPTIMISATION}

[consensus]
n = ${CONSENSUS_N}
mm10db = True
sgrnascorer2 = True
chopchop = True

[input]
exon-sequences = ${INPUT}
offtarget-sites = ${ISSL_PATH}
$( [[ -n "$GFF_PATH" ]] && echo "gff-annotation = ${GFF_PATH}" )
bowtie2-index = ${BOWTIE_PREFIX_PATH}
batch-size = ${PAGE_LEN}

[output]
dir = ${RUNROOT}/
filename = guides.txt
delimiter = ,

[offtargetscore]
enabled = True
binary = ${ISSLSCORE_BIN}
method = ${OFFTARGET_MODE}
threads = ${ISSL_THREADS}
page-length = ${PAGE_LEN}
score-threshold = ${OFFTARGET_THRESHOLD}
max-distance = 4

[sgrnascorer2]
model = ${SG_MODEL}
score-threshold = 0

[bowtie2]
binary = ${BOWTIE2_BIN}
threads = ${BOWTIE_THREADS}
page-length = ${PAGE_LEN}

[rnafold]
binary = ${RNAFOLD_BIN}
threads = ${RNAFOLD_THREADS}
page-length = ${PAGE_LEN}
low_energy_threshold = ${LOW_E}
high_energy_threshold = ${HIGH_E}
EOF
echo "# Wrote config: $CFG"

# ---------- run crackling ----------
rm -f "${LOGPREFIX}-"*.log "${LOGPREFIX}-"*.errlog "${LOGPREFIX}-guides.txt" || true
Crackling -c "$CFG"

# ---------- summary CSV (passedOffTargetScore == 1) ----------
GUIDES="${LOGPREFIX}-guides.txt"; [[ -s "$GUIDES" ]] || die "Guides file missing/empty: $GUIDES"
SUMMARY="$RUNROOT/${JOBTAG}_final_guides_summary.csv"
python3 - "$GUIDES" "$SUMMARY" <<'PY'
import csv, sys
src, dst = sys.argv[1], sys.argv[2]
with open(src, newline='') as f, open(dst, 'w', newline='') as g:
    r = csv.DictReader(f)
    cols = ["spacer","pam","strand","bowtieChr","bowtieStart","bowtieEnd",
            "sgrnascorer2score","acceptedByMm10db",
            "mitOfftargetscore","cfdOfftargetscore","consensusCount"]
    w = csv.DictWriter(g, fieldnames=cols); w.writeheader()
    for row in r:
        if row.get("passedOffTargetScore","0") != "1": continue
        seq = (row.get("seq","") or "").strip()
        w.writerow({
            "spacer": seq[:20], "pam": seq[20:], "strand": row.get("strand",""),
            "bowtieChr": row.get("bowtieChr",""), "bowtieStart": row.get("bowtieStart",""),
            "bowtieEnd": row.get("bowtieEnd",""),
            "sgrnascorer2score": row.get("sgrnascorer2score",""),
            "acceptedByMm10db": row.get("acceptedByMm10db",""),
            "mitOfftargetscore": row.get("mitOfftargetscore",""),
            "cfdOfftargetscore": row.get("cfdOfftargetscore",""),
            "consensusCount": row.get("consensusCount",""),
        })
print("Wrote", dst)
PY

set +x
echo
echo "=== DONE ==="
echo "Run root         : $RUNROOT"
echo "Master run log   : $MASTER_LOG"
echo "Crackling logs   : ${LOGPREFIX}-*.log / *.errlog"
echo "Guides (full)    : $GUIDES"
echo "Guides (summary) : $SUMMARY"
echo
echo "Don't forget to cite: Bradford, Chappell, Perrin. The CRISPR Journal (2022) doi:10.1089/crispr.2021.0102"
