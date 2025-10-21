#!/usr/bin/env bash
set -euo pipefail

# Config (override via env)
ASSETS_DIR="${ASSETS_DIR:-$HOME/genomes/personalised}"  # contains *.idx folders
CRACKLING_INPUT="${CRACKLING_INPUT:-$HOME/crackling_inputs/CXCL10_exon1.fa}"
CRACKLING_OUTDIR="${CRACKLING_OUTDIR:-$HOME/cracklingOutput}"
THREADS="${THREADS:-24}"
OFFTARGET_MODE="${OFFTARGET_MODE:-avg}"
OFFTARGET_THRESHOLD="${OFFTARGET_THRESHOLD:-65}"

need(){ command -v "$1" >/dev/null 2>&1 || { echo "Missing: $1" >&2; exit 2; }; }
need bash; need Crackling
[[ -x "$HOME/runcrackling.sh" ]] || { echo "runcrackling.sh not found/executable in ~" >&2; exit 2; }

log(){ echo "[$(date +%F\ %T)] $*"; }

shopt -s nullglob
idxs=( "$ASSETS_DIR"/*.idx )
((${#idxs[@]})) || { echo "No *.idx under $ASSETS_DIR" >&2; exit 2; }

for d in "${idxs[@]}"; do
  issl="$d/offtargets.issl"
  btprefix="$d/consensus"
  [[ -s "$issl" ]] || { echo "Missing $issl — skip $d" >&2; continue; }
  [[ -s "${btprefix}.1.bt2" ]] || { echo "Missing Bowtie2 index in $d — skip" >&2; continue; }

  job="CR_$(basename "$d" .idx)"
  log "Running Crackling: $job"
  bash "$HOME/runcrackling.sh" \
    --input "$CRACKLING_INPUT" \
    --outdir "$CRACKLING_OUTDIR" \
    --genome-dir "$d" \
    --bowtie-prefix consensus \
    --issl "$(basename "$issl")" \
    --threads "$THREADS" \
    --offtarget-mode "$OFFTARGET_MODE" \
    --offtarget-threshold "$OFFTARGET_THRESHOLD" \
    --jobtag "$job"
done

log "Done."
