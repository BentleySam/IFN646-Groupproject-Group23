#!/usr/bin/env bash
set -euo pipefail

GUIDES="${1:-}"    # CSV or txt (one 23-mer per line)
FASTA_DIR="${2:-}" # dir with *.fa / *.fasta / *.fa.gz / *.fasta.gz
[[ -n "$GUIDES" && -n "$FASTA_DIR" ]] || { echo "Usage: $(basename "$0") <guides.csv|guides.txt> <fasta_dir>"; exit 1; }
[[ -d "$FASTA_DIR" ]] || { echo "No such dir: $FASTA_DIR" >&2; exit 2; }

# build list of 23-mers (spacer+PAM)
tmp="$(mktemp)"; trap 'rm -f "$tmp"' EXIT
if [[ "$GUIDES" == *.csv ]]; then
  python3 - "$GUIDES" "$tmp" <<'PY'
import sys,csv
src,out=sys.argv[1],sys.argv[2]
n=0; kept=0
with open(src,newline='') as f, open(out,'w') as g:
    r=csv.DictReader(f)
    for row in r:
        n+=1
        s=(row.get('spacer','')+row.get('pam','')).upper().strip()
        if len(s)==23 and all(c in "ACGTN" for c in s):
            g.write(s+'\n'); kept+=1
print(f"# read {n} rows; kept {kept} guides of length 23", file=sys.stderr)
PY
else
  # plain list: keep only exact 23-mers
  awk 'BEGIN{IGNORECASE=1} length($0)==23 && $0 ~ /^[ACGTN]+$/ {print toupper($0)}' "$GUIDES" > "$tmp"
fi

python3 - "$tmp" "$FASTA_DIR" <<'PY'
import sys,glob,gzip,io,os
guides=open(sys.argv[1]).read().split()
fa_dir=sys.argv[2]

# reverse complement (uppercase)
rc_tbl=str.maketrans('ACGTN','TGCAN')
def revcomp(s): return s.translate(rc_tbl)[::-1]

# expand fasta patterns
fas = sorted(
    glob.glob(os.path.join(fa_dir, "*.fa")) +
    glob.glob(os.path.join(fa_dir, "*.fasta")) +
    glob.glob(os.path.join(fa_dir, "*.fa.gz")) +
    glob.glob(os.path.join(fa_dir, "*.fasta.gz"))
)

print("fasta,guide23,present,strand")
for fa in fas:
    # open gz or plain
    if fa.endswith(".gz"):
        fh = gzip.open(fa, "rt")
    else:
        fh = open(fa, "r")

    try:
        # read full seq uppercase (simple, OK for your machine)
        seq = []
        for line in fh:
            if line.startswith(">"): continue
            seq.append(line.strip().upper())
        s = "".join(seq)
    finally:
        fh.close()

    for g in guides:
        if g in s:
            print(f"{os.path.basename(fa)},{g},YES,+")
        else:
            gr = revcomp(g)
            if gr in s:
                print(f"{os.path.basename(fa)},{g},YES,-")
            else:
                print(f"{os.path.basename(fa)},{g},NO,.")
PY
