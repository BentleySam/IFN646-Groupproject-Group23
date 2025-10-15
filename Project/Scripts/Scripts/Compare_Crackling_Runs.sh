#!/usr/bin/env bash
set -euo pipefail
BASE="${1:-}"; SHIFT_DIR="${2:-}"

if [[ -z "$BASE" || -z "$SHIFT_DIR" ]]; then
  echo "Usage: $(basename "$0") <base_guides_summary.csv> <personalized_output_dir>"
  exit 1
fi
[[ -s "$BASE" ]] || { echo "Base CSV not found: $BASE" >&2; exit 2; }

shopt -s nullglob
csvs=( "$SHIFT_DIR"/*/*_final_guides_summary.csv )
((${#csvs[@]})) || { echo "No personalized summaries in $SHIFT_DIR" >&2; exit 0; }

python3 - "$BASE" "${csvs[@]}" <<'PY'
import sys, csv, os, math
base = sys.argv[1]
pers = sys.argv[2:]

def load(path):
    m={}
    with open(path, newline='') as f:
        r=csv.DictReader(f)
        for row in r:
            key=(row.get('spacer','')+row.get('pam','')).upper()
            m[key]=row
    return m

def fnum(x):
    if x in ("","?") or x is None: return math.nan
    try: return float(x)
    except: return math.nan

base_map = load(base)

for p in pers:
    new_map = load(p)
    out = os.path.splitext(p)[0] + "_delta.tsv"
    with open(out, "w", newline='') as g:
        w=csv.writer(g, delimiter='\t')
        w.writerow(["guide23","in_new","delta_MIT","delta_CFD","base_MIT","new_MIT","base_CFD","new_CFD"])
        for k,b in base_map.items():
            n = new_map.get(k)
            bMIT, bCFD = fnum(b.get("mitOfftargetscore")), fnum(b.get("cfdOfftargetscore"))
            if n:
                nMIT, nCFD = fnum(n.get("mitOfftargetscore")), fnum(n.get("cfdOfftargetscore"))
                dMIT = "" if (math.isnan(bMIT) or math.isnan(nMIT)) else nMIT-bMIT
                dCFD = "" if (math.isnan(bCFD) or math.isnan(nCFD)) else nCFD-bCFD
                w.writerow([k,"YES",dMIT,dCFD,
                            "" if math.isnan(bMIT) else bMIT, "" if math.isnan(nMIT) else nMIT,
                            "" if math.isnan(bCFD) else bCFD, "" if math.isnan(nCFD) else nCFD])
            else:
                w.writerow([k,"NO","","",
                            "" if math.isnan(bMIT) else bMIT, "",
                            "" if math.isnan(bCFD) else bCFD, ""])
    print("Wrote", out)
PY
