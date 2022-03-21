#!/bin/sh
# Put country specific datasets in output/country
# depends: jq
set -exou pipefail
COUNTRIES="$(jq -r '.[]' data/countries.json | tr '\n' ' ')"
for c in $COUNTRIES; do
    python3 analysis.py --country=$c
    python3 analysis.py --country=$c --age
    python3 analysis.py --country=$c --age --group
    python3 analysis.py --country=$c --gender
    python3 analysis.py --country=$c --gender --group
done
