#!/bin/sh
# Put country specific datasets in output/country
# depends: jq
set -exou pipefail

COUNTRIES="$(jq -r '.[]' data/countries.json | tr '\n' ' ')"
for c in $COUNTRIES; do
    python analysis.py --country=$c
    python analysis.py --country=$c --age
    python analysis.py --country=$c --age --group
    python analysis.py --country=$c --gender
    python analysis.py --country=$c --gender --group
done
