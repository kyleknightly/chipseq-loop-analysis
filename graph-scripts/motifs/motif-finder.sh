#!/bin/bash

IN_DIR="/mnt/altnas/work/Kyle.Knightly/component_networks"
MINDER="/mnt/altnas/work/Kyle.Knightly/contact-network/graph-scripts/motifs/mfinder1.21/mfinder"

for file in "$IN_DIR"/*component*.txt; do
    base=$(basename "$file" .txt)
    echo "Running mfinder on $base"

    "$MINDER" "$file" -s 4 -r 15 -nd -omem
done
