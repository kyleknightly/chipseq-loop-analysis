#!/bin/bash

input_dir="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/K562/beds/chip-tracks"
output_dir="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/K562/beds/filtered"

mkdir -p "$output_dir"

for file in "$input_dir"/*.bed "$input_dir"/*.narrowPeak; do
  [ -e "$file" ] || continue  # Skip if no match

  base=$(basename "$file")

  # Count non-empty lines
  total_lines=$(grep -cv '^$' "$file")

  # Compute top 10% (rounded up)
  top_n=$(awk -v n="$total_lines" 'BEGIN { printf("%d", (n * 0.10 + 0.5)) }')

  # Sort and keep top 10% of lines
  sort -k7,7nr "$file" | head -n "$top_n" > "$output_dir/$base"
done