#!/bin/bash

# Set the directory containing the BED files
INDIR="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/filtered-tracks"

# Set the directory to store the output files
OUTDIR="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/merged-filtered-tracks"

# Create the output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Move into the directory with the BED files
cd "$INDIR" || exit 1

# 1) Identify all unique PROTNAMEs by splitting on the first dash (-).
#    We list all .bed files, cut off everything after the first '-', then sort uniquely.
for PROTNAME in $(ls *.bed | awk -F'-' '{print $1}' | sort | uniq); do
    
    echo "Now processing: $PROTNAME"

    # 2) Gather all BED files that match this PROTNAME
    FILES=( "${PROTNAME}-"*.bed )

    # 3) If no files match (edge case), continue
    if [ ${#FILES[@]} -eq 0 ]; then
        echo "  No files found for $PROTNAME, skipping."
        continue
    fi

    # 4) Concatenate all those files, then sort by chromosome and coordinate.
    #    -k1,1V  sorts chr1 < chr2 < chr10
    #    -k2,2n  sorts numerically by start
    #    (Optionally) -k3,3n for tie-breaking on end.
    #    We put the result in a temporary file.
    COMBINED_FILE=$(mktemp)
    cat "${FILES[@]}" | sort -k1,1V -k2,2n > "$COMBINED_FILE"

    # 5) Merge intervals across the combined, sorted file
    bedtools merge -i "$COMBINED_FILE" > "${OUTDIR}/${PROTNAME}-merged.bed"

    # 6) Clean up
    rm "$COMBINED_FILE"

    
done
