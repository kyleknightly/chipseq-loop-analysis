#!/bin/bash

# Set the directory containing the BED files
BED_DIR="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/new-filtered"

# Set the directory to store the output files
OUTPUT_DIR="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/new-merged-filtered"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Move into the directory with the BED files
cd "$BED_DIR" || exit 1

# Loop over each unique PROTNAME
for PROTNAME in $(ls *.bed | awk -F'-' '{print $1}' | sort | uniq); do
    # Find all BED files with this PROTNAME
    FILES=""
    for FILE in ${PROTNAME}-*.bed; do
        CLEAN_FILE="${FILE%.bed}-clean.bed"
        
        
        # Sort the file, keep only the first 3 columns, and ensure correct tab delimiters
        sort -k1,1 -k2,2n "$FILE" | cut -f1-3 --output-delimiter=$'\t' > "$CLEAN_FILE"
        
        # Add the cleaned and sorted file to the list of files for multiinter
        FILES="$FILES $CLEAN_FILE"
    done
    # awk 'NF != 3' CTCF-human_HepG2_ENCFF119PKI-clean.bed
    # head -n 10 CTCF-human_HepG2_ENCFF119PKI-clean.bed

    # If there are multiple files, run bedtools multiinter
    if [[ $(echo "$FILES" | wc -w) -gt 1 ]]; then
        bedtools multiinter -i $FILES > "${OUTPUT_DIR}/${PROTNAME}-multiinter.bed"
    else
        # If there's only one file, copy it directly as the output
        cp $FILES "${OUTPUT_DIR}/${PROTNAME}-multiinter.bed"
    fi

    # Clean up the temporary cleaned files
    rm $FILES
done