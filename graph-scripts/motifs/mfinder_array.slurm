#!/bin/bash
#SBATCH --job-name=mfinder_array
#SBATCH --output=slurm_%A_%a.out
#SBATCH --error=slurm_%A_%a.err
#SBATCH --array=0-22
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=weka2
#SBATCH --nodelist=voltpve17  # optional if static binary works anywhere

# Variables
FILE_LIST="/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/sorted-loops/loop_file_list.txt"
LOOP_DIR="/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/sorted-loops"
MINDER="/mnt/altnas/work/Kyle.Knightly/contact-network/graph-scripts/motifs/mfinder1.21/mfinder"

# Read file list into array
mapfile -t FILES < "$FILE_LIST"
INPUT_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"
BASE=$(basename "$INPUT_FILE" .txt)

echo "Running mfinder on $INPUT_FILE"

"$MINDER" "$LOOP_DIR/$INPUT_FILE" -s 6 -r 15 -nd -omem -f "$LOOP_DIR/${BASE}"
