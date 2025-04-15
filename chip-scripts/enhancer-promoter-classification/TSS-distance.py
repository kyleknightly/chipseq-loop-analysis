

import os
import pandas as pd
import numpy as np
from glob import glob

# === CONFIGURATION ===
tss_bed_path = "/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/TSSs_ENCFF246WDH.bed"         # BED file of TSSs
chipseq_folder = "/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/merged-filtered-tracks"  # Folder containing ChIP-seq BEDs (one per protein)
output_csv = "average_distances_to_tss.csv"

# Load TSSs
tss_df = pd.read_csv(tss_bed_path, sep="\t", header=None, engine="python")
tss_df = tss_df[[0, 1, 2]]
tss_df.columns = ["chrom", "start", "end"]

print(tss_df)
# We'll use the TSS center
tss_df["tss"] = ((tss_df["start"] + tss_df["end"]) // 2).astype(int)

# Group TSSs by chromosome for fast lookup
tss_by_chr = {chrom: tss_df[tss_df["chrom"] == chrom]["tss"].values for chrom in tss_df["chrom"].unique()}

results = []

# Loop through each protein BED file
for bed_file in glob(os.path.join(chipseq_folder, "*.bed")):
    protein_name = os.path.basename(bed_file).split("-")[0]

    # Load ChIP peaks
    chip_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end"])

    distances = []

    for _, row in chip_df.iterrows():
        chrom = row["chrom"]
        peak_center = (row["start"] + row["end"]) // 2
        # print(chrom)
        # print(tss_by_chr)
        if chrom in tss_by_chr:
            tss_positions = tss_by_chr[chrom]
            nearest_distance = np.min(np.abs(tss_positions - peak_center))
            distances.append(nearest_distance)
    print(distances)
    avg_distance = np.mean(distances) if distances else float("nan")
    results.append({"protein": protein_name, "average_distance_to_tss": avg_distance})

print(results)
# Save results
df_out = pd.DataFrame(results)
df_out.to_csv(output_csv, index=False)
