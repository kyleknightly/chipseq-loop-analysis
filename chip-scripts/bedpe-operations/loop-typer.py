import pandas as pd
import subprocess
import os
import sys
from itertools import product

# === File paths ===
loop_bedpe = "/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/beds/hepg2-loops.bedpe"
CCRE_bed = "/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/beds/ENCFF546MZK_ENCFF732PJK_ENCFF795ONN_ENCFF357NFO.bed"
output_file = "/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/beds/hepg2-ccre-loop-types.bedpe"

# Temp files
tmp1 = "anchor1.tmp.bed"
tmp2 = "anchor2.tmp.bed"
annot1_ccre = "annot1_ccre.tmp.bed"
annot2_ccre = "annot2_ccre.tmp.bed"

# === Load loops and write anchor bed files ===
loops = pd.read_csv(loop_bedpe, sep="\t", header=None,
                   names=["chr1", "start1", "end1", "chr2", "start2", "end2"])

loops[["chr1", "start1", "end1"]].to_csv(tmp1, sep="\t", header=False, index=False)
loops[["chr2", "start2", "end2"]].to_csv(tmp2, sep="\t", header=False, index=False)

# === Run bedtools intersections with cCRE file ===
bt = "/gpfs0/apps/x86_64/anaconda3/bin/bedtools"

# Intersect anchors with cCRE file and keep both anchor coords and cCRE type
subprocess.run(f"{bt} intersect -a {tmp1} -b {CCRE_bed} -wa -wb > {annot1_ccre}", shell=True, check=True)
subprocess.run(f"{bt} intersect -a {tmp2} -b {CCRE_bed} -wa -wb > {annot2_ccre}", shell=True, check=True)

# === Load annotated anchors ===
def load_ccre_annot(file):
    """Load bedtools intersect output with cCRE annotations"""
    if os.path.getsize(file) == 0:
        # Return empty dataframe with correct columns if file is empty
        return pd.DataFrame(columns=["chr", "start", "end", "type"])
    
    # FIXED: Correct column parsing for your cCRE format
    # Anchor cols (3) + cCRE cols: chr, start, end, name, score, strand, thick_start, thick_end, rgb, type, classification
    df = pd.read_csv(file, sep="\t", header=None, 
                     names=["chr", "start", "end", "ccre_chr", "ccre_start", "ccre_end", 
                           "ccre_name", "ccre_score", "ccre_strand", "ccre_thick_start", 
                           "ccre_thick_end", "ccre_rgb", "type", "classification"])
    
    # Keep only anchor coordinates and cCRE type
    df = df[["chr", "start", "end", "type"]].drop_duplicates()
    return df

a1_ccre = load_ccre_annot(annot1_ccre)
a2_ccre = load_ccre_annot(annot2_ccre)

print(f"🔎 Found {len(a1_ccre)} anchor1 annotations")
print(f"🔎 Found {len(a2_ccre)} anchor2 annotations")

# === Merge types for each anchor ===
def merge_anchor_types(a_ccre, prefix):
    """Group by anchor coordinates and collect all cCRE types for that anchor"""
    if a_ccre.empty:
        return pd.DataFrame(columns=[f"{prefix}chr", f"{prefix}start", f"{prefix}end", f"{prefix}types"]), 0
    
    grouped = a_ccre.groupby(["chr", "start", "end"])["type"].agg(lambda x: list(set(x))).reset_index()
    grouped.columns = [f"{prefix}chr", f"{prefix}start", f"{prefix}end", f"{prefix}types"]
    
    # Count multi-type anchors
    multi_type_count = grouped[grouped[f"{prefix}types"].apply(lambda x: len(x) > 1)].shape[0]
    return grouped, multi_type_count

a1, a1_multi_count = merge_anchor_types(a1_ccre, "anchor1_")
a2, a2_multi_count = merge_anchor_types(a2_ccre, "anchor2_")

print(f"🔎 Found {a1_multi_count} anchor1s with multiple cCRE types")
print(f"🔎 Found {a2_multi_count} anchor2s with multiple cCRE types")

# Show the types we found
if not a1_ccre.empty:
    print(f"🔎 Anchor1 cCRE types found: {sorted(a1_ccre['type'].unique())}")
if not a2_ccre.empty:
    print(f"🔎 Anchor2 cCRE types found: {sorted(a2_ccre['type'].unique())}")

# === Annotate loops ===
annotated = loops.merge(
    a1, left_on=["chr1", "start1", "end1"], right_on=["anchor1_chr", "anchor1_start", "anchor1_end"], how="left"
).merge(
    a2, left_on=["chr2", "start2", "end2"], right_on=["anchor2_chr", "anchor2_start", "anchor2_end"], how="left"
)

# === Drop loops with unannotated anchors ===
missing = annotated[annotated["anchor1_types"].isna() | annotated["anchor2_types"].isna()]
if not missing.empty:
    print(f"⚠️ Dropping {len(missing)} loops with unannotated anchors")
    annotated = annotated.dropna(subset=["anchor1_types", "anchor2_types"])

if annotated.empty:
    print("❌ No loops remain after annotation! Check if coordinates match between files.")
    sys.exit(1)

# === Expand all combinations of anchor1_types × anchor2_types ===
expanded_rows = []
for _, row in annotated.iterrows():
    for t1, t2 in product(row["anchor1_types"], row["anchor2_types"]):
        expanded_rows.append([
            row["chr1"], row["start1"], row["end1"], t1,
            row["chr2"], row["start2"], row["end2"], t2
        ])

expanded_df = pd.DataFrame(expanded_rows,
                          columns=["chr1", "start1", "end1", "type1", "chr2", "start2", "end2", "type2"])

# === Show summary of loop types ===
if not expanded_df.empty:
    type_combinations = expanded_df.groupby(['type1', 'type2']).size().sort_values(ascending=False)
    print(f"\n🔎 Loop type combinations found:")
    for (t1, t2), count in type_combinations.items():
        print(f"   {t1} - {t2}: {count}")

# === Save final expanded output ===
expanded_df.to_csv(output_file, sep="\t", index=False, header=False)
print(f"✅ Done! Saved to {output_file} with {len(expanded_df)} rows.")

# === Clean up ===
for f in [tmp1, tmp2, annot1_ccre, annot2_ccre]:
    if os.path.exists(f):
        os.remove(f)