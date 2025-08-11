"""
This takes an enrichment matrix and an average TSS file and creates a scatter plot
"""

import pandas as pd
import matplotlib.pyplot as plt


distance_df = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/enhancer-promoter/avg_tss_dist.csv")  # or pass this as an argument
distance_dict = dict(zip(distance_df["protein"], distance_df["average_distance_to_tss"]))

cis_enrichment_matrix = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/merged-filtered/tss-ordered/noctcf_cis_enrichment_matrix", sep="\t",index_col=0)
trans_enrichment_matrix = pd.read_csv("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/merged-filtered/tss-ordered/noctcf_trans_enrichment_matrix", sep="\t",index_col=0)

cis_enrichment_dict = cis_enrichment_matrix.sum().to_dict()
trans_enrichment_dict = trans_enrichment_matrix.sum().to_dict()
ratio_dict = {
    k: trans_enrichment_dict[k] / cis_enrichment_dict[k]
    for k in cis_enrichment_dict.keys() & trans_enrichment_dict.keys()
}

plots = {'noctcf-cis-tss-enrichment-correlation': cis_enrichment_dict,
            'noctcf-trans-tss-enrichment-correlation': trans_enrichment_dict,
            'noctcf-trans-over-cis-tss-enrichment-correlation': ratio_dict}

for name, enrichment_dict in plots.items():
    common_proteins = set(enrichment_dict) & set(distance_dict)

    # Extract data
    x = [enrichment_dict[p] for p in common_proteins]         # enrichment values
    y = [distance_dict[p] for p in common_proteins]    # average TSS distances
    labels = list(common_proteins)

    # Plot
    plt.figure(figsize=(16, 12))
    plt.scatter(x, y, s=1)
    # plt.xlim(0, 1.05*max(max(cis_enrichment_dict.values()), max(trans_enrichment_dict.values())))
    plt.ylim(0, 1.05*max(distance_dict.values()))

    # Add labels
    for i, label in enumerate(labels):
        plt.text(x[i], y[i], label, fontsize=8, ha='right', va='bottom')

    # Axis labels
    plt.xlabel("Enrichment (Column Sum)")
    plt.ylabel("Average Distance to TSS")
    plt.title("Protein Enrichment vs. TSS Distance")
    plt.savefig(name)
