import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load CSV
input_csv = "/mnt/altnas/work/Kyle.Knightly/residue_type_counts.csv"
df = pd.read_csv(input_csv, index_col=0)

# Identify residue type columns
residue_types = ['nonpolar', 'polar_uncharged', 'polar_positive', 'polar_negative']
df = df[[col for col in residue_types if col in df.columns]]  # drop unknown if not needed

# Compute proportions
df_prop = df.div(df.sum(axis=1), axis=0).fillna(0)

# Order proteins (if you have a custom order, insert here)
protein_order = list(df_prop.index)
df_prop = df_prop.loc[protein_order]

# Plotting
available_proteins = df_prop.index.tolist()
x_spacing = 10
x_pos = np.arange(len(available_proteins)) * x_spacing
bar_width = 10
colors = sns.color_palette("Paired", n_colors=len(residue_types))

# Create the figure
fig, ax = plt.subplots(figsize=(len(available_proteins) * 0.3, 4))

# Build stacked bars
bottom = np.zeros(len(available_proteins))
bars = []

for i, res_type in enumerate(residue_types):
    values = df_prop[res_type].values
    bar = ax.bar(x_pos, values, bottom=bottom, label=res_type.replace("_", " ").title(),
                 color=colors[i], alpha=0.9, edgecolor='white', linewidth=0.2, width=bar_width)
    bars.append(bar)
    bottom += values

# Labeling and styling
ax.set_xlabel("Proteins", fontsize=12, fontweight='bold')
ax.set_ylabel("Proportion of Residues", fontsize=12, fontweight='bold')
ax.set_title(f"Amino Acid Residue Composition per Protein\n(Stacked Proportions)",
             fontsize=14, fontweight='bold')
ax.set_xticks(x_pos)
ax.set_xticklabels(available_proteins, rotation=90, ha='center', fontsize=8)
ax.legend(title="Residue Type", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
ax.grid(False)
ax.set_ylim(0, 1.0)

# Add percentage labels for larger slices
if len(available_proteins) <= 20:
    bottom = np.zeros(len(available_proteins))
    for i, res_type in enumerate(residue_types):
        values = df_prop[res_type].values
        for j, value in enumerate(values):
            if value > 0.05:
                ax.text(j * x_spacing, bottom[j] + value / 2, f"{value:.0%}",
                        ha='center', va='center', fontweight='bold', fontsize=6)
        bottom += values

plt.tight_layout()
plt.savefig("residue_type_distribution_styled.png", dpi=300, bbox_inches="tight")
plt.show()

# Print proportions for reference
print("\n" + "="*80)
print("RESIDUE TYPE PROPORTIONS BY PROTEIN")
print("="*80)
print(df_prop.round(3))
