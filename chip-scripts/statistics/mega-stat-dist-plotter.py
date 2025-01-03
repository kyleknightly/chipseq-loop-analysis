import pandas as pd
import ast
import matplotlib.pyplot as plt
import numpy as np

file = '/mnt/altnas/work/Kyle.Knightly/med-extended-shuffle-STRING-group-stats.tsv'
df = pd.read_csv(file, sep='\t', header = 0)
df = df.rename(columns={'Unnamed: 0': 'prots'})
df['prots']=df['prots'].apply(ast.literal_eval)
df['prot'] = df['prots'].apply(lambda x: x[0])
df['shufSNRs']=df['shufSNRs'].apply(ast.literal_eval)
df['shufcophs']=df['shufcophs'].apply(ast.literal_eval)
# print(df)

# List of proteins to iterate over
proteins = ['CTCF', 'RAD21', 'SMC3', 'STAG1']
# ['FOSL1', 'FOSL2', 'JUNB', 'JUN']
# ['NFYA', 'NFYB', 'NFYC']


# Create a figure with 2 columns and len(proteins) rows
fig, axs = plt.subplots(len(proteins), 2, figsize=(12, 6 * len(proteins)))

# Loop through each protein and create two plots (SNR and coph) for each
for i, prot in enumerate(proteins):
    print(prot)
    df_prot = df[df['prot'] == prot].iloc[0]  # Filter the dataframe for the current protein

    # Extract values for plotting
    SNR = df_prot['SNR']
    coph = df_prot['coph']
    shufSNRs = df_prot['shufSNRs']
    shufcophs = df_prot['shufcophs']

    # Plot for SNR vs shufSNRs in the first column of the current row
    axs[i, 0].hist(shufSNRs, bins=50, alpha=0.7, color='blue', label='Shuffled SNRs')
    axs[i, 0].axvline(SNR, color='red', linestyle='dashed', linewidth=2, label=f'Actual SNR ({SNR})')
    axs[i, 0].set_title(f'SNR Distribution for {prot}')
    axs[i, 0].set_xlabel('SNR')
    axs[i, 0].set_ylabel('Frequency')
    axs[i, 0].legend()

    # Plot for coph vs shufcophs in the second column of the current row
    axs[i, 1].hist(shufcophs, bins=50, alpha=0.7, color='green', label='Shuffled cophs')
    axs[i, 1].axvline(coph, color='red', linestyle='dashed', linewidth=2, label=f'Actual coph ({coph})')
    axs[i, 1].set_title(f'coph Distribution for {prot}')
    axs[i, 1].set_xlabel('coph')
    axs[i, 1].set_ylabel('Frequency')
    axs[i, 1].legend()

# Adjust layout to prevent overlapping
plt.tight_layout()

plt.savefig('cohesin-med-stat-distributions.png')
