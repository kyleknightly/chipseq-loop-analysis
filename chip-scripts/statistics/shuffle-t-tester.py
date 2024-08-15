"""
This takes a heatmap and a file in the format of the output of interactions-list-maker
For each group it calculates SNRs and cophenetic distances, and shuffles the matrix labels
n times to create a random sample to compare it to for p values
"""

import pandas as pd
from scipy.cluster.hierarchy import linkage, cophenet, leaves_list
from scipy.spatial.distance import pdist, squareform
import scipy.stats as stats
from scipy.stats import ttest_ind
import numpy as np
import ast
import os
import random

# Load the heatmap data from a TSV file
heatmap_path = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/same-anchor-interactions/geq250-geq50-anchor-enrichments-pseudo-anchor-contact-matrix.tsv'
df = pd.read_csv(heatmap_path, sep='\t', index_col=0)
#print(len(df))
# Create linkage matrix
Z = linkage(df.values, method='ward')
ordered_index = leaves_list(Z)
protein_list = df.index.tolist()
# Calculate cophenetic distance matrix
coph_dists = cophenet(Z, pdist(df.values))[1]
#print(cophenet(Z, pdist(df.values))[0])
#print(coph_dists)
coph_dist_matrix = squareform(coph_dists)
coph_dist_df = pd.DataFrame(coph_dist_matrix, index=df.index, columns=df.index)

ordered_df = df.iloc[ordered_index, ordered_index]


def coph_distance(protein1, protein2):
    try:
        return coph_dist_df.loc[protein1, protein2]
    except KeyError:
        return np.nan


def coph(list):
    dists = []
    for i in range(len(list)):
        for j in range(i + 1, len(list)):
            protein1 = list[i]
            protein2 = list[j]
            cophenetic_distance = coph_distance(protein1, protein2)
            dists.append(cophenetic_distance)
    return np.mean(dists)
    
    
def SNR(prots): 
    #print(prots)
    prots_to_remove = [protein for protein in prots if protein not in df.index]
    for protein in prots_to_remove:
        prots.remove(protein)
    #print(prots)
    signal = df.loc[prots, prots]
    np.fill_diagonal(signal.values, np.nan)
    signal_average = signal.mean(axis=1, skipna=True)

    N_average = {}
    SNRs = {}
    
    for prot in prots:
        row = df.loc[prot]
        N_row = row.drop(index=prots)
        average_value = N_row.mean()
        N_average[prot] = average_value
        SNRs[prot] = signal_average[prot] / N_average[prot]
    return np.mean(list(SNRs.values()))


list_file = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/STRING/groups/strong-interaction-lists.tsv'
list_df = pd.read_csv(list_file, sep='\t', header=0)
list_df['associated_proteins'] = list_df['associated_proteins'].apply(ast.literal_eval)
lists = [[row['protein']] + row['associated_proteins'] for index, row in list_df.iterrows()]
print(lists)
stats_df = pd.DataFrame(columns = ['coph', 'SNR'])
print(stats_df)
for protein_list in lists:
    #print(protein_list)
    s_n_r = SNR(protein_list)
    avg_cophenetic_distance = coph(protein_list)
    
    # Assign the calculated values to the DataFrame
    stats_df.loc[str(protein_list)] = [avg_cophenetic_distance, s_n_r]

stats_df['shufcophs'] = [[] for _ in range(len(stats_df))]
stats_df['shufSNRs'] = [[] for _ in range(len(stats_df))]

print(ordered_df.index)
idtoprot = {n: ordered_df.index[n] for n in range(len(ordered_df.index))}
prottoid = {ordered_df.index[n]:n for n in range(len(ordered_df.index))}

"""SHUFFLING"""
for i in range(100):
    rand_ints = random.sample(range(len(ordered_df.index)), len(ordered_df.index))
    for index, row in stats_df.iterrows():
        # print("----")
        # print(index)
        current_nums = [prottoid[prot] for prot in ast.literal_eval(index)]
        new_nums = [rand_ints[num] for num in current_nums]
        new_prots = [idtoprot[num] for num in new_nums]
        # print(new_prots)
        row['shufSNRs'].append(SNR(new_prots))
        row['shufcophs'].append(coph(new_prots))

"""
cophs is a list of pairwise cophs
shufcophs is a list of lists of pairwise cophs
SNR is a groupwise SNR
shufSNRs is a list of groupwise SNRs
"""

# Add a new column 'avg_shufcophs' for the average of the lists in 'shufcophs'
stats_df['avg_shufcophs'] = stats_df['shufcophs'].apply(np.mean)

# Add a new column 'avg_shufSNRs' for the average of the lists in 'shufSNRs'
stats_df['avg_shufSNRs'] = stats_df['shufSNRs'].apply(np.mean)

#stats_df.to_csv('spotcheck.tsv', sep='\t', index = True)

"""T TESTS"""
for index, row in stats_df.iterrows():
    t_coph, p_coph = stats.ttest_1samp(row['shufcophs'], row['coph'])
    
    # Adjust p-value based on the t-value
    if t_coph > 0:
        p_coph = p_coph / 2
    else:
        p_coph = 1 - (p_coph / 2)
    
    # Assign results back to the DataFrame
    stats_df.at[index, 't_coph'] = t_coph
    stats_df.at[index, 'p_coph'] = p_coph
    
    # Calculate t_snr and p_snr and assign back to the DataFrame
    t_snr, p_snr = stats.ttest_1samp(row['shufSNRs'], row['SNR'])

    if t_snr < 0:
        p_snr = p_snr / 2
    else:
        p_snr = 1 - (p_snr / 2)

    stats_df.at[index, 't_snr'] = t_snr
    stats_df.at[index, 'p_snr'] = p_snr

#print(stats_df)



final_df = stats_df[['coph', 'avg_shufcophs', 't_coph', 'p_coph', 'SNR', 'avg_shufSNRs', 't_snr', 'p_snr']].copy()
final_df['proteins'] = stats_df.index
final_df['coph_enrichment'] = final_df.apply(lambda row: row['coph']/row['avg_shufcophs'], axis=1)
final_df['SNR_enrichment'] = final_df.apply(lambda row: row['SNR']/row['avg_shufSNRs'], axis=1)
final_df = final_df[['proteins', 'coph', 'avg_shufcophs', 'coph_enrichment', 't_coph', 'p_coph', 'SNR', 'avg_shufSNRs', 'SNR_enrichment', 't_snr', 'p_snr']]
final_df=final_df.sort_values(by='SNR', ascending=False)
final_df.to_csv('shuffle-STRING-group-stats.tsv', sep='\t', index=False)
#stats_df.to_csv('extended-shuffle-STRING-group-stats.tsv', sep='\t', index=True)