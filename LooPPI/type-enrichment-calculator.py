#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import seaborn as sns
from tqdm.auto import tqdm
tqdm.pandas(desc="Computing enrichments")
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, optimal_leaf_ordering

indir = '/mnt/altnas/work/Kyle.Knightly/looppi/out/'
annot_loop_file = indir + 'annotated_loops.bedpe'
trans_count_file = indir + 'trans-contacts.tsv'
prot_end_props_file = indir + 'loop_end_proportions.tsv'
type_end_props_file = indir + 'end_type_proportions.tsv'
prot_type_props_file = indir + 'protein_type_proportions.tsv'
loop_type_props_file = indir + 'loop_type_proportions.tsv'


# Loading files
trans_count_long = pd.read_csv(trans_count_file, sep='\t', index_col = False, names = ['p1', 'p2', 'count'])
prot_end_props = pd.read_csv(prot_end_props_file, sep='\t', index_col=0, names = ['count', 'prop'])
prots = prot_end_props.index.tolist()
type_end_props = pd.read_csv(type_end_props_file, sep='\t', index_col=0, names = ['count', 'prop'])
types = type_end_props.index.tolist()
prot_type_props = pd.read_csv(prot_type_props_file, sep='\t', index_col=0, names = ['type', 'count', 'type_total', 'prop'])
loop_type_props = pd.read_csv(loop_type_props_file, sep='\t', index_col = False, names = ['t1', 't2', 'count', 'prop'])

prot_type_props = prot_type_props.rename_axis("protein")
prot_type_props = (prot_type_props
             .reset_index()
             .set_index(['protein','type']))['prop']

# loop_type_props = loop_type_props.set_index(['t1','t2'])['prop']
loop_type_counts = loop_type_props.set_index(['t1','t2'])['count']

annot_loops = pd.read_csv(annot_loop_file, sep='\t', names = ['c1', 's1', 'e1', 'id1', 't1', 'p1', 'c2', 's2', 'e2', 'id2', 't2', 'p2'], index_col=False)
loops = len(annot_loops)
# for col in ["p1", "p2"]:
#     annot_loops[col] = annot_loops[col].apply(
#         lambda x: [] if pd.isna(x) or x == ""
#         else [t.strip() for t in str(x).split("|") if t.strip()]
#     )


def expected(p1, p2):
    expected = 0
    for (t1, t2), count in loop_type_counts.items():
        if p1!=p2:
            prob=prot_type_props.loc[(p1,t1)]*prot_type_props.loc[(p2,t2)]+prot_type_props.loc[(p2,t1)]*prot_type_props.loc[(p1,t2)]
        else:
            prob=prot_type_props.loc[(p1,t1)]*prot_type_props.loc[(p2,t2)]
        expected+=(prob*count)
    return(expected)


def _enrich(r):
    e = expected(r['p1'], r['p2'])
    return np.nan if (e is None or e == 0) else r['count'] / e
print(trans_count_long)
trans_count_long['enrichment'] = trans_count_long.progress_apply(_enrich, axis=1)
print(trans_count_long)
trans_matrix_upper = (trans_count_long.pivot_table(index='p1', columns='p2',
                                        values='enrichment', aggfunc='sum', fill_value=0)
           .reindex(index=prots, columns=prots, fill_value=0))
trans_matrix = trans_matrix_upper.add(trans_matrix_upper.T, fill_value=0)
np.fill_diagonal(trans_matrix.values, np.diag(trans_matrix_upper.values))

print(trans_matrix)
trans_matrix.to_csv("trans_enrichment_matrix.tsv",sep="\t",index=True, header=True)