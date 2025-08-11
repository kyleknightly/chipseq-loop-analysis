import networkx as nx
import pandas as pd

G = nx.DiGraph()
skipcounter = 0
df = pd.read_csv('/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/beds/hepg2-dELS-loop-types.bedpe',
                sep="\t", header=None, 
                names = ['chr1', 'start1', 'end1', 'type1', 'chr2', 'start2', 'end2', 'type2'])


for _, row in df.iterrows():
    node1 = (row['chr1'], row['start1'], row['end1'])
    node2 = (row['chr2'], row['start2'], row['end2'])
    t1, t2 = row['type1'], row['type2']

    # enhancer → promoter
    if 'ELS' in t1 and 'PLS' in t2:
        G.add_edge(node1, node2)
    elif 'PLS' in t1 and 'ELS' in t2:
        G.add_edge(node2, node1)
    # optionally: promoter → promoter or enhancer → enhancer
    else:
        skipcounter+=1  # skip or G.add_edge(node1, node2) if desired

print(skipcounter)
print(skipcounter/len(df))

nx.write_gpickle(G, 'directed-dEP-graph.gpickle')