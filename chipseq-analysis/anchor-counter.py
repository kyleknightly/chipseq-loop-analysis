import pandas as pd
file = '/mnt/altnas/work/Kyle.Knightly/hepg2-ctcfanchors-intallbin10_sortedmerged_wchip3_catsort450bp-wawb.bed'
df = pd.read_csv(file, sep = '\t', header = None, names = ['chr', 'start', 'end', 'peak_chr', 'peak_start', 'peak_end'])
print(df)
unique_count = len(df[['chr', 'start', 'end']].drop_duplicates())
print(unique_count)

# file = '/mnt/altnas/work/Kyle.Knightly/S5_HepG2.attributes_anchorAPA_pass.bed_sieveML_only_good.bedpe'

# df = pd.read_csv(file, sep='\t', header = 0)
# print(df)