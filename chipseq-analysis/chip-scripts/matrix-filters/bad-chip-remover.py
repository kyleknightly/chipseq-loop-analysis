"""
Takes a matrix (heatmap of protein interactions of some sort) and removes proteins
with <peak_filter peaks in it's file.

This iterates through a chipdir, which holds the peaks for each protein
It expects the names to be "PROTEIN-whatever.bed"
"""

import os
import pandas as pd

peak_filter=250
chipdir = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/chipbin10beds_merged'
peak_counts = {}
for bed_filename in os.listdir(chipdir):
    if bed_filename.endswith('.bed'): #verify filetype
        bed_filepath = os.path.join(chipdir, bed_filename) #get path
        with open(bed_filepath, 'r') as bed_file: #open
            num = sum(1 for line in bed_file)
            peak_counts[bed_filename.split('-')[0]] = num

# remove = []
# for key, value in peak_counts.items():
#     if value <= peak_filter:
#         remove.append(key)

# print(remove)
# print(len(remove))
remove = ['RBAK', 'NPAS2', 'ZNF484', 'HOXA7', 'ZNF17', 'ZNF660', 'LRRFIP1', 'ZNF530', 'ZNF512', 'ZNF490', 'ZNF260', 'ZBTB24', 'NCOA5', 'ZNF707', 'ZNF326', 'GLMP', 'PAF1', 'ZMAT5', 'ZNF432', 'MATR3', 'ZNF10', 'ATRX', 'ZNF382', 'HMGA1', 'GPBP1L1', 'SCMH1', 'CENPX', 'ZNF496', 'GRHL1', 'ZNF562', 'MTA3', 'NR3C1', 'MLLT10', 'ZUFSP', 'ZNF146', 'SETDB1', 'NFE2', 'DDIT3', 'ZNF7', 'MAFG', 'ZC3H4', 'ZNF329', 'NFIB', 'ZBTB1', 'ZFP14', 'ZNF577', 'JDP2', 'ZNF75D', 'ATAD3A', 'CSRNP2', 'ZNF549', 'ZNF841', 'ZSCAN22', 'ZNF597', 'AKAP8L', 'CERS6', 'ZNF180', 'ZNF720', 'PRDM4', 'KLF16', 'PLSCR1', 'ZFP62', 'ZFAT', 'ZNF850', 'ZNF672', 'ZNF737', 'TEAD2', 'GPN1', 'ZBED5', 'PREB', 'TERF1', 'REXO4', 'MIER2', 'MTF2', 'ZKSCAN5', 'CBX2', 'ZNF20', 'CENPBD1', 'ZNF251', 'CREM', 'ELK4', 'FLYWCH1', 'PIN1', 'ZNF296', 'ZBTB49', 'AKNA', 'ZFP36L2', 'NFE2L1', 'ZNF778', 'IRF1', 'ZNF101', 'ZFYVE20', 'MTF1', 'NR1H2', 'ZNF121', 'ZNF576', 'SSRP1', 'THYN1', 'HSF2', 'ZNF460', 'ZFP64', 'ADNP', 'ZNF629', 'ZNF564', 'XBP1', 'AFF4', 'ZNF646', 'RORA', 'TRIM24', 'ZNF865', 'ZNF3', 'ATF4', 'ZNF317', 'ZBTB44', 'DBP', 'ISX', 'ETS1', 'ZNF468', 'PAWR', 'KLF13', 'ZNF780A', 'BHLHA15', 'MLXIP', 'DMTF1', 'ZNF160', 'ZNF550', 'MZF1', 'ZNF367', 'ZNF30', 'SP110', 'PBX3', 'TBP', 'ZNF513', 'ZC3H13', 'ZBTB46', 'STAT5B', 'ZBTB8A', 'RARG', 'KIAA2018', 'FUBP1', 'GZF1', 'FOXO4', 'TEF', 'TP53', 'ZNF615', 'ZSCAN30', 'CC2D1A', 'PRMT3', 'ZNF775', 'ZNF12', 'ZNF33A', 'JUNB', 'ZBTB37', 'ONECUT2', 'RNF219', 'TCF3', 'FOSL1', 'ZNF703', 'KDM1A', 'ZNF652', 'ZBTB3', 'AKAP8', 'WIZ', 'THAP8', 'BATF2', 'TIGD3', 'MYRF', 'ZHX2', 'MYNN', 'ZNF713', 'ZNF485', 'ZBTB33', 'ATF6', 'ZSCAN12', 'ATF5', 'ZNF234', 'ZNF674', 'ZFP41', 'ZNF253', 'ZNF768', 'RFX3', 'ZNF800', 'ZNF431', 'TSC22D1', 'HOXA9', 'SMAD9', 'ZBTB40', 'ZNF773', 'CENPT', 'TRAFD1', 'ZNF569', 'ZBTB7A', 'KLF9', 'ETV6', 'ZNF527', 'EEA1', 'DNMT1', 'HES4', 'ZNF276', 'ZNF280B', 'ZNF235', 'HIC2', 'ZNF512B', 'ZFP1', 'SOX18', 'ZBTB4', 'FUBP3', 'IRX3', 'ZNF589', 'ZNF790', 'ZNF558', 'REL', 'ZSCAN20', 'ZNF777', 'ZNF263', 'ZBTB21', 'USF2', 'NAIF1', 'TOE1', 'SNAPC4', 'EED', 'ZNF526', 'ZNF740', 'ZBTB14', 'GMEB2', 'E2F5', 'ZNF552', 'ZNF138', 'ZCCHC11', 'ZNF839', 'ZNF616', 'PAX8', 'PBX2', 'NRF1', 'ZFP37', 'ZNF441', 'NFATC3', 'ZBTB10', 'ZHX3', 'CREBL2', 'RFXANK', 'CAMTA2', 'GTF3A', 'E2F2', 'MAF1', 'ZNF571', 'ZNF318', 'LBX2', 'KLF15', 'ZNF83', 'ZNF256', 'ZNF827', 'CHCHD3', 'CREB3', 'NFAT5', 'THAP4', 'ZNF678', 'JARID2', 'DZIP1', 'ZXDC', 'ZNF879', 'ZNF25', 'ZNF446', 'ZNF383', 'RELA', 'ZNF770', 'PRDM15', 'ZNF572', 'ZNF354B', 'ZNF18', 'ZBTB42', 'SATB2', 'SNAPC5', 'SP140L', 'ZNF557', 'NRL', 'ZNF891', 'ZNF619', 'SMAD1', 'HINFP', 'ATF1', 'ZNF232', 'ZNF781', 'MTERF4', 'ZNF548', 'ZNF343', 'ZNF761', 'KLF11', 'TSC22D2', 'ZBTB2', 'SAFB2', 'AHR', 'ZNF766', 'NR0B2', 'ZNF776', 'ZNF275', 'ZNF451', 'ZNF563', 'E2F1', 'TIGD6', 'ZNF143', 'ZNF335', 'BRF2', 'NKX3', 'ZMYM3', 'ZNF281', 'PRDM10', 'ZNF33B', 'ZNF786', 'THAP7', 'ZNF749', 'ZNF34', 'ZNF724P', 'DNMT3B', 'ZNF608', 'ZNF181', 'JRK', 'ELF4', 'ZNF44', 'ZC3H8', 'ZNF782', 'ZNF510', 'ZNF337', 'BRCA1', 'ZNF142', 'MTERF2', 'ZFP82', 'ZNF224', 'SNAI1', 'ZNF697', 'ZNF333', 'TMF1', 'ZNF567', 'ZNF225', 'GLI4', 'ZNF704', 'SP2', 'ZZZ3', 'ELK1', 'ZIK1', 'IKZF4', 'SRY', 'ARHGAP35', 'IRF5', 'ZNF546', 'JUND', 'THRB', 'ESRRA', 'NR2F1', 'ZMYM2', 'HNF4G', 'ZNF816', 'ZBTB26', 'SRF', 'ZSCAN25', 'ZNF746', 'ZBTB39', 'ZFHX3', 'MEF2A', 'HHEX', 'NCOA1', 'SKIL', 'JUN', 'ZBTB34', 'SIX4', 'FOXQ1', 'ONECUT1', 'SIX1', 'STAT6', 'SMYD3', 'ATF3', 'ZSCAN29', 'ZNF264', 'MED13', 'HMG20B', 'ZNF878', 'NFIA', 'CEBPD', 'MED8', 'MLX', 'TCF12', 'ZBTB43', 'ZNF136', 'ARID2', 'SNAPC2', 'KDM4B', 'BRD4', 'ZHX1', 'ZNF709', 'TOPORS', 'HOXD1', 'FBXL19', 'ZFP36L1', 'ETV4', 'ZNF670', 'ZSCAN5A', 'ZNF362', 'TEAD4']

matrix = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/pseudo-trans-contacts.tsv'
df = pd.read_csv(matrix, sep='\t', header=0, index_col=0)
remove_rows = [item for item in remove if item in df.index]
remove_columns = [item for item in remove if item in df.columns]

# Remove the specified rows and columns
df = df.drop(index=remove_rows, columns=remove_columns)

name = os.path.basename(matrix)
df.to_csv('subset-' + name, sep='\t', index=True)
