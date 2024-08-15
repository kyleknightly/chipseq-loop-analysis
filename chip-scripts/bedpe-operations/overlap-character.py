"""
This file takes in a 6 row bed file with anchors and peaks, (Ragini's file of choice)
adds row 7 with the number of peaks at each anchor
adds row 8 with the protein that caused each peak

I usually don't use this anymore and instead use paired-anchor-prot-finder.py,
but if you have Ragini's aformentioned file format, use this
"""

import csv
from collections import defaultdict
import os

#Add row 7 with anchor count

def anchor_count(inbed, outbed):
    interval_counts = defaultdict(int)
    bed_lines = []
    with open(inbed, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            # Assuming the first interval is a combination of the first three columns (chrom, start, end)
            interval = (row[0], row[1], row[2])
            interval_counts[interval] += 1
            bed_lines.append(row)

    # Write the modified lines to the new BED file with the 7th column
    with open(outbed, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        for row in bed_lines:
            interval = (row[0], row[1], row[2])
            count = interval_counts[interval]
            writer.writerow(row + [count])

#anchor_count("finalbeds/hepg2.nonctcf-wawb-withallchipbin10mergedsorted.bed","hepg2.nonctcf-wawb-withallchipbin10mergedsorted-counts.bed")
#anchor_count("finalbeds/hepg2.ctcf-wawb-withallchipbin10mergedsorted.bed","hepg2.ctcf-wawb-withallchipbin10mergedsorted-counts.bed")
#anchor_count("/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/all-anchor-peaks.bed", '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/all-anchor-peaks-counts.bed')


#Creates a lookup dictionary

def lookup_dic(chipdir):
    interval_origins = defaultdict(list) #the dic
    for bed_filename in os.listdir(chipdir):
        if bed_filename.endswith('.bed'): #verify filetype
            bed_filepath = os.path.join(chipdir, bed_filename) #get path
            with open(bed_filepath, 'r') as bed_file: #open
                reader = csv.reader(bed_file, delimiter='\t') #make csv
                for row in reader: #parse
                    interval = (row[0], row[1], row[2])  # (ch2, start2, end2) in origin's inbed var
                    interval_origins[interval].append(bed_filename)

                    
    return interval_origins

# test_dict = lookup_dic("chipseq-analysis/chipbin10beds_merged")
# for i, (key, value) in enumerate(test_dict.items()):
#     if i < 10:
#         print(f'{key}: {value}')
#     else:
#         break

def origin(inbed, outbed, chipdir):
    interval_origins = lookup_dic(chipdir)
    with open(inbed, 'r') as infile, open(outbed, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        n=0
        for row in reader:
            interval = (row[3], row[4], row[5])  # (ch2, start2, end2)
            #print(interval)
            #print(interval_origins[interval])
            if interval not in interval_origins.keys():
                n+=1
                print(n)
                #print(interval)
                #print(int(interval[2])-int(interval[1]))
            else: 
                origin = interval_origins[interval]
                originsplit = [item.split('-')[0] for item in origin] #gets rid of the file details and just leaves the protein name
                writer.writerow(row + [originsplit[0]])
                # Right now theres an issue where col8 will have 2 proteins IFF col1-7 are identical.
                # Basically, if there was already an exact col1-7 we should pop
                # *there are 2 elements in the list IFF col1-7 are duplicated, as far as I know
                if len(originsplit)>1:
                    interval_origins[interval].pop(0)

#origin("chipseq-analysis/hepg2.ctcf-wawb-withallchipbin10mergedsorted-counts.bed", "chipseq-analysis/hepg2.ctcf-wawb-withallchipbin10mergedsorted-counts-origins.bed", "chipseq-analysis/chipbin10beds_merged")

origin('/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/all-anchor-peaks-counts.bed', '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/all-anchor-peaks-origins.bed', '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/chipbin10beds_merged')