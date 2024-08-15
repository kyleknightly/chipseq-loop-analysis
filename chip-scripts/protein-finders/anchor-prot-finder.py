"""
Iterating through a directory of chipseq data and a list of anchors, this finds all anchors
with a given protein and creates a new bed file"""

from collections import defaultdict
import os
import csv

chipdir = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/chipbin10beds_merged'
anchors = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/hepg2-anchors.bed'

def lookup_dic(chipdir):
    interval_origins = defaultdict(list) #the dic
    for bed_filename in os.listdir(chipdir):
        if bed_filename.endswith('.bed'): #verify filetype
            bed_filepath = os.path.join(chipdir, bed_filename) #get path
            with open(bed_filepath, 'r') as bed_file: #open
                reader = csv.reader(bed_file, delimiter='\t') #make csv
                for row in reader: #parse
                    interval = (row[0], row[1], row[2])  # (ch2, start2, end2) in origin's inbed var
                    interval_origins[interval].append(bed_filename.split('-')[0])
    #print(interval_origins)
    print("made lookup_dic")
    return interval_origins

def find_overlaps(interval, lookup_dict):
    overlap_values = []
    chrom, start, end = interval
    for (lookup_chrom, lookup_start, lookup_end), values in lookup_dict.items():
        if chrom == lookup_chrom and int(start) <= int(lookup_end) and int(end) >= int(lookup_start):
            overlap_values.extend(values)

    return list(set(overlap_values))  # Remove duplicates


def origin(inbed, outbed, chipdir):
    interval_prots = lookup_dic(chipdir)
    with open(inbed, 'r') as infile, open(outbed, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        seen = {} #makes this more efficient since often anchors are seen multiple times
        
        for row in reader:
            
            interval1 = (row[0], row[1], row[2])  # (ch1, start1, end1)
            if interval1 in seen:
                prots1 = seen[interval1]
                #print(interval1)
            else:
                prots1 = find_overlaps(interval1, interval_prots)
                seen[interval1] = prots1
            

            writer.writerow(row[:3] + [prots1])
            #print(row[:3] + [prots1])

            

origin(anchors, 'hepg2-anchor-TFs.bed',chipdir)