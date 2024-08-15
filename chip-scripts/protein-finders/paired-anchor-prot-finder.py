"""
This creates our paired-anchor-TF file

This takes a loop list and a directory of chipseq tracks and finds overlaps in each anchor
It adds cols after each anchor with a list of the proteins in that anchor
"""


from collections import defaultdict
import os
import csv

chipdir = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/extended-set/chipbin10beds_merged'
paired_anchors = '/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/restricted-set/paired-anchors/hepg2-loops.bed'

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
    with open(inbed, 'r') as infile, open(outbed, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter=' ')
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
            
            interval2 = (row[3], row[4], row[5])  # (ch2, start2, end2)

            if interval2 in seen:
                prots2 = seen[interval2]
                #print(interval2)
            else:
                prots2 = find_overlaps(interval2, interval_prots)
                seen[interval2] = prots2

            writer.writerow(row[:3] + [prots1] + row[3:] + [prots2])

            

origin(paired_anchors, 'paired-anchor-TFs.bed',chipdir)