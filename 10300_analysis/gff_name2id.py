#!/usr/bin/python

'''
This program parses information from fasta files and gff files for the location,
sequence and functional information for annotated gene models and RxLRs.
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument('--gff',required=True,type=str,help='Gff file of predicyted gene models')


conf = ap.parse_args()



with open(conf.gff) as f:
    gene_lines = f.readlines()


for line in gene_lines:
    line = line.rstrip()
    if line.startswith('#'):
        continue
    split_line = line.split("\t")
    if 'gene' in split_line[2]:
        # Col9 = split_line[8]
        Gene_split_line = split_line
    elif 'transcript' in split_line[2]:
        Col9 = split_line[8]
        Name = Col9.split(";")[2]
        ID = Name.replace('Name=', '')
        Transcript_split_line = split_line
        Gene_split_line[8] = 'ID=' + ID
        Transcript_split_line[8] = 'ID=' + ID + '.t1;Parent=' + ID
        print("\t".join(Gene_split_line))
        print("\t".join(Transcript_split_line))
