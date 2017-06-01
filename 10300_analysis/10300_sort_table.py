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
from operator import itemgetter
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--gene_table',required=True,type=str,help='gene table as output by 10300_gene_tables.py')
conf = ap.parse_args()

with open(conf.gene_table) as f:
    table_lines = f.readlines()

fpkm_set = set()
fpkm_dict = defaultdict(list)

for line in table_lines[1:]:
    line = line.rstrip()
    split_line = line.split("\t")
    fpkm = split_line[19]
    fpkm_dict[fpkm].append(line)
    # print fpkm
    fpkm_set.add(fpkm)

fpkm_sorted = sorted(fpkm_set, key=float, reverse=True)

print "\t".join(["rank", table_lines[0].rstrip()])

i = 1
for fpkm in fpkm_sorted:
    outlines = fpkm_dict[fpkm]
    for line in outlines:
        rank = i
        # print str(rank) + "\t" + str(fpkm) + "\t" + line
        print "\t".join([str(rank), line])
    i += len(outlines)
