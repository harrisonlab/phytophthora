#!/usr/bin/python

'''
Make a 3-way venn diagram from a vcf file showin SNPs present in isolates.
'''

import sys,argparse
from collections import defaultdict
from sets import Set
import numpy as np
from operator import itemgetter

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument('--vcf',required=True,type=str,help='input vcf file')
ap.add_argument('--gene_list',required=True,type=str,help='list of gene_ids to extract SNPs from')


conf = ap.parse_args()

with open(conf.vcf) as f:
    vcf_lines = f.readlines()

with open(conf.gene_list) as f:
    gene_lines = f.readlines()

#-----------------------------------------------------
# Step 2
# make a set of gene names
#-----------------------------------------------------

gene_set = set()
for line in gene_lines:
    gene = line.rstrip()
    gene_set.add(gene)


#-----------------------------------------------------
# Step 3
# extract the lines for the required genes
#-----------------------------------------------------

for line in vcf_lines:
    line = line.rstrip()
    if line.startswith('#'):
        print line
        continue
    split_line = line.split()
    info_col = split_line[7]
    gene_id = info_col.split('|')[4]
    transcript_id = info_col.split('|')[6]
    if gene_id in gene_set or transcript_id in gene_set:
        # print gene_id
        print line
