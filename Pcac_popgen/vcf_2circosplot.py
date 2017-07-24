#!/usr/bin/python

'''
Program to convert annotated vcf files of SNPs
into a text file that can be used to generate a scatterpolot
'''

import sys,argparse
from collections import defaultdict
from sets import Set

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument('--inp_vcf',required=True,type=str,help='input vcf file')
conf = ap.parse_args()

with open(conf.inp_vcf) as f:
    inp_lines = f.readlines()

#-----------------------------------------------------
# Step 2
# print locations of SNPs as circos scatterplot values
#-----------------------------------------------------

contig_list = []

for line in inp_lines:
    line = line.strip("\n")
    if line.startswith("#"):
        continue
    split_line = line.split()
    contig = split_line[0]
    SNP_bp = split_line[1]
    print("\t".join([contig, SNP_bp, SNP_bp ,"1"]))
