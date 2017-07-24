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
# ap.add_argument('--gff_windows',required=False,type=str,help='gff file giving genomic windows over which SNPs will be recorded')
# ap.add_argument('--per_X_bp',required=False, default=1000, type=float, help='The number of bp over which you want to report the number of features eg. SNPs per 1000bp')


conf = ap.parse_args() #sys.argv
# multiplier = conf.per_X_bp
#
# if conf.gff_windows:
#     windows = conf.conf.gff_windows

with open(conf.inp_vcf) as f:
    inp_lines = f.readlines()


#-----------------------------------------------------
# Step 2
# Order gff features by input contigs and gene start
#-----------------------------------------------------

contig_list = []
# gene_start_dict = defaultdict(list)
# features_dict = defaultdict(list)
# conversion_dict = defaultdict(list)

for line in inp_lines:
    line = line.strip("\n")
    if line.startswith("#"):
        # print line
        continue
    split_line = line.split()
    contig = split_line[0]
    SNP_bp = split_line[1]
    print("\t".join([contig, SNP_bp, SNP_bp ,"1"]))
