#!/usr/bin/python

'''
Program to flag SNPs in a reference genome that show 1/1 SNPs from its own illumina
reads, indicating that these SNPs represent errors within its own genome (a 1/0)
SNP may reflect its diploid nature.
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
ap.add_argument('--ref_isolate',required=True,type=str,help='The isolate id in the input vcf that should be investigated')
conf = ap.parse_args()

ref_isolate = conf.ref_isolate

with open(conf.inp_vcf) as f:
    inp_lines = f.readlines()


#-----------------------------------------------------
# Step 2
#
#-----------------------------------------------------

for line in inp_lines:
    line = line.rstrip()
    if line.startswith('##'):
        continue
    elif line.startswith('#CHROM'):
        split_line = line.split("\t")
        ref_index = split_line.index(ref_isolate)
        # print isolate_index
        continue
    split_line = line.split("\t")
    contig = split_line[0]
    location = split_line[1]
    ref_SNP = split_line[ref_index]
    GT = ref_SNP.split(':')[0]
    if GT == '1/1' or GT == '1|1':
        print "\t".join([contig, location, ref_SNP])
        # print "\t".join([contig, location])
