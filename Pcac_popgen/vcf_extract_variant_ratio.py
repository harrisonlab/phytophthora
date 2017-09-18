#!/usr/bin/python

'''
Extract the variant ratio for 0/1 SNPs in a user defined isolate
'''

import sys,argparse
from collections import defaultdict
from sets import Set
import numpy as np

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

print("\t".join(["GT", "AD", "GQ", "AD_ratio"]))

for line in inp_lines:
    line = line.rstrip()
    if line.startswith('##'):
        continue
    elif line.startswith('#CHROM'):
        split_line = line.split("\t")
        ref_index = split_line.index(ref_isolate)
        continue
    split_line = line.split("\t")
    contig = split_line[0]
    location = split_line[1]
    ref_SNP = split_line[ref_index]
    GT = ref_SNP.split(':')[0]
    AD = ref_SNP.split(':')[1] # Also known as the allele depth, this annotation gives the unfiltered count of reads that support a given allele for an individual sample.
    DP = ref_SNP.split(':')[2] # At the sample level (FORMAT), the DP value is the count of reads that passed the caller's internal quality control metrics (such as MAPQ > 17, for example). At the site level (INFO), the DP value is the unfiltered depth over all samples.
    GQ = ref_SNP.split(':')[3] # Genotype qulaity
    if GT == '0/1' or GT == '0|1':
        AD_split = AD.split(",")
        # AD_ratio = np.divide(float(AD_split[0]), float(AD_split[1]))
        AD_list = [float(i) for i in AD_split]
        AD_ratio = np.divide(min(AD_list), max(AD_list))
        AD_ratio = AD_ratio.round(2)
        print("\t".join([GT, AD, GQ, str(AD_ratio)]))
        # DP_split = DP.split(",")
        # DP_ratio = np.divide(int(DP_split[0]), int(DP_split[1]))
        # DP_ratio = DP_ratio.round(2)
        # print("\t".join([AD, GT, str(AD_ratio), str(DP_ratio)]))
