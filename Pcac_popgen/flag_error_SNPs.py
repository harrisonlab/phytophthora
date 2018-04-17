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
ap.add_argument('--errors',required=True,type=str,help='tsv file of contig and location of error SNPs')
ap.add_argument('--filtered',required=True,type=str,help='filtered vcf file with error SNP lines removed')
ap.add_argument('--ploidy',required=True,type=str,help='ploidy of genotypes (GT) in the vcf file: haploid or diploid')


conf = ap.parse_args()

ref_isolate = conf.ref_isolate
f_errors = conf.errors
f_filtered = conf.filtered
ploidy = conf.ploidy
if not any(ploidy != x for x in ['haploid', 'diploid']):
    sys.exit("ploidy error")

with open(conf.inp_vcf) as f:
    inp_lines = f.readlines()


#-----------------------------------------------------
# Step 2
#
#-----------------------------------------------------

errors_out = []
filtered_out = []

for line in inp_lines:
    line = line.rstrip()
    if line.startswith('##'):
        filtered_out.append(line)
        continue
    elif line.startswith('#CHROM'):
        split_line = line.split("\t")
        ref_index = split_line.index(ref_isolate)
        filtered_out.append(line)
        continue
    split_line = line.split("\t")
    contig = split_line[0]
    location = split_line[1]
    ref_SNP = split_line[ref_index]
    GT = ref_SNP.split(':')[0]
    if ploidy == 'diploid':
        if GT == '1/1' or GT == '1|1':
            errors_out.append("\t".join([contig, location]))
        else:
            filtered_out.append(line)
    elif ploidy == 'haploid':
        if GT == '1':
            errors_out.append("\t".join([contig, location]))
        else:
            filtered_out.append(line)
    else:
        print "ploidy error"
        exit()

f = open(f_errors,"w")
f.write("\n".join(errors_out) + "\n")
f.close()

f = open(f_filtered,"w")
f.write("\n".join(filtered_out) + "\n")
f.close()
