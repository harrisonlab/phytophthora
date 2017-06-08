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
ap.add_argument('--exp_rank',required=True,type=str,help='Gff file format')
ap.add_argument('--domain',required=True,type=str,help='Gff file format')


conf = ap.parse_args()


domain_dict = defaultdict(list)

with open(conf.exp_rank) as f:
    exp_rank_lines = f.readlines()

with open(conf.domain) as f:
    domain_lines = f.readlines()

for line in domain_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    domain = split_line[0]
    orthogroup = split_line[1]
    gene = split_line[2]
    if gene.startswith('Pcac'):
        transcript = gene.split('|')[1]
        domain_dict[transcript].append("\t".join([orthogroup, domain]))
    else:
        continue



for line in exp_rank_lines:
    line = line.rstrip()
    split_line = line.split()
    rank = split_line[0]
    transcript = split_line[1]
    domain_info = "\t".join(domain_dict[transcript])
    print "\t".join([rank, transcript, domain_info])
