#!/usr/bin/python

'''
Extraction of upstream regions from gene starts. Distance upstream can be set.
'''

import argparse
from collections import defaultdict



# -----------------------------------------------------
# Step 1
# Import variables, load input files & create set of genes
# If using a different number of files, arguments & appending to list of genes
# will need to be changed
# -----------------------------------------------------

ap = argparse.ArgumentParser()
ap.add_argument('--fasta', required=True, type=str, help='Fasta file of assembly contigs')
ap.add_argument('--gff', required=True, type=str, help='Gene model Gff')
conf = ap.parse_args()

# -----------------------------------------------------
# Step 2
#
# -----------------------------------------------------


with open(conf.fasta) as f:
    gff_lines = f.readlines()

with open(conf.gff) as f:
    gff_lines = f.readlines()
