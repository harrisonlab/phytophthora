#!/usr/bin/python

'''
This script uses text files of differentially expressed genes to create a
table for plotting venn diagrams.
'''

import argparse
from collections import defaultdict
import numpy
import csv

# -----------------------------------------------------
# Step 1
# Import variables, load input files & create set of genes
# If using a different number of files, arguments & appending to list of genes
# will need to be changed
# -----------------------------------------------------

ap = argparse.ArgumentParser()
ap.add_argument('--input', required=True, type=str, help='File containing predicted fpkm values for libraries')
# ap.add_argument('--conditions', required=True, type=str, nargs='+', help='list of conditions relating to input files')
#ap.add_argument('--out_dir', required=True, type=str, help='the tsv file where the count table is output to')
conf = ap.parse_args()

#-----------------------------------------------------
# Step 2
# Define classes
#-----------------------------------------------------

class fpkm_obj(object):
    """A gene identified as differentially expressed in one isolate.
    Attributes:
        gene_name: A string representing the gene name.
        conditions_tested: List of conditions that the gene was tested for DE.
        conditions_positive: A binary list of integers representing whether the
            gene tested positive for each condition listed in conditions tested.
    """

    def __init__(self, gene_name, treatments_all = [], fpkm_values = [], treatments_grouped = [], fpkm_grouped = []):
        """Return a fpkm_obj whose name is *gene_name*"""
        self.gene_name = gene_name
        self.treatments_all = treatments_all
        self.fpkm_values = fpkm_values
        self.treatments_grouped = treatments_grouped
        self.fpkm_grouped = fpkm_grouped

    def fill_obj(self, gene_line, treatment_list):
        """Reset conditions_tested to a list of conditions"""
        self.treatments_all = treatment_list
        split_line = gene_line.split()
        self.gene_name = split_line[0]
        self.fpkm_values = split_line[1:]
        fpkm_dict = defaultdict(list)
        for treatment, fpkm in zip(treatment_list, split_line[1:]):
            fpkm_dict[treatment].append(float(fpkm))
        # print fpkm_dict
        treatments_grouped = fpkm_dict.keys()
        self.treatments_grouped = treatments_grouped
        fpkm_list = []
        for key in treatments_grouped:
            fpkm_by_treatment = fpkm_dict[key]
            # print key
            # print fpkm_by_treatment
            # print numpy.mean(fpkm_by_treatment)
            if 5 < numpy.mean(fpkm_by_treatment):
                fpkm_list.append('1')
            else:
                fpkm_list.append('0')
        self.fpkm_grouped = fpkm_list
        # print self.fpkm_grouped

#-----------------------------------------------------
# Step 3
# Build evidence of expression gene table for each condition
#-----------------------------------------------------

gene_dict = defaultdict(list)

with open(conf.input) as f:
    lines = f.readlines()
    first_line = lines[0].rstrip()
    conditions_list = first_line.split("\t")
    # print conditions_list
    for line in lines[1:]:
        split_line = line.split()
        gene_name = split_line[0]
        # print gene_name
        gene_obj = fpkm_obj(gene_name)
        gene_obj.fill_obj(line, conditions_list)
        gene_dict[gene_name] = gene_obj
        # print gene_obj

gene_list = sorted(gene_dict.keys(), key=lambda x:int(x[1:]))
treatments = [x.replace(' h', '_h').replace(' - ', '_') for x in gene_obj.treatments_grouped]
print "\t".join(['Gene'] + treatments)

for gene in gene_list:
    outline = gene_dict[gene].fpkm_grouped
    print "\t".join([gene] + outline)
