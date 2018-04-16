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
ap.add_argument('--input', required=True, type=str, nargs='+', help='list of genes files containing differentially expressed genes at timepoints')
ap.add_argument('--conditions', required=True, type=str, nargs='+', help='list of conditions relating to input files')
#ap.add_argument('--out_dir', required=True, type=str, help='the tsv file where the count table is output to')
conf = ap.parse_args()

input_list = []
input_list = conf.input
conditions_list = []
conditions_list = conf.conditions

#-----------------------------------------------------
# Step 2
# Define classes
#-----------------------------------------------------

class DEG_obj(object):
    """A gene identified as differentially expressed in one isolate.
    Attributes:
        gene_name: A string representing the gene name.
        conditions_tested: List of conditions that the gene was tested for DE.
        conditions_positive: A binary list of integers representing whether the
            gene tested positive for each condition listed in conditions tested.
    """

    def __init__(self, gene_name, conditions_tested = [], conditions_positive = []):
        """Return a DEG_objs whose name is *gene_name*, with no conditions or positive results"""
        self.gene_name = gene_name
        self.conditions_tested = conditions_tested
        self.conditions_positive = conditions_positive

    def set_conditions(self, c_list):
        """Reset conditions_tested to a list of conditions"""
        self.conditions_tested = c_list

    def set_positive(self, c_list):
        """Reset conditions_positive to a binary list reflecting which conditions were positive"""
        p_list = []
        for x in c_list:
            p_list.append('0')
        self.conditions_positive = p_list

    def add_positive(self, condition):
        """Edit the binary list """
        c_list = self.conditions_tested
        p_list = self.conditions_positive
        i = c_list.index(condition)
        p_list[i] = '1'
        self.conditions_positive = p_list

#-----------------------------------------------------
# Step 3
# Build DEG gene table for each condition
#-----------------------------------------------------

gene_dict = defaultdict(list)
for file, condition in zip(input_list, conditions_list):
    with open(file) as f:
        lines = f.readlines()
        for x in lines:
            gene_name = x.rstrip()
            if gene_dict[gene_name]:
                gene_dict[gene_name].add_positive(condition)
            else:
                gene_obj = DEG_obj(gene_name)
                gene_obj.set_conditions(conditions_list)
                gene_obj.set_positive(conditions_list)
                gene_obj.add_positive(condition)
                gene_dict[gene_name] = gene_obj

gene_list = sorted(gene_dict.keys(), key=lambda x:int(x[1:]))
print "\t".join(['Gene'] + conditions_list)
for gene in gene_list:
    outline = gene_dict[gene].conditions_positive
    print "\t".join([gene] + outline)
