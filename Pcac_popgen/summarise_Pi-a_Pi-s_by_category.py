#!/usr/bin/python

'''
Program to summarise Pi and Ps scores by effector categories as described in a
supplied annotation table.
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import os
import argparse
import re
import math
import numpy as np
from sets import Set
from collections import defaultdict
from collections import Counter
from operator import itemgetter


ap = argparse.ArgumentParser()

ap.add_argument('--annotation_table',required=True,type=str,help='P414 annotation table')
ap.add_argument('--Pcac_Fa_Pi_stats',required=True,type=str,help='P. cactorum ex. strawberry Pi(a)/Pi(s) stats')
ap.add_argument('--Pcac_Md_Pi_stats',required=True,type=str,help='P. cactorum ex. apple Pi(a)/Pi(s) stats')

conf = ap.parse_args()


#-----------------------------------------------------
# Step X
# Define classes
#-----------------------------------------------------

class VariantsObj(object):
    """
        A dictionary of GeneObj and methods to summarise popgen stats
    """
    def __init__(self):
        """set the object, defining the dictionary."""
        dict = defaultdict(list)
        self.gene_family_dict = dict
        self.gene_family_dict['Genes'] = PiObj()
        self.gene_family_dict['RxLR'] = PiObj()
        self.gene_family_dict['CRN'] = PiObj()
        self.gene_family_dict['CAZY'] = PiObj()
        self.gene_family_dict['Apoplastic:Cutinase'] = PiObj()
        self.gene_family_dict['Apoplastic:Glucanase inhibitor'] = PiObj()
        self.gene_family_dict['Apoplastic:NLP'] = PiObj()
        self.gene_family_dict['Apoplastic:Phytotoxin'] = PiObj()
        self.gene_family_dict['Apoplastic:Protease inhibitor (cathepsin)'] = PiObj()
        self.gene_family_dict['Apoplastic:Protease inhibitor (cystatin-like)'] = PiObj()
        self.gene_family_dict['Apoplastic:Protease inhibitor (Kazal-type)'] = PiObj()
        self.gene_family_dict['MAMP:Elicitin'] = PiObj()
        self.gene_family_dict['MAMP:Transglutaminase'] = PiObj()

class PiObj(object):
    """
        An object containing Popgen stats and effector status for a given gene.
    """
    def __init__(self):
        """set the object, defining the dictionary."""
        # dict = defaultdict(lambda: GeneObj())
        # self.gene_dict = dict
        dict = defaultdict(str)
        self.gene_dict = dict
    def add_gene(self, gene_id, Pi):
        # print Pi
        self.gene_dict[gene_id] = str(Pi)
        # print self.gene_dict[gene_id]
    def get_Pi(self):
        "Get Pi(a)/Pi(s) stats for this gene subset"
        subset_Pi = []
        for key in self.gene_dict.keys():
            print key
            # self.gene_dict[key].
            Pi = self.gene_dict[key]
            # print Pi
            if any(str(Pi) not in x for x in ['Inf', 'NaN']):
                continue
            print Pi
            subset_Pi.append(Pi)
        return subset_Pi
    def get_num_genes(self):
        "Get Pi(a)/Pi(s) stats for this gene subset"
        subset_Pi = []
        for key in self.gene_dict.keys():
            # self.gene_dict[key].
            subset_Pi.append(self.gene_dict[key])
        return len(subset_Pi)

#-----------------------------------------------------
# Step X
# Read input files
#-----------------------------------------------------

with open(conf.Pcac_Fa_Pi_stats) as f:
    Pc_Fa_lines = f.readlines()

with open(conf.Pcac_Md_Pi_stats) as f:
    Pc_Md_lines = f.readlines()

# with open(conf.Pcac_Md_Pi_stats) as f:
#     Pc_Md_lines = f.readlines()

with open(conf.annotation_table) as f:
    a_lines = f.readlines()


#-----------------------------------------------------
# Step X
# Load Pi information
#-----------------------------------------------------

Pc_Fa_dict = defaultdict(str)
for line in Pc_Fa_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[1]
    gene_id = gene_id.replace('ID=', '').replace(';', '')
    Pi = split_line[2]
    # print Pi
    Pc_Fa_dict[gene_id] = str(Pi)

Pc_Md_dict = defaultdict(str)
for line in Pc_Md_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[1]
    gene_id = gene_id.replace('ID=', '').replace(';', '')
    Pi = split_line[2]
    Pc_Md_dict[gene_id] = str(Pi)

#-----------------------------------------------------
# Step X
# Identify variants
#-----------------------------------------------------

Pc_Fa_GenesObj = VariantsObj()
Pc_Md_GenesObj = VariantsObj()

# SNP_genes_obj = TaxaVarObj()

for line in a_lines[1:]:
    line = line.rstrip()
    split_line = line.split('\t')
    gene_id = split_line[0]
    gene_id = gene_id.split('.')[0]
    secreted_col = split_line[10]
    RxLR_col = split_line[13]
    CRN_col = split_line[16]
    CAZY_col = split_line[18]
    apoplast_col = split_line[26]
    score_Pc_Fa = Pc_Fa_dict[gene_id]
    score_Pc_Md = Pc_Md_dict[gene_id]
    # print score_Pc_Fa
    Pc_Fa_GenesObj.gene_family_dict['Genes'].add_gene(gene_id, score_Pc_Fa)
    Pc_Md_GenesObj.gene_family_dict['Genes'].add_gene(gene_id, score_Pc_Md)
    if RxLR_col:
        Pc_Fa_GenesObj.gene_family_dict['RxLR'].add_gene(gene_id, score_Pc_Fa)
        Pc_Md_GenesObj.gene_family_dict['RxLR'].add_gene(gene_id, score_Pc_Md)
    if CRN_col:
        Pc_Fa_GenesObj.gene_family_dict['CRN'].add_gene(gene_id, score_Pc_Fa)
        Pc_Md_GenesObj.gene_family_dict['CRN'].add_gene(gene_id, score_Pc_Md)
    if all([secreted_col, CAZY_col]):
        Pc_Fa_GenesObj.gene_family_dict['CAZY'].add_gene(gene_id, score_Pc_Fa)
        Pc_Md_GenesObj.gene_family_dict['CAZY'].add_gene(gene_id, score_Pc_Md)
    if all([secreted_col, apoplast_col]):
        Pc_Fa_GenesObj.gene_family_dict[apoplast_col].add_gene(gene_id, score_Pc_Fa)
        Pc_Md_GenesObj.gene_family_dict[apoplast_col].add_gene(gene_id, score_Pc_Md)

family_list = [
    'Genes',
    'RxLR',
    'CRN',
    'CAZY',
    'Apoplastic:Cutinase',
    'Apoplastic:Glucanase inhibitor',
    'Apoplastic:NLP',
    'Apoplastic:Phytotoxin',
    'Apoplastic:Protease inhibitor (cathepsin)',
    'Apoplastic:Protease inhibitor (cystatin-like)',
    'Apoplastic:Protease inhibitor (Kazal-type)',
    'MAMP:Elicitin',
    'MAMP:Transglutaminase']

for Obj in [Pc_Fa_GenesObj, Pc_Md_GenesObj]:
    for family in family_list:
        total_genes = Obj.gene_family_dict[family].get_num_genes()
        genes_Pi = len(Obj.gene_family_dict[family].get_Pi())
        print("\t".join([family, str(total_genes), str(genes_Pi)]))
        print Obj.gene_family_dict[family].get_Pi()
