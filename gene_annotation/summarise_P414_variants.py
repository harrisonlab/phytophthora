#!/usr/bin/python

'''
This program is used to summarise variant information contained within the
P414 annotation table. Showing numbers of non-synonymous SNPs, indels and SVs
by different gene families.
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

conf = ap.parse_args()


#-----------------------------------------------------
# Step X
# Define classes
#-----------------------------------------------------

class FamilyVarObj(object):
    """
        An object containing TaxaVarObj's for different gene families.
    """
    def __init__(self):
        """set the object, defining the dictionary."""
        dict = defaultdict(list)
        self.gene_family_dict = dict
    def add_family(self, family):
        "Add a TaxaVarObj to the for a given gene family"
        self.gene_family_dict[family] = TaxaVarObj()



class TaxaVarObj(object):
    """
        An object containing information on number of variants by different
        gene families.
    """
    def __init__(self):
        """set the object, defining the dictionary."""
        dict = defaultdict(list)
        self.variant_dict = dict
    def add_variant(self, gene_id, taxa):
        "Add a variant to the object"
        for strain in taxa.split(','):
            if any([strain == x for x in ['371', 'SCRP370', 'SCRP376']]):
                self.variant_dict['P.idaei'].append(gene_id)
            if any([strain == x for x in ['62471', 'P295', 'R36_14']]):
                self.variant_dict['P.cactorum_apple'].append(gene_id)
            if any([strain == x for x in ['12420', '15_13', '15_7', '2003_3', '4032', '404', '415', '416', 'PC13_15', '414', '4040', '11-40', 'P421']]):
                self.variant_dict['P.cactorum_strawberry'].append(gene_id)
            if any([strain == x for x in ['17-21']]):
                self.variant_dict['P.cactorum_strawberry_LR'].append(gene_id)
    def count_variants(self, group):
        "Count the number of variants present for a species"
        SNP_count = len(self.variant_dict[group])
        return(SNP_count)
    def count_gene_variants(self, group):
        "Count the number of variants present for a species"
        SNP_count = len(set(self.variant_dict[group]))
        return(SNP_count)
        # for a in ['P.cactorum_strawberry', 'P.cactorum_apple', 'P.cactorum_strawberry_LR', 'P.idaei']:
        #     SNP_count = len(self.variant_dict[a])
        #     print("\t".join([a, str(SNP_count)]))
            # PcFa_count = len(self.variant_dict['P.cactorum_strawberry'])
            # PcMd_count = len(self.variant_dict['P.cactorum_apple'])
            # PcFaLR_count = len(self.variant_dict['P.cactorum_strawberry_LR'])
            # Pi_count = len(self.variant_dict['P.idaei'])
            # print ''


#-----------------------------------------------------
# Step X
# Additional methods
#-----------------------------------------------------

# def add_protseq(obj_dict, protein_lines):
#     """Add protein sequence data to each gene"""
#     for line in protein_lines:
#         line = line.rstrip()
#         if line.startswith('>'):
#             transcript_id = line.replace('>', '')
#         else:
#             obj_dict[transcript_id].protein_seq += line


#-----------------------------------------------------
# Step X
# Read input files
#-----------------------------------------------------

with open(conf.annotation_table) as f:
    a_lines = f.readlines()



#-----------------------------------------------------
# Step X
# Identify variants
#-----------------------------------------------------

GenesObj = FamilyVarObj()
GenesObj.add_family('Genes')
GenesObj.add_family('RxLR')
GenesObj.add_family('CRN')
# SNP_genes_obj = TaxaVarObj()

for line in a_lines[1:]:
    line = line.rstrip()
    split_line = line.split('\t')
    gene_id = split_line[0]
    RxLR_col = split_line[13]
    CRN_col = split_line[16]
    # print CRN_col
    SNP_col = split_line[29]
    if SNP_col:
        for SNP_variants in SNP_col.split('|'):
            split_col = SNP_variants.split(':')
            # print split_col
            effect = split_col[0]
            taxa = split_col[1]
            GenesObj.gene_family_dict['Genes'].add_variant(gene_id, taxa)
            if RxLR_col:
                GenesObj.gene_family_dict['RxLR'].add_variant(gene_id, taxa)
            if CRN_col:
                GenesObj.gene_family_dict['CRN'].add_variant(gene_id, taxa)


# print(" ".join(SNP_obj.variant_dict['P.idaei']))
# print("SNP gene variants:")


#-----------------------------------------------------
# Step X
# Summarise variants
#-----------------------------------------------------

header_line = ['Group', 'Total non-syn SNPs', 'Genes with non-syn SNPs', 'RxLRs with non-syn SNPs', 'CRNs with non-syn SNPs']
print("\t".join(header_line))
for a in ['P.cactorum_strawberry', 'P.cactorum_apple', 'P.cactorum_strawberry_LR', 'P.idaei']:
    outline = [a]
    outline.append(str(GenesObj.gene_family_dict['Genes'].count_variants(a)))
    outline.append(str(GenesObj.gene_family_dict['Genes'].count_gene_variants(a)))
    outline.append(str(GenesObj.gene_family_dict['RxLR'].count_gene_variants(a)))
    outline.append(str(GenesObj.gene_family_dict['CRN'].count_gene_variants(a)))
    print("\t".join(outline))
    # len(self.variant_dict[a])
    # print("\t".join([a, str(SNP_count)]))

    # SNP_variants =
    # print(SNP_col)
