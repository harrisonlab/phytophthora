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
        self.gene_family_dict['Genes'] = TaxaVarObj()
        self.gene_family_dict['RxLR'] = TaxaVarObj()
        self.gene_family_dict['CRN'] = TaxaVarObj()
        self.gene_family_dict['CAZY'] = TaxaVarObj()
        self.gene_family_dict['Apoplastic:Cutinase'] = TaxaVarObj()
        self.gene_family_dict['Apoplastic:Glucanase inhibitor'] = TaxaVarObj()
        self.gene_family_dict['Apoplastic:NLP'] = TaxaVarObj()
        self.gene_family_dict['Apoplastic:Phytotoxin'] = TaxaVarObj()
        self.gene_family_dict['Apoplastic:Protease inhibitor (cathepsin)'] = TaxaVarObj()
        self.gene_family_dict['Apoplastic:Protease inhibitor (cystatin-like)'] = TaxaVarObj()
        self.gene_family_dict['Apoplastic:Protease inhibitor (Kazal-type)'] = TaxaVarObj()
        self.gene_family_dict['MAMP:Elicitin'] = TaxaVarObj()
        self.gene_family_dict['MAMP:Transglutaminase'] = TaxaVarObj()
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
    def get_all_variants(self):
        "Count the number of variants present for all species"
        list1 = self.variant_dict['P.idaei']
        list2 = self.variant_dict['P.cactorum_apple']
        list3 = self.variant_dict['P.cactorum_strawberry']
        list4 = self.variant_dict['P.cactorum_strawberry_LR']
        total_list = list1 + list2 + list3 + list4
        return(total_list)
        # num_genes = len(set(total_list))
        # return(num_genes)



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

# TotalGenes = FamilyVarObj()
# AllVarObj = FamilyVarObj()
SnpObj = FamilyVarObj()
InDelObj = FamilyVarObj()
SvObj = FamilyVarObj()

# SNP_genes_obj = TaxaVarObj()

for line in a_lines[1:]:
    line = line.rstrip()
    split_line = line.split('\t')
    gene_id = split_line[0]
    # if gene_id == 'g26864.t1':
    #     print line
    secreted_col = split_line[10]
    RxLR_col = split_line[13]
    CRN_col = split_line[16]
    CAZY_col = split_line[18]
    apoplast_col = split_line[26]
    SNP_col = split_line[29]
    if SNP_col:
        for SNP_variants in SNP_col.split('|'):
            split_col = SNP_variants.split(':')
            # print split_col
            effect = split_col[0]
            SNP_taxa = split_col[1]
            SnpObj.gene_family_dict['Genes'].add_variant(gene_id, SNP_taxa)
            if RxLR_col:
                SnpObj.gene_family_dict['RxLR'].add_variant(gene_id, SNP_taxa)
            if CRN_col:
                SnpObj.gene_family_dict['CRN'].add_variant(gene_id, SNP_taxa)
            if all([secreted_col, CAZY_col]):
                SnpObj.gene_family_dict['CAZY'].add_variant(gene_id, SNP_taxa)
            if all([secreted_col, apoplast_col]):
                SnpObj.gene_family_dict[apoplast_col].add_variant(gene_id, SNP_taxa)
    InDel_col = split_line[30]
    if InDel_col:
        for InDel_variants in InDel_col.split('|'):
            split_col = InDel_variants.split(':')
            # print split_col
            effect = split_col[0]
            indel_taxa = split_col[1]
            InDelObj.gene_family_dict['Genes'].add_variant(gene_id, indel_taxa)
            if RxLR_col:
                InDelObj.gene_family_dict['RxLR'].add_variant(gene_id, indel_taxa)
            if CRN_col:
                InDelObj.gene_family_dict['CRN'].add_variant(gene_id, indel_taxa)
            if all([secreted_col, CAZY_col]):
                InDelObj.gene_family_dict['CAZY'].add_variant(gene_id, indel_taxa)
            if all([secreted_col, apoplast_col]):
                InDelObj.gene_family_dict[apoplast_col].add_variant(gene_id, indel_taxa)
    SV_col = split_line[31]
    if SV_col:
        for SV_variants in SV_col.split('|'):
            split_col = SV_variants.split(':')
            effect = split_col[0]
            sv_taxa = split_col[1]
            SvObj.gene_family_dict['Genes'].add_variant(gene_id, sv_taxa)
            if RxLR_col:
                SvObj.gene_family_dict['RxLR'].add_variant(gene_id, sv_taxa)
            if CRN_col:
                SvObj.gene_family_dict['CRN'].add_variant(gene_id, sv_taxa)
            if all([secreted_col, CAZY_col]):
                SvObj.gene_family_dict['CAZY'].add_variant(gene_id, sv_taxa)
            if all([secreted_col, apoplast_col]):
                SvObj.gene_family_dict[apoplast_col].add_variant(gene_id, sv_taxa)



# print(" ".join(SNP_obj.variant_dict['P.idaei']))
# print("SNP gene variants:")


#-----------------------------------------------------
# Step X
# Summarise variants
#-----------------------------------------------------

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
header_line = ['Group', 'Total non-syn SNPs']
header_line.extend(family_list)

for variant_type in [SnpObj, InDelObj, SvObj]:
    print("\t".join(header_line))
    for a in ['P.cactorum_strawberry', 'P.cactorum_apple', 'P.cactorum_strawberry_LR', 'P.idaei']:
        outline = [a]
        outline.append(str(variant_type.gene_family_dict['Genes'].count_variants(a)))
        for family in family_list:
            outline.append(str(variant_type.gene_family_dict[family].count_gene_variants(a)))
        print("\t".join(outline))


header_line = ['Group']
header_line.extend(family_list)
print("\t".join(header_line))
outline = ['total_vairants']
# outline.append(str(variant_type.gene_family_dict['Genes'].count_variants(a)))
for family in family_list:
    SNP_vars = SnpObj.gene_family_dict[family].get_all_variants()
    indel_vars = InDelObj.gene_family_dict[family].get_all_variants()
    sv_vars = SvObj.gene_family_dict[family].get_all_variants()
    total_vars = set(SNP_vars + indel_vars + sv_vars)
    outline.append(str(len(total_vars)))
print("\t".join(outline))

for a in ['P.cactorum_strawberry', 'P.cactorum_apple', 'P.cactorum_strawberry_LR', 'P.idaei']:
    outline = [a]
    for family in family_list:
        SNP_vars = SnpObj.gene_family_dict[family].variant_dict[a]
        indel_vars = InDelObj.gene_family_dict[family].variant_dict[a]
        sv_vars = SvObj.gene_family_dict[family].variant_dict[a]
        total_vars = set(SNP_vars + indel_vars + sv_vars)
        outline.append(str(len(total_vars)))
    print("\t".join(outline))



# SNP_vars = SnpObj.gene_family_dict['Apoplastic:Protease inhibitor (Kazal-type)'].get_all_variants()
# indel_vars = InDelObj.gene_family_dict['Apoplastic:Protease inhibitor (Kazal-type)'].get_all_variants()
# sv_vars = SvObj.gene_family_dict['Apoplastic:Protease inhibitor (Kazal-type)'].get_all_variants()
# total_vars = set(SNP_vars + indel_vars + sv_vars)
# print(" ".join(total_vars))
# print(str(len(total_vars)))
    # # Output InDel variants
    # outline.append(str(InDelObj.gene_family_dict['Genes'].count_variants(a)))
    # for family in family_list:
    #     outline.append(str(InDelObj.gene_family_dict[family].count_gene_variants(a)))
    # # Output structural variants
    # outline.append(str(SvObj.gene_family_dict['Genes'].count_variants(a)))
    # for family in family_list:
    #     outline.append(str(SvObj.gene_family_dict[family].count_gene_variants(a)))
    # print("\t".join(outline))
