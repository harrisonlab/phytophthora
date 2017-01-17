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

ap = argparse.ArgumentParser()
ap.add_argument('--add_headers',required=True,type=str,help='List of ORFs to build the table for')
ap.add_argument('--gene_gff',required=True,type=str,help='Gff file of predicyted gene models')
ap.add_argument('--gene_fasta',required=True,type=str,help='amino acid sequence of predicted proteins')
ap.add_argument('--SigP2',required=True,type=str,help='fasta file of genes testing positive for signal peptide using SigP2.0')
# ap.add_argument('--SigP4',required=True,type=str,help='fasta file of genes testing positive for signal peptide using SigP4.1')
ap.add_argument('--phobius',required=True,type=str,help='txt file ofheaders from gene testing positive for signal peptide using phobius')
ap.add_argument('--RxLR_motif',required=True,type=str,help='fasta file of genes testing positive for RxLR-EER motifs')
ap.add_argument('--RxLR_Hmm',required=True,type=str,help='fasta file of genes testing positive for RxLR-EER domains using an hmm model')
ap.add_argument('--RxLR_WY',required=True,type=str,help='fasta file of genes testing positive for WY domains using an hmm model')
ap.add_argument('--CRN_LFLAK_DWL',required=True,type=str,help='fasta file of genes testing positive for LFLAK domains using an hmm model')


conf = ap.parse_args()

with open(conf.add_headers) as f:
    these_ORFs_lines = f.readlines()

with open(conf.gene_gff) as f:
    gene_lines = f.readlines()

with open(conf.gene_fasta) as f:
    prot_lines = f.readlines()

with open(conf.SigP2) as f:
    sigP2_lines = f.readlines()

# with open(conf.SigP4) as f:
#     sigP4_lines = f.readlines()

with open(conf.phobius) as f:
    phobius_lines = f.readlines()

with open(conf.RxLR_motif) as f:
    RxLR_motif_lines = f.readlines()

with open(conf.RxLR_Hmm) as f:
    RxLR_hmm_lines = f.readlines()

with open(conf.RxLR_WY) as f:
    RxLR_WY_lines = f.readlines()

with open(conf.CRN_LFLAK_DWL) as f:
    CRN_lines = f.readlines()

#-----------------------------------------------------
# Load Additional ORF headers into a set
#-----------------------------------------------------

these_ORFs_set = Set()
for line in these_ORFs_lines:
    header = line.rstrip()
    these_ORFs_set.add(header)

#-----------------------------------------------------
# Load protein sequence data into a dictionary
#-----------------------------------------------------

prot_dict = defaultdict(list)
for line in prot_lines:
    line = line.rstrip()
    if line.startswith('>'):
        header = line.replace('>', '')
    else:
        prot_dict[header] += line

#-----------------------------------------------------
# Load signalP2.0 files into a set
#-----------------------------------------------------

SigP2_set = Set()
for line in sigP2_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        header = split_line[0].replace('>', '')
        SigP2_set.add(header)

# #-----------------------------------------------------
# # Load signalP4.0 files into a set
# #-----------------------------------------------------
#
# SigP4_set = Set()
# for line in sigP4_lines:
#     line = line.rstrip()
#     if line.startswith('>'):
#         split_line = line.split()
#         header = split_line[0].replace('>', '')
#         SigP4_set.add(header)

#-----------------------------------------------------
# Load phobius files into a set
#-----------------------------------------------------

phobius_set = Set()
for line in phobius_lines:
    header = line.rstrip()
    phobius_set.add(header)

#-----------------------------------------------------
# Load RxLR motif +ve proteins into a set
#-----------------------------------------------------

RxLR_motif_set = Set()
for line in RxLR_motif_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        header = split_line[0].replace('>', '')
        RxLR_motif_set.add(header)

#-----------------------------------------------------
# Load RxLR hmm +ve proteins into a set
#-----------------------------------------------------

RxLR_hmm_set = Set()
for line in RxLR_hmm_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        header = split_line[0].replace('>', '')
        RxLR_hmm_set.add(header)

#-----------------------------------------------------
# Load RxLR hmm +ve proteins into a set
#-----------------------------------------------------

RxLR_WY_set = Set()
for line in RxLR_WY_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        header = split_line[0].replace('>', '')
        RxLR_WY_set.add(header)

#-----------------------------------------------------
# Load CRN LFLAK hmm +ve proteins into a set
#-----------------------------------------------------

CRN_set = Set()
for line in CRN_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        header = split_line[0].replace('>', '')
        CRN_set.add(header)

#-----------------------------------------------------
# Store genes and their ortholog groups in a dictionary
#-----------------------------------------------------

# organism_name = conf.ortho_name
# ortho_dict = defaultdict(list)
# for line in ortho_lines:
#     line = line.rstrip()
#     split_line = line.split()
#     orthogroup = split_line[0]
#     orthogroup = orthogroup.replace(":", "")
#     genes_in_group = [ x for x in split_line if organism_name in x ]
#     for gene in genes_in_group:
#         gene = gene.replace(organism_name, '').replace('|', '')
#         ortho_dict[gene] = orthogroup


#-----------------------------------------------------
# Step 3
# Itterate through genes in file, identifying if
# they ahve associated information
#-----------------------------------------------------
transcript_lines = []
# print these_ORFs_set

for line in gene_lines:
    line = line.rstrip()
    if line.startswith('#'):
        continue
    split_line = line.split()
    if 'transcript' in split_line[2] or 'mRNA' in split_line[2]:
        # Identify gene id
        split_col9 = split_line[8].split(';')
        transcript_id = "".join([ x for x in split_col9 if 'ID' in x ])
        transcript_id = transcript_id.replace('ID=', '')
        transcript_name = "".join([ x for x in split_col9 if 'Name' in x ])
        transcript_name = transcript_name.replace('Name=', '')
        # print transcript_id + "\t" + transcript_name
        if transcript_name not in these_ORFs_set:
            continue
        # Set defaults
        sigP2 = ''
        # sigP4 = ''
        phobius = ''
        RxLR_motif = ''
        RxLR_hmm = ''
        WY_hmm = ''
        CRN_LFLAK = ''
        CRN_DWL = ''
        orthogroup = 'NA'
        prot_seq = ''


        if transcript_name in SigP2_set:
            sigP2 = 'Yes'
        # if transcript_name in SigP4_set:
        #     sigP4 = 'Yes'
        if transcript_name in phobius_set:
            phobius = 'Yes'
        if transcript_name in RxLR_motif_set:
            RxLR_motif = 'Yes'
        if transcript_name in RxLR_hmm_set:
            RxLR_hmm = 'Yes'
        if transcript_name in RxLR_WY_set:
            WY_hmm = 'Yes'
        if transcript_name in CRN_set:
            CRN_LFLAK = 'Yes'
            CRN_DWL = 'Yes'
        # if ortho_dict[transcript_name]:
        #     orthogroup = ortho_dict[transcript_name]

        prot_seq = "".join(prot_dict[transcript_name])
        outline = [transcript_name, sigP2, phobius ,RxLR_motif, RxLR_hmm, WY_hmm, CRN_LFLAK, CRN_DWL, orthogroup]
        outline.extend(split_line)
        outline.append(prot_seq)
        print "\t".join(outline)
