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
ap.add_argument('--gene_gff',required=True,type=str,help='Gff file of predicyted gene models')
ap.add_argument('--gene_fasta',required=True,type=str,help='amino acid sequence of predicted proteins')
ap.add_argument('--SigP2',required=True,type=str,help='fasta file of genes testing positive for signal peptide using SigP2.0')
# ap.add_argument('--SigP4',required=True,type=str,help='fasta file of genes testing positive for signal peptide using SigP4.1')
ap.add_argument('--phobius',required=True,type=str,help='txt file ofheaders from gene testing positive for signal peptide using phobius')
ap.add_argument('--RxLR_motif',required=True,type=str,help='fasta file of genes testing positive for RxLR-EER motifs')
ap.add_argument('--RxLR_Hmm',required=True,type=str,help='fasta file of genes testing positive for RxLR-EER domains using an hmm model')
ap.add_argument('--RxLR_WY',required=True,type=str,help='fasta file of genes testing positive for WY domains using an hmm model')
ap.add_argument('--CRN_LFLAK',required=True,type=str,help='fasta file of genes testing positive for LFLAK domains using an hmm model')
ap.add_argument('--CRN_DWL',required=True,type=str,help='fasta file of genes testing positive for DWL domains using an hmm model')

conf = ap.parse_args()



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

with open(conf.CRN_LFLAK) as f:
    CRN_LFLAK_lines = f.readlines()

with open(conf.CRN_DWL) as f:
    CRN_DWL_lines = f.readlines()

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

CRN_LFLAK_set = Set()
for line in CRN_LFLAK_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        header = split_line[0].replace('>', '')
        CRN_LFLAK_set.add(header)

#-----------------------------------------------------
# Load CRN DWL hmm +ve proteins into a set
#-----------------------------------------------------

CRN_DWL_set = Set()
for line in CRN_DWL_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        header = split_line[0].replace('>', '')
        CRN_DWL_set.add(header)

#-----------------------------------------------------
# Step 3
# Itterate through genes in file, identifying if
# they ahve associated information
#-----------------------------------------------------

for line in gene_lines:
    line = line.rstrip()
    split_line = line.split()
    if 'transcript' in split_line[2]:
        # Set defaults
        sigP2 = ''
        sigP4 = ''
        phobius = ''
        RxLR_motif = ''
        RxLR_hmm = ''
        WY_hmm = ''
        CRN_LFLAK = ''
        CRN_DWL = ''
        orthogroup = ''
        prot_seq = ''
        # Identify gene id
        transcript_id = split_line[8]

        if transcript_id in SigP2_set:
            sigP2 = 'Yes'
        # if transcript_id in SigP4_set:
        #     sigP4 = 'Yes'
        if transcript_id in phobius_set:
            phobius = 'Yes'
        if transcript_id in RxLR_motif_set:
            RxLR_motif = 'Yes'
        if transcript_id in RxLR_hmm_set:
            RxLR_hmm = 'Yes'
        if transcript_id in RxLR_WY_set:
            WY_hmm = 'Yes'
        if transcript_id in CRN_LFLAK_set:
            CRN_LFLAK = 'Yes'
        if transcript_id in CRN_DWL_set:
            CRN_DWL = 'Yes'

        prot_seq = "".join(prot_dict[transcript_id])

        outline = [transcript_id, sigP2, sigP4, phobius ,RxLR_motif, RxLR_hmm, WY_hmm, CRN_LFLAK, CRN_DWL, orthogroup]
        outline.extend(split_line)
        outline.append(prot_seq)
        print "\t".join(outline)
