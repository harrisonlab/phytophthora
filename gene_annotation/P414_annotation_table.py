#!/usr/bin/python

'''
This program is used to build information on all the genes predicted in
an annotated genome. These commands take information on location of genes
& suppliment this information with information on interproscan domains
and swissprot annotations.
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

# ap.add_argument('--genome',required=True,type=str,help='A fasta file of the assembled contigs')
ap.add_argument('--protein_fasta',required=True,type=str,help='A fasta file of predicted protein sequences')
ap.add_argument('--genes_gff',required=True,type=str,help='A gff file of the genes')
ap.add_argument('--conversion_log',required=True,type=str,help='A file detailing gene names prior to renaming')
ap.add_argument('--ORF_gff',required=True,type=str,help='A gff file for ORFs from which the two versions of ORF gene names can be extracted')
ap.add_argument('--SigP2',required=True,type=str,help='A file containing a list of signal-peptide containing genes')
ap.add_argument('--SigP3',required=True,type=str,help='A file containing a list of signal-peptide containing genes')
ap.add_argument('--SigP4',required=True,type=str,help='A file containing a list of signal-peptide containing genes')
ap.add_argument('--phobius',required=True,type=str,help='A fasta file containing signal-peptide containing genes')
ap.add_argument('--TM_list',required=True,type=str,help='txt file of headers from gene testing positive for tranmembrane proteins by TMHMM')
# # ap.add_argument('--GPI_list',required=True,type=str,help='txt file of headers from gene testing positive for GPI anchors as identified by GPI-SOM')
ap.add_argument('--RxLR_regex',required=True,type=str,help='A file containing results of RxLR regex searches')
ap.add_argument('--RxLR_hmm',required=True,type=str,help='A file containing results of RxLR hmm searches')
ap.add_argument('--CRN_dwl',required=True,type=str,help='A file containing results of DWL hmm searches')
ap.add_argument('--CRN_lflak',required=True,type=str,help='A file containing results of LFLAK hmm searches')
ap.add_argument('--EffP_list',required=True,type=str,help='A file containing results of effectorP')
ap.add_argument('--CAZY',required=True,type=str,help='A file containing results of CAZY')
ap.add_argument('--PhiHits',required=True,type=str,help='BLAST results of Phibase proteins vs predicted gene models')
# ap.add_argument('--ToxinHits',required=True,type=str,help='BLAST results of CDC toxins vs predicted gene models')
ap.add_argument('--InterPro',required=True,type=str,help='The Interproscan functional annotation .tsv file')
ap.add_argument('--Swissprot',required=True,type=str,help='A parsed table of BLAST results against the Swissprot database. Note - must have been parsed with swissprot_parser.py')
ap.add_argument('--TFs',required=True,type=str,help='Tab seperated of putative transcription factors and their domains as identified by interpro2TFs.py')
ap.add_argument('--DEGs',required=True,type=str,nargs='+',help='A list of files containing lists of DEGs')
ap.add_argument('--fpkm',required=True,type=str,help='File containing FPKM values of genes under all conditions')
ap.add_argument('--SNPs',required=True,type=str,help='Annotated vcf file, detailing SNPs and effects on genes')
ap.add_argument('--InDels',required=True,type=str,help='Annotated vcf file, detailing InDels and effects on genes')
ap.add_argument('--SVs',required=True,type=str,help='Annotated vcf file, detailing SVs and effects on genes')
ap.add_argument('--phasing',required=True,type=str,help='A text file giving phasing information on SNP data')
ap.add_argument('--orthogroups', required=True,type=str,help='A file containing results of orthology analysis')
ap.add_argument('--strain_id',required=True,type=str,help='The identifier of this strain as used in the orthology analysis')
ap.add_argument('--OrthoMCL_all',required=True,type=str,nargs='+',help='The identifiers of all strains used in the orthology analysis')

conf = ap.parse_args()



#-----------------------------------------------------
# Step 2
# Define classes
#-----------------------------------------------------

def add_protseq(obj_dict, protein_lines):
    """Add protein sequence data to each gene"""
    for line in protein_lines:
        line = line.rstrip()
        if line.startswith('>'):
            transcript_id = line.replace('>', '')
        else:
            obj_dict[transcript_id].protein_seq += line
def summarise_orthogroup(content_str):
    """"""
    content_list = content_str.split(",")
    content_list = [x.split('|')[0] for x in content_list]
    content_set = set(content_list)
    crown_rot_isolates = ['Pc_CR1', 'Pc_CR2', 'Pc_CR3', 'Pc_CR4', 'Pc_CR5', 'Pc_CR6', 'Pc_CR7', 'Pc_CR8', 'Pc_CR9', 'Pc_CR10', 'Pc_CR11', 'Pc_CR12', 'Pc_CR13', 'Pc_LR1']
    leather_rot_isolates = ['Pc_LR2']
    apple_isolates = ['Pc_MD1', 'Pc_MD2', 'Pc_MD3']
    raspberry_isolates = ['Pi_RI1', 'Pi_RI2', 'Pi_RI3']
    cr = 0
    lr = 0
    md = 0
    ri = 0
    # print content_set
    for i in content_set:
        if i in crown_rot_isolates:
            cr += 1
        elif i in leather_rot_isolates:
            lr += 1
        elif i in apple_isolates:
            md  += 1
        elif i in raspberry_isolates:
            ri  += 1
    summarised_orthogroup = "".join(["CR(", str(cr), ")LR(", str(lr), ")Md(", str(md), ")Ri(", str(ri), ")"])
    return(summarised_orthogroup)

def get_orthogroup_expansion_status(orthogroup_content_dict):
    """"""
    # content_list = content_str.split(",")
    # content_list = [x.split('|')[0] for x in content_list]
    # content_set = set(content_list)
    # crown_rot_isolates = ['Pc_CR1', 'Pc_CR2', 'Pc_CR3', 'Pc_CR4', 'Pc_CR5', 'Pc_CR6', 'Pc_CR7', 'Pc_CR8', 'Pc_CR9', 'Pc_CR10', 'Pc_CR11', 'Pc_CR12', 'Pc_CR13', 'Pc_LR1']
    crown_rot_isolates = ['Pc_CR1', 'Pc_CR2', 'Pc_CR3', 'Pc_CR4', 'Pc_CR5', 'Pc_CR6', 'Pc_CR7', 'Pc_CR8', 'Pc_CR9', 'Pc_CR10', 'Pc_CR11', 'Pc_CR12', 'Pc_CR13']
    # leather_rot_isolates = ['Pc_LR2']
    apple_isolates = ['Pc_MD1', 'Pc_MD2', 'Pc_MD3']
    raspberry_isolates = ['Pi_RI1', 'Pi_RI2', 'Pi_RI3']
    cr = []
    # lr = []
    md = []
    ri = []
    for x in crown_rot_isolates:
        counts = int(orthogroup_content_dict[x])
        cr.append(counts)
    for x in apple_isolates:
        counts = int(orthogroup_content_dict[x])
        md.append(counts)
    for x in raspberry_isolates:
        counts = int(orthogroup_content_dict[x])
        ri.append(counts)
    if len(cr) == 0:
        cr = [0]
    if len(md) == 0:
        md = [0]
    if len(ri) == 0:
        ri = [0]
    expansion_status = ""
    if min(cr) > max(md) and min(cr) > max(ri):
        expansion_status = "crown rot expanded"
    elif max(cr) < min(md) and max(cr) < min(ri):
        expansion_status = "crown rot loss"
    elif min(md) > max(cr) and min(md) > max(ri):
        expansion_status = "apple expanded"
    elif max(md) < min(cr) and max(md) < min(ri):
        expansion_status = "apple loss"
    # print orthogroup_content_dict
    # print cr
    # print md
    # print ri
    # print expansion_status
    # exit()
    return(expansion_status)


def get_orthogroup_expansion_status2(orthogroup_content_dict):
    """"""
    clade_A_isolates = ['Pi_RI1', 'Pi_RI2', 'Pi_RI3']
    # clade_B_isolates = []
    clade_C_isolates = ['Pc_CR1', 'Pc_CR2', 'Pc_CR3', 'Pc_CR4', 'Pc_CR5', 'Pc_CR6', 'Pc_CR7', 'Pc_CR8', 'Pc_CR9', 'Pc_CR10', 'Pc_CR11', 'Pc_CR12', 'Pc_CR13','Pc_LR2', 'Pc_MD1', 'Pc_MD2', 'Pc_MD3']
    clade_D_isolates = ['Pc_LR2', 'Pc_MD1', 'Pc_MD2', 'Pc_MD3']
    clade_E_isolates = ['Pc_MD1', 'Pc_MD2', 'Pc_MD3']
    clade_F_isolates = ['Pc_CR1', 'Pc_CR2', 'Pc_CR3', 'Pc_CR4', 'Pc_CR5', 'Pc_CR6', 'Pc_CR7', 'Pc_CR8', 'Pc_CR9', 'Pc_CR10', 'Pc_CR11', 'Pc_CR12', 'Pc_CR13']
    clade_G_isolates = ['Pc_CR1', 'Pc_CR2', 'Pc_CR3', 'Pc_CR4', 'Pc_CR5', 'Pc_CR6', 'Pc_CR7', 'Pc_CR8', 'Pc_CR9', 'Pc_CR10', 'Pc_CR11', 'Pc_CR12']
    clade_H_isolates = ['Pc_CR13']
    clade_I_isolates = ['Pc_LR2']

    c_A = []
    # c_B = []
    c_C = []
    c_D = []
    c_E = []
    c_F = []
    c_G = []
    c_H = []
    c_I = []
    for x in clade_A_isolates:
        counts = int(orthogroup_content_dict[x])
        c_A.append(counts)
    # for x in clade_B_isolates:
    #     counts = int(orthogroup_content_dict[x])
    #     c_B.append(counts)
    for x in clade_C_isolates:
        counts = int(orthogroup_content_dict[x])
        c_C.append(counts)
    for x in clade_D_isolates:
        counts = int(orthogroup_content_dict[x])
        c_D.append(counts)
    for x in clade_E_isolates:
        counts = int(orthogroup_content_dict[x])
        c_E.append(counts)
    for x in clade_F_isolates:
        counts = int(orthogroup_content_dict[x])
        c_F.append(counts)
    for x in clade_G_isolates:
        counts = int(orthogroup_content_dict[x])
        c_G.append(counts)
    for x in clade_H_isolates:
        counts = int(orthogroup_content_dict[x])
        c_H.append(counts)
    for x in clade_I_isolates:
        counts = int(orthogroup_content_dict[x])
        c_I.append(counts)

    if len(c_A) == 0:
        C_A = [0]
    # if len(c_B) == 0:
    #     c_B = [0]
    if len(c_C) == 0:
        c_C = [0]
    if len(c_D) == 0:
        c_D = [0]
    if len(c_E) == 0:
        c_E = [0]
    if len(c_F) == 0:
        c_F = [0]
    if len(c_G) == 0:
        c_G = [0]
    if len(c_H) == 0:
        c_H = [0]
    if len(c_I) == 0:
        c_I = [0]

    expansion_status = []
    # if min(c_A) > max(c_C):
    #     expansion_status.append('A_gain')
    # if min(c_C) > max(c_A):
    #     expansion_status.append('A_loss')
    if min(c_G) > max(c_A + c_D + c_H):
        expansion_status.append('G_gain')
    if max(c_G) < min(c_A + c_D + c_H):
        expansion_status.append('G_loss')
    if min(c_H) > max(c_A + c_D + c_G):
        expansion_status.append('H_gain')
    if max(c_H) < min(c_A + c_D + c_G):
        expansion_status.append('H_loss')
    if min(c_E) > max(c_A + c_F + c_I):
        expansion_status.append('E_gain')
    if max(c_E) < min(c_A + c_F + c_I):
        expansion_status.append('E_loss')
    if min(c_I) > max(c_A + c_F + c_E):
        expansion_status.append('I_gain')
    if max(c_I) < min(c_A + c_F + c_E):
        expansion_status.append('I_loss')
    if min(c_F) > max(c_A + c_D):
        expansion_status.append('F_gain')
    if max(c_F) < min(c_A + c_D):
        expansion_status.append('F_loss')
    if min(c_D) > max(c_A + c_F):
        expansion_status.append('D_gain')
    if max(c_D) < min(c_A + c_F):
        expansion_status.append('D_loss')
    if min(c_C) > max(c_A):
        expansion_status.append('A_gain or C_loss')
    if max(c_C) < min(c_A):
        expansion_status.append('A_loss or C_gain')
    return(";".join(expansion_status))

class ConvObj(object):
    """A dictionary of old and current gene names, with methods to aid conversion
    Attributes:
        conversion_dict: A string representing the gene name.
    """
    def __init__(self):
        """set the object, defining the dictionary."""
        dict = defaultdict(list)
        self.conversion_dict = dict
    def set_dict(self, lines):
        """ Create the conversion dictionary from input lines """
        for line in lines:
            line = line.rstrip()
            split_line = line.split("\t")
            old_gene_id = split_line[0]
            new_gene_id = split_line[2]
            conv_dict = self.conversion_dict
            conv_dict[old_gene_id] = new_gene_id
            self.conversion_dict = conv_dict
    def set_ORF_dict(self, lines):
        """ Create a conversion dictionary from ORF gff input lines """
        for line in lines:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            split_line = line.split("\t")
            if split_line[2] == 'gene':
                continue
            col9 = split_line[8]
            split_col9 = col9.split(";")
            old_gene_id = split_col9[2].replace("Name=", "")
            new_gene_id = split_col9[0].replace("ID=", "")
            # conv_dict = self.conversion_dict
            # conv_dict[old_gene_id] = new_gene_id
            # self.conversion_dict = conv_dict
            # print str(old_gene_id) + "\t" + str(new_gene_id)
            self.conversion_dict[old_gene_id] = new_gene_id
        # print self.conversion_dict
    def convert_transcript(self, old_transcript_id):
        """  """
        old_gene_id = re.sub(r".t.*", "" , old_transcript_id)
        if self.conversion_dict[old_gene_id]:
            gene_id = self.conversion_dict[old_gene_id]
            transcript_id = old_transcript_id.replace(old_gene_id, gene_id)
        else:
            transcript_id = "not present"
        # print str(old_gene_id) + "\t" + str(gene_id)
        # print "\t".join([old_transcript_id, old_gene_id, gene_id, transcript_id])
        return transcript_id
    def convert_ORF_transcript(self, old_transcript_id):
        """  """
        transcript_id = self.conversion_dict[old_transcript_id]
        return transcript_id



class Expression_obj(object):
    """Information on a genes fpkm under various conditions and whether a gene
        is differentially expressed under different conditions:
    """
    def __init__(self, gene_id):
        """Return a Expression_obj whose name is *gene_id*"""
        self.gene_id = gene_id
        self.DEG_conditions = set()
        # self.fpkm_dict = defaultdict(list)
        self.condition_list = []
        self.fpkm_list = []

    def add_DEG(self, DEG_condition):
        """Return a Expression_obj whose name is *gene_id*"""
        self.DEG_conditions.add(DEG_condition)

    def add_fpkm(self, conditions_list, fpkm_list):
        """Return a Expression_obj whose name is *gene_id*"""
        for condition, fpkm in zip(conditions_list, fpkm_list):
            self.condition_list.append(condition)
            fpkm = int(np.round_(float(fpkm),  decimals=0))
            self.fpkm_list.append(str(fpkm))


class Annot_obj(object):
    # """A gene identified as differentially expressed in one isolate.
    # Attributes:
    #     transcript_id: A string representing the gene name.
    #     conditions_tested: List of conditions that the gene was tested for DE.
    #     conditions_positive: A binary list of integers representing whether the
    #         gene tested positive for each condition listed in conditions tested.
    # """

    def __init__(self):
        """Return a Annot_obj whose name is *transcript_id*"""
        self.transcript_id = ''
        self.old_transcript = ''
        self.contig = ''
        self.start = ''
        self.stop = ''
        self.strand = ''
        self.protein_seq = ''
        self.sigp2 = ''
        self.sigp3 = ''
        self.sigp4 = ''
        self.phobius = ''
        self.transmem = ''
        self.secreted = ''
        self.rxlr_regex = ''
        self.rxlr_deer = ''
        self.rxlr_hmm = ''
        self.putative_rxlr = ''
        self.dwl_hmm = ''
        self.lflak_hmm = ''
        self.putative_crn = ''
        self.effp = ''
        self.cazy = ''
        self.interpro = set()
        self.swissprot = ''
        self.phi = ''
        self.toxin = ''
        self.orthogroup_id = ''
        self.content_counts = ''
        self.content_str = ''
        self.transcriptionfactor = set()
        self.DEG_conditions = ''
        self.fpkm_conditions = ''
        self.fpkm = ''
        self.deg = ''
        self.exp_time = ''
        self.SNP_info = []
        self.SNP_phase = defaultdict(int)
        self.InDel_info = []
        self.SV_info = []
        self.variation_level = set()
        self.ipr_effectors = set()

    def set_conditions(self, gff_elements):
        """Reset conditions_tested to a list of conditions"""
        self.contig = gff_elements[0]
        self.start = gff_elements[3]
        self.stop = gff_elements[4]
        self.strand = gff_elements[6]
        gene_features = gff_elements[8].split(';')
        transcript_id = gene_features[0]
        transcript_id = transcript_id.replace('ID=', '')
        self.transcript_id = transcript_id

    def set_prev_id(self, old_transcript):
        """Reset conditions_tested to a list of conditions"""
        self.old_transcript = old_transcript
    def add_sigp2(self):
        """Add SignalP information"""
        self.sigp2 = 'Yes'
        self.secreted = 'Yes'
    def add_sigp3(self):
        """Add SignalP information"""
        self.sigp3 = 'Yes'
        self.secreted = 'Yes'
    def add_sigp4(self):
        """Add SignalP information"""
        self.sigp4 = 'Yes'
        self.secreted = 'Yes'
    def add_phobius(self):
        """Add signal peptide information"""
        self.phobius = 'Yes'
        # self.secreted = 'Yes'
    def add_transmem(self):
        """Add transmembrane protein information"""
        # print self.transcript_id
        if self.secreted == 'Yes':
            self.transmem = 'Yes'
            self.secreted = ''
    def add_rxlr_regex(self, line):
        """Add RxLR information"""
        self.rxlr_regex = 'Yes'
        if 'EER_motif_start' in line:
            self.rxlr_deer = 'Yes'
            if any('Yes' in x for x in ([self.sigp2, self.sigp3, self.sigp4])):
                self.putative_rxlr = 'RxLR'
    def add_rxlr_hmm(self):
        """Add RxLR information"""
        self.rxlr_hmm = 'Yes'
        if any('Yes' in x for x in ([self.sigp2, self.sigp3, self.sigp4])):
            self.putative_rxlr = 'RxLR'
    def add_dwl_hmm(self):
        """Add DWL information"""
        self.dwl_hmm = 'Yes'
    def add_lflak_hmm(self):
        """Add LFLAK information"""
        self.lflak_hmm = 'Yes'
        if self.dwl_hmm == 'Yes':
            self.putative_crn = 'CRN'
    def add_effp(self):
        """Add EffectorP information"""
        if any('Yes' in x for x in ([self.sigp2, self.sigp3, self.sigp4])):
            self.effp = 'Yes'
    def add_cazy(self, line):
        """Add CAZY information"""
        split_line = line.split()
        self.cazy = 'CAZY:' + split_line[0].replace('.hmm', '')
    def add_interpro(self, line):
        """Add InterPro information"""
        split_line = line.split("\t")
        interpro_columns = []
        index_list = [4, 5, 11, 12]
        for x in index_list:
            if len(split_line) > x:
                interpro_columns.append(split_line[x])
        ipr_line = ";".join(interpro_columns)
        self.interpro.add(ipr_line)
    def add_swissprot(self, line):
        """Add swissprot information"""
        split_line = line.split("\t")
        self.swissprot = "|".join(itemgetter(14, 12, 13)(split_line))
    def add_phi(self, line):
        line = line.rstrip()
        split_line = line.split()
        self.phi = split_line[1]
    def add_toxin(self, line):
        line = line.rstrip()
        split_line = line.split()
        self.toxin = split_line[1]
    def add_orthogroup(self, orthogroup_id, content_counts, content_str):
        """Add swissprot information"""
        self.orthogroup_id = orthogroup_id
        self.content_counts = content_counts
        self.content_str = content_str
    def add_TF(self, line):
        """Add swissprot information"""
        split_line = line.split("\t")
        TF_function = split_line[2]
        self.transcriptionfactor.add(TF_function)

    def add_SNP(self, line, isolate_list):
        """Add information on non-synonymous SNPs"""
        species_isolates = ['371','SCRP370','SCRP376']
        pathotype_isolates = ['62471','P295','R36_14']
        leather_rot = ['17-21', '11-40']
        population_isolates = ['415','15_7','15_13', '404', '12420', '12420', '2003_3', '4032', 'PC13_15', '414', '4040', '416', '415']
        split_line = line.split("\t")
        effect_col = split_line[7]
        SNP_effect = effect_col.split('|')[1]
        presence_list = []

        for isolate, isolate_col in zip(isolate_list, split_line[9:]):
            if any(x in isolate_col for x in (['1/1', '0/1', '1/0'])):
                presence_list.append(isolate)
        self.SNP_info.append(":".join([SNP_effect, ",".join(presence_list)]))
        if any(x in presence_list for x in (species_isolates)):
            self.variation_level.add('species')
        if any(x in presence_list for x in (pathotype_isolates)):
            self.variation_level.add('pathotype')
        if any(x in presence_list for x in (population_isolates)):
            self.variation_level.add('population')


    def add_InDel(self, line, isolate_list):
        """Add information on non-synonymous InDels"""
        species_isolates = ['371','SCRP370','SCRP376']
        pathotype_isolates = ['62471','P295','R36_14']
        population_isolates = ['415','15_7','15_13', '404', '12420', '12420', '2003_3', '4032', 'PC13_15', '414', '4040', '416']
        split_line = line.split("\t")
        effect_col = split_line[7]
        InDel_effect = effect_col.split('|')[1]
        presence_list = []
        for isolate, isolate_col in zip(isolate_list, split_line[9:]):
            if any(x in isolate_col for x in (['1/1', '0/1', '1/0'])):
                isolate = isolate.replace('_vs_414_aligned_sorted.bam', '')
                presence_list.append(isolate)
        if len(presence_list) > 0:
            self.InDel_info.append(":".join([InDel_effect, ",".join(presence_list)]))
        if any(x in presence_list for x in (species_isolates)):
            self.variation_level.add('species')
        if any(x in presence_list for x in (pathotype_isolates)):
            self.variation_level.add('pathotype')
        if any(x in presence_list for x in (population_isolates)):
            self.variation_level.add('population')

    def add_SV(self, line, isolate_list):
        """Add information on non-synonymous SVs"""
        species_isolates = ['371','SCRP370','SCRP376']
        pathotype_isolates = ['62471','P295','R36_14']
        population_isolates = ['415','15_7','15_13', '404', '12420', '12420', '2003_3', '4032', 'PC13_15', '414', '4040', '416']
        split_line = line.split("\t")
        effect_col = split_line[7]
        SV_effect = effect_col.split('|')[1]
        presence_list = []
        for isolate, isolate_col in zip(isolate_list, split_line[9:]):
            if any(x in isolate_col for x in (['1/1', '0/1', '1/0'])):
                isolate = isolate.replace('_vs_414_aligned_sorted.bam', '')
                presence_list.append(isolate)
        # print self.transcript_id + "\t" + ",".join(presence_list)
        self.SV_info.append(":".join([SV_effect, ",".join(presence_list)]))
        if any(x in presence_list for x in (species_isolates)):
            self.variation_level.add('species')
        if any(x in presence_list for x in (pathotype_isolates)):
            self.variation_level.add('pathotype')
        if any(x in presence_list for x in (population_isolates)):
            self.variation_level.add('population')

    def summarise_phase(self):
        """"""
        outline = []
        for state in self.SNP_phase.keys():
            count = self.SNP_phase[state]
            outline.append(state + "(" + str(count) + ")")
        return(";".join(outline))



    def add_expr(self, exp_obj):
        """"""
        self.DEG_conditions = ";".join(exp_obj.DEG_conditions)
        self.fpkm_conditions = "\t".join(exp_obj.condition_list)
        self.fpkm = "\t".join(exp_obj.fpkm_list)
        self.treatments = []
        self.mean_fpkm = []
        # self.FC_conditions =['Mycelium vs Emily 12 hpi',
        #     'Mycelium vs Emily 48 hpi',
        #     'Emily 12 hpi vs Emily 48 hpi',
        #     'Mycelium vs Fenella 12 hpi',
        #     'Mycelium vs Fenella 48 hpi',
        #     'Fenella 12 hpi vs Fenella 48 hpi']
        self.FC_values = []
        self.FC_values_adj = []


    def add_LFC(self):
        "Add LFC information based upon fpkm data already in the dictionary"
        fpkm_dict = defaultdict(list)
        for condition, fpkm in zip(self.fpkm_conditions.split("\t"), self.fpkm.split("\t")):
            fpkm_dict[condition].append(int(fpkm))
        fpkm_list = []
        condition_list = []
        for condition in fpkm_dict.keys():
            # condition_list.append(condition)
            fpkm_mean = np.mean(fpkm_dict[condition])
            fpkm_dict[condition].append(fpkm_mean)
            self.treatments.append(condition)
            self.mean_fpkm.append(str(np.round_(fpkm_mean,  decimals=1)))
            # int(np.round_(float(fpkm),  decimals=0)))
        # print "\t".join(condition_list)
        # print "\t".join(fpkm_list)
        LFC_conditions = [
        ('Mycelium - Mycelium', 'Emily - 12 hours'),
        ('Mycelium - Mycelium', 'Emily - 48 hours'),
        ('Emily - 12 hours', 'Emily - 48 hours'),
        ('Mycelium - Mycelium', 'Fenella - 12 hours'),
        ('Mycelium - Mycelium', 'Fenella - 48 hours'),
        ('Fenella - 12 hours', 'Fenella - 48 hours'),
        ]
        for x, y in LFC_conditions:
            x_fpkm_mean = fpkm_dict[x][-1]
            y_fpkm_mean = fpkm_dict[y][-1]
            if y_fpkm_mean == x_fpkm_mean:
                FC = 1
                # FC_adj = 1
                FC_adj = 0
            # elif any(num == 0.0 for num in [x, y]):
            #     LFC = 0
            elif y_fpkm_mean >= x_fpkm_mean:
                FC = np.divide(y_fpkm_mean, x_fpkm_mean)
                # if x_fpkm_mean < 2:
                #     x_adj = 2
                # else:
                #     x_adj = x_fpkm_mean
                # if y_fpkm_mean < 2:
                #     y_adj = 2
                # else:
                #     y_adj = y_fpkm_mean
                # # FC_adj = np.divide(y_adj, x_adj)
                # FC_adj = np.log2(y_adj) - np.log2(x_adj)
            else:
                FC = np.negative(np.divide(x_fpkm_mean, y_fpkm_mean))
                # if x_fpkm_mean < 2:
                #     x_adj = 2
                # else:
                #     x_adj = x_fpkm_mean
                # if y_fpkm_mean < 2:
                #     y_adj = 2
                # else:
                #     y_adj = y_fpkm_mean
                # # FC_adj = np.negative(np.divide(x_adj, y_adj))
                # FC_adj = np.log2(y_adj) - np.log2(x_adj)
            if x_fpkm_mean < 1:
                x_adj = 1
            else:
                x_adj = x_fpkm_mean
            if y_fpkm_mean < 1:
                y_adj = 1
            else:
                y_adj = y_fpkm_mean
            FC_adj = np.log2(y_adj) - np.log2(x_adj)
            self.FC_values.append(str(np.round_(FC,  decimals=1)))
            self.FC_values_adj.append(str(np.round_(FC_adj,  decimals=1)))
            # print "\t".join([x, str(np.round_(x_fpkm_mean,  decimals=1)), y, str(np.round_(y_fpkm_mean,  decimals=1)), 'FC', str(np.round_(FC,  decimals=1)), 'FC_adj', str(np.round_(FC_adj,  decimals=1))])

        # print self.DEG_conditions
        # print self.fpkm_conditions
        # print self.fpkm

            # self.fpkm_dict[condition].append(fpkm)
    def add_early_late(self):
        FC_conditions = ['Mycelium vs Emily 12 hpi',
            'Mycelium vs Emily 48 hpi',
            'Emily 12 hpi vs Emily 48 hpi',
            'Mycelium vs Fenella 12 hpi',
            'Mycelium vs Fenella 48 hpi',
            'Fenella 12 hpi vs Fenella 48 hpi']
        conditions = self.DEG_conditions.split("\t")

        FC_values_adj = self.FC_values_adj

        if all(abs(float(x)) < 2 for x in FC_values_adj):
            exp_time = "constant"
        elif all(float(x) > 0 for x in [FC_values_adj[0], FC_values_adj[3]]) and all(float(x) < -2 for x in [FC_values_adj[2], FC_values_adj[5]]):
            exp_time = "early expressed"
        elif all(float(x) > 0 for x in [FC_values_adj[1], FC_values_adj[4]]) and all(float(x) > 2 for x in [FC_values_adj[2], FC_values_adj[5]]):
            exp_time = "late expressed"
        elif all(float(x) >=2 for x in FC_values_adj[0:2] + FC_values_adj[3:5]):
            exp_time = "upregulated in planta (all)"
        elif all(float(x) <= -2 for x in FC_values_adj[0:2] + FC_values_adj[3:5]):
            exp_time = "downregulated in planta (all)"
        elif any(float(x) >=2 for x in FC_values_adj[0:2] + FC_values_adj[3:5]) and all(float(x) >=0 for x in FC_values_adj[0:2] + FC_values_adj[3:5]):
            exp_time = "upregulated in planta (some)"
        elif any(float(x) <= -2 for x in FC_values_adj[0:2] + FC_values_adj[3:5]) and all(float(x) <=0 for x in FC_values_adj[0:2] + FC_values_adj[3:5]):
            exp_time = "downregulated in planta (some)"
        else:
            exp_time = ""
        # print "\t".join([self.transcript_id, exp_time])
        # print self.treatments
        # print self.mean_fpkm
        # print FC_values_adj
        self.exp_time = exp_time
        if any(abs(float(x)) >= 2 for x in FC_values_adj) and self.DEG_conditions != '':
            self.deg = 'DEG'


    def ipr2effectors(self):
        """"""
        ipr_dict = {
        'IPR002200' : 'MAMP:Elicitin',
        'IPR032048' : 'MAMP:Transglutaminase',
        'IPR001314' : 'Apoplastic:Glucanase inhibitor',
        'PF09461' : 'Apoplastic:Phytotoxin',
        'PF05630' : 'Apoplastic:NLP',
        'IPR008701' : 'Apoplastic:NLP',
        'PF01083' : 'Apoplastic:Cutinase',
        'IPR002350' : 'Apoplastic:Protease inhibitor (Kazal-type)',
        'IPR013201' : 'Apoplastic:Protease inhibitor (cathepsin)',
        'IPR000010' : 'Apoplastic:Protease inhibitor (cystatin-like)',
        'IPR027214' : 'Apoplastic:Protease inhibitor (cystatin-like)',
        'IPR018073' : 'Apoplastic:Protease inhibitor (cystatin-like)',
        'IPR020381' : 'Apoplastic:Protease inhibitor (cystatin-like)'
        }
        for ipr_line in self.interpro:
            for x in ipr_line.split(';'):
                if x in ipr_dict:
                    self.ipr_effectors.add(ipr_dict[x])


#-----------------------------------------------------
# Step X
# Read input files
#-----------------------------------------------------

# with open(conf.genome) as f:
#     contig_lines = f.readlines()
with open(conf.genes_gff) as f:
    gff_lines = f.readlines()
with open(conf.protein_fasta) as f:
    protein_lines = f.readlines()
with open(conf.conversion_log) as f:
    conversion_lines = f.readlines()
with open(conf.ORF_gff) as f:
    ORF_gff_lines = f.readlines()
with open(conf.SigP2) as f:
    sigP2_lines = f.readlines()
with open(conf.SigP3) as f:
    sigP3_lines = f.readlines()
with open(conf.SigP4) as f:
    sigP4_lines = f.readlines()
with open(conf.phobius) as f:
    phobius_lines = f.readlines()
with open(conf.TM_list) as f:
    TM_lines = f.readlines()
with open(conf.RxLR_regex) as f:
    rxlr_regex_lines = f.readlines()
with open(conf.RxLR_hmm) as f:
    rxlr_hmm_lines = f.readlines()
with open(conf.CRN_dwl) as f:
    crn_dwl_lines = f.readlines()
with open(conf.CRN_lflak) as f:
    crn_lflak_lines = f.readlines()
with open(conf.EffP_list) as f:
    effP_lines = f.readlines()
with open(conf.TFs) as f:
    TF_lines = f.readlines()
with open(conf.CAZY) as f:
    cazy_lines = f.readlines()
with open(conf.PhiHits) as f:
    phibase_lines = f.readlines()
# with open(conf.ToxinHits) as f:
#     toxin_lines = f.readlines()
with open(conf.InterPro) as f:
    InterPro_lines = f.readlines()
with open(conf.Swissprot) as f:
    swissprot_lines = f.readlines()
with open(conf.orthogroups) as f:
    orthogroup_lines = f.readlines()
with open(conf.SNPs) as f:
    SNP_lines = f.readlines()
with open(conf.InDels) as f:
    InDel_lines = f.readlines()
with open(conf.SVs) as f:
    SV_lines = f.readlines()
with open(conf.phasing) as f:
    phase_lines = f.readlines()

#-----------------------------------------------------
# Step X
# Create an object to hold expression information
#-----------------------------------------------------


expression_dict = defaultdict(list)
for file in conf.DEGs:
    with open(file) as f:
        DEG_lines = f.readlines()
    condition = os.path.basename(file)
    condition = condition.replace('_DEGs.txt', '')
    for line in DEG_lines:
        gene_id = line.rstrip()
        if expression_dict[gene_id]:
            exp_obj = expression_dict[gene_id][0]
        else:
            exp_obj = Expression_obj(gene_id)
        # exp_obj.gene_id = gene_id
        exp_obj.add_DEG(condition)
        expression_dict[gene_id].append(exp_obj)

with open(conf.fpkm) as f:
    fpkm_lines = f.readlines()

conditions_list = fpkm_lines[0].rstrip().split("\t")
for line in fpkm_lines[1:]:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    fpkm_list = split_line[1:]
    if expression_dict[gene_id]:
        expression_dict[gene_id][0].add_fpkm(conditions_list, fpkm_list)
    else:
        exp_obj = Expression_obj(gene_id)
        exp_obj.add_fpkm(conditions_list, fpkm_list)
        expression_dict[gene_id].append(exp_obj)


gene_dict = defaultdict(list)


#-----------------------------------------------------
# Step X
# Create gene conversion object
#-----------------------------------------------------

conversion_obj = ConvObj()
conversion_obj.set_dict(conversion_lines)

ORF_conversion_obj = ConvObj()
ORF_conversion_obj.set_ORF_dict(ORF_gff_lines)

#-----------------------------------------------------
# Step X
# Create an annotation object for each gene
#-----------------------------------------------------


for line in gff_lines:
    if "gff-version" in line:
        continue
    if line.startswith('#'):
        continue
    line = line.rstrip()
    split_line = line.split("\t")
    if 'mRNA' in split_line[2]:
        gene_features = split_line[8].split(';')
        transcript_id = gene_features[0]
        transcript_id = transcript_id.replace('ID=', '')
        gene_obj = Annot_obj()
        gene_obj.set_conditions(split_line)
        gene_dict[transcript_id] = gene_obj


#-----------------------------------------------------
# Step X
# Add annotation information to gene objects
#-----------------------------------------------------

add_protseq(gene_dict, protein_lines)

for line in sigP2_lines:
    line = line.rstrip()
    split_line = line.split()
    old_transcript_id = split_line[0]
    if '.t' not in old_transcript_id:
        old_transcript_id = ORF_conversion_obj.convert_ORF_transcript(old_transcript_id)
    transcript_id = conversion_obj.convert_transcript(old_transcript_id)
    if transcript_id == 'not present' and 'ORF' in old_transcript_id:
        continue
    gene_dict[transcript_id].add_sigp2()

for line in sigP3_lines:
    line = line.rstrip()
    split_line = line.split()
    old_transcript_id = split_line[0]
    if '.t' not in old_transcript_id:
        old_transcript_id = ORF_conversion_obj.convert_ORF_transcript(old_transcript_id)
    transcript_id = conversion_obj.convert_transcript(old_transcript_id)
    if transcript_id == 'not present' and 'ORF' in old_transcript_id:
        continue
    gene_dict[transcript_id].add_sigp3()

for line in sigP4_lines:
    line = line.rstrip()
    split_line = line.split()
    old_transcript_id = split_line[0]
    if '.t' not in old_transcript_id:
        old_transcript_id = ORF_conversion_obj.convert_ORF_transcript(old_transcript_id)
    transcript_id = conversion_obj.convert_transcript(old_transcript_id)
    if transcript_id == 'not present' and 'ORF' in old_transcript_id:
        continue
    gene_dict[transcript_id].add_sigp4()

for line in phobius_lines:
    line = line.rstrip()
    split_line = line.split()
    old_transcript_id = split_line[0]
    if '.t' not in old_transcript_id:
        old_transcript_id = ORF_conversion_obj.convert_ORF_transcript(old_transcript_id)
    transcript_id = conversion_obj.convert_transcript(old_transcript_id)
    if transcript_id == 'not present' and 'ORF' in old_transcript_id:
        continue
    gene_dict[transcript_id].add_phobius()

for line in TM_lines:
    line = line.rstrip()
    old_transcript_id = line.split("\t")[0]
    # print(old_transcript_id)
    transcript_id = conversion_obj.convert_transcript(old_transcript_id)
    # gene_dict[old_transcript_id].add_transmem()
    # print transcript_id
    gene_dict[transcript_id].add_transmem()

for line in rxlr_regex_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        old_transcript_id = split_line[0].replace('>', '')
        if '.t' not in old_transcript_id:
            old_transcript_id = ORF_conversion_obj.convert_ORF_transcript(old_transcript_id)
        transcript_id = conversion_obj.convert_transcript(old_transcript_id)
        if transcript_id == 'not present' and 'ORF' in old_transcript_id:
            continue
        gene_dict[transcript_id].add_rxlr_regex(line)

for line in rxlr_hmm_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        old_transcript_id = split_line[0].replace('>', '')
        if '.t' not in old_transcript_id:
            old_transcript_id = ORF_conversion_obj.convert_ORF_transcript(old_transcript_id)
        transcript_id = conversion_obj.convert_transcript(old_transcript_id)
        if transcript_id == 'not present' and 'ORF' in old_transcript_id:
            continue
        gene_dict[transcript_id].add_rxlr_hmm()

for line in crn_dwl_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        old_transcript_id = split_line[0].replace('>', '')
        if '.t' not in old_transcript_id:
            old_transcript_id = ORF_conversion_obj.convert_ORF_transcript(old_transcript_id)
        transcript_id = conversion_obj.convert_transcript(old_transcript_id)
        if transcript_id == 'not present' and 'ORF' in old_transcript_id:
            continue
        gene_dict[transcript_id].add_dwl_hmm()

for line in crn_lflak_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        old_transcript_id = split_line[0].replace('>', '')
        if '.t' not in old_transcript_id:
            old_transcript_id = ORF_conversion_obj.convert_ORF_transcript(old_transcript_id)
        transcript_id = conversion_obj.convert_transcript(old_transcript_id)
        if transcript_id == 'not present' and 'ORF' in old_transcript_id:
            continue
        gene_dict[transcript_id].add_lflak_hmm()

for line in effP_lines:
    line = line.rstrip()
    if 'probability' in line:
        continue
    old_transcript_id = line
    transcript_id = conversion_obj.convert_transcript(old_transcript_id)
    # transcript_id = line
    gene_dict[transcript_id].add_effp()

for line in TF_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    transcript_id = split_line[0]
    gene_dict[transcript_id].add_TF(line)

for line in phibase_lines:
    line = line.rstrip()
    transcript_id = line.split("\t")[0]
    gene_dict[transcript_id].add_phi(line)

# for line in toxin_lines:
#     line = line.rstrip()
#     transcript_id = line.split("\t")[0]
#     gene_dict[transcript_id].add_toxin(line)

cazy_dict = defaultdict(list)
for line in cazy_lines:
    if line.startswith('#'):
        continue
    line = line.rstrip()
    split_line = line.split()
    transcript_id = split_line[2]
    gene_dict[transcript_id].add_cazy(line)

for line in swissprot_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    transcript_id = split_line[0]
    gene_dict[transcript_id].add_swissprot(line)

interpro_set =  Set([])
interpro_dict = defaultdict(list)

for line in InterPro_lines:
    line = line.rstrip()
    split_line = line.split()
    transcript_id = split_line[0]
    # print transcript_id
    gene_dict[transcript_id].add_interpro(line)

isolate_list = []
for line in SNP_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    if line.startswith('##'):
        continue
    elif line.startswith('#'):
        header = line
        isolate_list = split_line[9:]
        # print isolate_list
        continue
    if 'transcript' in split_line[7] and any(x in split_line[7] for x in ['missense_variant', 'stop_gained', 'stop_lost', 'start_lost']):
        transcript_id = split_line[7].split("|")[6]
        # print transcript_id
        # print gene_dict[transcript_id]
        gene_dict[transcript_id].add_SNP(line, isolate_list)

#
for line in phase_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    phase_info = split_line[0]
    transcript_id = split_line[1]
    gene_dict[transcript_id].SNP_phase[phase_info] += 1

isolate_list = []
for line in InDel_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    if line.startswith('##'):
        continue
    elif line.startswith('#'):
        header = line
        isolate_list = split_line[9:]
        # print isolate_list
        continue
    if 'transcript' in split_line[7] and any(x in split_line[7] for x in ['frameshift_variant', 'disruptive_inframe_deletion', 'disruptive_inframe_insertion', 'inframe_deletion', 'inframe_insertion',  'stop_lost', 'start_lost', 'stop_gained', 'transcript_ablation', 'exon_loss_variant']):
        transcript_id = split_line[7].split("|")[6]
        gene_dict[transcript_id].add_InDel(line, isolate_list)

isolate_list = []
for line in SV_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    if line.startswith('##'):
        continue
    elif line.startswith('#'):
        header = line
        isolate_list = split_line[9:]
        # print isolate_list
        continue
    if 'transcript' in split_line[7] and any(x in split_line[7] for x in ['frameshift_variant', 'disruptive_inframe_deletion', 'disruptive_inframe_insertion', 'inframe_deletion', 'inframe_insertion',  'stop_lost', 'start_lost', 'stop_gained', 'transcript_ablation', 'exon_loss_variant']):
        transcript_id = split_line[7].split("|")[6]
        gene_dict[transcript_id].add_SV(line, isolate_list)

strain_id = conf.strain_id + "|"
all_isolates = conf.OrthoMCL_all
orthogroup_dict = defaultdict(lambda: "\t")
orthogroup_content_dict = defaultdict(list)
for line in orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    orthogroup_contents = []
    orthogroup_content_dict.clear()
    for isolate in all_isolates:
        num_genes = line.count((isolate + "|"))
        orthogroup_contents.append(str(isolate) + "(" + str(num_genes) + ")")
        content_counts = ":".join(orthogroup_contents)
        orthogroup_content_dict[isolate] = num_genes
    for transcript_id in split_line[1:]:
        content_str = ",".join(split_line[1:])
        if strain_id in transcript_id:
            transcript_id = transcript_id.replace(strain_id, '')
            gene_dict[transcript_id].add_orthogroup(orthogroup_id, content_counts, content_str)
    # print orthogroup_id
    expansion_status = get_orthogroup_expansion_status(orthogroup_content_dict)
    expansion_status2 = get_orthogroup_expansion_status2(orthogroup_content_dict)
    content_str = ",".join(split_line[1:])
    summary = summarise_orthogroup(content_str)
    orthogroup_dict[orthogroup_id] = "\t".join([summary, expansion_status, expansion_status2])


gene_list = sorted(gene_dict.keys(), key=lambda x:int(x[1:].split('.')[0]))

for transcript_id in gene_list:
    gene_dict[transcript_id].ipr2effectors()
    gene_id = re.sub(r".t.*", "" , transcript_id)
    # print gene_id
    if expression_dict[gene_id]:
        # print 'badgers'
        # print "\t".join([transcript_id, gene_id])
        gene_dict[transcript_id].add_expr(expression_dict[gene_id][0])
        gene_dict[transcript_id].add_LFC()
        gene_dict[transcript_id].add_early_late()

FC_conditions = ['Mycelium vs Emily 12 hpi',
    'Mycelium vs Emily 48 hpi',
    'Emily 12 hpi vs Emily 48 hpi',
    'Mycelium vs Fenella 12 hpi',
    'Mycelium vs Fenella 48 hpi',
    'Fenella 12 hpi vs Fenella 48 hpi']

print "\t".join([
    'transcript_id',
    'contig',
    'start',
    'stop',
    'strand',
    'sigp2',
    'sigp3',
    'sigp4',
    'phobius',
    'transmem',
    'secreted',
    'rxlr_regex',
    'rxlr_hmm',
    'putative_rxlr',
    'dwl_hmm',
    'lflak_hmm',
    'putative_crn',
    'effp',
    'cazy',
    'transcriptionfactor',
    'orthogroup_id',
    'content_counts',
    'content_str',
    'expansion_status',
    'phi',
    'ipr_effectors',
    'variation_level',
    'SNP_phase',
    'SNP_info',
    'InDel_info',
    'SV_info',
    'swissprot',
    'interpro',
    # "\t".join(conditions_list),
    "\t".join(gene_dict[gene_list[0]].treatments),
    # "\t".join(FC_conditions),
    "\t".join(FC_conditions),
    'DEG',
    'DEG profile',
    # 'DEG_conditions',
    'protein sequence'
    ])

for transcript_id in gene_list:
    gene_obj = gene_dict[transcript_id]
    print "\t".join([
        gene_obj.transcript_id,
        gene_obj.contig,
        gene_obj.start,
        gene_obj.stop,
        gene_obj.strand,
        gene_obj.sigp2,
        gene_obj.sigp3,
        gene_obj.sigp4,
        gene_obj.phobius,
        gene_obj.transmem,
        gene_obj.secreted,
        gene_obj.rxlr_regex,
        gene_obj.rxlr_hmm,
        gene_obj.putative_rxlr,
        gene_obj.dwl_hmm,
        gene_obj.lflak_hmm,
        gene_obj.putative_crn,
        gene_obj.effp,
        gene_obj.cazy,
        ";".join(gene_obj.transcriptionfactor),
        gene_obj.orthogroup_id,
        gene_obj.content_counts,
        orthogroup_dict[gene_obj.orthogroup_id],
        gene_obj.phi,
        '|'.join(gene_obj.ipr_effectors),
        ",".join(gene_obj.variation_level),
        gene_obj.summarise_phase(),
        "|".join(gene_obj.SNP_info),
        "|".join(gene_obj.InDel_info),
        "|".join(gene_obj.SV_info),
        gene_obj.swissprot,
        "|".join(gene_obj.interpro),
        "\t".join(gene_obj.mean_fpkm),
        # "\t".join(gene_obj.FC_values),
        "\t".join(gene_obj.FC_values_adj),
        # gene_obj.fpkm,
        gene_obj.deg,
        gene_obj.exp_time,
        # gene_obj.DEG_conditions,
        gene_obj.protein_seq
        ])
