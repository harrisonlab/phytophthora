#!/usr/bin/python
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

ap.add_argument('--orthogroups', required=True,type=str,help='A file containing results of orthology analysis')
ap.add_argument('--strain_id',required=True,type=str,help='The identifier of this strain as used in the orthology analysis')
ap.add_argument('--OrthoMCL_all',required=True,type=str,nargs='+',help='The identifiers of all strains used in the orthology analysis')

conf = ap.parse_args()


with open(conf.orthogroups) as f:
    orthogroup_lines = f.readlines()

#-----------------------------------------------------
# Step 2
# Define classes
#-----------------------------------------------------


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

    total_set = set(c_A + c_C)
    # print total_set
    # print len(total_set)
    if len(total_set) == 1:
        expansion_status.append('equal')

    return(";".join(expansion_status))


#--------
#
#--------

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
    expansion_status2 = get_orthogroup_expansion_status2(orthogroup_content_dict)
    content_str = ",".join(split_line[1:])
    outline = "\t".join([orthogroup_id, expansion_status2, content_counts])
    print outline
