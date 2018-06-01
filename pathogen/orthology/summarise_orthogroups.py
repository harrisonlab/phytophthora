#!/usr/bin/python

'''
This program is used to build information on all the genes predicted in
an annotated genome. These commands take information on location of genes
& suppliment this information with information on interproscan domains
and swissprot annotations.
'''

# summarise orthogroup (accepts orthogroup_contents and returns a string containing summarised path groups)
# CR(), LR(), Md(), Pi()


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import os
import argparse
import re
import numpy as np
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--orthogroups', required=True,type=str,help='A file containing results of orthology analysis')
conf = ap.parse_args()

with open(conf.orthogroups) as f:
    orthogroup_lines = f.readlines()

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
    # return(summarised_orthogroup)
    print "\t".join([orthogroup_id, summarised_orthogroup])



for line in orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    content_str = ",".join(split_line[1:])
    summarise_orthogroup(content_str)
    # summarised_orthogroup = summarise_orthogroup(content_str)
    # print "\t".join([orthogroup_id, summarised_orthogroup])
