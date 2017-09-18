#!/usr/bin/python

'''
Make a 3-way venn diagram from a vcf file showin SNPs present in isolates.
'''

import sys,argparse
from collections import defaultdict
from sets import Set
import numpy as np
from operator import itemgetter

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument('--vcf',required=True,type=str,help='input vcf file')
ap.add_argument('--g1_name',required=True,type=str,help='name to use for group1 isolates')
ap.add_argument('--g1_isolates',required=True,nargs='+',type=str,help='list of group1 isolates')
ap.add_argument('--g2_name',required=True,type=str,help='name to use for group2 isolates')
ap.add_argument('--g2_isolates',required=True,nargs='+',type=str,help='list of group2 isolates')
ap.add_argument('--g3_name',required=True,type=str,help='name to use for group3 isolates')
ap.add_argument('--g3_isolates',required=True,nargs='+',type=str,help='list of group3 isolates')
ap.add_argument('--prefix',required=True,type=str,help='output_file prefix')


conf = ap.parse_args()

g1_name = conf.g1_name
g1_isolates = conf.g1_isolates
g2_name = conf.g2_name
g2_isolates = conf.g2_isolates
g3_name = conf.g3_name
g3_isolates = conf.g3_isolates

f_prefix = conf.prefix

with open(conf.vcf) as f:
    inp_lines = f.readlines()


#-----------------------------------------------------
# Step 2
#
#-----------------------------------------------------

GT_lines = ["\t".join([g1_name, g2_name, g3_name])]
# print("\t".join(["GT_AD_GQ_AD_ratio"]))
# print g1_isolates
# print "\t".join(g1_isolates)
# uncharacterised_SNPs = 0
for line in inp_lines:
    line = line.rstrip()
    if line.startswith('##'):
        continue
    elif line.startswith('#CHROM'):
        split_line = line.split("\t")
        # g1_cols
        # find columns for group1 isolates
        g1_index = [split_line.index(x) for x in g1_isolates]
        print "group1 name:\t" + g1_name
        print "isolates:\t" + "\t".join(g1_isolates)
        print "indexes:\t" + "\t".join(str(x) for x in g1_index) + "\n"
        # find columns for group2 isoaltes
        g2_index = [split_line.index(x) for x in g2_isolates]
        print "group2 name:\t" + g2_name
        print "isolates:\t" + "\t".join(g2_isolates)
        print "indexes:\t" + "\t".join(str(x) for x in g2_index) + "\n"
        # find columns for group2 isoaltes
        g3_index = [split_line.index(x) for x in g3_isolates]
        print "group3 name:\t" + g3_name
        print "isolates:\t" + "\t".join(g3_isolates)
        print "indexes:\t" + "\t".join(str(x) for x in g3_index) + "\n"
        continue
    split_line = line.split("\t")
    contig = split_line[0]
    location = split_line[1]
    info_col = split_line[7]
    # gene_id = info_col.split('|')[4]
    transcript_id = info_col.split('|')[6]

    group1_GT_cols = [split_line[x] for x in g1_index]
    group1_GTs = [x.split(":")[0] for x in group1_GT_cols]
    # print group1_GTs
    if "0/1" in group1_GTs:
        group1_GT = "AB"
    elif "1/1" in group1_GTs and "0/0" in group1_GTs:
        group1_GT = "AB"
    elif all(["0/0" == x for x in group1_GTs]):
        group1_GT = "AA"
        # print group1_GTs
    elif all(["1/1" == x for x in group1_GTs]):
        group1_GT = "BB"
        # print group1_GTs
    else: # remaining sites may represent polymporhic sites or sites containing missing data
        # uncharacterised_SNPs += 1
        group1_GT = "?"

    group2_GT_cols = [split_line[x] for x in g2_index]
    group2_GTs = [x.split(":")[0] for x in group2_GT_cols]
    if "0/1" in group2_GTs:
        group2_GT = "AB"
    elif "1/1" in group2_GTs and "0/0" in group2_GTs:
        group2_GT = "AB"
    elif all(["0/0" == x for x in group2_GTs]):
        group2_GT = "AA"
    elif all(["1/1" == x for x in group2_GTs]):
        group2_GT = "BB"
    else: # remaining sites may represent polymporhic sites or sites containing missing data
        # uncharacterised_SNPs += 1
        group2_GT = "?"

    group3_GT_cols = [split_line[x] for x in g3_index]
    group3_GTs = [x.split(":")[0] for x in group3_GT_cols]
    if "0/1" in group3_GTs:
        group3_GT = "AB"
    elif "1/1" in group3_GTs and "0/0" in group3_GTs:
        group3_GT = "AB"
    elif all(["0/0" == x for x in group3_GTs]):
        group3_GT = "AA"
    elif all(["1/1" == x for x in group3_GTs]):
        group3_GT = "BB"
    else: # remaining sites may represent polymporhic sites or sites containing missing data
        # uncharacterised_SNPs += 1
        group3_GT = "?"
    GT_lines.append("\t".join([transcript_id, group1_GT, group2_GT, group3_GT]))

out_genotypes = f_prefix + "_genotype_by_group.txt"
f = open(out_genotypes ,"w")
f.write("\n".join(GT_lines))
f.close()




#-----------------------------------------------------
# Step 2
#
#-----------------------------------------------------
#
# out_genotypes = f_prefix + "_genotype_by_group.txt"
# with open(out_genotypes) as f:
#     GT_lines = f.readlines()
out_dict = defaultdict(list)
monomophic_lines = [x for x in GT_lines[1:] if not "?" in x]

A = 0
B = 0
C = 0
AB = 0
AC = 0
BC = 0
ABC = 0

for line in monomophic_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    transcript_id = split_line[0]
    GT_combo = "_".join(split_line[1:])
    # print GT_combo
    if GT_combo == "AB_AA_AA":
        A += 1
        out_dict['A'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "BB_AA_AA":
        A += 1
        out_dict['A'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AA_AB_AA":
        B += 1
        out_dict['B'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AA_BB_AA":
        B += 1
        out_dict['B'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AA_AA_AB":
        C += 1
        out_dict['C'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AA_AA_BB":
        C += 1
        out_dict['C'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AB_AB_AA":
        AB += 1
        out_dict['AB'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "BB_AB_AA":
        AB += 1
        out_dict['AB'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AB_BB_AA":
        AB += 1
        out_dict['AB'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "BB_BB_AA":
        AB += 1
        out_dict['AB'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AB_AA_AB":
        AC += 1
        out_dict['AC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "BB_AA_AB":
        AC += 1
        out_dict['AC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AB_AA_BB":
        AC += 1
        out_dict['AC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "BB_AA_BB":
        AC += 1
        out_dict['AC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AA_AB_AB":
        BC += 1
        out_dict['BC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AA_BB_BB":
        BC += 1
        out_dict['BC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AA_BB_AB":
        BC += 1
        out_dict['BC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AA_AB_BB":
        BC += 1
        out_dict['BC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AB_AB_AB":
        ABC += 1
        out_dict['ABC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AB_BB_BB":
        ABC += 1
        out_dict['ABC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AB_AB_BB":
        ABC += 1
        out_dict['ABC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "AB_BB_AB":
        ABC += 1
        out_dict['ABC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "BB_AB_BB":
        ABC += 1
        out_dict['ABC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "BB_BB_AB":
        ABC += 1
        out_dict['ABC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "BB_AB_AB":
        ABC += 1
        out_dict['ABC'].append(":".join([transcript_id, GT_combo]))
    elif GT_combo == "BB_BB_BB":
        ABC += 1
        out_dict['ABC'].append(":".join([transcript_id, GT_combo]))

# print "".join(monomophic_lines)
labelA = g1_name
labelB = g2_name
labelC = g3_name

# print A

from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_unweighted, venn3_circles
plt.figure(figsize=(4,4))
# values for each sector in order of A, B, AB, C, AC, BC, ABC
v = venn3_unweighted(subsets=(A, B, AB, C, AC, BC, ABC), set_labels = (labelA, labelB, labelC))
# modify a colour of a specified sector
# v.get_patch_by_id('100').set_alpha(1.0)
# v.get_patch_by_id('100').set_color('white')
# replace a sector with specified text
# v.get_label_by_id('100').set_text('Unknown')
# v.get_label_by_id('A').set_text('Set "A"')
# formatting of lines drawn around the circles, scaled to circles of size A, B, AB, C, AC, BC, ABC
c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
# modify group1 line weight
# c[0].set_lw(1.0)
# modify group1 line style
# c[0].set_ls('dotted')
plt.title("Number of monomorphic SNPs")
# plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
#              ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
#              arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))

# plt.show()
venn_pdf = f_prefix + "_venn.pdf"
plt.savefig(venn_pdf, bbox_inches='tight')


out_lines = []
for venn_sector in ['A', 'B', 'C', 'AB', 'AC', 'BC', 'ABC']:
    # print "Genes in Sector:\t" +  venn_sector
    sector_lines = out_dict[venn_sector]
    for sector_info in sector_lines:
        split_sector_info = sector_info.split(":")
        split_sector_info.append(venn_sector)
        out_lines.append("\t".join(split_sector_info))


out_venn_sectors = f_prefix + "_genes_by_venn_sector.txt"
f = open(out_venn_sectors ,"w")
f.write("\n".join(out_lines))
f.close()
