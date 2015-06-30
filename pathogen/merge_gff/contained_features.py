#!/usr/bin/python

'''
This tool will open a gffutils database and identify gene features from a
particular source. It will then count the number of genes from other sources
contained within this gene.
'''

import sys,argparse
import gffutils
from collections import defaultdict
from itertools import chain

#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

ap.add_argument('--inp',required=True,type=str,help='databases to input')
ap.add_argument('--out_txt',required=True,type=str,help='Name of output text document of IDs')
ap.add_argument('--out_db',required=True,type=str,help='Name of output database containing the merged unique features')
ap.add_argument('--gff',required=False,action='store_true',help='Will print the output in .gff3 format to stdout if set')
ap.add_argument('--A',required=True,type=str,help='Genes that may encapsulate other genes are from source A')
ap.add_argument('--B',required=True,type=str,help='Source of genes that may be encapsulated by genes from A')
conf = ap.parse_args() #sys.argv

f = conf.inp
o_txt = open(conf.out_txt, 'w')
a = conf.A
b = conf.B
# out_lines = []
out_lines1 = []
out_lines2 = []
used_ids = []
d = defaultdict(list)
s = conf.gff
merge_count = 0


#######################################
#          Open the gffutils          #
#               database              #
#                                     #
#######################################

db = gffutils.FeatureDB(f)

#######################################
#        Extract all genes from the   #
#            specified source         #
#                                     #
#######################################

# Initialise counts of genes from source A (aug) and B (atg)
aug = 0
atg = 0
genes = db.features_of_type('gene')
for gene in genes:
    # If from source specified at stdin
    if a in gene.source:
        aug += 1
        # Identify the strand and find all genes intersecting this gene
        # on this strand
        strand = gene.strand
        list = db.region(region = gene, strand=strand, featuretype = 'gene')
        # Initialise a count of intersecting genes
        c = 0
        # Initiase variable to store attributes of genes from source B
        B_ID = []
        B_note = []
        # For each intersected gene:
        for contained_feature in list:
            # out_lines.append(str(contained_feature))
            # print(contained_feature)
            # If not from source A:
            if not a in contained_feature.source:
                # Count the number of B genes intersected
                c += 1
                # Extract their attributes for transfer to a merged feature
                if 'ID' in contained_feature.attributes:
                    B_ID.extend(contained_feature.attributes['ID'])
                if 'Note' in contained_feature.attributes:
                    B_note.extend(contained_feature.attributes['Note'])
            # If from source A:
            else:
                # Store the IDs of genes used for intersection
                # out_lines.append(str(contained_feature))
                A_ID = "".join(contained_feature.attributes['ID'])
        # out_lines.append("number of ORF_fragments in this gene:\t" + str(c))

        # If a gene from A contains 1 or more genes from B:
        if c >= 1:
            # Count the number of features with intersects
            merge_count += 1
            # Modify the current gene from A to contain attributes of B
            gene['Note'] = B_note
            gene['Name'] = B_ID
            # Add the modified gene to a list for output (ID and db feature)
            out_lines1.append(str(A_ID))
            out_lines2.append(str(gene))
            # Add the IDs of the used gene and intersected genes to a list
            # to record the genes that have been merged
            used_ids.extend(gene['ID'])
            used_ids.extend(B_ID)
            # for sub-features in the gene:
            for child in db.children(gene):
                # Add these lines to the list for output to a db
                out_lines2.append(str(child))

    # Count the number of genes from source B
    elif b in gene.source:
        atg += 1



#######################################
#        Extract genes that had       #
#          no intersections or        #
#          were not intersected       #
#######################################

# make a dictionary of the gene ids used to create merged features
for k in used_ids:
    d[k].append("True")


# Initialise counts for merged and unmerged genes
found = 0
unmatched = 0

# Go through the all the genes in the database again
genes = db.features_of_type('gene')
# For each gene:
for gene in genes:
    # Check if the Id is present in the dictionary of merged genes
    this_id = "".join(gene.attributes['ID'])
    if this_id in d:
        found += 1
    # If gene wasn't intersect by a gene in A:
    else:
        unmatched += 1
        # Append it and its sub-features to a list of lines to be input
        # into a new db
        out_lines2.append(str(gene))
        for child in db.children(gene):
            out_lines2.append(str(child))


#######################################
#            Print output             #
#                                     #
#                                     #
#######################################

# Output IDs of A genes which have been used to merge:
for line in out_lines1:
    o_txt.write(line + "\n")

# Create a db from lines of merged genes, their sub-features
# and from genes that did not intersect with any other genes
# (along with their sub-features)
db_lines = "\n".join(out_lines2)
db3 = gffutils.create_db(
	db_lines,
	from_string=True,
	dbfn=conf.out_db,
	force=True,
	keep_order=False,
	sort_attribute_values='merge',
	merge_strategy='merge',
	id_spec=['ID']
	)


# if the gff switch was set:
if s == True:
    # Print a gff output of the new database to stdout
    print("##gff-version 3")
    for feature in db3.all_features():
        print(feature)
# If it wasn't set
else:
    # Print summary statistics
    print("The total number of Augustus genes are:\t" + str(aug))
    print("The total number of atg genes are:\t" + str(atg))
    print("Of these, this many were merged:\t" + str(found))
    print("Into this many features:\t" + str(merge_count))
    print("And this many remain unmerged:\t" + str(unmatched))
    # Count the number of genes following merging
    # (Including genes without intersections)
    merge_num = db3.count_features_of_type(featuretype='gene')
    print("The final dataset contains the following number of features:\t" + str(merge_num))


quit()
