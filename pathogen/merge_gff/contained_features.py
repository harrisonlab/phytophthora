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
#ap.add_argument('--db',required=True,type=str,help='output db file')
#ap.add_argument('--source',required=False,type=str,help='(Optional) String to replace the source column - aids merging.')
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
#regex = conf.str
#out_f = open(conf.out, 'w')


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

aug = 0
atg = 0
genes = db.features_of_type('gene')
for gene in genes:
    if a in gene.source:
        aug += 1
        strand = gene.strand
        # print (strand)
        list = db.region(region = gene, strand=strand, featuretype = 'gene')
        c = 0
        # start = []
        # end = []
        # start = ""
        # end = ""
        # B_ID = ""
        B_ID = []
        B_note = []
        B_name = []
        for contained_feature in list:
            # out_lines.append(str(contained_feature))
            # print(contained_feature)
            if not a in contained_feature.source:
#                print(contained_feature)
                c += 1
                if 'ID' in contained_feature.attributes:
                    # B_ID = "".join(contained_feature.attributes['ID'])
                    B_ID.extend(contained_feature.attributes['ID'])
                if 'Note' in contained_feature.attributes:
                    B_note.extend(contained_feature.attributes['Note'])
                # if 'Name' in contained_feature.attributes:
                #     B_name.extend(contained_feature.attributes['Name'])
                # start.append(contained_feature.start)
                # end.append(contained_feature.end)
            else:
                # out_lines.append(str(contained_feature))
#                print(contained_feature)
                A_ID = "".join(contained_feature.attributes['ID'])
                # print("badgers\t" + A_ID)
        # out_lines.append("number of ORF_fragments in this gene:\t" + str(c))
        # print("number of ORF_fragments in this gene:\t" + str(c))
        # if c > 1:
        #     o_txt.write(A_ID + "\n")

        if c >= 1:
            # out_lines2.append("A consensus of these gene features may look like:")
            # print(str(B_ID))
            # gene.attributes['ID']
            merge_count += 1
            gene['Note'] = B_note
            # B_name.append(B_ID)
            gene['Name'] = B_ID
            out_lines1.append(str(A_ID))
            out_lines2.append(str(gene))
            used_ids.extend(gene['ID'])
            used_ids.extend(B_ID)

            for child in db.children(gene):
                # print(child)
                out_lines2.append(str(child))
            # a.start
            # a.end
            # start.sort()
            # end.sort()
            # print("start:\t" + str(min(start)))
            # print("end:\t" + str(max(end)))


    elif b in gene.source:
        atg += 1



#######################################
#        Identify genes from another  #
#            source contained         #
#          within these genes         #
#######################################

# print(str("\n".join(out_lines)))
# print("\n".join(out_lines2))
#


# print("Used IDs are:\t" + str(used_ids))

for k in used_ids:
    d[k].append("True")

# print(d.keys())

genes = db.features_of_type('gene')
found = 0
unmatched = 0
for gene in genes:
    this_id = "".join(gene.attributes['ID'])
    if this_id in d:
        # print(gene)
        # print("found gene:\t" + this_id)
        found += 1
    else:
        # print("unmatched gene\t\t" + this_id)
        #print(gene)
        unmatched += 1
        out_lines2.append(str(gene))
        for child in db.children(gene):
            out_lines2.append(str(child))

for line in out_lines1:
    o_txt.write(line + "\n")

db_lines = "\n".join(out_lines2)
# print(db_lines)


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

merge_num = db3.count_features_of_type(featuretype='gene')

if not s == True:
    print("The total number of Augustus genes are:\t" + str(aug))
    print("The total number of atg genes are:\t" + str(atg))
    print("Of these, this many were merged:\t" + str(found))
    print("Into this many features:\t" + str(merge_count))
    print("And this many remain unmerged:\t" + str(unmatched))
    print("The final dataset contains the following number of features:\t" + str(merge_num))

if s == True:
    print("##gff-version 3")
    for feature in db3.all_features():
        print(feature)
# print("\n".join(out_lines1))
# print(len(out_lines1))

quit()
