#!/usr/bin/python

'''
This tool will open a gffutils database and identify gene features that overlap
one another. These features will be merged when creating a new database.
'''

import sys,argparse
import gffutils
from itertools import chain

#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

ap.add_argument('--inp',required=True,type=str,help='databases to input')
ap.add_argument('--id',required=True,type=str,help='The ID for newly merged features will start with this string')
ap.add_argument('--source',required=False,type=str,default='merge_db_features.py',help='A string describing the source of newly created features')
ap.add_argument('--gff',required=False,action='store_true',help='Will print the output in .gff3 format to stdout if set')
ap.add_argument('--out',required=True,type=str,help='The name of the output database.')

conf = ap.parse_args() #sys.argv

f = conf.inp
o = open(conf.out, 'w')
s = conf.gff

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
d = {}
i = 0
out_lines = []

genes = db.features_of_type('gene')
#print(db.count_features_of_type('gene'))
# merged = db.merge(genes, ignore_strand=False)
# for feature in merged:
#       print(feature)
for gene in genes:
    ID = "".join(gene.attributes['ID'])
#    print(ID)
#    print(d.keys())
    if not ID in d.keys():
#        print("bagers")
        strand = gene.strand
        overlaps = db.region(region = gene, strand=strand, featuretype = 'gene')
        list = []
        for feature in overlaps:
            x = "".join(feature.attributes['ID'])
            d[x] = ""
            list.append(x)
            # print("badger")
            # print(feature)
        overlaps = db.region(region = gene, strand=strand, featuretype = 'gene')
        merged = db.merge(overlaps, ignore_strand=False)
        for new in merged:
            i += 1
            new.source = str(conf.source)
            new.attributes['ID'] = [str(conf.id) + "_" + str(i)]
            new.attributes['Note'] = ['merged(' + " ".join(list) + ')']
            out_lines.append(new)
            if s == True:
                print(new)


db3 = gffutils.create_db(
	out_lines,
	from_string=True,
	dbfn=conf.out,
	force=True,
	keep_order=False,
	sort_attribute_values='merge',
#	transform=transform_func,
	merge_strategy='merge',
	id_spec=['ID']
	)
#write(out_lines, "\n")

print(db3.count_features_of_type('gene'))




#print(db.count_features_of_type('gene'))
