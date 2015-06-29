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
#        Serarch for overlapping      #
#                features             #
#                                     #
#######################################

d = {}
i = 0
out_lines = []

genes = db.features_of_type('gene')
# For each gene:
for gene in genes:
    ID = "".join(gene.attributes['ID'])
    # See if that gene has been previously merged
    if not ID in d.keys():
        # Otherwise, look for overlapping features on the same strand
        strand = gene.strand
        overlaps = db.region(region = gene, strand=strand, featuretype = 'gene')
        list = []
        # For each overlapping feature, add this to a dictionary of
        # all the merged features.
        # Also create a list of the current overlapping features
        for feature in overlaps:
            x = "".join(feature.attributes['ID'])
            d[x] = ""
            list.append(x)
        # Identify overlaps a second time, this time to perfrom merging.
        overlaps = db.region(region = gene, strand=strand, featuretype = 'gene')
        merged = db.merge(overlaps, ignore_strand=False)
        # For each newly created merged feature:
        for new in merged:
            # Add an ID, with a name set from stdin + itterator number.
            # Add a source, set from stdin
            # Add a note containing the IDs of the merged features.
            i += 1
            new.source = str(conf.source)
            new.attributes['ID'] = [str(conf.id) + "_" + str(i)]
            new.attributes['Note'] = ['merged(' + " ".join(list) + ')']
            out_lines.append(new)


#######################################
#        Create a database of merged  #
#                features             #
#                                     #
#######################################


merged_db = gffutils.create_db(
	out_lines,
	from_string=True,
	dbfn=conf.out,
	force=True,
	keep_order=False,
	sort_attribute_values='merge',
	merge_strategy='merge',
	id_spec=['ID']
	)

# If the switch was set on stdin, print a gff file to stdout
if s == True:
    print("##gff-version 3")
    for line in out_lines:
        print(str(line))
