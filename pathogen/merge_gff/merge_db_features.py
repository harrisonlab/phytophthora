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
#ap.add_argument('--out',required=True,type=str,help='Name of output text document of IDs')
conf = ap.parse_args() #sys.argv

f = conf.inp
#o = open(conf.out, 'w')



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

genes = db.features_of_type('gene')
print(db.count_features_of_type('gene'))
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
        list = db.region(region = gene, strand=strand, featuretype = 'gene')
        for feature in list:
            x = "".join(feature.attributes['ID'])
            d[x] = ""
            # print("badger")
            # print(feature)
        list = db.region(region = gene, strand=strand, featuretype = 'gene')
        merged = db.merge(list, ignore_strand=False)
        for new in merged:
            print(new)

#        print(merged)



#print(db.count_features_of_type('gene'))
