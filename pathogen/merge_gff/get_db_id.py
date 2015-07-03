#!/usr/bin/python

'''
This tool will extract all the IDs in from a gffutils database
'''

import sys,argparse
import gffutils
from itertools import chain


#######################################
#           Load sys. args.           #
#                                     #
#                                     #
#######################################


ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--db',required=True,type=str,help='The gffutils database to be searched')
ap.add_argument('--type',required=True,type=str,help='The feature type to extract IDs from')
ap.add_argument('--out',required=True,type=str,help='The name of the output text file')

conf = ap.parse_args() #sys.argv
f = conf.db
t = conf.type
o = open(conf.out, 'w')

#######################################
#          Open the gffutils          #
#               database              #
#                                     #
#######################################

db = gffutils.FeatureDB(f)

#######################################
#          Open the gffutils          #
#               database              #
#                                     #
#######################################


num_feat = db.count_features_of_type(t)

print("Searching for features of type:\t" + t)
print("Number of features:\t" + str(num_feat))
print("Extracting to:\t" + conf.out)

features = db.features_of_type(t, limit=None, strand=None, order_by=None, reverse=False, completely_within=False)


for feature in features:
    ID = "".join(feature.attributes['ID'])
    o.write(str(ID) + "\n")
