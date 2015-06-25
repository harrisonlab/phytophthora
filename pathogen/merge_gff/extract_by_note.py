#!/usr/bin/python

'''
This tool can be used to search the notes attribute of a gff feature in a 
gffutils database. It will search for a given string and extract the feature 
as a gff file.
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
ap.add_argument('--str',required=True,type=str,nargs='+',help='The string or list of strings to be searched for in the attributes collumn')
ap.add_argument('--out',required=True,type=str,help='The name of the output gff file')

conf = ap.parse_args() #sys.argv
f = conf.db
regex = conf.str
out_f = open(conf.out, 'w')


#######################################
#          Open the gffutils          #
#               database              #
#                                     #
#######################################


db = gffutils.FeatureDB(f)

#######################################
#        Search for the provided      #
#               strings &             #
#         output to a gff file        #
#######################################

transcripts = db.count_features_of_type('transcript')

result = []
for feature in db.features_of_type('transcript', limit=None, strand=None, order_by=None, reverse=False, completely_within=False):
	if all(word in str(feature.attributes['Note']) for word in regex):
 		result.append(feature)
		out_f.write(str(feature) + "\n")
		parents = db.parents(feature, level=None, featuretype=None, order_by=None, reverse=False, completely_within=False, limit=None)
		for higher_feat in parents:
			print(str(higher_feat))
		children = db.children(higher_feat, level=None, featuretype=None, order_by=None, reverse=False, completely_within=False, limit=None)
		for lower_feat in children:
			print(str(lower_feat))

print("Number of Transcripts present:")
print(transcripts)
print("Number of features with both WY and RxLR annotations:")
print(len(result))		

quit()