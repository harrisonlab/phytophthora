#!/usr/bin/python

'''
This tool can be used to search the notes attribute of a gff feature in a
gffutils database. It will search for a given string and extract the feature
as a gff file.
'''

import sys,argparse
import gffutils
from itertools import chain
from collections import defaultdict

#######################################
#           Load sys. args.           #
#                                     #
#                                     #
#######################################


ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--db',required=True,type=str,help='The gffutils database to be searched')
ap.add_argument('--str',required=True,type=str,nargs='+',help='The string or list of strings to be searched for in the attributes collumn')
ap.add_argument('--out',required=True,type=str,help='The name of the output gff file')
ap.add_argument('--type',required=True,type=str,nargs='+',help='features of type will be searched for provided IDs')

conf = ap.parse_args() #sys.argv
f = conf.db
regex = conf.str
out_f = open(conf.out, 'w')
out_f.write("##gff-version 3\n")
type_list= conf.type
d = defaultdict(list)

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

for t in type_list:
	transcripts = db.count_features_of_type(t)
	print("Number of features of type \'" + t + "\':\t")
	print(transcripts)

print("looking for the following regex:")
print(regex)


ID_list = []
for t in type_list:
	print("looking in features of type:\t" + t)
	# result = []
	for feature in db.features_of_type(t, limit=None, strand=None, order_by=None, reverse=False, completely_within=False):
		if 'Note' in feature.attributes:
			# if all(word in str(feature.attributes['Note']) for word in regex):
			k = "".join(feature['ID'])
			# print(feature.attributes['Note'])
			if not k in d:
					for this_regex in regex:
						if all(word in str(feature.attributes['Note']) for word in this_regex):
					# if all(word in str(feature.attributes['Note']) for word in regex):
					# if (word in "\t".join(feature.attributes['Note']) for word in regex):
					# ID_list.append(feature['ID'])
					# 	result.append(feature.ID)
						# higher_feat = feature
							if t == 'gene':
								# print(feature)
								out_f.write(str(feature) + "\n")
								k = "".join(feature['ID'])
								d[k].append("True")
								children = db.children(feature, level=None, featuretype=None, order_by=None, reverse=False, completely_within=False, limit=None)
							else:
								parents = db.parents(feature, level=None, featuretype=None, order_by=None, reverse=False, completely_within=False, limit=None)
								for higher_feat in parents:
									out_f.write(str(higher_feat) + "\n")
									# print(higher_feat)
									k = "".join(higher_feat['ID'])
									d[k].append("True")
									children = db.children(higher_feat, level=None, featuretype=None, order_by=None, reverse=False, completely_within=False, limit=None)
							for lower_feat in children:
								out_f.write(str(lower_feat) + "\n")
								# print(lower_feat)
								k = "".join(lower_feat['ID'])
								d[k].append("True")
							break

# print(ID_list)

# for ID in ID_list:


# print("Number of features with both WY and RxLR annotations:")
# print(len(result))

quit()
