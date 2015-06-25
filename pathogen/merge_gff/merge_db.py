#!/usr/bin/python

'''
Combine two gffutils databases into a single database
This will merge duplicate features, but preserve the unique
notes from each database.
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
ap.add_argument('--inp',required=True,type=str,nargs='+',help='list of databases to input')
ap.add_argument('--db',required=True,type=str,help='output db file')
ap.add_argument('--source',required=False,type=str,help='(Optional) String to replace the source column - aids merging.')
conf = ap.parse_args() #sys.argv


#######################################
#         Load input databases        #
#         into a single string        #
#                                     #
#######################################

all_gff = []
for f in conf.inp:
	print("Parsing input database:\n" + f)
	loaded_db = gffutils.FeatureDB(f).all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)
	for feature in loaded_db:
 		all_gff.append(str(feature))

gff_str = "\n".join(all_gff)


#######################################
#        Create the output database   #
#               from string           #
#                                     #
#######################################

new_source = conf.source
def transform_func(x):
   # Replace the source column with text from stdin
	if new_source:
		x.source = new_source
	return x
# 	if not 'ID' in x.attributes:
#  		new_id = "".join(x.attributes['Parent']) + "." + str(x.featuretype)
# 		if id_dict.has_key(new_id):
# 			id_dict[new_id] += 1
# 		else:
# 			id_dict[new_id] = 1
# 		x.attributes['ID'] = new_id + str(id_dict[new_id])
# 	return x

out_db = conf.db
print("Creting the merged database:\n" + out_db)
db3 = gffutils.create_db(
	gff_str,
	from_string=True,
	dbfn=out_db,
	force=True,
	keep_order=False, 
	sort_attribute_values='merge',
	transform=transform_func, 
	merge_strategy='merge',
	id_spec=['ID']
	)

quit()

	