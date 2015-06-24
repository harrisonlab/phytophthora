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

out_db = conf.db
print("Creting the merged database:\n" + out_db)
db3 = gffutils.create_db(gff_str, from_string=True, dbfn=out_db, force=True, keep_order=False, 
	sort_attribute_values='merge', merge_strategy='merge', id_spec=['ID'])

quit()

	