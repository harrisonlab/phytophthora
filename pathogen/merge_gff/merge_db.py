#!/usr/bin/python

'''
Combine two gffutils databases into a single database
This will merge duplicate features, but preserve the unique
notes from each database
'''

import sys,argparse
import gffutils
from itertools import chain

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,nargs='+',help='list of databases to input')
ap.add_argument('--db',required=True,type=str,help='output db file')
conf = ap.parse_args() #sys.argv

#-------------------------
# Make a combined gff database
#-------------------------
# 
# feat_list = []
# 
# in_list = conf.inp
# for infile in in_list:
# 	print(infile)
# 	db = gffutils.FeatureDB(infile)
# 	print(db)
# 	gff_feat = db.all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)
# 	print(gff_feat)
# 	feat_list.append(gff_feat)
# 
# print(feat_list)

# db1 = conf.inp[0]
#db1 = gffutils.FeatureDB(conf.inp[0])
# #db1 = gffutils.FeatureDB(in_list[0])
#print(db1)
# db2 = conf.inp[1]
# db2 = gffutils.FeatureDB(conf.inp[1])
# #db2 = gffutils.FeatureDB(in_list[1])
# print(db2)
out_db = conf.db
print(out_db)
# 
# 
# 
#gff_feat1 = db1.all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)
# print(gff_feat1)
# gff_feat2 = db2.all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)
# print(gff_feat1)

feat_list = []
for infile in conf.inp:
	gff_feat = gffutils.FeatureDB(infile).all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)
	feat_list.append(gff_feat)
# gff_feat1 = gffutils.FeatureDB(conf.inp[0]).all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)
# gff_feat2 = gffutils.FeatureDB(conf.inp[1]).all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)



# for gff in feat_list:
# 	for feature in gff:
# 		print(feature)

# all_gff = chain(gff_feat1, gff_feat2)
# all_gff = chain(feat_list)
# for feature in feat_list:
# 	print(str(feature))

all_gff = []
for gff in chain(feat_list):
 	for feature in gff:
 		all_gff.append(str(feature))
# print(all_gff)
gff_str = "\n".join(all_gff)
#print(gff_str)

db3 = gffutils.create_db(gff_str, from_string=True, dbfn=out_db, force=True, keep_order=False, \
sort_attribute_values='merge', \
merge_strategy='merge', \
id_spec=['ID'])



#-------------------------
# Perform analysis of database
#-------------------------

all_gff_feat = db3.all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)

# for line in all_gff_feat:
# 	print(str(line))

db3.count_features_of_type('transcript')
result = []
for i in db3.features_of_type('transcript', limit=None, strand=None, order_by=None, reverse=False, completely_within=False):
	if 'RxLR_motif' in str(i.attributes['Note']) and 'WY_hmmer' in str(i.attributes['Note']) :
		result.append(i)
print("Number of features with both WY and RxLR annotations:")
print(len(result))		

quit()	
	