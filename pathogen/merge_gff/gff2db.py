#!/usr/bin/python

'''

This script builds a gffutils database from one or more gff files.
Any duplicate gff features will be combined into a single feature.
Unique notes in any of the duplicated gff features will be preserved
in the final feature. 

'''

import gffutils
from itertools import chain

#-------------------------
# Make temporary databases
#-------------------------

## Iteration for input file 1
#-------------------------
gff_fn1 = 'tmp1.gff'

db1 = gffutils.create_db(gff_fn1, dbfn='test1.db', force=True, keep_order=True, \
merge_strategy='merge', sort_attribute_values=True)

gff_feat1 = db1.all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)

## Iteration for input file 2
#-------------------------
gff_fn2 = 'tmp2.gff'

db2 = gffutils.create_db(gff_fn2, dbfn='test2.db', force=True, keep_order=True, \
merge_strategy='merge', sort_attribute_values=True)

gff_feat2 = db2.all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)

#-------------------------
# Make a combined gff database
#-------------------------

all_gff = chain(gff_feat1, gff_feat2)

db3 = gffutils.create_db(all_gff, dbfn='test3.db', force=True, keep_order=False, \
sort_attribute_values='merge', \
merge_strategy='merge', \
id_spec=['ID'])

all_gff_feat = db3.all_features(limit=None, strand=None, featuretype=None, order_by=None, reverse=False, completely_within=False)

#-------------------------
# Perform analysis of database
#-------------------------

db3.count_features_of_type('transcript')
result = []
for i in db3.features_of_type('transcript', limit=None, strand=None, order_by=None, reverse=False, completely_within=False):
	if 'RxLR_motif' in str(i.attributes['Note']) and 'WY_hmmer' in str(i.attributes['Note']) :
		result.append(i)
len(result)		

quit()