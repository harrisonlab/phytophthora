#!/usr/bin/python

'''
make a gffutils db from a gff file
'''

import sys,argparse
import gffutils
i = 0

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--inp',required=True,type=str,help='input gff file')
ap.add_argument('--db',required=True,type=str,help='output db file')
conf = ap.parse_args() #sys.argv

#db = gffutils.create_db(conf.inp, force=True, dbfn=conf.db)


def transform_func(x):
    # adds some text to the end of transcript IDs
	if not 'ID' in x.attributes:
# 		i += 1
 		new_id = "".join(x.attributes['Parent']) + "." + str(x.featuretype)
# 		print(new_id)
		x.attributes['ID'] = new_id
#	x.attributes['gene_name'] = 'spoons'
	return x


db1 = gffutils.create_db(
	conf.inp,
	dbfn=conf.db,
	force=True,
	keep_order=True, 
	merge_strategy='merge',
	transform=transform_func, 
	sort_attribute_values=True
	)


