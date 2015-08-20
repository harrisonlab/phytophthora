#!/usr/bin/python

'''
This tool can be used to transfer a note attribute from a transcript
gff feature in a gffutils database to a parent feature (gene). parents
of a transcript will receieve notes from their child transcripts.
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
ap.add_argument('--in_db',required=True,type=str,help='The gffutils database to be searched')
ap.add_argument('--out_db',required=True,type=str,help='The name of the output gff file')
conf = ap.parse_args() #sys.argv


#######################################
#      Define the transformation      #
#          functions used to          #
#          edit the database          #
#######################################

i = 0
def transform_func(x):
    if x.featuretype == 'gene':
        children = in_db.children(x, level=None, featuretype='transcript', order_by=None, reverse=False, completely_within=False, limit=None)
        for lower_feat in children:
            if 'Note' in lower_feat.attributes:
                new_notes = lower_feat.attributes['Note']
                if 'Note' in x.attributes:
                    x.attributes['Note'].extend(new_notes)
                else:
                    x.attributes['Note'] = new_notes
            global i
            i += 1
    return x



#######################################
#        Create output database       #
#                                     #
#                                     #
#######################################

f = conf.in_db
in_db = gffutils.FeatureDB(f)

# for feature in in_db.all_:

gffutils.create_db(
    in_db,
    dbfn=conf.out_db,
    force=True,
    keep_order=True,
    merge_strategy='merge',
    transform=transform_func,
    sort_attribute_values=False
    )

print("Number of features with notes added:\t" + str(i))
