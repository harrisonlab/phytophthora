#!/usr/bin/python

'''
This tool can be used to add a note attribute to a gff feature in a 
gffutils database. Gff features matching a given ID will have notes 
added to them.
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
ap.add_argument('--id_file',required=True,type=str,help='A text file containing gene IDs or transcript IDs for annotation')
ap.add_argument('--str',required=True,type=str,help='Text to add to the notes section for each matched feature.')
conf = ap.parse_args() #sys.argv


dict = {}
with open(conf.id_file) as f:
    for line in f:
    	id = line.strip()
    	dict[id] = None

#######################################
#      Define the transformation      #
#          functions used to          #
#          edit the database          #
#######################################

new_note = conf.str
def transform_func(x):
   # Replace the source column with text from stdin
	if not 'Note' in in_db[id].attributes:
#		in_db[id].attributes['Note'] += "," + new_note
		x.attributes['Note'] = new_note
	else:
#		in_db[id].attributes['Note'] = "monkeys"
		x.attributes['Note'].append(new_note)
	print(str(x))
	return x

#######################################
#        Search for the provided      #
#               strings &             #
#         output to a gff file        #
#######################################

f = conf.in_db
in_db = gffutils.FeatureDB(f)


# d = {}
# new_note = conf.str
# print("note to add:\t" + new_note)
# with open(conf.id_file) as f:
#     for line in f:
#     	id = line.strip()
#     	print(id)
#     	print(in_db[id].attributes['ID'])
#     	print(in_db[id].attributes.items())
#     	in_db[id].start = 23
#     	print(in_db[id].start)
#     	print(in_db[id].attributes.items())    	
# #    	in_db[id].attributes['Parent'].append('monkeys')
# #    	print(in_db[id].attributes['Note'])
# #    	print(in_db[id].attributes['Parent'])
# #    	in_db[id].attributes['ID'].append('monkeys')
# #    	print(str(in_db[id].attributes['ID'].append('monkeys')))
# # #   	print(in_db[id])
# #     	if in_db[id].attributes['ID']:
# #     		in_db[id].attributes['ID'] += "," + new_note
# #     	else:
# #     		in_db[id].attributes['ID'] = "monkeys"
# #     		#new_note
# #    	print(in_db[id])


gffutils.create_db(
	in_db,
	dbfn=conf.out_db,
	force=True,
	keep_order=True, 
	merge_strategy='merge',
	transform=transform_func, 
	sort_attribute_values=True
	)