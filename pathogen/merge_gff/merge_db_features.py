#!/usr/bin/python

'''
This tool will open a gffutils database and identify gene features that overlap
one another. These features will be merged when creating a new database.
'''

import sys,argparse
import gffutils
import re
from itertools import chain

#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

ap.add_argument('--inp',required=True,type=str,help='databases to input')
ap.add_argument('--id',required=True,type=str,help='The ID for newly merged features will start with this string')
ap.add_argument('--source',required=False,type=str,default='merge_db_features.py',help='A string describing the source of newly created features')
ap.add_argument('--gff',required=False,action='store_true',help='Will print the output in .gff3 format to stdout if set')
ap.add_argument('--out',required=True,type=str,help='The name of the output database.')

conf = ap.parse_args() #sys.argv

f = conf.inp
o = open(conf.out, 'w')
s = conf.gff

#######################################
#          Open the gffutils          #
#               database              #
#                                     #
#######################################

db = gffutils.FeatureDB(f)

#######################################
#        Create a function to         #
#      Serarch for overlapping        #
#             features                #
#######################################

def merge_func(func_db):
    genes = func_db.features_of_type('gene')
    d = {}
    i = 0
    out_lines = []
    # For each gene:
    for gene in genes:
        ID = "".join(gene.attributes['ID'])
        # See if that gene has been previously merged
        if not ID in d.keys():
            # Otherwise, look for overlapping features on the same strand
            strand = gene.strand
            overlaps = func_db.region(region = gene, strand=strand, featuretype = 'gene')
            list = []
            # For each overlapping feature, add this to a dictionary of
            # all the merged features.
            for feature in overlaps:
                x = "".join(feature.attributes['ID'])
                d[x] = ""
                # Create a list of the current overlapping features
                # A list may have already been created in the notes section
                # From the first itteration of the function, if this is the
                # case, then take the list of IDs from the notes attribute
                # and add them to the new list of IDs.
                dont_append = False
                prev_notes = []
#                if 'Note' in gene.attributes:
#                    prev_notes = gene.attributes['Note']
                if 'Note' in feature.attributes:
                    note_list = []
                    note_list = feature.attributes['Note']
                    for note in note_list:
                        m = re.match(r"merged\[(.*)\]", note)
                        if m and m.group(1):
                            prev_ids = m.group(1)
                            for prev_id in prev_ids.split(" "):
                                list.append(prev_id)
                                dont_append = True
                        else:
                            prev_notes.append(note)
                if dont_append == False:
                    list.append(x)
                # for transcripts in feature.children()
            # Identify overlaps a second time, this time to perfrom merging.
            overlaps = func_db.region(region = gene, strand=strand, featuretype = 'gene')
            merged = func_db.merge(overlaps, ignore_strand=False)
            # For each newly created merged feature:
            for new in merged:
                # Add an ID, with a name set from stdin + itterator number.
                # Add a source, set from stdin
                # Add a note containing the IDs of the merged features.
                i += 1
                out_list = " ".join(sorted(set(list)))
                new.source = str(conf.source)
                this_id = str(conf.id) + "_" + str(i)
                new.attributes['ID'] = [this_id]
                prev_notes.append(str('merged[' + out_list + ']'))
                new.attributes['Note'] = prev_notes
                out_lines.append(new)
                # Identify which genes were merged to make new features.
                # Maintain the old gene features as child feature transcripts
                # Within the newly created genes.
                seqid = new.seqid
                start = new.start
                stop = new.stop
                merged_features = func_db.region(region=(seqid, start, stop), strand=strand, featuretype = 'gene')
                for transcript in merged_features:
                    # print(transcript)
                    transcript.featuretype = 'transcript'
                    transcript.attributes['parents'] = [this_id]
                    out_lines.append(transcript)

    return(out_lines)


#######################################
#        Create a database of merged  #
#                features             #
#                                     #
#######################################

# Perfrom a first itteration of merging features
tmp_lines =  merge_func(db)
db2 = gffutils.create_db(
	tmp_lines,
	from_string=True,
	dbfn=':memory:',
	force=True,
	keep_order=False,
	sort_attribute_values='merge',
	merge_strategy='merge',
	id_spec=['ID']
	)

# Perform a second itteration of merging
out_lines =  merge_func(db2)
merged_db = gffutils.create_db(
	out_lines,
	from_string=True,
	dbfn=conf.out,
	force=True,
	keep_order=False,
	sort_attribute_values='merge',
	merge_strategy='merge',
	id_spec=['ID']
	)

# If the switch was set on stdin, print a gff file to stdout
if s == True:
    print("##gff-version 3")
    for line in out_lines:
        print(str(line))
