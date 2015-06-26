#!/usr/bin/python

'''
This tool will open a gffutils database and identify gene features from a
particular source. It will then count the number of genes from other sources
contained within this gene.
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

ap.add_argument('--inp',required=True,type=str,help='databases to input')
ap.add_argument('--out',required=True,type=str,help='Name of output text document of IDs')
ap.add_argument('--A',required=True,type=str,help='Genes that may encapsulate other genes are from source A')
ap.add_argument('--B',required=True,type=str,help='Source of genes that may be encapsulated by genes from A')
#ap.add_argument('--db',required=True,type=str,help='output db file')
#ap.add_argument('--source',required=False,type=str,help='(Optional) String to replace the source column - aids merging.')
conf = ap.parse_args() #sys.argv

f = conf.inp
o = open(conf.out, 'w')
a = conf.A
b = conf.B
#regex = conf.str
#out_f = open(conf.out, 'w')


#######################################
#          Open the gffutils          #
#               database              #
#                                     #
#######################################

db = gffutils.FeatureDB(f)

#######################################
#        Extract all genes from the   #
#            specified source         #
#                                     #
#######################################

aug = 0
atg = 0
genes = db.features_of_type('gene')
for gene in genes:
    if a in gene.source:
        aug += 1
        strand = gene.strand
        # print (strand)
        list = db.region(region = gene, strand=strand, featuretype = 'gene')
        c = 0
        start = []
        end = []
        for contained_feature in list:
#            print(contained_feature)
            if not a in contained_feature.source:
                print(contained_feature)
                c += 1
                # start.append(contained_feature.start)
                # end.append(contained_feature.end)
            else:
                print(contained_feature)
                ID = "".join(contained_feature.attributes['ID'])
        print("number of ORF_fragments in this gene:\t" + str(c))
        if c > 1:
            o.write(ID + "\n")
        # if c > 1:
        #     print("A consensus of these gene features may look like:")
        #     start.sort()
        #     end.sort()
        #     print("start:\t" + str(min(start)))
        #     print("end:\t" + str(max(end)))


    elif b in gene.source:
        atg += 1



#######################################
#        Identify genes from another  #
#            source contained         #
#          within these genes         #
#######################################

print("The total number of Augustus genes are:\t" + str(aug))
print("The total number of atg genes are:\t" + str(atg))
