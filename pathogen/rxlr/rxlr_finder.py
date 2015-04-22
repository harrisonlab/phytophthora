#!/usr/bin/python
'''
This script searches for RxLR motifs within amino acid sequences.
It annotates the position of the RxLR sequence within the genes.
It also notes the presence  of WL motifs and DEER domains and
their positions.

The program usage is:
rxlr_finder.py [protein_file.fa] > RxLR_outfile.fa
'''

# import modules incl biopython
from Bio import SeqIO
from Bio import Motif
import sys
import re
from os import path

motif="R.LR"
#expression=

# Open fasta file
filename = sys.argv[1]
with open(filename) as file:


# For fasta accesion
	for rec in SeqIO.parse(file,"fasta"):
		# print rec.id
		seq = str(rec.seq)
# Search within the sequence for RxLR.
		match = re.search(r".{35,100}R.LR", seq)
		if match:
			print rec.description ,
			
# Note position of RxLR
			motifPos = (len(match.group()) - 4) 
			print "\t--RxLR_start: " + str(motifPos) ,

# if the sequence has an RxLR, then
	# Search for a [D]EER motif after the RxLR.
	# Note the position of the EER motif and the variant that matched the query.
			hitDEER = re.search(r"(^.*?)(D?E*EER)", str(rec.seq))
			if hitDEER: 
# 				print hitDEER.group(0)
# 				print str(len(hitDEER.group(0)))
				print "\t--" + hitDEER.group(2) + "_start: " + str(len(hitDEER.group(1))) ,
				
				#print str(len(hitDEER.group(1)))
# 				print hitDEER.group(2)
# 				print str(len(hitDEER.group(2)))
				#lenDEER = 
			#motifPos = (len(match.group()) - 4) 

	
	# Search for WL motifs within the sequence following the RxLR
	# Note the number of WL motifs and the position of each of these.
#			hitWY = [m.start() for m in re.finditer('WY', str(rec.seq))]
# 			hits = m in re.finditer(
#			if hitWY: print "\t--WY_start: " + str(hitWY).strip('[]') ,
			#m.start() for m in hitWL,
	# Retrieve the accession header 
	# Modify the header to include fields of
	#'--RxLR_start: $pos --EER_var: $str --EER_start: $pos --WL_pos: $pos;$pos;$pos'
	# print the modified header and sequence
			print ""
			print str(rec.seq)