#!/bin/bash

#INFILE=analysis/blast_homology/P.cactorum/404/404_PHI_36_accessions_homologs.csv
#INFILE=$1
#echo "Infile is: $INFILE"

for INFILE in $@; do
#	echo "$INFILE"
	head -n1 $INFILE | cut -f1 > "$INFILE"_present.csv
	head -n1 $INFILE | cut -f1 > "$INFILE"_absent.csv
	while read line; do
#	for COL_COL in $(cat $INFILE | cut -f1 | tail -n+2); do
#	for COL_COL in $(cat $INFILE | cut -f1,1020 | tail -n+2); do
#		echo "$COL_COL"
#		LINE= split -p'\t' $line
		ID=$(printf $line | cut -d' '  -f1)
		HIT=$(echo $line | cut -d' ' -f1020)
#		echo "ID is: $ID"
#		echo "Hit is: $HIT"
		if [ "$HIT" != "0" ]; then
#			echo "hit"
			printf "$ID" >> "$INFILE"_present.csv
			printf "\n" >> "$INFILE"_present.csv
		else
#			echo "absent"
			printf "$ID" >> "$INFILE"_absent.csv
			printf "\n" >> "$INFILE"_absent.csv
		fi
	done<$INFILE
done 
exit