#!/bin/bash

#INFILE=analysis/blast_homology/P.cactorum/404/404_PHI_36_accessions_homologs.csv
#INFILE=$1
#echo "Infile is: $INFILE"

for INFILE in $@; do
#	echo "$INFILE"
	head -n1 $INFILE | cut -f1 > "$INFILE"_present.csv
	head -n1 $INFILE | cut -f1 > "$INFILE"_absent.csv
	head -n1 $INFILE | cut -f1 > "$INFILE"_differentials.csv	
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
			printf "$ID""\t1\n" >> "$INFILE"_differentials.csv
		else
#			echo "absent"
			printf "$ID" >> "$INFILE"_absent.csv
			printf "\n" >> "$INFILE"_absent.csv
			printf "$ID""\t0\n" >> "$INFILE"_differentials.csv
		fi
	done<$INFILE
done 

NUMBER=1
cut -f1 "$1" > tmp_tab"$NUMBER".csv

for INFILE in $@; do
	NEXT_NUM=$((NUMBER+1))
	paste -d '\t' tmp_tab"$NUMBER".csv <(cut -f2 "$INFILE"_differentials.csv) > tmp_tab"$NEXT_NUM".csv
	NUMBER=$((NUMBER+1))
done
mv tmp_tab"$NEXT_NUM".csv presence_tab.csv
rm tmp_tab*


# while read line; do
# 	if [ $(printf "$line" | grep "\t0\t0\t0\t0") ]; then
# 		printf "$line""\n" >> absent_all.csv
# 	elif $(printf "$line" | grep "\t1\t1\t1\t1"); then
# 		printf "$line""\n" >> present_all.csv
# 	elif $(printf "$line" | grep "\t?\t?\t?\t?"); then
# 		printf "$line""\n" >> differentials.csv
# 	else
# 		echo "problem line:"
# 		echo "$line"
# 	fi	
# done;<presence_tab.csv

grep -P '\s0\s0\s0\s0' presence_tab.csv > absent_all.csv
grep -P '\s1\s1\s1\s1' presence_tab.csv > present_all.csv
grep -vP '\s0\s0\s0\s0' presence_tab.csv | grep -vP '\s1\s1\s1\s1' > differentials.csv
exit