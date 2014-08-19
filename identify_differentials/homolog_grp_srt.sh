#!/bin/bash

# Commands to process the presence_tabs.csv file of blast_differentials.sh script. 
# Collects blast_self data from any of the initial blast_pipe.sh output files 
# and orders the presence_tabs.csv file by these homolog groups. The groupings 
# are collapsed, so that they are more readable and then the file is split into 
# queries that are present in all genomes, absent in all genomes or are differential.

$USAGE="homolog_grp_srt.sh <presence_tab.csv> <any_blast_pipe_outfile.csv>";

echo "$USAGE"
PRESENCE_TAB=$1
HOMOLOG_GRP_TAB=$2
echo "homolog_grp_srt.sh is running using:"
echo "$PRESENCE_TAB"
echo "$HOMOLOG_GRP_TAB"

#	IMPORTANT
#
#	Edit this next line -- cut -f 2-<final_column_of_homolog_grps_in_$2> $HOMOLOG_GRP_TAB
#
paste -d '\t' $PRESENCE_TAB <( cut -f 2-1121 $HOMOLOG_GRP_TAB | tail -n +2) > presence_tab_grps.csv

sort -k7 presence_tab_grps.csv > presence_tab_grps_srt.csv

sed -e 's/[\t]-//g' presence_tab_grps_srt.csv > presence_tab_grps_srt_ed.csv

grep -P '\s0\s0\s0' presence_tab_grps_srt_ed.csv > absent_all.csv		# Edit this line before running
grep -P '\s1\s1\s1' presence_tab_grps_srt_ed.csv > present_all.csv		# Edit this line before running
grep -vP '\s0\s0\s0' presence_tab_grps_srt_ed.csv | grep -vP '\s1\s1\s1' > differentials.csv	# Edit this line before running
