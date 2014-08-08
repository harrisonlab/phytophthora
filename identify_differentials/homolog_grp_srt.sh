#!/bin/bash

#INFILE=analysis/blast_homology/P.fragariae/JHVZ02/JHVZ02_appended_rxlr_sp.fa_homologs.csv


paste -d '\t' presence_tab.csv <(cut -f 2-1121 analysis/blast_homology/P.fragariae/JHVZ02/JHVZ02_appended_rxlr_sp.fa_homologs.csv) > presence_tab_grps.csv

sort -k7 presence_tab_grps.csv > presence_tab_grps_srt.csv

sed -e 's/[\t]-//g' presence_tab_grps_srt.csv > presence_tab_grps_srt_ed.csv

grep -P '\s0\s0\s0' presence_tab_grps_srt_ed.csv > absent_all.csv		# Edit this line before running
grep -P '\s1\s1\s1' presence_tab_grps_srt_ed.csv > present_all.csv		# Edit this line before running
grep -vP '\s0\s0\s0' presence_tab_grps_srt_ed.csv | grep -vP '\s1\s1\s1' > differentials.csv	# Edit this line before running
