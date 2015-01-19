for directory in $(ls -d analysis/rxlr_unmasked/P.*/*); do
	STRAIN=$(echo $directory | rev | cut -d '/' -f1 | rev)
	cat $directory/"$STRAIN"_F_atg_ORF.gff | sed 's/ORF_/ORF_F_/g' > $directory/"$STRAIN"_F_atg_ORF_ed.gff
	cat $directory/"$STRAIN"_R_atg_ORF.gff | sed 's/_RC//g' | sed 's/ORF_/ORF_R_/g' > $directory/"$STRAIN"_R_atg_ORF_ed.gff
	cat $directory/"$STRAIN"_F_atg_ORF_ed.gff $directory/"$STRAIN"_R_atg_ORF_ed.gff > $directory/"$STRAIN"_ORF_ed.gff
done

# gff files supplied but they are not currently of use as the feature IDs do not match the names of the genes in the input fasta.

set -- 10300 404 414 JHVZ02 371
for a; do 
	shift
	for b; do
	printf "%s - %s\n" "$a" "$b"
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/sub_inparanoid.sh analysis/rxlr_unmasked/P.*/$a/"$a"_sp_rxlr.fa analysis/rxlr_unmasked/P.*/$b/"$b"_sp_rxlr.fa analysis/rxlr_unmasked/P.*/$a/"$a"_ORF_ed.gff analysis/rxlr_unmasked/P.*/$b/"$b"_ORF_ed.gff
	done 
done


mkdir analysis/inparanoid/summary_tables
cat analysis/inparanoid/*/*_seqs.txt | cut -f1 | sort | uniq | less
cat analysis/inparanoid/*/*_seqs.txt | cut -f1 | sort | uniq > analysis/inparanoid/summary_tables/all_genes.txt

/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_orthology_tab.pl analysis/inparanoid/summary_tables/all_genes.txt analysis/inparanoid/*/sqltable.* > analysis/inparanoid/summary_tables/rxlr_orthology_tab.csv

cat analysis/inparanoid/10300-404/*_seqs.txt analysis/inparanoid/10300-414/*_seqs.txt analysis/inparanoid/404-414/*_seqs.txt | cut -f1 | sort | uniq > analysis/inparanoid/summary_tables/all_P.cac_RxLR.txt
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_orthology_tab.pl analysis/inparanoid/summary_tables/all_P.cac_RxLR.txt analysis/inparanoid/10300-404/sqltable.10300-404 analysis/inparanoid/10300-414/sqltable.10300-414 analysis/inparanoid/404-414/sqltable.404-414 > analysis/inparanoid/summary_tables/rxlr_orthology_P.cact.csv

# 404 = 1024 genes, 1016 ortholog groups - 8 genes different
# 414 = 1048 genes, 1036 ortholog groups - 12 genes different
# 10300 = 1144 genes, 1054 ortholog groups - 90 genes different

# Identify which genes are present as inparalogs.
# 404 = 1x3copy 2x2copy
# 414 = 1x3copy, 4x2copy
# 10300 = 45x2 copy
for FILE in $(ls analysis/inparanoid/summary_tables/orthogroups4/*.txt); do 
	cat $FILE | cut -f1 -d '|' | uniq -dc; printf "$FILE\n"
done | grep -B1 '414' | less

# identify the number of genes that are represented by orthogroups from 10300
cat analysis/inparanoid/summary_tables/new_final_tab.csv analysis/rxlr_unmasked/P.cactorum/10300/10300_sp_rxlr.fa | grep -o '10300.*' | cut -f1 | sed s'/\.txt//' | sed s'/_singleton_insl.*bp//' | sed 's/_Pcac10300//' | sed 's/_contig//' | sed 's/ //' | tail -n+2 | sort | uniq -d | wc -l
# 1054 = this number will not take account of genes that are duplicated within an orthogroup.
# identify the number of genes that are not present in the orthogroups in the outfile.
cat analysis/inparanoid/summary_tables/new_final_tab.csv analysis/rxlr_unmasked/P.cactorum/10300/10300_sp_rxlr.fa | grep -o '10300.*' | cut -f1 | sed s'/\.txt//' | sed s'/_singleton_insl.*bp//' | sed 's/_Pcac10300//' | sed 's/_contig//' | sed 's/ //' | tail -n+2 | sort | uniq -u | wc -l
# 90 genes different

# Understanding the difference:
# 404 = 1x3copy, 2x2copy ~ 4 genes missing from analysis & 4 genes overpredicted 
# 414 = 1x3copy, 4x2copy ~ 6 genes missing from analysis and 6 genes overpredicted
# 10300 = 45x2 copy ~ 45 genes missing from analysis and 45 genes overpredicted


