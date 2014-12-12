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
