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