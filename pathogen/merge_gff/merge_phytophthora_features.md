

# ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
# $ProgDir/make_gff_database.py --inp tmp1.gff --db test1.db
# $ProgDir/make_gff_database.py --inp tmp2.gff --db test2.db
# $ProgDir/merge_db.py --inp test1.db test2.db --db merge3.db --source monkeys
# $ProgDir/extract_by_note.py --db merge3.db --str RxLR WY --out utils_out1.gff
# 

# Merge effector features from Augustus predictions

```shell
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	Gff1=analysis/hmmer/WY/P.infestans/T30-4_unmasked/T30-4_unmasked_aug_WY_hmmer.gff
	WY_Names=analysis/hmmer/WY/P.infestans/T30-4_unmasked/T30-4_unmasked_aug_WY_hmmer_names.txt
	Gff2=analysis/sigP_rxlr/P.infestans/T30-4_unmasked/T30-4_unmasked_aug_RxLR_finder.gff
	RXLR_Names=analysis/sigP_rxlr/P.infestans/T30-4_unmasked/T30-4_unmasked_aug_RxLR_finder_names.txt


	$ProgDir/make_gff_database.py --inp $Gff1 --db test3.db
	$ProgDir/note2db.py --in_db test3.db --out_db test3_notes.db --id_file $WY_Names --str WY_domains
	$ProgDir/make_gff_database.py --inp $Gff2 --db test4.db
	$ProgDir/note2db.py --in_db test4.db --out_db test4_notes.db --id_file $RXLR_Names --str RxLR_motif
	$ProgDir/merge_db.py --inp test3_notes.db test4_notes.db --db merge4.db --source Augustus
	$ProgDir/extract_by_note.py --db merge4.db --str RxLR_motif WY_domain --out utils_out2.gff
```
