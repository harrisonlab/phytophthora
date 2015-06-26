

# Merge effector features from Augustus predictions

<!--
# ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
# $ProgDir/make_gff_database.py --inp tmp1.gff --db test1.db
# $ProgDir/make_gff_database.py --inp tmp2.gff --db test2.db
# $ProgDir/merge_db.py --inp test1.db test2.db --db merge3.db --source monkeys
# $ProgDir/extract_by_note.py --db merge3.db --str RxLR WY --out utils_out1.gff
#
 -->

This is working in the unmasked P. infestans genome T30-4

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

# Merge effector features from both Augustus and atg.pl predictions

This is working in the unmasked P. cactorum genome 404

```shell
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	GffAug=analysis/hmmer/WY/P.cactorum/404/404_aug_WY_hmmer.gff
	AugNames=analysis/hmmer/WY/P.cactorum/404/404_aug_WY_hmmer_names.txt
	GffAtg=analysis/hmmer/WY/P.cactorum/404/404_ORF_WY_hmmer.gff3
	AtgNames=analysis/hmmer/WY/P.cactorum/404/404_ORF_WY_hmmer_names.txt

	$ProgDir/make_gff_database.py --inp $GffAug --db aug_gff.db
	$ProgDir/note2db.py --in_db aug_gff.db --out_db aug_gff_notes.db --id_file $AugNames --str Augustus_effector
	$ProgDir/make_gff_database.py --inp $GffAtg --db atg_gff.db
	$ProgDir/note2db.py --in_db atg_gff.db --out_db atg_gff_notes.db --id_file $AtgNames --str atg_effector --attribute Name
	$ProgDir/merge_db.py --inp aug_gff_notes.db atg_gff_notes.db --db merge_aug_atg_effectors.db
	$ProgDir/extract_by_note.py --db merge_aug_atg_effectors.db --str Augustus_effector atg.pl_effector --out aug_atg_effectors.gff

```
next identify how many atg.pl transcripts are within augustus genes.

The program contained_features.py was written to print the augustus features containing
atg.pl features

```shell
	$ProgDir/contained_features.py --inp merge_aug_atg_effectors.db --out aug_genes_containing_atg.txt --A WY --B atg
```

This output a text file contianing a list of IDs that could be used to identify
augustus genes containing putative effectors predicted from atg.pl

These genes were annotated as follows:

```shell
	$ProgDir/note2db.py --in_db aug_gff_notes.db --out_db aug_gff_suppoted.db --id_file aug_genes_containing_atg.txt --str atg_effector --attribute ID
```

All genes overlapping one another in a gff database were combined using the
following commands:

```shell
	$ProgDir/merge_db_features.py --inp atg_gff.db
```
