

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
	$ProgDir/merge_db_features.py --inp atg_gff.db --id atg_WY --source WY_ORFs --out WY_ORF_merged.db
```

The number of atg.pl WY features contained within WY augustus genes was identified
the following commands:
```shell
	$ProgDir/get_db_id.py --db WY_ORF_merged.db --type gene --out WY_ORF_merged_id.txt
	$ProgDir/note2db.py --in_db WY_ORF_merged.db --out_db WY_ORF_merged_notes.db --id_file WY_ORF_merged_id.txt --str WY_ORFs --attribute ID
	$ProgDir/merge_db.py --inp aug_gff_notes.db WY_ORF_merged_notes.db --db merge_aug_atg_effectors2.db
	$ProgDir/contained_features.py --inp merge_aug_atg_effectors2.db --out_db final_out.db --A WY_hmmer --B WY_ORFs --out_gff tmpX.gff
```



# Merge 404 ORF RxLR motif effectors with Augustus genes
This process was repeated on ORF models containing RxLRs in all Augustus gene
models using the following commands:

## Set variables
```shell
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	# Infiles
	GffAug=gene_pred/augustus_unmasked/P.cactorum/404/404_augustus_preds.gtf
	ORF_RxLRs=analysis/rxlr_atg_unmasked/P.cactorum/404/404_ORF_sp_rxlr.gff3
	# Outfiles
	AugDB=404_aug.db
	OrfDB=404_ORF.db
	OrfMerged=404_ORF_rxlr_merged.db
	MergedDB=404_Aug_ORF_merged.db
	FinalDB=404_Aug_ORF_RxLR.db
	FinalGff=404_Aug_ORF_RxLR.gff
```

## Make a db of aug genes and effector ORFs
```shell
	$ProgDir/make_gff_database.py --inp $GffAug --db $AugDB
	$ProgDir/make_gff_database.py --inp $ORF_RxLRs --db $OrfDB
```

## Merge all features in the effector ORF database
```shell
	$ProgDir/merge_db_features.py --inp $OrfDB --id ORF_RxLR --source atg_RxLR --out $OrfMerged
```

## Merge the effector ORF database and the augustus database together
```shell
	$ProgDir/merge_db.py --inp $AugDB $OrfMerged --db $MergedDB
```

## Identify all effector ORFs contained within Augustus genes
```shell
	$ProgDir/contained_features.py --inp $MergedDB --out_db $FinalDB --A AUGUSTUS --B atg_RxLR --out_gff $FinalGff
```

```shell
	$ProgDir/get_db_id.py --db WY_ORF_merged.db --type gene --out WY_ORF_merged_id.txt
	$ProgDir/note2db.py --in_db WY_ORF_merged.db --out_db WY_ORF_merged_notes.db --id_file WY_ORF_merged_id.txt --str WY_ORFs --attribute ID
	$ProgDir/merge_db.py --inp aug_gff_notes.db WY_ORF_merged_notes.db --db merge_aug_atg_effectors2.db
	$ProgDir/contained_features.py --inp merge_aug_atg_effectors2.db --out_db final_out.db --A WY_hmmer --B WY_ORFs --out_gff tmpX.gff
```
The total number of Augustus genes are:	16452
The total number of atg genes are:	967
Of these, this many were merged:	783
Into this many features:	390
And this many remain unmerged:	16636
The final dataset contains the following number of features:	17026



# Merge 414 ORF RxLR motif effectors and WY domain effectors with Augustus genes
This process was repeated on ORF models containing RxLRs in all Augustus gene
models using the following commands:

## Set variables
```shell
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	# Infiles
	GffAug=gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf
	ORF_RxLRs=analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3
	ORF_WYs=analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3
	# Outfiles
	AugDB=414_aug.db
	WyDB=414_WY.db
	WyID=WY_id.txt
	WyDB_mod=414_WY_note.db
	RxlrDB=414_rxlr.db
	RxlrID=rxlr_id.txt
	RxlrDB_mod=414_rxlr_note.db
	Rxlr_Wy_DB=414_rxlr_WY.db

	OrfMerged=414_rxlr_WY_merged.db
	MergedDB=414_Aug_ORF_merged.db
	FinalDB=414_Aug_ORF.db
	FinalGff=414_Aug_ORF.gff
```

## Make a db of aug genes and effector ORFs
```shell
	$ProgDir/make_gff_database.py --inp $GffAug --db $AugDB
	$ProgDir/make_gff_database.py --inp $ORF_WYs --db $WyDB
	$ProgDir/make_gff_database.py --inp $ORF_RxLRs --db $RxlrDB
```

## Get the IDs of all the genes in the RxLR and WY ORF databases and add notes
```shell
	$ProgDir/get_db_id.py --db $WyDB --type gene --out $WyID
	$ProgDir/note2db.py --in_db $WyDB --out_db $WyDB_mod --id_file $WyID --str ORF_WY_hmm --attribute ID
	$ProgDir/get_db_id.py --db $RxlrDB --type gene --out $RxlrID
	$ProgDir/note2db.py --in_db $RxlrDB --out_db $RxlrDB_mod --id_file $RxlrID --str ORF_RxLR_atg --attribute ID
```

## Merge the RxLR effector and WY effector databases together
```shell
	$ProgDir/merge_db.py --inp $WyDB_mod $RxlrDB_mod --db $Rxlr_Wy_DB
```


## Merge all features in the Rxlr and WY effector ORF database
```shell
	$ProgDir/merge_db_features.py --inp $Rxlr_Wy_DB --id ORF_RxLR --source ORF_RxLR --out $OrfMerged
```

## Merge the effector ORF database and the augustus database together
```shell
	$ProgDir/merge_db.py --inp $AugDB $OrfMerged --db $MergedDB
```

## Identify all effector ORFs contained within Augustus genes
```shell
	$ProgDir/contained_features.py --inp $MergedDB --out_db $FinalDB --A AUGUSTUS --B RF_RxLR --out_gff $FinalGff
```
The total number of Augustus genes are:	16889
The total number of atg genes are:	1146
Of these, this many were merged:	875
Into this many features:	435
And this many remain unmerged:	17160
The final dataset contains the following number of features:	17595

## Finally, add notes to the db from effectors predicted from augustus genes

```shell
	RxLR_aug=analysis/sigP_rxlr/P.cactorum/414/414_aug_RxLR_finder_names.txt
	WY_aug=analysis/hmmer/WY/P.cactorum/414/414_aug_WY_hmmer_names.txt
	FinalDB_mod1=414_Aug_ORF_mod1.db
	FinalDB_mod2=414_Aug_ORF_full_rxlr.db

	$ProgDir/note2db.py --in_db $FinalDB --out_db $FinalDB_mod1 --id_file $RxLR_aug --str Aug_RxLR --attribute ID
	$ProgDir/note2db.py --in_db $FinalDB_mod1 --out_db $FinalDB_mod2 --id_file $WY_aug --str Aug_WY_hmm --attribute ID

	$ProgDir/extract_by_note.py --db $FinalDB_mod2 --str Aug_RxLR Aug_WY_hmm ORF_WY_hmm ORF_RxLR_atg --out 414_effectors.gff --type gene transcript
```




# Building into a pipeline
```shell
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	qsub qsub_effectors2genemod.sh gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3 analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3 analysis/sigP_rxlr/P.cactorum/414/414_aug_RxLR_finder_names.txt analysis/hmmer/WY/P.cactorum/414/414_aug_WY_hmmer_names.txt
```


<!-- $ProgDir/effectors2genemodels.sh \
../gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf \
../analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3 \
../analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3 \
../analysis/sigP_rxlr/P.cactorum/414/414_aug_RxLR_finder_names.txt \
../analysis/hmmer/WY/P.cactorum/414/414_aug_WY_hmmer_names.txt \
414_Aug_ORF_full_rxlr.db 414_effectors.gff \
/ 2>&1 | tee logfile.txt -->

##414
```shell
	screen -a
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	ProjDir=/home/groups/harrisonlab/project_files/idris
	WorkDir=$ProjDir/analysis/databases/P.cactorum/414
	mkdir -p $WorkDir
	cd $WorkDir

	$ProgDir/effectors2genemodels.sh \
	$ProjDir/gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf \
	$ProjDir/analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3 \
	$ProjDir/analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3 \
	$ProjDir/analysis/sigP_rxlr/P.cactorum/414/414_aug_RxLR_finder_names.txt \
	$ProjDir/analysis/hmmer/WY/P.cactorum/414/414_aug_WY_hmmer_names.txt \
	414_Aug_ORF_full_rxlr.db \
	414_effectors.gff \
	2>&1 | tee 404_logfile.txt
```

The total number of Augustus genes are: 16889
The total number of atg genes are:      1146
Of these, this many were merged:        875
Into this many features:        435
And this many remain unmerged:  17160
The final dataset contains the following number of features:    17595
No. IDs in infile:      115
Number of features with notes added:    115

```shell
	cat *_effectors.gff | grep 'gene' | wc -l
```
Including this many effectors:	1148

##404
```shell
	screen -a
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	ProjDir=/home/groups/harrisonlab/project_files/idris
	WorkDir=$ProjDir/analysis/databases/P.cactorum/404
	mkdir -p $WorkDir
	cd $WorkDir

	$ProgDir/effectors2genemodels.sh \
	$ProjDir/gene_pred/augustus_unmasked/P.cactorum/404/404_augustus_preds.gtf \
	$ProjDir/analysis/rxlr_atg_unmasked/P.cactorum/404/404_ORF_sp_rxlr.gff3 \
	$ProjDir/analysis/hmmer/WY/P.cactorum/404/404_ORF_WY_hmmer.gff3 \
	$ProjDir/analysis/sigP_rxlr/P.cactorum/404/404_aug_RxLR_finder_names.txt \
	$ProjDir/analysis/hmmer/WY/P.cactorum/404/404_aug_WY_hmmer_names.txt \
	404_Aug_ORF_full_rxlr.db \
	404_effectors.gff \
	2>&1 | tee 404_logfile.txt
```
The total number of Augustus genes are: 16452
The total number of atg genes are:      1112
Of these, this many were merged:        855
Into this many features:        426
And this many remain unmerged:  16709
The final dataset contains the following number of features:    17135

```shell
	cat *_effectors.gff | grep 'gene' | wc -l
```
Including this many effectors:	1118


##10300
```shell
	screen -a
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	ProjDir=/home/groups/harrisonlab/project_files/idris
	WorkDir=$ProjDir/analysis/databases/P.cactorum/10300
	mkdir -p $WorkDir
	cd $WorkDir

	$ProgDir/effectors2genemodels.sh \
	$ProjDir/gene_pred/augustus_unmasked/P.cactorum/10300/10300_augustus_preds.gtf \
	$ProjDir/analysis/rxlr_atg_unmasked/P.cactorum/10300/10300_ORF_sp_rxlr.gff3 \
	$ProjDir/analysis/hmmer/WY/P.cactorum/10300/10300_ORF_WY_hmmer.gff3 \
	$ProjDir/analysis/sigP_rxlr/P.cactorum/10300/10300_aug_RxLR_finder_names.txt \
	$ProjDir/analysis/hmmer/WY/P.cactorum/10300/10300_aug_WY_hmmer_names.txt \
	10300_Aug_ORF_full_rxlr.db \
	10300_effectors.gff \
	2>&1 | tee 10300_logfile.txt
```

The total number of Augustus genes are: 17539
The total number of atg genes are:      1223
Of these, this many were merged:        967
Into this many features:        481
And this many remain unmerged:  17795
The final dataset contains the following number of features:    18276

```shell
	cat *_effectors.gff | grep 'gene' | wc -l
```
Including this many effectors:	1228

##371
```shell
	screen -a
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	ProjDir=/home/groups/harrisonlab/project_files/idris
	WorkDir=$ProjDir/analysis/databases/P.ideai/371
	mkdir -p $WorkDir
	cd $WorkDir

	$ProgDir/effectors2genemodels.sh \
	$ProjDir/gene_pred/augustus_unmasked/P.ideai/371/371_augustus_preds.gtf \
	$ProjDir/analysis/rxlr_atg_unmasked/P.ideai/371/371_ORF_sp_rxlr.gff3 \
	$ProjDir/analysis/hmmer/WY/P.ideai/371/371_ORF_WY_hmmer.gff3 \
	$ProjDir/analysis/sigP_rxlr/P.ideai/371/371_aug_RxLR_finder_names.txt \
	$ProjDir/analysis/hmmer/WY/P.ideai/371/371_aug_WY_hmmer_names.txt \
	371_Aug_ORF_full_rxlr.db \
	371_effectors.gff \
	2>&1 | tee 371_logfile.txt
```

The total number of Augustus genes are: 15675
The total number of atg genes are:      1093
Of these, this many were merged:        781
Into this many features:        388
And this many remain unmerged:  15987
The final dataset contains the following number of features:    16375

```shell
	cat *_effectors.gff | grep 'gene' | wc -l
```
Including this many effectors:	1092

##T30-4
```shell
	screen -a
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	ProjDir=/home/groups/harrisonlab/project_files/idris
	WorkDir=$ProjDir/analysis/databases/P.infestans/T30-4
	mkdir -p $WorkDir
	cd $WorkDir

	$ProgDir/effectors2genemodels.sh \
	$ProjDir/gene_pred/augustus_unmasked/P.infestans/T30-4_unmasked/T30-4_augustus_preds.gtf \
	$ProjDir/analysis/rxlr_atg_unmasked/P.infestans/T30-4_unmasked/T30-4_ORF_sp_rxlr.gff3 \
	$ProjDir/analysis/hmmer/WY/P.infestans/T30-4_unmasked/T30-4_unmasked_ORF_WY_hmmer.gff3 \
	$ProjDir/analysis/sigP_rxlr/P.infestans/T30-4_unmasked/T30-4_unmasked_aug_RxLR_finder_names.txt \
	$ProjDir/analysis/hmmer/WY/P.infestans/T30-4_unmasked/T30-4_unmasked_aug_WY_hmmer_names.txt \
	T30-4_Aug_ORF_full_rxlr.db \
	T30-4_effectors.gff \
	2>&1 | tee T30-4_logfile.txt
```

The total number of Augustus genes are: 38789
The total number of atg genes are:      3400
Of these, this many were merged:        2192
Into this many features:        1093
And this many remain unmerged:  39997
The final dataset contains the following number of features:    41090

```shell
	cat *_effectors.gff | grep 'gene' | wc -l
```
Including this many effectors:
