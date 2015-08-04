#!/bin/bash

set -e
set -u


ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
# Infiles
GffAug=$1 # ../gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf
ORF_RxLRs=$2 # ../analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3
ORF_WYs=$3 # ../analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3
RxLR_aug=$4 # ../analysis/sigP_rxlr/P.cactorum/414/414_aug_RxLR_finder_names.txt
WY_aug=$5 # ../analysis/hmmer/WY/P.cactorum/414/414_aug_WY_hmmer_names.txt
ORF_EER_ID=$7
Aug_EER_ID=$8

# Outfiles
AugDB=tmp1.db # 414_aug.db
WyDB=tmp2.db # 414_WY.db
WyID=tmp3.db # WY_id.txt
WyDB_mod=tmp4.db # 414_WY_note.db
RxlrDB=tmp5.db # 414_rxlr.db
RxlrID=tmp6.db # rxlr_id.txt
RxlrDB_mod1=tmp7A.db # 414_rxlr_note.db
RxlrDB_mod2=tmp7B.db
Rxlr_Wy_DB=tmp8.db # 414_rxlr_WY.db


OrfMerged=tmp9.db # 414_rxlr_WY_merged.db
MergedDB=tmp10A.db # 414_Aug_ORF_merged.db
MergedDB_mod=tmp10B.db
FinalDB=tmp11.db # 414_Aug_ORF.db
FinalGff=tmp12.db # 414_Aug_ORF.gff



FinalDB_mod1=tmp13.db # 414_Aug_ORF_mod1.db
FinalDB_mod2=$9 # 414_Aug_ORF_full_rxlr.db
FinalGff_mod2=$10 # 414_effectors.gff

$ProgDir/make_gff_database.py --inp $GffAug --db $AugDB
$ProgDir/make_gff_database.py --inp $ORF_WYs --db $WyDB
$ProgDir/make_gff_database.py --inp $ORF_RxLRs --db $RxlrDB


## Get the IDs of all the genes in the RxLR and WY ORF databases and add notes

$ProgDir/get_db_id.py --db $WyDB --type gene --out $WyID
$ProgDir/note2db.py --in_db $WyDB --out_db $WyDB_mod --id_file $WyID --str ORF_WY_hmm --attribute ID
$ProgDir/get_db_id.py --db $RxlrDB --type gene --out $RxlrID
$ProgDir/note2db.py --in_db $RxlrDB --out_db $RxlrDB_mod1 --id_file $RxlrID --str ORF_RxLR_atg --attribute ID
$ProgDir/note2db.py --in_db $RxlrDB_mod1 --out_db $RxlrDB_mod2 --id_file $ORF_EER_ID --str ORF_RxLR_EER --attribute ID


## Merge the RxLR effector and WY effector databases together

$ProgDir/merge_db.py --inp $WyDB_mod $RxlrDB_mod2 --db $Rxlr_Wy_DB



## Merge all features in the Rxlr and WY effector ORF database

$ProgDir/merge_db_features.py --inp $Rxlr_Wy_DB --id ORF_RxLR --source ORF_RxLR --out $OrfMerged


## Merge the effector ORF database and the augustus database together

$ProgDir/merge_db.py --inp $AugDB $OrfMerged --db $MergedDB


## Identify all effector ORFs contained within Augustus genes

$ProgDir/contained_features.py --inp $MergedDB --out_db $FinalDB --A AUGUSTUS --B RF_RxLR --out_gff $FinalGff

# The total number of Augustus genes are:	16889
# The total number of atg genes are:	1146
# Of these, this many were merged:	875
# Into this many features:	435
# And this many remain unmerged:	17160
# The final dataset contains the following number of features:	17595

## Finally, add notes to the db from effectors predicted from augustus genes


$ProgDir/note2db.py --in_db $FinalDB --out_db $FinalDB_mod1 --id_file $RxLR_aug --str Aug_RxLR --attribute ID
$ProgDir/note2db.py --in_db $FinalDB_mod1 --out_db $FinalDB_mod2 --id_file $Aug_EER_ID --str Aug_RxLR_EER --attribute ID
$ProgDir/note2db.py --in_db $FinalDB_mod2 --out_db $Rxlr_Wy_DB --id_file $WY_aug --str Aug_WY_hmm --attribute ID
$ProgDir/extract_by_note.py --db $FinalDB_mod2 --str Aug_RxLR Aug_WY_hmm ORF_WY_hmm ORF_RxLR_atg --out $FinalGff_mod2 --type gene transcript
