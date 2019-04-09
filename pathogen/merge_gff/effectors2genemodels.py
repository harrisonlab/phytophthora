#!/usr/bin/python

'''
This script will run A series of python scripts to combine evidence from
a number of sources with a database of gene models
'''
import os
import sys,argparse
import gffutils
from collections import defaultdict
from itertools import chain

#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

ap.add_argument('--GffAug',required=True,type=str,help='augustus_preds.gtf')
ap.add_argument('--ORF_RxLRs',required=True,type=str,help='ORF_sp_rxlr.gff3')
ap.add_argument('--ORF_EER_IDs',required=True,type=str,help='IDs of ORF fragments containing EER domains')
ap.add_argument('--ORF_WYs',required=True,type=str,help='ORF_WY_hmmer.gff3')
ap.add_argument('--RxLR_aug',required=True,type=str,help='aug_RxLR_finder_names.txt')
ap.add_argument('--aug_EER_IDs',required=True,type=str,help='IDs of Augustus genes containing EER domains')
ap.add_argument('--WY_aug',required=True,type=str,help='aug_WY_hmmer_names.txt')

ap.add_argument('--out_db',required=True,type=str,help='Name of output database containing the merged unique features')
ap.add_argument('--out_gff',required=False,type=str,help='Name of output gff document of final database')
conf = ap.parse_args() #sys.argv

ProgDir = ~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff

# Infiles
# GffAug=gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf
# ORF_RxLRs=analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3
# ORF_WYs=analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3
# RxLR_aug=analysis/sigP_rxlr/P.cactorum/414/414_aug_RxLR_finder_names.txt
# WY_aug=analysis/hmmer/WY/P.cactorum/414/414_aug_WY_hmmer_names.txt
GffAug = conf.GffAug
ORF_RxLRs = conf.ORF_RxLRs
ORF_WYs = conf.ORF_WYs
RxLR_aug = conf.RxLR_aug
WY_aug = conf.WY_aug
ORF_EER_ID = conf.ORF_EER_IDs
Aug_EER_ID = conf.ORF_EER_IDs

# Outfiles
# AugDB=414_aug.db
# WyDB=414_WY.db
# WyID=WY_id.txt
# WyDB_mod=414_WY_note.db
# RxlrDB=414_rxlr.db
# RxlrID=rxlr_id.txt
# RxlrDB_mod=414_rxlr_note.db
# Rxlr_Wy_DB=414_rxlr_WY.db
AugDB=tmp1.db
WyDB=tmp2.db
WyID=tmp3.db
WyDB_mod=tmp4.db
RxlrDB=tmp5.db
RxlrID=tmp6.txt
RxlrDB_mod1=tmp7A.db
RxlrDB_mod2=tmp7B.db
Rxlr_Wy_DB=tmp8.db

OrfMerged=tmp9.db
MergedDB=tmp10.db
FinalDB=tmp11.db
FinalGff=tmp12.gff
FinalDB_mod1A=tmp13A.db
FinalDB_mod1B=tmp13B.db
FinalDB_mod2=conf.out_db
FinalGff_mod2=conf.out_gff


$ProgDir/make_gff_database.py --inp $GffAug --db $AugDB
$ProgDir/make_gff_database.py --inp $ORF_WYs --db $WyDB
$ProgDir/make_gff_database.py --inp $ORF_RxLRs --db $RxlrDB


## Get the IDs of all the genes in the RxLR and WY ORF databases and add notes
$ProgDir/get_db_id.py --db $WyDB --type gene --out $WyID
$ProgDir/note2db.py --in_db $WyDB --out_db $WyDB_mod --id_file $WyID --str ORF_WY_hmm --attribute ID
$ProgDir/get_db_id.py --db $RxlrDB --type gene --out $RxlrID
$ProgDir/note2db.py --in_db $RxlrDB --out_db $RxlrDB_mod1 --id_file $ORF_EER_ID --str ORF_RxLR_EER --attribute ID
$ProgDir/note2db.py --in_db $RxlrDB_mod1 --out_db $RxlrDB_mod2 --id_file $RxlrID --str ORF_RxLR_atg --attribute ID

## Merge the RxLR effector and WY effector databases together
$ProgDir/merge_db.py --inp $WyDB_mod $RxlrDB_mod --db $Rxlr_Wy_DB


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

$ProgDir/note2db.py --in_db $FinalDB --out_db $FinalDB_mod1A --id_file $RxLR_aug --str Aug_RxLR --attribute ID
$ProgDir/note2db.py --in_db $FinalDB_mod1A --out_db $FinalDB_mod1B --id_file $Aug_EER_ID --str Aug_RxLR_EER --attribute ID
$ProgDir/note2db.py --in_db $FinalDB_mod1B --out_db $FinalDB_mod2 --id_file $WY_aug --str Aug_WY_hmm --attribute ID
