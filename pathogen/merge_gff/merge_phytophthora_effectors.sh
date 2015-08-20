#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G

# #$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace

set -u
set -e
#set -o pipefail

# These commands were performed to merge gff evidence of Putative RxLR effectors
# from fusarium spp. isolates.

# Merge Fus2 ORF and Augustus RxLR motif effectors, WY hmm model effectors and RxLR hmm model effectors

## Convert ORF predictions into correct gff3 format

ORF_Gff=$1
Aug_Gff=$2


Organism=$(echo $Aug_Gff | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Aug_Gff | rev | cut -f2 -d '/' | rev)
WorkDir=$TMPDIR/merge_effectors
mkdir -p $WorkDir

CurPath=$PWD
OutDir=$CurPath/analysis/database/$Organism/$Strain
mkdir -p $OutDir


echo "$Organism"
echo "$Strain"

ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
# ORF_Gff=gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF.gff
ORF_Gff_mod=$WorkDir/"$Strain"_ORF_mod.gff
$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod


## Extract names of effectors from each source

### Augustus genes identified as putative effectors

# Extracting RxLR Regex genes

echo "Parsing files"

Aug_Regex_RxLR_FA=analysis/sigP_rxlr/$Organism/$Strain/"$Strain"_aug_RxLR_EER_regex.fa
# Aug_Regex_RxLR=analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/"$Strain"_Aug_RxLR_regex_names.txt
# Aug_Regex_RxLR_EER=analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/"$Strain"_Aug_RxLR_EER_regex_names.txt
Aug_Regex_RxLR=$WorkDir/"$Strain"_Aug_RxLR_regex_names.txt
Aug_Regex_RxLR_EER=$WorkDir/"$Strain"_Aug_RxLR_EER_regex_names.txt
cat $Aug_Regex_RxLR_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $Aug_Regex_RxLR
cat $Aug_Regex_RxLR_FA | grep '>' | grep 'EER_motif' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $Aug_Regex_RxLR_EER

echo "1"

# Extracting WY hmm domain containing genes


Aug_hmm_WY_FA=analysis/hmmer/WY/$Organism/$Strain/"$Strain"_aug_WY_hmmer_out.fa

# Aug_hmm_WY=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain/"$Strain"_Aug_WY_hmmer_names.txt
Aug_hmm_WY=$WorkDir/"$Strain"_Aug_WY_hmmer_names.txt
ls $Aug_hmm_WY_FA
echo $Aug_hmm_WY
cat $Aug_hmm_WY_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $Aug_hmm_WY

echo "5"

Aug_hmm_CRN_FA=analysis/hmmer/CRN/$Organism/$Strain/"$Strain"_Aug_CRN_hmmer.fa
# Aug_hmm_WY=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain/"$Strain"_Aug_WY_hmmer_names.txt
Aug_hmm_CRN=$WorkDir/"$Strain"_Aug_CRN_hmmer_names.txt
cat $Aug_hmm_CRN_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $Aug_hmm_CRN


echo "2"
# Extracting RxLR hmm domain containing genes


Aug_hmm_RxLR_FA=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"*_Aug_RxLR_hmmer.fa
# Aug_hmm_RxLR=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_Aug_RxLR_hmmer_names.txt
Aug_hmm_RxLR=$WorkDir/"$Strain"_Aug_RxLR_hmmer_names.txt
cat $Aug_hmm_RxLR_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $Aug_hmm_RxLR
echo "3"

### ORF fragments identified as putative effectors

# Extracting RxLR containing ORFs

echo "4"
ORF_Regex_RxLR_FA=analysis/rxlr_atg_unmasked/$Organism/$Strain/"$Strain"_ORF_RxLR_EER_regex.fa
# ORF_Regex_RxLR=analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/"$Strain"_ORF_RxLR_regex_names.txt
# ORF_Regex_RxLR_EER=analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/"$Strain"_ORF_RxLR_EER_regex_names.txt
ORF_Regex_RxLR=$WorkDir/"$Strain"_ORF_RxLR_regex_names.txt
ORF_Regex_RxLR_EER=$WorkDir/"$Strain"_ORF_RxLR_EER_regex_names.txt
cat $ORF_Regex_RxLR_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $ORF_Regex_RxLR
cat $ORF_Regex_RxLR_FA | grep '>' | grep 'EER_motif' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $ORF_Regex_RxLR_EER


# Extracting WY hmm domain containing ORFs
echo "5"

ORF_hmm_WY_FA=analysis/hmmer/WY/$Organism/$Strain/"$Strain"_ORF_WY_hmmer.fa
# ORF_hmm_WY=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain/"$Strain"_ORF_WY_hmmer_names.txt
ORF_hmm_WY=$WorkDir/"$Strain"_ORF_WY_hmmer_names.txt
cat $ORF_hmm_WY_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $ORF_hmm_WY

echo "5"

ORF_hmm_CRN_FA=analysis/hmmer/CRN/$Organism/$Strain/"$Strain"_ORF_CRN_hmmer.fa
# ORF_hmm_WY=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain/"$Strain"_ORF_WY_hmmer_names.txt
ORF_hmm_CRN=$WorkDir/"$Strain"_ORF_CRN_hmmer_names.txt
cat $ORF_hmm_CRN_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $ORF_hmm_CRN

# Extracting RxLR hmm domain containing ORFs
echo "6"

ORF_hmm_RxLR_FA=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_ORF_RxLR_hmmer.fa
# ORF_hmm_RxLR=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_ORF_RxLR_hmmer_names.txt
ORF_hmm_RxLR=$WorkDir/"$Strain"_ORF_RxLR_hmmer_names.txt
cat $ORF_hmm_RxLR_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $ORF_hmm_RxLR


echo "setting infile variables"

## Set variables

ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff

# Infiles
# Aug_Gff=gene_pred/augustus/$Organism/$Strain/"$Strain"_augustus_preds.gtf
# ORF_Gff_mod=gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF_mod.gff

Aug_Regex_RxLR=$WorkDir/"$Strain"_Aug_RxLR_regex_names.txt
Aug_Regex_RxLR_EER=$WorkDir/"$Strain"_Aug_RxLR_EER_regex_names.txt
Aug_hmm_WY=$WorkDir/"$Strain"_Aug_WY_hmmer_names.txt
Aug_hmm_CRN=$WorkDir/"$Strain"_Aug_CRN_hmmer_names.txt
Aug_hmm_RxLR=$WorkDir/"$Strain"_Aug_RxLR_hmmer_names.txt
Aug_Mimp_1500=analysis/mimps/$Organism/$Strain/"$Strain"_mimps_intersected_Aug_genes_names.txt

ORF_Regex_RxLR=$WorkDir/"$Strain"_ORF_RxLR_regex_names.txt
ORF_Regex_RxLR_EER=$WorkDir/"$Strain"_ORF_RxLR_EER_regex_names.txt
ORF_hmm_WY=$WorkDir/"$Strain"_ORF_WY_hmmer_names.txt
ORF_hmm_CRN=$WorkDir/"$Strain"_ORF_CRN_hmmer_names.txt
ORF_hmm_RxLR=$WorkDir/"$Strain"_ORF_RxLR_hmmer_names.txt
ORF_Mimp_1500=analysis/mimps/$Organism/$Strain/"$Strain"_mimps_intersected_ORF_genes_names.txt

# Outfiles

echo "Setting outfile variables"
AugDB=$WorkDir/"$Strain"_Aug.db
OrfDB=$WorkDir/"$Strain"_ORF.db

Aug_RxLR_DB=$WorkDir/"$Strain"_Aug_RxLR_rxlr.db
Aug_RxLR_EER_DB=$WorkDir/"$Strain"_Aug_RxLR_EER_rxlr.db
Aug_hmm_WY_DB=$WorkDir/"$Strain"_Aug_RxLR_EER_WY_rxlr.db
Aug_hmm_CRN_DB=$WorkDir/"$Strain"_Aug_RxLR_EER_WY_CRN_rxlr.db
Aug_hmm_RxLR_DB=$WorkDir/"$Strain"_Aug_RxLR_EER_WY_CRN_RxLR_rxlr.db
# Aug_mimp_DB=$WorkDir/"$Strain"_Aug_RxLR_EER_WY_RxLR_CRN_mimps_rxlr.db
Aug_named_parents_DB=$WorkDir/"$Strain"_Aug_annotated_parents.db

ORF_RxLR_DB=$WorkDir/"$Strain"_ORF_RxLR_rxlr.db
ORF_RxLR_DB_EER=$WorkDir/"$Strain"_ORF_RxLR_EER_rxlr.db
ORF_hmm_WY_DB=$WorkDir/"$Strain"_ORF_RxLR_EER_WY_rxlr.db
ORF_hmm_CRN_DB=$WorkDir/"$Strain"_ORF_RxLR_EER_WY_CRN_rxlr.db
ORF_hmm_RxLR_DB=$WorkDir/"$Strain"_ORF_RxLR_EER_WY_CRN_RxLR_rxlr.db
# ORF_mimp_DB=$WorkDir/"$Strain"_ORF_RxLR_EER_WY_CRN_RxLR_mimps_rxlr.db
ORF_named_parents_DB=$WorkDir/"$Strain"_ORF_annotated_parents.db

ORF_Out_Gff=$WorkDir/"$Strain"_ORF_effectors.gff

Orf_combined_DB=$WorkDir/"$Strain"_ORF_combined.db
Orf_merged_DB=$WorkDir/"$Strain"_ORF_merged.db

Combined_DB=$WorkDir/"$Strain"_Aug_ORF_combined.db
Merged_DB=$WorkDir/"$Strain"_Aug_ORF_merged.db
Merged_Gff=$WorkDir/"$Strain"_Aug_ORF_merged.gff
Effectors_Gff=$WorkDir/"$Strain"_effectors.gff
Effectors_high_conf_Gff=$WorkDir/"$Strain"_high_conf_effectors.gff

## Make a db of aug genes and effector ORFs

$ProgDir/make_gff_database.py --inp $Aug_Gff --db $AugDB
$ProgDir/make_gff_database.py --inp $ORF_Gff_mod --db $OrfDB
# $ProgDir/make_gff_database.py --inp $ORF_WYs --db $WyDB
# $ProgDir/make_gff_database.py --inp $ORF_RxLRs --db $RxlrDB




## Get the IDs of all the genes in the RxLR and WY ORF databases and add notes
echo "Step 1"
echo "Get the IDs of all the genes in the RxLR and WY ORF databases and add notes"
# Add notes to Augustus genes with effector evidence


# $ProgDir/get_db_id.py --db $WyDB --type gene --out $WyID

$ProgDir/note2db.py --in_db $AugDB --out_db $Aug_RxLR_DB --id_file $Aug_Regex_RxLR --str Aug_RxLR_motif --attribute ID
$ProgDir/note2db.py --in_db $Aug_RxLR_DB --out_db $Aug_RxLR_EER_DB --id_file $Aug_Regex_RxLR_EER --str Aug_RxLR_EER_motif --attribute ID
$ProgDir/note2db.py --in_db $Aug_RxLR_EER_DB --out_db $Aug_hmm_WY_DB --id_file $Aug_hmm_WY --str Aug_WY_hmm --attribute ID
$ProgDir/note2db.py --in_db $Aug_hmm_WY_DB --out_db $Aug_hmm_CRN_DB --id_file $Aug_hmm_CRN --str Aug_CRN_hmm --attribute ID
$ProgDir/note2db.py --in_db $Aug_hmm_CRN_DB --out_db $Aug_hmm_RxLR_DB --id_file $Aug_hmm_RxLR --str Aug_RxLR_hmm --attribute ID
$ProgDir/note2db.py --in_db $Aug_hmm_RxLR_DB --out_db $Aug_mimp_DB --id_file $Aug_Mimp_1500 --str Aug_mimp_intersect --attribute ID
$ProgDir/notes2parents.py --in_db $Aug_mimp_DB --out_db $Aug_named_parents_DB


# Add notes to ORF fragments with effector evidence

echo "Step 2"
echo "Add notes to ORF fragments with effector evidence"

$ProgDir/note2db.py --in_db $OrfDB --out_db $ORF_RxLR_DB --id_file $ORF_Regex_RxLR --str ORF_RxLR_motif --attribute Name
$ProgDir/note2db.py --in_db $ORF_RxLR_DB --out_db $ORF_RxLR_DB_EER --id_file $ORF_Regex_RxLR_EER --str ORF_RxLR_EER_motif --attribute Name
$ProgDir/note2db.py --in_db $ORF_RxLR_DB_EER --out_db $ORF_hmm_WY_DB --id_file $ORF_hmm_WY --str ORF_WY_hmm --attribute Name
$ProgDir/note2db.py --in_db $ORF_hmm_WY_DB --out_db $ORF_hmm_CRN_DB --id_file $ORF_hmm_CRN --str ORF_CRN_hmm --attribute Name
$ProgDir/note2db.py --in_db $ORF_hmm_CRN_DB --out_db $ORF_hmm_RxLR_DB --id_file $ORF_hmm_RxLR --str ORF_RxLR_hmm --attribute Name
$ProgDir/note2db.py --in_db $ORF_hmm_RxLR_DB --out_db $ORF_mimp_DB --id_file $ORF_Mimp_1500 --str ORF_mimp_intersect --attribute ID
$ProgDir/notes2parents.py --in_db $ORF_mimp_DB --out_db $ORF_named_parents_DB


$ProgDir/extract_by_note.py --db $ORF_named_parents_DB --str ORF_RxLR_motif ORF_RxLR_EER_motif ORF_WY_hmm ORF_CRN_hmm ORF_RxLR_hmm --out $ORF_Out_Gff --type gene transcript
$ProgDir/make_gff_database.py --inp $ORF_Out_Gff --db $Orf_combined_DB


## Merge all features in the ORF database

echo "Step 3"
echo "Merge all features in the ORF database"

$ProgDir/merge_db_features.py --inp $Orf_combined_DB --id ORF_finder --source ORF_finder --out $Orf_merged_DB


## Merge the effector ORF database and the augustus database together
echo "Step 4"
echo "Merge the effector ORF database and the augustus database together"

$ProgDir/merge_db.py --inp $Aug_named_parents_DB $Orf_merged_DB --db $Combined_DB


## Identify all effector ORFs contained within Augustus genes
echo "Step 5"
echo "Identify all effector ORFs contained within Augustus genes"

$ProgDir/contained_features.py --inp $Combined_DB --out_db $Merged_DB --A AUGUSTUS --B ORF_finder --out_gff $Merged_Gff

#
# The total number of Augustus genes are: 11055
# The total number of atg genes are:      732
# Of these, this many were merged:        662
# Into this many features:        326
# And this many remain unmerged:  11125
# The final dataset contains the following number of features:    11451
#
# The total number of Augustus genes are:	11055
# The total number of atg genes are:	732
# Of these, this many were merged:	662
# Into this many features:	326
# And this many remain unmerged:	11125
# The final dataset contains the following number of features:	11451
echo "Building a final gff file"

echo "extracting all effectors:"
$ProgDir/extract_by_note.py --db $Merged_DB --str \
Aug_RxLR_motif  \
Aug_RxLR_EER_motif \
Aug_WY_hmm \
Aug_CRN_hmm \
Aug_RxLR_hmm \
# Aug_mimp_intersect \
ORF_RxLR_motif \
ORF_RxLR_EER_motif \
ORF_WY_hmm \
ORF_CRN_hmm \
ORF_RxLR_hmm \
# ORF_mimp_intersect \
--out $Effectors_Gff --type gene transcript

for i in 1; do
  printf "Aug_RxLR_motif:\t"
  cat $WorkDir/Fus2_effectors.gff | grep 'Aug_RxLR_motif' | wc -l
  printf "Aug_RxLR_EER_motif:\t"
  cat $WorkDir/Fus2_effectors.gff | grep 'Aug_RxLR_EER_motif' | wc -l
  printf "Aug_WY_hmm:\t"
  cat $WorkDir/Fus2_effectors.gff | grep 'Aug_WY_hmm' | wc -l
  printf "Aug_RxLR_hmm:\t"
  cat $WorkDir/Fus2_effectors.gff | grep 'Aug_RxLR_hmm' | wc -l
  printf "Aug_mimp_intersect:\t"
  cat $WorkDir/Fus2_effectors.gff | grep 'Aug_mimp_intersect' | wc -l
  printf "ORF_RxLR_motif:\t"
  cat $WorkDir/Fus2_effectors.gff | grep 'ORF_RxLR_motif' | wc -l
  printf "ORF_RxLR_EER_motif:\t"
  cat $WorkDir/Fus2_effectors.gff | grep 'ORF_RxLR_EER_motif' | wc -l
  printf "ORF_WY_hmm:\t"
  cat $WorkDir/Fus2_effectors.gff | grep 'ORF_WY_hmm' | wc -l
  printf "ORF_WY_hmm:\t"
  cat $WorkDir/Fus2_effectors.gff | grep 'ORF_RxLR_hmm' | wc -l
  printf "ORF_mimp_intersect:\t"
  cat $WorkDir/Fus2_effectors.gff | grep 'ORF_mimp_intersect' | wc -l
done


echo "extracting high confidence effectors:"
$ProgDir/extract_by_note.py --db $Merged_DB --str \
Aug_RxLR_motif  \
Aug_RxLR_EER_motif \
Aug_WY_hmm \
Aug_RxLR_hmm \
# Aug_mimp_intersect \
# ORF_RxLR_motif \
ORF_RxLR_EER_motif \
ORF_WY_hmm \
ORF_RxLR_hmm \
ORF_mimp_intersect \
--out $Effectors_high_conf_Gff --type gene transcript

for i in 1; do
  printf "Aug_RxLR_motif:\t"
  cat $Effectors_high_conf_Gff  | grep 'Aug_RxLR_motif' | wc -l
  printf "Aug_RxLR_EER_motif:\t"
  cat $Effectors_high_conf_Gff  | grep 'Aug_RxLR_EER_motif' | wc -l
  printf "Aug_WY_hmm:\t"
  cat $Effectors_high_conf_Gff  | grep 'Aug_WY_hmm' | wc -l
  printf "Aug_RxLR_hmm:\t"
  cat $Effectors_high_conf_Gff  | grep 'Aug_RxLR_hmm' | wc -l
  printf "Aug_mimp_intersect:\t"
  cat $Effectors_high_conf_Gff  | grep 'Aug_mimp_intersect' | wc -l
  printf "ORF_RxLR_motif:\t"
  cat $Effectors_high_conf_Gff  | grep 'ORF_RxLR_motif' | wc -l
  printf "ORF_RxLR_EER_motif:\t"
  cat $Effectors_high_conf_Gff  | grep 'ORF_RxLR_EER_motif' | wc -l
  printf "ORF_WY_hmm:\t"
  cat $Effectors_high_conf_Gff  | grep 'ORF_WY_hmm' | wc -l
  printf "ORF_WY_hmm:\t"
  cat $Effectors_high_conf_Gff  | grep 'ORF_RxLR_hmm' | wc -l
  printf "ORF_mimp_intersect:\t"
  cat $Effectors_high_conf_Gff  | grep 'ORF_mimp_intersect' | wc -l
done



echo "Final step"
echo "Copying files"

cp $WorkDir/. $OutDir/.


#
# Aug_RxLR_motif:	13
# Aug_RxLR_EER_motif:	1
# Aug_WY_hmm:	0
# Aug_RxLR_hmm:	2
# Aug_mimp_intersect:	2
# ORF_RxLR_motif:	564
# ORF_RxLR_EER_motif:	16
# ORF_WY_hmm:	4
# ORF_WY_hmm:	30
# ORF_mimp_intersect:	124
