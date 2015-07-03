#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace
 set -u
 set -e

IN1=$1
IN2=$2
IN3=$3
IN4=$4
IN5=$5

GFF_AUG=$(basename "$1") # ../gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf
ORF_RXLR=$(basename "$2") # ../analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3
ORF_WY=$(basename "$3") # ../analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3
RXLR_AUG=$(basename "$4") # ../analysis/sigP_rxlr/P.cactorum/414/414_aug_RxLR_finder_names.txt
WY_AUG=$(basename "$5") # ../analysis/hmmer/WY/P.cactorum/414/414_aug_WY_hmmer_names.txt

ORGANISM=$(echo $IN1 | rev | cut -d "/" -f3 | rev)
STRAIN=$(echo $IN1 | rev | cut -d "/" -f2 | rev)

FINAL_DB="$STRAIN"_Aug_ORF_full_rxlr.db # 414_Aug_ORF_full_rxlr.db
FINAL_GFF="$STRAIN"_effectors.gff  # 414_effectors.gff
LOGFILE="$STRAIN"_logfile.txt

PROG_DIR=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
CUR_PATH=$PWD

WORK_DIR=$TMPDIR/effectors_"$STRAIN"
mkdir -p $WORK_DIR
cd $WORK_DIR
cp $CUR_PATH/$IN1 $GFF_AUG
cp $CUR_PATH/$IN2 $ORF_RXLR
cp $CUR_PATH/$IN3 $ORF_WY
cp $CUR_PATH/$IN4 $RXLR_AUG
cp $CUR_PATH/$IN5 $WY_AUG

$PROG_DIR/effectors2genemodels.sh $GFF_AUG $ORF_RXLR $ORF_WY $RXLR_AUG $WY_AUG \
$FINAL_DB $FINAL_GFF / 2>&1 | tee $LOGFILE

OUTDIR=$CUR_PATH/analysis/databases/$ORGANISM/$STRAIN

mkdir -p $OUTDIR
cp $FINAL_DB $OUTDIR/$FINAL_DB
cp $FINAL_GFF $OUTDIR/$FINAL_DB
cp logfile.txt $OUTDIR/$LOGFILE
