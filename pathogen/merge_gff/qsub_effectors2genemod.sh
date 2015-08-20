#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G
#$ -l h=blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace
 set -u
 set -e
 set -o pipefail

IN1=gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf
IN2=analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3
IN3=analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3
IN4=analysis/sigP_rxlr/P.cactorum/414/414_aug_RxLR_finder_names.txt
IN5=analysis/hmmer/WY/P.cactorum/414/414_aug_WY_hmmer_names.txt


# IN1=$1
# IN2=$2
# IN3=$3
# IN4=$4
# IN5=$5

GFF_AUG=$(basename "$1") # ../gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf
GFF_ORF_RXLR=$(basename "$2") # ../analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3
GFF_ORF_WY=$(basename "$3") # ../analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3

TXT_RXLR_AUG=$(basename "$4") # ../analysis/sigP_rxlr/P.cactorum/414/414_aug_RxLR_finder_names.txt
TXT_WY_AUG=$(basename "$5") # ../analysis/hmmer/WY/P.cactorum/414/414_aug_WY_hmmer_names.txt

IN6=analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3
IN7=analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3
GFF_ORF_RXLR_HMM=$(basename $6)
TXT_RXLR_AUG_HMM=$(basename $7)

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
cp $CUR_PATH/$IN2 $GFF_ORF_RXLR
cp $CUR_PATH/$IN3 $GFF_ORF_WY
cp $CUR_PATH/$IN4 $TXT_RXLR_AUG
cp $CUR_PATH/$IN5 $TXT_WY_AUG

$PROG_DIR/effectors2genemodels.sh $GFF_AUG $GFF_ORF_RXLR $GFF_ORF_WY $TXT_RXLR_AUG $TXT_WY_AUG \
$FINAL_DB $FINAL_GFF / 2>&1 | tee $LOGFILE

OUTDIR=$CUR_PATH/analysis/databases/$ORGANISM/$STRAIN

mkdir -p $OUTDIR
cp $FINAL_DB $OUTDIR/$FINAL_DB
cp $FINAL_GFF $OUTDIR/$FINAL_DB
cp logfile.txt $OUTDIR/$LOGFILE
