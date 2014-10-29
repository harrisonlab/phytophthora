#!/usr/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=1G

# Commands used to identify putative crinkler genes in phytophthora genomes


INFILE=$1
MOTIF1=$2
MOTIF2=$3
MOTIF3=$4

ORGANISM=$(echo $INFILE | rev | cut -d "/" -f4 | rev)
STRAIN=$(echo $INFILE | rev | cut -d "/" -f3 | rev)
IN_FASTA=$(echo $INFILE | rev | cut -d "/" -f1 | rev)

CUR_PATH=$PWD
WORK_DIR=$TMPDIR/"$STRAIN""_motifsearch"
mkdir -p $WORK_DIR
cd $WORK_DIR

cp -r $CUR_PATH/$INFILE $IN_FASTA


/home/armita/git_repos/emr_repos/tools/pathogen/crinkler/find_crinkler.pl $IN_FASTA $MOTIF1 $MOTIF2 $MOTIF3

# cat findmotif_$MOTIF1.fa findmotif_$MOTIF3.fa | grep '>' | cut -f1 | sort | uniq -d > findmotif_"$MOTIF1"_"$MOTIF3".txt
# cat findmotif_$MOTIF2.fa findmotif_$MOTIF3.fa | grep '>' | cut -f1 | sort | uniq -d > findmotif_"$MOTIF2"_"$MOTIF3".txt  

rm $IN_FASTA

mkdir -p $CUR_PATH/analysis/crinkler/$ORGANISM/$STRAIN
cp * $CUR_PATH/analysis/crinkler/$ORGANISM/$STRAIN/.

rm -r $TMPDIR

exit

