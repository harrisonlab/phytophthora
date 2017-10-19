#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=5.9G

#	This script Runs genome assembly via Spades on trimmed sequence data for
#	P. cactorum isolate 10300. Assembly is performed at a range of hash lengths.
#	Sequence data for five illumina runs are used.



#---	Step 1		---
# 		Set Variables
#----------------------


CurPath=$PWD
TrimPath=qc_dna/paired/P.cactorum/10300
MatePath=qc_dna/mate-paired/P.cactorum/10300


Lib1F=$(ls qc_dna/paired/P.cactorum/10300/F/*.fq.gz | grep -v 'Pcact10300_S2_L001_R1_001.fastq.gz' | head -n1 | tail -n1)
Lib1R=$(ls qc_dna/paired/P.cactorum/10300/R/*.fq.gz | grep -v 'Pcact10300_S2_L001_R2_001.fastq.gz'| head -n1 | tail -n1)
# Lib2InsLgth=1000
Lib2F=$(ls qc_dna/paired/P.cactorum/10300/F/*.fq.gz | grep -v 'Pcact10300_S2_L001_R1_001.fastq.gz'| head -n2 | tail -n1)
Lib2R=$(ls qc_dna/paired/P.cactorum/10300/R/*.fq.gz | grep -v 'Pcact10300_S2_L001_R2_001.fastq.gz'| head -n2 | tail -n1)
# Lib3InsLgth=1000
Lib3F=$(ls qc_dna/paired/P.cactorum/10300/F/*.fq.gz | grep -v 'Pcact10300_S2_L001_R1_001.fastq.gz' | head -n3 | tail -n1)
Lib3R=$(ls qc_dna/paired/P.cactorum/10300/R/*.fq.gz | grep -v 'Pcact10300_S2_L001_R2_001.fastq.gz' | head -n3 | tail -n1)
# Lib4InsLgth=300
Lib4F=$(ls qc_dna/paired/P.cactorum/10300/F/*.fq.gz | grep -v 'Pcact10300_S2_L001_R1_001.fastq.gz' | head -n4 | tail -n1)
Lib4R=$(ls qc_dna/paired/P.cactorum/10300/R/*.fq.gz | grep -v 'Pcact10300_S2_L001_R2_001.fastq.gz' | head -n4 | tail -n1)
# Lib5InsLgth=5000
# Lib5F=$(ls qc_dna/mate-paired/P.cactorum/10300/F/*.fq.gz | grep -v 'rev' | head -n1 | tail -n1)
# Lib5R=$(ls qc_dna/mate-paired/P.cactorum/10300/R/*.fq.gz | grep -v 'rev' | head -n1 | tail -n1)
Lib5F=$(ls qc_dna/mate-paired/P.cactorum/10300/nextclip_no_rev/P.cactorum_10300_trim_clip_D_R1.fastq.gz)
Lib5R=$(ls qc_dna/mate-paired/P.cactorum/10300/nextclip_no_rev/P.cactorum_10300_trim_clip_D_R2.fastq.gz)

Strain=10300
Organism=P.cactorum
WorkDir=$TMPDIR/"$Strain"_assembly
OutDir=$CurPath/assembly/spades2/$Organism/$Strain
Cutoff="auto"
AssemblyName="$Strain"_spades



#---	Step 2		---
# 		Copy data onto
#		Worker node
#----------------------

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$Lib1F Lib1_1.fq.gz
cp $CurPath/$Lib1R Lib1_2.fq.gz

cp $CurPath/$Lib2F Lib2_1.fq.gz
cp $CurPath/$Lib2R Lib2_2.fq.gz

cp $CurPath/$Lib3F Lib3_1.fq.gz
cp $CurPath/$Lib3R Lib3_2.fq.gz

cp $CurPath/$Lib4F Lib4_1.fq.gz
cp $CurPath/$Lib4R Lib4_2.fq.gz

cp $CurPath/$Lib5F Lib5_1.fq.gz
cp $CurPath/$Lib5R Lib5_2.fq.gz

#---	Step 3		---
# 		Assemble
#----------------------

# abyss-pe k=64 -np 16 -j 16 name=P.cac_10300 lib='pe1 pe2 pe3 pe4 mp5' pe1='Lib1_1.fq.gz Lib1_2.fq.gz' pe2='Lib2_1.fq.gz Lib2_2.fq.gz' pe3='Lib3_1.fq.gz Lib3_2.fq.gz' pe4='Lib4_1.fq.gz Lib4_2.fq.gz' mp5='Lib5_1.fq.gz Lib5_2.fq.gz'

spades.py \
    -k 21,33,45 \
    -m 200 \
    --phred-offset 33 \
    --pe1-1 Lib1_1.fq.gz \
    --pe1-2 Lib1_2.fq.gz \
    --pe2-1 Lib2_1.fq.gz \
    --pe2-2 Lib2_2.fq.gz \
    --pe3-1 Lib3_1.fq.gz \
    --pe3-2 Lib3_2.fq.gz \
    --pe4-1 Lib4_1.fq.gz \
    --pe4-2 Lib4_2.fq.gz \
    --mp1-1 Lib5_1.fq.gz \
    --mp1-2 Lib5_2.fq.gz \
    -t 16  \
    -o $WorkDir/. \
    --cov-cutoff "$Cutoff"

    # spades.py \
    #     -k 21,33,45,57 \
    #     -m 200 \
    #     --phred-offset 33 \
    #     --pe1-1 Lib1_1.fq.gz \
    #     --pe1-2 Lib1_2.fq.gz \
    #     --pe2-1 Lib2_1.fq.gz \
    #     --pe2-2 Lib2_2.fq.gz \
    #     --pe3-1 Lib3_1.fq.gz \
    #     --pe3-2 Lib3_2.fq.gz \
    #     --pe4-1 Lib4_1.fq.gz \
    #     --pe4-2 Lib4_2.fq.gz \
    #     --mp1-fr \
    #     --mp1-1 Lib5_1.fq.gz \
    #     --mp1-2 Lib5_2.fq.gz \
    #     -t 16  \
    #     -o $WorkDir/. \
    #     --cov-cutoff "$Cutoff"

echo "Filtering contigs smaller than 500bp"
mkdir -p $WorkDir/filtered_contigs
FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
$FilterDir/filter_abyss_contigs.py $WorkDir/scaffolds.fasta 500 > $WorkDir/filtered_contigs/contigs_min_500bp.fasta



#---	Step 4		---
# 		Cleanup
#----------------------

rm Lib1_1.fq.gz
rm Lib1_2.fq.gz

rm Lib2_1.fq.gz
rm Lib2_2.fq.gz

rm Lib3_1.fq.gz
rm Lib3_2.fq.gz

rm Lib4_1.fq.gz
rm Lib4_2.fq.gz

rm Lib5_1.fq.gz
rm Lib5_2.fq.gz

mkdir -p $OutDir
cp -r $WorkDir/* $OutDir/.
