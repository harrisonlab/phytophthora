#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=4G

set -u
set -e
set -o pipefail

#	abyss_10300_assembly.sh
#
#	This script Runs genome assembly via Abyss on trimmed sequence data for 
#	P. cactorum isolate 10300. Assembly is performed at a range of hash lengths.
#	Sequence data for five illumina runs are used.



#---	Step 1		---
# 		Set Variables
#----------------------


CurPath=$PWD
TrimPath=qc_dna/paired/P.cactorum/10300
MatePath=qc_dna/mate-paired/P.cactorum/10300

# MinHash=31
# MaxHash=81
# HashStep=10
# GenomeSz=70
# ExpCov=70
# MinCov=20
# Lib1InsLgth=300
# Lib1F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane4_300bp_R1_trim.fq.gz
# Lib1R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane4_300bp_R2_trim.fq.gz
# # Lib2InsLgth=1000
# Lib2F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane5_1Kb_R1_trim.fq.gz
# Lib2R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane5_1Kb_R2_trim.fq.gz
# # Lib3InsLgth=1000
# Lib3F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane3_1Kb_R1_trim.fq.gz  
# Lib3R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane3_1Kb_R2_trim.fq.gz  
# # Lib4InsLgth=300
# Lib4F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane4_300bp_R1_trim.fq.gz
# Lib4R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane4_300bp_R2_trim.fq.gz
# # Lib5InsLgth=5000
# Lib5F=$CurPath/$MatePath/F/Pcact10300_S2_L001_R1_001_trim_rev.fq.gz
# Lib5R=$CurPath/$MatePath/R/Pcact10300_S2_L001_R2_001_trim_rev.fq.gz


# MinHash=31
# MaxHash=31
# HashStep=2
# GenomeSz=70
# ExpCov=70
# MinCov=20
# Lib1InsLgth=300
Lib1F=$CurPath/tmp_F_1.fq.gz
Lib1R=$CurPath/tmp_R_1.fq.gz
Lib2InsLgth=1000
Lib2F=$CurPath/tmp_F_2.fq.gz
Lib2R=$CurPath/tmp_R_2.fq.gz
Lib3InsLgth=1000
Lib3F=$CurPath/tmp_F_3.fq.gz  
Lib3R=$CurPath/tmp_R_3.fq.gz
Lib4InsLgth=300
Lib4F=$CurPath/tmp_F_4.fq.gz
Lib4R=$CurPath/tmp_R_4.fq.gz
Lib5InsLgth=5000
Lib5F=$CurPath/tmp_F_5.fq.gz
Lib5R=$CurPath/tmp_R_5.fq.gz


# Strain=$(printf $TrimPath | rev | cut -f1 -d '/' | rev)
# Organism=$(printf $TrimPath | rev | cut -f2 -d '/' | rev)
Strain=10300
Organism=P.cactorum
WorkDir=$TMPDIR/"$Strain"_assembly
OutDir=$CurPath/assembly/abyss/$Organism/$Strain
AssemblyName="$Strain"_abyss



#---	Step 2		---
# 		Copy data onto
#		Worker node
#----------------------

mkdir -p $WorkDir
cd $WorkDir

cp $Lib1F Lib1_1.fq.gz
cp $Lib1R Lib1_2.fq.gz

cp $Lib2F Lib2_1.fq.gz
cp $Lib2R Lib2_2.fq.gz

cp $Lib3F Lib3_1.fq.gz
cp $Lib3R Lib3_2.fq.gz

cp $Lib4F Lib4_1.fq.gz
cp $Lib4R Lib4_2.fq.gz

cp $Lib5F Lib5_1.fq.gz
cp $Lib5R Lib5_2.fq.gz

#---	Step 3		---
# 		Assemble
#----------------------

# $VelvetPath/velveth 10300_assembly $MinHash,$MaxHash,$HashStep -fastq -shortPaired -separate Lib1_F.fq.gz Lib1_R.fq.gz -shortPaired2 -separate Lib2_F.fq.gz Lib2_R.fq.gz -shortPaired3 -separate  Lib3_F.fq.gz Lib3_R.fq.gz -shortPaired4 -separate  Lib4_F.fq.gz Lib4_R.fq.gz -shortPaired5 -separate Lib5_F.fq.gz Lib5_R.fq.gz
# for Directory in $(ls -d $WorkDir/*/); do
# 	cd $Directory
# 	$VelvetPath/velvetg . -exp_cov $ExpCov -ins_length $Lib1InsLgth -ins_length2 $Lib2InsLgth -ins_length3 $Lib3InsLgth -ins_length4 $Lib4InsLgth -ins_length5 $Lib5InsLgth -shortMatePaired yes -min_contig_lgth 500
# 	process_contigs.pl -i $Directory/contigs.fa -o $Directory
# 	cd $WorkDir
# done

abyss-pe k=64 -np 16 -j 16 name=P.cac_10300 lib='pe1 pe2 pe3 pe4 mp5' pe1='Lib1_1.fq.gz Lib1_2.fq.gz' pe2='Lib2_1.fq.gz Lib2_2.fq.gz' pe3='Lib3_1.fq.gz Lib3_2.fq.gz' pe4='Lib4_1.fq.gz Lib4_2.fq.gz' mp5='Lib5_1.fq.gz Lib5_2.fq.gz'



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

rm -r $TMPDIR
