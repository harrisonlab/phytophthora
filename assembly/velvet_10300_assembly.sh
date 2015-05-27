#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=4G


#	velvet_10300_assembly.sh
#
#	This script Runs genome assembly via velvet on trimmed sequence data for 
#	P. cactorum isolate 10300. Assembly is performed at a range of hash lengths.
#	Sequence data for five illumina runs are used. This requires velvet to be 
#	compiled with 5 categories. As such a local install of velvet is used which
#	has been compiled accordingly.



#---	Step 1		---
# 		Set Variables
#----------------------

# CurPath=$PWD
# VelvetPath=/home/armita/prog/velvet_1.2.08
# TrimPath=qc_dna/paired/P.cactorum/10300
# MatePath=qc_dna/mate-paired/P.cactorum/10300
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/velvet_1.2.08
# 
# MinHash=41
# MaxHash=81
# HashStep=2
# GenomeSz=70
# ExpCov=70
# MinCov=20
# Lib1InsLgth=300
# Lib1F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane4_300bp_R1_trim.fq.gz
# Lib1R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane4_300bp_R2_trim.fq.gz
# Lib2InsLgth=1000
# Lib2F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane5_1Kb_R1_trim.fq.gz
# Lib2R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane5_1Kb_R2_trim.fq.gz
# Lib3InsLgth=1000
# Lib3F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane3_1Kb_R1_trim.fq.gz  
# Lib3R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane3_1Kb_R2_trim.fq.gz  
# Lib4InsLgth=300
# Lib4F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane4_300bp_R1_trim.fq.gz
# Lib4R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane4_300bp_R2_trim.fq.gz
# Lib5InsLgth=5000
# Lib5F=$CurPath/$MatePath/F/Pcact10300_S2_L001_R1_001_trim_rev.fq.gz
# Lib5R=$CurPath/$MatePath/R/Pcact10300_S2_L001_R2_001_trim_rev.fq.gz

CurPath=$PWD
VelvetPath=/home/armita/prog/velvet_1.2.08
TrimPath=qc_dna/paired/P.cactorum/10300
MatePath=qc_dna/mate-paired/P.cactorum/10300
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/velvet_1.2.08

MinHash=41
MaxHash=81
HashStep=2
GenomeSz=70
ExpCov=70
MinCov=20
Lib1InsLgth=300
Lib1F=$CurPath/tmp_F_1.fq.gz
Lib1R=$CurPath/tmp_F_1.fq.gz
Lib2InsLgth=1000
Lib2F=$CurPath/tmp_F_2.fq.gz
Lib2R=$CurPath/tmp_F_2.fq.gz
Lib3InsLgth=1000
Lib3F=$CurPath/tmp_F_3.fq.gz  
Lib3R=$CurPath/tmp_F_3.fq.gz
Lib4InsLgth=300
Lib4F=$CurPath/tmp_F_4.fq.gz
Lib4R=$CurPath/tmp_F_4.fq.gz
Lib5InsLgth=5000
Lib5F=$CurPath/tmp_F_5.fq.gz
Lib5R=$CurPath/tmp_F_5.fq.gz

Strain=$(printf $TrimPath | rev | cut -f1 -d '/' | rev)
Organism=$(printf $TrimPath | rev | cut -f2 -d '/' | rev)
WorkDir=$TMPDIR/"$Strain"_assembly
AssemblyName="$Strain"_velvet



#---	Step 2		---
# 		Copy data onto
#		Worker node
#----------------------

mkdir -p $WorkDir
cd $WorkDir

cp $Lib1F Lib1_F.fq.gz
cp $Lib1R Lib1_R.fq.gz

cp $Lib2F Lib2_F.fq.gz
cp $Lib2R Lib2_R.fq.gz

cp $Lib3F Lib3_F.fq.gz
cp $Lib3R Lib3_R.fq.gz

cp $Lib4F Lib4_F.fq.gz
cp $Lib4R Lib4_R.fq.gz

cp $Lib5F Lib5_F.fq.gz
cp $Lib5R Lib5_R.fq.gz

#---	Step 3		---
# 		Assemble
#----------------------

$VelvetPath/velveth 10300_assembly 31,33,2 -fastq -shortPaired -separate Lib1_F.fq.gz Lib1_R.fq.gz -shortPaired2 -separate Lib2_F.fq.gz Lib2_R.fq.gz -shortPaired3 -separate  Lib3_F.fq.gz Lib3_R.fq.gz -shortPaired4 -separate  Lib4_F.fq.gz Lib4_R.fq.gz -shortPaired5 -separate Lib5_F.fq.gz Lib5_R.fq.gz
for Directory in $(ls -d $WorkDir/*/); do
	cd $Directory
	$VelvetPath/velvetg . -exp_cov $ExpCov -ins_length $Lib1InsLgth -ins_length2 $Lib2InsLgth -ins_length3 $Lib3InsLgth -ins_length4 $Lib4InsLgth -ins_length5 $Lib5InsLgth -shortMatePaired yes -min_contig_lgth 500
	process_contigs.pl -i $Directory/contigs.fa -o $Directory
	cd $WorkDir
done



#---	Step 4		---
# 		Cleanup
#----------------------

rm Lib1_F.fq.gz
rm Lib1_R.fq.gz

rm Lib2_F.fq.gz
rm Lib2_R.fq.gz

rm Lib3_F.fq.gz
rm Lib3_R.fq.gz

rm Lib4_F.fq.gz
rm Lib4_R.fq.gz

rm Lib5_F.fq.gz
rm Lib5_R.fq.gz

mkdir -p $CurPath/assembly/velvet/$Organism/$Strain/
cp -r $WorkDir/* $CurPath/assembly/velvet/$Organism/$Strain/.

rm -r $TMPDIR
