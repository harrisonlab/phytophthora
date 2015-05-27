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

CurPath=$PWD
ProgDir=/home/armita/prog/velvet_1.2.08
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
Lib1F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane4_300bp_R1_trim.fq.gz
Lib1R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane4_300bp_R2_trim.fq.gz
Lib2InsLgth=1000
Lib2F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane5_1Kb_R1_trim.fq.gz
Lib2R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane5_1Kb_R2_trim.fq.gz
Lib3InsLgth=1000
Lib3F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane3_1Kb_R1_trim.fq.gz  
Lib3R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane3_1Kb_R2_trim.fq.gz  
Lib4InsLgth=300
Lib4F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane4_300bp_R1_trim.fq.gz
Lib4R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane4_300bp_R2_trim.fq.gz
Lib5InsLgth=5000
Lib5F=$CurPath/$MatePath/F/Pcact10300_S2_L001_R1_001_trim_rev.fq.gz
Lib5R=$CurPath/$MatePath/R/Pcact10300_S2_L001_R2_001_trim_rev.fq.gz


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

$ProgDir/velveth 10300_assembly 31,41,2 -fastq -shortPaired -separate Lib1_F.fq.gz Lib1_R.fq.gz -shortPaired2 -separate Lib2_F.fq.gz Lib2_R.fq.gz -shortPaired3 -separate  Lib3_F.fq.gz Lib3_R.fq.gz -shortPaired4 -separate  Lib4_F.fq.gz Lib4_R.fq.gz -shortPaired5 -separate Lib5_F.fq.gz Lib5_R.fq.gz
for Directory in $(ls -d $WorkDir/*/); do
	cd $Directory
	$ProgDir/velvetg . -exp_cov $ExpCov -cov_cutoff $CovCutoff -ins_length $Lib1InsLgth -ins_lgth2 $Lib2InsLgth -ins_lgth3 $Lib3InsLgth  -ins_lgth4 $Lib4InsLgth -ins_lgth5 $Lib5InsLgth -shortMatePaired yes -min_contig_lgth 500
	process_contigs.pl -i $Directory/contigs.fa -o $Directory
	cd $WorkDir
done

mkdir -p $CurPath/assembly/velvet/$Organism/$Strain/
cp -r $WorkDir/* $CurPath/assembly/velvet/$Organism/$Strain/.

#---	Step 4		---
# 		Cleanup
#----------------------

rm -r $TMPDIR
