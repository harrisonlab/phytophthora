#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l virtual_free=4G

set -u
set -e

# P.10300 Dipspades assembly

#---	Step 1		---
# 		Set Variables
#----------------------

Strain=10300
Organism=P.cactorum

CurPath=/home/groups/harrisonlab/project_files/idris
TrimPath=qc_dna/paired/P.cactorum/10300
MatePath=qc_dna/mate-paired/P.cactorum/10300

# AssemblyName="$Strain"_dip-spades
AssemblyName="$Strain"_dip-spades_no_rev
WorkDir=$TMPDIR/"$Strain"_dip-spades
OutDir=$CurPath/assembly/dip-spades/$Organism/$Strain


Lib1F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane4_300bp_R1_trim.fq.gz
Lib1R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane4_300bp_R2_trim.fq.gz

Lib2F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane5_1Kb_R1_trim.fq.gz
Lib2R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane5_1Kb_R2_trim.fq.gz

Lib3F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane3_1Kb_R1_trim.fq.gz  
Lib3R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane3_1Kb_R2_trim.fq.gz  

Lib4F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane4_300bp_R1_trim.fq.gz
Lib4R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane4_300bp_R2_trim.fq.gz

# Lib5F=$CurPath/$MatePath/nextclip/P.cactorum_10300_trim_rev_clip_D_R1.fastq.gz
# Lib5R=$CurPath/$MatePath/nextclip/P.cactorum_10300_trim_rev_clip_D_R2.fastq.gz
Lib5F=$CurPath/$MatePath/nextclip_no_rev/P.cactorum_10300_trim_clip_D_R1.fastq.gz
Lib5R=$CurPath/$MatePath/nextclip_no_rev/P.cactorum_10300_trim_clip_D_R2.fastq.gz

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

echo "Running dip-spades" 2>&1 |  tee -a $WorkDir/"$Organism"_"$Strain"_dipspades.log
dipspades.py --pe1-1 Lib1_1.fq.gz --pe1-2 Lib1_2.fq.gz \
--pe2-1 Lib2_1.fq.gz --pe2-2 Lib2_2.fq.gz \
--pe3-1 Lib3_1.fq.gz --pe3-2 Lib3_2.fq.gz \
--pe4-1 Lib4_1.fq.gz --pe4-2 Lib4_2.fq.gz \
--mp1-1 Lib5_1.fq.gz --mp1-2 Lib5_2.fq.gz \
-t 16 -m 96 -k 21,33,55,77 --careful -o $AssemblyName 2>&1 |  tee -a $WorkDir/"$Organism"_"$Strain"_dipspades.log

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
cp -r $WorkDir/$AssemblyName/* $OutDir/.
cp $WorkDir/"$Organism"_"$Strain"_dipspades.log $OutDir/.
rm -r $WorkDir


