#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h=blacklace06.blacklace
###

#Input: Just the name for Output results files. The script needs to be run from the directory containing all the Output files from sub_bwa_mem.sh for the samples of interest (but not any others!!)

#Output: Genotype (*.gt) and VCF (*.vcf) files with structural variants.

Output=$1

#First, obtain command line arguments with: 1) all the alignments
# 2) discordant alignments, 3) split alignments to be used as input
# by Lumpy

#1)
Bams=()
for b in *rg.bam; do Bams+=("$b"); done;
Bam=$(IFS=, ; echo "${Bams[*]}")

#2)
Discs=()
for b in *discordants.bam; do Discs+=("$b"); done;
discordant=$(IFS=, ; echo "${Discs[*]}")

#3)
Splits=()
for b in *splitters.bam; do Splits+=("$b"); done;
Splitters=$(IFS=, ; echo "${Splits[*]}")

Out=${Output}.gt

#Detect structural variants
lumpyexpress -B $Bam -S $Splitters -D $discordant -o $Out

#Call genotypes
# svtyper=/home/sobczm/bin/svtyper
svtyper -B $Bam -S $Splitters -i $Out >${Output}.vcf
