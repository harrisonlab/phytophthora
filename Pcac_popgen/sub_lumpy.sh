#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h=blacklace06.blacklace
###

#Input: Just the name for Output results files. The script needs to be run from the directory containing all the Output files from sub_bwa_mem.sh for the samples of interest (but not any others!!)

#Output: Genotype (*.gt) and VCF (*.vcf) files with structural variants.

Prefix=$1
Strain=$2

#First, obtain command line arguments with: 1) all the alignments
# 2) discordant alignments, 3) split alignments to be used as input
# by Lumpy

#1)
Bams=()
for b in $(ls *_sorted.bam | grep -v -e 'splitters' -e 'discordant' | grep "${Strain}_"); do
  Bams+=("$b");
done;
Bam=$(IFS=, ; echo "${Bams[*]}")

#2)
Discs=()
for b in $(ls *discordants_sorted.bam | grep "${Strain}_"); do
  Discs+=("$b");
done;
Discordant=$(IFS=, ; echo "${Discs[*]}")

#3)
Splits=()
for b in $(ls *splitters_sorted.bam | grep "${Strain}_"); do
  Splits+=("$b");
done;
Splitters=$(IFS=, ; echo "${Splits[*]}")

# Out="$Prefix".gt

echo "Strain: $Strain"
echo "Bams: $Bams"
echo "Discordant: $Discordant"
echo "Splitters: $Splitters"

#Detect structural variants
lumpyexpress -B $Bam -S $Splitters -D $Discordant -o "$Prefix"_express.vcf

#Call genotypes
# svtyper=/home/sobczm/bin/svtyper
svtyper -B $Bam -S $Splitters -i "$Prefix"_express.vcf > "$Prefix".vcf
