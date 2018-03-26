#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=5.9G
#$ -l h=blacklace01.blacklace

# # Prefix="PcacP414"
# Prefix="Pi_SCRP370"
# Assembly=repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa
# F_reads=tmp_concat_dir/SCRP370_F_reads.fq.gz
# R_reads=tmp_concat_dir/SCRP370_R_reads.fq.gz
# # ControlBAM=PcacP414.bam
# OutDir=analysis/popgen/indel_calling/svaba

Usage="sub_svaba <prefix> <reference_assembly.fa> <directory_containing_BAM_Files> <output_directory>"
echo "$Usage"

Prefix=$1
Assembly=$2
BamDir=$3
OutDir=$4

CWD=$PWD

WorkDir="$TMPDIR"
# WorkDir=/tmp/svaba
mkdir -p $WorkDir

cp -r $Assembly $WorkDir/.
cp $BamDir/*.bam $WorkDir/.
cp $BamDir/*.bam.bai $WorkDir/.
# Cb=$(basename "$ControlBAM")
As=$(basename "$Assembly")

cd $WorkDir
# Output=${Fr%.*}.bam

# Align the data
bwa index $As
# bwa mem -t 16 $As $Fr $Rr | samtools view -S -b - > "$Prefix".bam

# samtools sort -@ 16 -o "$Prefix".sorted.bam "$Prefix".bam
# samtools index "$Prefix".sorted.bam

# samtools sort -@ 16 -o control.sorted.bam $Cb
# samtools index control.sorted.bam

# Control=normal.bam
# TargetRegion=22
# svaba -t "$Prefix".bam -n $Control -k $TargetRegion -G ref.fa -a test_id -p -4

# svaba run -t "$Prefix".bam -G $As -a "$Prefix"_sv -p 16
# svaba run -t "$Prefix".sorted.bam -n control.sorted.bam -G $As -a "$Prefix"_sv -p 24
BamFiles=$(ls *.bam | sed -E "s/$/ -t /g" | tr -d '\n' | sed -E "s/ -t $/ /g")
svaba run -t $BamFiles -G $As -a "$Prefix"_sv -p 16

rm $As*
rm $BamFiles
rm *.bam.bai
mkdir -p $OutDir
cp -r $WorkDir/* $CWD/$OutDir/.
# wget "https://data.broadinstitute.org/snowman/dbsnp_indel.vcf" ## get a DBSNP known indel file
# DBSNP=dbsnp_indel.vcf
# CORES=8 ## set any number of cores
# REF=/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta
# ## -a is any string you like, which gives the run a unique ID
# svaba run -t $TUM_BAM -n $NORM_BAM -p $CORES -D $DBSNP -a somatic_run -G $REF
