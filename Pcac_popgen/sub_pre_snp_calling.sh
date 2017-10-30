#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace
##############################################
# Prep mappings from Bowtie2 for SNP calling
### Remove multimapping reads, discordant reads. PCR and optical duplicates, and
### add read group and sample name to each mapped read (preferably, the shortest ID possible)
#INPUT:
# 1st argument: input SAM file with your mappings
# 2nd argument: sample name (Prefix) to be used to identify it in the future
#OUTPUT:
# Indexed BAM file with suffix "nodup_rg" to be fed into SNP calling with GATK.
#############################################

InputSam=$1
Prefix=$2
# Filename=$(basename "$InputSam")
OutDir=$(dirname "$InputSam")
Name="${Filename%.*}"

### Output folder
CurDir=$PWD
WorkDir="$TMPDIR"

### Prep
mkdir -p $WorkDir
cp $InputSam $WorkDir/.
cd $WorkDir

### Get rid of multimapping reads by filtering out on the XS:i: tag
grep -v "XS:i" *.sam > tmp.sam
# mv tmp.sam $Filename
samtools view --threads 16 -bS -o $Name.bam tmp.sam
samtools sort --threads 16 -o "$Name"_nomulti_sorted.bam $Name.bam
samtools index --threads 16 "$Name"_nomulti_sorted.bam

### Keep only reads with "paired reads" and "properly paired reads" flags.
samtools view -b -h -f 3 -o "$Name"_nomulti_proper.bam "$Name"_nomulti_sorted.bam
### Sort for downstream analyses
samtools sort "$Name"_nomulti_proper.bam "$Name"_nomulti_proper_sorted
samtools index "$Name"_nomulti_proper_sorted.bam

### Remove PCR and optical duplicates
ProgDir=/home/sobczm/bin/picard-tools-2.5.0
java -jar $ProgDir/picard.jar MarkDuplicates \
  INPUT="$Name"_nomulti_proper_sorted.bam \
  OUTPUT="$Name"_nomulti_proper_sorted_nodup.bam \
  METRICS_FILE="$Name"_nomulti_proper_sorted_nodup.txt \
  REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE MAX_RECORDS_IN_RAM=500000000 \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT
### Add group and sample name (Prefix)
java -jar $ProgDir/picard.jar AddOrReplaceReadGroups \
  INPUT="$Name"_nomulti_proper_sorted_nodup.bam \
  OUTPUT="$Name"_nomulti_proper_sorted_nodup_rg.bam \
  SORT_ORDER=coordinate CREATE_INDEX=true RGID=$Prefix  RGSM=$Prefix \
  RGPL=Illumina RGLB=library RGPU=barcode VALIDATION_STRINGENCY=LENIENT
  samtools index "$Name"_nomulti_proper_sorted_nodup_rg.bam

### Cleanup
mv "$Name"_nomulti_proper_sorted_nodup.txt $CurDir/$OutDir/.
mv "$Name"_nomulti_proper_sorted_nodup_rg.bam $CurDir/$OutDir/.
mv "$Name"_nomulti_proper_sorted_nodup_rg.bam.bai $CurDir/$OutDir/.
rm -rf $WorkDir
