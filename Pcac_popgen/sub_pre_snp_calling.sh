#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace
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
Filename=$(basename "$InputSam")
OutDir=$(dirname "$InputSam")
Name="${Filename%.*}"

### Output folder
CurDir=$PWD
WorkDir="$TMPDIR"

### Prep
mkdir -p $WorkDir
cp $InputSam $WorkDir/.
cd $WorkDir


### If the provided file is a .bam file then convert it into sam format.

Extension=$(echo $Filename | rev | cut -f1 -d '.' | rev)
if [[ $Extension == "bam" ]]; then
  echo ".bam file extension found"
  echo "converting to sam"
  samtools view --threads 8 -o $Name.sam $CurDir/$InputSam
  echo "Reheadering the bam file using the original sam headers"
  samtools view -H $CurDir/$InputSam > header.sam
  cat header.sam $Name.sam > ${Name}_reheader.sam
  Filename=${Name}_reheader.sam
fi

### Get rid of multimapping reads by filtering out on the XS:i: tag
echo "Removing reads with the XS:i flag (multimapping reads)"
grep -v "XS:i" $Filename > tmp.sam
# mv tmp.sam $Filename
echo "Converting sam to bam format"
samtools view --threads 8 -bS -o $Name.bam tmp.sam

echo "Sorting and indexing the bam file"
samtools sort --threads 8 -o "$Name"_nomulti_sorted.bam $Name.bam
samtools index -@ 8 "$Name"_nomulti_sorted.bam "$Name"_nomulti_sorted.bam.index

### Keep only reads with "paired reads" and "properly paired reads" flags.
echo "Filtering reads to retain only concordant paired reads"
samtools view -b -h -f 3 -o "$Name"_nomulti_proper.bam "$Name"_nomulti_sorted.bam
### Sort for downstream analyses
"Sorting and indexing concordant, paired reads bam file"
samtools sort --threads 8 -o "$Name"_nomulti_proper_sorted.bam "$Name"_nomulti_proper.bam
samtools index -@ 8 "$Name"_nomulti_proper_sorted.bam "$Name"_nomulti_proper_sorted.bam.index

### Remove PCR and optical duplicates
echo "Using picard tools to remove duplicate reads"
ProgDir=/home/sobczm/bin/picard-tools-2.5.0
java -jar $ProgDir/picard.jar MarkDuplicates \
  INPUT="$Name"_nomulti_proper_sorted.bam \
  OUTPUT="$Name"_nomulti_proper_sorted_nodup.bam \
  METRICS_FILE="$Name"_nomulti_proper_sorted_nodup.txt \
  REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE MAX_RECORDS_IN_RAM=500000000 \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT
### Add group and sample name (Prefix)
echo "Adding the given sample information to the BAM file: $Prefix"
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
