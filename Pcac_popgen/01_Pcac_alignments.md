
# 1. Alignment of Pcac raw reads vs the 414 genome

Alignment of reads from a single run:

```bash
  Reference=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa)
  for StrainPath in $(ls -d qc_dna/paired/P.*/* | grep -v -e '10300' -e '404' -e '414' -e '415' -e '416' -e 'PC13_15' -e '2003_3' | grep -e 'P.cactorum' -e 'P.idaei'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
    R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_414
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
  done
```

Alignment of reads from multiple sequencing runs:

For isolates with two runs of data:

```bash
  Reference=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa)
  for StrainPath in $(ls -d qc_dna/paired/P*/* | grep -e 'P.cactorum' -e 'P.idaei' | grep -w -e '2003_3' -e '415' -e '416' -e 'PC13_15'); do
    echo $StrainPath
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    F1_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1);
    R1_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1);
    F2_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n2 | tail -n1);
    R2_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n2 | tail -n1);
    echo $F1_Read
    echo $R1_Read
    echo $F2_Read
    echo $R2_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_414
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie_2lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $OutDir
  done
```

for isolates with three runs of data:

```bash
  Reference=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa)
  for StrainPath in $(ls -d qc_dna/paired/P*/* | grep -e 'P.cactorum' -e 'P.idaei' | grep -w -e '404' -e '414'); do
      echo $StrainPath
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
      Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
      Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
      echo $Strain
      echo $Organism
      F1_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1);
      R1_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1);
      F2_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n2 | tail -n1);
      R2_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n2 | tail -n1);
      F3_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n3 | tail -n1);
      R3_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n3 | tail -n1);
      echo $F1_Read
      echo $R1_Read
      echo $F2_Read
      echo $R2_Read
      echo $F3_Read
      echo $R3_Read
      OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_414
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
      qsub $ProgDir/bowtie/sub_bowtie_3lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $F3_Read $R3_Read $OutDir
    done
```


# 2. Pre SNP calling cleanup


## 2.1 Rename input mapping files in each folder by prefixing with the strain ID

```bash
  for File in $(ls analysis/genome_alignment/bowtie/*/*/vs_414/414_v2_contigs_unmasked.fa_aligned.sam | grep -v '10300'); do
    Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
    echo $Strain
    echo $Organism
    OutDir=analysis/popgen/$Organism/$Strain
    CurDir=$PWD
    mkdir -p $OutDir
    cd $OutDir
    cp -s $CurDir/$File "$Strain"_vs_414_v2_aligned.sam
    cd $CurDir
  done
```

## 2.2 Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used:
qsub $ProgDir/sub_pre_snp_calling.sh <INPUT SAM FILE> <SAMPLE_ID>

```bash
  for Sam in $(ls $PWD/analysis/popgen/*/*/*_vs_414_v2_aligned.sam | grep -v '10300' | tail -n+2); do
    Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Sam | rev | cut -f3 -d '/' | rev)
    CurDir=$PWD
    OutDir=$(dirname $Sam)
    cd $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
    qsub $ProgDir/sub_pre_snp_calling.sh $Sam $Strain
    cd $CurDir
  done
```

<!-- ##Copy outputs from cleanup to alignment folder

```bash
for Strain in A4 Bc1 Bc23 Nov27 Nov5 Nov71 Nov77 Nov9 ONT3 SCRP245_v2 SCRP249 SCRP324 SCRP333
do
    Bam="$Strain"_95m_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.bam
    rgBam="$Strain"_95m_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam
    Bai="$Strain"_95m_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam.bai
    Txt="$Strain"_95m_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.txt
    Directory=analysis/genome_alignment/bowtie/*/$Strain/vs_Bc16_unmasked_max1200/
    mv $Bam $Directory
    mv $rgBam $Directory
    mv $Bai $Directory
    mv $Txt $Directory
done
Strain=Bc16
Bam="$Strain"_95m_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.bam
rgBam="$Strain"_95m_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam
Bai="$Strain"_95m_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam.bai
Txt="$Strain"_95m_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.txt
Directory=analysis/genome_alignment/bowtie/*/$Strain/vs_Bc16_unmasked_max1200_SNP/
mv $Bam $Directory
mv $rgBam $Directory
mv $Bai $Directory
mv $Txt $Directory
``` -->

# 3. Run SNP calling

#Runs a SNP calling script from Maria in order to be able to draw up a phylogeny
To change in each analysis:

```bash
input=/home/groups/harrisonlab/project_files/phytophthora_fragariae/analysis/genome_alignment/bowtie
reference=repeat_masked/P.fragariae/Bc16/filtered_contigs_repmask/95m_contigs_unmasked.fa

filename=$(basename "$reference")
output="${filename%.*}.dict"
```

##Prepare genome reference indexes required by GATK

```bash
Reference=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa)
OutDir=$(dirname $Reference)
mkdir -p $OutDir
ProgDir=/home/sobczm/bin/picard-tools-2.5.0
java -jar $ProgDir/picard.jar CreateSequenceDictionary R=$Reference O=$OutDir/414_v2_contigs_unmasked.dict
samtools faidx $Reference
```

###Copy index file to same folder as BAM alignments

```bash
Reference=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa)
for AlignDir in $(ls -d analysis/popgen/P.*/*/); do
    Index="$Reference".dict
    Directory=analysis/genome_alignment/bowtie/*/$Strain/vs_Bc16_unmasked_max1200/
    cp $Index $AlignDir/.
done
```

Move to the directory where the output of SNP calling should be placed. Then
Start SNP calling with GATK.
The submission script required need to be custom-prepared for each analysis,
depending on what samples are being analysed. See inside the submission script
below:

```bash
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling
mkdir -p $OutDir
cd $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
qsub $ProgDir/sub_SNP_calling_multithreaded.sh
cd $CurDir
```

```bash

ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
qsub $ProgDir/sub_SNP_calling_multithreaded.sh
```
