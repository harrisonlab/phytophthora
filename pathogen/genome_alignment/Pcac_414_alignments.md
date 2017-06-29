
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

The multiple alignments of 10300 reads were aligned using the following commands:

```bash
  OutDir=analysis/genome_alignment/bowtie/P.cactorum/10300/vs_414
  mkdir -p $OutDir
  samtools merge -f $OutDir/414_contigs_unmasked.fa_aligned_sorted.bam analysis/genome_alignment/bowtie/P.cactorum/10300/*_vs_414/414_contigs_unmasked.fa_aligned_sorted.bam
```


Alignment of reads from P.idaei:

```bash
  Reference=$(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w '414')
  for StrainPath in $(ls -d qc_dna/paired/P.idaei/*); do
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_414
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
  done
```
<!--
Bowtie is not appropriate for aligning Pacbio reads - too much memory overhead.

Alignment of Pacbio reads vs the P414 genome. PacBio reads maxed out the RAM of
the node and therefore the wrapper script was modified to allow submission to
a node with greater RAM.

```bash
  Reference=$(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w '414')
  for Pacbio_Reads in $(ls -d raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq); do
    Organism=$(echo $Pacbio_Reads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Pacbio_Reads | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/genome_alignment/bowtie/$Organism/"$Strain"_pacbio/vs_414
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub -h $ProgDir/bowtie/sub_bowtie_unpaired.sh $Reference $Pacbio_Reads $OutDir
    # The job ID was identified as:
    JobID=$(qstat | grep 'sub_bowtie' | tail -n1 | cut -f1 -d ' ')
    qalter -N "Pacbio_bowtie" $JobID
    qalter -l h=blacklace01.blacklace $JobID
    qalter -l virtual_free=5G $JobID
    qalter -h U $JobID
  done
``` -->
