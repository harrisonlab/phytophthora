
# 1. Alignment of Pcac raw reads vs the 414 genome

Alignment of reads from a single run:

```bash
  Reference=$(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w '414')
  for StrainPath in $(ls -d qc_dna/paired/P.cactorum/* | grep -v '10300' | grep -v '404' | grep -v '414' | grep -v '411'); do
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

Alignment of reads from multiple sequencing runs:

```bash
  Reference=$(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w '414')
  for StrainPath in $(ls -d qc_dna/paired/P.cactorum/* | grep -v -w -e '404' -e 414); do
    echo $StrainPath
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    F1_Read=$(ls $StrainPath/F/*.fq.gz | head -n1);
    R1_Read=$(ls $StrainPath/R/*.fq.gz | head -n1);
    F2_Read=$(ls $StrainPath/F/*.fq.gz | tail -n1);
    R2_Read=$(ls $StrainPath/R/*.fq.gz | tail -n1);
    echo $F1_Read
    echo $R1_Read
    echo $F2_Read
    echo $R2_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_414
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie_2lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $OutDir
  done
  for StrainPath in $(ls -d qc_dna/paired/P.cactorum/* | grep -w -e '10300'); do
    echo $StrainPath
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    F1_Read=$(ls $StrainPath/F/*300bp*.fq.gz | head -n1);
    R1_Read=$(ls $StrainPath/R/*300bp*.fq.gz | head -n1);
    F2_Read=$(ls $StrainPath/F/*300bp*.fq.gz | tail -n1);
    R2_Read=$(ls $StrainPath/R/*300bp*.fq.gz | tail -n1);
    echo $F1_Read
    echo $R1_Read
    echo $F2_Read
    echo $R2_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/300bp_vs_414
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie_2lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $OutDir
    F1_Read=$(ls $StrainPath/F/*1Kb*.fq.gz | head -n1);
    R1_Read=$(ls $StrainPath/R/*1Kb*.fq.gz | head -n1);
    F2_Read=$(ls $StrainPath/F/*1Kb*.fq.gz | tail -n1);
    R2_Read=$(ls $StrainPath/R/*1Kb*.fq.gz | tail -n1);
    echo $F1_Read
    echo $R1_Read
    echo $F2_Read
    echo $R2_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/1Kb_vs_414
    qsub $ProgDir/bowtie/sub_bowtie_2lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $OutDir
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
