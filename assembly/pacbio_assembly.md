
## Data extraction


for P.cactorum data:
```bash
  cd /home/groups/harrisonlab/project_files/idris
  RawDatDir=/home/harrir/projects/pacbio_test/p_cact
  mkdir -p raw_dna/pacbio/P.cactorum/414
  cp -r $RawDatDir/A07_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/B07_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/G06_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/H06_1 raw_dna/pacbio/P.cactorum/414/.
  OutDir=raw_dna/pacbio/P.cactorum/414/extracted
  mkdir -p $OutDir
  cat raw_dna/pacbio/P.cactorum/414/*/Analysis_Results/*.subreads.fastq | gzip -cf > $OutDir/concatenated_pacbio.fastq.gz
  #
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/pacbio/Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage.tar.gz
  mkdir -p raw_dna/pacbio/P.cactorum/414
  cp -r $RawDat raw_dna/pacbio/P.cactorum/414/.
  cd raw_dna/pacbio/P.cactorum/414
  tar -zxvf Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage.tar.gz
  # Data for 414 was contained in E02_1, F02_1 and G02_1
  cat \
  Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage/E02_1/Analysis_Results/*.subreads.fastq \
  Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage/F02_1/Analysis_Results/*.subreads.fastq \
  Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage/G02_1/Analysis_Results/*.subreads.fastq \
  | gzip -cf > extracted/concatenated_pacbio_extra_coverage.fastq.gz
  rm Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage.tar.gz
  rm -r Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage
```

Pacbio coverage was determined using:


Data quality was visualised once again following trimming:

```bash
for RawData in $(ls raw_dna/pacbio/P.cactorum/414/extracted/*q.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
GenomeSz=65
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

86.44 coverage was obtained from pacbio sequencing
plus 76.49 illumina sequencing

<!--
for P. fragariae data (commands for tom to run)

```bash
  cd /home/groups/harrisonlab/project_files/phytophthora_fragariae
  RawDatDir=/home/harrir/projects/pacbio_test/p_frag
  mkdir -p raw_dna/pacbio/P.fragariae/Bc16
  cp -r $RawDatDir/C07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  cp -r $RawDatDir/D07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  cp -r $RawDatDir/E07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  cp -r $RawDatDir/F07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  OutDir=raw_dna/pacbio/P.fragariae/Bc16/extracted
  mkdir -p $OutDir
  cat raw_dna/pacbio/P.fragariae/Bc16/*/Analysis_Results/*.subreads.fastq > $OutDir/concatenated_pacbio.fastq
``` -->

<!-- Sequencing depth and possible contamination was identified using kmer counting
kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

```bash
  for Reads in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
    echo $Reads
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/kmc_kmer_counting.sh $Reads
  done
``` -->

## Assembly


### Canu assembly

<!--
```bash
  Reads=$(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq.gz)
  Run1=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq.gz)
  Run2=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz)
  Reads=raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_both.fastq.gz
  cat $Run1 $Run2 > $Reads
  GenomeSz="65m"
  Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir=assembly/canu/$Organism/"$Strain"_3
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
```

```bash
Run1=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq.gz)
GenomeSz="75m"
Strain=$(echo $Run1 | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Run1 | rev | cut -f4 -d '/' | rev)
Prefix="$Strain"_canu
OutDir=assembly/canu/$Organism/"$Strain"_run1only
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/submit_canu.sh $Run1 $GenomeSz $Prefix $OutDir

  Run2=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz)
  GenomeSz="75m"
  Strain=$(echo $Run2 | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Run2 | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir=assembly/canu/$Organism/"$Strain"_run2only
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu.sh $Run2 $GenomeSz $Prefix $OutDir
```
-->

```bash
  Run1=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq.gz)
  Run2=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz)
  GenomeSz="75m"
  Strain=$(echo $Run1 | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Run1 | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir=assembly/canu/$Organism/"$Strain"_modified_script
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu_2lib.sh $Run1 $Run2 $GenomeSz $Prefix $OutDir
```


```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/canu/*/*/*.contigs.fasta | grep -v 'old' | grep -w '414'); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

Assemblies were polished using Pilon

```bash
  for Assembly in $(ls assembly/canu/*/*/*.contigs.fasta | grep -v 'old' | grep -w '414'); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
    TrimF2_Read=$(ls $IlluminaDir/F/414_run2_F_trim.fq.gz);
    TrimR2_Read=$(ls $IlluminaDir/R/414_run2_R_trim.fq.gz);
    TrimF3_Read=$(ls $IlluminaDir/F/414_170210_F_trim.fq.gz);
    TrimR3_Read=$(ls $IlluminaDir/R/414_170210_R_trim.fq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    echo $TrimF3_Read
    echo $TrimR3_Read
    OutDir=assembly/canu/$Organism/$Strain/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon_3_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
  done
```

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/canu/*/*/polished/pilon.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/canu/$Organism/$Strain/polished
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

### Spades Assembly

For P. cactorum

```bash
  # for PacBioDat in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
for PacBioDat in $(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz); do
echo $StrainPath
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
TrimF2_Read=$(ls $IlluminaDir/F/414_run2_F_trim.fq.gz);
TrimR2_Read=$(ls $IlluminaDir/R/414_run2_R_trim.fq.gz);
TrimF3_Read=$(ls $IlluminaDir/F/414_170210_F_trim.fq.gz);
TrimR3_Read=$(ls $IlluminaDir/R/414_170210_R_trim.fq.gz);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
# OutDir=assembly/spades_pacbio/$Organism/$Strain
OutDir=assembly/spades_pacbio/$Organism/"$Strain"_auto_cuttoff
qsub $ProgDir/subSpades_3lib_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
done
```

Contigs shorter than 500bp were removed from the assembly

```bash
for Contigs in $(ls assembly/spades_pacbio/*/*/contigs.fasta); do
AssemblyDir=$(dirname $Contigs)
mkdir $AssemblyDir/filtered_contigs
FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
$FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs/contigs_min_500bp.fasta
done
```

Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades_pacbio/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```



## Merging pacbio and hybrid assemblies

```bash
  # for PacBioAssembly in $(ls assembly/canu/*/*/polished/pilon.fasta); do
    # Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    # Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
Jobs=$(qstat | grep 'sub_pilon' | wc -l)
while [ $Jobs -gt 0 ]; do
printf "."
sleep 5m
Jobs=$(qstat | grep 'sub_pilon' | wc -l)
done		
for PacBioAssembly in $(ls assembly/canu/P.cactorum/414/414_canu.contigs.fasta); do
Organism=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $PacBioAssembly | rev | cut -f2 -d '/' | rev)
# HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
HybridAssembly=$(ls assembly/spades_pacbio/P.cactorum/414/contigs.fasta)
OutDir=assembly/merged_canu_spades/$Organism/$Strain
AnchorLength=500000
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
done
```

This merged assembly was polished using Pilon

```bash
for Assembly in $(ls -d assembly/merged_canu_spades/P.*/*/filtered_contigs/contigs_min_500bp_renamed.fasta | grep -w -e '414_v2'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/$Organism/414)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
TrimF2_Read=$(ls $IlluminaDir/F/414_run2_F_trim.fq.gz);
TrimR2_Read=$(ls $IlluminaDir/R/414_run2_R_trim.fq.gz);
TrimF3_Read=$(ls $IlluminaDir/F/414_170210_F_trim.fq.gz);
TrimR3_Read=$(ls $IlluminaDir/R/414_170210_R_trim.fq.gz);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
OutDir=assembly/merged_canu_spades/$Organism/$Strain/polished
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon_3_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
done
```

Contigs were renamed in accordance with ncbi recomendations.

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/*/*/polished/pilon.fasta | grep '414_v2'); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/merged_canu_spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -e '414_v2'); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

checking using busco

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -e '414_v2' | grep -v 'polished'); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
  OutDir=gene_pred/busco/$Organism/$Strain/assembly
  qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt | grep '414_v2'); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

```
P.cactorum	414_v2	278	2	23	303
```


Contigs were identified that had blast hits to non-phytophthora genomes

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -e '414_v2' | grep -v 'polished'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Exclude_db="bact,virus,hsref"
Exclude_db="paenibacillus"
Good_db="phytoph"
AssemblyDir=$(dirname $Assembly)
OutDir=$AssemblyDir/../deconseq_Paen
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
done
```

Results were summarised using the commands:

```bash
# for File in $(ls assembly/spades/P.*/*/deconseq/log.txt); do
for File in $(ls assembly/merged_canu_spades/P.*/*/deconseq_Paen/log.txt | grep '414_v2'); do
  Name=$(echo $File | rev | cut -f3 -d '/' | rev);
  Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
  Both=$(cat $File |cut -f2 | head -n2 | tail -n1);
  Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
  printf "$Name\t$Good\t$Both\t$Bad\n";
done
```

```
  414_v2	186	0	0
```

# Preliminary analysis

## Checking MiSeq coverage against P414 contigs

```bash
for Assembly in $(ls assembly/merged_canu_spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -e '414_v2' | grep -v 'polished'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
# IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
IlluminaDir=$(ls -d qc_dna/paired/$Organism/414)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz | head -n1 | tail -n1);
TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz | head -n1 | tail -n1);
TrimF2_Read=$(ls $IlluminaDir/F/*trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $IlluminaDir/R/*trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $IlluminaDir/F/*trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $IlluminaDir/R/*trim.fq.gz | head -n3 | tail -n1);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
InDir=$(dirname $Assembly)
OutDir=$InDir/aligned_MiSeq
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie_3lib.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
done
```

## Checking PacBio coverage against P414 contigs

The accuracy of PacBio assembly pipelines is currently unknown. To help identify
regions that may have been missassembled the pacbio reads were aligned back to
the assembled genome. Coverage was determined using bedtools genomecov and
regions with low coverage flagged using a python script flag_low_coverage.py.
These low coverage regions were visually inspected using IGV.

```bash
# for Assembly in $(ls assembly/merged_canu_spades/P.cactorum/414/filtered_contigs/contigs_min_500bp_renamed.fasta); do
for Assembly in $(ls assembly/canu/P.cactorum/*/filtered_contigs/contigs_min_500bp_renamed.fasta); do
  Reads=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq)
  OutDir=analysis/genome_alignment/bwa/P.cactorum/414/vs_414
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
  #
  # AlignedBam=$OutDir/414_contigs_renamed.fasta_aligned_sorted.bam.gz
  # CoverageTxt=$OutDir/414_bp_genome_cov.txt
  # bedtools genomecov -max 5 -bga -d -ibam $AlignedBam -g $Assembly > $CoverageTxt
  #
  # Threshold=5
  # FlaggedRegions=$OutDir/414_flagged_regions.txt
  # $ProgDir/flag_low_coverage.py --genomecov $CoverageTxt --min $Threshold > $FlaggedRegions
  done
```

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/merged_canu_spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -w -e '414_v2' | grep -v 'polished'); do
qsub $ProgDir/rep_modeling.sh $BestAss
qsub $ProgDir/transposonPSI.sh $BestAss
done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
  for RepDir in $(ls -d repeat_masked/P.*/*/filtered_contigs_repmask | grep -w -e '414_v2'); do
    Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
    RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
    TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    printf "$Organism\t$Strain\n"
    printf "The number of bases masked by RepeatMasker:\t"
    sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The number of bases masked by TransposonPSI:\t"
    sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The total number of masked bases are:\t"
    cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    echo
  done
```
```
  P.cactorum	414_v2
  The number of bases masked by RepeatMasker:	15459761
  The number of bases masked by TransposonPSI:	4939043
  The total number of masked bases are:	16794791
```


The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

  for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep '414_v2'); do
    OutDir=$(dirname $File)
    TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
    echo "$OutFile"
    bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
    echo "Number of masked bases:"
    cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
  done
  # The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
  for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep '414_v2'); do
    OutDir=$(dirname $File)
    TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
    echo "$OutFile"
    bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
  done
```

# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

## Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  for Genome in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
    echo $Genome;
    qsub $ProgDir/sub_cegma.sh $Genome dna;
  done
```

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/*/*/*_dna_cegma.completeness_report | grep -w -e '414_v2'); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
```

#Gene prediction

Gene prediction was performed for the P. cactorum genome. Two gene prediction
approaches were used:

Gene prediction using Braker1
Prediction of all putative ORFs in the genome using the ORF finder (atg.pl)
approach.


## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to Fusarium genomes.



#### Aligning published RNAseq data

* qc of RNA seq data was performed as part of sequencing the 10300 genome:

<!-- ```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNA in $(ls qc_rna/raw_rna/genbank/*/*/*_trim.fq.gz); do
      Timepoint=$(echo $RNA | rev | cut -f1 -d '/' | rev | sed 's/_trim.*//g')
      echo "$Timepoint"
      OutDir=alignment/$Organism/$Strain/$Timepoint
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
      qsub $ProgDir/tophat_alignment_unpaired.sh $Assembly $RNA $OutDir
    done
  done
``` -->
```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNA in $(ls qc_rna/raw_rna/genbank/*/*/*_trim.fq.gz); do
      Timepoint=$(echo $RNA | rev | cut -f1 -d '/' | rev | sed 's/_trim.*//g')
      echo "$Timepoint"
      Prefix="$Tiimepoint"
      OutDir=alignment/star/$Organism/"$Strain"/$Timepoint/$Prefix
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
      qsub $ProgDir/sub_star_unpaired.sh $Assembly $RNA $OutDir
    done
  done
```


### Aligning in house RNAseq data

make symbolic links to timecourse data

```bash
# Create symbolic links for all F read files
for File in $(ls /home/groups/harrisonlab/raw_data/raw_seq/fragaria/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07/*/*.gz | grep 'R1.fastq.gz'); do
  # echo $File;
  Sample=$(echo $File | rev | cut -d '/' -f2 | rev)
  echo "$Sample"
  OutDir=qc_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07/$Sample/F
  mkdir -p "$OutDir"
  cp -s $File $OutDir/.
done
# Create symbolic links for all R read files
for File in $(ls /home/groups/harrisonlab/raw_data/raw_seq/fragaria/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07/*/*.gz | grep 'R2.fastq.gz'); do
  # echo $File;
  Sample=$(echo $File | rev | cut -d '/' -f2 | rev)
  echo "$Sample"
  OutDir=qc_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07/$Sample/R
  mkdir -p "$OutDir"
  cp -s $File $OutDir/.
done
```


Perform qc of RNAseq timecourse data
```bash
  for FilePath in $(ls -d raw_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07/*); do
    echo $FilePath
    FileNum=$(ls $FilePath/F/*.gz | wc -l)
    for num in $(seq 1 $FileNum); do
      FileF=$(ls $FilePath/F/*.gz | head -n $num | tail -n1)
      FileR=$(ls $FilePath/R/*.gz | head -n $num | tail -n1)
      echo $FileF
      echo $FileR
      Jobs=$(qstat | grep 'rna_qc' | grep 'qw' | wc -l)
      while [ $Jobs -gt 16 ]; do
        sleep 5m
        printf "."
        Jobs=$(qstat | grep 'rna_qc' | grep 'qw' | wc -l)
      done		
      printf "\n"
      IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
      qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
    done
  done
```

<!-- Data quality was visualised using fastqc:
```bash
for RawData in $(ls qc_rna/paired/*/*/*/*.fq.gz); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
Jobs=$(qstat | grep 'fastqc' | wc -l)
while [ $Jobs -gt 20 ]; do
sleep 5m
printf "."
Jobs=$(qstat | grep 'fastqc' | wc -l)
done
sleep 1m
printf "\n"
echo $RawData;
qsub $ProgDir/run_fastqc.sh $RawData
done
``` -->

#### Aligning


```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNADir in $(ls -d qc_rna/paired/Transcriptome*/*); do
      FileNum=$(ls $RNADir/F/*_trim.fq.gz | wc -l)
      for num in $(seq 1 $FileNum); do
        while [ $Jobs -gt 1 ]; do
          sleep 1m
          printf "."
          Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        done
        printf "\n"
        FileF=$(ls $RNADir/F/*_trim.fq.gz | head -n $num | tail -n1)
        FileR=$(ls $RNADir/R/*_trim.fq.gz | head -n $num | tail -n1)
        echo $FileF
        echo $FileR
        Prefix=$(echo $FileF | rev | cut -f1 -d '/' | rev | sed "s/_R.*_trim.fq.gz//g")
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
        echo "$Timepoint"
        OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
        qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
      done
    done
  done
```

Alignment stats were collected:

```bash
for File in $(ls alignment/star/P.cactorum/414_v2/*/*/star_aligmentLog.final.out); do
Sample=$(echo $File | rev | cut -f2 -d '/' | rev);
ReadNumU=$(cat $File | grep 'Uniquely' | grep 'number' | cut -f2);
ReadPercU=$(cat $File | grep 'Uniquely' | grep '%' | cut -f2);
ReadNumM=$(cat $File | grep 'multiple' | grep 'Number' | cut -f2);
ReadPercM=$(cat $File | grep 'multiple' | grep '%' | cut -f2);
echo -e "$Sample""\t""$ReadNumU""\t""$ReadPercU""\t""$ReadNumM""\t""$ReadPercM";  
done
```


Alignments were concatenated prior to gene prediction

```bash

BamFiles=$(ls alignment/star/P.cactorum/414_v2/Sample_*/*/star_aligmentAligned.sortedByCoord.out.bam | grep -v -e 'PRO1467_S1_' -e 'PRO1467_S2_' -e 'PRO1467_S3_' -e 'PRO1467_S10_' -e 'PRO1467_S11_' -e 'PRO1467_S12_' | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=alignment/star/P.cactorum/414_v2/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/12_24hrs_concatenated.bam $BamFiles

BamFiles=$(ls alignment/star/P.cactorum/414_v2/*/*/star_aligmentAligned.sortedByCoord.out.bam | grep -v -e 'PRO1467_S1_' -e 'PRO1467_S2_' -e 'PRO1467_S3_' -e 'PRO1467_S10_' -e 'PRO1467_S11_' -e 'PRO1467_S12_' | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=alignment/star/P.cactorum/414_v2/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/12_24hrs_published_concatenated.bam $BamFiles
```



#### Braker prediction

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
    	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    	echo "$Organism - $Strain"
    	AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/12_24hrs_published_concatenated.bam)
    	OutDir=gene_pred/braker/$Organism/"$Strain"_braker_pacbio
    	GeneModelName="$Organism"_"$Strain"_braker_pacbio
    	rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_pacbio
    	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
    	qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  	done
```


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
    mkdir -p $OutDir
    AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/12_24hrs_published_concatenated.bam)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
  done
```

Secondly, genes were predicted using CodingQuary:

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred/final_genes/$Organism/$Strain
		CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
		qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
	done
```

Some gene models were manually removed. Note - this step has been inserted after
an initial prediction of the gene models and feedback has been received by ncbi.

```bash
# The following gene was removed from predicted gene models as a result of recomendation by ncbi
# CUFF_1141_1_145 - alternatively spliced gene with no stop codon (remember that '_' need to be replaced with '.')
# g7251 - overlapping exon with an RxLR gene
for BrakerGff in $(ls gene_pred/braker/P.*/*/*/augustus.gff3 | grep '414_v2'); do
  Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g')
  Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  BrakerGffEdits=$(echo $BrakerGff | sed 's/.gff3/_edits.gff3/g')
  cat $BrakerGff | grep -v -w 'CUFF.1141.1.145' | grep -v -w 'g7251' > $BrakerGffEdits
  CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
  CodingQuaryGffEdits=$(echo $CodingQuaryGff | sed 's/.gff3/_edits.gff3/g')
  PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
  PGNGffEdits=$(echo $PGNGff | sed 's/.gff3/_edits.gff3/g')
  cat $CodingQuaryGff | grep -v -w 'CUFF.1141.1.145' | grep -v -w 'g7251' > $CodingQuaryGffEdits
  cat $PGNGff | grep -v -w 'CUFF.1141.1.145' | grep -v -w 'g7251' > $PGNGffEdits
done
```


Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
for BrakerGff in $(ls gene_pred/braker/P.*/*/*/augustus_edits.gff3 | grep '414_v2'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/"$Strain"_contigs_softmasked.fa)
CodingQuaryGff=$(ls gene_pred/codingquary/$Organism/$Strain/out/PredictedPass_edits.gff3)
PGNGff=$(ls gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass_edits.gff3)
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/final_genes/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
# -
# This section is edited
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuary_unspliced.gff3
$ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuary_unspliced.gff3 > $FinalDir/final_genes_CodingQuary.gff3
# -
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

GffBraker=$FinalDir/final_genes_Braker.gff3
GffQuary=$FinalDir/final_genes_CodingQuary.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended
done
```

The final number of genes per isolate was observed using:
```bash
  for DirPath in $(ls -d gene_pred/final_genes/P.*/*/final | grep '414_v2'); do
    echo $DirPath;
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
    echo "";
  done
```

```
gene_pred/final_genes/P.cactorum/414_v2/final
24892
4912
29804
```

<!--
# Assessing gene space in predicted transcriptomes

```bash
  for Transcriptome in $(ls gene_pred/final_genes/*/*/final/final_genes_combined.gene.fasta | grep '414_v2'); do
    Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    # BuscoDB="Fungal"
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
    OutDir=gene_pred/busco/$Organism/$Strain/genes
    qsub $ProgDir/sub_busco2.sh $Transcriptome $BuscoDB $OutDir
  done
```
```bash
  for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt | grep '414_v2'); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

```
P.cactorum	414_v2	277	4	22	303
``` -->

## Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing
Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
  	qsub $ProgDir/run_ORF_finder.sh $Assembly
  done
```


The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
	ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
  for OrfGff in $(ls gene_pred/ORF_finder/P.*/*/*_ORF.gff | grep -v 'atg' | grep -w -e '414_v2'); do
    echo "$OrfGff"
  	OrfGffMod=$(echo $OrfGff | sed 's/.gff/.gff3/g')
  	$ProgDir/gff_corrector.pl $OrfGff > $OrfGffMod
  done
```
<!--
#Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	cd /home/groups/harrisonlab/project_files/idris
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
Strain=$(echo $Proteome | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
echo $Strain
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteome $InterProRaw
done
```


## B) SwissProt

```bash
  for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/swissprot/$Organism/$Strain
    SwissDbDir=../../uniprot/swissprot
    SwissDbName=uniprot_sprot
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
    qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
  done
```

## C) Summarising annotation in annotation table

```bash
  for GeneGff in $(ls gene_pred/final_genes/*/*/*/final_genes_appended.gff3 | grep '414_v2'); do
    Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
    Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked.fa)
    InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
    SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
    OutDir=gene_pred/annotation/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/annotation_tables
    $ProgDir/build_annot_tab.py --genome $Assembly --genes_gff $GeneGff --InterPro $InterPro --Swissprot $SwissProt > $OutDir/"$Strain"_annotation.tsv
  done
```
 -->

#Genomic analysis

## RxLR genes

Putative RxLR genes were identified within Augustus gene models using a number
of approaches:

 * A) From Augustus gene models - Signal peptide & RxLR motif  
 * B) From Augustus gene models - Hmm evidence of WY domains  
 * C) From Augustus gene models - Hmm evidence of RxLR effectors
 * D) From Augustus gene models - Hmm evidence of CRN effectors  
 * E) From ORF fragments - Signal peptide & RxLR motif  
 * F) From ORF fragments - Hmm evidence of WY domains  
 * G) From ORF fragments - Hmm evidence of RxLR effectors


 ### A) From Augustus gene models - Signal peptide & RxLR motif

 Required programs:
  * SigP
  * biopython

#### A.1) Signal peptide prediction using SignalP 2.0

 Proteins that were predicted to contain signal peptides were identified using
 the following commands:

 ```bash
for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/final_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_final
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_final_*); do
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
while [ $Jobs -gt 20 ]; do
sleep 1
printf "."
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
done
printf "\n"
echo $File
qsub $ProgDir/pred_sigP.sh $File
qsub $ProgDir/pred_sigP.sh $File signalp-3.0
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
 ```

 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

 ```bash
for SplitDir in $(ls -d gene_pred/final_split/P.*/* | grep '414_v2'); do
Strain=$(echo $SplitDir | cut -d '/' -f4)
Organism=$(echo $SplitDir | cut -d '/' -f3)
echo "$Organism - $Strain"
for SigpDir in $(ls -d gene_pred/final_sig* | cut -f2 -d'/'); do
echo "$SigpDir"
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
for GRP in $(ls -l $SplitDir/*_final_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.aa";
InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp_neg.aa";  
InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.tab";
InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.txt";
done
cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
done
done
```


#### B.2) Prediction using Phobius

 Secreted proteins were also predicted using Phobius

 ```bash
for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/phobius/$Organism/$Strain
mkdir -p $OutDir
phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
$ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
cat $OutDir/"$Strain"_phobius.fa | grep '>' | cut -f1 | sed 's/>//g' > $OutDir/"$Strain"_phobius_headers.txt
done
 ```


Secreted proteins from different sources were combined into a single file:

```bash
for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/combined_sigP/$Organism/$Strain
mkdir -p $OutDir
echo "The following number of sequences were predicted as secreted:"
cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | wc -l
echo "This represented the following number of unique genes:"
cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
done
```

```
  P.cactorum - 414_v2
  The following number of sequences were predicted as secreted:
  11279
  This represented the following number of unique genes:
  3689
```


Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
qsub $ProgDir/submit_TMHMM.sh $Proteome
done
```
<-
Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt | grep '414_v2'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
NonTmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $NonTmHeaders
SigP=$(ls gene_pred/final_sigP/$Organism/$Strain/*_sp.aa | grep -v 'neg')
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $NonTmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
echo "Number of SigP proteins:"
cat $SigP | grep '>' | wc -l
echo "Number without transmembrane domains:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
echo "Number of gene models:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l

# A text file was also made containing headers of proteins testing +ve
PosFile=$(ls gene_pred/trans_mem/$Organism/$Strain/"$Strain"_TM_genes_pos.txt)
TmHeaders=$(echo $PosFile | sed 's/.txt/_headers.txt/g')
cat $PosFile | cut -f1 > $TmHeaders

done
```

Proteins containing GPI anchors were also removed using GPIsom


These proteins were identified through submitting the combined protein file to
the webserver at: http://gpi.unibe.ch

An output directory was made to download the file to:
"GPI anchored (C&N-term signal) (SignalP):"

```bash
  for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/trans_mem/$Organism/$Strain/GPIsom
    mkdir -p $OutDir
  done
```

Results were pasted into the file:

```bash
  nano gene_pred/trans_mem/P.cactorum/414_v2/GPIsom/GPI_pos.fa
```

Those proteins with GPI anchors were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/GPIsom/GPI_pos.fa | grep '414_v2'); do
Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
TmHeaders=$(echo "$File" | sed 's/.fa/.txt/g')
cat $File | grep '>' | cut -f1 -d ' ' | sed 's/>//g' > $TmHeaders
SigP=$(ls gene_pred/final_sigP/$Organism/$Strain/*_sp_no_trans_mem.aa)
SigPHeaders=gene_pred/final_sigP/$Organism/$Strain/"$Strain"_sp_no_trans_mem_headers.txt
cat $SigP | grep '>' | cut -f1 | sed 's/>//g'> $SigPHeaders
GoodHeaders=$(echo "$File" | sed 's/_pos.fa/_neg.txt/g')
cat $SigPHeaders | grep -v -f $TmHeaders > $GoodHeaders
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# cat $SigP | grep -v -A1 -f $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $GoodHeaders  > $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa
echo "Number of SigP proteins:"
cat $SigP | grep '>' | wc -l
echo "Number with GPI anchors in entire proteome:"
cat $TmHeaders | wc -l
echo "Number without GPI anchors:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa | grep '>' | wc -l
echo "Number of gene models:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l
done
```

```
  P.cactorum - 414_v2
  Number of SigP proteins:
  1796
  Number with GPI anchors in entire proteome:
  560
  Number without GPI anchors:
  1540
  Number of gene models:
  1536
```



### C.i) From Augustus gene models - Effector-like structure identification using EffectorP

Required programs:
 * EffectorP.py

```bash
for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
BaseName="$Organism"_"$Strain"_EffectorP
OutDir=analysis/effectorP/$Organism/$Strain
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
done
```

Those genes that were predicted as secreted and tested positive by effectorP
were identified:

```bash
for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt | grep '414_v2'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
cat $File | grep 'Effector' | cut -f1 > $Headers
printf "EffectorP headers:\t"
cat $Headers | wc -l
Secretome=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_secreted.fa)
OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
printf "Secreted effectorP headers:\t"
cat $OutFileHeaders | wc -l
Gff=$(ls gene_pred/final_genes/$Organism/$Strain/*/final_genes_appended.gff3)
EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
done
```
```
P.cactorum - 414_v2
EffectorP headers:	18326
Secreted effectorP headers:	966
```


## C.ii) SSCP

Small secreted cysteine rich proteins were identified within secretomes. These
proteins may be identified by EffectorP, but this approach allows direct control
over what constitutes a SSCP.

```bash

for Secretome in $(ls gene_pred/combined_sigP/*/*/*_secreted.fa | grep -v 'all_secreted' |  grep '414_v2'); do
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/sscp/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_all_results.fa
cat $OutDir/"$Strain"_sscp_all_results.fa | grep 'Yes' > $OutDir/"$Strain"_sscp.fa
printf "number of SSC-rich genes:\t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' | cut -f1 -d '.' | sort | uniq | wc -l
printf "Number of effectors predicted by EffectorP:\t"
EffectorP=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP_secreted_headers.txt)
cat $EffectorP | wc -l
printf "Number of SSCPs predicted by both effectorP and this approach: \t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' > $OutDir/"$Strain"_sscp_headers.txt
cat $OutDir/"$Strain"_sscp_headers.txt $EffectorP | cut -f1 | sort | uniq -d | wc -l
echo ""
done > tmp.txt
cat tmp.txt
```

```
P.cactorum - 414_v2
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	435
number of SSC-rich genes:	434
Number of effectors predicted by EffectorP:	966
Number of SSCPs predicted by both effectorP and this approach: 	255
```


## D) CAZY proteins


Cazy HMMs may need an older version of HMMer to run

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/CAZY/$Organism/$Strain
mkdir -p $OutDir
Prefix="$Strain"_CAZY
CazyHmm=../../dbCAN/dbCAN-fam-HMMs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
done
```


The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep '414_v2'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $File)
echo "$Organism - $Strain"
ProgDir=/home/groups/harrisonlab/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
printf "number of CAZY genes identified:\t"
cat $CazyHeaders | wc -l
Gff=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
CazyGff=$OutDir/"$Strain"_CAZY.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

SecretedProts=$(ls gene_pred/combined_sigP/$Organism/$Strain/*_all_secreted.fa)
SecretedHeaders=$(echo $SecretedProts | sed 's/.fa/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
printf "number of Secreted CAZY genes identified:\t"
cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
done > tmp.txt
```
```
P.cactorum - 414_v2
number of CAZY genes identified:	679
number of Secreted CAZY genes identified:	398
```

Note - the CAZY genes identified may need further filtering based on e value and
cuttoff length - see below:

Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage

* For fungi, use E-value < 1e-17 and coverage > 0.45

* The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)


# D) Prediction of RxLRs from - Augustus/CodingQuary gene models

 The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

 The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash

for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep '414_v2'); do
Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
Proteome=$(ls gene_pred/final_genes/$Organism/$Strain/*/final_genes_combined.pep.fasta)
Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
mkdir -p $OutDir;
printf "\nstrain: $Strain\tspecies: $Organism\n";
printf "the total number of SigP gene is:\t";
cat $Secretome | grep '>' | wc -l;
printf "the number of unique SigP gene is:\t";
cat $Secretome | grep '>' | cut -f1 | tr -d ' '| sort | uniq | wc -l;

printf "the number of SigP-RxLR genes are:\t";
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_all_secreted_RxLR_regex.fa;
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $OutDir/"$Strain"_RxLR_regex.txt
cat $OutDir/"$Strain"_RxLR_regex.txt | wc -l

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa


printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $OutDir/"$Strain"_RxLR_EER_regex.txt
cat $OutDir/"$Strain"_RxLR_EER_regex.txt | wc -l
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa

printf "\n"

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_regex.txt
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_EER_regex.txt

cat $Gff | grep -w -f $OutDir/"$Strain"_RxLR_regex.txt> $OutDir/"$Strain"_RxLR_regex.gff3
cat $Gff | grep -w -f $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.gff3
done > tmp.txt
```

```
  strain: 414_v2	species: P.cactorum
  the total number of SigP gene is:	11279
  the number of unique SigP gene is:	3689
  the number of SigP-RxLR genes are:	336
  the number of SigP-RxLR-EER genes are:	148
```


### E) From gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search predicted proteins. These were run with the following commands:


```bash
for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep '414_v2'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=$(ls /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm)
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_WY_hmmer.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_WY_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_WY_hmmer_headers.txt
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
echo "Total genes with WY domains"
cat $OutDir/$Headers | sort | uniq | wc -l
done
```

P.cactorum 414_v2
Initial search space (Z):              11279  [actual number of targets]
Domain search space  (domZ):             407  [number of targets reported over threshold]


### G) From Secreted gene models - Hmm evidence of RxLR effectors

```bash
  for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
    HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
    mkdir -p $OutDir
    HmmResults="$Strain"_RxLR_hmmer.txt
    hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
    echo "$Organism $Strain"
    cat $OutDir/$HmmResults | grep 'Initial search space'
    cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
    HmmFasta="$Strain"_RxLR_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
    Headers="$Strain"_RxLR_hmmer_headers.txt
    cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' | sort | uniq > $OutDir/$Headers
    Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
    cat $Gff | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
  done
```

```
  P.cactorum 414_v2
  Initial search space (Z):              29806  [actual number of targets]
  Domain search space  (domZ):             152  [number of targets reported over threshold]
```


### F) Combining RxLRs from Regex and hmm searches


The total RxLRs are

```bash
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_RxLR_EER_regex.txt | grep '414_v2'); do
Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
Proteome=$(ls gene_pred/final_genes/$Organism/$Strain/*/final_genes_combined.pep.fasta)
HmmRxLR=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer_headers.txt
echo "$Organism - $Strain"
echo "Number of RxLRs identified by Regex:"
cat $RegexRxLR | sort | uniq | wc -l
echo "Number of RxLRs identified by Hmm:"
cat $HmmRxLR | sort | uniq | wc -l
echo "Number of RxLRs in combined dataset:"
cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
# echo "Number of RxLRs in both datasets:"
# cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
echo ""
# echo "Extracting RxLRs from datasets"
OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
mkdir -p $OutDir
cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_RxLR_headers.txt
Gff=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
cat $Gff | grep -w -f $OutDir/"$Strain"_total_RxLR_headers.txt > $OutDir/"$Strain"_total_RxLR.gff
echo "Number of genes in the extracted gff file:"
cat $OutDir/"$Strain"_total_RxLR.gff | grep -w 'gene' | wc -l
done
```

```
P.cactorum - 414_v2
Number of RxLRs identified by Regex:
148
Number of RxLRs identified by Hmm:
152
Number of RxLRs in combined dataset:
181

Number of genes in the extracted gff file:
181
```
<!--
72 of 137 RxLRs were in the predicted 797 effectorP genes.


```bash
  cat analysis/RxLR_effectors/combined_evidence/P.idaei/SCRP370/SCRP370_total_RxLR_headers.txt $OutFileHeaders | sort | cut -f1 -d '.' | uniq -d | wc -l
```
-->

### D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:


```bash
  HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
  LFLAK_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm)
  DWL_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm)
  for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep '414_v2'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
    mkdir -p $OutDir
    echo "$Organism - $Strain"
    # Run hmm searches LFLAK domains
    CrinklerProts_LFLAK=$OutDir/"$Strain"_pub_CRN_LFLAK_hmm.txt
    hmmsearch -T0 $LFLAK_hmm $Proteome > $CrinklerProts_LFLAK
    cat $CrinklerProts_LFLAK | grep 'Initial search space'
    cat $CrinklerProts_LFLAK | grep 'number of targets reported over threshold'
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
    $ProgDir/hmmer2fasta.pl $CrinklerProts_LFLAK $Proteome > $OutDir/"$Strain"_pub_CRN_LFLAK_hmm.fa
    # Run hmm searches DWL domains
    CrinklerProts_DWL=$OutDir/"$Strain"_pub_CRN_DWL_hmm.txt
    hmmsearch -T0 $DWL_hmm $Proteome > $CrinklerProts_DWL
    cat $CrinklerProts_DWL | grep 'Initial search space'
    cat $CrinklerProts_DWL | grep 'number of targets reported over threshold'
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
    $ProgDir/hmmer2fasta.pl $CrinklerProts_DWL $Proteome > $OutDir/"$Strain"_pub_CRN_DWL_hmm.fa
    # Identify the genes detected in both models
    cat $OutDir/"$Strain"_pub_CRN_LFLAK_hmm.fa $OutDir/"$Strain"_pub_CRN_DWL_hmm.fa | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $OutDir/"$Strain"_pub_CRN_LFLAK_DWL.txt
    cat $OutDir/"$Strain"_pub_CRN_LFLAK_DWL.txt | wc -l
  done
```

```
  P.cactorum - 414_v2
  Initial search space (Z):              29806  [actual number of targets]
  Domain search space  (domZ):             225  [number of targets reported over threshold]
  Initial search space (Z):              29806  [actual number of targets]
  Domain search space  (domZ):             188  [number of targets reported over threshold]
  171
```
<!--
23 of 92 P.idaei CRNs were in the effectorP dataset.
```bash
cat analysis/CRN_effectors/hmmer_CRN/P.idaei/SCRP370/SCRP370_pub_CRN_LFLAK_DWL.txt $OutFileHeaders | sort | cut -f1 -d '.' | uniq -d | wc -l
``` -->

Extract gff annotations for Crinklers:

```bash
  for CRNlist in $(ls analysis/CRN_effectors/hmmer_CRN/*/*/*_pub_CRN_LFLAK_DWL.txt | grep -e '414_v2'); do
    Strain=$(echo $CRNlist | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $CRNlist | rev | cut -f3 -d '/' | rev)
    OutName=$(echo $CRNlist | sed 's/.txt/.gff/g')
    echo "$Organism - $Strain"
    Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
    cat $CRNlist | sed -r 's/\.t.$//g' > tmp.txt
    cat $Gff | grep -w -f tmp.txt > $OutName
    rm tmp.txt
  done
```


### E) From ORF gene models - Signal peptide & RxLR motif

Required programs:
 * SigP
 * Phobius
 * biopython


#### E.1) Prediction using SignalP
Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep '414_v2'); do
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
SplitDir=gene_pred/ORF_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_ORF_preds
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_ORF_preds_*); do
Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
while [ $Jobs -gt 6 ]; do
sleep 1
printf "."
Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
done		
printf "\n"
echo $File
qsub $ProgDir/pred_sigP.sh $File
qsub $ProgDir/pred_sigP.sh $File signalp-3.0
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```


 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

```bash
for SplitDir in $(ls -d gene_pred/ORF_split/*/* | grep '414_v2'); do
Strain=$(echo $SplitDir | cut -d '/' -f4)
Organism=$(echo $SplitDir | cut -d '/' -f3)
echo "$Organism - $Strain"
for SigpDir in $(ls -d gene_pred/ORF_sig* | cut -f2 -d'/'); do
echo "$SigpDir"
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
for GRP in $(ls -l $SplitDir/*_ORF_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.aa";
InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp_neg.aa";
InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.tab";
InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.txt";
done
cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
done
done
```

#### E.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
	for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep '414_v2'); do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/phobius/$Organism/$Strain
    mkdir -p $OutDir
  	phobius.pl $Proteome > $OutDir/"$Strain"_phobius_ORF.txt
  	cat $OutDir/"$Strain"_phobius_ORF.txt | grep -B1 'SIGNAL' | grep 'ID' | sed s'/ID.*g/g/g' > $OutDir/"$Strain"_phobius_headers_ORF.txt
  done
```

Secreted proteins from different sources were combined into a single file:

```bash
for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep '414_v2'); do
Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/combined_sigP_ORF/$Organism/$Strain
mkdir -p $OutDir
echo "The following number of sequences were predicted as secreted:"
# cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa > $OutDir/"$Strain"_all_secreted.fa
cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | tr -d '>' | tr -d ' ' | sed "s/HMM_score\t/HMM_score=\t/g" > $OutDir/"$Strain"_all_secreted_headers.txt
cat $OutDir/"$Strain"_all_secreted_headers.txt | wc -l
echo "This represented the following number of unique genes:"
# cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
done
```


#### E.3) Prediction of RxLRs


Names of ORFs containing signal peptides were extracted from fasta files. This
included information on the position and hmm score of RxLRs.

```bash
  for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep '414_v2'); do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/ORF_sigP/$Organism/$Strain
  	SigP_headers=$OutDir/"$Strain"_ORF_sp_names.txt
  	cat $Proteome | grep '>' | sed -r 's/>//g' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | sed 's/--//g' > $SigP_headers
  done
```

Due to the nature of predicting ORFs, some features overlapped with one another.
A single ORF was selected from each set of overlapped ORFs. This was was
selected on the basis of its SignalP Hmm score. Biopython was used to identify
overlaps and identify the ORF with the best signalP score.

```bash
for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff3 | grep -v -e 'atg' | grep '414_v2'); do
Organism=$(echo $ORF_Gff | rev |  cut -d '/' -f3 | rev) ;
Strain=$(echo $ORF_Gff | rev | cut -d '/' -f2 | rev);
OutDir=$(ls -d gene_pred/combined_sigP_ORF/$Organism/$Strain)
echo "$Organism - $Strain"
# SigP_fasta=$(ls $OutDir/"$Strain"_all_secreted.fa)
SigP_headers=$(ls $OutDir/"$Strain"_all_secreted_headers.txt)
ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)

SigP_Gff=$OutDir/"$Strain"_all_secreted_unmerged.gff
SigP_Merged_Gff=$OutDir/"$Strain"_all_secreted_merged.gff
SigP_Merged_txt=$OutDir/"$Strain"_all_secreted_merged.txt
SigP_Merged_AA=$OutDir/"$Strain"_all_secreted_merged.aa

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $SigP_headers $ORF_Gff SigP Name > $SigP_Gff
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $SigP_Gff --db sigP_ORF.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF.db --id sigP_ORF --out sigP_ORF_merged.db --gff > $SigP_Merged_Gff
cat $SigP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d'=' | rev > $SigP_Merged_txt
# $ProgDir/extract_from_fasta.py --fasta $SigP_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
done
```



The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep '414_v2'); do
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
mkdir -p $OutDir;
printf "\nstrain: $Strain\tspecies: $Organism\n";
printf "the number of SigP gene is:\t";
cat $Secretome | grep '>' | cut -f1 | sort | uniq | wc -l
printf "the number of SigP-RxLR genes are:\t";
$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa;
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt
cat $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l
printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' '> $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt  $SigP_Merged_Gff 	RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
SigP_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt  $SigP_Gff	RxLR_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt  $SigP_Gff	RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.gff
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt $SigP_Merged_Gff RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex.gff
RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_regex_merged.gff
RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_regex_merged.txt
RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_regex_merged.aa
RxLR_EER_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.gff
RxLR_EER_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.txt
RxLR_EER_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff --db sigP_ORF_RxLR.db
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.gff --db sigP_ORF_RxLR_EER.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR.db --id sigP_ORF_RxLR --out sigP_ORF_RxLR_merged.db --gff > $RxLR_Merged_Gff
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR_EER.db --id sigP_ORF_RxLR_EER --out sigP_ORF_RxLR_EER_merged.db --gff > $RxLR_EER_Merged_Gff
cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
cat $RxLR_EER_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_EER_Merged_txt
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_EER_Merged_txt > $RxLR_EER_Merged_AA
printf "Merged RxLR regex proteins:\t"
cat $RxLR_Merged_AA | grep '>' | wc -l
printf "Merged RxLR-EER regex proteins:\t"
cat $RxLR_EER_Merged_AA | grep '>' | wc -l
printf "\n"
done
```

```
strain: 414_v2	species: P.cactorum
the number of SigP gene is:	54092
the number of SigP-RxLR genes are:	1924
the number of SigP-RxLR-EER genes are:	223
Merged RxLR regex proteins:	1511
Merged RxLR-EER regex proteins:	196
```

Quantification of ORF RxLRs was also performed.

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w -e '414_v2'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
InGff=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain*/*_ORF_RxLR_regex_merged.gff)
InGff_features=$(echo $InGff | sed 's/.gff/_features.gff/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/add_ORF_features.pl $InGff $Assembly >> $InGff_features
for BamFile in $(ls ../../../..//sobczm/popgen/rnaseq/vesca_*/pcac/star_aligmentAligned.sortedByCoord.out.bam); do
OutDir=alignment/star/$Organism/$Strain/featureCounts_maria_alignment_ORF_RxLR
mkdir -p $OutDir
Prefix=$(echo $BamFile | rev | cut -f3 -d '/' | rev | sed 's/vesca_//g')
Jobs=$(qstat | grep 'sub_fea' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
  sleep 1m
  printf "."
  Jobs=$(qstat | grep 'sub_fea' | grep 'qw'| wc -l)
done
printf "\n"
echo $Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
# echo "$BamFile $Gff $OutDir $Prefix"
qsub $ProgDir/sub_featureCounts.sh $BamFile $InGff_features $OutDir $Prefix
done
done
```

Those RxLRs with evidence of expression in planta (fpkm >5) were included in
ORF RxLR gene models:

```bash
  Organism="P.cactorum"
  Strain="414_v2"
  OutDir=$(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain)
  for FeatureCounts in $(ls alignment/star/$Organism/$Strain/featureCounts_maria_alignment_ORF_RxLR/*_featurecounts.txt); do
    OutFile=$(echo $FeatureCounts | sed 's/.txt/_fpkm_>5.txt/g')
  cat $FeatureCounts | tail -n+3 | cut -f1,7 | grep -v -e '\s0$' -e '\s1$' -e '\s2$' -e '\s3$' -e '\s4$' > $OutFile
  cat $OutFile
  done | cut -f1 | sort | uniq > $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers.txt
```



```bash
echo "Number of secreted ORFs"
cat gene_pred/ORF_signalp-4.1/P.cactorum/414_v2/414_v2_aug_sp.aa gene_pred/ORF_sigP/P.cactorum/414_v2/414_v2_aug_sp.aa | grep '>' | cut -f1 | sort | uniq | wc -l
```


### E5) From ORF gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search ORFs predicted with atg.pl. These were run with the following commands:


```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/"$Strain"_all_secreted_merged.aa | grep '414_v2'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=$(ls /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm)
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_ORF_WY_hmmer.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_ORF_WY_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_ORF_WY_hmmer_headers.txt
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
SigP_Merged_Gff=$(ls gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_merged.gff)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_WY_hmmer.gff
done
```

P.cactorum 414_v2
Initial search space (Z):              21882  [actual number of targets]
Domain search space  (domZ):             109  [number of targets reported over threshold]



### E6) From ORF gene models - Hmm evidence of RxLR effectors

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep '414_v2'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=$(ls /home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm)
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_ORF_RxLR_hmmer_unmerged.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_ORF_RxLR_hmmer_headers_unmerged.txt
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
SigP_Gff=$(ls gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_RxLR_hmmer_unmerged.gff3
RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.gff
RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.txt
RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_hmmer_unmerged.gff3 --db sigP_ORF_RxLR_hmm.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR_hmm.db --id sigP_ORF_RxLR_hmm --out sigP_ORF_RxLR_hmm_merged.db --gff > $RxLR_Merged_Gff
cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
printf "Merged RxLR-EER Hmm proteins:\t"
cat $RxLR_Merged_AA | grep '>' | wc -l
done
```
```
  P.cactorum 414_v2
  Initial search space (Z):              65593  [actual number of targets]
  Domain search space  (domZ):             498  [number of targets reported over threshold]
  Merged RxLR-EER Hmm proteins:	154
```


### E7) Combining RxLRs from Regex and hmm searches


The total RxLRs are

```bash
  for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_ORF_RxLR_EER_regex_merged.txt | grep '414_v2'); do
    Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
    Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF.gff3)
    Proteome=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
    HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_ORF_RxLR_hmm_merged.txt)
    echo "$Organism - $Strain"
    echo "Number of RxLRs identified by Regex:"
    cat $RegexRxLR | sort | uniq | wc -l
    echo "Number of RxLRs identified by Hmm:"
    cat $HmmRxLR | sort | uniq | wc -l
    echo "Number of RxLRs in combined dataset:"
    cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
    # echo "Number of RxLRs in both datasets:"
    # cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
    echo ""
    # echo "Extracting RxLRs from datasets"
    OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
    mkdir -p $OutDir
    cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
    # cat $Gff | grep -w -f $OutDir/"$Strain"_total_ORF_RxLR_headers.txt > $OutDir/"$Strain"_total_ORF_RxLR.gff
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_total_ORF_RxLR_headers.txt $Gff ORF_RxLR Name Augustus > $OutDir/"$Strain"_total_ORF_RxLR.gff
    echo "Number of genes in the extracted gff file:"
    cat $OutDir/"$Strain"_total_ORF_RxLR.gff | grep -w 'gene' | wc -l
  done
```

```
P.cactorum - 414_v2
Number of RxLRs identified by Regex:
196
Number of RxLRs identified by Hmm:
154
Number of RxLRs in combined dataset:
214
Number of genes in the extracted gff file:
214
```

## 4.2.c Analysis of RxLR effectors - merger of Augustus / published genes with ORFs

Intersection between the coodinates of putative RxLRs from gene models and ORFs
were identified to determine the total number of RxLRs predicted in these
genomes.

The RxLR effectors from both Gene models and ORF finding approaches were
combined into a single file.

This step was complicated by the inconsistency in downloaded gff files for gene
models.


```bash
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/*/* | grep '414_v2'); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$(ls gene_pred/final_genes/$Species/$Strain/final/final_genes_appended.gff3)
AugRxLR=$(ls $MergeDir/"$Strain"_total_RxLR.gff)
AugTxt=$(ls $MergeDir/"$Strain"_total_RxLR_headers.txt)
AugFa=$(ls gene_pred/final_genes/"$Species"/"$Strain"/final/final_genes_combined.pep.fasta)

ORFGff=$(ls $MergeDir/"$Strain"_total_ORF_RxLR.gff)
ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
ORFsTxt=$(ls $MergeDir/"$Strain"_total_ORF_RxLR_headers.txt)

ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_RxLR_EER_motif_hmm.gff
AugInORFs=$MergeDir/"$Strain"_AugInORFs_RxLR_EER_motif_hmm.gff
ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_RxLR_EER_motif_hmm.gff
# RxLR ORFs may be located in predicted gene models not identified as RxLR
# proteins - manual inspection suggests that many of the original gene models
# appear have RNAseq data aligning accross their entire length
ORFsUniqNoGenes=$MergeDir/"$Strain"_ORFsUniq_excl_genes_RxLR_EER_motif_hmm.gff
# ORFsUniqInGenes=$MergeDir/"$Strain"_ORFsUniq_in_genes_RxLR_EER_motif_hmm.gff
AugUniq=$MergeDir/"$Strain"_Aug_Uniq_RxLR_EER_motif_hmm.gff
TotalRxLRsTxt=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.txt
TotalRxLRsGff=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.gff

bedtools intersect -wa -u -a $ORFGff -b $AugRxLR > $ORFsInAug
bedtools intersect -wa -u -a $AugRxLR -b $ORFGff > $AugInORFs
bedtools intersect -v -wa -a $ORFGff -b $AugRxLR > $ORFsUniq
bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniqNoGenes
bedtools intersect -v -wa -a $AugRxLR -b $ORFGff > $AugUniq

echo "$Species - $Strain"
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l
echo "The number of RxLRs unique to ORF models and not intersecting nonRxLR gene models:"
cat $ORFsUniqNoGenes | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
echo "The total number of putative RxLRs are:"
cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalRxLRsTxt
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
# cat $ORFsUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f3 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
cat $ORFsUniqNoGenes | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f3 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
cat $TotalRxLRsTxt | wc -l

# Later stages needed lists including ORF 'IDs' rather than names
cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm_ID.txt
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm_ID.txt
# cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm_ID.txt
cat $ORFsUniqNoGenes | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm_ID.txt


# cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalRxLRsTxt > $TotalRxLRsGff
cat $AugInORFs $AugUniq $ORFsUniqNoGenes | grep -w -f $TotalRxLRsTxt > $TotalRxLRsGff

RxLRsFa=$MergeDir/"$Strain"_final_RxLR_EER.fa
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/unwrap_fasta.py --inp_fasta $AugFa | grep -A1 -w -f $AugTxt | grep -v -E '^--$' > $RxLRsFa
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalRxLRsTxt > $RxLRsFa
# echo "$Strain"
$ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalRxLRsTxt >> $RxLRsFa
# echo "$Strain"
echo "The number of sequences extracted is"
cat $RxLRsFa | grep '>' | wc -l
done
```

```
P.cactorum - 414_v2
The number of ORF RxLRs overlapping Augustus RxLRs:
159
The number of Augustus RxLRs overlapping ORF RxLRs:
159
The number of RxLRs unique to ORF models:
55
The number of RxLRs unique to ORF models and not intersecting nonRxLR gene models:
20
The number of RxLRs unique to Augustus models:
22
The total number of putative RxLRs are:
201
```


### H) From ORF gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in ORF gene models. This was done with the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/*/*/*.aa_cat.fa | grep '414_v2'); do
# Setting variables
Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
mkdir -p $OutDir
# Hmmer variables
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
# Searches for LFLAK domain
LFLAK_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm
HmmResultsLFLAK="$Strain"_ORF_CRN_LFLAK_unmerged_hmmer.txt
hmmsearch -T 0 $LFLAK_hmm $Proteome > $OutDir/$HmmResultsLFLAK
echo "Searching for LFLAK domains in: $Organism $Strain"
cat $OutDir/$HmmResultsLFLAK | grep 'Initial search space'
cat $OutDir/$HmmResultsLFLAK | grep 'number of targets reported over threshold'
HmmFastaLFLAK="$Strain"_ORF_CRN_LFLAK_unmerged_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResultsLFLAK $Proteome > $OutDir/$HmmFastaLFLAK
# Searches for DWL domain
DWL_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm
HmmResultsDWL="$Strain"_ORF_CRN_DWL_unmerged_hmmer.txt
hmmsearch -T 0 $DWL_hmm $Proteome > $OutDir/$HmmResultsDWL
echo "Searching for DWL domains in: $Organism $Strain"
cat $OutDir/$HmmResultsDWL | grep 'Initial search space'
cat $OutDir/$HmmResultsDWL | grep 'number of targets reported over threshold'
HmmFastaDWL="$Strain"_ORF_CRN_DWL_unmerged_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResultsDWL $Proteome > $OutDir/$HmmFastaDWL
# Identify ORFs found by both models
CommonHeaders=$OutDir/"$Strain"_ORF_CRN_DWL_LFLAK_unmerged_headers.txt
cat $OutDir/$HmmFastaLFLAK $OutDir/$HmmFastaDWL | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $CommonHeaders
echo "The number of CRNs common to both models are:"
cat $CommonHeaders | wc -l
# The sequences will be merged based upon the strength of their DWL domain score
# For this reason headers as they appear in the DWL fasta file were extracted
Headers=$OutDir/"$Strain"_CRN_hmmer_unmerged_headers.txt
cat $OutDir/$HmmFastaDWL | grep '>' | grep -w -f $CommonHeaders | tr -d '>' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | tr -d '-' | sed 's/hmm_score/HMM_score/g' > $Headers
# As we are dealing with JGI and Broad sequences, some features need formatting:
ORF_Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/*_ORF.gff3)
# Gff features were extracted for each header
CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $Headers $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
# Gff features were merged based upon the DWL hmm score
DbDir=analysis/databases/$Organism/$Strain
mkdir -p $DbDir
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $CRN_unmerged_Gff --db $DbDir/CRN_ORF.db
CRN_Merged_Gff=$OutDir/"$Strain"_CRN_merged_hmmer.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp $DbDir/CRN_ORF.db --id LFLAK_DWL_CRN --out $DbDir/CRN_ORF_merged.db --gff > $CRN_Merged_Gff
# Final results are reported:
echo "Number of CRN ORFs after merging:"
cat $CRN_Merged_Gff | grep 'gene' | wc -l
done
```

```
Searching for LFLAK domains in: P.cactorum 414_v2
Initial search space (Z):             548526  [actual number of targets]
Domain search space  (domZ):             342  [number of targets reported over threshold]
The number of CRNs common to both models are:
223
Number of CRN ORFs after merging:
155
```


Extract crinklers from published gene models


```bash
  for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/*/* | grep '414_v2'); do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    AugGff=$(ls gene_pred/final_genes/$Species/$Strain/final/final_genes_appended.gff3)
    AugCRN=$(ls $MergeDir/"$Strain"_pub_CRN_LFLAK_DWL.gff)
    AugFa=$(ls gene_pred/final_genes/"$Species"/"$Strain"/final/final_genes_combined.pep.fasta)
    ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
    ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer.gff3
    ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_CRN_hmmer.gff
    AugInORFs=$MergeDir/"$Strain"_AugInORFs_CRN_hmmer.gff
    ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_CRN_hmmer.gff
    # RxLR ORFs may be located in predicted gene models not identified as RxLR
    # proteins - manual inspection suggests that many of the original gene models
    # appear have RNAseq data aligning accross their entire length
    ORFsUniqNoGenes=$MergeDir/"$Strain"_ORFsUniq_excl_genes_CRN_hmmer.gff
    AugUniq=$MergeDir/"$Strain"_Aug_Uniq_CRN_hmmer.bed
    TotalCRNsTxt=$MergeDir/"$Strain"_final_CRN.txt
    TotalCRNsGff=$MergeDir/"$Strain"_final_CRN.gff
    TotalCRNsHeaders=$MergeDir/"$Strain"_Total_CRN_headers.txt
    bedtools intersect -wa -u -a $ORFGff -b $AugCRN > $ORFsInAug
    bedtools intersect -wa -u -a $AugCRN -b $ORFGff > $AugInORFs
    bedtools intersect -v -wa -a $ORFGff -b $AugCRN > $ORFsUniq
    bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniqNoGenes
    bedtools intersect -v -wa -a $AugCRN -b $ORFGff > $AugUniq
    echo "$Species - $Strain"

    echo "The number of ORF CRNs overlapping Augustus CRNs:"
    cat $ORFsInAug | grep -w -e 'transcript' -e 'mRNA'  | sort | uniq | wc -l
    echo "The number of Augustus CRNs overlapping ORF CRNs:"
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA' | sort | uniq | wc -l
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '='  | sort | uniq > $TotalCRNsTxt
    echo "The number of CRNs unique to ORF models:"
    cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' | sort | uniq | wc -l
    # cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' | sort | uniq >> $TotalCRNsTxt
    echo "The number of RxLRs unique to ORF models and not intersecting non-CRN gene models:"
    cat $ORFsUniqNoGenes | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l   
    cat $ORFsUniqNoGenes | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' | sort | uniq >> $TotalCRNsTxt
    echo "The number of CRNs unique to Augustus models:"
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | sort | uniq | wc -l
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | sort | uniq >> $TotalCRNsTxt
    echo "The total number of crinklers is:"
    cat $TotalCRNsTxt | wc -l

    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '='  | sort | uniq > $MergeDir/"$Strain"_final_CRN_ID.txt
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | sort | uniq >> $MergeDir/"$Strain"_final_CRN_ID.txt
    # cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f2 -d ';' | cut -f2 -d '=' | sort | uniq >> $MergeDir/"$Strain"_final_CRN_ID.txt
    cat $ORFsUniqNoGenes | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f2 -d ';' | cut -f2 -d '=' | sort | uniq >> $MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm_ID.txt

    # cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalCRNsTxt > $TotalCRNsGff
    cat $AugInORFs $AugUniq $ORFsUniqNoGenes | grep -w -f $TotalCRNsTxt > $TotalCRNsGff

    CRNsFa=$MergeDir/"$Strain"_final_CRN.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalCRNsTxt > $CRNsFa
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalCRNsTxt >> $CRNsFa
    echo "The number of sequences extracted is"
    cat $CRNsFa | grep '>' | wc -l

  done
```

Multiple Augustus genes overlapped ORFs where alternative transcripts were predicted

```
P.cactorum - 414_v2
The number of ORF CRNs overlapping Augustus CRNs:
151
The number of Augustus CRNs overlapping ORF CRNs:
154
The number of CRNs unique to ORF models:
4
The number of RxLRs unique to ORF models and not intersecting nonRxLR gene models:
3
The number of CRNs unique to Augustus models:
18
The total number of crinklers is:
175
The number of sequences extracted is
175
```


## C.ii) SSCP

Small secreted cysteine rich proteins were identified within secretomes. These
proteins may be identified by EffectorP, but this approach allows direct control
over what constitutes a SSCP.

```bash

for Secretome in $(ls gene_pred/combined_sigP/*/*/*_secreted.fa | grep -v 'all_secreted' |  grep '414_v2'); do
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/sscp/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_all_results.fa
cat $OutDir/"$Strain"_sscp_all_results.fa | grep 'Yes' > $OutDir/"$Strain"_sscp.fa
printf "number of SSC-rich genes:\t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' | cut -f1 -d '.' | sort | uniq | wc -l
printf "Number of effectors predicted by EffectorP:\t"
EffectorP=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP_secreted_headers.txt)
cat $EffectorP | wc -l
printf "Number of SSCPs predicted by both effectorP and this approach: \t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' > $OutDir/"$Strain"_sscp_headers.txt
cat $OutDir/"$Strain"_sscp_headers.txt $EffectorP | cut -f1 | sort | uniq -d | wc -l
echo ""
done > tmp.txt
cat tmp.txt
```


```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep '414_v2'); do
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"

OutDir=analysis/sscp_ORF/$Organism/$Strain
mkdir -p $OutDir
printf "the number of SigP gene is:\t";
cat $Secretome | grep '>' | cut -f1 | sort | uniq | wc -l
printf "the number of SSCP genes are:\t";
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_ORF.fa

cat $OutDir/"$Strain"_sscp_ORF.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/"$Strain"_sscp_ORF_headers_unmerged.txt
cat $OutDir/"$Strain"_sscp_ORF_headers_unmerged.txt | tr -d ' ' | sort | uniq | wc -l

SigP_Gff=$(ls gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_sscp_ORF_headers_unmerged.txt  $SigP_Gff	sscp_filter.py Name Augustus > $OutDir/"$Strain"_sscp_ORF_unmerged.gff

SSCP_Merged_Gff=$OutDir/"$Strain"_ORF_sscp_merged.gff
SSCP_Merged_txt=$OutDir/"$Strain"_ORF_sscp_merged.txt
SSCP_Merged_AA=$OutDir/"$Strain"_ORF_sscp_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_sscp_ORF_unmerged.gff --db sigP_ORF_sscp.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_sscp.db --id sigP_ORF_sscp --out sigP_ORF_sscp_merged.db --gff > $SSCP_Merged_Gff
cat $SSCP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $SSCP_Merged_txt

ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $SSCP_Merged_txt > $SSCP_Merged_AA
printf "Merged SSCP proteins:\t"
cat $SSCP_Merged_AA | grep '>' | wc -l
printf "\n"
done
```



# Making a combined file of Braker, Coding quary genes as well as additional ORF effector candidates


A gff file containing the combined Braker and CodingQuary genes as well as the
additional CRN and RxLR genes predicted by ORF analysis was made.

```bash
GeneGff=$(ls gene_pred/final_genes/P.cactorum/414_v2/final/final_genes_appended.gff3)
GffOrfRxLR=$(ls analysis/RxLR_effectors/combined_evidence/P.cactorum/414_v2/414_v2_ORFsUniq_excl_genes_RxLR_EER_motif_hmm.gff)
GffOrfCRN=$(ls analysis/CRN_effectors/hmmer_CRN/P.cactorum/414_v2/414_v2_ORFsUniq_excl_genes_CRN_hmmer.gff)
Assembly=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_softmasked_repeatmasker_TPSI_appended.fa )
OutDir=gene_pred/annotation/P.cactorum/414_v2
mkdir -p $OutDir
cat $GeneGff > $OutDir/414_v2_genes_incl_ORFeffectors.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/add_ORF_features.pl $GffOrfRxLR $Assembly >> $OutDir/414_v2_genes_incl_ORFeffectors.gff3
$ProgDir/add_ORF_features.pl $GffOrfCRN $Assembly >> $OutDir/414_v2_genes_incl_ORFeffectors.gff3
# Make gene models from gff files.
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
Assembly=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_softmasked_repeatmasker_TPSI_appended.fa)
$ProgDir/gff2fasta.pl $Assembly $OutDir/414_v2_genes_incl_ORFeffectors.gff3 $OutDir/414_v2_genes_incl_ORFeffectors
```


In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed
<!-- - CUFF_969_2_184 was identified as a duplicated gene -->


```bash
GffAppended=$(ls gene_pred/annotation/P.cactorum/414_v2/414_v2_genes_incl_ORFeffectors.gff3)
OutDir=gene_pred/final_ncbi/P.cactorum/414_v2/final_ncbi
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/remove_dup_features.py --inp_gff $GffAppended
# cat $GffAppended | grep -v -w 'CUFF_969_2_184' > $OutDir/414_v2_genes_incl_ORFeffectors_filtered.gff3
cat $GffAppended > $OutDir/414_v2_genes_incl_ORFeffectors_filtered.gff3

GffRenamed=$OutDir/414_v2_genes_incl_ORFeffectors_renamed.gff3
RenameLog=$OutDir/414_v2_genes_incl_ORFeffectors_renamed.log
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/gff_rename_genes.py --inp_gff $OutDir/414_v2_genes_incl_ORFeffectors_filtered.gff3 --conversion_log $RenameLog > $GffRenamed

Assembly=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_softmasked_repeatmasker_TPSI_appended.fa )
$ProgDir/gff2fasta.pl $Assembly $GffRenamed $OutDir/414_v2_genes_incl_ORFeffectors_renamed

# The proteins fasta file contains * instead of Xs for stop codons, these should
# be changed
sed -i 's/\*/X/g' $OutDir/414_v2_genes_incl_ORFeffectors_renamed.pep.fasta
```


# Re-running Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
  # mkdir -p gene_pred/final_genes/P.cactorum/414_v2/final_ncbi
  # gene_pred/final_ncbi/P.cactorum/414_v2/414_v2_genes_incl_ORFeffectors_* gene_pred/final_genes/P.cactorum/414_v2/final_ncbi/.
	screen -a
	cd /home/groups/harrisonlab/project_files/idris
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/final_ncbi/*/*/final_ncbi/414_v2_genes_incl_ORFeffectors_renamed.pep.fasta | grep '414_v2'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteome in $(ls gene_pred/final_ncbi/*/*/final_ncbi/414_v2_genes_incl_ORFeffectors_renamed.pep.fasta | grep '414_v2'); do
Strain=$(echo $Proteome | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
echo $Strain
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteome $InterProRaw
done
```


## B) SwissProt

```bash
  for Proteome in $(ls gene_pred/final_ncbi/*/*/final_ncbi/414_v2_genes_incl_ORFeffectors_renamed.pep.fasta | grep '414_v2'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/swissprot/$Organism/$Strain
    SwissDbDir=../../uniprot/swissprot
    SwissDbName=uniprot_sprot
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
    qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
  done
```


# Quantifying expression of predicted genes

## Quantification of gene models


Maria performed alignment of RNAseq data vs the P.cac genome, removing any reads
that aligned to the F. annanassa transcriptome.



Quantification of these genes was performed using featureCounts program as part
of the Subreads package.

```bash
  Gff=$(ls gene_pred/final_ncbi/P.cactorum/414_v2/final_ncbi/414_v2_genes_incl_ORFeffectors_renamed.gff3 | grep '414_v2')
  # for BamFile in $(ls alignment/star/P.cactorum/414_v2/*/*/star_aligmentAligned.sortedByCoord.out.bam); do
  for BamFile in $(ls ../../../..//sobczm/popgen/rnaseq/vesca_*/pcac/star_aligmentAligned.sortedByCoord.out.bam); do
    OutDir=$(dirname $BamFile)
    OutDir=alignment/star/P.cactorum/414_v2/featureCounts_maria_alignment
    Prefix=$(echo $BamFile | rev | cut -f3 -d '/' | rev | sed 's/vesca_//g')
    Jobs=$(qstat | grep 'sub_fea' | grep 'qw'| wc -l)
    while [ $Jobs -gt 1 ]; do
      sleep 1m
      printf "."
      Jobs=$(qstat | grep 'sub_fea' | grep 'qw'| wc -l)
    done
    printf "\n"
    echo $Prefix
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_featureCounts.sh $BamFile $Gff $OutDir $Prefix
  done
```

A file was created with columns referring to experimental treatments:

```bash
  OutDir=alignment/star/P.cactorum/414_v2/DeSeq
  mkdir -p $OutDir
  # make file in excel and copy accross
  # Parse the file if it was made in windows:
  #cat $OutDir/P.cactorum_RNAseq_design.txt | tr -d '\r' | sed 's/Timepoint/Timepoint\n/g' | sed "s/hours/hours\n/g" > $OutDir/P.cactorum_RNAseq_design_parsed.txt
  # Parse the file and remove technical replicate information:
  cat $OutDir/P.cactorum_RNAseq_design.txt | tr -d '\r' | sed 's/Timepoint/Timepoint\n/g' | sed "s/hours/hours\n/g" | sed 's/_L...//g' | cut -f1-4,6- | uniq | grep -v '\-\-' > $OutDir/P.cactorum_RNAseq_design_parsed.txt

  # Edit header lines of feature coutn files to ensure they have the treament name rather than file name
  OutDir=alignment/star/P.cactorum/414_v2/DeSeq
  for File in $(ls alignment/star/P.cactorum/414_v2/featureCounts_maria_alignment/*_featurecounts.txt); do
    echo $File;
    cp $File $OutDir/.;
  done
  for File in $(ls $OutDir/*_featurecounts.txt); do
    Prefix=$(echo $File | rev | cut -f1 -d '/' | rev | sed 's/_featurecounts.txt//g' | sed "s/_totRNA_S.*_L/_L/g")
    sed -ie "s/star_aligmentAligned.sortedByCoord.out.bam/$Prefix/g" $File
  done
```

DeSeq commands as used in R are documented in:
RNAseq/P414/DESeq_analysis.md


# Building summary tables of all data


For P. cactorum

```bash
for GeneGff in $(ls gene_pred/final_ncbi/P.cactorum/414_v2/final_ncbi/414_v2_genes_incl_ORFeffectors_renamed.gff3 | grep '414_v2'); do
  Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
  Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked.fa)
  InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
  SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
  OutDir=gene_pred/annotation/$Organism/$Strain
  mkdir -p $OutDir
	# GeneFasta=$(ls gene_pred/annotation/P.cactorum/414_v2/414_v2_genes_incl_ORFeffectors.pep.fasta)
  GeneFasta=$(ls gene_pred/final_ncbi/P.cactorum/414_v2/final_ncbi/414_v2_genes_incl_ORFeffectors_renamed.cds.fasta)
	SigP2=$(ls gene_pred/final_sigP/$Organism/$Strain/*_aug_sp.aa)
	SigP4=$(ls gene_pred/final_signalp-4.1/$Organism/$Strain/*_aug_sp.aa)
  TMHMM_headers=$(ls gene_pred/trans_mem/$Organism/$Strain/*_TM_genes_pos_headers.txt)
  GPI_headers=$(ls gene_pred/trans_mem/$Organism/$Strain/GPIsom/GPI_pos.txt)
	PhobiusTxt=$(ls analysis/phobius/$Organism/$Strain/*_phobius_headers.txt)
	#RxLR_Motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/*_RxLR_EER_regex.fa | grep -v 'ORF')
	#RxLR_Hmm=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer.fa | grep -v 'ORF')
	#RxLR_WY=$(ls analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain/*_WY_hmmer_headers.txt | grep -v 'ORF')
  RxLR_total=$(ls analysis/RxLR_effectors/combined_evidence/$Organism/$Strain/*_Total_RxLR_EER_motif_hmm_ID.txt)
	#CRN_LFLAK=$(ls analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/*_pub_CRN_LFLAK_hmm.fa | grep -v 'ORF')
	#CRN_DWL=$(ls analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/*_pub_CRN_DWL_hmm.fa | grep -v 'ORF')
  CRN_total=$(ls analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/*_final_CRN_ID.txt)
#	OrthoName=Pcac
#	OrthoFile=$(ls analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_Pinf_Ppar_Pcap_Psoj_orthogroups.txt)
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/gene_annotation
  DEG_Files=$(ls alignment/star/P.cactorum/414_v2/DeSeq/*_vs_*.txt  | grep -v -e 'up' -e 'down' | sed -e "s/$/ /g" | tr -d "\n")
	# $ProgDir/pacbio_anntoation_tables.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP2 $SigP2 --SigP4 $SigP4 --phobius $PhobiusTxt --RxLR_motif $RxLR_Motif --RxLR_Hmm $RxLR_Hmm --RxLR_WY $RxLR_WY --RxLR_total $RxLR_total --CRN_LFLAK $CRN_LFLAK --CRN_DWL $CRN_DWL --CRN_total $CRN_total --DEG_files $DEG_Files  > $OutDir/414_v2_gene_table_incl_exp.tsv
  # NormCount=$(ls alignment/star/P.cactorum/414_v2/DeSeq/normalised_counts.txt)
  RawCount=$(ls alignment/star/P.cactorum/414_v2/DeSeq/raw_counts.txt)
  FPKM=$(ls alignment/star/P.cactorum/414_v2/DeSeq/fpkm_counts.txt)
  # File showing gene name conversions
  ConversionLog=$(ls gene_pred/final_ncbi/P.cactorum/414_v2/final_ncbi/414_v2_genes_incl_ORFeffectors_renamed.log)
  VcfFiles=$(ls analysis/popgen/SNP_calling/*_vs_P414/*_no_indels.recode_gene.vcf)
  $ProgDir/pacbio_anntoation_tables.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --SigP2 $SigP2 --SigP4 $SigP4 --phobius $PhobiusTxt --trans_mem $TMHMM_headers --GPI_anchor $GPI_headers --RxLR_total $RxLR_total --CRN_total $CRN_total --DEG_files $DEG_Files --raw_counts $RawCount --fpkm $FPKM --InterPro $InterPro --Swissprot $SwissProt --SNP $VcfFiles --gene_conversion $ConversionLog > $OutDir/414_v2_gene_table_incl_exp.tsv
done
```
