
## Data extraction

```bash
  cd /home/groups/harrisonlab/project_files/idris
  RawDatDir=/home/harrir/projects/pacbio_test/p_cact
  mkdir -p raw_dna/pacbio/P.cactorum/414
  cp -r $RawDatDir/A07_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/B07_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/G07_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/H07_1 raw_dna/pacbio/P.cactorum/414/.
  OutDir=raw_dna/pacbio/P.cactorum/414/extracted
  mkdir -p $OutDir
  cat raw_dna/pacbio/P.cactorum/414/*/Analysis_Results/*.subreads.fastq > $OutDir/concatenated_pacbio.fastq

```

## Assembly


### Canu assembly

```bash
  Reads=$(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq)
  GenomeSz="65m"
  Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir="assembly/canu/$Organism/$Strain"
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
```


### Spades Assembly

```bash
for PacBioDat in $(ls raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/Fus2_pacbio.fastq); do
Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/s_6_1_sequence_trim.fq.gz);
TrimR1_Read=$(ls $IlluminaDir/R/s_6_2_sequence_trim.fq.gz);
TrimF2_Read=$(ls $IlluminaDir/F/FUS2_S2_L001_R1_001_trim.fq.gz);
TrimR2_Read=$(ls $IlluminaDir/R/FUS2_S2_L001_R2_001_trim.fq.gz);
OutDir=assembly/spades_pacbio/$Organism/"$Strain"
echo $TrimR1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
qsub $ProgDir/subSpades_2lib_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir 50
done
```
