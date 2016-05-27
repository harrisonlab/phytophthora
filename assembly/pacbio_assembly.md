
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
  cat raw_dna/pacbio/P.cactorum/414/*/Analysis_Results/*.subreads.fastq > $OutDir/concatenated_pacbio.fastq
```
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

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/canu/*/*/*_canu.contigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/canu/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

### Spades Assembly

<!-- For P. fragariae

```bash
for PacBioDat in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/Bc16_S1_L001_R1_001_trim.fq.gz);
TrimR1_Read=$(ls $IlluminaDir/R/Bc16_S1_L001_R2_001_trim.fq.gz);
TrimF2_Read=$(ls $IlluminaDir/F/Bc16_S2_L001_R1_001_160129_trim.fq.gz);
TrimR2_Read=$(ls $IlluminaDir/R/Bc16_S2_L001_R2_001_160129_trim.fq.gz);
OutDir=assembly/spades_pacbio/$Organism/"$Strain"
echo $TrimR1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
qsub $ProgDir/subSpades_2lib_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir 50
done
``` -->

For P. cactorum

```bash
  for PacBioDat in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
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
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    OutDir=assembly/spades/$Organism/$Strain
    qsub $ProgDir/subSpades_2lib_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir 28
  done
```

Contigs shorter than 500bp were renomed from the assembly

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



Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/spades_pacbio/*/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/canu/*/*/*_canu.contigs.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/"$Strain"_canu_contigs_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```



# Analysis of preliminary assemblies

## Preparing data


```bash
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir

  Fus2_pacbio_canu=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_canu/test/Fus2_pacbio_canu.fa
  Fus2_pacbio_spades=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_spades/test/Fus2_pacbio_spades.fa
  Fus2_pacbio_merged=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/test/Fus2_pacbio_merged.fa

  mkdir -p $(dirname $Fus2_pacbio_canu)
  mkdir -p $(dirname $Fus2_pacbio_spades)
  mkdir -p $(dirname $Fus2_pacbio_merged)
  cp /home/harrir/projects/pacbio_test/fus2/fus2-auto/fus2.contigs.fasta $Fus2_pacbio_canu
  cp /home/harrir/projects/pacbio_test/spades/FUS2/scaffolds.fasta $Fus2_pacbio_spades
  cp /home/harrir/projects/pacbio_test/hybrid_merge/fus2/merged.fasta $Fus2_pacbio_merged
```

## Repeatmasking

```bash
Fus2_pacbio_canu=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_canu/test/Fus2_pacbio_canu.fa
Fus2_pacbio_spades=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_spades/test/Fus2_pacbio_spades.fa
Fus2_pacbio_merged=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/test/Fus2_pacbio_merged.fa
for Assembly in $(ls $Fus2_pacbio_spades $Fus2_pacbio_merged $Fus2_pacbio_canu); do
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  qsub $ProgDir/rep_modeling.sh $Assembly
  qsub $ProgDir/transposonPSI.sh $Assembly
done
```

## Blast searches of LS region genes vs FoC

Some preliminary commands were used to analyse Pacbio assemblies Richard had generated

Blast searches were performed against these assembled contigs to identify which
contigs contained blast homologs from known FoL LS genes.

Headers of LS genes from FoL had previously been extracted. These were used to
extract the relevant proteins from fasta files.

```bash
  ProtFastaUnparsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta
  ProtFasta=analysis/FoL_ls_genes/4287_proteins_renamed.fasta
  cat $ProtFastaUnparsed | sed -r "s/^>.*FOXG/>FOXG/g" > $ProtFasta
  for File in $(ls analysis/FoL_ls_genes/chr_*_gene_headers.txt); do
    Chr=$(echo $File | rev | cut -f3 -d '_' | rev);
    echo $File;
    echo "extracting proteins associated with chromosome: $Chr";
    OutFile=$(echo $File | sed 's/_headers.txt/.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $ProtFasta --headers $File | grep -v -P '^$' > $OutFile
  done
```

```bash
  Fus2_pacbio_canu=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_canu/test/Fus2_pacbio_canu.fa
  Fus2_pacbio_spades=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_spades/test/Fus2_pacbio_spades.fa
  Fus2_pacbio_merged=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged/test/Fus2_pacbio_merged.fa
  for Genome in $(ls $Fus2_pacbio_spades $Fus2_pacbio_merged $Fus2_pacbio_canu); do
    for Proteome in $(ls analysis/FoL_ls_genes/chr_*_gene.aa); do
      Strain=$(echo $Genome | rev | cut -f3 -d '/' | rev)
      Organism=$(echo $Genome | rev | cut -f4 -d '/' | rev)
      Chr=$(echo $Proteome | rev | cut -f2 -d '_' | rev);
      echo "$Organism - $Strain - $Chr"
      OutDir=analysis/blast_homology/$Organism/$Strain
      ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
      qsub $ProgDir/run_blast2csv.sh $Proteome protein $Genome $OutDir
    done
  done
```

Convert top blast hits into gff annotations

```bash
  for BlastHitsCsv in $(ls analysis/blast_homology/*/Fus2_pacbio_*/*_chr_*_gene.aa_hits.csv); do
    Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
    Chr=$(echo $BlastHitsCsv | rev | cut -f3 -d '_' | rev);
    Column2=Chr"$Chr"_gene_homolog
    NumHits=1
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
  done
```
