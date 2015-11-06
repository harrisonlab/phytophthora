
Braker_prediction.md

This document details commands used for an initial run of Braker to train and
predict gene models.

Braker is a pipeline based upon Augustus and GeneMark-est.



## Align RNAseq reads to the genome:

RNAseq data from pre-infection life stages of P. cactorum were used to train gene
models to 10300. These data were from NCBI bioproject PRJNA242795

### 1) QC

Perform qc of RNAseq timecourse data. These reads are not actually paired reads
but this is irrelivent for processing usin fast-mcf.

```bash
  FileF=raw_rna/genbank/P.cactorum/10300/SRR1206032.fastq
  FileR=raw_rna/genbank/P.cactorum/10300/SRR1206033.fastq
  IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
  qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
```

## 2) Align reads vs. Fus2 genome
Alignments of RNAseq reads were made against the Fus2 Genome using tophat:

## 2.1) Alignment

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
Genome=assembly/abyss/P.cactorum/10300/10300_abyss_53/10300_abyss-scaffolds_500bp_renamed.fa
FileF=qc_rna/raw_rna/genbank/P.cactorum/F/SRR1206032_trim.fq.gz
FileR=qc_rna/raw_rna/genbank/P.cactorum/R/SRR1206033_trim.fq.gz
OutDir=alignment/P.cactorum/10300
qsub $ProgDir/tophat_alignment.sh $Genome $FileF $FileR $OutDir
```


# Run Braker1

screen -a
set variables
```bash
  qlogin
  WorkDir=/tmp/braker
  ProjDir=/home/groups/harrisonlab/project_files/idris
  Assembly=$ProjDir/assembly/abyss/P.cactorum/10300/10300_abyss_53/10300_abyss-scaffolds_500bp_renamed.fa
  OutDir=$ProjDir/gene_pred/braker/P.cactorum/10300
```

move to working directory
```
mkdir -p $WorkDir
cd $WorkDir
```


```bash
  braker.pl \
    --cores 16 \
    --genome=$Assembly \
    --GENEMARK_PATH=/home/armita/prog/genemark/gm_et_linux_64/gmes_petap \
    --BAMTOOLS_PATH=/home/armita/prog/bamtools/bamtools/bin \
    --species=P.cactorum \
    --bam=$ProjDir/alignment/P.cactorum/10300/accepted_hits.bam
```

```bash
  mkdir -p $OutDir
  cp -r braker/* $OutDir/.

  rm -r $WorkDir
```
