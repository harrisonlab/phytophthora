
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

## 2) Align reads vs. publihsed genomes
Alignments of RNAseq reads were made against the published Genomes using tophat:

## 2.1) Alignment

```bash
  Pinf=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
  Ppar=assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra_310.i2.scaffolds.genome.parsed.fa
  Pcap=assembly/external_group/P.capsici/LT1534/dna/Phyca11_unmasked_genomic_scaffolds.fasta
  Psoj=assembly/external_group/P.sojae/67593/dna/Phytophthora_sojae.ASM14975v1.26.dna.genome.parsed.fa
  for Genome in $Pinf $Ppar $Pcap $Psoj; do
    Strain=$(echo $Genome| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Genome | rev | cut -d '/' -f4 | rev)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    FileF=qc_rna/raw_rna/genbank/P.cactorum/F/SRR1206032_trim.fq.gz
    FileR=qc_rna/raw_rna/genbank/P.cactorum/R/SRR1206033_trim.fq.gz
    OutDir=alignment/$Organism/$Strain
    qsub $ProgDir/tophat_alignment.sh $Genome $FileF $FileR $OutDir
  done
```


# Run Braker1

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
  for Assembly in $Pinf $Ppar $Pcap $Psoj; do
    # Jobs=$(qstat | grep 'tophat_ali' | wc -l)
    # while [ $Jobs -gt 0 ]; do
    #   sleep 10
    #   printf "."
    #   Jobs=$(qstat | grep 'tophat_ali' | wc -l)
    # done
    printf "\n"
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    OutDir=gene_pred/braker/$Organism/$Strain
    AcceptedHits=alignment/$Organism/$Strain/accepted_hits.bam
    GeneModelName="$Organism"_"$Strain"_braker
    qsub $ProgDir/sub_braker.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```
