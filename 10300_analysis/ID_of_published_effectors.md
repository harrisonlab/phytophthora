# ID of published effectors

This file documents the identification of previously reported RxLR and Crinkler
effectors in the 10300 genome.


The documentation details commands used to identify
* RxLR effectors from the 10300 transcriptome - Chen et al 2014
* CRN effectors from the 10300 transcriptome - Chen et al 2014
* RxLR effectors from the P. infestans, P. parisitica, P. capsici and P. sojae genomes
* Previously characterised AVR genes


## 1 RxLR effectors from the 10300 transcriptome - Chen et al 2014

### 1.1 preperation of query fasta
The 96 putative RxLR proteins identified within the 10300 transcriptome were
downloaded in an excel table from:

```
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4289400/
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4289400/bin/12864_2014_6857_MOESM9_ESM.xls
```

This table was converted to comma seperated format and uploaded to the cluster at
the following location:
```bash
  cd /home/groups/harrisonlab/project_files/idris
  mkdir -p analysis/blast_homology/oomycete_avr_genes/
  # File saved as: analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.csv
  sed -i analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.csv 's/^M/\n/'g
```

A query file was built from the RxLRs identified from transcriptome sequencing
in Chen et al 2014.
```bash
	cat analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.csv | grep -v 'Additional' | grep -v 'Index' | grep -v -e '^,' | cut -f 2,5 -d ',' | sed -r 's/^/>/g' | sed 's/,/\n/g' > analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.fa
```

### 1.2 Performing BLAST searches

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Assembly=repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
	Query=analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.fa
	qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	BlastHits=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.csv
	HitsGff=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.gff
	Column2=Chen_RxLR
	NumHits=1
	$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
```

### 1.3 Characterising BLAST hits

Genes overlapped by BLAST hits were identified:

```bash
HitsGff=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.gff
Proteins=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.gff
ORFs=gene_pred/ORF_finder/P.cactorum/10300/10300_ORF_corrected.gff3
OutDir=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_analysis
mkdir -p $OutDir
BrakerIntersect=$OutDir/10300_chen_et_al_2014_RxLR_BrakerIntersect.gff
BrakerNoIntersect=$OutDir/10300_chen_et_al_2014_RxLR_BrakerNoIntersect.gff
ORFIntersect=$OutDir/10300_chen_et_al_2014_RxLR_ORFIntersect.gff

RxLRs=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_Total_RxLR_EER_motif_hmm_headers.txt
AugGff=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_Aug_RxLR_EER_motif_hmm.gff
ORFGff=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_ORF_RxLR_EER_motif_hmm.gff
ChenMissingRxLRs=$OutDir/10300_chen_et_al_2014_RxLR_MissingRxLRs.gff
ChenSupportedRxLRs=$OutDir/10300_chen_et_al_2014_RxLR_SupportedRxLRs.gff

echo "The following number of blast hits intersected Braker gene models:"
bedtools intersect -wa -u -a $HitsGff -b $Proteins > $BrakerIntersect
cat $BrakerIntersect | wc -l
echo "The following Blast hits did not intersect Braker gene models:"
bedtools intersect -v -a $HitsGff -b $Proteins > $BrakerNoIntersect
cat $BrakerNoIntersect | wc -l
echo "The following number of blast hits intersected ORF gene models:"
bedtools intersect -wa -u -a $BrakerNoIntersect -b $ORFs > $ORFIntersect
cat $ORFIntersect | wc -l
echo "The following Chen et al 2014 RxLRs were not found in this study:"
bedtools intersect -v -a $HitsGff -b $AugGff $ORFGff > $ChenMissingRxLRs
cat $ChenMissingRxLRs | wc -l
echo "The following Chen et al 2014 RxLRs were also found in this study:"
bedtools intersect -wa -u -a $HitsGff -b $AugGff $ORFGff > $ChenSupportedRxLRs
cat $ChenSupportedRxLRs | wc -l
```
Of the 94 queries all 77 showed overlap to genes predicted in Braker1 gene
models and the remaining 17 overlapped predicted ORFs. 47 of these 94 genes were
predicted as RxLRs in this study. The remaining 47 query genes that did not
overlap genes identified as RxLRs in this study were manually inspected. It was
noted that of the 22 of the missing genes were identified in Chen et al as lower
confidence RxLRs as they lacked RxLR and EER motifs and were not complete
transcripts. Of the missing Chen et al RxLRs, 3 lacked the RxLR motif completely and 10
had non-typical RxLR sequences. This was also true for the EER motif, with 4 of
the missing RxLRs lacking the EER motif completely and 8 possessing an EER motif
that did not conform to the [ED]+[ED][KR] motif used in this study.


## 4 Previously characterised AVR genes

## 4.2 Performing BLAST searches

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  Assembly=repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
  Query=analysis/blast_homology/oomycete_avr_genes/appended_oomycete_avr_cds.fasta
  qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  BlastHits=analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_avr_cds.fasta_homologs.csv
  HitsGff=analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_avr_cds.fasta_homologs.gff
  Column2=Avr_homolog
  NumHits=1
  $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
```


-->
