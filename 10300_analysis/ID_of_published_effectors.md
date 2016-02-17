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
  cat analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.csv | grep -v 'Additional' | grep -v 'Index' | grep -v -e '^,' | cut -f 2,7 -d ',' | sed -r 's/^/>/g' | sed 's/,/\n/g' > analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR_prots.fa
```

### 1.2 Performing BLAST searches

#### 1.2.a)

The chen et al 2014 RxLRs were BLASTed against the 10300 genome to identify
locations of RxLR genes that may not have been detected in this study.

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

Bedtools was used to identify which of the Chen et al RxLRs intersected with
RxLRs predicted in our study:

```bash
  HitsGff=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.gff
  RxLRGff=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_Total_RxLR_EER_motif_hmm.gff
  OutDir=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR
  ChenMissingRxLRs=$OutDir/10300_chen_et_al_2014_RxLR_MissingRxLRs_2.gff
  echo "The following Chen et al 2014 RxLRs were not found in this study:"
  bedtools intersect -v -s -a $HitsGff -b $RxLRGff > $ChenMissingRxLRs
  cat $ChenMissingRxLRs | wc -l
```

The genes that missing RxLRs intersected were identified:

```bash
  AugGff=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus_extracted.gff
  AdditionalRxLRs=$OutDir/10300_chen_et_al_2014_RxLR_additional_genes.gff
  ChenNoAugRxLRs=$OutDir/10300_chen_et_al_2014_RxLR_noAug.gff
  echo "The following additional genes have been identified as RxLRs through BLAST homology to Chen et al 2014 RxLRs:"
  bedtools intersect -s -a $AugGff -b $ChenMissingRxLRs > $AdditionalRxLRs
  cat $AdditionalRxLRs | grep -w 'gene' | wc -l
  echo "The following Chen et al 2014 RxLRs have still not been identified:"
  bedtools intersect -v -s -a $ChenMissingRxLRs  -b $AugGff > $ChenNoAugRxLRs
  cat $ChenNoAugRxLRs | wc -l
```

This identified an additional 17 genes, leaving 30 remaining.

Finally, ORFs intersecting Chen et al RxLRs were identified:

```bash
  OrfSigPGff=gene_pred/ORF_sigP/P.cactorum/10300/10300_ORF_sp_merged.gff
  ChenNoAugRxLRs=$OutDir/10300_chen_et_al_2014_RxLR_noAug.gff
  AdditionalRxLRs=$OutDir/10300_chen_et_al_2014_RxLR_additional_genes.gff
  ChenNoSigPORFRxLRs=$OutDir/10300_chen_et_al_2014_RxLR_noSigPORF.gff
  echo "The following additional genes have been identified as RxLRs through BLAST homology to Chen et al 2014 RxLRs:"
  bedtools intersect -s -wo -a $OrfSigPGff -b $ChenNoAugRxLRs >> $AdditionalRxLRs
  cat $AdditionalRxLRs | grep -w 'gene' | wc -l
  echo "The following Chen et al 2014 RxLRs have still not been identified:"
  bedtools intersect -v -s -a $ChenNoAugRxLRs  -b $OrfSigPGff > $ChenNoSigPORFRxLRs
  cat $ChenNoSigPORFRxLRs | wc -l
```

This identified an additional 4 genes, leaving 26 remaining.

<!-- ```bash
  OrfGff=gene_pred/ORF_finder/P.cactorum/10300/10300_ORF_corrected.gff3
  ChenNoSigPORFRxLRs=$OutDir/10300_chen_et_al_2014_RxLR_noSigPORF.gff
  AdditionalRxLRsORF=$OutDir/10300_chen_et_al_2014_RxLR_additional_ORFs.gff
  echo "The following additional genes have been identified as RxLRs through BLAST homology to Chen et al 2014 RxLRs:"
  bedtools intersect -s -a $OrfGff -b $ChenNoSigPORFRxLRs > $AdditionalRxLRsORF
  cat $AdditionalRxLRsORF | grep -w 'gene' | wc -l
  echo "The following Chen et al 2014 RxLRs have still not been identified:"
  bedtools intersect -v -s -a $ChenNoAugRxLRs  -b $OrfSigPGff > $ChenNoSigPORFRxLRs
  cat $ChenNoSigPORFRxLRs | wc -l
``` -->


#### 1.2.b)

Furthermore, 10300 predicted proteins and RxLRs were BLASTed against P. infestans
predicted proteins to identify homology to P.infestans RxLRs.

```bash
PinfProts=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
OutDir=analysis/blast_homology/oomycete_avr_genes
PinfRxLRheaders=$OutDir/P.infestans_published_RxLR_headers.txt
PinfRxLRs=$OutDir/P.infestans_published_RxLR.fa
mkdir -p $OutDir
cat $PinfProts | grep -i 'RxLR'  | cut -f1 -d ' ' | tr -d '>' > $PinfRxLRheaders
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $PinfProts --headers $PinfRxLRheaders > $PinfRxLRs
```

```bash
PredProts=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
PredRxLRs=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_Total_RxLR_EER_motif_hmm_headers.fa
PinfProts=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
DatabaseDir=analysis/databases/P.infestans/T30-4
OutDir=analysis/blast_homology/P.cactorum/10300

cat $Proteins | grep -i 'RxLR' | cut -f1 -d ' ' | tr -d '>' > $Pinf_RxLRs
makeblastdb -in $PinfProts -input_type fasta -dbtype prot -out $DatabaseDir/PinfProt.db
makeblastdb -in $PinfRxLRs -input_type fasta -dbtype prot -out $DatabaseDir/PinfRxLR.db

mkdir -p $OutDir
blastp \
-db $DatabaseDir/PinfProt.db \
-query $PredProts \
-out $OutDir/10300_prots_vs_Pinf_prots.tbl \
-evalue 1e-5 \
-outfmt 6 \
-num_threads 16 \
-num_alignments 1

# blastx \
# -db $OutDir/PinfRxLR.db \
# -query $PredProts \
# -out $OutDir/10300_prots_vs_Pinf_RxLRs.tbl \
# -evalue 1e-5 \
# -outfmt 6 \
# -num_threads 16 \
# -num_alignments 1

cat $OutDir/10300_prots_vs_Pinf_prots.tbl | grep -w -f $Pinf_RxLR_headers > $OutDir/10300_prots_vs_Pinf_prots_homologs.tbl
cat $OutDir/10300_prots_vs_Pinf_prots_homologs.tbl | wc -l
```

These results reported by Chen et al could not be replicated. Either when
performing blast searches against the proteome (67 chen et al RxLRs have hits)
or when performing searches against 486 proteins with RxLR putative functions
according to the published proteome (81 chen et al RxLRs have hits). The
Gene models may be different from those described in chen et al. The genes not
were investigated using the following commands:

```bash
cat $OutDir/10300_chen_et_al_2014_RxLR.fa_RxLR_homologs.tbl analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.fa | cut -f1 | grep -E -o 'PcRXLR\w+' | sort | uniq -u
cat $Proteins | grep 'PITG_23092'
```
It was noted that PcRXLR86 was absent, the corresponding gene identified as its
homolog in chen et al 2014 was identified within the P.infestans proteins and
was noted not have an RxLR function:
```
>PITG_23092T0 pep:known supercontig:ASM14294v1:supercont1.51:943348:943863:1 gene:PITG_23092 transcript:PITG_23092T0 description:"fucker"
```

#### 1.2.c)

To compare our results to Chen et al 2014. The RxLRs predicted in their
study were BLASTed against the P.infestans genome and Proteome as well. In their
study they report 89 of 94 of their RxLRs have BLAST homology to P. infestans
RxLRs. That was based on P. infestans gene models containing a total of 480
RxLRs.


This was also performed against the P.infestans genome:
```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Assembly=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
	Query=analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.fa
	qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
```
Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	BlastHits=analysis/blast_homology/P.infestans/T30-4/T30-4_chen_et_al_2014_RxLR.fa_homologs.csv
	HitsGff=analysis/blast_homology/P.infestans/T30-4/T30-4_chen_et_al_2014_RxLR.fa_homologs.gff
	Column2=Chen_RxLR
	NumHits=1
	$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
```

Blast searches were also performed against P. infestans predicted proteins.

```bash
Query=analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.fa
Proteins=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
BLASTProteins=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all_parsed.fa
cat $Proteins | cut -f1 -d ' ' > $BLASTProteins
OutDir=analysis/blast_homology/oomycete_avr_genes
Pinf_RxLR_headers=$OutDir/P.infestans_published_RxLR_headers.txt
Pinf_RxLR_fasta=$OutDir/P.infestans_published_RxLR_fasta.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Proteins --headers $Pinf_RxLR_headers > $Pinf_RxLR_fasta
cat $Proteins | grep -i 'RxLR' | cut -f1 -d ' ' | tr -d '>' > $Pinf_RxLRs
makeblastdb -in $Pinf_RxLR_fasta -input_type fasta -dbtype prot -out tmp.db
#Strain=$(echo $Proteins | rev | cut -f3 -d '/' | rev)
#Organism=$(echo $Proteins | rev | cut -f4 -d '/' | rev)
#OutDir=$ProjDir/gene_pred/augustus/$Species/$Strain/swissplot

mkdir -p $OutDir
blastx \
-db tmp.db \
-query $Query \
-out $OutDir/10300_chen_et_al_2014_RxLR.fa_protin_homologs.tbl \
-evalue 1e-5 \
-outfmt 6 \
-num_threads 16 \
-num_alignments 1
cat $OutDir/10300_chen_et_al_2014_RxLR.fa_protin_homologs.tbl | grep -w -f $Pinf_RxLR_headers > $OutDir/10300_chen_et_al_2014_RxLR.fa_RxLR_homologs.tbl
cat $OutDir/10300_chen_et_al_2014_RxLR.fa_RxLR_homologs.tbl | wc -l
```

These results reported by Chen et al could not be replicated. Either when
performing blast searches against the proteome (67 chen et al RxLRs have hits)
or when performing searches against 486 proteins with RxLR putative functions
according to the published proteome (81 chen et al RxLRs have hits). The
Gene models may be different from those described in chen et al. The genes not
were investigated using the following commands:

```bash
cat $OutDir/10300_chen_et_al_2014_RxLR.fa_RxLR_homologs.tbl analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.fa | cut -f1 | grep -E -o 'PcRXLR\w+' | sort | uniq -u
cat $Proteins | grep 'PITG_23092'
```
It was noted that PcRXLR86 was absent, the corresponding gene identified as its
homolog in chen et al 2014 was identified within the P.infestans proteins and
was noted not have an RxLR function:
```
>PITG_23092T0 pep:known supercontig:ASM14294v1:supercont1.51:943348:943863:1 gene:PITG_23092 transcript:PITG_23092T0 description:"fucker"
```

### 1.3 Characterising BLAST hits

Genes overlapped by BLAST hits were identified:

```bash
  PcacAugGff=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus_extracted.gff
  # PparPubGff=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_transcripts.gtf
  # PinfPubGff=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts.gtf
  # PinfPubGffParsed=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts_parsed.gtf
  # cat $PinfPubGff | sed 's/Supercontig_/supercont/g' > $PinfPubGffParsed
  # PcapPubGff=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_genes.gff
  # PsojPubGff=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_genes_20110401.gff

  # for Proteins in $PcacAugGff $PinfPubGffParsed $PparPubGff $PcapPubGff $PsojPubGff; do
  for Proteins in $PcacAugGff; do
    Species=$(echo "$Proteins" | rev | cut -f4 -d '/' | rev)
    Strain=$(echo "$Proteins" | rev | cut -f3 -d '/' | rev)
    HitsGff=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.gff
    ORFs=gene_pred/ORF_finder/P.cactorum/10300/10300_ORF_corrected.gff3
    OutDir=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_analysis
    mkdir -p $OutDir
    BrakerIntersect=$OutDir/10300_chen_et_al_2014_RxLR_BrakerIntersect.gff
    BrakerNoIntersect=$OutDir/10300_chen_et_al_2014_RxLR_BrakerNoIntersect.gff
    ORFIntersect=$OutDir/10300_chen_et_al_2014_RxLR_ORFIntersect.gff
    RxlrProteins=$OutDir/"$Strain"_chen_et_al_2014_RxLR_BrakerIntersect.txt
    RxlrORFs=$OutDir/"$Strain"_chen_et_al_2014_RxLR_ORFIntersect.txt

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
    if [ $Species == 'P.infestans' ] || [ $Species == 'P.parisitica' ]; then
    bedtools intersect -wao -a $HitsGff -b $Proteins | grep -w 'exon' | cut -f9,18 | cut -f1,3 -d ";" | tr -d '"' | tr -d ';' | sed 's/ID=//g' | sed 's/transcript_id //g' | sed 's/transcriptId //g' |sed 's/ /\t/g' > $RxlrProteins
    elif [ $Species == 'P.capsici' ] || [ $Species == 'P.sojae' ]; then
    bedtools intersect -wao -a $HitsGff -b $Proteins | grep -w 'exon' | cut -f9,18 | cut -f1,2 -d ";" | tr -d '"' | tr -d ';' | sed 's/ID=//g' | sed 's/name //g' | sed 's/ /\t/g' > $RxlrProteins
    else
      bedtools intersect -wao -a $HitsGff -b $Proteins | grep -w 'transcript' | cut -f9,18 | tr -d '"' | tr -d ';' | sed 's/ID=//g' > $RxlrProteins
    fi
    bedtools intersect -wao -a $BrakerNoIntersect -b $ORFs | grep -w 'transcript' | cut -f9,18 | cut -d ';' -f1,4 | tr -d '"' | sed 's/;/\t/g' | sed 's/ID=//g' | sed 's/Name=//g' > $RxlrORFs
  done
```

```
  The following number of blast hits intersected Braker gene models:
  77
  The following Blast hits did not intersect Braker gene models:
  17
  The following number of blast hits intersected ORF gene models:
  17
  The following Chen et al 2014 RxLRs were not found in this study:
  47
  The following Chen et al 2014 RxLRs were also found in this study:
  47
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



The genes relating to each RxLR Blast hit were extracted.

```bash
  Pcac_pep=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
  # Pinf_pep=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
  # Ppar_pep=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins.pep.all.fa
  # Pcap_pep=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins.fasta
  # Psoj_pep=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401.aa.fasta
  # for Proteins in $Pcac_pep $Pinf_pep $Ppar_pep $Pcap_pep $Psoj_pep; do
  for Proteins in $Pcac_pep; do
    Species=$(echo "$Proteins" | rev | cut -f4 -d '/' | rev)
    Strain=$(echo "$Proteins" | rev | cut -f3 -d '/' | rev)
    echo "$Species - $Strain"
    OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_chen_et_al_2014_analysis
    RxlrProteins=$(ls $OutDir/"$Strain"_chen_et_al_2014_RxLR_BrakerIntersect.txt)
    while read Line; do
      RxlrName=$(echo $Line | cut -f1 -d ' ' )
      HitName=$(echo $Line | cut -f2 -d ' ')
      echo "$RxlrName - $HitName"
      if [ $Species == 'P.cactorum' ] || [ $Species == 'P.infestans' ] || [ $Species == 'P.parisitica' ]; then
        echo "$HitName" > tmp.txt
      else
        cat $Proteins | grep -w "$HitName" | tr -d '>' > tmp.txt
      fi
      OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_chen_et_al_2014_RxLR/fasta
      mkdir -p $OutDir
      OutFasta=$(echo ""$Species"_"$Strain"_""$RxlrName""_homolog.fa" | sed 's/_BlastHit_1//g')
      ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
      $ProgDir/extract_from_fasta.py --fasta $Proteins --headers tmp.txt > $OutDir/$OutFasta
    done < $RxlrProteins
    OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_chen_et_al_2014_analysis
    echo ""
    ORFs=$(ls gene_pred/ORF_finder/$Species/$Strain/"$Strain".aa_cat.fa)
    RxlrORFs=$(ls $OutDir/"$Strain"_chen_et_al_2014_RxLR_ORFIntersect.txt )
    for RxlrName in $(cat $RxlrORFs | cut -f1 | sort | uniq); do
      cat $RxlrORFs | grep "$RxlrName" | cut -f2 > tmp.txt
      echo "$RxlrName - ORF hits"
      OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_chen_et_al_2014_RxLR_analysis/fasta
      mkdir -p $OutDir
      OutFasta=$(echo ""$Species"_"$Strain"_""$RxlrName""_homolog.fa" | sed 's/_BlastHit_1//g')
      ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
      $ProgDir/extract_from_fasta.py --fasta $ORFs --headers tmp.txt > $OutDir/$OutFasta
      sed -E "s/^>/>$Species|/g" -i $OutDir/$OutFasta
    done
    echo ""
  done
```




## 2 CRN effectors from the 10300 transcriptome - Chen et al 2014


### 2.1 preperation of query fasta
The 64 putative CRN proteins identified within the 10300 transcriptome were
downloaded in an excel table from:

```
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4289400/
http://static-content.springer.com/esm/art%3A10.1186%2F1471-2164-15-980/MediaObjects/12864_2014_6857_MOESM11_ESM.xls
```

This table was converted to comma seperated format and uploaded to the cluster at
the following location:
```bash
  cd /home/groups/harrisonlab/project_files/idris
  mkdir -p analysis/blast_homology/oomycete_avr_genes/
  # File saved as: analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_CRN.csv
  sed -i 's/^M/\n/g' analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_CRN.csv
  # Note you will need to edit the '^M' in the above text when you put it into the command line by pressing ctrl+v then ctrl+M
```

A query file was built from the CRNs identified from transcriptome sequencing
in Chen et al 2014.
```bash
	cat analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_CRN.csv | grep -v 'Additional' | grep -v 'Index' | grep -v -e '^,' | cut -f 2,5 -d ',' | sed -r 's/^/>/g' | sed 's/,/\n/g' > analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_CRN.fa
```

### 2.2 Performing BLAST searches

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Assembly=repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
	Query=analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_CRN.fa
	qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	BlastHits=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_CRN.fa_homologs.csv
	HitsGff=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_CRN.fa_homologs.gff
	Column2=Chen_CRN
	NumHits=1
	$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
```

### 2.3 Characterising BLAST hits

Genes overlapped by BLAST hits were identified:

```bash
  PcacAugGff=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus_extracted.gff
  # PparPubGff=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_transcripts.gtf
  # PinfPubGff=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts.gtf
  # PinfPubGffParsed=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts_parsed.gtf
  # cat $PinfPubGff | sed 's/Supercontig_/supercont/g' > $PinfPubGffParsed
  # PcapPubGff=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_genes.gff
  # PsojPubGff=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_genes_20110401.gff

  # for Proteins in $PcacAugGff $PinfPubGffParsed $PparPubGff $PcapPubGff $PsojPubGff; do
  for Proteins in $PcacAugGff; do
    Species=$(echo "$Proteins" | rev | cut -f4 -d '/' | rev)
    Strain=$(echo "$Proteins" | rev | cut -f3 -d '/' | rev)
    HitsGff=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_CRN.fa_homologs.gff
    ORFs=gene_pred/ORF_finder/P.cactorum/10300/10300_ORF_corrected.gff3
    OutDir=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_analysis
    mkdir -p $OutDir
    BrakerIntersect=$OutDir/10300_chen_et_al_2014_CRN_BrakerIntersect.gff
    BrakerNoIntersect=$OutDir/10300_chen_et_al_2014_CRN_BrakerNoIntersect.gff
    ORFIntersect=$OutDir/10300_chen_et_al_2014_CRN_ORFIntersect.gff    
    CrnProteins=$OutDir/"$Strain"_chen_et_al_2014_CRN_BrakerIntersect.txt
    CrnORFs=$OutDir/"$Strain"_chen_et_al_2014_CRN_ORFIntersect.txt

    AugGff=analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_Aug_CRN_hmmer.gff3
    ORFGff=analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORF_CRN_hmmer.gff3
    ChenMissingCRNs=$OutDir/10300_chen_et_al_2014_CRN_MissingCRNs.gff
    ChenSupportedCRNs=$OutDir/10300_chen_et_al_2014_CRN_SupportedCRNs.gff

    echo "The following number of blast hits intersected Braker gene models:"
    bedtools intersect -wa -u -a $HitsGff -b $Proteins > $BrakerIntersect
    cat $BrakerIntersect | wc -l
    echo "The following Blast hits did not intersect Braker gene models:"
    bedtools intersect -v -a $HitsGff -b $Proteins > $BrakerNoIntersect
    cat $BrakerNoIntersect | wc -l
    echo "The following number of blast hits intersected ORF gene models:"
    bedtools intersect -wa -u -a $BrakerNoIntersect -b $ORFs > $ORFIntersect
    cat $ORFIntersect | wc -l
    echo "The following Chen et al 2014 CRNs were not found in this study:"
    bedtools intersect -v -a $HitsGff -b $AugGff $ORFGff > $ChenMissingCRNs
    cat $ChenMissingCRNs | wc -l
    echo "The following Chen et al 2014 CRNs were also found in this study:"
    bedtools intersect -wa -u -a $HitsGff -b $AugGff $ORFGff > $ChenSupportedCRNs
    cat $ChenSupportedCRNs | wc -l
    if [ $Species == 'P.infestans' ] || [ $Species == 'P.parisitica' ]; then
      bedtools intersect -wao -a $HitsGff -b $Proteins | grep -w 'exon' | cut -f9,18 | cut -f1,3 -d ";" | tr -d '"' | tr -d ';' | sed 's/ID=//g' | sed 's/transcript_id //g' | sed 's/transcriptId //g' |sed 's/ /\t/g' > $CrnProteins
    elif [ $Species == 'P.capsici' ] || [ $Species == 'P.sojae' ]; then
      bedtools intersect -wao -a $HitsGff -b $Proteins | grep -w 'exon' | cut -f9,18 | cut -f1,2 -d ";" | tr -d '"' | tr -d ';' | sed 's/ID=//g' | sed 's/name //g' | sed 's/ /\t/g' > $CrnProteins
    else
      bedtools intersect -wao -a $HitsGff -b $Proteins | grep -w 'transcript' | cut -f9,18 | tr -d '"' | tr -d ';' | sed 's/ID=//g' > $CrnProteins
    fi
    bedtools intersect -wao -a $BrakerNoIntersect -b $ORFs | grep -w 'transcript' | cut -f9,18 | cut -d ';' -f1,4 | tr -d '"' | sed 's/;/\t/g' | sed 's/ID=//g' | sed 's/Name=//g' > $CrnORFs
  done
```

```
  The following number of blast hits intersected Braker gene models:
  43
  The following Blast hits did not intersect Braker gene models:
  21
  The following number of blast hits intersected ORF gene models:
  21
  The following Chen et al 2014 CRNs were not found in this study:
  30
  The following Chen et al 2014 CRNs were also found in this study:
  34
```

Of the 64 queries 43 showed overlap to genes predicted in Braker1 gene
models and the remaining 21 overlapped predicted ORFs. 34 of these 64 genes were
predicted as CRNs in this study. The remaining 30 query genes that did not
overlap genes identified as CRNs in this study were manually inspected.


The genes relating to each CRN BLast hit were extracted.

```bash
  Pcac_pep=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
  # Pinf_pep=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
  # Ppar_pep=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins.pep.all.fa
  # Pcap_pep=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins.fasta
  # Psoj_pep=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401.aa.fasta
  # for Proteins in $Pcac_pep $Pinf_pep $Ppar_pep $Pcap_pep $Psoj_pep; do
  for Proteins in $Pcac_pep; do
    Species=$(echo "$Proteins" | rev | cut -f4 -d '/' | rev)
    Strain=$(echo "$Proteins" | rev | cut -f3 -d '/' | rev)
    echo "$Species - $Strain"
    OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_chen_et_al_2014_analysis
    CrnProteins=$(ls $OutDir/"$Strain"_chen_et_al_2014_CRN_BrakerIntersect.txt)
    while read Line; do
      CrnName=$(echo $Line | cut -f1 -d ' ' )
      HitName=$(echo $Line | cut -f2 -d ' ')
      echo "$CrnName - $HitName"
      if [ $Species == 'P.cactorum' ] || [ $Species == 'P.infestans' ] || [ $Species == 'P.parisitica' ]; then
        echo "$HitName" > tmp.txt
      else
        cat $Proteins | grep -w "$HitName" | tr -d '>' > tmp.txt
      fi
      OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_chen_et_al_2014_CRN/fasta
      mkdir -p $OutDir
      OutFasta=$(echo ""$Species"_"$Strain"_""$CrnName""_homolog.fa" | sed 's/_BlastHit_1//g')
      ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
      $ProgDir/extract_from_fasta.py --fasta $Proteins --headers tmp.txt > $OutDir/$OutFasta
    done < $CrnProteins
    OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_chen_et_al_2014_analysis
    echo ""
    ORFs=$(ls gene_pred/ORF_finder/$Species/$Strain/"$Strain".aa_cat.fa)
    CrnORFs=$(ls $OutDir/"$Strain"_chen_et_al_2014_CRN_ORFIntersect.txt)
    for CrnName in $(cat $CrnORFs | cut -f1 | sort | uniq); do
      cat $CrnORFs | grep "$CrnName" | cut -f2 > tmp.txt
      echo "$CrnName - ORF hits"
      OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_chen_et_al_2014_CRN_analysis/fasta
      mkdir -p $OutDir
      OutFasta=$(echo ""$Species"_"$Strain"_""$CrnName""_homolog.fa" | sed 's/_BlastHit_1//g')
      ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
      $ProgDir/extract_from_fasta.py --fasta $ORFs --headers tmp.txt > $OutDir/$OutFasta
      sed -E "s/^>/>$Species|/g" -i $OutDir/$OutFasta
    done
    echo ""
  done
```

## 4 Previously characterised AVR genes

## 4.2 Performing BLAST searches

BLast searches were performed against each genome:

```bash
  Pcac_ass=repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
  Pinf_ass=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.fa
  Ppar_ass=assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra-310_2_supercontigs.fasta
  Pcap_ass=assembly/external_group/P.capsici/LT1534/dna/Phyca11_unmasked_genomic_scaffolds.fasta
  Psoj_ass=assembly/external_group/P.sojae/P6497/dna/Physo3_AssemblyScaffolds.fasta
  for Assembly in $Pcac_ass $Pinf_ass $Ppar_ass $Pcap_ass $Psoj_ass; do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=analysis/blast_homology/oomycete_avr_genes/appended_oomycete_avr_cds.fasta
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
  for BlastHits in $(ls analysis/blast_homology/*/T30-4/*_appended_oomycete_avr_cds.fasta_homologs.csv); do
    Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)  
    Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    HitsGff=analysis/blast_homology/$Organism/$Strain/"$Strain"_appended_oomycete_avr_cds.fasta_homologs.gff
    Column2=Avr_homolog
    NumHits=1
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
  done
```



### 4.3 Characterising BLAST hits

Genes overlapped by BLAST hits were identified:

```bash
  PcacAugGff=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus_extracted.gff
  PparPubGff=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_transcripts.gtf
  PinfPubGff=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts.gtf
  PinfPubGffParsed=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts_parsed.gtf
  cat $PinfPubGff | sed 's/Supercontig_/supercont/g' > $PinfPubGffParsed
  PcapPubGff=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_genes.gff
  PsojPubGff=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_genes_20110401.gff

  for Proteins in $PcacAugGff $PinfPubGffParsed $PparPubGff $PcapPubGff $PsojPubGff; do
    Species=$(echo "$Proteins" | rev | cut -f4 -d '/' | rev)
    Strain=$(echo "$Proteins" | rev | cut -f3 -d '/' | rev)
    HitsGff=$(ls analysis/blast_homology/$Species/$Strain/"$Strain"_appended_oomycete_avr_cds.fasta_homologs.gff)
    ORFs=$(ls gene_pred/ORF_finder/$Species/$Strain/"$Strain"_ORF_corrected.gff3)
    OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_appended_oomycete_avr_analysis
    mkdir -p $OutDir
    BrakerIntersect=$OutDir/"$Strain"_appended_oomycete_avr_BrakerIntersect.gff
    BrakerNoIntersect=$OutDir/"$Strain"_appended_oomycete_avr_BrakerNoIntersect.gff
    ORFIntersect=$OutDir/"$Strain"_appended_oomycete_avr_ORFIntersect.gff
    AvrProteins=$OutDir/"$Strain"_appended_oomycete_avr_BrakerIntersect.txt
    AvrORFs=$OutDir/"$Strain"_appended_oomycete_avr_ORFIntersect.txt

    echo "The following number of Avr genes had blast hits:"
    cat $HitsGff | wc -l
    echo "The following number of blast hits intersected Braker gene models:"
    bedtools intersect -wa -u -a $HitsGff -b $Proteins > $BrakerIntersect
    cat $BrakerIntersect | wc -l
    echo "The following Blast hits did not intersect Braker gene models:"
    bedtools intersect -v -a $HitsGff -b $Proteins > $BrakerNoIntersect
    cat $BrakerNoIntersect | wc -l
    echo "The following number of blast hits intersected ORF gene models:"
    bedtools intersect -wa -u -a $BrakerNoIntersect -b $ORFs > $ORFIntersect
    cat $ORFIntersect | wc -l
    if [ $Species == 'P.infestans' ] || [ $Species == 'P.parisitica' ]; then
    bedtools intersect -wao -a $HitsGff -b $Proteins | grep -w 'exon' | cut -f9,18 | cut -f1,3 -d ";" | tr -d '"' | tr -d ';' | sed 's/ID=//g' | sed 's/transcript_id //g' | sed 's/transcriptId //g' |sed 's/ /\t/g' > $AvrProteins
    elif [ $Species == 'P.capsici' ] || [ $Species == 'P.sojae' ]; then
    bedtools intersect -wao -a $HitsGff -b $Proteins | grep -w 'exon' | cut -f9,18 | cut -f1,2 -d ";" | tr -d '"' | tr -d ';' | sed 's/ID=//g' | sed 's/name //g' | sed 's/ /\t/g' > $AvrProteins
    else
      bedtools intersect -wao -a $HitsGff -b $Proteins | grep -w 'transcript' | cut -f9,18 | tr -d '"' | tr -d ';' | sed 's/ID=//g' > $AvrProteins
    fi
    bedtools intersect -wao -a $BrakerNoIntersect -b $ORFs | grep -w 'transcript' | cut -f9,18 | cut -d ';' -f1,4 | tr -d '"' | sed 's/;/\t/g' | sed 's/ID=//g' | sed 's/Name=//g' > $AvrORFs
  done
```

Of the 57 oomycete avr genes searched for, 36 had blast hits within the
P. cactorum genome. 34 of these intersected Braker gene models and the remaining
2 intersected ORF gene models.

The genes relating to each Avr BLast hit were extracted.

```bash
  Pcac_pep=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
  Pinf_pep=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
  Ppar_pep=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins.pep.all.fa
  Pcap_pep=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins.fasta
  Psoj_pep=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401.aa.fasta
  for Proteins in $Pcac_pep $Pinf_pep $Ppar_pep $Pcap_pep $Psoj_pep; do
    Species=$(echo "$Proteins" | rev | cut -f4 -d '/' | rev)
    Strain=$(echo "$Proteins" | rev | cut -f3 -d '/' | rev)
    echo "$Species - $Strain"
    OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_appended_oomycete_avr_analysis
    AvrProteins=$(ls $OutDir/"$Strain"_appended_oomycete_avr_BrakerIntersect.txt)
    while read Line; do
      AvrName=$(echo $Line | cut -f1 -d ' ' )
      HitName=$(echo $Line | cut -f2 -d ' ')
      echo "$AvrName - $HitName"
      if [ $Species == 'P.cactorum' ] || [ $Species == 'P.infestans' ] || [ $Species == 'P.parisitica' ]; then
        echo "$HitName" > tmp.txt
      else
        cat $Proteins | grep -w "$HitName" | tr -d '>' > tmp.txt
      fi
      OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_appended_oomycete_avr_analysis/fasta
      mkdir -p $OutDir
      OutFasta=$(echo ""$Species"_"$Strain"_""$AvrName""_homolog.fa" | sed 's/_BlastHit_1//g')
      ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
      $ProgDir/extract_from_fasta.py --fasta $Proteins --headers tmp.txt > $OutDir/$OutFasta
    done < $AvrProteins
    OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_appended_oomycete_avr_analysis
    echo ""
    ORFs=$(ls gene_pred/ORF_finder/$Species/$Strain/"$Strain".aa_cat.fa)
    AvrORFs=$(ls $OutDir/"$Strain"_appended_oomycete_avr_ORFIntersect.txt)
    for AvrName in $(cat $AvrORFs | cut -f1 | sort | uniq); do
      cat $AvrORFs | grep "$AvrName" | cut -f2 > tmp.txt
      echo "$AvrName - ORF hits"
      OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_appended_oomycete_avr_analysis/fasta
      mkdir -p $OutDir
      OutFasta=$(echo ""$Species"_"$Strain"_""$AvrName""_homolog.fa" | sed 's/_BlastHit_1//g')
      ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
      $ProgDir/extract_from_fasta.py --fasta $ORFs --headers tmp.txt > $OutDir/$OutFasta
      sed -E "s/^>/>$Species|/g" -i $OutDir/$OutFasta
    done
    echo ""
  done
```

Individual fasta accessions for Strains were concatentated into a single file:

```bash
  for AvrGene in $(ls analysis/blast_homology/P.*/*/*_appended_oomycete_avr_analysis/fasta/*_homolog.fa | rev | cut -d '/' -f1 | rev | cut -f3- -d '_' | sort | uniq); do
    echo $AvrGene
    OutDir=analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_avr_analysis/concatenated_results
    mkdir -p $OutDir
    ConcatList=$(ls analysis/blast_homology/P.*/*/*_appended_oomycete_avr_analysis/fasta/*$AvrGene)
    cat $ConcatList > $OutDir/$AvrGene
  done
```


#PcF
PcF was not identified in the Braker gene models although 3 ORF models
overlapped the location of the top BLAST hit. To further investigate
interproscan annotations were searched for putative PcF-like proteins:


```bash
cat gene_pred/interproscan/P.cactorum/10300/10300_interproscan.tsv | grep 'PcF'
cat gene_pred/braker/P.cactorum/10300/P.cactorum/augustus_extracted.gff | grep 'g19537'
cat analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_Pinf_Ppar_Pcap_Psoj_orthogroups.txt | grep 'Pcac|g19537' | sed 's/ /\n/g'
```

```
  orthogroup1045:
  Pcac|g19537.t1
  Pcac|g19539.t1
  Pinf|PITG_17290T0
  Pinf|PITG_22670T0
  Pinf|PITG_19001T0
  Pcap|124555
  Pcap|53545
  Pcap|102642
  Pcap|133596
  Pcap|124564
  Pcap|124560
  Pcap|63419
  Pcap|130489
  Pcap|124296
  Pcap|73399
```
