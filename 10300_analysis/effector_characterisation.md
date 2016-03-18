## 4.2.a Analysis of RxLR effectors - Augustus gene models

Due to RxLR effectors being predicted from a number of sources the number of
unique RxLRs were identified from motif and Hmm searches within gene models.
221 RxLR effectors were predicted in total from the P.cactorum genome. Of these,
90 were shared between both datasets.
<!--
```bash
for InDir in $(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/*/10300); do
Strain=$(echo "$InDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$InDir" | rev | cut -f2 -d '/' | rev)
Source="pub"
if [ $Strain == '10300' ]; then
Source="Aug"
fi
RxLR_motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Species/$Strain/"$Strain"_"$Source"_RxLR_EER_regex.txt)
RxLR_hmm=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Species/$Strain/"$Strain"_"$Source"_RxLR_hmmer_headers.txt)
WY_hmm=$(ls analysis/RxLR_effectors/hmmer_WY/$Species/$Strain/"$Strain"_"$Source"_WY_hmmer_headers.txt)
echo "$Species - $Strain"
echo "Total number of RxLRs in predicted genes:"
cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq | wc -l
echo "Total number of RxLRs shared between prediction sources:"
cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
OutDir=analysis/RxLR_effectors/combined_evidence/$Species/$Strain
mkdir -p $OutDir
cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq > $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt
echo "The number of combined RxLR containing proteins containing WY domains are:"
cat $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt $WY_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
echo ""
GeneModels=$(ls gene_pred/braker/$Species/$Strain/*/augustus_extracted.gff)
cat $GeneModels | grep -w -f $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt > $OutDir/"$Strain"_Aug_RxLR_EER_motif_hmm.gff3
done
``` -->

```bash
  PcacAugGff=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus_extracted.gff
  PparPubGff=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_transcripts.gtf
  PinfPubGff=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts.gff3
  PcapPubGff=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_genes.gff
  PsojPubGff=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_genes_20110401.gff
  for GeneGff in $PcacAugGff $PparPubGff $PinfPubGff $PsojPubGff; do
  # for GeneGff in $PparPubGff; do
  echo "$GeneGff"
  Strain=$(echo "$GeneGff" | rev | cut -f3 -d '/' | rev)
  Species=$(echo "$GeneGff" | rev | cut -f4 -d '/' | rev)
  InDir=$(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/$Species/$Strain)
  Source="pub"
  if [ $Strain == '10300' ]; then
    Source="Aug"
  fi
  RxLR_motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Species/$Strain/"$Strain"_"$Source"_RxLR_EER_regex.txt)
  RxLR_hmm=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Species/$Strain/"$Strain"_"$Source"_RxLR_hmmer_headers.txt)
  WY_hmm=$(ls analysis/RxLR_effectors/hmmer_WY/$Species/$Strain/"$Strain"_"$Source"_WY_hmmer_headers.txt)
  echo "$Species - $Strain"
  echo "Total number of RxLRs in predicted genes:"
  cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq | wc -l
  echo "Total number of RxLRs shared between prediction sources:"
  cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
  OutDir=analysis/RxLR_effectors/combined_evidence/$Species/$Strain
  mkdir -p $OutDir
  cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq | rev | cut -f2 -d '|' | rev > $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt
  echo "The number of combined RxLR containing proteins containing WY domains are:"
  cat $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt $WY_hmm | cut -f1 -d' ' | rev | cut -f2 -d '|' | rev | sort | uniq -d | wc -l
  echo ""
  cat $GeneGff | grep -w -f $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt > $OutDir/"$Strain"_Aug_RxLR_EER_motif_hmm.gff
  done

  # For P. capsici
  for GeneGff in $PcapPubGff $PsojPubGff; do
    echo "$GeneGff"
    Strain=$(echo "$GeneGff" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$GeneGff" | rev | cut -f4 -d '/' | rev)
    InDir=$(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/$Species/$Strain)
    Source="pub"
    RxLR_motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Species/$Strain/"$Strain"_"$Source"_RxLR_EER_regex.txt)
    RxLR_hmm=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Species/$Strain/"$Strain"_"$Source"_RxLR_hmmer_headers.txt)
    WY_hmm=$(ls analysis/RxLR_effectors/hmmer_WY/$Species/$Strain/"$Strain"_"$Source"_WY_hmmer_headers.txt)
    echo "$Species - $Strain"
    echo "Total number of RxLRs in predicted genes:"
    cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq | wc -l
    echo "Total number of RxLRs shared between prediction sources:"
    cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
    OutDir=analysis/RxLR_effectors/combined_evidence/$Species/$Strain
    mkdir -p $OutDir
    cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq | sed -e 's/|$//g' | rev | cut -f1 -d '|' | rev > $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt
    echo "The number of combined RxLR containing proteins containing WY domains are:"
    cat $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt $WY_hmm | cut -f1 -d' ' | rev | cut -f1 -d '|' | rev | sort | uniq -d | wc -l
    echo ""
    cat $GeneGff | grep -w -f $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt > $OutDir/"$Strain"_Aug_RxLR_EER_motif_hmm.gff
  done
```

```
P.cactorum - 10300
Total number of RxLRs in predicted genes:
145
Total number of RxLRs shared between prediction sources:
88
The number of combined RxLR containing proteins containing WY domains are:
69

P.capsici - LT1534
Total number of RxLRs in predicted genes:
125
Total number of RxLRs shared between prediction sources:
67
The number of combined RxLR containing proteins containing WY domains are:
42

P.infestans - T30-4
Total number of RxLRs in predicted genes:
398
Total number of RxLRs shared between prediction sources:
249
The number of combined RxLR containing proteins containing WY domains are:
157

P.parisitica - 310
Total number of RxLRs in predicted genes:
267
Total number of RxLRs shared between prediction sources:
156
The number of combined RxLR containing proteins containing WY domains are:
108

P.sojae - P6497
Total number of RxLRs in predicted genes:
344
Total number of RxLRs shared between prediction sources:
176
The number of combined RxLR containing proteins containing WY domains are:
132
```

## 4.2.a Analysis of RxLR effectors - ORF gene models

Due to RxLR effectors being predicted from a number of sources the number of
unique RxLRs were identified from motif and Hmm searches within predicted ORFs.

This is a similar analysis to that performed above:

```bash
  for InDir in $(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/*/* | grep -v 'P.infestans'); do
  Strain=$(echo "$InDir" | rev | cut -f1 -d '/' | rev)
  Species=$(echo "$InDir" | rev | cut -f2 -d '/' | rev)
  RxLR_motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Species/$Strain/"$Strain"_ORF_RxLR_EER_regex.txt)
  RxLR_hmm=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Species/$Strain/"$Strain"_ORF_RxLR_hmmer_headers.txt)
  WY_hmm=$(ls analysis/RxLR_effectors/hmmer_WY/$Species/$Strain/"$Strain"_ORF_WY_hmmer_headers.txt)
  echo "$Species - $Strain"
  echo "Total number of RxLRs in predicted ORFs:"
  cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq | wc -l
  echo "Total number of RxLRs shared between prediction sources:"
  cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
  OutDir=analysis/RxLR_effectors/combined_evidence/$Species/$Strain
  mkdir -p $OutDir
  cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq > $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm_headers.txt
  echo "The number of combined RxLR containing ORFs containing WY domains are:"
  cat $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm_headers.txt $WY_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
  echo ""
  GeneModels=$(ls gene_pred/ORF_finder/$Species/$Strain/"$Strain"_ORF_corrected.gff3)
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
  $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm_headers.txt $GeneModels RxLR_EER_combined Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm.gff
  done
  # For P. infestans
  for InDir in $(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/P.infestans/T30-4); do
  Strain=$(echo "$InDir" | rev | cut -f1 -d '/' | rev)
  Species=$(echo "$InDir" | rev | cut -f2 -d '/' | rev)
  RxLR_motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Species/$Strain/"$Strain"_ORF_RxLR_EER_regex.txt)
  RxLR_hmm=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Species/$Strain/"$Strain"_ORF_RxLR_hmmer_headers.txt)
  WY_hmm=$(ls analysis/RxLR_effectors/hmmer_WY/$Species/$Strain/"$Strain"_ORF_WY_hmmer_headers.txt)
  echo "$Species - $Strain"
  echo "Total number of RxLRs in predicted ORFs:"
  cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq | wc -l
  echo "Total number of RxLRs shared between prediction sources:"
  cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
  OutDir=analysis/RxLR_effectors/combined_evidence/$Species/$Strain
  mkdir -p $OutDir
  cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq > $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm_headers.txt
  echo "The number of combined RxLR containing ORFs containing WY domains are:"
  cat $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm_headers.txt $WY_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
  echo ""
  GeneModels=$(ls gene_pred/ORF_finder/$Species/$Strain/"$Strain"_ORF_corrected.gff3)
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
  $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm_headers.txt $GeneModels RxLR_EER_combined Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm.gff
  cat $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm.gff | sed 's/^supercont/Supercontig_/' | sed -e 's/_dna.*\tRxLR_EER_combined/\tRxLR_EER_combined/' > $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm_mod.gff
  done
```

```
P.cactorum - 10300
Total number of RxLRs in predicted ORFs:
192
Total number of RxLRs shared between prediction sources:
121
The number of combined RxLR containing ORFs containing WY domains are:
74

P.capsici - LT1534
Total number of RxLRs in predicted ORFs:
250
Total number of RxLRs shared between prediction sources:
152
The number of combined RxLR containing ORFs containing WY domains are:
76

P.parisitica - 310
Total number of RxLRs in predicted ORFs:
317
Total number of RxLRs shared between prediction sources:
204
The number of combined RxLR containing ORFs containing WY domains are:
114

P.sojae - P6497
Total number of RxLRs in predicted ORFs:
317
Total number of RxLRs shared between prediction sources:
187
The number of combined RxLR containing ORFs containing WY domains are:
117

P.infestans - T30-4
Total number of RxLRs in predicted ORFs:
437
Total number of RxLRs shared between prediction sources:
251
The number of combined RxLR containing ORFs containing WY domains are:
148

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
# for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/P.*/*); do
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/P.*/P6497); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$MergeDir/"$Strain"_Aug_RxLR_EER_motif_hmm.gff
ORFGff=$MergeDir/"$Strain"_ORF_RxLR_EER_motif_hmm.gff
WY_Aug_hmm=$(ls analysis/RxLR_effectors/hmmer_WY/$Species/$Strain/"$Strain"_*_WY_hmmer_headers.txt | grep -v 'ORF')
WY_ORF_hmm=$(ls analysis/RxLR_effectors/hmmer_WY/$Species/$Strain/"$Strain"_ORF_WY_hmmer_headers.txt)
if [ $Species == P.infestans ]; then
ORFGff=$MergeDir/"$Strain"_ORF_RxLR_EER_motif_hmm_mod.gff
fi
ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_RxLR_EER_motif_hmm.gff
AugInORFs=$MergeDir/"$Strain"_AugInORFs_RxLR_EER_motif_hmm.gff
ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_RxLR_EER_motif_hmm.gff
AugUniq=$MergeDir/"$Strain"_Aug_Uniq_RxLR_EER_motif_hmm.gff
TotalRxLRsTxt=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.txt
TotalRxLRsGff=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.gff
TotalRxLRsWYTxt=$MergeDir/"$Strain"_Total_RxLR_EER_WY_motif_hmm.txt
TotalRxLRsWYGff=$MergeDir/"$Strain"_Total_RxLR_EER_WY_motif_hmm.gff
# TotalRxLRsHeaders=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm_headers.txt
bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
# bedtools intersect -wb -u -a $ORFGff -b $AugGff > $AugInORFs
bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq
echo "$Species - $Strain"
if [ $Strain == "10300" ]; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' | wc -l
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' >> $TotalRxLRsTxt
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep -w 'gene' | wc -l
echo "The total number of putative RxLRs are:"
cat $AugInORFs | grep -w 'gene' | cut -f9 > $TotalRxLRsTxt
cat $AugUniq | grep -w 'gene' | cut -f9 >> $TotalRxLRsTxt
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f4 -d '=' >> $TotalRxLRsTxt
cat $TotalRxLRsTxt | wc -l
echo "The number of these RxLRs containing WY domains are:"
cat $TotalRxLRsTxt $WY_Aug_hmm $WY_ORF_hmm | cut -f1 -d ' ' | rev | cut -f2 -d '|' | rev | sort | uniq -d > $TotalRxLRsWYTxt
cat $TotalRxLRsWYTxt | wc -l
elif [ $Strain == "T30-4" ]; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' | wc -l
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' >> $TotalRxLRsTxt
echo "The number of RxLRs unique to Augustus models:"
cat $AugInORFs | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq > Augtmp.txt
cat $AugUniq | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | grep -v -f Augtmp.txt | wc -l
echo "The total number of putative RxLRs are:"
cat $AugInORFs | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq > $TotalRxLRsTxt
cat $AugUniq | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | grep -v -f Augtmp.txt >> $TotalRxLRsTxt
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f4 -d '=' >> $TotalRxLRsTxt
cat $TotalRxLRsTxt | wc -l
echo "The number of these RxLRs containing WY domains are:"
cat $TotalRxLRsTxt $WY_Aug_hmm $WY_ORF_hmm | cut -f1 -d ' ' | rev | cut -f2 -d '|' | rev | sort | uniq -d > $TotalRxLRsWYTxt
cat $TotalRxLRsWYTxt | wc -l
elif [ $Strain == "310" ]; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' | wc -l
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' >> $TotalRxLRsTxt
echo "The number of RxLRs unique to Augustus models:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | cut -f2 -d '"'> Augtmp.txt
cat $AugUniq | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | grep -v -f Augtmp.txt | wc -l
echo "The total number of putative RxLRs are:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq|  cut -f2 -d '"' > $TotalRxLRsTxt
cat $AugUniq | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | cut -f2 -d '"' | grep -v -f Augtmp.txt >> $TotalRxLRsTxt
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f4 -d '=' >> $TotalRxLRsTxt
cat $TotalRxLRsTxt | wc -l
echo "The number of these RxLRs containing WY domains are:"
cat $TotalRxLRsTxt $WY_Aug_hmm $WY_ORF_hmm | cut -f1 -d ' ' | rev | cut -f2 -d '|' | rev | sort | uniq -d > $TotalRxLRsWYTxt
cat $TotalRxLRsWYTxt | wc -l
elif [ $Strain == "LT1534" ]; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' | wc -l
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' >> $TotalRxLRsTxt
echo "The number of RxLRs unique to Augustus models:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq > Augtmp.txt
cat $AugUniq | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq | grep -v -f Augtmp.txt | wc -l
echo "The total number of putative RxLRs are:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq > $TotalRxLRsTxt
cat $AugUniq | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq | grep -v -f Augtmp.txt >> $TotalRxLRsTxt
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f4 -d '=' >> $TotalRxLRsTxt
cat $TotalRxLRsTxt | wc -l
echo "The number of these RxLRs containing WY domains are:"
cat $TotalRxLRsTxt $WY_Aug_hmm $WY_ORF_hmm | cut -f1 -d ' ' | sed -e 's/|$//g'| rev | cut -f1 -d '|' | rev | sort | uniq -d > $TotalRxLRsWYTxt
cat $TotalRxLRsWYTxt | wc -l
elif [ $Strain == "P6497" ]; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' | wc -l
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' >> $TotalRxLRsTxt
echo "The number of RxLRs unique to Augustus models:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq > Augtmp.txt
cat $AugUniq | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq | grep -v -f Augtmp.txt | wc -l
echo "The total number of putative RxLRs are:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq > $TotalRxLRsTxt
cat $AugUniq | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq | grep -v -f Augtmp.txt >> $TotalRxLRsTxt
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f4 -d '=' >> $TotalRxLRsTxt
cat $TotalRxLRsTxt | wc -l
echo "The number of these RxLRs containing WY domains are:"
cat $TotalRxLRsTxt $WY_Aug_hmm $WY_ORF_hmm | cut -f1 -d ' ' | sed -e 's/|$//g'| rev | cut -f1 -d '|' | rev | sort | uniq -d > $TotalRxLRsWYTxt
cat $TotalRxLRsWYTxt | wc -l
fi
cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalRxLRsTxt > $TotalRxLRsGff
cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalRxLRsWYTxt > $TotalRxLRsWYGff
done
```


```
  P.cactorum - 10300
  The number of ORF RxLRs overlapping Augustus RxLRs:
  128
  The number of Augustus RxLRs overlapping ORF RxLRs:
  128
  The number of RxLRs unique to ORF models:
  64
  The number of RxLRs unique to Augustus models:
  17
  The total number of putative RxLRs are:
  209
  The number of these RxLRs containing WY domains are:
  79
  P.capsici - LT1534
  The number of ORF RxLRs overlapping Augustus RxLRs:
  110
  The number of Augustus RxLRs overlapping ORF RxLRs:
  110
  The number of RxLRs unique to ORF models:
  140
  The number of RxLRs unique to Augustus models:
  15
  The total number of putative RxLRs are:
  265
  The number of these RxLRs containing WY domains are:
  82
  P.infestans - T30-4
  The number of ORF RxLRs overlapping Augustus RxLRs:
  372
  The number of Augustus RxLRs overlapping ORF RxLRs:
  372
  The number of RxLRs unique to ORF models:
  65
  The number of RxLRs unique to Augustus models:
  26
  The total number of putative RxLRs are:
  463
  The number of these RxLRs containing WY domains are:
  157
  P.parisitica - 310
  The number of ORF RxLRs overlapping Augustus RxLRs:
  235
  The number of Augustus RxLRs overlapping ORF RxLRs:
  235
  The number of RxLRs unique to ORF models:
  82
  The number of RxLRs unique to Augustus models:
  32
  The total number of putative RxLRs are:
  349
  The number of these RxLRs containing WY domains are:
  124
  P.sojae - P6497
  The number of ORF RxLRs overlapping Augustus RxLRs:
  276
  The number of Augustus RxLRs overlapping ORF RxLRs:
  276
  The number of RxLRs unique to ORF models:
  41
  The number of RxLRs unique to Augustus models:
  68
  The total number of putative RxLRs are:
  385
  The number of these RxLRs containing WY domains are:
  133
```

Fasta sequences for RxLRs were extracted for each isolate

```bash
  PcacFa=$(ls gene_pred/braker/*/10300/*/augustus.aa)
  PparFa=$(ls assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins.pep.all.fa)
  PinfFa=$(ls assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all_parsed.fa)
  PcapFa=$(ls assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins.fasta)
  PsojFa=$(ls assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401.aa.fasta)

  for AugFa in $PcacFa $PparFa $PinfFa $PcapFa $PsojFa; do
  # for AugFa in $PinfFa; do
    Strain=$(echo "$AugFa" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$AugFa" | rev | cut -f4 -d '/' | rev)
    ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
    MergeDir=analysis/RxLR_effectors/combined_evidence/$Species/$Strain
    TotalRxLRsTxt=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.txt
    RxLRsFa=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm_headers.fa
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
    $ProgDir/unwrap_fasta.py --inp_fasta $AugFa | grep -A1 -w -f $TotalRxLRsTxt | grep -v -E '^--$' > $RxLRsFa
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalRxLRsTxt >> $RxLRsFa
    echo "$Strain"
    echo "The number of sequences extracted is"
    cat $RxLRsFa | grep '>' | wc -l
  done
```

```
  10300
    The number of sequences extracted is
    209
  310
    The number of sequences extracted is
    349
  T30-4
    The number of sequences extracted is
    463
  LT1534
    The number of sequences extracted is
    265
  P6497
    The number of sequences extracted is
    131
```


## 4.2.d Expression of P.cactorum 10300 RxLR genes

Expression data of 10300 was used to provide expression support for P. cactorum
RxLR genes.

This was done by intersecting the location of RxLR with the RNAseq data aligned
to the 10300 genome.

```bash
MergeDir=$(ls -d analysis/RxLR_effectors/combined_evidence/P.*/10300)
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
TotalRxLRsGff=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm_headers.gff

cufflinks -o tmp -p 16 -G gene_pred/braker/P.cactorum/10300/P.cactorum/augustus_extracted.gff alignment/P.cactorum/10300/accepted_hits.bam
bedtools intersect -s -u -a tmp/transcripts.gtf -b $TotalRxLRsGff  > analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/"$Strain"_Total_RxLR_EER_motif_hmm_expressed.gtf
echo "The top 20 expressed RxLRs are:"
cat analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/"$Strain"_Total_RxLR_EER_motif_hmm_expressed.gtf | grep -w 'transcript' | sort -r -n -k 14 -t '"' | cut -f2,14 -d '"' --output-delimite " - " |  head -n 20
echo "The total number of RxLRs was:"
cat analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/"$Strain"_Total_RxLR_EER_motif_hmm_expressed.gtf | grep -w 'transcript' | wc -l
echo "The number of RxLRs with 1x coverage or greater was:"
cat analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/"$Strain"_Total_RxLR_EER_motif_hmm_expressed.gtf | grep -w 'transcript' | sort -r -n -k 14 -t '"' | cut -f14 -d '"' | grep -v -E '0\.' | wc -l
echo "The number of RxLRs with 0x coverage was:"
cat analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/"$Strain"_Total_RxLR_EER_motif_hmm_expressed.gtf | grep -w 'transcript' | sort -r -n -k 14 -t '"' | cut -f14 -d '"' | grep -E '0\.00' | wc -l
```

```
  The top 20 expressed RxLRs are:
  g11594 - 29.889056 - Pcac unique - orthogroup8833
  g11370 - 29.864802 - Pcac unique - orthogroup8833
  g12117 - 5.227685
  g18158 - 5.131975
  g8566 - 4.951625
  g17215 - 4.929356
  g14216 - 4.866843
  g14208 - 4.818298
  g1876 - 4.194974
  g8314 - 2.084139
  g5978 - 2.010123
  g11506 - 1.970131
  g19634 - 1.808004
  g5288 - 1.629739 - clade 1 orthogroup
  g7929 - 1.503567
  g8468 - 1.450752
  g3902 - 1.238576
  g15964 - 0.988131
  g12572 - 0.933141
  g3433 - 0.900869
  The total number of RxLRs was:
  160
  The number of RxLRs with 1x coverage or greater was:
  17
  The number of RxLRs with 0x coverage was:
  81
```


The expression of the P. cactorum unique RxLR g1610 was investigated. It was
found to have no evidence of expression.

The P. cactorum uniuq ortholog group of 2 RxLRs g11370 and g11594 were both
expressed at 30x coverage.


```bash
  cat analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/"$Strain"_Total_RxLR_EER_motif_hmm_expressed.gtf | grep -w 'g1610' | less
  cat analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/Group1_RxLR_Orthogroups_hits.txt
  cat analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/"$Strain"_Total_RxLR_EER_motif_hmm_expressed.gtf | grep -w 'g11370' | less
```

## 4.2.f Functional annotation of 10300 RxLRs

Interproscan annotations and swissprot similarities were identified for 10300
RxLRs. Non of the RxLRs carried an interproscan annotation or BLAST
homology >1e-100 to a swissprot protein.

```bash
  MergeDir=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300
  Strain=10300
  for Gene in $(cat $MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm_expressed.gtf | grep -w 'transcript' | sort -r -n -k 14 -t '"' | cut -f4,14 -d '"' --output-delimite " - " | cut -f1 -d '.'); do     
    echo $Gene;
    cat gene_pred/interproscan/P.cactorum/10300/10300_interproscan.tsv | grep '$Gene';
    cat gene_pred/swissprot/P.cactorum/10300/swissprot_v2015_10_hits.tbl | grep '$Gene';
    echo "";
  done
```

## 4.3.a Analysis of Crinkler effectors - merger of Augustus / published genes with ORFs


Intersection between the coodinates of putative CRNs from gene models and ORFs
were identified to determine the total number of CRNs predicted in these
genomes.

The CRN effectors from both Gene models and ORF finding approaches were
combined into a single file.

Again, this step was complicated by the inconsistency in downloaded gff files for gene
models.


Extract crinklers from published gene models

```bash
  PcacAugGff=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus_extracted.gff
  PparPubGff=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_transcripts.gtf
  PinfPubGff=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts.gff3
  PcapPubGff=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_genes.gff
  PsojPubGff=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_genes_20110401.gff
  for GeneGff in $PcacAugGff; do
    echo "$GeneGff"
    Strain=$(echo "$GeneGff" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$GeneGff" | rev | cut -f4 -d '/' | rev)
    CrnDir=$(ls -d analysis/CRN_effectors/hmmer_CRN/$Species/$Strain)
    Source="pred"
    CRN_hmm_txt=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.txt
    CRN_hmm_txt_mod=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL_mod.txt
    CRN_hmm_gff=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.gff
    cat $CRN_hmm_txt | sed -E 's/\.t.+$//g' > $CRN_hmm_txt_mod
    echo "$Species - $Strain"
    cat $GeneGff | grep -w -f $CRN_hmm_txt_mod > $CRN_hmm_gff
    cat $CRN_hmm_gff | cut -f2 -d '"' | sort | uniq | wc -l
  done
  # For P. capsici & P. sojae
  for GeneGff in $PcapPubGff $PsojPubGff; do
    echo "$GeneGff"
    Strain=$(echo "$GeneGff" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$GeneGff" | rev | cut -f4 -d '/' | rev)
    CrnDir=$(ls -d analysis/CRN_effectors/hmmer_CRN/$Species/$Strain)
    Source="pub"
    CRN_hmm_txt=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.txt
    CRN_hmm_txt_mod=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL_mod.txt
    CRN_hmm_gff=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.gff
    cat $CRN_hmm_txt | cut -f4 -d '|' > $CRN_hmm_txt_mod
    echo "$Species - $Strain"
    cat $GeneGff | grep -w -f $CRN_hmm_txt_mod > $CRN_hmm_gff
    rm $CRN_hmm_txt_mod
    cat $CRN_hmm_gff | cut -f2 -d '"' | sort | uniq | wc -l
  done
  # For P. infestans & P. parisitica
  for GeneGff in $PinfPubGff $PparPubGff; do
    echo "$GeneGff"
    Strain=$(echo "$GeneGff" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$GeneGff" | rev | cut -f4 -d '/' | rev)
    CrnDir=$(ls -d analysis/CRN_effectors/hmmer_CRN/$Species/$Strain)
    Source="pub"
    CRN_hmm_txt=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.txt
    CRN_hmm_gff=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.gff
    echo "$Species - $Strain"
    cat $GeneGff | grep -w -f $CRN_hmm_txt > $CRN_hmm_gff
    cat $CRN_hmm_gff | grep 'exon' | cut -f9 | cut -f2 -d ':' | sort | uniq | wc -l
  done
```


```bash
  for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/*/10300 | grep -v -e '67593' -e 'masked'); do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    AugGff=$MergeDir/"$Strain"_pub_CRN_LFLAK_DWL.gff
    ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer.gff3
    if [ $Species == P.infestans ]; then
      cat $ORFGff | sed 's/^supercont/Supercontig_/' | sed -e 's/_dna.*\tCRN_HMM/\tCRN_HMM/' > $MergeDir/"$Strain"_CRN_merged_hmmer_mod.gff3
      ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer_mod.gff3
    fi
    ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_CRN_hmmer.bed
    AugInORFs=$MergeDir/"$Strain"_AugInORFs_CRN_hmmer.bed
    ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_CRN_hmmer.bed
    AugUniq=$MergeDir/"$Strain"_Aug_Uniq_CRN_hmmer.bed
    TotalCRNsTxt=$MergeDir/"$Strain"_Total_CRN.txt
    TotalCRNsGff=$MergeDir/"$Strain"_Total_CRN.gff
    TotalCRNsHeaders=$MergeDir/"$Strain"_Total_CRN_headers.txt
    bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
    bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
    bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
    bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq
    echo "$Species - $Strain"
    if [ $Strain == "10300" ]; then
      echo "The number of ORF CRNs overlapping Augustus CRNs:"
      cat $ORFsInAug | grep -w 'gene' | wc -l
      cat $ORFsInAug | grep -w 'gene' > $TotalCRNsTxt
      echo "The number of Augustus CRNs overlapping ORF CRNs:"
      cat $AugInORFs | grep -w 'gene' | wc -l
      cat $AugInORFs | grep -w 'gene' >> $TotalCRNsTxt
      echo "The number of CRNs unique to ORF models:"
      cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' | wc -l
      cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' >> $TotalCRNsTxt
      echo "The number of CRNs unique to Augustus models:"
      cat $AugUniq | grep -w 'gene' | wc -l
      cat $AugUniq | grep -w 'gene' >> $TotalCRNsTxt
    elif [ $Strain == "T30-4" ]; then
      echo "The number of ORF CRNs overlapping Augustus CRNs:"
      cat $ORFsInAug | grep -w 'gene' | wc -l
      cat $ORFsInAug | grep -w 'gene' > $TotalCRNsTxt
      echo "The number of Augustus CRNs overlapping ORF CRNs:"
      cat $AugInORFs | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | wc -l
      cat $AugInORFs | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq >> $TotalCRNsTxt
      echo "The number of CRNs unique to ORF models:"
      cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' | wc -l
      cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' >> $TotalCRNsTxt
      echo "The number of CRNs unique to Augustus models:"
      cat $AugInORFs | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq > Augtmp.txt
      cat $AugUniq | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | grep -v -f Augtmp.txt | wc -l
      cat $AugUniq | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | grep -v -f Augtmp.txt >> $TotalCRNsTxt
    elif [ $Strain == "310" ]; then
      echo "The number of ORF CRNs overlapping Augustus CRNs:"
      cat $ORFsInAug | grep -w 'gene' | wc -l
      cat $ORFsInAug | grep -w 'gene' > $TotalCRNsTxt
      echo "The number of Augustus CRNs overlapping ORF CRNs:"
      cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | wc -l
      cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq >> $TotalCRNsTxt
      echo "The number of CRNs unique to ORF models:"
      cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' | wc -l
      cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' >> $TotalCRNsTxt
      echo "The number of CRNs unique to Augustus models:"
      cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq > Augtmp.txt
      cat $AugUniq | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | grep -v -f Augtmp.txt | wc -l
      cat $AugUniq | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | grep -v -f Augtmp.txt >> $TotalCRNsTxt
    elif [ $Strain == "LT1534" ]; then
      echo "The number of ORF CRNs overlapping Augustus CRNs:"
      cat $ORFsInAug | grep -w 'gene' | wc -l
      cat $ORFsInAug | grep -w 'gene' > $TotalCRNsTxt
      echo "The number of Augustus CRNs overlapping ORF CRNs:"
      cat $AugInORFs | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | wc -l
      cat $AugInORFs | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq >> $TotalCRNsTxt
      echo "The number of CRNs unique to ORF models:"
      cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' | wc -l
      cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' >> $TotalCRNsTxt
      echo "The number of CRNs unique to Augustus models:"
      cat $AugInORFs  | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | sed -e 's/^ //g' > Augtmp.txt
      cat $AugUniq | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | grep -v -f Augtmp.txt | wc -l
      cat $AugUniq | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | grep -v -f Augtmp.txt >> $TotalCRNsTxt
      echo ""
    elif [ $Strain == "P6497" ]; then
      echo "The number of ORF CRNs overlapping Augustus CRNs:"
      cat $ORFsInAug | grep -w 'gene' | wc -l
      cat $ORFsInAug | grep -w 'gene' > $TotalCRNsTxt
      echo "The number of Augustus CRNs overlapping ORF CRNs:"
      cat $AugInORFs | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | wc -l
      cat $AugInORFs | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq >> $TotalCRNsTxt
      echo "The number of CRNs unique to ORF models:"
      cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' | wc -l
      cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f5 -d '=' >> $TotalCRNsTxt
      echo "The number of CRNs unique to Augustus models:"
      cat $AugInORFs  | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | sed -e 's/^ //g' > Augtmp.txt
      cat $AugUniq | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | grep -v -f Augtmp.txt | wc -l
      cat $AugUniq | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | grep -v -f Augtmp.txt >> $TotalCRNsTxt
      echo ""
    fi
    cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalCRNsTxt > $TotalCRNsGff
    echo "The total number of CRNs are:"
    if [ $Strain == "10300" ]; then
      cat $TotalCRNsGff | grep 'AUGUSTUS' | cut -f9 > $TotalCRNsHeaders
      cat $TotalCRNsGff | grep -o -E -e 'Name=.*$' -e 'name ".*";' | sed -e 's/^Name=//g' | sed 's/^name "//g' | sed 's/";$//g' | sort | uniq >> $TotalCRNsHeaders
      cat $TotalCRNsHeaders | wc -l
    elif [ $Strain == "310" ]; then
      cat $TotalCRNsGff | grep -o -E -e 'Name=.*$' -e 'gene_id ".*"; t' | sed -e 's/^Name=//g' | sed 's/^gene_id "//g' | sed 's/"; t$//g' | sort | uniq > $TotalCRNsHeaders
      cat $TotalCRNsHeaders | wc -l
    elif [ $Strain == "T30-4" ]; then
      cat $TotalCRNsGff | grep -o -E -e 'Name=.*$' -e 'Parent=.*T.$' | sed -e 's/^Name=//g' | sed 's/Parent=//g' | sed 's/";$//g' | sort | uniq > $TotalCRNsHeaders
      cat $TotalCRNsHeaders | wc -l
    elif [ $Strain == "LT1534" ]; then
      cat $TotalCRNsGff | grep -o -E -e 'Name=.*$' -e 'name ".*";' | sed -e 's/^Name=//g' | sed 's/^name "//g' | sed 's/";$//g' | sort | uniq > $TotalCRNsHeaders
      cat $TotalCRNsHeaders | wc -l
    elif [ $Strain == "P6497" ]; then
      cat $TotalCRNsGff | grep -o -E -e 'Name=.*$' -e 'name ".*";' | sed -e 's/^Name=//g' | sed 's/^name "//g' | sed 's/";$//g' | sort | uniq > $TotalCRNsHeaders
      cat $TotalCRNsHeaders | wc -l
    fi
  done
```

```
P.cactorum - 10300
The number of ORF CRNs overlapping Augustus CRNs:
55
The number of Augustus CRNs overlapping ORF CRNs:
55
The number of CRNs unique to ORF models:
3
The number of CRNs unique to Augustus models:
12
The total number of CRNs are:
70

P.capsici - LT1534
The number of ORF CRNs overlapping Augustus CRNs:
72
The number of Augustus CRNs overlapping ORF CRNs:
71
The number of CRNs unique to ORF models:
32
The number of CRNs unique to Augustus models:
11
The total number of CRNs are:
114

P.infestans - T30-4
The number of ORF CRNs overlapping Augustus CRNs:
157
The number of Augustus CRNs overlapping ORF CRNs:
157
The number of CRNs unique to ORF models:
98
The number of CRNs unique to Augustus models:
10
The total number of CRNs are:
265

P.parisitica - 310
The number of ORF CRNs overlapping Augustus CRNs:
22
The number of Augustus CRNs overlapping ORF CRNs:
22
The number of CRNs unique to ORF models:
4
The number of CRNs unique to Augustus models:
9
The total number of CRNs are:
35

P.sojae - P6497
The number of ORF CRNs overlapping Augustus CRNs:
96
The number of Augustus CRNs overlapping ORF CRNs:
89
The number of CRNs unique to ORF models:
51
The number of CRNs unique to Augustus models:
19
The total number of CRNs are:
159
```

<!-- ```
P.cactorum - 10300
The number of ORF CRNs overlapping Augustus CRNs:
77
The number of Augustus CRNs overlapping ORF CRNs:
77
The number of CRNs unique to ORF models:
24
The number of CRNs unique to Augustus models:
15
The total number of CRNs are:
116
P.capsici - LT1534
The number of ORF CRNs overlapping Augustus CRNs:
90
The number of Augustus CRNs overlapping ORF CRNs:
88
The number of CRNs unique to ORF models:
69
The number of CRNs unique to Augustus models:
18

The total number of CRNs are:
175
P.infestans - T30-4
The number of ORF CRNs overlapping Augustus CRNs:
177
The number of Augustus CRNs overlapping ORF CRNs:
177
The number of CRNs unique to ORF models:
190
The number of CRNs unique to Augustus models:
10
The total number of CRNs are:
377
P.parisitica - 310
The number of ORF CRNs overlapping Augustus CRNs:
36
The number of Augustus CRNs overlapping ORF CRNs:
36
The number of CRNs unique to ORF models:
17
The number of CRNs unique to Augustus models:
6
The total number of CRNs are:
59
P.sojae - P6497
The number of ORF CRNs overlapping Augustus CRNs:
142
The number of Augustus CRNs overlapping ORF CRNs:
133
The number of CRNs unique to ORF models:
70
The number of CRNs unique to Augustus models:
19
The total number of CRNs are:
222
``` -->

Fasta sequences for CRNs were extracted from each isolate

```bash
  PcacFa=$(ls gene_pred/braker/*/10300/*/augustus.aa)
  PparFa=$(ls assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins.pep.all.fa)
  PinfFa=$(ls assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all_parsed.fa)
  PcapFa=$(ls assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins.fasta)
  PsojFa=$(ls assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401.aa.fasta)

  for AugFa in $PcacFa $PparFa $PinfFa $PcapFa $PsojFa; do
  # for AugFa in $PinfFa; do
    Strain=$(echo "$AugFa" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$AugFa" | rev | cut -f4 -d '/' | rev)
    ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
    MergeDir=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain
    TotalCRNsHeaders=$MergeDir/"$Strain"_Total_CRN_headers.txt
    CRNsFa=$MergeDir/"$Strain"_Total_CRN.fa
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
    $ProgDir/unwrap_fasta.py --inp_fasta $AugFa | grep -A1 -w -f $TotalCRNsHeaders | grep -v -E '^--$' > $CRNsFa
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalCRNsHeaders >> $CRNsFa
    echo "$Strain"
    echo "The number of sequences extracted is"
    cat $CRNsFa | grep '>' | wc -l
  done
```



## 4.3.b Expression of P.cactorum 10300 Crinkler genes



Expression data of 10300 was used to provide expression support for P. cactorum
CRN genes.

This was done by intersecting the location of CRN with the RNAseq data aligned
to the 10300 genome.

```bash
  MergeDir=$(ls -d analysis/CRN_effectors/hmmer_CRN/P.*/10300)
  Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
  TotalCRNsGff=$MergeDir/"$Strain"_Total_CRN.gff
  cufflinks -o tmp -p 16 -G gene_pred/braker/P.cactorum/10300/P.cactorum/augustus_extracted.gff alignment/P.cactorum/10300/accepted_hits.bam
  bedtools intersect -s -u -a tmp/transcripts.gtf -b $TotalCRNsGff  > $MergeDir/"$Strain"_Total_CRN_expressed.gtf
  echo "The top 20 expressed CRNs are:"
  cat $MergeDir/"$Strain"_Total_CRN_expressed.gtf | grep -w 'transcript' | sort -r -n -k 14 -t '"' | cut -f4,14 -d '"' --output-delimite " - " |  head -n 20
  echo "The total number of CRNs was:"
  cat $MergeDir/"$Strain"_Total_CRN_expressed.gtf | grep -w 'transcript' | wc -l
  echo "The number of CRNs with 1x coverage or greater was:"
  cat $MergeDir/"$Strain"_Total_CRN_expressed.gtf | grep -w 'transcript' | sort -r -n -k 14 -t '"' | cut -f14 -d '"' | grep -v -E '0\.' | wc -l
  echo "The number of CRNs with 0x coverage was:"
  cat $MergeDir/"$Strain"_Total_CRN_expressed.gtf | grep -w 'transcript' | sort -r -n -k 14 -t '"' | cut -f14 -d '"' | grep -E '0\.00' | wc -l
```

```
  The top 20 expressed CRNs are:
  g1786.t1 - 781.976786
  g17786.t1 - 383.018279
  g18873.t1 - 373.339419
  g13219.t1 - 346.771781
  g14548.t1 - 266.755802
  g19529.t1 - 254.736537
  g14980.t1 - 247.070946
  g778.t1 - 243.816930
  g15212.t1 - 142.053576
  g8090.t1 - 96.736487
  g15485.t1 - 87.810172
  g15471.t1 - 60.885705
  g14356.t1 - 56.727299
  g17134.t1 - 28.375429
  g10998.t1 - 26.094800
  g14894.t1 - 19.551693
  g10526.t1 - 18.019939
  g445.t1 - 16.065962
  g20243.t1 - 14.698188
  g11419.t1 - 13.194537
  The total number of CRNs was:
  67
  The number of CRNs with 1x coverage or greater was:
  42
  The number of CRNs with 0x coverage was:
  11
```

The ortholog groups that these genes belonged to were investigated:

```bash
  for Gene in $(cat $MergeDir/"$Strain"_Total_CRN_expressed.gtf | grep -w 'transcript' | sort -r -n -k 14 -t '"' | cut -f4,14 -d '"' --output-delimite " - " |  head -n 20 | cut -f1 -d '.'); do
    echo $Gene;
    cat analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_Pinf_Ppar_Pcap_Psoj_orthogroups.txt | grep -w "$Gene";
    echo "";
  done
```

Expression of Group 1 unique Crinklers was determined. These genes showed no evidence of expression:

```bash
  cat analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_CRN/Group1_CRN_Orthogroups_hits.txt
  cat $MergeDir/"$Strain"_Total_CRN_expressed.gtf | grep -w 'g16875'
  cat $MergeDir/"$Strain"_Total_CRN_expressed.gtf | grep -w 'g3179'
```

## 4.3.c Functional annotation of 10300 Crinklers

Interproscan annotations and swissprot similarities were identified for 10300
crinklers. Non of the crinklers carried an interproscan annotation or BLAST
homology >1e-100 to a swissprot protein.

```bash
  for Gene in $(cat $MergeDir/"$Strain"_Total_CRN_expressed.gtf | grep -w 'transcript' | sort -r -n -k 14 -t '"' | cut -f4,14 -d '"' --output-delimite " - " | cut -f1 -d '.'); do     
    echo $Gene;
    cat gene_pred/interproscan/P.cactorum/10300/10300_interproscan.tsv | grep '$Gene';
    cat gene_pred/swissprot/P.cactorum/10300/swissprot_v2015_10_hits.tbl  | grep '$Gene';
    echo "";
  done
```
