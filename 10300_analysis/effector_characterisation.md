## 4. 2 Ananlysis of RxLR effectors

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
  for GeneGff in $PcacAugGff $PparPubGff $PinfPubGff $PcapPubGff $PsojPubGff; do
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
    cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq | rev | cut -f1 -d '|' | rev > $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt
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

A similar analysis was performed with ORF predicted RxLRs:

```bash
for InDir in $(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*); do
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
```

```
P.cactorum - 10300
Total number of RxLRs in predicted ORFs:
194
Total number of RxLRs shared between prediction sources:
121
The number of combined RxLR containing ORFs containing WY domains are:
74

P.capsici - LT1534
Total number of RxLRs in predicted ORFs:
260
Total number of RxLRs shared between prediction sources:
159
The number of combined RxLR containing ORFs containing WY domains are:
81

P.infestans - T30-4
Total number of RxLRs in predicted ORFs:
447
Total number of RxLRs shared between prediction sources:
254
The number of combined RxLR containing ORFs containing WY domains are:
151

P.parisitica - 310
Total number of RxLRs in predicted ORFs:
342
Total number of RxLRs shared between prediction sources:
211
The number of combined RxLR containing ORFs containing WY domains are:
123

P.sojae - 67593
Total number of RxLRs in predicted ORFs:
322
Total number of RxLRs shared between prediction sources:
188
The number of combined RxLR containing ORFs containing WY domains are:
118
```


The RxLR effectors from both Gene models and ORF finding approaches were
combined into a single file.

```bash
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/P.*/* | grep -v '67593'); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$MergeDir/"$Strain"_Aug_RxLR_EER_motif_hmm.gff
ORFGff=$MergeDir/"$Strain"_ORF_RxLR_EER_motif_hmm.gff
ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_RxLR_EER_motif_hmm.gff
AugInORFs=$MergeDir/"$Strain"_AugInORFs_RxLR_EER_motif_hmm.gff
ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_RxLR_EER_motif_hmm.gff
AugUniq=$MergeDir/"$Strain"_Aug_Uniq_RxLR_EER_motif_hmm.gff
bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
bedtools intersect -wb -u -a $ORFGff -b $AugGff > $AugInORFs
bedtools intersect -v -a $ORFGff -b $AugGff > $ORFsUniq
bedtools intersect -v -a $AugGff -b $ORFGff > $AugUniq
echo "$Species - $Strain"
if $Strain == "10300"; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep -w 'gene' | wc -l
fi
if $Strain == "T30-4"; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | wc -l
fi
if $Strain == "310"; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep 'exon' | rev | cut -f1 -d ';' | rev | sort | uniq | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep 'exon' | rev | cut -f1 -d ';' | rev | sort | uniq | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep 'exon' | rev | cut -f1 -d ';' | rev | sort | uniq | wc -l
fi
if $Strain == "LT1534"; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep 'transcriptId' | rev | cut -f1 -d ':' | rev | sort | uniq | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep 'transcriptId' | rev | cut -f1 -d ':' | rev | sort | uniq | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep 'transcriptId' | rev | cut -f1 -d ':' | rev | sort | uniq | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep 'transcriptId' | rev | cut -f1 -d ':' | rev | sort | uniq | wc -l
fi
if $Strain == "P6497"; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep -w 'gene' | wc -l
fi
done
```


```
P.cactorum - 10300
The number of ORF RxLRs overlapping Augustus RxLRs:
128
The number of Augustus RxLRs overlapping ORF RxLRs:
128
The number of RxLRs unique to ORF models:
66
The number of RxLRs unique to Augustus models:
17

P.capsici - LT1534
Total number of RxLRs in predicted ORFs:
260
Total number of RxLRs shared between prediction sources:
159
The number of combined RxLR containing ORFs containing WY domains are:
81

P.infestans - T30-4
Total number of RxLRs in predicted ORFs:
447
Total number of RxLRs shared between prediction sources:
254
The number of combined RxLR containing ORFs containing WY domains are:
151

P.parisitica - 310
Total number of RxLRs in predicted ORFs:
342
Total number of RxLRs shared between prediction sources:
211
The number of combined RxLR containing ORFs containing WY domains are:
123

P.sojae - 67593
Total number of RxLRs in predicted ORFs:
322
Total number of RxLRs shared between prediction sources:
188
The number of combined RxLR containing ORFs containing WY domains are:
118

```

```bash
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/P.infestans/T30-4); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$MergeDir/10300_Aug_RxLR_EER_motif_hmm.gff3
ORFGff=$MergeDir/10300_ORF_RxLR_EER_motif_hmm.gff
ORFsInAug=$MergeDir/10300_ORFsInAug_RxLR_EER_motif_hmm.gff
AugInORFs=$MergeDir/10300_AugInORFs_RxLR_EER_motif_hmm.gff
ORFsUniq=$MergeDir/10300_ORFsUniq_RxLR_EER_motif_hmm.gff
AugUniq=$MergeDir/10300_Aug_Uniq_RxLR_EER_motif_hmm.gff
bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
bedtools intersect -wb -u -a $ORFGff -b $AugGff > $AugInORFs
bedtools intersect -v -a $ORFGff -b $AugGff > $ORFsUniq
bedtools intersect -v -a $AugGff -b $ORFGff > $AugUniq
echo "$Species - $Strain"
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep -w 'gene' | wc -l
done
```
