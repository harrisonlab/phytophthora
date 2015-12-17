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
  for GeneGff in $PcapPubGff; do
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
  for GeneGff in $PcapPubGff; do
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
## 4.2.c Analysis of RxLR effectors - merger of Augustus / published genes with ORFs

Intersection between the coodinates of putative RxLRs from gene models and ORFs
were identified to determine the total number of RxLRs predicted in these
genomes.

The RxLR effectors from both Gene models and ORF finding approaches were
combined into a single file.

This step was complicated by the inconsistency in downloaded gff files for gene
models.


```bash
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/P.*/T30-4 | grep -v '67593'); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$MergeDir/"$Strain"_Aug_RxLR_EER_motif_hmm.gff
ORFGff=$MergeDir/"$Strain"_ORF_RxLR_EER_motif_hmm.gff
if [ $Species == P.infestans ]; then
  ORFGff=$MergeDir/"$Strain"_ORF_RxLR_EER_motif_hmm_mod.gff
fi
ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_RxLR_EER_motif_hmm.gff
AugInORFs=$MergeDir/"$Strain"_AugInORFs_RxLR_EER_motif_hmm.gff
ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_RxLR_EER_motif_hmm.gff
AugUniq=$MergeDir/"$Strain"_Aug_Uniq_RxLR_EER_motif_hmm.gff
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
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep -w 'gene' | wc -l
elif [ $Strain == "T30-4" ]; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugInORFs | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq > Augtmp.txt
cat $AugUniq | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | grep -v -f Augtmp.txt | wc -l
elif [ $Strain == "310" ]; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | wc -l
elif [ $Strain == "LT1534" ]; then
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to Augustus models:"
cat $AugInORFs  | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | sed -e 's/^ //g' > Augtmp.txt
cat $AugUniq | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | grep -v -f Augtmp.txt | wc -l
echo ""
elif [ $Strain == "P6497" ]; then
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
P.infestans - T30-4
The number of ORF RxLRs overlapping Augustus RxLRs:
380
The number of Augustus RxLRs overlapping ORF RxLRs:
380
The number of RxLRs unique to ORF models:
67
The number of RxLRs unique to Augustus models:
18
P.capsici - LT1534
The number of ORF RxLRs overlapping Augustus RxLRs:
114
The number of Augustus RxLRs overlapping ORF RxLRs:
114
The number of RxLRs unique to ORF models:
146
The number of RxLRs unique to Augustus models:
11
```

This analysis didn't work for published genomes whose gff files did not match up
with fasta genome files - different names for contigs.
```
P. parisitica
Gff - Supercontig_2.1
Fa - >7000000185249020
P. infestans
Gff - Supercontig_1.4892
Fa - >supercont1.4921_dna:supercontig_supercontig:ASM14294v1:supercont1.4921:1:139:1
P. capsici
Gff - PHYCAscaffold_1
Fa - >PHYCAscaffold_1
P. sojae
Gff - scaffold_1
Fa - >scaffold_1
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
  # For P. capsici & P. sojae
  for GeneGff in $PcapPubGff $PsojPubGff; do
    echo "$GeneGff"
    Strain=$(echo "$GeneGff" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$GeneGff" | rev | cut -f4 -d '/' | rev)
    CrnDir=$(ls -d analysis/CRN_effectors/hmmer_CRN/$Species/$Strain)
    Source="pub"
    CRN_hmm_fa=$(ls analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_hmmer_out.fa)
    CRN_hmm_txt=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_hmmer.txt
    CRN_hmm_gff=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_hmmer.gff
    echo "$Species - $Strain"
    cat $CRN_hmm_fa | grep '>' | cut -f1 | tr -d '>' | cut -f4 -d '|' > $CRN_hmm_txt
    cat $GeneGff | grep -w -f $CRN_hmm_txt > $CRN_hmm_gff
  done
  # For P. infestans & P. parisitica
  for GeneGff in $PinfPubGff $PparPubGff; do
    echo "$GeneGff"
    Strain=$(echo "$GeneGff" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$GeneGff" | rev | cut -f4 -d '/' | rev)
    CrnDir=$(ls -d analysis/CRN_effectors/hmmer_CRN/$Species/$Strain)
    Source="pub"
    CRN_hmm_fa=$(ls analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_hmmer_out.fa)
    CRN_hmm_txt=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_hmmer.txt
    CRN_hmm_gff=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_hmmer.gff
    echo "$Species - $Strain"
    cat $CRN_hmm_fa | grep '>' | cut -f1 | tr -d '>' | cut -f4 -d '|' > $CRN_hmm_txt
    cat $GeneGff | grep -w -f $CRN_hmm_txt > $CRN_hmm_gff
  done
```

```bash
for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/P.infestans/* | grep -v -e '67593' -e 'masked'); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$MergeDir/"$Strain"_pub_CRN_hmmer.gff
if [ $Species == "P.cactorum" ]; then
  AugGff=$MergeDir/"$Strain"_Aug_CRN_hmmer.gff3
fi
ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer.gff3
if [ $Species == P.infestans ]; then
  cat $ORFGff | sed 's/^supercont/Supercontig_/' | sed -e 's/_dna.*\tCRN_HMM/\tCRN_HMM/' > $MergeDir/"$Strain"_CRN_merged_hmmer_mod.gff3
  ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer_mod.gff3
fi
ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_CRN_hmmer.bed
AugInORFs=$MergeDir/"$Strain"_AugInORFs_CRN_hmmer.bed
ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_CRN_hmmer.bed
AugUniq=$MergeDir/"$Strain"_Aug_Uniq_CRN_hmmer.bed
bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
# bedtools intersect -wb -u -a $ORFGff -b $AugGff > $AugInORFs
bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq
echo "$Species - $Strain"
if [ $Strain == "10300" ]; then
echo "The number of ORF CRNs overlapping Augustus CRNs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus CRNs overlapping ORF CRNs:"
cat $AugInORFs | grep -w 'gene' | wc -l
echo "The number of CRNs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of CRNs unique to Augustus models:"
cat $AugUniq | grep -w 'gene' | wc -l
elif [ $Strain == "T30-4" ]; then
echo "The number of ORF CRNs overlapping Augustus CRNs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus CRNs overlapping ORF CRNs:"
cat $AugInORFs | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | wc -l
echo "The number of CRNs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of CRNs unique to Augustus models:"
cat $AugInORFs | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq > Augtmp.txt
cat $AugUniq | grep -w 'exon' | rev | cut -f2 -d ':' | rev | sort | uniq | grep -v -f Augtmp.txt | wc -l
elif [ $Strain == "310" ]; then
echo "The number of ORF CRNs overlapping Augustus CRNs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus CRNs overlapping ORF CRNs:"
cat $AugInORFs | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | wc -l
echo "The number of CRNs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of CRNs unique to Augustus models:"
cat $AugUniq | grep -w 'exon' | cut -f9 | cut -f2 -d ';' | sort | uniq | wc -l
elif [ $Strain == "LT1534" ]; then
echo "The number of ORF CRNs overlapping Augustus CRNs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus CRNs overlapping ORF CRNs:"
cat $AugInORFs | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | wc -l
echo "The number of CRNs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of CRNs unique to Augustus models:"
cat $AugInORFs  | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | sed -e 's/^ //g' > Augtmp.txt
cat $AugUniq | grep -w 'exon' | grep 'transcriptId' | rev | cut -f1 -d ';' | rev | sort | uniq | grep -v -f Augtmp.txt | wc -l
echo ""
elif [ $Strain == "P6497" ]; then
echo "The number of ORF CRNs overlapping Augustus CRNs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus CRNs overlapping ORF CRNs:"
cat $AugInORFs | grep -w 'gene' | wc -l
echo "The number of CRNs unique to ORF models:"
cat $ORFsUniq | grep -w 'gene' | wc -l
echo "The number of CRNs unique to Augustus models:"
cat $AugUniq | grep -w 'gene' | wc -l
fi
done
```


```
P.cactorum - 10300
The number of ORF CRNs overlapping Augustus CRNs:
85
The number of Augustus CRNs overlapping ORF CRNs:
85
The number of CRNs unique to ORF models:
30
The number of CRNs unique to Augustus models:
7
P.capsici - LT1534
The number of ORF CRNs overlapping Augustus CRNs:
101
The number of Augustus CRNs overlapping ORF CRNs:
98
The number of CRNs unique to ORF models:
80
The number of CRNs unique to Augustus models:
8
```
