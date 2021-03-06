
This document details how hmm models were trained to LFLAK and DWL domains of
Phytophthora Crinklers.


Please find a Powerpoint presentation attached, which gives an overview of the
rational and process behind the analysis.


Alignments of Phytophthora Crinklers can be found in this directory:

```bash
  ProjDir=/home/groups/harrisonlab/project_files/idris
  HmmDir=$ProjDir/analysis/CRN_effectors/hmmer_models
  mkdir -p $HmmDir
  PhytCRN=$HmmDir/Pinf_Pram_Psoj_Pcap_CRN.fa
  PhytCRN_LFLAK=$HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.fa
  PhytCRN_DWL=$HmmDir/Pinf_Pram_Psoj_Pcap_DWL.fa
```

Alignments were made in Geneious. Inital alignments were made using the MAFFT
algorithm, before being adjusted by eye. Proteins were removed that were
considered to be psudogenes or that did not carry that particular domain.


## Alignment conversion

Alignments were exported from geneious in fasta format. These alignments were
converted into stockholm alignment format using the wesbite:
http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_stockholm.php

The alignments can be found at:
```bash
  PhytCRN_LFLAK_aln=$HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.stockholm
  PhytCRN_DWL_aln=$HmmDir/Pinf_Pram_Psoj_Pcap_DWL.stockholm
```


## Building hmm models

Hmm models were built for the LFAK and DWL domains using hmmer 3.0.

```bash
  LFLAK_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm
  hmmbuild -n Pinf_Pram_Psoj_Pcap_LFLAK --amino $LFLAK_hmm $PhytCRN_LFLAK_aln
  DWL_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm
  hmmbuild -n Pinf_Pram_Psoj_Pcap_DWL --amino $DWL_hmm $PhytCRN_DWL_aln
```

## Validating hmm models

### Validation within P. infestans published gene models:

Crinklers were identified among the published P. infestans proteins:
```bash
  PinfProts=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all_parsed.fa
  OutDir=analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4

  PinfProts_LFLAK=$OutDir/T30-4_pub_CRN_LFLAK_hmm.txt
  hmmsearch -T0 $LFLAK_hmm $PinfProts > $PinfProts_LFLAK
  cat $PinfProts_LFLAK | grep 'Initial search space'
  cat $PinfProts_LFLAK | grep 'number of targets reported over threshold'
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
  $ProgDir/hmmer2fasta.pl $PinfProts_LFLAK $PinfProts > $OutDir/T30-4_pub_CRN_LFLAK_hmm.fa

  PinfProts_DWL=$OutDir/T30-4_pub_CRN_DWL_hmm.txt
  hmmsearch -T0 $DWL_hmm $PinfProts > $PinfProts_DWL
  cat $PinfProts_DWL | grep 'Initial search space'
  cat $PinfProts_DWL | grep 'number of targets reported over threshold'
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
  $ProgDir/hmmer2fasta.pl $PinfProts_DWL $PinfProts > $OutDir/T30-4_pub_CRN_DWL_hmm.fa

  cat $OutDir/T30-4_pub_CRN_LFLAK_hmm.fa $OutDir/T30-4_pub_CRN_DWL_hmm.fa | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $OutDir/T30-4_pub_CRN_LFLAK_DWL.txt
  cat $OutDir/T30-4_pub_CRN_LFLAK_DWL.txt | grep '>' | wc -l
```

The LFLAK and DWL hmm models identified 193 and 168 proteins respectively. Of
these, 167 proteins were common to both.

Gff annotations were extracted for these putative crinklers:

```bash
  PinfPubGff=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts.gff3
  for GeneGff in $PinfPubGff; do
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


### Validation within P. infestans ORFs:

Crinklers were identified among the predicted P. infestans ORFs:
<!--
```bash
  PinfORFs=gene_pred/ORF_finder/P.infestans/T30-4/T30-4.aa_cat.fa

  PinfORFs_LFLAK=tmp/Pinf_ORF_CRN_LFLAK_hmm.txt
  hmmsearch -T0 $LFLAK_hmm $PinfORFs > $PinfORFs_LFLAK
  cat $PinfORFs_LFLAK | grep 'Initial search space'
  cat $PinfORFs_LFLAK | grep 'number of targets reported over threshold'
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
  $ProgDir/hmmer2fasta.pl $PinfORFs_LFLAK $PinfORFs > tmp/Pinf_ORF_CRN_LFLAK_hmm.fa

  PinfORFs_DWL=tmp/Pinf_ORF_CRN_DWL_hmm.txt
  hmmsearch -T0 $DWL_hmm $PinfORFs > $PinfORFs_DWL
  cat $PinfORFs_DWL | grep 'Initial search space'
  cat $PinfORFs_DWL | grep 'number of targets reported over threshold'
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
  $ProgDir/hmmer2fasta.pl $PinfORFs_DWL $PinfORFs > tmp/Pinf_ORF_CRN_DWL_hmm.fa

  CRN_ORF_list=tmp/Pinf_ORF_CRN_LFLAK_DWL.txt
  cat tmp/Pinf_ORF_CRN_LFLAK_hmm.fa tmp/Pinf_ORF_CRN_DWL_hmm.fa | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $CRN_ORF_list
  cat $CRN_ORF_list | wc -l
``` -->

The LFLAK and DWL hmm models identified 568 and 760 ORFs respectively. Of
these, 373 ORFs were common to both.
<!--
The 373 ORFs identified in both CRN hmm models were were extracted from gff
annotations:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
  Col2=atg_CRN
  GeneModels=gene_pred/ORF_finder/P.infestans/T30-4/T30-4_ORF_corrected.gff3
  # OutFile=analysis/CRN/P.infestans/T30-4/T30-4_ORF_LxLFLAK_HVLVVVP.gff3
  OutFile=tmp/Pinf_ORF_CRN_LFLAK_DWL.gff3
  $ProgDir/gene_list_to_gff.pl $CRN_ORF_list $GeneModels $Col2 Name > $OutFile
``` -->

```bash
  for Proteome in $(ls gene_pred/ORF_finder/P.infestans/T30-4/T30-4.aa_cat.fa); do
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
    Headers="$Strain"_CRN_hmmer_unmerged_headers.txt
    cat $OutDir/$HmmFastaDWL | grep '>' | grep -w -f $CommonHeaders | tr -d '>' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | tr -d '-' | sed 's/hmm_score/HMM_score/g' > $OutDir/$Headers
    # As we are dealing with JGI and Broad sequences, some headers need formatting:
    cat $OutDir/$Headers | sed 's/:/_a_/g' | sed 's/supercont1./supercont1_b_/g' | sed 's/Supercontig_2./Supercontig_c_/g' > tmp.txt
    # As we are dealing with JGI and Broad sequences, some features need formatting:
    ORF_Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/*_ORF_corrected.gff3)
    cat $ORF_Gff | sed 's/:/_a_/g' | sed 's/supercont1./supercont1_b_/g' | sed 's/Supercontig_2./Supercontig_c_/g' > tmp.gff
    # Gff features were extracted for each header
    CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl tmp.txt tmp.gff CRN_HMM Name > $CRN_unmerged_Gff
    # Gff features were merged based upon the DWL hmm score
    DbDir=analysis/databases/$Organism/$Strain
    mkdir -p $DbDir
    ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $CRN_unmerged_Gff --db $DbDir/CRN_ORF.db
    CRN_Merged_Gff=$OutDir/"$Strain"_CRN_merged_hmmer.gff3
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp $DbDir/CRN_ORF.db --id LFLAK_DWL_CRN --out $DbDir/CRN_ORF_merged.db --gff > $CRN_Merged_Gff
    # As we are dealing with JGI and Broad sequences, some features need formatting:
    sed -i 's/_a_/:/g' $CRN_Merged_Gff
    sed -i 's/supercont1_b_/supercont1./g' $CRN_Merged_Gff
    sed -i 's/Supercontig_c_/Supercontig_2./g' $CRN_Merged_Gff
    # Final results are reported:
    echo "Number of CRN ORFs after merging:"
    cat $CRN_Merged_Gff | grep 'gene' | wc -l
    # Temporary files containing the reformatted headers and features were deleted
    rm tmp.txt
    rm tmp.gff
  done
```

### Intersecting Crinklers predicted from gene models and ORFs

bedtools intersect was used to identify the number of ORF CRNs intersecting CRNs
from published gene models.

```bash
  for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4 ); do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    AugGff=$MergeDir/"$Strain"_pub_CRN_LFLAK_DWL.gff
    ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer.gff3
    cat $ORFGff | sed 's/^supercont/Supercontig_/' | sed -e 's/_dna.*\tCRN_HMM/\tCRN_HMM/' > $MergeDir/"$Strain"_CRN_merged_hmmer_mod.gff3
    ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer_mod.gff3
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
    cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalCRNsTxt > $TotalCRNsGff
    echo "The total number of CRNs are:"
    cat $TotalCRNsGff | grep -o -E -e 'Name=.*$' -e 'Parent=.*T.$' | sed -e 's/^Name=//g' | sed 's/Parent=//g' | sed 's/";$//g' | sort | uniq > $TotalCRNsHeaders
    cat $TotalCRNsHeaders | wc -l
  done
```
```
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
```

Fasta sequences for CRNs were extracted

```bash
  PinfFa=$(ls assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all_parsed.fa)
  for AugFa in $PinfFa; do
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



265 Crinklers were identified in the genome. This may still represent a false
 positives being predicted. The Crinklers predicted in Haas 2009 were compared
 to the Crinklers predicted in this study.


 ```bash
   OutDir=analysis/Pinf_effectors
   PinfCRNGff=$OutDir/PinfCRN.gff
   PinfPredCRN=analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_Total_CRN.gff
   PinfIntersectCRN=$OutDir/Pinf_ORF_CRN_LFLAK_DWL_intersected.gff
   PinfIntersectCRN_headers=$OutDir/Pinf_ORF_CRN_LFLAK_DWL_intersected.txt
   echo "Of the 194 published CRNs, the number we identified were:"
   bedtools intersect -s -wo -a $PinfPredCRN -b $PinfCRNGff > $PinfIntersectCRN
   cat $PinfIntersectCRN | grep 'exon' | cut -f2 -d '"' | sort | uniq > $PinfIntersectCRN_headers
   cat $PinfIntersectCRN_headers | wc -l
   echo "The number missing were:"
   bedtools intersect -s -v -a $PinfCRNGff -b $PinfPredCRN | grep -w 'exon' | cut -f9 | cut -f2 -d '"' | sort | uniq | grep -v -f $PinfIntersectCRN_headers > $OutDir/T30-4_pub_CRN_missing.txt
   cat $OutDir/T30-4_pub_CRN_missing.txt | wc -l
```

```bash
  echo "Of the 159 published CRN genes that have BLAST homologs to our CRNs, the following number were from published gene models:"
  cat $PinfIntersectCRN | grep -v 'CRN_HMM' | grep 'exon' | cut -f4 -d '"' | sort | uniq | wc -l
  echo "and the following number were from ORFs:"
  cat $PinfIntersectCRN | grep 'CRN_HMM' | grep 'exon' | cut -f4 -d '"' | sort | uniq | wc -l
```

The 35 Crinklers that were predicted in published models, but missing from CRNs
in this study were identified and fasta accessions extracted.
Study of these proteins in geneious showed that one of them contained an LQLFLAK
and none of them contained a canonical DWL domain (HVLVVVP). Some LIVQVP
domain or DWL domain (HVLVVVP sequence). The closest was the presence of an
LFVAQ sequence in one and RQVVVP in five sequences. 15 of the sequences had a
very similar start.


```bash
  PinfProt=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
  PinfProtMod=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all_mod.fa
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
  $ProgDir/unwrap_fasta.py --inp_fasta $PinfProt > $PinfProtMod

  echo "The following Crinklers were not found:"
  cat $OutDir/T30-4_pub_CRN_missing.txt | wc -l
  cat $PinfProtMod | grep -w -A1 -f $OutDir/T30-4_pub_CRN_missing.txt | grep -v -E '^--' > $OutDir/T30-4_pub_CRN_LFLAK_DWL_missing.fa
```
