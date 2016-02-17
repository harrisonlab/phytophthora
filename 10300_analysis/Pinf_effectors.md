
P. infestans is known to have 194 Crinklers and 486 RxLRs. These were identified in the
published proteins:

```bash
  OutDir=analysis/Pinf_effectors
  mkdir -p $OutDir
  PinfProt=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
  PinfProtMod=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all_mod.fa
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
  $ProgDir/unwrap_fasta.py --inp_fasta $PinfProt > $PinfProtMod
  PinfCRN=$OutDir/PinfCRN.txt
  PinfCRNFa=$OutDir/PinfCRN.fa
  cat $PinfProt | grep '>' | grep 'CRN' | cut -f1  -d ' ' | tr -d '>' > $PinfCRN
  cat $PinfProtMod | grep -A1 -E '^>.*CRN' | grep -v -E '^--' > $PinfCRNFa
  echo "Number of CRN effectors extracted:"
  cat $PinfCRN | wc -l
  cat $PinfCRNFa | grep '>' | wc -l
  PinfRxLR=$OutDir/PinfRxLR.txt
  PinfRxLRFa=$OutDir/PinfRxLR.fa
  cat $PinfProt | grep '>' | grep -i 'rxlr' | cut -f1  -d ' ' | tr -d '>' > $PinfRxLR
  cat $PinfProtMod | grep -i -A1 -E '>.*rxlr' | grep -v -E '^--' > $PinfRxLRFa
  echo "Number of RxLR effectors extracted:"
  cat $PinfRxLR | wc -l
  cat $PinfRxLRFa | grep '>' | wc -l
```

The gff features for these proteins were extracted


```bash
  PinfGff=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts.gtf
  PinfGff_mod=$OutDir/PinfGff_mod.gff
  cat $PinfGff | sed 's/ PI_T30-4_FINAL_CALLGENES_4//g' > $PinfGff_mod
  PinfCRNGff=$OutDir/PinfCRN.gff
  cat $PinfGff_mod | grep -w -f $PinfCRN > $PinfCRNGff
  echo "Number of CRN effectors extracted:"
  cat $PinfCRNGff | grep 'exon' | cut -f4 -d '"' | sort | uniq | wc -l
  PinfRxLRGff=$OutDir/PinfRxLR.gff
  cat $PinfGff_mod | grep -w -f $PinfRxLR > $PinfRxLRGff
  echo "Number of RxLR effectors extracted:"
  cat $PinfRxLRGff | grep 'exon' | cut -f4 -d '"' | sort | uniq | wc -l
```

The location previously identified CRNs and RxLRs were intersected with the
CRNs and RxLRs from this study.

```bash
  PinfPredCRN=analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_Total_CRN.gff
  PinfIntersectCRN=$OutDir/PinfCRN_intersected.gff
  echo "Of the 194 published CRNs, the number we identified were:"
  bedtools intersect -s -wo -a $PinfPredCRN -b $PinfCRNGff > $PinfIntersectCRN
  cat $PinfIntersectCRN  | grep 'exon' | cut -f4 -d '"' | sort | uniq | wc -l
  bedtools intersect -s -v -a $PinfCRNGff -b $PinfPredCRN | grep -w 'exon' | cut -f9 | cut -f2 -d '"' > $OutDir/T30-4_pub_CRN_missing.txt

  PinfPredRxLR=analysis/RxLR_effectors/combined_evidence/P.infestans/T30-4/T30-4_Total_RxLR_EER_motif_hmm.gff
  PinfIntersectRxLR=$OutDir/PinfRxLR_intersected.gff
  echo "Of the 486 published RxLRs, the number we identified were:"
  bedtools intersect -s -wb -a $PinfPredRxLR -b $PinfRxLRGff > $PinfIntersectRxLR
  cat $PinfIntersectRxLR | grep 'exon' | cut -f4 -d '"' | sort | uniq | wc -l
```


It seemed that a large number of false positives were being predicted for
Crinklers. The number of false positives in published gene models was
investigated.

```bash
  echo "Of the 175 published CRN genes that ahve BLAST homologs to our CRNs, the following number were from published gene models:"
  cat $PinfIntersectCRN | grep -v 'CRN_HMM' | grep 'exon' | cut -f4 -d '"' | sort | uniq | wc -l
  echo "and the following number were from ORFs:"
  cat $PinfIntersectCRN | grep -v 'CRN_HMM' | grep 'exon' | cut -f4 -d '"' | sort | uniq | wc -l
```

The Crinklers that were predicted in published models, but missing from CRNs
in this study were identified and fasta accessions extracted.
Study of these proteins in geneious showed that none of them contained a LFLAK
domain or DWL domain (HVLVVVP sequence). The closest was the presence of an
LFVAQ sequence in one and RQVVVP in five sequences. 15 of the sequences had a
very similar start.

CRNs that were predicted but were not present in P. infestans were identified
and their protein sequence extracted.

```bash
  echo "The following Crinklers were not found:"
  bedtools intersect -s -v -a $PinfCRNGff -b $PinfPredCRN | grep -w 'exon' | cut -f9 | cut -f2 -d '"' > $OutDir/T30-4_pub_CRN_missing.txt
  cat $OutDir/T30-4_pub_CRN_missing.txt | wc -l
  cat $PinfProtMod | grep -w -A1 -f $OutDir/T30-4_pub_CRN_missing.txt > $OutDir/T30-4_pub_CRN_missing.fa

  PinfCRN_pred_fa=analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_Total_CRN.fa
  bedtools intersect -s -v -a $PinfPredCRN -b $PinfCRNGff | grep -w -e 'exon' -e 'transcript' | grep -v -E 'PI_T30-4_FINAL_CALLGENES_4.*transcript' | cut -f9 | sed -E 's/.*Parent=//g' | sed -E 's/.*Name=//g' | sort | uniq > $OutDir/T30-4_total_CRN_false_pos.txt
  cat $OutDir/T30-4_total_CRN_false_pos.txt | wc -l
  cat $PinfCRN_pred_fa | grep -w -A1 -f $OutDir/T30-4_total_CRN_false_pos.txt | grep -v -E '^--' > $OutDir/T30-4_total_CRN_false_pos.fa
  cat $OutDir/T30-4_total_CRN_false_pos.txt | wc -l
```





```bash
cat $PinfCRNFa | grep -E -e 'VLV.VP|VLVVVP|VDDFLP|VDPRLP|VPVYRP|VTLKYP|VLVELP|VCLDGP|VLVALP|VVDQTP|VEQLLP|VLVVAP|VKIGLP|VEGVGP|VDAVKP|VLVDCP|VKTKIP|VILGNP|VEGRGP' | wc -l
cat $PinfCRNFa |  grep -o -E 'VLV..P' | sort | uniq -c | sort -r -n | less
cat $PinfCRNFa | grep -E 'VLV..P' | sort | wc -l
cat analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_Total_CRN.fa | grep -o -E 'VLV..P' | sort | uniq -c | sort -r -n | less
cat analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_Total_CRN.fa | grep -E 'VLV..P' | sort | wc -l
```

Published genes that did not contain a VL...P motif were identified

```bash
  PinfCRN_no_HVLVVVP=$OutDir/PinfCRN_no_HVLVVVP.fa
  echo "The number of published CRNs not containing a VLV..P motif were:"
  cat $PinfCRNFa |  grep -v -E 'VLV..P' | grep -B1 -E '^M' | grep '>' | wc -l
  cat $PinfCRNFa |  grep -v -E 'VLV..P' | grep -B1 -E '^M' | grep -v -E '^--' > $PinfCRN_no_HVLVVVP
```


<!--

A hmm model was built for the DWL domain of Phytophthora crinklers.

The DWL domain sequences published by Haas et 2009 was used to create the model.
An alignment was made using MAFFT in geneious and exported as a fasta.

This was converted to stockholm format using the wesite:
http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_stockholm.php


```bash
  DWL_Align=/home/armita/git_repos/emr_repos/tools/pathogen/crinkler/Haas_et_al_2009_hmm/C-domain/DWL.stockholm
  hmmbuild -n DWL_ex_Haas_2009 --amino DWL_ex_Haas_2009.hmm $DWL_Align
```


```bash

hmmsearch -T0 DWL_ex_Haas_2009.hmm analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_Total_CRN.fa > tmp_Pinf_CRN_DWL.txt
cat tmp_Pinf_CRN_DWL.txt | grep 'Initial search space'
cat tmp_Pinf_CRN_DWL.txt | grep 'number of targets reported over threshold'
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
$ProgDir/hmmer2fasta.pl tmp_Pinf_CRN_DWL.txt analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_Total_CRN.fa > tmp_Pinf_CRN_DWL.fa

hmmsearch -T0 DWL_ex_Haas_2009.hmm analysis/CRN_effectors/hmmer_CRN/P.capsici/LT1534/LT1534_Total_CRN.fa > tmp_Pcap_CRN_DWL.txt
cat tmp_Pcap_CRN_DWL.txt | grep 'Initial search space'
cat tmp_Pcap_CRN_DWL.txt | grep 'number of targets reported over threshold'
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
$ProgDir/hmmer2fasta.pl tmp_Pcap_CRN_DWL.txt analysis/CRN_effectors/hmmer_CRN/P.capsici/LT1534/LT1534_Total_CRN.fa > tmp_Pcap_CRN_DWL.fa

``` -->
