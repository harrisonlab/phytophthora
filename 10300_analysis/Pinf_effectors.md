
P. infestans is known to have 194 Crinklers and 486 RxLRs. These were identified in the
published proteins:

```bash
  OutDir=analysis/Pinf_effectors
  mkdir -p $OutDir
  PinfProt=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
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
  echo "Of the 175 CRN genes identified in this study the following number were form published gene models:"
  cat $PinfIntersectCRN | grep -v 'CRN_HMM' | grep 'exon' | cut -f4 -d '"' | sort | uniq | wc -l
  echo "and the following number were from ORFs:"
  cat $PinfIntersectCRN | grep -v 'CRN_HMM' | grep 'exon' | cut -f4 -d '"' | sort | uniq | wc -l
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
