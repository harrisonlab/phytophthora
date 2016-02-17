P. capsici is reported to have 237 Crinklers according to Stam et al 2013. These
were identified in the published proteins:

## All Pcap predicted crinklers

```bash
  OutDir=analysis/Pcap_effectors
  mkdir -p $OutDir

  DownloadedCRNtab=$OutDir/Pcap_sup.tab1_CRN.txt
  PubCRN_fa=$OutDir/Pcap_CRN.fa
  cat $DownloadedCRNtab | sed -E 's/\r/\n/g' | tail -n+2 | cut -f1,15 | sed -E 's/^/>/g' | grep 'CRN' |  sed 's/\t/\n/g' > $PubCRN_fa
  cat $PubCRN_fa | grep '>' | wc -l
  echo "Within these crinklers, the following number contained X's in their sequence:"
  cat $PubCRN_fa | grep -v '>' | grep 'X' |  wc -l
```

36 of the 237 genes contained X's and may not represent valid gene models as
they carry X's within their sequence.


Published Crinklers were from isolate Pcap11, these crinklers were Blasted
against the LT1534 genome.

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  Assembly=assembly/external_group/P.capsici/LT1534/dna/Phyca11_unmasked_genomic_scaffolds.fasta
  Query=$PubCRN_fa
  qsub $ProgDir/blast_pipe.sh $Query protein $Assembly
```

Manual inspection of the BLAST results found that most Stam et al CRNs had 100%
alignment length in the top BLAST hit. Those that didn't typically had many X's
in the Crinkler sequence. The top BLAST hit of each query was found to have 100%
identity in the region it aligned to. A single query (CRN75_29) did not have a
blast hit against the genome. In this query 22/239 aa were not X's.

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	BlastHits=analysis/blast_homology/P.capsici/LT1534/LT1534_Pcap_CRN.fa_homologs.csv
	HitsGff=analysis/blast_homology/P.capsici/LT1534/LT1534_Pcap_CRN.fa_homologs.gff
	Column2=Stam_et_al_CRN
	NumHits=1
	$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
```

Blast hits of Stam et al CRNs from the Pcap11 genome were intersected with
predicted Crinklers from the LT1534 genome:

```bash
  PcapCRNGff=$HitsGff
  PcapPredCRN=analysis/CRN_effectors/hmmer_CRN/P.capsici/LT1534/LT1534_Total_CRN.gff
  PcapIntersectCRN=$OutDir/PcapCRN_intersected.gff
  echo "Of the 237 published CRNs, the number we identified were:"
  bedtools intersect -s -wo -a $PcapPredCRN -b $PcapCRNGff > $PcapIntersectCRN
  cat $PcapIntersectCRN | cut -f 18 | cut -f2 -d '"' | cut -f2 -d '=' | sort -V | uniq | wc -l
  cat $PcapIntersectCRN | cut -f 9 | rev | cut -f1 -d ';' | rev | sed 's/Name=//g' | sed -E 's/^ //g' | sort -V | uniq | wc -l
  bedtools intersect -s -v -a $PcapCRNGff -b $PcapPredCRN | cut -f9 | cut -f2 -d '"' | cut -f2 -d '=' | sort -V | less
```

Of the 237 Stam et al Crinklers, 142 of them intersect with the 175 crinklers
we predicted.


## Full length P.cap predicted crinklers

84 crinklers carried HVLVVVP and LFAK motifs and were termed full length
crinklers. These were extracted as well:

```bash
  echo "Full length crinklers were also extracted:"
  PubCRN_fa_FL=$OutDir/Pcap_CRN_full_lgth.fa
  cat $DownloadedCRNtab | sed -E 's/\r/\n/g' | tail -n+2 | cut -f1,7,15 | grep -w -i 'Y' | cut -f 1,3 | sed -E 's/^/>/g' | sed 's/\t/\n/g' > $PubCRN_fa_FL
  cat $PubCRN_fa_FL | grep '>' | wc -l
  echo "Within these crinklers, the following number contained X's in their sequence:"
  cat $PubCRN_fa_FL | grep -v '>' | grep 'X' |  wc -l
```

19 of the 84 proteins contained X's and may not represent valid gene models.

These full length Crinklers were Blasted against the LT1534 genome.

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  Assembly=assembly/external_group/P.capsici/LT1534/dna/Phyca11_unmasked_genomic_scaffolds.fasta
  Query=$PubCRN_fa_FL
  qsub $ProgDir/blast_pipe.sh $Query protein $Assembly
```

<!-- Manual inspection of the BLAST results found that most Stam et al CRNs had 100%
alignment length in the top BLAST hit. Those that didn't typically had many X's
in the Crinkler sequence. The top BLAST hit of each query was found to have 100%
identity in the region it aligned to. A single query (CRN75_29) did not have a
blast hit against the genome. In this query 22/239 aa were not X's. -->

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	BlastHits=analysis/blast_homology/P.capsici/LT1534/LT1534_Pcap_CRN_full_lgth.fa_homologs.csv
	HitsGff=analysis/blast_homology/P.capsici/LT1534/LT1534_Pcap_CRN_full_lgth.fa_homologs.gff
	Column2=Stam_et_al_CRN
	NumHits=1
	$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
```

Blast hits of Stam et al CRNs from the Pcap11 genome were intersected with
predicted Crinklers from the LT1534 genome:

```bash
  PcapCRNGff=$HitsGff
  PcapPredCRN=analysis/CRN_effectors/hmmer_CRN/P.capsici/LT1534/LT1534_Total_CRN.gff
  PcapIntersectCRN_FL=$OutDir/PcapCRN_intersected.gff
  echo "Of the 237 published CRNs, the number we identified were:"
  bedtools intersect -s -wo -a $PcapPredCRN -b $PcapCRNGff > $PcapIntersectCRN_FL
  cat $PcapIntersectCRN_FL | cut -f 18 | cut -f2 -d '"' | cut -f2 -d '=' | sort -V | uniq | wc -l
  cat $PcapIntersectCRN_FL | cut -f 9 | rev | cut -f1 -d ';' | rev | sed 's/Name=//g' | sed -E 's/^ //g' | sort -V | uniq | wc -l
```

Of the 83 full length Stam et al Crinklers, 71 of them intersect with the 175 crinklers
we predicted. The 13 that were missing were all 'X-rich'

<!--
  PcapProt=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
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
``` -->
