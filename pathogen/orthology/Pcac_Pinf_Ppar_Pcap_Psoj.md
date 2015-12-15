# For a comparison between 3 Phytophthora clade 1 isolates (P.cac 10300, P.inf T30-4, P.par 310), a clade 2 isolate (P.cap LT1534) & a clade 7 isolate (P.soj 67593)


```bash
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  IsolateAbrv=Pcac_Pinf_Ppar_Pcap_Psoj
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## Format fasta files

### for P.cac 10300
```bash
  Taxon_code=Pcac
  Fasta_file=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.inf T30-4
```bash
  Taxon_code=Pinf
  Fasta_file=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.par 310
```bash
  Taxon_code=Ppar
  Fasta_file=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins.pep.all.fa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.cap LT1534
```bash
  Taxon_code=Pcap
  Fasta_file=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins.fasta
  Id_field=3
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.soj 67593
```bash
  Taxon_code=Psoj
  Fasta_file=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401.aa.fasta
  Id_field=3
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


## Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## Perform an all-vs-all blast of the proteins

```bash
  BlastDB=$WorkDir/blastall/$IsolateAbrv.db

  makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
  BlastOut=$WorkDir/all-vs-all_results.tsv
  mkdir -p $WorkDir/splitfiles

  SplitDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  $SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $WorkDir/splitfiles --out_base goodProteins

  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology  
  for File in $(find $WorkDir/splitfiles); do
    Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
      sleep 10
      printf "."
      Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
    done
    printf "\n"
    echo $File
    BlastOut=$(echo $File | sed 's/.fa/.tab/g')
    qsub $ProgDir/blast_500.sh $BlastDB $File $BlastOut
  done
```

## Merge the all-vs-all blast results  
```bash  
  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts
```

## Plot venn diagrams:

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/venn_diagrams
  $ProgDir/ven_diag_5_way.R --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs

```
  [1] "Pcac (8814)"
  [1] 567
  [1] 118
  [1] "Pcap (7646)"
  [1] 333
  [1] 52
  [1] "Pinf (8335)"
  [1] 562
  [1] 100
  [1] "Ppar (8987)"
  [1] 695
  [1] 79
  [1] "Psoj (9156)"
  [1] 883
  [1] 388
```

# Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.


### P. cactotum unique gene families

The genes unique to P.cactorum were identified within the orthology analysis.

First variables were set:
```bash
  WorkDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj
  PcacUniqDir=$WorkDir/Pcac_unique
  Orthogroups=$WorkDir/Pcac_Pinf_Ppar_Pcap_Psoj_orthogroups.txt
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  Braker_genes=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
  Uniq_Pcac_groups=$PcacUniqDir/Pcac_uniq_orthogroups.txt
  mkdir -p $PcacUniqDir
```

Orthologroups only containing P.cactorum 10300 genes were extracted:

```bash
  cat $Orthogroups | grep -v 'Pinf' | grep -v 'Ppar' | grep -v 'Pcap' | grep -v 'Psoj' > $Uniq_Pcac_groups
  echo "The number of orthogroups unique to P. cactorum are:"
  cat $Uniq_Pcac_groups | wc -l
  echo "The following number genes are contained in these orthogorups:"
  cat $Uniq_Pcac_groups | grep -o 'Pcac|' | wc -l  
```

```
  The number of orthogroups unique to P. cactorum are:
  126
  The following number genes are contained in these orthogorups:
  347
```

### P. cactorum unique RxLR families

P. cactorum strain 10300 RxLR genes were parsed to the same format as the gene
names used in the analysis:

```bash
  RxLR_Names_10300=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_Aug_RxLR_EER_motif_hmm_headers.txt
  WorkDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj
  RxLR_Dir=$WorkDir/Pcac_RxLR
  Orthogroups=$WorkDir/Pcac_Pinf_Ppar_Pcap_Psoj_orthogroups.txt
  RxLR_ID_10300=$RxLR_Dir/10300_aug_RxLR_EER_IDs.txt
  mkdir -p $RxLR_Dir
  cat $RxLR_Names_10300 | sed 's/g/Pcac|g/g' > $RxLR_ID_10300
```

Ortholog groups containing RxLR proteins were identified using the following
commands:
```bash
  echo "The number of RxLRs searched for is:"
  cat $RxLR_ID_10300 | wc -l
  echo "Of these, the following number were found in orthogroups:"
  RxLR_Orthogroup_hits_10300=$RxLR_Dir/Pcac_RxLR_Orthogroups_hits.txt
  cat $Orthogroups | grep -o -w -f $RxLR_ID_10300 > $RxLR_Orthogroup_hits_10300
  cat $RxLR_Orthogroup_hits_10300 | wc -l
  echo "These were distributed through the following number of Orthogroups:"
  RxLR_Orthogroup_10300=$RxLR_Dir/Pcac_RxLR_Orthogroups.txt
  cat $Orthogroups | grep -w -f $RxLR_ID_10300 > $RxLR_Orthogroup_10300
  cat $RxLR_Orthogroup_10300 | wc -l
  echo "The following RxLRs were found in Pcac unique orthogroups:"
  RxLR_Pcac_uniq_groups=$RxLR_Dir/Pcac_RxLR_Orthogroups_hits.txt
  cat $RxLR_Orthogroup_10300 | grep -v 'Pinf' | grep -v 'Ppar' | grep -v 'Pcap' | grep -v 'Psoj' > $RxLR_Pcac_uniq_groups
  cat $RxLR_Pcac_uniq_groups | wc -l
  echo "The following RxLRs were found in Group1 unique orthogroups:"
  RxLR_Group1_uniq_groups=$RxLR_Dir/Group1_RxLR_Orthogroups_hits.txt
  cat $RxLR_Orthogroup_10300 | grep -v 'Pcap' | grep -v 'Psoj' > $RxLR_Group1_uniq_groups
  cat $RxLR_Group1_uniq_groups | wc -l
```

```
  The number of RxLRs searched for is:
  145
  Of these, the following number were found in orthogroups:
  145
  These were distributed through the following number of Orthogroups:
  79
  The following RxLRs were found in Pcac unique orthogroups:
  2
  The following RxLRs were found in Group1 unique orthogroups:
  17
```

The P.cactorum RxLR genes that were not found in orthogroups were identified:

```bash
  RxLR_10300_uniq=$RxLR_Dir/Pcac_unique_RxLRs.txt
  cat $RxLR_ID_10300 | grep -v -w -f $RxLR_Orthogroup_hits_10300 | tr -d 'Pcac|' > $RxLR_10300_uniq
  echo "The number of P.cac 10300 unique RxLRs are:"
  cat $RxLR_10300_uniq | wc -l
  RxLR_Seq_10300=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_Aug_RxLR_EER_motif_hmm.fa
  Braker_genes=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
  RxLR_10300_uniq_fa=$RxLR_Dir/Pcac_unique_RxLRs.fa
  cat $Braker_genes | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -w -A1 -f $RxLR_10300_uniq | grep -E -v '^--' > $RxLR_10300_uniq_fa
```

```
  The number of P.cac 10300 unique RxLRs are:
  0
```

```bash
  mkdir -p analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/orthogroups_fasta_Pcac_RxLR
  ~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL/orthoMCLgroups2fasta.py --orthogroups analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/Pcac_RxLR_Orthogroups.txt --fasta analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/goodProteins/goodProteins.fasta --out_dir analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/orthogroups_fasta_Pcac_RxLR
  mkdir -p analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/orthogroups_fasta_clade1_RxLR
  ~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL/orthoMCLgroups2fasta.py --orthogroups analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/Group1_RxLR_Orthogroups_hits.txt --fasta analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/goodProteins/goodProteins.fasta --out_dir analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/orthogroups_fasta_clade1_RxLR
```






<!--
### P. cactorum unique RxLR families

P. cactorum strain 10300 RxLR genes were parsed to the same format as the gene
names used in the analysis:

```bash
  RxLR_Motif_Names_10300=analysis/RxLR_effectors/RxLR_EER_regex_finder/P.cactorum/10300/10300_Aug_RxLR_EER_regex.txt
  RxLR_HMM_Names_10300=analysis/RxLR_effectors/hmmer_RxLR/P.cactorum/10300/10300_Aug_RxLR_hmmer_headers.txt
  RxLR_Dir=$WorkDir/Pcac_RxLR
  Orthogroups=$WorkDir/"$IsolateAbrv"_orthogroups.txt
  RxLR_ID_10300=$RxLR_Dir/10300_Aug_RxLR_EER_IDs.txt
  mkdir -p $RxLR_Dir
  cat $RxLR_Motif_Names_10300 $RxLR_HMM_Names_10300 | sort | uniq | sed 's/g/Pcac|g/g' > $RxLR_ID_10300
```

Ortholog groups containing RxLR proteins were identified using the following
commands:
```bash
  RxLR_Orthogroup_10300=$RxLR_Dir/Pcac_RxLR_Orthogroups.txt
  cat $Orthogroups | grep -w -f $RxLR_ID_10300 > $RxLR_Orthogroup_10300
  RxLR_Orthogroup_hits_10300=$RxLR_Dir/Pcac_RxLR_Orthogroups_hits.txt
  cat $Orthogroups | grep -o -w -f $RxLR_ID_10300 > $RxLR_Orthogroup_hits_10300
```

None of the 145 predicted RxLRs were grouped into orthogroups, grouped with genes
from other taxa or with inparalogs within the genome. The 145 putative RxLRs
distributed through 79 orthogroups.

Orthogroup 1 was a large gene family of 1824 genes. This included 359 P. cactorum
gene including 2 P. cactorum RxLRs.
```bash
  cat $Orthogroups | grep -w 'orthogroup1' | sed 's/ /\n/g' | sort | grep -v 'orthogroup1' | wc -l
  Orthogroup1_Pcac_ID=$RxLR_Dir/Pcac_RxLR_Orthogroup1_IDs.txt
  cat $Orthogroups | grep -w 'orthogroup1' | sed 's/ /\n/g' | sort | grep 'Pcac' > $Orthogroup1_Pcac_ID.txt
  cat $Orthogroups | grep -w 'orthogroup1' | sed 's/ /\n/g' | sort | grep -v 'orthogroup1' | cut -f1 -d'|' | uniq -c
  cat $Orthogroup1_Pcac_ID.txt $RxLR_Orthogroup_hits_10300 | sed -r 's/\.t.//g' | sort | uniq | wc -l
  # cat assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa | grep -i 'rxlr' | cut -f1 -d ' ' | tr -d '>' > $RxLR_Dir/Pinf_RxLR_pub_RxLR_ID.txt
  # Pinf_RxLR_ID=analysis/RxLR_effectors/RxLR_EER_regex_finder/P.infestans/T30-4/T30-4_pub_RxLR_EER_regex.txt
  # cat $Pinf_RxLR_ID | cut -f1 -d ' ' > $RxLR_Dir/Pinf_RxLR_EER_motif_IDs.txt
  # cat $RxLR_Dir/Pinf_RxLR_EER_motif_IDs.txt $RxLR_Dir/Pinf_RxLR_Orthogroup1_IDs.txt | sort | uniq -d | wc -l
  # cat $Orthogroups | grep -w 'orthogroup1' | sed 's/ /\n/g' | sort | grep 'Pinf' | sed 's/Pinf|//g' > $RxLR_Dir/Pinf_RxLR_Orthogroup1_IDs.txt
  # cat $RxLR_Dir/Pinf_RxLR_pub_RxLR_ID.txt $RxLR_Dir/Pinf_RxLR_Orthogroup1_IDs.txt | sort | uniq -d | wc -l
```
```
  359 Pcac
  376 Pcap
  357 Pinf
  380 Ppar
  352 Psoj
```

The distibution of publihsed P.infestans RxLR genes throughout the ortholog groups was
identified. The 486 RxLR genes were present in 112 ortholog groups.

```bash
  Pinf_RxLR_ID=analysis/RxLR_effectors/RxLR_EER_regex_finder/P.infestans/T30-4/T30-4_pub_RxLR_EER_regex.txt
  Pinf_RxLR_Groups=analysis/RxLR_effectors/RxLR_EER_regex_finder/P.infestans/T30-4/T30-4_pub_RxLR_EER_regex_orthogroups.txt
  cat $Orthogroups | grep -w -f $RxLR_Dir/Pinf_RxLR_pub_RxLR_ID.txt > $Pinf_RxLR_Groups
  cat $Pinf_RxLR_Groups | sed 's/ /\n/g' | sort | grep -v 'orthogroup' | cut -f1 -d'|' | uniq -c
```

```bash
  245 Pcac
  159 Pcap
  559 Pinf
  420 Ppar
  172 Psoj
``` -->
