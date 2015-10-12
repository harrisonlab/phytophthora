# For a comparison between P.cac 10300 and P.inf T30-4, P.par 310, P.ram 163428 & P.soj 67593


```bash
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  IsolateAbrv=Pcac_Pinf_Ppar_Pram_Psoj
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
  Fasta_file=gene_pred/augustus_unmasked/P.cactorum/10300/10300_augustus_preds.aa
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

### for P.ram T3163428
```bash
  Taxon_code=Pram
  Fasta_file=assembly/external_group/P.ramorum/164328/pep/Phytophthora_ramorum.ASM14973v1.26.pep.all.fa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.soj 67593
```bash
  Taxon_code=Psoj
  Fasta_file=assembly/external_group/P.sojae/67593/pep/Phytophthora_sojae.ASM14975v1.26.pep.all.fa
  Id_field=1
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
    Jobs=$(qstat | grep 'blast_500' | wc -l)
    while [ $Jobs -gt 32 ]; do
      sleep 10
      printf "."
      Jobs=$(qstat | grep 'blast_500' | wc -l)
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

The following additional information was also provided. THe format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs

```
  [1] "Pcac (8580)"
  [1] 1529
  [1] 63
  [1] "Pinf (7969)"
  [1] 724
  [1] 157
  [1] "Ppar (8518)"
  [1] 822
  [1] 114
  [1] "Pram (6682)"
  [1] 115
  [1] 57
  [1] "Psoj (7461)"
  [1] 645
  [1] 152
```

# Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.

### P. cactorum unique RxLR families

P. cactorum strain 10300 RxLR genes were parsed to the same format as the gene
names used in the analysis:

```bash
  RxLR_Names_10300=analysis/sigP_rxlr/P.cactorum/10300/10300_aug_RxLR_EER_regex.txt
  RxLR_Dir=analysis/orthology/orthomcl/Pcac_Pinf_Pram_Psoj/Pcac_RxLR
  Orthogroups=analysis/orthology/orthomcl/Pcac_Pinf_Pram_Psoj/Pcac_Pinf_Pram_Psoj_orthogroups.txt
  RxLR_ID_10300=$RxLR_Dir/10300_aug_RxLR_EER_IDs.txt
  mkdir -p $RxLR_Dir
  cat $RxLR_Names_10300 | sed 's/g/Pcac|g/g' > $RxLR_ID_10300
```

Ortholog groups containing RxLR proteins were identified using the following
commands:
```bash
  RxLR_Orthogroup_10300=$RxLR_Dir/Pcac_RxLR_Orthogroups.txt
  cat $Orthogroups | grep -w -f $RxLR_ID_10300 > $RxLR_Orthogroup_10300
  RxLR_Orthogroup_hits_10300=$RxLR_Dir/Pcac_RxLR_Orthogroups_hits.txt
  cat $Orthogroups | grep -o -w -f $RxLR_ID_10300 > $RxLR_Orthogroup_hits_10300
```

All 61 predicted RxLRs were found to have orthologs in other taxa. The 61 RxLRs
were distributed through 41 orthogroups.

Orthogroup 8 was a large gene family of 298 genes. This included 46 P. cactorum
gene including 11 P. cactorum RxLRs.
```bash
  cat $Orthogroups | grep -w 'orthogroup8' | sed 's/ /\n/g' | sort | wc -l
  Orthogroup8_Pcac_ID=$RxLR_Dir/Pcac_RxLR_Orthogroups8_IDs.txt
  cat $Orthogroups | grep -w 'orthogroup8' | sed 's/ /\n/g' | sort | grep 'Pcac' > $Orthogroup8_Pcac_ID.txt
  cat $Orthogroups | grep -w 'orthogroup8' | sed 's/ /\n/g' | sort | cut -f1 -d'|' | uniq -c
```
```
  1 orthogroup8:
  46 Pcac
  100 Pinf
  73 Pram
  74 Psoj
```
