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

## 1) Format fasta files

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


## 2) Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## 3.1) Perform an all-vs-all blast of the proteins

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

## 3.2) Merge the all-vs-all blast results  
```bash  
  MergeHits=$WorkDir/"$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## 4) Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```

## 5) Plot venn diagrams:

```bash
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
  $ProgDir/Pcac_Pinf_Ppar_Pcap_Psoj_venn.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name
total number of orthogroups
number of unique singleton genes
number of unique groups of inparalogs

```
  [1] "Pcac"
  [1] 12631
  [1] 1564
  [1] 413
  [1] "Pcap"
  [1] 11158
  [1] 1147
  [1] 375
  [1] "Psoj"
  [1] 13822
  [1] 2688
  [1] 1002
  [1] "Pinf"
  [1] 11206
  [1] 1230
  [1] 258
  [1] "Ppar"
  [1] 12849
  [1] 1755
  [1] 315
```

# 6) Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.


### 6.1 ) P. cactotum unique gene families

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
  118
  The following number genes are contained in these orthogorups:
  328
```

### 6.2.a) P. cactorum unique RxLR families

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
  RxLR_Pcac_uniq_groups=$RxLR_Dir/Pcac_uniq_RxLR_Orthogroups_hits.txt
  cat $RxLR_Orthogroup_10300 | grep -v 'Pinf' | grep -v 'Ppar' | grep -v 'Pcap' | grep -v 'Psoj' > $RxLR_Pcac_uniq_groups
  cat $RxLR_Pcac_uniq_groups | wc -l
  echo "These orthogroups contain the folloing number of RxLRs:"
  cat $RxLR_Pcac_uniq_groups | grep -w -o -f $RxLR_ID_10300 | wc -l
  echo "The following RxLRs were found in Group1 unique orthogroups:"
  RxLR_Group1_uniq_groups=$RxLR_Dir/Group1_RxLR_Orthogroups_hits.txt
  cat $RxLR_Orthogroup_10300 | grep -v 'Pcap' | grep -v 'Psoj' > $RxLR_Group1_uniq_groups
  cat $RxLR_Group1_uniq_groups | wc -l
  echo "These orthogroups contain the folloing number of RxLRs:"
  cat $RxLR_Group1_uniq_groups | grep -w -o -f $RxLR_ID_10300 | wc -l
```

```
  The number of RxLRs searched for is:
  145
  Of these, the following number were found in orthogroups:
  144
  These were distributed through the following number of Orthogroups:
  76
  The following RxLRs were found in Pcac unique orthogroups:
  2
  These orthogroups contain the folloing number of RxLRs:
  3
  The following RxLRs were found in Group1 unique orthogroups:
  15
  These orthogroups contain the folloing number of RxLRs:
  20

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
  1
```

#### 6.2.b) Extracting fasta files for orthogroups containing P. cactorum putative RxLRs

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogrouTxt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/Pcac_RxLR_Orthogroups.txt
  GoodProt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/goodProteins/goodProteins.fasta
  OutDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/orthogroups_fasta_Pcac_RxLR
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogrouTxt --fasta $GoodProt --out_dir $OutDir
```

#### 6.2.c) Extracting fasta files for Clade 1 orthogroups containing P. cactorum putative RxLRs
```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogrouTxt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/Group1_RxLR_Orthogroups_hits.txt
  GoodProt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/goodProteins/goodProteins.fasta
  OutDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_RxLR/orthogroups_fasta_clade1_RxLR
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogrouTxt --fasta $GoodProt --out_dir $OutDir
```

### 6.3.a) P. cactorum unique Crinkler families

P. cactorum strain 10300 crinkler genes were parsed to the same format as the gene
names used in the analysis:

```bash
  CRN_Names_10300=analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_Aug_CRN_hmmer_headers.txt
  WorkDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj
  CRN_Dir=$WorkDir/Pcac_CRN
  Orthogroups=$WorkDir/Pcac_Pinf_Ppar_Pcap_Psoj_orthogroups.txt
  CRN_ID_10300=$CRN_Dir/10300_CRN_hmmer_IDs.txt
  mkdir -p $CRN_Dir
  cat $CRN_Names_10300 | sed 's/g/Pcac|g/g' > $CRN_ID_10300
```

Ortholog groups containing CRN proteins were identified using the following
commands:
```bash
  echo "The number of CRNs searched for is:"
  cat $CRN_ID_10300 | wc -l
  echo "Of these, the following number were found in orthogroups:"
  CRN_Orthogroup_hits_10300=$CRN_Dir/Pcac_CRN_Orthogroups_hits.txt
  cat $Orthogroups | grep -o -w -f $CRN_ID_10300 > $CRN_Orthogroup_hits_10300
  cat $CRN_Orthogroup_hits_10300 | wc -l
  echo "These were distributed through the following number of Orthogroups:"
  CRN_Orthogroup_10300=$CRN_Dir/Pcac_CRN_Orthogroups.txt
  cat $Orthogroups | grep -w -f $CRN_ID_10300 > $CRN_Orthogroup_10300
  cat $CRN_Orthogroup_10300 | wc -l
  echo "The following CRNs were found in Pcac unique orthogroups:"
  CRN_Pcac_uniq_groups=$CRN_Dir/Pcac_uniq_CRN_Orthogroups_hits.txt
  cat $CRN_Orthogroup_10300 | grep -v 'Pinf' | grep -v 'Ppar' | grep -v 'Pcap' | grep -v 'Psoj' > $CRN_Pcac_uniq_groups
  cat $CRN_Pcac_uniq_groups | wc -l
  echo "The following CRNs were found in Group1 unique orthogroups:"
  cat $CRN_Pcac_uniq_groups | grep -w -o -f $CRN_ID_10300 | wc -l
  echo "The following CRNs were found in 10300 unique orthogroups:"
  CRN_Group1_uniq_groups=$CRN_Dir/Group1_CRN_Orthogroups_hits.txt
  cat $CRN_Orthogroup_10300 | grep -v 'Pcap' | grep -v 'Psoj' > $CRN_Group1_uniq_groups
  cat $CRN_Group1_uniq_groups | wc -l
  echo "The following CRNs were found in Group1 unique orthogroups:"
  cat $CRN_Group1_uniq_groups | grep -w -o -f $CRN_ID_10300 | wc -l
```

```
  The number of CRNs searched for is:
  92
  Of these, the following number were found in orthogroups:
  91
  These were distributed through the following number of Orthogroups:
  27
  The following CRNs were found in Pcac unique orthogroups:
  0
  The following CRNs were found in Group1 unique orthogroups:
  2
```

The P.cactorum CRN genes that were not found in orthogroups were identified:

```bash
  CRN_10300_uniq=$CRN_Dir/Pcac_unique_CRNs.txt
  cat $CRN_ID_10300 | grep -v -w -f $CRN_Orthogroup_hits_10300 | tr -d 'Pcac|' > $CRN_10300_uniq
  echo "The number of P.cac 10300 unique CRNs are:"
  cat $CRN_10300_uniq | wc -l
  CRN_Seq_10300=analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_aug_CRN_hmmer_out.fa
  Braker_genes=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
  CRN_10300_uniq_fa=$CRN_Dir/Pcac_unique_CRNs.fa
  cat $Braker_genes | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -w -A1 -f $CRN_10300_uniq | grep -E -v '^--' > $CRN_10300_uniq_fa
```

```
  The number of P.cac 10300 unique CRNs are:
  1
```


#### 6.3.b) Extracting fasta files for orthogroups containing P. cactorum putative CRNs

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogrouTxt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_CRN/Pcac_CRN_Orthogroups.txt
  GoodProt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/goodProteins/goodProteins.fasta
  OutDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_CRN/orthogroups_fasta_Pcac_CRN
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogrouTxt --fasta $GoodProt --out_dir $OutDir
```

#### 6.2.c) Extracting fasta files for Clade 1 orthogroups containing P. cactorum putative CRNs
```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogrouTxt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_CRN/Group1_CRN_Orthogroups_hits.txt
  GoodProt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/goodProteins/goodProteins.fasta
  OutDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj/Pcac_CRN/orthogroups_fasta_clade1_CRN
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogrouTxt --fasta $GoodProt --out_dir $OutDir
```

#### 6.3 Determining function of orthogroups

Lists of genes from P.cactorum unique genes, clade1 orthogorups and the largest
shared gene families were identified.


```bash
WorkDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj
InterProFile=gene_pred/interproscan/P.cactorum/10300/10300_interproscan.tsv
```

#### 6.3.a) Function of the 10 largest orthogroups

```bash
  FuncDir=$WorkDir/Ortholog_function/10_largest
  mkdir -p $FuncDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
  for num in $(seq 1 10); do
    cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | head -n$num | tail -n1 | sed -E 's/ /\n/g' | grep 'Pcac' | sed 's/Pcac|//g' > $FuncDir/orthogroup"$num"_gene_list.txt
    $ProgDir/feature_extractor.py --inp_txt $FuncDir/orthogroup"$num"_gene_list.txt --inp_interpro $InterProFile > $FuncDir/orthogroup"$num"_functions.txt
  done
```
#### 6.3.a) Function of core orthogroups

```bash
FuncDir=$WorkDir/Ortholog_function/common
mkdir -p $FuncDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
echo "The following number of Orthogroups were common to all species specific:"
groups_txt=$FuncDir/common_orthogroups.txt
cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Pcac'| grep 'Pinf' | grep 'Ppar' | grep 'Pcap' | grep 'Psoj' > $groups_txt
num_groups=$(cat $groups_txt | wc -l)
echo $num_groups
for num in $(seq 1 $num_groups); do
  Orthogroup=$(cat $groups_txt | head -n$num | tail -n1 | cut -f1 -d ':')
  echo "$Orthogroup" > $FuncDir/"$Orthogroup"_functions.txt
  echo "Extracting functions for: $Orthogroup"
  cat $Group1_uniq_groups| head -n$num | tail -n1 | sed -E 's/ /\n/g' | grep 'Pcac' | sed 's/Pcac|//g' > $FuncDir/"$Orthogroup"_gene_list.txt
  $ProgDir/feature_extractor.py --inp_txt $FuncDir/"$Orthogroup"_gene_list.txt --inp_interpro $InterProFile >> $FuncDir/"$Orthogroup"_functions.txt
done
```


#### 6.3.b) Function of the Clade1 unique orthogroups

```bash
  FuncDir=$WorkDir/Ortholog_function/clade1_unique
  mkdir -p $FuncDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
  echo "The following number of Orthogroups were Clade1 specific:"
  Group1_uniq_groups=$FuncDir/Group1_Orthogroups.txt
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | grep 'Pcac' | grep -v 'Pcap' | grep -v 'Psoj' > $Group1_uniq_groups
  cat $Group1_uniq_groups | wc -l
  for num in $(seq 1 404); do
    Orthogroup=$(cat $Group1_uniq_groups | head -n$num | tail -n1 | cut -f1 -d ':')
    echo "$Orthogroup" > $FuncDir/"$Orthogroup"_functions.txt
    echo "Extracting functions for: $Orthogroup"
    cat $Group1_uniq_groups| head -n$num | tail -n1 | sed -E 's/ /\n/g' | grep 'Pcac' | sed 's/Pcac|//g' > $FuncDir/"$Orthogroup"_gene_list.txt
    $ProgDir/feature_extractor.py --inp_txt $FuncDir/"$Orthogroup"_gene_list.txt --inp_interpro $InterProFile >> $FuncDir/"$Orthogroup"_functions.txt
  done
```
 The following commands were used to concatenate the individual results files
 into a single sumarry csv file. This can be imported into Excel.

 ```bash
   InString=""
  for num in $(ls $FuncDir/orthogroup*_functions.txt | sed "s&$FuncDir/orthogroup&&g" | sed 's/_functions.txt//g' | sort -n); do
    ThisFile=$(ls $FuncDir/orthogroup"$num"_functions.txt);
    InString="$InString $ThisFile";
  done
  paste -d ',' $InString > $FuncDir/Group1_functions.txt
 ```



# Run TribeMCL on orthoMCL data


## Set variable

```bash
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  IsolateAbrv=Pcac_Pinf_Ppar_Pcap_Psoj
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir/tribeMCL
```
## Run TribeMCL

```bash
  cat $WorkDir/"$IsolateAbrv"_blast.tab | cut -f1,2,11 | clusterx --method mcl -p inflation=1.5 - > $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/orthology/tribeMCl
  $ProgDir/tribeMCLgroups2tab.py --orthogroups $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL.txt --out_tab $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL_orthogroups.tab --out_txt $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL_orthogroups.txt
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
  $ProgDir/Pcac_Pinf_Ppar_Pcap_Psoj_venn.r --inp $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL_orthogroups.tab --out $WorkDir/tribeMCL/tribeMCL_orthogroups.pdf
```

```
  [1] "Pcac"
  [1] 0
  [1] 733
  [1] "Pcap"
  [1] 0
  [1] 357
  [1] "Psoj"
  [1] 0
  [1] 1243
  [1] "Pinf"
  [1] 0
  [1] 598
  [1] "Ppar"
  [1] 0
  [1] 648
```
