# For a comparison between RxLR compliments of 3 Phytophthora clade 1 isolates (P.cac 10300, P.inf T30-4, P.par 310), a clade 2 isolate (P.cap LT1534) & a clade 7 isolate (P.soj 67593)


```bash
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  IsolateAbrv=Pcac_Pinf_Ppar_Pcap_Psoj_RxLR
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
  Fasta_file=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_Total_RxLR_EER_motif_hmm_headers.fa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.inf T30-4
```bash
  Taxon_code=Pinf
  Unparsed_Fasta_file=analysis/RxLR_effectors/combined_evidence/P.infestans/T30-4/T30-4_Total_RxLR_EER_motif_hmm_headers.fa
  Fasta_file=$WorkDir/T30-4_Total_RxLR_EER_motif_hmm_headers_parsed.fa
  Id_field=1
  cat $Unparsed_Fasta_file | sed -e 's/supercont1..*_dna:supercontig_supercontig:ASM14294v1://g' | sed 's/supercont1./supercont1_/g' | tr ':' '-' > $Fasta_file
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.par 310
```bash
  Taxon_code=Ppar
  Unparsed_Fasta_file=analysis/RxLR_effectors/combined_evidence/P.parisitica/310/310_Total_RxLR_EER_motif_hmm_headers.fa
  Fasta_file=$WorkDir/310_Total_RxLR_EER_motif_hmm_headers_parsed.fa
  cat $Unparsed_Fasta_file | sed 's/Supercontig_2./Supercontig_2_/g' > $Fasta_file
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.cap LT1534
```bash
  Taxon_code=Pcap
  Unparsed_Fasta_file=analysis/RxLR_effectors/combined_evidence/P.capsici/LT1534/LT1534_Total_RxLR_EER_motif_hmm_headers.fa
  Fasta_file=$WorkDir/LT1534_Total_RxLR_EER_motif_hmm_headers_parsed.fa
  cat $Unparsed_Fasta_file | sed 's/^>jgi|Phyca11|/>/g' > $Fasta_file
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.soj 67593
```bash
  Taxon_code=Psoj
  Unparsed_Fasta_file=analysis/RxLR_effectors/combined_evidence/P.sojae/P6497/P6497_Total_RxLR_EER_motif_hmm_headers.fa
  Fasta_file=$WorkDir/P6497_Total_RxLR_EER_motif_hmm_headers_parsed.fa
  cat $Unparsed_Fasta_file | sed 's/^>jgi|Physo3|/>/g' > $Fasta_file
  Id_field=1
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
      sleep 5
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

<!-- ```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/venn_diagrams
  $ProgDir/ven_diag_5_way.R --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
``` -->

```bash
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
  $ProgDir/Pcac_Pinf_Ppar_Pcap_Psoj_venn.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs

At an inflation value of 1:
```
[1] "Pcac"
[1] 26
[1] 3
[1] "Pcap"
[1] 20
[1] 4
[1] "Psoj"
[1] 25
[1] 15
[1] "Pinf"
[1] 34
[1] 4
[1] "Ppar"
[1] 72
[1] 19
NULL
```

At an inflation value of 5:

```
[1] "Pcac"
[1] 26
[1] 5
[1] "Pcap"
[1] 20
[1] 15
[1] "Psoj"
[1] 25
[1] 35
[1] "Pinf"
[1] 34
[1] 8
[1] "Ppar"
[1] 72
[1] 33
NULL
```


# 6) Downstream analysis

## 6.1) Identification of Pseudogenes

Some predicted ORFs carried X's in their amino acid sequence. These ORFs are not
likely to be real genes. The number of Pseudogenes were identified using the
following commands:

```bash
  cat $GoodProts | grep -v '>' | grep 'X' | wc -l
  cat $GoodProts | grep '>' |  wc -l
```

This identified that 32 of 1701 RxLRs in the analysis carried X's in their amino
acid sequence


### 6.2.a ) P. cactotum unique gene families

The genes unique to P.cactorum were identified within the orthology analysis.

First variables were set:
```bash
  IsolateAbrv=Pcac_Pinf_Ppar_Pcap_Psoj_RxLR
  WorkDir=analysis/orthology/orthomcl/"$IsolateAbrv"
  PcacUniqDir=$WorkDir/Pcac_unique
  Orthogroups=$WorkDir/"$IsolateAbrv"_orthogroups.txt
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
  echo "The following number of genes are contained in these orthogorups:"
  cat $Uniq_Pcac_groups | grep -o 'Pcac|' | wc -l  
```

```
  The number of orthogroups unique to P. cactorum are:
  3
  The following number of genes are contained in these orthogorups:
  6
```

#### 6.2.b)

```bash


$WorkDir/Pcac_Pinf_Ppar_Pcap_Psoj_RxLR_orthogroups.txt | head -n1 | sed -E 's/ /\n/g' | grep 'Pcac' | sed 's/Pcac|//g' | wc-l
$WorkDir/Pcac_Pinf_Ppar_Pcap_Psoj_RxLR_orthogroups.txt | head -n1 | grep -o -E 'P\w*\|' | sort | uniq -c
RxLRDir=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300
$WorkDir/Pcac_Pinf_Ppar_Pcap_Psoj_RxLR_orthogroups.txt | head -n1 | grep -o -E 'Pcac\|\w*' | sed 's/Pcac|//g' | grep -f $RxLRDir/10300_Total_RxLR_EER_WY_motif_hmm_headers.txt | wc -l

```

The largest RxLR group contained 141 proteins, 19 of which were from P.cactorum.
Of these 19, 15 contained WY domains.

'''
  19 Pcac|
  16 Pcap|
  29 Pinf|
  38 Ppar|
  39 Psoj|
'''

#### 6.3) Extracting fasta files for all orthogroups

```bash
  IsolateAbrv=Pcac_Pinf_Ppar_Pcap_Psoj_RxLR
  WorkDir=analysis/orthology/orthomcl/"$IsolateAbrv"
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogroupTxt=$WorkDir/orthogroups_fasta/"$IsolateAbrv"_Orthogroups.txt
  GoodProt=$WorkDir/goodProteins/goodProteins.fasta
  OutDir=$WorkDir/orthogroups_fasta_inflation_5
  mkdir -p $OutDir
  # cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | cut -f1 -d ' ' > $OrthogroupTxt
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $WorkDir/"$IsolateAbrv"_orthogroups.txt --fasta $GoodProt --out_dir $OutDir
```

# Run TribeMCL on orthoMCL data


## Run TribeMCL

```bash
  mkdir -p $WorkDir/tribeMCL
  cat $WorkDir/"$IsolateAbrv"_blast.tab | cut -f1,2,11 | clusterx --method mcl -p inflation=5 - > $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/orthology/tribeMCl
  $ProgDir/tribeMCLgroups2tab.py --orthogroups $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL.txt --out_tab $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL_orthogroups.tab --out_txt $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL_orthogroups.txt
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
  $ProgDir/Pcac_Pinf_Ppar_Pcap_Psoj_venn.r --inp $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL_orthogroups.tab --out $WorkDir/tribeMCL/tribeMCL_orthogroups.pdf
```

```
  [1] "Pcac"
  [1] 0
  [1] 51
  [1] "Pcap"
  [1] 0
  [1] 35
  [1] "Psoj"
  [1] 0
  [1] 91
  [1] "Pinf"
  [1] 0
  [1] 34
  [1] "Ppar"
  [1] 0
  [1] 31
```
