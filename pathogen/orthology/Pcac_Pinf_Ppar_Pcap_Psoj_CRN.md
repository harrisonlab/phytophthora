# For a comparison between CRN compliments of 3 Phytophthora clade 1 isolates (P.cac 10300, P.inf T30-4, P.par 310), a clade 2 isolate (P.cap LT1534) & a clade 7 isolate (P.soj 67593)


```bash
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  IsolateAbrv=Pcac_Pinf_Ppar_Pcap_Psoj_CRN_publication
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins
```

## 1) Format fasta files

### for P.cac 10300



```bash
OutDir=$WorkDir
mkdir -p $OutDir
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
cat $AnnotTab | grep -w 'CRN' | cut -f1,26 | sed "s/^/>/g" | sed "s/\t/\n/g" | sed "s/X$//g" > $OutDir/Pcac_CRN.fa
```

```bash
  Taxon_code=Pcac
  Fasta_file=$OutDir/Pcac_CRN.fa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.inf T30-4
```bash
  Taxon_code=Pinf
  Unparsed_Fasta_file=analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_Total_CRN.fa
  Fasta_file=$WorkDir/T30-4_Total_CRN_parsed.fa
  cat $Unparsed_Fasta_file | sed -e 's/supercont1..*_dna:supercontig_supercontig:ASM14294v1://g' | sed 's/supercont1./supercont1_/g' | tr ':' '-' > $Fasta_file
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.par 310
```bash
  Taxon_code=Ppar
  Unparsed_Fasta_file=analysis/CRN_effectors/hmmer_CRN/P.parisitica/310/310_Total_CRN.fa
  Fasta_file=$WorkDir/310_Total_CRN_parsed.fa
  cat $Unparsed_Fasta_file | sed 's/Supercontig_2./Supercontig_2_/g' > $Fasta_file
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.cap LT1534
```bash
  Taxon_code=Pcap
  Unparsed_Fasta_file=analysis/CRN_effectors/hmmer_CRN/P.capsici/LT1534/LT1534_Total_CRN.fa
  Fasta_file=$WorkDir/LT1534_Total_CRN_parsed.fa
  cat $Unparsed_Fasta_file | sed 's/^>jgi|Phyca11|/>/g' > $Fasta_file
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P.soj P6497
```bash
  Taxon_code=Psoj
  Unparsed_Fasta_file=analysis/CRN_effectors/hmmer_CRN/P.sojae/P6497/P6497_Total_CRN.fa
  Fasta_file=$WorkDir/P6497_Total_CRN_parsed.fa
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
  MergeHits=$WorkDir/"$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```

Also try using orthofinder

```bash
qlogin -pe smp 16

ProjDir=/home/groups/harrisonlab/project_files/idris
cd $ProjDir
IsolateAbrv=Pcac_Pinf_Ppar_Pcap_Psoj_CRN_publication
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv

#16 threads used
orthofinder -f $WorkDir/formatted -t 16
```

```
Orthogroup statistics:
   Statistics_PerSpecies.csv   Statistics_Overall.csv   Orthogroups_SpeciesOverlaps.csv

OrthoFinder assigned 547 genes (84.2% of total) to 42 orthogroups. Fifty percent of all genes were in orthogroups
with 19 or more genes (G50 was 19) and were contained in the largest 11 orthogroups (O50 was 11). There were 8
orthogroups with all species present and 0 of these consisted entirely of single-copy genes.

When publishing work that uses OrthoFinder please cite:
    D.M. Emms & S. Kelly (2015), OrthoFinder: solving fundamental biases in whole genome comparisons
    dramatically improves orthogroup inference accuracy, Genome Biology 16:157.
```

## 5) Plot venn diagrams:

Orthofinder output:

```bash
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  $ProgDir/orthoMCLgroups2tab.py $GoodProts $WorkDir/formatted/Results_Jan15/Orthogroups.txt > $WorkDir/formatted/Results_Jan15/"$IsolateAbrv"_orthogroups.tab
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
  $ProgDir/Pcac_Pinf_Ppar_Pcap_Psoj_venn.r --inp $WorkDir/formatted/Results_Jan15/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/formatted/Results_Jan15/"$IsolateAbrv"_orthogroups.pdf
```


Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs

```
  [1] "Pcac"
  [1] 41
  [1] 0
  [1] 0
  [1] "Pcap"
  [1] 49
  [1] 0
  [1] 0
  [1] "Psoj"
  [1] 83
  [1] 0
  [1] 0
  [1] "Pinf"
  [1] 52
  [1] 0
  [1] 0
  [1] "Ppar"
  [1] 19
  [1] 0
  [1] 0
```


## Extracting fasta files for orthogroups:

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogrouTxt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN_publication/formatted/Results_Jan15/Orthogroups.txt
  GoodProt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN_publication/goodProteins/goodProteins.fasta
  OutDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN_publication/formatted/Results_Jan15/orthogroups_fasta_Pcac
  #
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogrouTxt --fasta $GoodProt --out_dir $OutDir
```

```bash
ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
OrthogrouTxt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN_publication/*orthogroups.txt
GoodProt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN_publication/goodProteins/goodProteins.fasta
OutDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN_publication/orthogroups_fasta_Pcac
mkdir -p $OutDir
$ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogrouTxt --fasta $GoodProt --out_dir $OutDir
```

# 6) Downstream analysis


## 6.0) Identifying Signal peptides in Crinklers



Signal peptides were identified in the predicted Crinklers:

```bash
  for Crinklers in $(ls analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/formatted/*.fasta); do
    echo "$Crinklers"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    Strain=$(echo $Crinklers | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Crinklers | rev | cut -f3 -d '/' | rev)
    qsub $ProgDir/pred_sigP.sh $Crinklers
  done
```



```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogroupTxt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/Pcac_Pinf_Ppar_Pcap_Psoj_CRN_orthogroups.txt
  GoodProt=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/goodProteins/goodProteins.fasta
  OutDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/fasta/all_orthogroups
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogroupTxt --fasta $GoodProt --out_dir $OutDir > $OutDir/extractionlog.txt
  for File in $(ls $OutDir/orthogroup*.fa); do
    cat $File | grep '>' | tr -d '> '
  done > $OutDir/orthogroup_genes.txt
  cat $GoodProt | grep '>' | tr -d '> ' | grep -v -f $OutDir/orthogroup_genes.txt > $OutDir/singleton_genes.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $GoodProt --headers $OutDir/singleton_genes.txt > $OutDir/singleton_genes.fa
  echo "The numbe of singleton genes extracted is:"
  cat $OutDir/singleton_genes.fa | grep '>' | wc -l
  cat $OutDir/singleton_genes.fa | grep -A1 '>Pcac' > $OutDir/Pcac_singletons.fa
  echo "the number of singleton genes in P.cactorum is:"
  cat $OutDir/Pcac_singletons.fa | grep '>' | wc -l
```

```bash
  cat gene_pred/orthomcl_sigP/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/formatted/split/P*.fasta_sp.aa | grep '>' | cut -f1 | sed 's/>//g' | tr -d ' ' > tmp.txt
  # cat analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/Pcac_Pinf_Ppar_Pcap_Psoj_CRN_orthogroups.txt | grep -f tmp.txt | less -S
  for File in $(ls -v analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/fasta/all_orthogroups/*.fa); do
  echo $(basename $File);
  echo "total genes:"
  cat $File | grep '>' | tr -d '>' | cut -f1 -d '|' | sort | uniq -c | sort -r -n
  echo "SigP genes:"
  cat $File | grep -f tmp.txt | tr -d '>' | cut -f1 -d '|' | sort | uniq -c | sort -r -n
  cat $File | grep -f tmp.txt;
  echo ""
  done > gene_pred/orthomcl_sigP/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/formatted/split/CRNs_containing_sigP.txt
```

## 6.1) Identification of Pseudogenes

Some predicted ORFs carried X's in their amino acid sequence. These ORFs are not
likely to be real genes. The number of Pseudogenes were identified using the
following commands:

```bash
  cat $GoodProts | grep -v '>' | grep 'X' | wc -l
  cat $GoodProts | grep '>' |  wc -l
```
<!--
This identified that 7 of 843 CRNs in the analysis carried X's in their amino
acid sequence. These genes were:
```
Pcac|g12904.t1 Pcac|g13271.t1 Pcac|g15212.t1
Pcac|g15627.t1 Pcac|g17041.t1 Pcac|g17676.t1 504842|fgenesh2_kg.PHYCAscaffold_10_#_2_#_4097787:1
``` -->




### 6.2 ) P. cactotum unique gene families

The genes unique to P.cactorum were identified within the orthology analysis.

First variables were set:
```bash
  IsolateAbrv=Pcac_Pinf_Ppar_Pcap_Psoj_CRN
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

#### 6.2) Extracting fasta files for all orthogroups

```bash
  IsolateAbrv=Pcac_Pinf_Ppar_Pcap_Psoj_CRN
  WorkDir=analysis/orthology/orthomcl/"$IsolateAbrv"
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  # OrthogroupTxt=$WorkDir/orthogroups_fasta/"$IsolateAbrv"_Orthogroups.txt
  GoodProt=$WorkDir/goodProteins/goodProteins.fasta
  OutDir=$WorkDir/orthogroups_fasta
  mkdir -p $OutDir
  # cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | cut -f1 -d ' ' > $OrthogroupTxt
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $WorkDir/"$IsolateAbrv"_orthogroups.txt --fasta $GoodProt --out_dir $OutDir
  GenesInOrthogroups=$OutDir/"$IsolateAbrv"_orthogroup_genes.txt
  GenesNotInOrthogroupsTxt=$OutDir/"$IsolateAbrv"_non-orthogroup_genes.txt
  GenesNotInOrthogroupsFasta=$OutDir/"$IsolateAbrv"_non-orthogroup_genes.fa
  cat $WorkDir/"$IsolateAbrv"_orthogroups.txt | sed 's/ /\n/g' | grep -v 'orthogroup' > $GenesInOrthogroups
  cat $WorkDir/goodProteins/goodProteins.fasta $GenesInOrthogroups | grep '|' | tr -d '>' | cut -f1 -d '.' | sort | uniq -u > $GenesNotInOrthogroupsTxt
  cat $WorkDir/goodProteins/goodProteins.fasta | grep -w -A1 -f $GenesNotInOrthogroupsTxt | grep -v -E '^--'> $GenesNotInOrthogroupsFasta
```

The names of Crinklers in each orthogroup (& P.cactorum orthogroups) were
extracted using the commands:

```bash
  for File in $(ls -v analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/fasta/all_orthogroups/*.fa | grep -v 'Pcac_singletons'); do
    Orthogroup=$(basename $File | sed 's/.fa//g');
    for Gene in $(cat $File | grep '>' | sed 's/>//g'); do
      printf "$Orthogroup\t$Gene\n";
    done;
  done > analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/fasta/all_orthogroups/CRN_headers_by_orthogroup.txt
  for File in $(ls -v analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/fasta/all_orthogroups/*.fa | grep -v 'Pcac_singletons'); do
    Orthogroup=$(basename $File | sed 's/.fa//g');
    for Gene in $(cat $File | grep 'Pcac' | sed 's/>Pcac|//g'); do
      printf "$Orthogroup\t$Gene\n";
    done;
  done > analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_CRN/fasta/all_orthogroups/Pcac_CRN_headers_by_orthogroup.txt
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
  [1] 1
  [1] "Pcap"
  [1] 0
  [1] 1
  [1] "Psoj"
  [1] 0
  [1] 4
  [1] "Pinf"
  [1] 0
  [1] 2
  [1] "Ppar"
  [1] 0
  [1] 0
```


## Identify the number of incomplete crinklers


```bash
OutDir=analysis/gene_tables/P.cactorum/10300/subset
cat $OutDir/10300_gene_table_CRN.tsv | cut -f26 | grep "X$" | wc -l
cat $OutDir/10300_gene_table_CRN.tsv | cut -f26 | grep -v "X$" | wc -l
```
