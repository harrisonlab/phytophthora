
# Run TribeMCL on Pcactorum Crinklers


## Set variable

```bash
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  IsolateAbrv=Pcac_CRN_publication
  WorkDir=analysis/orthology/tribemcl/$IsolateAbrv
  mkdir -p $WorkDir
```



## 1) Format fasta files

## Investigating Crinklers:

```bash
OutDir=$WorkDir/formatted
mkdir -p $OutDir
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
cat $AnnotTab | grep -w 'CRN' | cut -f1,26 | sed "s/^/>/g" | sed "s/\t/\n/g" | sed "s/X$//g" > $OutDir/Pcac.fa
```

```bash
  Taxon_code=Pcac
  Fasta_file=$OutDir/Pcac.fa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


## 2) Filter proteins into good and poor sets.
<!--
```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
``` -->


## 3.1) Perform an all-vs-all blast of the proteins

```bash
  Good_proteins_file=$WorkDir/formatted/Pcac.fasta
  BlastDB=$WorkDir/blastall/$IsolateAbrv.db

  makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
  BlastOut=$WorkDir/all-vs-all_results.tsv
  mkdir -p $WorkDir/splitfiles

  SplitDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  $SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $WorkDir/splitfiles --out_base goodProteins

  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology  
  for File in $(ls $WorkDir/splitfiles/goodProteins_*.fa); do
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

## Run TribeMCL

```bash
  mkdir -p $WorkDir/tribeMCL
  cat $WorkDir/"$IsolateAbrv"_blast.tab | cut -f1,2,11 | clusterx --method mcl -p inflation=1.5 - > $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/orthology/tribeMCl
  $ProgDir/tribeMCLgroups2tab.py --orthogroups $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL.txt --out_tab $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL_orthogroups.tab --out_txt $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL_orthogroups.txt
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
  $ProgDir/Pcac_Pinf_Ppar_Pcap_Psoj_venn.r --inp $WorkDir/tribeMCL/"$IsolateAbrv"_tribeMCL_orthogroups.tab --out $WorkDir/tribeMCL/tribeMCL_orthogroups.pdf
```
