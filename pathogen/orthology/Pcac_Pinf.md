# Orthology analysis between P. cactorum and P. idaei
 Including 13 P. cactorum crown rot, 2 P. cactorum leather rot isolates, 3 P. cactorum apple isolates and three P. idaei raspberry isolates


```bash
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  IsolateAbrv=Pcac_Pinf_publication
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins
```

## 1.a) Format fasta files - P. cactorum crown rot


### for P. cactorum crown rot isolate P414

```bash
  Taxon_code=Pc_CR1
  Fasta_file=$(ls /data/scratch/armita/idris/gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P. cactorum crown rot isolate 12-420

```bash
  Taxon_code=Pc_CR2
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/12420/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P. cactorum crown rot isolate 15-13

```bash
  Taxon_code=Pc_CR3
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/15_13/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


### for P. cactorum crown rot isolate 15-7

```bash
  Taxon_code=Pc_CR4
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/15_7/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


### for P. cactorum crown rot isolate 2003-3

```bash
  Taxon_code=Pc_CR5
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/2003_3/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


### for P. cactorum crown rot isolate 4032

```bash
  Taxon_code=Pc_CR6
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/4032/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


### for P. cactorum crown rot isolate 4040

```bash
  Taxon_code=Pc_CR7
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/4040/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


### for P. cactorum crown rot isolate P404

```bash
  Taxon_code=Pc_CR8
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/404/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


### for P. cactorum crown rot isolate P415

```bash
  Taxon_code=Pc_CR9
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/415/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


### for P. cactorum crown rot isolate P416

```bash
  Taxon_code=Pc_CR10
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/416/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


### for P. cactorum crown rot isolate P421

```bash
  Taxon_code=Pc_CR11
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/P421_v2/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```


### for P. cactorum crown rot isolate PC13/15

```bash
  Taxon_code=Pc_CR12
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/PC13_15/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P. cactorum crown rot isolate 10300

```bash
  Taxon_code=Pc_CR13
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```



## 1.b) Format fasta files - P. cactorum leather rot

### for P. cactorum leather rot isolate 11-40

```bash
  Taxon_code=Pc_LR1
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/11-40/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P. cactorum leather rot isolate 17-21

```bash
  Taxon_code=Pc_LR2
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/17-21/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```



## 1.c) Format fasta files - P. cactorum malus x domestica

### for  P. cactorum malus x domestica isolate 62471

```bash
  Taxon_code=Pc_MD1
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/62471/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for  P. cactorum malus x domestica isolate P295

```bash
  Taxon_code=Pc_MD2
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/P295/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for  P. cactorum malus x domestica isolate R36/14

```bash
  Taxon_code=Pc_MD3
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/R36_14/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```



## 1.d) Format fasta files - P. idaei Rubus idaeu

### for P. idaei Rubus idaeu isolate P371

```bash
  Taxon_code=Pi_RI1
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.idaei/371/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P. idaei Rubus idaeu isolate SCRP370

```bash
  Taxon_code=Pi_RI2
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.idaei/SCRP370/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for P. idaei Rubus idaeu isolate SCRP376

```bash
  Taxon_code=Pi_RI3
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.idaei/SCRP376/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
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
  MergeHits="$IsolateAbrv"_blast.tab
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

Orthomcl results


```bash
Orthogroups=$(ls analysis/orthology/orthomcl/Pcac_Pinf_publication/Pcac_Pinf_publication_orthogroups.txt)
echo "Number of orthgroups"
cat $Orthogroups | wc -l
echo "Number in all isoaltes:"
echo "Orthogroups present in all isoaltes"
cat $Orthogroups | grep -w -e "Pc_CR1" | grep -w -e "Pc_CR2" | grep -w -e "Pc_CR3" | grep -w -e "Pc_CR4" | grep -w -e "Pc_CR5" | grep -w -e "Pc_CR6" | grep -w -e "Pc_CR7" | grep -w -e "Pc_CR8" | grep -w -e "Pc_CR9" | grep -w -e "Pc_CR10" | grep -w -e "Pc_CR11" | grep -w -e "Pc_CR12" | grep -w -e "Pc_CR13" | grep -w -e "Pc_LR1" | grep -w -e "Pc_LR2" | grep -w -e "Pc_MD1" | grep -w -e "Pc_MD2" | grep -w -e "Pc_MD3" | grep -w -e "Pi_RI1" | grep -w -e "Pi_RI2" | grep -w -e "Pi_RI3" | wc -l
echo "orthogroups found in all cactorum ex strawberry but not in ex apple"
cat $Orthogroups | grep -w -e "Pc_CR1" | grep -w -e "Pc_CR2" | grep -w -e "Pc_CR3" | grep -w -e "Pc_CR4" | grep -w -e "Pc_CR5" | grep -w -e "Pc_CR6" | grep -w -e "Pc_CR7" | grep -w -e "Pc_CR8" | grep -w -e "Pc_CR9" | grep -w -e "Pc_CR10" | grep -w -e "Pc_CR11" | grep -w -e "Pc_CR12" | grep -w -e "Pc_CR13" | grep -v -w -e "Pc_MD1" -e "Pc_MD2" -e "Pc_MD3" | wc -l
echo "orthogroups found in all cactorum ex apple but not in ex strawberry"
cat $Orthogroups | grep -w -e "Pc_CR1" | grep -w -e "Pc_CR2" | grep -w -e "Pc_CR3" | grep -w -e "Pc_CR4" | grep -w -e "Pc_CR5" | grep -w -e "Pc_CR6" | grep -w -e "Pc_CR7" | grep -w -e "Pc_CR8" | grep -w -e "Pc_CR9" | grep -w -e "Pc_CR10" | grep -w -e "Pc_CR11" | grep -w -e "Pc_CR12" | grep -w -e "Pc_CR13" | grep -v -w -e "Pc_MD1" -e "Pc_MD2" -e "Pc_MD3" | wc -l

cat $Orthogroups | grep -e "Pc_MD1" | grep  -e "Pc_MD2" | grep -e "Pc_MD3" | grep -v -w -e "Pc_CR1" -e "Pc_CR2" -e "Pc_CR3" -e "Pc_CR4" -e "Pc_CR5" -e "Pc_CR6" -e "Pc_CR7" -e "Pc_CR8" -e "Pc_CR9" -e "Pc_CR10" -e "Pc_CR11" -e "Pc_CR12" -e "Pc_CR13" | wc -l


cat $Orthogroups | grep -w -e "Pc_CR1" -e "Pc_CR2" -e "Pc_CR3" -e "Pc_CR4" -e "Pc_CR5" -e "Pc_CR6" -e "Pc_CR7" -e "Pc_CR8" -e "Pc_CR9" -e "Pc_CR10" -e "Pc_CR11" -e "Pc_CR12" -e "Pc_CR13" -e "Pc_LR1" -e "Pc_LR2" -e "Pc_MD1" -e "Pc_MD2" -e "Pc_MD3" -e "Pi_RI1" -e "Pi_RI2" -e "Pi_RI3" | wc -l
```


Also try using orthofinder

```bash
qlogin -pe smp 24 -l virtual_free=1G -l h=blacklace07.blacklace

#16 threads used
ProjDir=/home/groups/harrisonlab/project_files/idris
cd $ProjDir
IsolateAbrv=Pcac_Pinf_publication
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
orthofinder -f $WorkDir/formatted -t 24 -a 24
```

orthofinder results:

```
OrthoFinder assigned 157038 genes (99.2% of total) to 14187 orthogroups. Fifty percent of all genes were in orthogroups
with 12 or more genes (G50 was 12) and were contained in the largest 6160 orthogroups (O50 was 6160). There were 10669
orthogroups with all species present and 10016 of these consisted entirely of single-copy genes.
```

output files are in:
```bash
ls $WorkDir/formatted/Results_Apr10
```

```bash
OutDir=$(ls -d analysis/orthology/orthomcl/Pcac_Pinf_publication)
Orthogroups=$(ls $OutDir/Pcac_Pinf_publication_orthogroups.txt)
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
$ProgDir/summarise_orthogroups.py --orthogroups $Orthogroups > $OutDir/summarised_orthogroups.tsv
```

```bash
OutDir=$(ls -d analysis/orthology/orthomcl/Pcac_Pinf_publication)
Orthogroups=$(ls $OutDir/Pcac_Pinf_publication_orthogroups.txt)
OrthoMclIds='Pc_CR1 Pc_CR2 Pc_CR3 Pc_CR4 Pc_CR5 Pc_CR6 Pc_CR7 Pc_CR8 Pc_CR9 Pc_CR10 Pc_CR11 Pc_CR12 Pc_CR13 Pc_LR1 Pc_LR2 Pc_MD1 Pc_MD2 Pc_MD3 Pi_RI1 Pi_RI2 Pi_RI3'
TrueNames='414 12420 15_13 15_7 2003_3 4032 4040 404 415 416 P421 PC13_15 10300 11-40 17-21 62471 P295 R36_14 371 SCRP370 SCRP376'
MissingIsolates='LV007'
mkdir -p $OutDir/gloome
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
$ProgDir/orthoMCL2Gloome.py --orthogroups $Orthogroups --isolate_list $OrthoMclIds --true_names $TrueNames --missing_isolates $MissingIsolates > $OutDir/gloome/gloome.fasta
```

Gloome was used to predict gain/loss of orthogroups within the tree
http://gloome.tau.ac.il/source/GLOOME.CoPAP.gainLoss.Manual.pdf


A tree file from my local computer was used as a reference
(Pcac_phylogeny.consensus.scored_geneious2.tre). The tree file contents needed to be edited to remove node labels.

The original tree was:

```
((SCRP370:1.0,(371:1.0,SCRP376:1.0)0.9:0.08245622453296164)1:1.6434064099824397,(LV007:1.0000000000000004,((17-21:1.0,(R36_14:1.0,(P295:1.0,62471:1.0)1:0.29696569109508353)1:0.32373562653903853)1:0.15949253028150512,(10300:1.0,(((2003_3:1.0,4032:1.0)0.54:0.03507127382577924,(P421:1.0,(414:1.0,PC13_15:1.0)0.48:0.027700345429534146)0.48:0.019840795963525615)0.61:0.03554532977083369,((4040:1.0,416:1.0)0.77:0.05715039539753075,(11-40:1.0,((12420:1.0,15_13:1.0)0.69:0.04700123987863236,(404:1.0,(15_7:1.0,415:1.0)0.45:0.015690963218303544)0.47:0.021387852567743337)0.52:0.02567263702771161)0.37:0.0046454517207186186)0.39:0.007662872745568983)0.88:0.07713662906987118)1:0.5668858014697449)1:1.0567984198056735):1.6434064099824397);
```
and it was edited to:
```
((SCRP370:1.0,(371:1.0,SCRP376:1.0):1.0):1.0,(LV007:1.0,((17-21:1.0,(R36_14:1.0,(P295:1.0,62471:1.0):1.0):1.0):1.0,(10300:1.0,(((2003_3:1.0,4032:1.0):1.0,(P421:1.0,(414:1.0,PC13_15:1.0):1.0):1.0):1.0,((4040:1.0,416:1.0):1.0,(11-40:1.0,((12420:1.0,15_13:1.0):1.0,(404:1.0,(15_7:1.0,415:1.0):1.0):1.0):1.0):1.0):1.0):1.0):1.0):1.0):1.0);
```

```bash
OutDir=analysis/orthology/orthomcl/Pcac_Pinf_publication/gloome
mkdir -p $OutDir
ParamFile=$OutDir/param.txt
Tree=tmp.nwk
printf "((SCRP370:1.0,(371:1.0,SCRP376:1.0):1.0):1.0,(LV007:1.0,((17-21:1.0,(R36_14:1.0,(P295:1.0,62471:1.0):1.0):1.0):1.0,(10300:1.0,(((2003_3:1.0,4032:1.0):1.0,(P421:1.0,(414:1.0,PC13_15:1.0):1.0):1.0):1.0,((4040:1.0,416:1.0):1.0,(11-40:1.0,((12420:1.0,15_13:1.0):1.0,(404:1.0,(15_7:1.0,415:1.0):1.0):1.0):1.0):1.0):1.0):1.0):1.0):1.0):1.0);" > $Tree
Fasta=$(ls $OutDir/gloome.fasta)
# printf \
# "_seqFile $Fasta
# _treeFile $Tree
# _gainLossDist 1
# _optimizationLevel low
# _gainDistributionType GENERAL_GAMMA_PLUS_INV
# _lossDistributionType GENERAL_GAMMA_PLUS_INV
# " \
# > $ParamFile
printf \
"_seqFile $Fasta
_treeFile $Tree
_gainLossDist 1
_optimizationLevel high
_minNumOfOnes 2
_minNumOfZeros 1
_gainDistributionType GENERAL_GAMMA_PLUS_INV
_lossDistributionType GENERAL_GAMMA_PLUS_INV
" \
> $ParamFile
gainLoss.VR01.266 $ParamFile
mv RESULTS $OutDir/RESULTS3
ls $OutDir/RESULTS2/TheTree.INodes.ph
# Analysis of the output tree identified the following nodes:
# N2 P. idaei clade
# N4 LV007 clade (no results)
# N5 P. cactorum clade
# N6 Leather rot + apple clade
# N7 apple clade
# N9 crown rot clade
# maximum parsimony
OutFile=$(ls $OutDir/RESULTS2/gainLossMP.2.00099.PerPosPerBranch.txt)
# Stochastic modeling
# OutFile=$(ls $OutDir/RESULTS2/ProbabilityPerPosPerBranch.txt)
cat $OutFile | grep -w 'N2' > $OutDir/N2_P.idaei_clade.txt
cat $OutFile | grep -w 'N4' > $OutDir/N4_LV007_clade.txt
cat $OutFile | grep -w 'N5' > $OutDir/N5_P.cactorum_clade.txt
cat $OutFile | grep -w 'N6' > $OutDir/N6_leather-rot_apple_clade.txt
cat $OutFile | grep -w 'N7' > $OutDir/N7_apple_clade.txt
cat $OutFile | grep -w 'N9' > $OutDir/N9_crown_clade.txt
for File in $(ls $OutDir/*_clade.txt); do
  Name=$(basename $File)
  Gains=$(cat $File | grep 'gain' | wc -l)
  Losses=$(cat $File | grep 'loss' | wc -l)
  Total=$(cat $File | grep -e 'gain' -e 'loss' | wc -l)
  printf "$Name\t$Gains\t$Losses\t$Total\n"
done
```

Maximum parsimony results:
```
N2_P.idaei_clade.txt	0	6409	6409
N4_LV007_clade.txt	0	0	0
N5_P.cactorum_clade.txt	0	3098	3098
N6_leather-rot_apple_clade.txt	0	696	696
N7_apple_clade.txt	0	507	507
N9_crown_clade.txt	0	1065	1065
```

stochastic modeling results:

```
N2_P.idaei_clade.txt	4652	8068	12720
N4_LV007_clade.txt	8394	4481	12875
N5_P.cactorum_clade.txt	5014	8387	13401
N6_leather-rot_apple_clade.txt	5243	7639	12882
N7_apple_clade.txt	4752	5236	9988
N9_crown_clade.txt	4357	8473	12830
```


```bash
Orthogroups=$(ls analysis/orthology/orthomcl/Pcac_Pinf_publication/Pcac_Pinf_publication_orthogroups.txt)
OrthoMclIds='Pc_CR1 Pc_CR2 Pc_CR3 Pc_CR4 Pc_CR5 Pc_CR6 Pc_CR7 Pc_CR8 Pc_CR9 Pc_CR10 Pc_CR11 Pc_CR12 Pc_CR13 Pc_LR1 Pc_LR2 Pc_MD1 Pc_MD2 Pc_MD3 Pi_RI1 Pi_RI2 Pi_RI3'
TrueNames='414 12420 15_13 15_7 2003_3 4032 4040 404 415 416 P421 PC13_15 10300 11-40 17-21 62471 P295 R36_14 371 SCRP370 SCRP376'
OutDir=analysis/orthology/orthomcl/Pcac_Pinf_publication/gloome
for CladeFile in $(ls $OutDir/*_clade.txt | grep 'N2'); do
  Prefix=$(basename $CladeFile .txt)
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  $ProgDir/gloome2orthogroups.py --orthogroups $Orthogroups --gloome $OutDir/N2_P.idaei_clade.txt --OrthoMCL_all $OrthoMclIds  > $OutDir/${Prefix}_orthogroup_summary.txt
  cat $OutDir/${Prefix}_orthogroup_summary.txt \
  | grep -e "Pc_CR(13)" -e "Pc_CR(0)" \
  | grep -e "Pc_MD(3)" -e "Pc_MD(0)" \
  | grep -e "Pci_RI(3)" -e "Pi_RI(0)" \
  > $OutDir/${Prefix}_orthogroup_summary_filtered.txt
done
```


<!--
## 5) Plot venn diagrams:

Orthofinder output:

```bash
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  $ProgDir/orthoMCLgroups2tab.py $GoodProts $WorkDir/formatted/Results_Apr10/Orthogroups.txt > $WorkDir/formatted/Results_Apr10/"$IsolateAbrv"_orthogroups.tab
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
  $ProgDir/Pcac_Pinf_Ppar_Pcap_Psoj_venn.r --inp $WorkDir/formatted/Results_Apr10/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/formatted/Results_Apr10/"$IsolateAbrv"_orthogroups.pdf
```

Orthomcl output:
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
``` -->

# 6) Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.


```bash
IsolateAbrv=Pcac_Pinf_publication
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
OrthogroupsTxt=$(ls $WorkDir/${IsolateAbrv}_orthogroups.txt)
echo "Number of orthogroups"
cat $OrthogroupsTxt | wc -l
echo "Number of clustered proteins:"
cat $OrthogroupsTxt | grep -o '|' | wc -l
echo "Total proteins:"
GoodProts=$(ls $WorkDir/goodProteins/goodProteins.fasta)
cat $GoodProts | grep '>' | wc -l
```

## Extracting protein fasta files for orthogroups:

```bash
  IsolateAbrv=Pcac_Pinf_publication
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogroupsTxt=$(ls $WorkDir/${IsolateAbrv}_orthogroups.txt)
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  OutDir=$WorkDir/orthogroups_fasta
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogroupsTxt --fasta $GoodProts --out_dir $OutDir
```

## Extracting nucleotide fasta files for orthogroups:

```bash
  mkdir $WorkDir/formatted_nuc
```

P414

```bash
  Taxon_code=Pc_CR1
  Fasta_file=$(ls /data/scratch/armita/idris/gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

12-420

```bash
  Taxon_code=Pc_CR2
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/12420/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

15-13

```bash
  Taxon_code=Pc_CR3
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/15_13/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

15-7

```bash
  Taxon_code=Pc_CR4
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/15_7/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

2003-3

```bash
  Taxon_code=Pc_CR5
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/2003_3/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

4032

```bash
  Taxon_code=Pc_CR6
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/4032/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

4040

```bash
  Taxon_code=Pc_CR7
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/4040/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

P404

```bash
  Taxon_code=Pc_CR8
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/404/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```


P415

```bash
  Taxon_code=Pc_CR9
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/415/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

P416

```bash
  Taxon_code=Pc_CR10
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/416/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

P421

```bash
  Taxon_code=Pc_CR11
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/P421_v2/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

PC13/15

```bash
  Taxon_code=Pc_CR12
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/PC13_15/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

10300

```bash
  Taxon_code=Pc_CR13
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

Format fasta files - P. cactorum leather rot

11-40

```bash
  Taxon_code=Pc_LR1
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/11-40/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

17-21

```bash
  Taxon_code=Pc_LR2
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/17-21/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

Format fasta files - P. cactorum malus x domestica

62471

```bash
  Taxon_code=Pc_MD1
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/62471/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

P295

```bash
  Taxon_code=Pc_MD2
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/P295/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

R36/14

```bash
  Taxon_code=Pc_MD3
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.cactorum/R36_14/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

Format fasta files - P. idaei Rubus idaeu

P371

```bash
  Taxon_code=Pi_RI1
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.idaei/371/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

SCRP370

```bash
  Taxon_code=Pi_RI2
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.idaei/SCRP370/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

SCRP376

```bash
  Taxon_code=Pi_RI3
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/idris/gene_pred/final_incl_ORF/P.idaei/SCRP376/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  cat $Fasta_file | sed "s/>/>${Taxon_code}\|/g" > $WorkDir/formatted_nuc/"$Taxon_code".fasta
```

```bash
  cat $WorkDir/formatted_nuc/*_*.fasta > $WorkDir/formatted_nuc/goodCDS.fasta

  IsolateAbrv=Pcac_Pinf_publication
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogroupsTxt=$(ls $WorkDir/${IsolateAbrv}_orthogroups.txt)
  GoodProts=$WorkDir/formatted_nuc/goodCDS.fasta
  OutDir=$WorkDir/orthogroups_fasta_nuc
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogroupsTxt --fasta $GoodProts --out_dir $OutDir
```






<!-- # 6) Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.


### 6.1 ) Clade unique gene families

The genes unique to A. tenuissima apple pathotypes were identified within the orthology analysis.

First variables were set:
```bash
  WorkDir=analysis/orthology/orthomcl/At_Aa_Ag_all_isolates
  Orthogroups=$(ls $WorkDir/formatted/Results_May04/Orthogroups.txt)
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  # Braker_genes=
```

#### 6.1.a ) Orthologroups only containing A. tenuissima genes were extracted:

```bash
for num in 1; do
  # AtenUniq
AtenUniqDir=$WorkDir/A.tenuissima_unique
Uniq_Aten_groups=$AtenUniqDir/A.tenuissima_unique.txt
mkdir -p $AtenUniqDir
cat $Orthogroups | grep -e 'At_1' | grep -e 'At_2' | grep -e 'At_3' | grep -e 'At_4' | grep -e 'At_5' | grep -e 'At_6' | grep -e 'At_7' | grep -e 'At_8' | grep -v -e 'Aa' -e 'Ag' > $Uniq_Aten_groups
echo "The number of orthogroups unique to A.tenuissima apple pathotype are:"
cat $Uniq_Aten_groups | wc -l
echo "The following number genes from isolate 648 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_1' | wc -l
echo "The following number genes from isolate 1082 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_2' | wc -l
echo "The following number genes from isolate 1164 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_3' | wc -l
echo "The following number genes from isolate 24350 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_4' | wc -l
echo "The following number genes from isolate 648 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_5' | wc -l
echo "The following number genes from isolate 743 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_6' | wc -l
echo "The following number genes from isolate 1166 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_7' | wc -l
echo "The following number genes from isolate 1177 are contained in these orthogorups:"
cat $Uniq_Aten_groups | grep -o -e 'At_8' | wc -l
done
```

```
The number of orthogroups unique to A.tenuissima apple pathotype are:
51
The following number genes from isolate 648 are contained in these orthogorups:
51
The following number genes from isolate 1082 are contained in these orthogorups:
51
The following number genes from isolate 1164 are contained in these orthogorups:
51
The following number genes from isolate 24350 are contained in these orthogorups:
52
The following number genes from isolate 648 are contained in these orthogorups:
51
The following number genes from isolate 743 are contained in these orthogorups:
51
The following number genes from isolate 1166 are contained in these orthogorups:
51
The following number genes from isolate 1177 are contained in these orthogorups:
51
```

#### 6.1.b ) Orthologroups only containing A. arborescens genes were extracted:

```bash
for num in 1; do
  # AarbUniq
  AarbUniqDir=$WorkDir/A.arborescens_unique
  Uniq_Aarb_groups=$AarbUniqDir/A.arborescens_unique.txt
  mkdir -p $AarbUniqDir
  cat $Orthogroups | grep -e 'Aa_1' | grep -e 'Aa_2' | grep -e 'Aa_3' | grep -v -e 'At' -e 'Ag' > $Uniq_Aarb_groups
  echo "The number of orthogroups unique to A.arborescens are:"
  cat $Uniq_Aarb_groups | wc -l
  echo "The following number genes from isolate 675 are contained in these orthogorups:"
  cat $Uniq_Aarb_groups | grep -o -e 'Aa_1' | wc -l
  echo "The following number genes from isolate 97.0013 are contained in these orthogorups:"
  cat $Uniq_Aarb_groups | grep -o -e 'Aa_2' | wc -l
  echo "The following number genes from isolate 97.0016 are contained in these orthogorups:"
  cat $Uniq_Aarb_groups | grep -o -e 'Aa_3' | wc -l
done
```

```
The number of orthogroups unique to A.arborescens are:
150
The following number genes from isolate 675 are contained in these orthogorups:
155
The following number genes from isolate 97.0013 are contained in these orthogorups:
154
The following number genes from isolate 97.0016 are contained in these orthogorups:
156
```

#### 6.1.c ) Orthologroups only containing A. gaisen genes were extracted:

```bash
for num in 1; do
  # AgaiPathUniq
  AgaiUniqDir=$WorkDir/A.gaisen_unique
  Uniq_Agai_groups=$AgaiUniqDir/A.gaisen_unique.txt
  mkdir -p $AgaiUniqDir
  cat $Orthogroups | grep -e 'Ag_1' | grep -v  -e 'At' -e 'Aa' > $Uniq_Agai_groups
  echo "The number of orthogroups unique to A.gaisen pear pathotype pathotype are:"
  cat $Uniq_Agai_groups | wc -l
  echo "The following number genes from isolate 650 are contained in these orthogorups:"
  cat $Uniq_Agai_groups | grep -o -e 'Ag_1' | wc -l
done
```

```
The number of orthogroups unique to A.gaisen pear pathotype pathotype are:
327
The following number genes from isolate 650 are contained in these orthogorups:
329
```

#### 6.1.d ) Orthologroups only containing A. tenuissima non pathotype genes were extracted:

```bash
for num in 1; do
  # AtenNonPathUniq
  AtenNonPathUniqDir=$WorkDir/A.tenuissima_non_pathotype_unique
  Uniq_AtenNonPath_groups=$AtenNonPathUniqDir/A.tenuissima_non_pathotype_unique.txt
  mkdir -p $AtenNonPathUniqDir
  cat $Orthogroups | grep -e 'At_1' | grep -e 'At_2' | grep -e 'At_3' | grep -e 'At_4' | grep -v -e 'At_5' -e 'At_6' -e 'At_7' -e 'At_8' -e 'Aa' -e 'Ag' > $Uniq_AtenNonPath_groups
  echo "The number of orthogroups unique to A.tenuissima non-apple pathotype are:"
  cat $Uniq_AtenNonPath_groups | wc -l
  echo "The following number genes from isolate 648 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_1' | wc -l
  echo "The following number genes from isolate 1082 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_2' | wc -l
  echo "The following number genes from isolate 1164 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_3' | wc -l
  echo "The following number genes from isolate 24350 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_4' | wc -l
done
```

```
The number of orthogroups unique to A.tenuissima non-apple pathotype are:
0
The following number genes from isolate 648 are contained in these orthogorups:
0
The following number genes from isolate 1082 are contained in these orthogorups:
0
The following number genes from isolate 1164 are contained in these orthogorups:
0
The following number genes from isolate 24350 are contained in these orthogorups:
0
```

#### 6.1.d ) Orthologroups only containing A. tenuissima apple pathotype genes were extracted:

```bash
for num in 1; do
  # AtenPathUniq
  AtenPathUniqDir=$WorkDir/A.tenuissima_apple_pathotype_unique
  Uniq_AtenPath_groups=$AtenPathUniqDir/A.tenuissima_apple_pathotype_unique.txt
  mkdir -p $AtenPathUniqDir
  cat $Orthogroups | grep -e 'At_5' | grep -e 'At_6' | grep -e 'At_7' | grep -e 'At_8' | grep -v -e 'At_1' -e 'At_2' -e 'At_3' -e 'At_4' -e 'Aa' -e 'Ag' > $Uniq_AtenPath_groups
  echo "The number of orthogroups unique to A.tenuissima apple pathotype are:"
  cat $Uniq_AtenPath_groups | wc -l
  echo "The following number genes from isolate 635 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_5' | wc -l
  echo "The following number genes from isolate 743 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_6' | wc -l
  echo "The following number genes from isolate 1166 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_7' | wc -l
  echo "The following number genes from isolate 1177 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_8' | wc -l
done
```

```
The number of orthogroups unique to A.tenuissima apple pathotype are:
49
The following number genes from isolate 635 are contained in these orthogorups:
56
The following number genes from isolate 743 are contained in these orthogorups:
57
The following number genes from isolate 1166 are contained in these orthogorups:
72
The following number genes from isolate 1177 are contained in these orthogorups:
55
```

```bash
cat $Uniq_AtenPath_groups | grep -o -P "At_7\|\S+" | cut -f2 -d '|' > tmp.txt
Isolate='1166'
AnnotTab=$(ls gene_pred/annotation/A.*/$Isolate/${Isolate}_annotation_ncbi.tsv)
OutFile=$(echo $AnnotTab | sed 's/_annotation_ncbi.tsv/_apple_pathotype_unique_orthogroups.tsv/g')
for Gene in $(cat tmp.txt); do
  cat $AnnotTab | grep -w "^$Gene"
done > $OutFile
```


#### 6.1.e ) Orthologroups only containing genes from apple and pear pathotypes were extracted:

```bash
for num in 1; do
  # PathUniq
  PathUniqDir=$WorkDir/Pathotype_unique
  Uniq_Path_groups=$PathUniqDir/Path_unique.txt
  mkdir -p $PathUniqDir
  cat $Orthogroups | grep -e 'Ag' | grep -e 'At_5' | grep -e 'At_6' | grep -e 'At_7' | grep -e 'At_8' | grep -v -e 'At_1' -e 'At_2' -e 'At_3' -e 'At_4' -e 'Aa' > $Uniq_Path_groups
  echo "The number of orthogroups unique to apple and pear pathotypes are:"
  cat $Uniq_Path_groups | wc -l
  echo "The following number genes from isolate 650 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'Ag_1' | wc -l
  echo "The following number genes from isolate 648 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_5' | wc -l
  echo "The following number genes from isolate 743 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_6' | wc -l
  echo "The following number genes from isolate 1166 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_7' | wc -l
  echo "The following number genes from isolate 1177 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_8' | wc -l
done
```

```
The number of orthogroups unique to apple and pear pathotypes are:
48
The following number genes from isolate 650 are contained in these orthogorups:
62
The following number genes from isolate 648 are contained in these orthogorups:
58
The following number genes from isolate 743 are contained in these orthogorups:
58
The following number genes from isolate 1166 are contained in these orthogorups:
54
The following number genes from isolate 1177 are contained in these orthogorups:
60
``` -->
