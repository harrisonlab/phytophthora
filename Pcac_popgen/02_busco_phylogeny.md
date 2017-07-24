


# 1 Find single copy busco genes in P.cactorum assemblies


Create a list of all BUSCO IDs

```bash
cd /home/groups/harrisonlab/project_files/idris

# pushd /home/sobczm/bin/BUSCO_v1.22/fungi/hmms
OutDir=analysis/popgen/busco_phylogeny
mkdir -p $OutDir
BuscoDb="eukaryota_odb9"
ls -1 /home/groups/harrisonlab/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt
```

For each busco gene create a folder and move all single copy busco hits from
each assembly to the folder.
Then create a fasta file containing all the aligned reads for each busco gene for
alignment later.

```bash
printf "" > analysis/popgen/busco_phylogeny/single_hits.txt
for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
echo $Busco
OutDir=analysis/popgen/busco_phylogeny/$Busco
mkdir -p $OutDir
for Fasta in $(ls gene_pred/busco/*/*/assembly/run_contigs_min_500bp*/single_copy_busco_sequences/$Busco*.fna | grep -v -w '414'); do
Strain=$(echo $Fasta | rev | cut -f5 -d '/' | rev)
Organism=$(echo $Fasta | rev | cut -f6 -d '/' | rev)
FileName=$(basename $Fasta)
cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
done
cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny/single_hits.txt
done
```


If all isolates have a single copy of a busco gene, move the appended fasta to
a new folder

```bash
  OutDir=analysis/popgen/busco_phylogeny/alignments
  mkdir -p $OutDir
  OrganismNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | cut -f2 | sort -nr | head -n1)
  for Busco in $(cat analysis/popgen/busco_phylogeny/all_buscos_*.txt); do
  echo $Busco
  HitNum=$(cat analysis/popgen/busco_phylogeny/single_hits.txt | grep "$Busco" | cut -f2)
  if [ $HitNum == $OrganismNum ]; then
    cp analysis/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta $OutDir/.
  fi
  done
```

Submit alignment for single copy busco genes with a hit in each organism


```bash
  AlignDir=analysis/popgen/busco_phylogeny/alignments
  CurDir=$PWD
  cd $AlignDir
  ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
  qsub $ProgDir/sub_mafft_alignment.sh
  cd $CurDir
```



```bash
# For closely related organisms (same species etc.): identify genes with high nucleotide diversity (Pi) and average number of pairwise differences, medium number of segregating sites
# (avoid alignments with low homology and lots of phylogenetically uninformative singletons).
# For analyses involving cross-species comparisons involving highly diverged sequences with high nucleotide diversity
# (e.g. 0.1<Pi<0.4), looking for genes with the lowest number of segregating sites.
AlignDir=analysis/popgen/busco_phylogeny/alignments
CurDir=$PWD
cd $AlignDir

# pip install dendropy --user
for Alignment in $(ls *aligned.fasta); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
python $ProgDir/calculate_nucleotide_diversity.py $Alignment
Busco=$(echo $Alignment | cut -f1 -d '_')
mv sequence_stats.txt "$Busco"_seqeunce_stats.txt
mv excel_stats.txt "$Busco"_excel_stats.txt
mkdir -p ../phylogeny
## Copy FASTA files of the aligments into a new directory
cp $Alignment ../phylogeny/.
done

cd $CurDir
```

Visually inspect the alignments of select genes (genes_selected_for_phylogeny.txt) to be used in
constructing the phylogenies and trim them as necessary in MEGA7.
Copy the relevant trimmed alignment FASTA files into

```bash
  # mkdir $CurDir/beast_runs/candidates/select/trimmed
```


##PartitionFinder (nucleotide sequence evolution model)

```bash
cd analysis/popgen/busco_phylogeny/phylogeny

config_template=/home/sobczm/bin/PartitionFinder1.1.1/partition_finder.cfg
ct=$(basename "$config_template")

mkdir NEXUS

# prepare directory for PartitionFinder run:
for f in $(ls *fasta); do
sed -i 's/:/_/g' $f
c="$(cat $f | awk 'NR%2==0' | awk '{print length($1)}' | head -1)"
p="${f%.fasta}.phy"
n="${f%.fasta}.NEXUS"
dir="${f%.fasta}"

mkdir $dir
cp $config_template $dir/.

# Substitute the name of the alignment file and the sequence length in the config file to become correct for the current run.
sed -i 's,^\(alignment = \).*,\1'"$p;"',' $dir/$ct
sed -i 's,^\(Gene1_pos1 = \).*,\1'"1-$c\\\3;"',' $dir/$ct
sed -i 's,^\(Gene1_pos2 = \).*,\1'"2-$c\\\3;"',' $dir/$ct
sed -i 's,^\(Gene1_pos3 = \).*,\1'"3-$c\\\3;"',' $dir/$ct

# Convert FASTA to phylip for the Partition Finder run
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
$ProgDir/fasta2phylip.pl $f>$p
mv $p $dir

# Convert FASTA to NEXUS for the BEAST run
$ProgDir/Fasta2Nexus.pl $f>$dir/$n

#Problems running PartitionFinder on the cluster. May have to be run locally on your Mac or Windows machine.
# qsub $ProgDir/sub_partition_finder.sh $dir
done
```

Partition finder wasnt run on the cluster. As such fasta alignment files were
downloaded to the local machine where partitionfinder was run
patritionfinder2 was downloaded from:
http://www.robertlanfear.com/partitionfinder/

and the anaconda libraries to support it were downloaded from:
https://www.continuum.io/downloads#macos


copy the fasta files and the partitionfinder config files to
your local computer

```bash
cd Users/armita/Downloads
scp -r cluster:/home/groups/harrisonlab/project_files/idris/analysis/popgen/busco_phylogeny/phylogeny .
```

Alignments were loaded into Geneious where they were visualised and manually sorted into
three categories:
* Good - All sequences present no trimming needed
* Trim - All sequences present short regions may need trimming from the beginning / end of the alignment before use in phylogenetics
* Bad - a region of one or more sequences is missing or the sequences / alignment is not appropriate for phylogenetics

These alignments were then exported from Geneious into the following folders:

```bash
cd Users/armita/Downloads/phylogeny
mkdir good_alignments
mkdir trim_alignments
mkdir bad_alignments
```

Alignments within the "good alignments" directory were taken forward for further
analysis

```bash
  for Dir in $(ls -d *_alignments); do
    for Alignment in $(ls $Dir/*_appended_aligned.phy); do
      Prefix=$(echo $Alignment | cut -f2 -d '/' | sed 's/.phy//g')
      echo $Prefix
      cp $Prefix/$Prefix.NEXUS $Dir/$Prefix/.
      cp -r $Prefix $Dir/.
      /Users/armita/anaconda2/bin/python ../partitionfinder-2.1.1/PartitionFinder.py $Dir/$Prefix --no-ml-tree --force-restart
    done
  done > log.txt
```


Upload partition models back to the cluster:

```bash
ClusterDir=/home/groups/harrisonlab/project_files/idris/analysis/popgen/busco_phylogeny/phylogeny
scp -r bad_alignments cluster:$ClusterDir/.
```


## Preparing to run BEAST


Using trimmed FASTA alignments and nucleotide substitution models identified with PartitionFinder:
create an XML input file using BEAUTi, with StarBeast template.

Prepare a 30 loci dataset, in addition to a 5 loci subset to compare convergence.

Run after qlogin into a worker node (BEAST does not find BEAGLE libraries when using qsub -
as the BEAST package is quite fiddly, may troubleshoot it later when necessary.

StarBeast settings used here:
* Substitution rate: default HKY
* Strict clock
* Species Tree Population Size: Linear with constant root
* Yule prior on species tree
* Chain length: 300 million (this may vary, change run convergence with Tracer during the run to establish the number of iterations required
* Tracer: /home/sobczm/bin/beast/Tracer_v1.6/bin/tracer
some runs may never converge)
* Store every: 10000

```bash

cd /home/groups/harrisonlab/project_files/idris


for File in $(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/analysis/best_scheme.txt); do
Busco=$(echo $File | cut -f6 -d '/' | cut -f1 -d '_')
Model=$(cat $File | grep -A1 'Best Model' | tail -n1 | cut -f2 -d '|')
printf "$Busco\t$Model\n"
done

# Edit NEXUS files:
for Nexus in $(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/*_appended_aligned.NEXUS); do
  sed -i -r "s/^.*_P\./P./g" $Nexus
  sed -i -r "s/_contig.*\t/\t/g" $Nexus
  sed -i -r "s/_NODE.*\t/\t/g" $Nexus
done

# OUtputs of partitionfinder were used to set models
# of DNA evolution in Beauti, as described on:
# http://www.robertlanfear.com/partitionfinder/faq/#toc-beast
# CHain length was modified from 10000000 to 500000000 as determined
# by a first run of beast where tracer reported the estimated sasmple size to be below 100 (3) - increase by 50 fold.

# Run Beauti
NexusFiles=$(ls analysis/popgen/busco_phylogeny/phylogeny/good_alignments/*_appended_aligned/*.NEXUS | sed -e 's/^/ -nex /g' | tr -d '\n')
OutFile=$(echo $Nexus | sed 's/.NEXUS/.xml/g')
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/beauti -template StarBeast.xml $NexusFiles




qlogin -pe smp 8
InXML=analysis/popgen/busco_phylogeny/phylogeny/Pcac_beauti_starBEAST2.xml
OutDir=$(dirname $InXML)"/BEAST4"
mkdir -p $OutDir
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/beast -threads 8 -prefix $OutDir $InXML > $OutDir/log.txt
# java -Djava.library.path="C:\Program Files (x86)\Common Files\libhmsbeagle-1.0" -jar "/BEAST175/lib/beast.jar"

#After the run, check convergence with Tracer, summarise the final tree with TreeAnnotator
for Tree in $(ls $OutDir/*.trees); do
BurnIn=10 # percentage of states to be considered as burnin
SumTree=$(echo $Tree | sed 's/.trees/_summary.tree/g')
ProgDir=/home/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/treeannotator -heights median -burnin $BurnIn $Tree $SumTree
done

#Visualise and beautify the final tree (suffix "summary") with FigTree
FigTree=/home/sobczm/bin/FigTree_v1.4.2/bin/figtree
$FigTree

```
