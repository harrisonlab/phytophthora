


# 1 Find single copy busco genes in P.cactorum assemblies

Firstly, busco results from P414 in /data/scratch/armita/idris
was copied into the idris project folder

```bash
mkdir -p gene_pred/busco/P.cactorum/414/assembly
# cp -r /data/scratch/armita/idris/assembly/merged_SMARTdenovo_spades/P.cactorum/414/filtered/run_filtered_contigs_renamed gene_pred/busco/P.cactorum/414/assembly/.
cp -r /data/scratch/armita/idris/gene_pred/busco/P.cactorum/414/assembly/run_414_contigs_unmasked gene_pred/busco/P.cactorum/414/assembly
```

Create a list of all BUSCO IDs

```bash
cd /projects/oldhome/groups/harrisonlab/project_files/idris

# pushd /projects/oldhome/sobczm/bin/BUSCO_v1.22/fungi/hmms
OutDir=analysis/popgen/busco_phylogeny_stramenopiles
mkdir -p $OutDir
# BuscoDb="eukaryota_odb9"
BuscoDb="alveolata_stramenophiles_ensembl"
ls -1 /projects/oldhome/groups/harrisonlab/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt
```

For each busco gene create a folder and move all single copy busco hits from
each assembly to the folder.
Then create a fasta file containing all the aligned reads for each busco gene for
alignment later.

```bash
printf "" > analysis/popgen/busco_phylogeny_stramenopiles/single_hits.txt
for Busco in $(cat analysis/popgen/busco_phylogeny_stramenopiles/all_buscos_*.txt); do
echo $Busco
OutDir=analysis/popgen/busco_phylogeny_stramenopiles/$Busco
mkdir -p $OutDir
for Fasta in $(ls gene_pred/busco/*/*/assembly/*/single_copy_busco_sequences/$Busco*.fna | grep -e 'P.cactorum' -e 'P.idaei' | grep -e 'contigs_unmasked' -e 'filtered_contigs_renamed' -e 'LV007' | grep -v -e '414_old' -e '414_v2'); do
Strain=$(echo $Fasta | rev | cut -f5 -d '/' | rev)
Organism=$(echo $Fasta | rev | cut -f6 -d '/' | rev)
FileName=$(basename $Fasta)
cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
done
cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny_stramenopiles/single_hits.txt
done
```

If all isolates have a single copy of a busco gene, move the appended fasta to
a new folder

```bash
  OutDir=analysis/popgen/busco_phylogeny_stramenopiles/alignments
  mkdir -p $OutDir
  OrganismNum=$(cat analysis/popgen/busco_phylogeny_stramenopiles/single_hits.txt | cut -f2 | sort -nr | head -n1)
  for Busco in $(cat analysis/popgen/busco_phylogeny_stramenopiles/all_buscos_*.txt); do
  echo $Busco
  HitNum=$(cat analysis/popgen/busco_phylogeny_stramenopiles/single_hits.txt | grep "$Busco" | cut -f2)
  if [ $HitNum == $OrganismNum ]; then
    # cp analysis/popgen/busco_phylogeny_stramenopiles/$Busco/"$Busco"_appended.fasta $OutDir/.
    cat analysis/popgen/busco_phylogeny_stramenopiles/$Busco/"$Busco"_appended.fasta \
    | sed "s/$Busco://g" \
    | sed "s/genome.ctg.fa://g" \
    | sed "s/_contigs_unmasked.fa//g" \
    | sed -E "s/:.*//g" \
    | sed "s/P.cactorum_//g" \
    | sed 's/Phytophthora_cactorum-//g' \
    | sed 's/.fa//g' \
    | tr '.,:' '_' \
    > $OutDir/"$Busco"_appended.fasta
  fi
  done
```

Submit alignment for single copy busco genes with a hit in each organism


```bash
  AlignDir=analysis/popgen/busco_phylogeny_stramenopiles/alignments
  CurDir=$PWD
  cd $AlignDir
  ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
  qsub $ProgDir/sub_mafft_alignment.sh
  cd $CurDir
```


```bash
  OutDir=analysis/popgen/busco_phylogeny_stramenopiles/trimmed_alignments
  mkdir -p $OutDir
  for Alignment in $(ls analysis/popgen/busco_phylogeny_stramenopiles/alignments/*_appended_aligned.fasta); do
    TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
    trimal -in $Alignment -out $OutDir/$TrimmedName -automated1
  done
```



```bash
for Alignment in $(ls analysis/popgen/busco_phylogeny_stramenopiles/trimmed_alignments/*aligned_trimmed.fasta); do
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
while [ $Jobs -gt 2 ]; do
sleep 2s
# printf "."
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Prefix
Prefix=$(basename $Alignment | cut -f1 -d '_')
OutDir=analysis/popgen/busco_phylogeny_stramenopiles/RAxML/$Prefix
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/phylogenetics
qsub $ProgDir/sub_RAxML.sh $Alignment $Prefix $OutDir
done
```

Run Astral to build a consensus phylogeny from a collective set of
"best phylogenies" from each BUSCO locus

* Note - "Recent versions of ASTRAL output a branch support value even without bootstrapping. Our analyses have revealed that this form of support is more reliable than bootstrapping (under the conditions we explored). Nevertheless, you may want to run bootstrapping as well."

Tutorial tips:
https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#running-with-unresolved-gene-trees


```bash
OutDir=analysis/popgen/busco_phylogeny_stramenopiles/ASTRAL
mkdir -p $OutDir
cat analysis/popgen/busco_phylogeny_stramenopiles/RAxML/*/RAxML_bestTree.* > $OutDir/Pcac_phylogeny.appended.tre
# InTree=$(ls /projects/oldhome/armita/prog/Astral/Astral/test_data/song_primates.424.gene.tre)
# -
# Trimm back branches that have less than 10% bootstrap support for each tree
# in the given file
# -
/projects/oldhome/armita/prog/newick_utilities/newick_utils/src/nw_ed $OutDir/Pcac_phylogeny.appended.tre 'i & b<=10' o > $OutDir/Pcac_phylogeny.appended.trimmed.tre
# -
# Calculate combined tree
# -
ProgDir=/projects/oldhome/armita/prog/Astral/Astral
java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -i $OutDir/Pcac_phylogeny.appended.tre -o $OutDir/Pcac_phylogeny.consensus.tre | tee 2> $OutDir/Pcac_phylogeny.consensus.log
java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -q $OutDir/Pcac_phylogeny.consensus.tre -i $OutDir/Pcac_phylogeny.appended.tre -o $OutDir/Pcac_phylogeny.consensus.scored.tre 2> $OutDir/Pcac_phylogeny.consensus.scored.log
```

GGtree was used to make a plot.

* Note- Tips can be found here: https://bioconnector.org/r-ggtree.html

The consensus tree was downloaded to my local machine

* Note - I had to import into geneious and export again in newick format to get around polytomy branches having no branch length.
* Terminal branch lengths are meanlingless from ASTRAL and should all be set to an arbitrary value. This will be done by geneious (set to 1), but it also introduces a branch length of 2 for one isolate that needs to be corrected with sed

```bash
# cat Pcac_phylogeny_stramenopiles.consensus.scored.geneious.tre | sed 's/:2/:1/g' > Pcac_phylogeny_stramenopiles.consensus.scored.geneious2.tre
cat /Users/armita/OneDrive\ -\ University\ of\ Greenwich/Armitage/Projects/Phytophthora\ EMR/Pcac_phylogeny.consensus.scored.tre.newick | sed 's/:2/:1/g' > /Users/armita/OneDrive\ -\ University\ of\ Greenwich/Armitage/Projects/Phytophthora\ EMR/Pcac_phylogeny.consensus.scored.tre.geneious_mod.newick
```


```r
# setwd("/Users/armita/Downloads/Pc/alignments/ASTRAL")
setwd("/Users/armita/OneDrive\ -\ University\ of\ Greenwich/Armitage/Projects/Phytophthora\ EMR")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)

# tree <- read.tree("Pcac_phylogeny.consensus.scored_geneious2.tre")
tree <- read.tree("Pcac_phylogeny.consensus.scored.tre.geneious_mod.newick")

tree$edge.length[tree$edge.length == 1] <- 0

mydata <- read.csv("/Users/armita/Downloads/Pc/alignments/ASTRAL/traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$label
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]


t <- ggtree(tree) # Core tree
t <- t + geom_treescale(offset=-0.5, fontsize = 3) # Add scalebar

# labels by values in another df
t <- t %<+% mydata
tips <- data.frame(t$data)
tips$label <- tips$Isolate
t <- t + geom_tiplab(data=tips, size=3, hjust=0, offset = +0.9, align=T, linetype = NULL)
tips2 <- data.frame(t$data)
tips2$label <- tips2$Host
t <- t + geom_tiplab(data=tips2, size=3, hjust=0, offset = +0.5, fontface = "italic", align=T)

# Format nodes by values
nodes <- data.frame(t$data)
#nodes <- nodes[!nodes$isTip,]
nodes$label <- as.numeric(nodes[nodes$label,])
as.numeric(nodes$label)
#nodes$label[nodes$label < 0.80] <- ''
nodes$support[nodes$isTip] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label > 0.80)] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label < 0.80)] <- 'unsupported'
nodes$support[(!nodes$isTip) & (nodes$label == '')] <- 'supported'
t <- t + aes(linetype=nodes$support)
nodes$label[nodes$label > 0.80] <- ''
t <- t + geom_nodelab(data=nodes, size=2, hjust=-0.05) # colours as defined by col2rgb

# Add in a further set of labels
# tree_mod <- tree
# tree_mod$tip.label <- mydata$Source
# t <- t + geom_tiplab(data=tree_mod, aes(label=label), size=2, offset = +1)

# Annotate a clade with a bar line
t <- t + geom_cladelabel(node=42, label='P. i', align=T, colour='black', offset=+1.75, fontface = "italic")
t <- t + geom_cladelabel(node=24, label='P. c', align=T, colour='black', offset=+1.75, fontface = "italic", fontface = "italic")

# Save as PDF and force a 'huge' size plot
ggsave("Pcac_phylogeny.pdf", width =20, height = 20, units = "cm", limitsize = FALSE)
````

The tree was also plotted using a matrix to display additional isoalte information.

```r

setwd("/Users/armita/OneDrive\ -\ University\ of\ Greenwich/Armitage/Projects/Phytophthora\ EMR")
​
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)

#===============================================================================
#       Import datasets
#===============================================================================

tree <- read.tree("/Users/armita/Google\ Drive/Pc\ comparative\ genomics/Figures/raw\ trees/Pcac_phylogeny.consensus.scored.tre.newick")
mydata <- read.csv("/Users/armita/Google\ Drive/Pc\ comparative\ genomics/Figures/raw\ trees/Phytophthora_traits2.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$label
​
## Remove terminal branch lengths
tree$edge.length[tree$edge.length == 1] <- 0
tree$edge.length[tree$edge.length == 2] <- 0


#===============================================================================
#       Prepare dataframe to rename isolates
#===============================================================================

mydata <- mydata[match(tree$tip.label,rownames(mydata)),]
​
​

#===============================================================================
#       Build the tree
#===============================================================================

# assign tree to variable, core tree
q <- ggtree(tree)
​
# add scalebar
q <- q + geom_treescale(offset=-0.5, fontsize = 3)
​
# Add bootstraps less than 100
d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label < 100,]
​
​
heatmap<-mydata[c(4,5,7,8,9)]
​
​
# Plot the tree and save to pdf
q2<- q +  geom_text2(data=d, aes(label=label), size=2.5, adj = -0.2, color="blue")
​
# create labels to rename isolates by values in another dataframe
q <- q %<+% mydata
tips <- data.frame(q$data)
tips$label <- tips$Isolate
​
q3 <- q2 %<+% mydata + geom_tiplab(data=tips, align=TRUE, linesize=0.5, offset=0.5, aes(fill=factor(Clade)),
              color="black",
              geom = "label",
              label.size = 0)

q4<-gheatmap(q3, heatmap, width=0.15, offset=1, font.size=3.5, colnames_angle=90, colnames_position='top',colnames_offset_y=1)  +
scale_fill_manual("",values=c('1'='#e9e2d0','2'='#dedede','Fragaria (crown)'='#440154FF','Fragaria (fruit)'='#443983','Malus'='#238A8DFF','Rubus'='#FDE725FF','Fagus'='#55C667FF','USA'='#0D0887','UK'='#7E03A7','Netherlands'='#CC4678','Norway'='#F1701F','Sweden'='#FEB131','Present'='#2b2a2a','Absent'='#fbf9fa','Premature stop codon'='#5c636e','More than a premature stop codon'='grey'),
  breaks=c('Fragaria (crown)','Fragaria (fruit)','Malus','Fagus','Rubus','UK','USA','Netherlands','Norway','Sweden','Present','Absent','Premature stop codon','More than a premature stop codon'),
  labels=c(bquote(italic('Fragaria')~'x'~italic('ananassa')~'(crown)'),bquote(italic('Fragaria')~'x'~italic('ananassa')~'(fruit)'),bquote(italic('Malus')~'x'~italic('domestica')),bquote(italic('Fagus sylvatica')),bquote(italic('Rubus idaeus')),bquote('U.K.'),bquote('U.S.A.'),bquote('Netherlands'),bquote('Norway'),bquote('Sweden'),bquote('Present'),bquote('Absent'),bquote('Single SNP, premature stop codon'),bquote('Multiple variations, premature stop codon')))+ #bquote, left hand justified
  theme(legend.position="bottom", legend.title=element_blank()) + guides(colour=FALSE)
pdf(file = "/Users/armita/Google\ Drive/Pc\ comparative\ genomics/Figures/raw\ trees/Phytophthora_tree.pdf",width=9,height=10)
q4
dev.off()
​

```





<!--


```bash
# For closely related organisms (same species etc.): identify genes with high nucleotide diversity (Pi) and average number of pairwise differences, medium number of segregating sites
# (avoid alignments with low homology and lots of phylogenetically uninformative singletons).
# For analyses involving cross-species comparisons involving highly diverged sequences with high nucleotide diversity
# (e.g. 0.1<Pi<0.4), looking for genes with the lowest number of segregating sites.
AlignDir=analysis/popgen/busco_phylogeny_stramenopiles/alignments
CurDir=$PWD
cd $AlignDir

# pip install dendropy --user
for Alignment in $(ls *aligned.fasta); do
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
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

Visually inspect the alignments of selected genes (genes_selected_for_phylogeny.txt) to be used in
constructing the phylogenies and trim them as necessary in MEGA7.
Copy the relevant trimmed alignment FASTA files into

```bash
  # mkdir $CurDir/beast_runs/candidates/select/trimmed
```


##PartitionFinder (nucleotide sequence evolution model)

```bash
cd analysis/popgen/busco_phylogeny_stramenopiles/phylogeny

config_template=/projects/oldhome/sobczm/bin/PartitionFinder1.1.1/partition_finder.cfg
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
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
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
scp -r cluster:/projects/oldhome/groups/harrisonlab/project_files/idris/analysis/popgen/busco_phylogeny_stramenopiles/phylogeny .
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
ClusterDir=/projects/oldhome/groups/harrisonlab/project_files/idris/analysis/popgen/busco_phylogeny_stramenopiles/phylogeny
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
* Tracer: /projects/oldhome/sobczm/bin/beast/Tracer_v1.6/bin/tracer
some runs may never converge)
* Store every: 10000

```bash

cd /projects/oldhome/groups/harrisonlab/project_files/idris


for File in $(ls analysis/popgen/busco_phylogeny_stramenopiles/phylogeny/good_alignments/*_appended_aligned/analysis/best_scheme.txt); do
Busco=$(echo $File | cut -f6 -d '/' | cut -f1 -d '_')
Model=$(cat $File | grep -A1 'Best Model' | tail -n1 | cut -f2 -d '|')
printf "$Busco\t$Model\n"
done

# Edit NEXUS files:
for Nexus in $(ls analysis/popgen/busco_phylogeny_stramenopiles/phylogeny/good_alignments/*_appended_aligned/*_appended_aligned.NEXUS); do
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
NexusFiles=$(ls analysis/popgen/busco_phylogeny_stramenopiles/phylogeny/good_alignments/*_appended_aligned/*.NEXUS | sed -e 's/^/ -nex /g' | tr -d '\n')
OutFile=$(echo $Nexus | sed 's/.NEXUS/.xml/g')
ProgDir=/projects/oldhome/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/beauti -template StarBeast.xml $NexusFiles




qlogin -pe smp 8
InXML=analysis/popgen/busco_phylogeny_stramenopiles/phylogeny/Pcac_beauti_starBEAST2.xml
OutDir=$(dirname $InXML)"/BEAST4"
mkdir -p $OutDir
ProgDir=/projects/oldhome/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/beast -threads 8 -prefix $OutDir $InXML > $OutDir/log.txt
# java -Djava.library.path="C:\Program Files (x86)\Common Files\libhmsbeagle-1.0" -jar "/BEAST175/lib/beast.jar"

#After the run, check convergence with Tracer, summarise the final tree with TreeAnnotator
for Tree in $(ls $OutDir/*.trees); do
BurnIn=10 # percentage of states to be considered as burnin
SumTree=$(echo $Tree | sed 's/.trees/_summary.tree/g')
ProgDir=/projects/oldhome/sobczm/bin/beast/BEASTv2.4.2/bin
$ProgDir/treeannotator -heights median -burnin $BurnIn $Tree $SumTree
done

#Visualise and beautify the final tree (suffix "summary") with FigTree
FigTree=/projects/oldhome/sobczm/bin/FigTree_v1.4.2/bin/figtree
$FigTree

```
-->
