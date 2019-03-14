# 04_popgenome_analysis.md

This documents the commands used in the analysis of SNP .vcf files generated in
01_Pcac_alignments.md


## Identify the population structure within populations

This is done using faststructure.

Structure was investigated within all P. cactorum isolates and within
the crown rot P. cactorum population.

Commands are documented in

```
popgenome_scripts/fastStructure.md
```



<!--
## Calculate Fst between two populations

An FST of zero indicates no divergence between populations, whereas an FST of one indicates complete isolation of populations. Negative Fst values should be effectively seen as zero values.



```bash
# Remove samples not of interest
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
# ExcludeList="11-40 12420 15_13 15_7 17-21 2003_3 371 4032 404 4040 414 415 416 62471 P295 P421 PC13_15 R36_14 SCRP370 SCRP376"
ExcludeList="11-40 17-21 371 SCRP370 SCRP376"
Prefix=apple_and_cown_rot_vs_414
OutDir=analysis/popgen/SNP_calling/fst/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
  $VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered

Pop1=analysis/popgen/SNP_calling/PcacFa_isolates.txt
printf "12420\n15_13\n15_7\n2003_3\n4032\n404\n4040\n414\n415\n416\nP421\nPC13_15" > $Pop1
Pop2=analysis/popgen/SNP_calling/PcacMd_isolates.txt
printf "62471\nP295\nR36_14" > $Pop2
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/"$Prefix"_filtered.recode.vcf --weir-fst-pop $Pop1 --weir-fst-pop $Pop2 --out $OutDir/${Prefix}_Fst
```

```
After filtering, kept 15 out of 15 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.83985
Weir and Cockerham weighted Fst estimate: 0.93549
After filtering, kept 26952 out of a possible 26952 Sites
```
-->

<!--
## Calculate linkage diseaquilibrium within a population

Linkage disequilibrium:
--hap-r2
--geno-r2
--geno-chisq

```bash

  Vcf=$(ls analysis/popgen/SNP_calling/1166_contigs_unmasked_filtered_no_errors_syn.vcf)
  ExcludeList="650 648 97.0013 97.0016"
  Prefix=tenuissima_vs_1166
  OutDir=analysis/popgen/SNP_calling/$Prefix
  mkdir -p $OutDir
  VcfLib=/home/sobczm/bin/vcflib/bin
  $VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
  VcfTools=/home/sobczm/bin/vcftools/bin
  $VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
  VcfTools=/home/sobczm/bin/vcftools/bin
  $VcfTools/vcftools --vcf $OutDir/"$Prefix"_filtered.recode.vcf --hap-r2 --ld-window-bp 50000 --out $OutDir/${Prefix}_ld_window_50kb
```
 -->



## Make a slimmed vcf file of crown rot and apple isolates

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="11-40 17-21"
Prefix=crown-rot_apple_isolates
OutDir=analysis/popgen/SNP_calling/summary_stats3/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```


https://cran.r-project.org/web/packages/PopGenome/PopGenome.pdf
http://popgenome.weebly.com/uploads/5/1/5/7/51573093/whole_genome_analyses_using_vcf_files.pdf


##Set inital variables

```bash
CurDir=/data/scratch/armita/idris
Prefix=crown-rot_apple_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats3/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen/popgenome_scripts
```

```
In order to calculate different statistics in Popgenome, the WorkDir has to be arranged in a particular way.
The WorkDir directory should contain two folders.
Folder No. 1: named "gff", contains GFF files for all the contigs output from the split_gff_contig.sh script
Folder No. 2: named "contigs", contains subfolders, each subfolder named with exact contig name and containing one individual contig FASTA file, also named with exact contig name, as output from vcf_to_fasta.py
```

##An example on how to create this directory structure

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats3/$Prefix
mkdir -p $WorkDir/all
mkdir -p $WorkDir/gff
mkdir -p $WorkDir/contigs
```

###copy the "gff" folder containing gff files

```bash
CurDir=/data/scratch/armita/idris
Prefix=crown-rot_apple_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats3/$Prefix
Gff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.gff3)

cd $WorkDir/gff
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/summary_stats
$ProgDir/split_gff_contig.sh $CurDir/$Gff
cd $CurDir
```

### Convert vcf files into fasta alignments

The vcf file needs to be converted into fasta alignment files. One file is
produced per reference contig

For the ploidy set <1|2|3>:
* 1 - haploid input
* 2 - diploid input, output as two separate phased haplotypes for each ind.
* 3 - diploid input, output as one sequence with ambiguity codes for each ind.

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats3/$Prefix

Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
Vcf=$(ls analysis/popgen/SNP_calling/summary_stats3/crown-rot_apple_isolates/crown-rot_apple_isolates_filtered.recode.vcf)
Ploidy=2

cd $WorkDir/contigs
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/summary_stats
python $ProgDir/vcf_to_fasta.py $CurDir/$Vcf $CurDir/$Reference $Ploidy
cd $CurDir
```


### This folder contains only contig FASTA files
### So create a new "contigs" directory to hold those files:

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats3/$Prefix
Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
cp $Reference $WorkDir/.
```

###Next step: in the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
CurDir=/data/scratch/armita/idris
Prefix=crown-rot_apple_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats3/$Prefix
cd $WorkDir/contigs
for File in $(ls *.fasta); do
    Prefix=${File%.fasta}
    mkdir $Prefix
    mv $File $Prefix/.
done
cd $CurDir
```

### Test if all contigs have a matching gff and remove any which do not

```bash
CurDir=/data/scratch/armita/idris
Prefix=crown-rot_apple_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats3/$Prefix
cd $WorkDir
for File in $(ls $PWD/contigs/*/*.fasta); do
    Prefix=$(basename "$File" .fasta)
    expected_gff="$PWD/gff/${Prefix}.gff"
    if [ ! -f "$expected_gff" ]; then
      echo "dirname $File"
      rm -rf $(dirname $File)
    fi
done
cd $CurDir
```

## Calculate diversity and nutrality stats


Nucleotide diversity stats were calculated on my local machine

the file structure was downloaded using the commands:

(on local computer)
```bash
  cd /Users/armita/Downloads/Pc
  scp -r cluster:/data/scratch/armita/idris/analysis/popgen/SNP_calling/summary_stats3 .
  cd summary_stats3
```

Commands used to do these analyses are documented in:
popgenome_scripts/calculate_nucleotide_diversity.R

### Pi

Outputs were analysed using:

```R
# install.packages("symbols")
library(ggplot2)
library("symbols")
Fa_Pi <- read.delim("/Users/armita/Downloads/Pc/summary_stats3/crown-rot_apple_isolates/genome_Pcac_Fa_Pi_per_gene_all.txt", header=FALSE)
Md_Pi <- read.delim("/Users/armita/Downloads/Pc/summary_stats3/crown-rot_apple_isolates/genome_Pcac_Md_Pi_per_gene_all.txt", header=FALSE)
Ri_Pi <- read.delim("/Users/armita/Downloads/Pc/summary_stats3/crown-rot_apple_isolates/genome_Pcac_Ri_Pi_per_gene_all.txt", header=FALSE)
df1 <- data.frame(Fa_Pi$V3)
df1$Pi_Md <- Md_Pi$V3
df1$Pi_Ri <- Ri_Pi$V3
colnames(df1) <- c('P.cac Fxa', 'P.cac Mxd', 'P.idaei Ri')
library(reshape)
df2 <- melt(df1)
df3 <- data.frame(df2[df2$value != 0, ])

p <- ggplot(df3, aes( x = variable, y = value)) +
  geom_boxplot()

p <- p + labs(x = '', y = "Nucleotide diversity (\u03C0)")
ggsave('Pcac_Pi_boxplot.jpg', p)


Fa_Pi <- read.delim("/Users/armita/Downloads/Pc/summary_stats3/crown-rot_apple_isolates/genome_Pcac_Fa_Pi_n_s_per_gene_all.txt", header=FALSE)
Md_Pi <- read.delim("/Users/armita/Downloads/Pc/summary_stats3/crown-rot_apple_isolates/genome_Pcac_Md_Pi_n_s_per_gene_all.txt", header=FALSE)
Ri_Pi <- read.delim("/Users/armita/Downloads/Pc/summary_stats3/crown-rot_apple_isolates/genome_Pcac_Ri_Pi_n_s_per_gene_all.txt", header=FALSE)
df1 <- data.frame(Fa_Pi$V3)
df1$Pi_Md <- Md_Pi$V3
df1$Pi_Ri <- Ri_Pi$V3
colnames(df1) <- c('P.cac Fxa', 'P.cac Mxd', 'P.idaei Ri')
library(reshape)
df2 <- melt(df1)
df3 <- data.frame(df2[df2$value != 0, ])
df3 <- na.omit(df3)
p <- ggplot(df3[which(df2$value>0)],, aes(x = variable, y = value)) +
  geom_boxplot()

p <- p + labs(x = '', y = "Nucleotide diversity (\u03C0)")
ggsave('Pcac_Pi_nonsyn_boxplot.jpg', p)
```


Calculate Pi for different groups of genes:

```R


```

### dxy pairwise nucelotide differences


```R
# install.packages("symbols")
library(ggplot2)
library("symbols")
pop1_vs_pop2 <- read.delim("~/Downloads/popstats/summary_stats/genome_pop1_vs_pop2_dxy_per_gene.txt", header=FALSE)
pop1_vs_pop3 <- read.delim("~/Downloads/popstats/summary_stats/genome_pop1_vs_pop3_dxy_per_gene.txt", header=FALSE)
pop2_vs_pop3 <- read.delim("~/Downloads/popstats/summary_stats/genome_pop2_vs_pop3_dxy_per_gene.txt", header=FALSE)
df1 <- data.frame(pop1_vs_pop2$V3)
df1$pop1_vs_pop3 <- pop1_vs_pop3$V3
df1$pop2_vs_pop3 <- pop2_vs_pop3$V3
colnames(df1) <- c('P.cac Fxa vs Mxd', 'P.cac Fxa vs P.idaei', 'P.cac Mxd vs P.idaei')
library(reshape)
df2 <- melt(df1)
df3 <- data.frame(df2[df2$value != 0, ])
p <- ggplot(df3, aes(x = variable, y = value)) +
  geom_boxplot()
p <- p + labs(x = '', y = "Pairwise diversity (dxy)")
ggsave('Pcac_dxy_boxplot.jpg', p)
```
<!--
The R script used below is custom-made for each run (see first few lines of it).
It requires custom definition of populations, and individual assignment to them.
The example below calculates nucleotide diversity within (Pc) and between (Pc) populations.
Other ProgDir (sub_calculate_neutrality_stats.sh) are used in analogous manner.


Calculate diversity stats on all P. cactorum isolates

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats
cd $WorkDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen/popgenome_scripts
qsub $ProgDir/sub_calculate_nucleotide_diversity.sh
cd $PWD
```



Calculate Fst between P. cactorum populations

```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen/popgenome_scripts
qsub $ProgDir/sub_calculate_fst.sh
```



```bash
ProgDir2=/home/armita/git_repos/emr_repos/ProgDir/phytophthora_fragariae/popgen_analysis/popgenome_ProgDir
qsub $ProgDir2/sub_calculate_nucleotide_diversity.sh
qsub $ProgDir2/sub_calculate_neutrality_stats.sh
qsub $ProgDir2/sub_calculate_fst.sh
# qsub $ProgDir2/sub_calculate_haplotype_based_stats.sh
``` -->


## Calculate Four Gamete Test


## Make a slimmed vcf file of crown rot isolates

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="11-40 17-21 P295 62471 R36_14 371 SCRP376 SCRP370 10300 LV007"
Prefix=crown-rot_isolates
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```

##create directory structure

```bash
CurDir=/data/scratch/armita/idris
Prefix=crown-rot_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $WorkDir/all
mkdir -p $WorkDir/gff
mkdir -p $WorkDir/contigs
```

###copy the "gff" folder containing gff files

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix
Gff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.gff3)

cd $WorkDir/gff
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen
$ProgDir/summary_stats/split_gff_contig.sh $CurDir/$Gff
cd $CurDir
```

### Convert vcf files into fasta alignments

The vcf file needs to be converted into fasta alignment files. One file is
produced per reference contig

For the ploidy set <1|2|3>:
* 1 - haploid input
* 2 - diploid input, output as two separate phased haplotypes for each ind.
* 3 - diploid input, output as one sequence with ambiguity codes for each ind.

```bash
Prefix=crown-rot_isolates
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix

Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
Vcf=$(ls analysis/popgen/SNP_calling/summary_stats/crown-rot_isolates/crown-rot_isolates_filtered.recode.vcf)
Ploidy=2

cd $WorkDir/contigs
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen
python $ProgDir/summary_stats/vcf_to_fasta.py $CurDir/$Vcf $CurDir/$Reference $Ploidy
cd $CurDir
```


###This folder contains only contig FASTA files
###So create a new "contigs" directory to hold those files:

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix
Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
cp $Reference $WorkDir/.
```

###Next step: in the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
CurDir=/data/scratch/armita/idris
Prefix=crown-rot_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix
cd $WorkDir/contigs
for File in $(ls *.fasta); do
    Prefix=${File%.fasta}
    mkdir $Prefix
    mv $File $Prefix/.
done
cd $CurDir
```

### Test if all contigs have a matching gff and remove any which do not

```bash
CurDir=/data/scratch/armita/idris
Prefix=crown-rot_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix
cd $WorkDir
for File in $(ls $PWD/contigs/*/*.fasta); do
    Prefix=$(basename "$File" .fasta)
    expected_gff="$PWD/gff/${Prefix}.gff"
    if [ ! -f "$expected_gff" ]; then
       rm -rf $(dirname $File)
    fi
done
cd $CurDir
```

```bash
  cd /data/scratch/armita/idris/analysis/popgen/SNP_calling/summary_stats/crown-rot_isolates
```

```R

setwd('/data/scratch/armita/idris/analysis/popgen/SNP_calling/summary_stats/crown-rot_isolates')
library("PopGenome")
library(ggplot2)

Pcac <- c(
"414_1", "15_13_1", "PC13_15_1", "404_1", "4040_1", "15_7_1", "2003_3_1", "416_1", "415_1", "4032_1", "P421_1", "12420_1",
"414_2", "15_13_2", "PC13_15_2", "404_2", "4040_2", "15_7_2", "2003_3_2", "416_2", "415_2", "4032_2", "P421_2", "12420_2"
)
populations <- list(Pcac)
population_names <- c("Pcac")
population_no <- length(populations)

interval <-  10000
jump_size <-  interval / 10

gff <- "gff"
all_folders <- list.dirs("contigs", full.names = FALSE)
#Remove the gff folder from PopGenome contig analysis
contig_folders <- all_folders[all_folders != "gff"]

###Loop through each contig-containing folder to calculate stats on each contig separately.
for (dir in contig_folders[contig_folders != ""])
{
  contig_folder <- paste("contigs/", dir, sep="")
  GENOME.class <- readData(contig_folder, gffpath=gff, include.unknown = TRUE)
  GENOME.class <- set.populations(GENOME.class, populations)

  GENOME.class.split <- splitting.data(GENOME.class, subsites="gene")
  GENOME.class.slide <- sliding.window.transform(GENOME.class,width=interval,jump=jump_size,type=2, whole.data=TRUE)
  #per gene
  GENOME.class.split <- recomb.stats(GENOME.class.split)
  fourgamete_split <- get.recomb(GENOME.class.split)
  #per interval
  GENOME.class.slide <- recomb.stats(GENOME.class.slide)
  fourgamete_slide <- get.recomb(GENOME.class.slide)
  ids <- length(GENOME.class.slide@region.names)
  xaxis <- seq(from = 1, to = ids, by = 1)

  #Loop over each population: print figure and table with raw data to file
  for (i in seq_along(population_names))
  {
    fgt <- unlist(fourgamete_split[i])
    file_table = paste(dir, "_", population_names[i], "_4GT_per_gene.txt", sep="")
    file_table2 = paste("genome_", population_names[i], "_4GT_per_gene.txt", sep="")
    current_gff <- paste(gff, "/", dir, ".gff", sep="")
    gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr=dir, feature=FALSE, extract.gene.names=TRUE)
    fgt_table <- as.data.frame(fourgamete_split[i])
    write.table(fgt_table, file=file_table, sep="\t",quote=FALSE, col.names=FALSE)
    write.table(fgt_table, file=file_table2, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)
  }

  for (i in seq_along(population_names))
  {
    fgt <- unlist(fourgamete_slide[i])
    #write table with raw data
    slide_table <- paste(dir, "_", population_names[i], "_4GT_per_sliding_window.txt", sep="")
    slide_table2 <- paste("genome_", population_names[i], "_4GT_per_sliding_window.txt", sep="")
    write.table(as.data.frame(fourgamete_slide[i]), file=slide_table, sep="\t",quote=FALSE, col.names=FALSE)
    write.table(as.data.frame(fourgamete_slide[i]), file=slide_table2, sep="\t",quote=FALSE, col.names=FALSE, append=TRUE)
  }

}

##Print files with combined results across the entire genome
for (i in seq_along(population_names))
{
  #Four gamete test
  file_table2 = paste("genome_", population_names[i], "_4GT_per_gene.txt", sep="")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i], "_4GT_per_gene.pdf", sep="")
  fgt_plot <- ggplot(x, aes(x=x[,2])) + geom_histogram(colour="black", fill="cornsilk") + xlab("Four gamete test") + ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[,2], n = 10))
  ggsave(file_hist, fgt_plot)

  file_table2 = paste("genome_", population_names[i], "_4GT_per_sliding_window.txt", sep="")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i], "_4GT_per_sliding_window.pdf", sep="")
  fgt_plot <- ggplot(x, aes(x=x[,2])) + geom_histogram(colour="black", fill="cornsilk") + xlab("Four gamete test") + ylab("Number of intervals") + scale_x_continuous(breaks = pretty(x[,2], n = 10))
  ggsave(file_hist, fgt_plot)
}

```

```bash
cat genome_Pcac_4GT_per_gene.txt | grep 'NA' | wc -l
cat genome_Pcac_4GT_per_gene.txt | grep -v 'NA' | grep "0$" | wc -l
cat genome_Pcac_4GT_per_gene.txt | grep -v 'NA' | grep "1$" | wc -l
cat genome_Pcac_4GT_per_gene.txt | grep -v 'NA' | grep "2$" | wc -l
```

Of 1388 tested genes, none showed evidence of recombination as determined from the four gamete test.

```
27428
1388
0
0
```

### Linkage disequilibrium:

Commands for calculating LD decay are found in:
```
popgenome_scripts/LD_decay.md
```

Within crown rot isolates:
* 1524 SNPs were present in the dataset
* Low LD decay was observed
* The LD decay model could not be fitted to the data
* Low LD decay indicates a lack of recombination within crown rot isolates
