# Calculating Neutrality stats

This includes calculating the number of segregating sites, Watterson's Theta, Tajima's D, Fu & Li's F* and Fu & Li'd D

## Comparison between Pcac crown rot, apple and P. idaei populations

Slim down the vcf to the three populations:

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="11-40 17-21 10300 LV007"
Prefix=core_isolates
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```


## create the directory structure

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality

mkdir -p $WorkDir/all
mkdir -p $WorkDir/gff
mkdir -p $WorkDir/contigs
```

###copy the "gff" folder containing gff files

```bash
CurDir=/data/scratch/armita/idris
Prefix=core_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality
Gff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.gff3)

cd $WorkDir/gff
ProgDir=/home/adamst/git_repos/scripts/popgen
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
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality


Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
Vcf=$(ls analysis/popgen/SNP_calling/summary_stats/$Prefix/${Prefix}_filtered.recode.vcf)
Ploidy=2

cd $WorkDir/contigs
ProgDir=/home/adamst/git_repos/scripts/popgen
python $ProgDir/summary_stats/vcf_to_fasta.py $CurDir/$Vcf $CurDir/$Reference $Ploidy
cd $CurDir
```


### This folder contains only contig FASTA files
### So create a new "contigs" directory to hold those files:

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality

Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
cp $Reference $WorkDir/.
```

###Next step: in the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality

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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality

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
Prefix=core_isolates
cd analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality
```

```R
library("PopGenome")
library("ggplot2")

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: number of segregating sites, Watterson's Theta, Tajima's D, Fu & Li's F*, Fu & Li'd D*
Pcac_Fa <- c(
"414_1", "15_13_1", "PC13_15_1", "404_1", "4040_1", "15_7_1", "2003_3_1", "416_1", "415_1", "4032_1", "P421_1", "12420_1",
"414_2", "15_13_2", "PC13_15_2", "404_2", "4040_2", "15_7_2", "2003_3_2", "416_2", "415_2", "4032_2", "P421_2", "12420_2"
)
Pcac_Md <- c(
  "P295_1", "62471_1", "R36_14_1",
  "P295_2", "62471_2", "R36_14_2"
)
Pi_ri <- c(
"371_1", "SCRP376_1", "SCRP370_1",
"371_2", "SCRP376_2", "SCRP370_2"
)
#Need to set argument diploid=TRUE if using diploid genomes in the below command:
populations <- list(Pcac_Fa, Pcac_Md, Pi_ri)
#Number of populations assigned above.
population_no <- length(populations)
population_names <- c("Pcac_Fa", "Pcac_Md", "Pi_ri")
#Given in the same order, as above.
#Interval and jump size used in the sliding window analysis
interval <-  10000
jump_size <-  interval / 10
#########################################################################

#Folder containing FASTA alignments in the current dir
gff <- "gff"
all_folders <- list.dirs("contigs", full.names = FALSE)
#Remove the gff folder from PopGenome contig analysis
contig_folders <- all_folders[all_folders != "gff"]

###Loop through each contig containing folder to calculate stats on each contig separately.
for (dir in contig_folders[contig_folders != ""]){
  contig_folder <- paste("contigs/", dir, sep = "")
  GENOME.class <- readData(contig_folder, gffpath = gff, include.unknown = TRUE)
  GENOME.class <- set.populations(GENOME.class, populations)

###############################################################################

#Calculate neutrality stats over genes
GENOME.class.split <- splitting.data(GENOME.class, subsites = "gene")
GENOME.class.split <- neutrality.stats(GENOME.class.split)
TajimaD <- GENOME.class.split@Tajima.D
TajimaD_d <- as.data.frame(TajimaD)

for (i in seq_along(population_names)){
  file_hist <- paste(dir, "_", population_names[i], "_TajimaD_per_gene.pdf",
  sep = "")
  Tajima_D_len <- length(TajimaD[, i])
  Tajima_D_na <- sum(sapply(TajimaD[, i], is.na))
  if (Tajima_D_len > Tajima_D_na){
  tajima_plot <- ggplot(TajimaD_d, aes(x = TajimaD_d[, i])) +
  geom_histogram(colour = "black", fill = "pink") + ggtitle(dir) +
  xlab("Tajima's D") + ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(TajimaD_d[, i], n = 10))
  ggsave(file_hist, tajima_plot)
}
  file_table <- paste(dir, "_", population_names[i], "_TajimaD_per_gene.txt",
  sep = "")
  file_table2 <- paste("genome_", population_names[i],
  "_TajimaD_per_gene_all.txt", sep = "")
  current_gff <- paste(gff, "/", dir, ".gff", sep = "")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
      feature = FALSE, extract.gene.names = TRUE)
  tajima_table <- cbind(gene_ids, TajimaD[, i])
  write.table(tajima_table, file = file_table, sep = "\t", quote = FALSE,
  col.names = FALSE)
  write.table(tajima_table, file = file_table2, sep = "\t", quote = FALSE,
  col.names = FALSE, append = TRUE)
}

FuLiF <- GENOME.class.split@Fu.Li.F
FuLiF_d <- as.data.frame(FuLiF)

for (i in seq_along(population_names)){
  file_hist <- paste(dir, "_", population_names[i], "_FuLiF_per_gene.pdf",
  sep = "")
  FuLiF_len <- length(FuLiF[, i])
  FuLiF_na <- sum(sapply(FuLiF[, i], is.na))
  if (FuLiF_len > FuLiF_na){
  FuLiF_plot <- ggplot(FuLiF_d, aes(x = FuLiF_d[, i])) +
  geom_histogram(colour = "black", fill = "maroon") + ggtitle(dir) +
  xlab("Fu & Li's F") + ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(FuLiF_d[, i], n = 10))
  ggsave(file_hist, FuLiF_plot)
}
  file_table <- paste(dir, "_", population_names[i], "_FuLiF_per_gene.txt",
  sep = "")
  file_table2 <- paste("genome_", population_names[i],
  "_FuLiF_per_gene_all.txt", sep = "")
  current_gff <- paste(gff, "/", dir, ".gff", sep = "")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
      feature = FALSE, extract.gene.names = TRUE)
  FuLiF_table <- cbind(gene_ids, FuLiF[, i])
  write.table(FuLiF_table, file = file_table, sep = "\t", quote = FALSE,
  col.names = FALSE)
  write.table(FuLiF_table, file = file_table2, sep = "\t", quote = FALSE,
  col.names = FALSE, append = TRUE)
}

FuLiD <- GENOME.class.split@Fu.Li.D
FuLiD_d <- as.data.frame(FuLiD)

for (i in seq_along(population_names)){
  file_hist <- paste(dir, "_", population_names[i], "_FuLiD_per_gene.pdf",
  sep = "")
  FuLiD_len <- length(FuLiD[, i])
  FuLiD_na <- sum(sapply(FuLiD[, i], is.na))
  if (FuLiD_len > FuLiD_na){
  FuLiD_plot <- ggplot(FuLiD_d, aes(x = FuLiD_d[, i])) +
  geom_histogram(colour = "black", fill = "orchid") + ggtitle(dir) +
  xlab("Fu & Li's D") + ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(FuLiD_d[, i], n = 10))
  ggsave(file_hist, FuLiD_plot)
}
  file_table <- paste(dir, "_", population_names[i], "_FuLiD_per_gene.txt",
  sep = "")
  file_table2 <- paste("genome_", population_names[i],
  "_FuLiD_per_gene_all.txt", sep = "")
  current_gff <- paste(gff, "/", dir, ".gff", sep = "")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
      feature = FALSE, extract.gene.names = TRUE)
  FuLiD_table <- cbind(gene_ids, FuLiD[, i])
  write.table(FuLiD_table, file = file_table, sep = "\t", quote = FALSE,
  col.names = FALSE)
  write.table(FuLiD_table, file = file_table2, sep = "\t", quote = FALSE,
  col.names = FALSE, append = TRUE)
}

#Theta per site
ThetaW <- GENOME.class.split@theta_Watterson / GENOME.class.split@n.sites
ThetaW_d <- as.data.frame(ThetaW)

for (i in seq_along(population_names)){
  file_hist <- paste(dir, "_", population_names[i], "_ThetaW_per_gene.pdf",
  sep = "")
  ThetaW_len <- length(ThetaW[, i])
  ThetaW_na <- sum(sapply(ThetaW[, i], is.na))
  if (ThetaW_len > ThetaW_na){
  ThetaW_plot <- ggplot(ThetaW_d, aes(x = ThetaW_d[, i])) +
  geom_histogram(colour = "black", fill = "orange") + ggtitle(dir) +
  xlab(expression(paste("Average Watterson's ", theta, " per site"))) +
  ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(ThetaW_d[, i], n = 10))
  ggsave(file_hist, ThetaW_plot)
}
  file_table <- paste(dir, "_", population_names[i], "_ThetaW_per_gene.txt",
  sep = "")
  file_table2 <- paste("genome_", population_names[i],
  "_ThetaW_per_gene_all.txt", sep = "")
  current_gff <- paste(gff, "/", dir, ".gff", sep = "")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
      feature = FALSE, extract.gene.names = TRUE)
  ThetaW_table <- cbind(gene_ids, ThetaW[, i])
  write.table(ThetaW_table, file = file_table, sep = "\t", quote = FALSE,
  col.names = FALSE)
  write.table(ThetaW_table, file = file_table2, sep = "\t", quote = FALSE,
  col.names = FALSE, append = TRUE)
}

#Average rate of segregating sites
no_seg <- GENOME.class.split@n.segregating.sites / GENOME.class.split@n.sites
no_seg_d <- as.data.frame(no_seg)

for (i in seq_along(population_names)){
  file_hist <- paste(dir, "_", population_names[i],
  "_no_seg_sites_per_gene.pdf", sep = "")
  no_seg_len <- length(no_seg[, i])
  no_seg_na <- sum(sapply(no_seg[, i], is.na))
  if (no_seg_len > no_seg_na){
  no_seg_plot <- ggplot(no_seg_d, aes(x = no_seg_d[, i])) +
  geom_histogram(colour = "black", fill = "gold") + ggtitle(dir) +
  xlab("Average frequency of segregating sites") + ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(no_seg_d[, i], n = 10))
  ggsave(file_hist, no_seg_plot)
}
  file_table <- paste(dir, "_", population_names[i],
  "_no_seg_sites_per_gene.txt", sep = "")
  file_table2 <- paste("genome_", population_names[i],
  "_no_seg_sites_per_gene_all.txt", sep = "")
  current_gff <- paste(gff, "/", dir, ".gff", sep = "")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
      feature = FALSE, extract.gene.names = TRUE)
  no_seg_table <- cbind(gene_ids, no_seg[, i])
  write.table(no_seg_table, file = file_table, sep = "\t", quote = FALSE,
  col.names = FALSE)
  write.table(no_seg_table, file = file_table2, sep = "\t", quote = FALSE,
  col.names = FALSE, append = TRUE)
}

#Calculate neutrality stats over sliding window
GENOME.class.slide <- sliding.window.transform(GENOME.class, width = interval,
    jump = jump_size, type = 2, whole.data = TRUE)
GENOME.class.slide <- neutrality.stats(GENOME.class.slide)

#Watterson's theta
#plot the results for all populations over one figure
ThetaW <- GENOME.class.slide@theta_Watterson / interval
ThetaW_d <- as.data.frame(ThetaW)
#x axis
ids <- length(GENOME.class.slide@region.names)
xaxis <- seq(from = 1, to = ids, by = 1)

#Plot individual populations
for (i in seq_along(population_names)){
  file_slide <- paste(dir, "_", population_names[i],
  "_ThetaW_per_sliding_window.pdf", sep = "")
  ThetaW_len <- length(ThetaW[, i])
  ThetaW_na <- sum(sapply(ThetaW[, i], is.na))
  if (ThetaW_len > ThetaW_na){
  slide_plot <- ggplot(ThetaW_d, aes(x = xaxis, y = ThetaW_d[, i])) +
  geom_smooth(colour = "black", fill = "orange") + ggtitle(dir) +
  xlab("Contig coordinate (kbp)") +
  ylab(expression(paste("Average Watterson's ", theta, " per site"))) +
  scale_x_continuous(breaks = pretty(xaxis, n = 10))
  ggsave(file_slide, slide_plot)
}
  #write table with raw data
  slide_table <- paste(dir, "_", population_names[i],
  "_ThetaW_per_sliding_window.txt", sep = "")
  theta_table <- paste(GENOME.class.slide@region.names, ThetaW_d[, i])
  write.table(theta_table, file = slide_table, sep = "\t", quote = FALSE,
  col.names = FALSE)
}

#Plot both populations for comparison
if (ThetaW_len > ThetaW_na){
title <- paste(dir, "Comparison of", population_names[1], "ver.",
population_names[2], sep = " ")
comp_slide_file <- paste(dir, "_ThetaW_sliding_window_comparison.pdf", sep = "")
slide_comparison <- ggplot(ThetaW_d, aes(x = xaxis)) +
geom_smooth(aes(y = ThetaW_d[, 1]), colour = "red4") +
geom_smooth(aes(y = ThetaW_d[, 2]), colour = "darkseagreen2") +
ggtitle(title) + xlab("Contig coordinate (kbp)") +
ylab(expression(paste("Average Watterson's ", theta, " per site"))) +
scale_x_continuous(breaks = pretty(xaxis, n = 10))
ggsave(comp_slide_file, slide_comparison)
}


#Rate of segregating sites
no_seg <- GENOME.class.slide@n.segregating.sites / interval
no_seg_d <- as.data.frame(no_seg)
#x axis
ids <- length(GENOME.class.slide@region.names)
xaxis <- seq(from = 1, to = ids, by = 1)

#Plot individual populations
for (i in seq_along(population_names)){
  file_slide <- paste(dir, "_", population_names[i],
  "_no_seg_sites_per_sliding_window.pdf", sep = "")
  no_seg_len <- length(no_seg[, i])
  no_seg_na <- sum(sapply(no_seg[, i], is.na))
  if (no_seg_len > no_seg_na){
  slide_plot <- ggplot(no_seg_d, aes(x = xaxis, y = no_seg_d[, i])) +
  geom_smooth(colour = "black", fill = "gold") + ggtitle(dir) +
  xlab("Contig coordinate (kbp)") +
  ylab("Average frequency of segregating sites") +
  scale_x_continuous(breaks = pretty(xaxis, n = 10))
  ggsave(file_slide, slide_plot)
}
  #write table with raw data
  slide_table <- paste(dir, "_", population_names[i],
  "_no_seg_sites_per_sliding_window.txt", sep = "")
  na_table <- paste(GENOME.class.slide@region.names, no_seg_d[, i])
  write.table(na_table, file = slide_table, sep = "\t", quote = FALSE,
  col.names = FALSE)
}

#Plot both populations for comparison
if (no_seg_len > no_seg_na){
title <- paste(dir, "Comparison of", population_names[1], "ver.",
population_names[2], sep = " ")
comp_slide_file <- paste(dir, "_no_seg_sites_window_comparison.pdf", sep = "")
slide_comparison <- ggplot(no_seg_d, aes(x = xaxis)) +
geom_smooth(aes(y = no_seg_d[, 1]), colour = "turquoise1") +
geom_smooth(aes(y = no_seg_d[, 2]), colour = "olivedrab2") +
ggtitle(title) + xlab("Contig coordinate (kbp)") +
ylab("Average frequency of segregating sites") +
scale_x_continuous(breaks = pretty(xaxis, n = 10))
ggsave(comp_slide_file, slide_comparison)
}

}

###Plot genome wide histograms
for (i in seq_along(population_names)){
  #Tajima's D table
  file_table2 <- paste("genome_", population_names[i],
  "_TajimaD_per_gene_all.txt", sep = "")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_",
  population_names[i], "_TajimaD_per_gene_all.pdf", sep = "")
  tajima_plot <- ggplot(x, aes(x = x[, 3])) +
  geom_histogram(colour = "black", fill = "pink") + xlab("Tajima's D") +
  ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[, 3], n = 10))
  ggsave(file_hist, tajima_plot)

  #Fu & Li's F* table
  file_table2 <- paste("genome_", population_names[i],
  "_FuLiF_per_gene_all.txt", sep = "")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i], "_FuLiF_per_gene_all.pdf",
  sep = "")
  FuLiF_plot <- ggplot(x, aes(x = x[, 3])) +
  geom_histogram(colour = "black", fill = "maroon") + xlab("Fu & Li's F") +
  ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[, 3], n = 10))
  ggsave(file_hist, FuLiF_plot)

  #Fu & Li's D* table
  file_table2 <- paste("genome_", population_names[i],
  "_FuLiD_per_gene_all.txt", sep = "")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i], "_FuLiD_per_gene_all.pdf",
  sep = "")
  FuLiD_plot <- ggplot(x, aes(x = x[, 3])) +
  geom_histogram(colour = "black", fill = "orchid") + xlab("Fu & Li's D") +
  ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[, 3], n = 10))
  ggsave(file_hist, FuLiD_plot)

  #Watterson's Theta table
  file_table2 <- paste("genome_", population_names[i],
  "_ThetaW_per_gene_all.txt", sep = "")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i],
  "_ThetaW_per_gene_all.pdf", sep = "")
  ThetaW_plot <- ggplot(x, aes(x = x[, 3])) +
  geom_histogram(colour = "black", fill = "orange") +
  xlab(expression(paste("Average Watterson's ", theta, " per site"))) +
  ylab("Number of genes") + scale_x_continuous(breaks = pretty(x[, 3], n = 10))
  ggsave(file_hist, ThetaW_plot)

  #number of seg sites table
  file_table2 <- paste("genome_", population_names[i],
  "_no_seg_sites_per_gene_all.txt", sep = "")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_", population_names[i],
  "_no_seg_sites_per_gene_all.pdf", sep = "")
  no_seg_plot <- ggplot(x, aes(x = x[, 3])) +
  geom_histogram(colour = "black", fill = "gold") +
  xlab("Average frequency of segregating sites") + ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(x[, 3], n = 10))
  ggsave(file_hist, no_seg_plot)
}
```
