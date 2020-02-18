# Calculating fst stats

This includes calculating the number of segregating sites, Watterson's Theta, Tajima's D, Fu & Li's F* and Fu & Li'd D


## Comparison between Pcac crown rot and apple populations

Slim down the vcf to the two populations (11-40 included in strawberry lineage):

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="17-21 10300 LV007 371 SCRP370 SCRP376"
Prefix=apple_vs_strawberry
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 1.0 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```

```
After filtering, kept 16 out of 16 Individuals
Outputting VCF file...
After filtering, kept 26963 out of a possible 329672 Sites
Run Time = 7.00 seconds
```

## create the directory structure

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

mkdir -p $WorkDir/all
mkdir -p $WorkDir/gff
mkdir -p $WorkDir/contigs
```

###copy the "gff" folder containing gff files

```bash
CurDir=/data/scratch/armita/idris
Prefix=apple_vs_strawberry
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst
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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst


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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
cp $Reference $WorkDir/.
```

###Next step: in the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

cd $WorkDir/contigs
for File in $(ls *.fasta); do
    Prefix=${File%.fasta}
    mkdir $Prefix
    mv $File $Prefix/.
done
cd $CurDir
```

Open an R session

```bash
cd $WorkDir
~/prog/R/R-3.2.2/bin/R
```

```R
options(download.file.method = "wget")
install.packages("PopGenome")
install.packages("ggplot2")
library("PopGenome")
library("ggplot2")

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: Fst all/pairwise per gene or sliding window or per SNP
#More than one population needs to be defined, of course!
#Output files either per contig or the entire genome (prefix genome_)
Pc_Fxa <- c("415_1","15_7_1","15_13_1","404_1","12420_1","2003_3_1","4032_1","PC13_15_1","414_1","4040_1","416_1", "P421_1", "11-40_1", "415_2","15_7_2","15_13_2","404_2","12420_2","2003_3_2","4032_2","PC13_15_2","414_2","4040_2","416_2", "P421_2", "11-40_2")
Pc_Mxd <- c("62471_1","P295_1","R36_14_1", "62471_2","P295_2","R36_14_2")
#In the output for pairwise FST, pop1, pop2 etc. refer to the order in which the populations have been listed here:
populations <- list(Pc_Fxa, Pc_Mxd)
#Number of populations assigned above.
population_no <- length(populations)
pairs <- choose(population_no, 2)
population_names <- c("Pc_Fxa", "Pc_Mxd")
#Given in the same order, as above.
#Interval and jump size used in the sliding window analysis
# The sliding window size was set to 5000bp as there were 26963
# Snps in a 66Mb genome meaning ~ 1 snp per 2.5 Kb.
interval <-  10000
jump_size <-  interval / 10
#########################################################################

#Folder containing FASTA alignments in the current dir
gff <- "gff"
all_folders <- list.dirs("contigs", full.names = FALSE)
#Remove the gff folder from PopGenome contig analysis
contig_folders <- all_folders[all_folders != "gff"]

###Loop through each contig containing folder to calculate stats on each contig separately.
contig_df <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("contig", "fst", "total_sites", "biallelic_sites"))))
for (dir in contig_folders[contig_folders != ""]){
  contig_folder <- paste("contigs/", dir, sep = "")
  GENOME.class <- readData(contig_folder, gffpath = FALSE, include.unknown = TRUE)
  # get.individuals(GENOME.class)
  GENOME.class <- set.populations(GENOME.class, populations, diploid = TRUE)

  # contig based analysis
  GENOME.class <- F_ST.stats(GENOME.class)
  GENOME.class@nuc.F_ST.pairwise
  #           contig_1.fasta
  # pop1/pop2      0.9000798
  contig_stats <- (c(dir, round(GENOME.class@nuc.F_ST.pairwise[1], 3), as.integer(GENOME.class@n.sites), as.integer(GENOME.class@n.biallelic.sites)))
  contig_stats_df <- as.data.frame(t(contig_stats))
  colnames(contig_stats_df) <- c("contig", "fst", "total_sites", "biallelic_sites")
  contig_df <- rbind(contig_df, contig_stats_df)
}
contig_df$fst <- as.numeric(as.character(contig_df$fst))
contig_df$total_sites <- as.numeric(as.character(contig_df$total_sites))
contig_df$biallelic_sites <- as.numeric(as.character(contig_df$biallelic_sites))
write.table(contig_df, file = "total_Fst_per_contig.tsv", sep = "\t", quote = FALSE,
col.names = TRUE)
```

```
contig               fst           total_sites      biallelic_sites
Length:142         Min.   :-0.1200   Min.   :  31782   Min.   :  2.00  
Class :character   1st Qu.: 0.8505   1st Qu.: 163568   1st Qu.: 50.25  
Mode  :character   Median : 0.8905   Median : 365662   Median :155.50  
                Mean   : 0.8755   Mean   : 454327   Mean   :189.86  
                3rd Qu.: 0.9197   3rd Qu.: 639063   3rd Qu.:288.00  
                Max.   : 1.0000   Max.   :1749730   Max.   :822.00  
```

The output fst file was downloaded to my local machine

```bash
cd /Users/armita/Downloads/Pc
mkdir -p Pc_popgen/apple_vs_strawberry
cd Pc_popgen/apple_vs_strawberry
scp cluster:/data/scratch/armita/idris/analysis/popgen/SNP_calling/summary_stats/apple_vs_strawberry/fst/total_Fst_per_contig.tsv .
```

This file was opened in R studio for plotting

```R
library("ggplot2")
df1 <- read.delim("~/Downloads/Pc/Pc_popgen/apple_vs_strawberry/total_Fst_per_contig.tsv")
df1$fst[df1$fst < 0] <- 0
# df1$comparison <- as.factor("apple vs strawberry isolates")
df1$comparison <- as.factor("a vs b")
p <- ggplot(df1, aes(x=comparison, y=fst)) +
  geom_boxplot()
p <- p + ylab(expression(italic("F")["ST"]))
p <- p + theme(
  axis.text.x = element_text(size=12, angle=-45),
  axis.text.y = element_text(size=12),
  axis.title.x = element_text(size=FALSE),
  axis.title.y = element_text(size=14),
        )
p <- p + scale_y_continuous(minor_breaks = seq(0 , 1, 0.05), breaks = seq(0, 1, 0.25), limits=c(0, 1))
ggsave("total_Fst_per_contig.pdf", p, width = 3, height = 3)
```
<!--
## Gene based analyses:

### Test if all contigs have a matching gff and remove any which do not

```bash
CurDir=/data/scratch/armita/idris
Prefix=apple_vs_strawberry
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

cd $WorkDir
for File in $(ls $PWD/contigs/*/*.fasta); do
    FilePrefix=$(basename "$File" .fasta)
    expected_gff="$PWD/gff/${FilePrefix}.gff"
    if [ ! -f "$expected_gff" ]; then
       rm -rf $(dirname $File)
    fi
done
cd $CurDir
```

Open an R session

```bash
cd $WorkDir
~/prog/R/R-3.2.2/bin/R
```

Generating Fst stats by gene or by sliding window:
```R
Open an R session

```bash
cd $WorkDir
~/prog/R/R-3.2.2/bin/R
```

```R
options(download.file.method = "wget")
install.packages("PopGenome")
install.packages("ggplot2")
library("PopGenome")
library("ggplot2")

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: Fst all/pairwise per gene or sliding window or per SNP
#More than one population needs to be defined, of course!
#Output files either per contig or the entire genome (prefix genome_)
Pc_Fxa <- c("415_1","15_7_1","15_13_1","404_1","12420_1","2003_3_1","4032_1","PC13_15_1","414_1","4040_1","416_1", "P421_1", "11-40_1", "415_2","15_7_2","15_13_2","404_2","12420_2","2003_3_2","4032_2","PC13_15_2","414_2","4040_2","416_2", "P421_2", "11-40_2")
Pc_Mxd <- c("62471_1","P295_1","R36_14_1", "62471_2","P295_2","R36_14_2")
#In the output for pairwise FST, pop1, pop2 etc. refer to the order in which the populations have been listed here:
populations <- list(Pc_Fxa, Pc_Mxd)
#Number of populations assigned above.
population_no <- length(populations)
pairs <- choose(population_no, 2)
population_names <- c("Pc_Fxa", "Pc_Mxd")
#Given in the same order, as above.
#Interval and jump size used in the sliding window analysis
# The sliding window size was set to 5000bp as there were 26963
# Snps in a 66Mb genome meaning ~ 1 snp per 2.5 Kb.
interval <-  10000
jump_size <-  interval / 10
#########################################################################

#Folder containing FASTA alignments in the current dir
gff <- "gff"
all_folders <- list.dirs("contigs", full.names = FALSE)
#Remove the gff folder from PopGenome contig analysis
contig_folders <- all_folders[all_folders != "gff"]

###Loop through each contig containing folder to calculate stats on each contig separately.
for (dir in contig_folders[contig_folders != ""][1]){
  contig_folder <- paste("contigs/", dir, sep = "")
  GENOME.class <- readData(contig_folder, gffpath = gff, include.unknown = TRUE)
  # get.individuals(GENOME.class)
  GENOME.class <- set.populations(GENOME.class, populations, diploid = TRUE)
  GENOME.class.genes.concat <- splitting.data(GENOME.class, subsites = "gene", whole.data = TRUE)
  GENOME.class.genes.concat <- F_ST.stats(GENOME.class.genes.concat)
  biallelic_list <- unlist(GENOME.class.genes.concat@region.data@biallelic.sites)
  positions=paste(biallelic_list, collapse=', ')

  positions
  # GENOME.class.genes.concat2 <- splitting.data(GENOME.class,positions=list(16238, 25926, 30952, 31113, 34093, 35760, 40438, 42933, 44390, 45772, 51760, 129374, 137075, 147119, 147723, 151657, 161261, 163280, 164015, 164190, 166494, 173228, 173651, 173842, 174664, 176566, 176567, 177631, 179371, 181668, 182451, 185332, 186469, 188990, 189972, 190687, 190814, 191447, 193772, 199753, 235102, 241566, 243136, 249380, 250316, 257327, 258498, 258884, 258912, 259084, 263373, 267399, 275149, 277381, 278599, 280790, 280979, 296012, 301213, 311413, 311893, 313376, 313588, 314883, 315742, 315921, 324617, 330096, 331514, 345495, 349452, 349531, 349995, 351243, 352324, 352639, 354139, 354387, 354449, 354780, 356452, 357584, 362723, 371049, 375038, 377217, 378462, 379628, 382100, 382313, 382886, 387606, 388050, 388198, 389383, 390839, 390844, 390941, 393195, 396871, 399309, 400593, 405856, 406970, 409596, 412020, 412037, 412712, 414681, 418443, 421168, 421723, 422471, 425956, 432157, 433199, 433299, 434634, 435894, 437562, 438164, 440261, 440307, 440599, 441067, 442305, 445678, 447268, 449609, 450971, 451835, 461231, 468082, 469744, 473085, 473170, 473499, 474415, 475794, 478424, 480493, 482556, 483537, 486513, 509195, 515803, 518478, 520552, 521147, 529380, 529717, 529747, 531225, 532550, 564360, 573681, 576785, 577277, 585853, 585902, 586294, 586708, 595603, 596233, 604684, 606069, 607522, 610734, 617519, 620519, 622048, 622742, 623695, 626191, 626196, 627133, 627329, 628128, 629294, 629400, 633202, 633893, 636286, 638561, 639016, 644042, 649248, 651452, 655133, 663111, 663113, 682328, 683380, 686008, 687864, 691424, 692183, 693104, 694288, 694444, 694745, 700065, 707218, 708662, 709203, 715066, 720221, 722651, 722840, 722944, 724332, 727064, 728676, 730234, 733875, 734112, 736308, 738571, 738618, 739714, 741401, 741744, 742842, 743810, 744087, 746064, 748347, 750414, 752800, 753964, 767669, 767922, 778469, 781316, 787357, 788886, 789865, 794533, 798638, 799192, 804321, 806031, 814040, 851391, 856950, 857681, 857682, 861807, 865262, 867335, 870372, 871102, 871491, 873495, 873952, 873968, 874667, 879646, 879915, 884786, 884995, 885141, 887757, 895674, 897702, 898787, 903404, 904067, 904280, 905108, 911596, 928412, 929361, 931398, 935305, 948832, 951616, 953422, 954681, 954725, 956216, 957915, 959780, 975690, 976452, 979676, 980343, 980350, 988338, 991767, 993094, 993756, 995257, 1004249, 1008177, 1008678, 1009365, 1018494, 1021958, 1026593, 1034068, 1035971, 1043681, 1043729, 1045596, 1048827, 1052075, 1055574, 1065335, 1065660, 1066610, 1066625, 1067960, 1068086, 1071780, 1073638, 1074382, 1075500, 1076210, 1077826, 1079868, 1079936, 1091086, 1096500, 1098014, 1102363, 1102822, 1104646, 1105558, 1107115, 1108696, 1113714, 1147451, 1148739, 1151103, 1154267, 1157647, 1158451, 1166972, 1167404, 1170479, 1170585, 1173951, 1175151, 1178251, 1191366, 1196781, 1198319, 1200146, 1201475, 1204709, 1206794, 1207854, 1209654, 1210627, 1210903, 1212865, 1213345, 1213692, 1215366, 1220158, 1230604, 1232123, 1234432, 1241038, 1249003, 1249286, 1254491, 1261268, 1262641, 1266603, 1267331, 1267442, 1270904, 1272976, 1283243, 1284438, 1286212, 1299054, 1299269, 1299819, 1299882, 1303394, 1306612, 1307272, 1310066, 1311493, 1315259, 1319037, 1323105, 1323880, 1326227, 1334188, 1371083, 1371399, 1371446, 1371679, 1371836, 1374207, 1376405, 1377152, 1380260, 1380292, 1410662, 1411428, 1411775, 1412564, 1416656, 1416797, 1421330, 1452141, 1454167, 1456451, 1456952, 1458073, 1486604, 1486956, 1488284, 1491483, 1500750, 1510235, 1512758, 1514822, 1515914, 1516914, 1524061, 1524611, 1526651, 1531627, 1533230, 1533975, 1534648, 1536921, 1537589, 1537931, 1538951, 1541175, 1543604, 1556512, 1556722, 1556738, 1557772, 1566587, 1566599, 1596942, 1602634, 1604377, 1618931, 1620576, 1622274, 1624556, 1624648, 1624974, 1626761, 1629127, 1629397, 1629556, 1630758, 1631884, 1632400, 1637155, 1657519, 1692043, 1693032, 1698650, 1716868, 1716927, 1718697, 1720116, 1720943, 1721944, 1721968, 1725656),type=2, whole.data=TRUE)
  # WHOLE <- concatenate.regions(GENOME.class.genes.concat2)
  # GENOME.class.genes.concat2 <- concatenate.regions(GENOME.class.genes.concat2)
  # GENOME.class.genes.concat2 <- F_ST.stats(GENOME.class.genes.concat2)
  # GENOME.class.genes.concat2@nuc.F_ST.pairwise

###############################################################################

#### Gene based analysis
    GENOME.class.split <- splitting.data(GENOME.class, subsites = "gene")
    GENOME.class.split <- F_ST.stats(GENOME.class.split)
    # get.F_ST(GENOME.class.split)

    FST_all <- GENOME.class.split@nuc.F_ST.vs.all
    FST_all_d <- as.data.frame(FST_all)
    FST_pairwise <- GENOME.class.split@nuc.F_ST.pairwise
    # Hudson_KST <- GENOME.class.split@Hudson.K_ST

    for (i in seq_along(population_names)){
      # file_hist <- paste(dir, "_", population_names[i], "_total_FST_per_gene.pdf", sep = "")
      # fst_plot <- ggplot(FST_all_d, aes(x = FST_all_d[, i])) +
      # geom_histogram(colour = "black", fill = "darkseagreen") + ggtitle(dir) +
      # xlab(expression(paste("Total FST per gene"))) + ylab("Number of genes") +
      # scale_x_continuous(breaks = pretty(FST_all_d[, i], n = 10))
      # ggsave(file_hist, fst_plot)
      file_table <- paste(dir, "_", population_names[i], "_total_FST_per_gene.txt",
      sep = "")
      file_table2 <- paste("genome_", population_names[i],
      "_total_FST_per_gene_all.txt", sep = "")
      current_gff <- paste(gff, "/", dir, ".gff", sep = "")
      gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
          feature = FALSE, extract.gene.names = TRUE)
      fst_table <- cbind(gene_ids, FST_all[, i])
      write.table(fst_table, file = file_table, sep = "\t", quote = FALSE,
      col.names = FALSE)
      write.table(fst_table, file = file_table2, sep = "\t", quote = FALSE,
      col.names = FALSE, append = TRUE)
    }

    for (i in seq(pairs)){
      FST_pairwise_d <- as.data.frame(as.vector(FST_pairwise[i, ]))
      labelling <- gsub("/", "_vs_", row.names(FST_pairwise)[i])
      file_hist <- paste(dir, "_pairwise_FST_per_gene_", labelling, ".pdf",
      sep = "")
      fst_plot <- ggplot(FST_pairwise_d, aes(x = FST_pairwise_d[, 1])) +
      geom_histogram(colour = "black", fill = "cadetblue") + ggtitle(dir) +
      xlab(expression(paste("Pairwise FST per gene"))) + ylab("Number of genes") +
      scale_x_continuous(breaks = pretty(FST_pairwise_d[, 1], n = 10))
      ggsave(file_hist, fst_plot)
      file_table <- paste(dir, "_pairwise_FST_per_gene_", labelling, ".txt",
      sep = "")
      file_table2 <- paste("genome_pairwise_FST_per_gene_", labelling, ".txt",
      sep = "")
      current_gff <- paste(gff, "/", dir, ".gff", sep = "")
      gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
          feature = FALSE, extract.gene.names = TRUE)
      fst_table <- cbind(gene_ids, FST_pairwise_d[, 1])
      write.table(fst_table, file = file_table, sep = "\t", quote = FALSE,
      col.names = FALSE, row.names = FALSE)
      write.table(fst_table, file = file_table2, sep = "\t", quote = FALSE,
      col.names = FALSE, row.names = FALSE, append = TRUE)
    }

    ### Plot Hudson KST
    Hudson_KST_d <- as.data.frame(as.vector(Hudson_KST[1, ]))
    file_hist <- paste(dir, "_Hudson_KST_per_gene", ".pdf", sep = "")
    fst_plot <- ggplot(Hudson_KST_d, aes(x = Hudson_KST_d[, 1])) +
    geom_histogram(colour = "black", fill = "springgreen") + ggtitle(dir) +
    xlab(expression(paste("Hudson KST per gene"))) + ylab("Number of genes") +
    scale_x_continuous(breaks = pretty(Hudson_KST_d[, 1], n = 10))
    ggsave(file_hist, fst_plot)
    file_table <- paste(dir, "_Hudson_KST_per_gene", ".txt", sep = "")
    file_table2 <- paste("genome_Hudson_KST_per_gene_all", ".txt", sep = "")
    current_gff <- paste(gff, "/", dir, ".gff", sep = "")
    gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
        feature = FALSE, extract.gene.names = TRUE)
    fst_table <- cbind(gene_ids, Hudson_KST_d[, 1])
    write.table(fst_table, file = file_table, sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE)
    write.table(fst_table, file = file_table2, sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE, append = TRUE)

    #### Sliding window analysis (interval)
    GENOME.class.slide <- sliding.window.transform(GENOME.class, width = interval,
        jump = jump_size, type = 2, whole.data = TRUE)
    GENOME.class.slide <- F_ST.stats(GENOME.class.slide)
    FST_all_slide <- GENOME.class.slide@nuc.F_ST.vs.all
    FST_pairwise_slide <- GENOME.class.slide@nuc.F_ST.pairwise
    FST_all_slide_d <- as.data.frame(FST_all_slide)
    FST_pairwise_slide_d <- as.data.frame(t(FST_pairwise_slide))
    fst_table <- cbind(GENOME.class.slide@region.names, FST_pairwise_slide_d[, 1])
    fst_table_d <- data.frame(fst_table)
    #x axis
    ids <- length(GENOME.class.slide@region.names)
    xaxis <- seq(from = 1, to = ids, by = 1)
    # record number of fixed, shared & polymorphic sites per interval
    GENOME.class.slide <- calc.fixed.shared(GENOME.class.slide)
    fixation_df <- data.frame()
    fixation_df <- data.frame(GENOME.class.slide@n.monomorphic.sites)
    fixation_df$shared <- GENOME.class.slide@n.shared.sites
    fixation_df$fixed <- GENOME.class.slide@n.fixed.sites
    fixation_df <- cbind(fixation_df, FST_all_slide_d[1])
    fixation_df <- cbind(fixation_df, FST_all_slide_d[2])
    fixation_df <- cbind(fixation_df, fst_table_d[2])
    colnames(fixation_df) <- c("monomorphic", "shared", "fixed", "Fst_pop1", "Fst_pop2", "Fst_pop1_vs_pop2")
    file_fixation_df <- paste(dir, "pairwise_fixation_and_FST_per_sliding_window2.txt", sep = "_")
    write.table(fixation_df, file = file_fixation_df, sep = "\t", quote = FALSE, col.names = TRUE)


    #Plot individual populations
    for (i in seq_along(population_names)){
      file_slide <- paste(dir, "_", population_names[i],
      "_total_FST_per_sliding_window.pdf", sep = "")
    #  slide_plot <- ggplot(FST_all_slide_d, aes(x = xaxis,
    #      y = FST_all_slide_d[, i])) +
    #      geom_smooth(colour = "black", fill = "plum") + ggtitle(dir) +
    #      xlab("Contig coordinate (kbp)") + ylab("Total FST per interval") +
    #      scale_x_continuous(breaks = pretty(xaxis, n = 10))
    #  ggsave(file_slide, slide_plot)
      #write table with raw data
      slide_table <- paste(dir, "_", population_names[i],
      "_total_FST_per_sliding_window.txt", sep = "")
      write.table(FST_all_slide[, i], file = slide_table, sep = "\t",
      quote = FALSE, col.names = FALSE)
    }

    #Plot pairwise FST
    for (i in seq(pairs)){
      FST_pairwise_slide_d <- as.data.frame(as.vector(FST_pairwise_slide[i, ]))
      labelling <- gsub("/", "_vs_", row.names(FST_pairwise_slide)[i])
      file_hist <- paste(dir, "_pairwise_FST_per_sliding_window_", labelling,
      ".pdf", sep = "")
      # slide_plot <- ggplot(FST_pairwise_slide_d, aes(x = xaxis,
      #     y = FST_pairwise_slide_d[, 1])) +
      #     geom_smooth(colour = "black", fill = "slateblue") +
      #     ggtitle(dir) + xlab("Contig coordinate (kbp)") +
      #     ylab("Pairwise FST per interval") +
      #     scale_x_continuous(breaks = pretty(xaxis, n = 10))
      # ggsave(file_slide, slide_plot)
      #write table with raw data
      slide_table <- paste(dir, "_pairwise_FST_per_sliding_window_", labelling,
      ".txt", sep = "")
      fst_table <- cbind(GENOME.class.slide@region.names, FST_pairwise_slide_d[, 1])
      write.table(fst_table, file = slide_table, sep = "\t", quote = FALSE,
      col.names = FALSE, row.names = FALSE)
    }

}

##Genomewide results
for (i in seq_along(population_names)){
#Total FST
file_table2 <- paste("genome_", population_names[i],
"_total_FST_per_gene_all.txt", sep = "")
x <- as.data.frame(read.delim(file_table2))
file_hist <- paste("genome_", population_names[i],
"_total_FST_per_gene_all.pdf", sep = "")
fst_plot <- ggplot(x, aes(x = x[, 3])) +
geom_histogram(colour = "black", fill = "darkseagreen") +
xlab(expression(paste("Total FST per gene"))) + ylab("Number of genes") +
scale_x_continuous(breaks = pretty(x[, 3], n = 10))
ggsave(file_hist, fst_plot)
}

#Hudson KST
file_table2 <- paste("genome_Hudson_KST_per_gene_all", ".txt", sep = "")
x <- as.data.frame(read.delim(file_table2))
file_hist <- paste("genome_Hudson_KST_per_gene_all", ".pdf", sep = "")
fst_plot <- ggplot(x, aes(x = x[, 2])) +
geom_histogram(colour = "black", fill = "springgreen") +
xlab(expression(paste("Hudson KST per gene"))) + ylab("Number of genes") +
scale_x_continuous(breaks = pretty(x[, 2], n = 10))
ggsave(file_hist, fst_plot)


for (i in seq(pairs)){
  #Pairwise FST
  labelling <- gsub("/", "_vs_", row.names(FST_pairwise)[i])
  file_table2 <- paste("genome_pairwise_FST_per_gene_", labelling, ".txt",
  sep = "")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_pairwise_FST_per_gene_", labelling, ".pdf",
  sep = "")
  fst_plot <- ggplot(x, aes(x = x[, 2])) +
  geom_histogram(colour = "black", fill = "cadetblue") +
  xlab(expression(paste("Pairwise FST per gene"))) + ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(x[, 2], n = 10))
  ggsave(file_hist, fst_plot)
}

```
 -->


## Comparison within all Pcac populations

Slim down the vcf to the two populations (11-40 included in strawberry lineage):

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="371 SCRP370 SCRP376"
Prefix=within_cactorum
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 1.0 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```

```
After filtering, kept 17 out of 17 Individuals
Outputting VCF file...
After filtering, kept 35428 out of a possible 329672 Sites
Run Time = 9.00 seconds
```

## create the directory structure

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

mkdir -p $WorkDir/all
mkdir -p $WorkDir/gff
mkdir -p $WorkDir/contigs
```

###copy the "gff" folder containing gff files

```bash
CurDir=/data/scratch/armita/idris
Prefix=within_cactorum
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst
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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst


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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
cp $Reference $WorkDir/.
```

###Next step: in the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

cd $WorkDir/contigs
for File in $(ls *.fasta); do
  Prefix=${File%.fasta}
  mkdir $Prefix
  mv $File $Prefix/.
done
cd $CurDir
```

Open an R session

```bash
cd $WorkDir
~/prog/R/R-3.2.2/bin/R
```

```R
options(download.file.method = "wget")
install.packages("PopGenome")
install.packages("ggplot2")
library("PopGenome")
library("ggplot2")

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: Fst all/pairwise per gene or sliding window or per SNP
#More than one population needs to be defined, of course!
#Output files either per contig or the entire genome (prefix genome_)
Pc_Fxa <- c("415_1","15_7_1","15_13_1","404_1","12420_1","2003_3_1","4032_1","PC13_15_1","414_1","4040_1","416_1", "P421_1", "11-40_1", "415_2","15_7_2","15_13_2","404_2","12420_2","2003_3_2","4032_2","PC13_15_2","414_2","4040_2","416_2", "P421_2", "11-40_2")
Pc_Mxd <- c("62471_1","P295_1","R36_14_1", "62471_2","P295_2","R36_14_2")
Pc_LR <- c("17-21_1", "17-21_2")
#In the output for pairwise FST, pop1, pop2 etc. refer to the order in which the populations have been listed here:
populations <- list(Pc_Fxa, Pc_Mxd, Pc_LR)
#Number of populations assigned above.
population_no <- length(populations)
pairs <- choose(population_no, 2)
population_names <- c("Pc_Fxa", "Pc_Mxd", "Pc_LR")
#Given in the same order, as above.
#Interval and jump size used in the sliding window analysis
# The sliding window size was set to 5000bp as there were 26963
# Snps in a 66Mb genome meaning ~ 1 snp per 2.5 Kb.
interval <-  10000
jump_size <-  interval / 10
#########################################################################

#Folder containing FASTA alignments in the current dir
gff <- "gff"
all_folders <- list.dirs("contigs", full.names = FALSE)
#Remove the gff folder from PopGenome contig analysis
contig_folders <- all_folders[all_folders != "gff"]

###Loop through each contig containing folder to calculate stats on each contig separately.
contig_df <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("contig", "fst_a_vs_b", "fst_a_vs_c", "fst_b_vs_c", "total_sites", "biallelic_sites"))))
for (dir in contig_folders[contig_folders != ""]){
contig_folder <- paste("contigs/", dir, sep = "")
GENOME.class <- readData(contig_folder, gffpath = FALSE, include.unknown = TRUE)
# get.individuals(GENOME.class)
GENOME.class <- set.populations(GENOME.class, populations, diploid = TRUE)

# contig based analysis
GENOME.class <- F_ST.stats(GENOME.class)
GENOME.class@nuc.F_ST.pairwise
#           contig_1.fasta
# pop1/pop2      0.9000798
contig_stats <- (c(dir, round(GENOME.class@nuc.F_ST.pairwise[1], 3), round(GENOME.class@nuc.F_ST.pairwise[2], 3), round(GENOME.class@nuc.F_ST.pairwise[3], 3), as.integer(GENOME.class@n.sites), as.integer(GENOME.class@n.biallelic.sites)))
contig_stats_df <- as.data.frame(t(contig_stats))
colnames(contig_stats_df) <- c("contig", "fst_a_vs_b", "fst_a_vs_c", "fst_b_vs_c", "total_sites", "biallelic_sites")
contig_df <- rbind(contig_df, contig_stats_df)
}
contig_df$fst_a_vs_b <- as.numeric(as.character(contig_df$fst_a_vs_b))
contig_df$fst_a_vs_c <- as.numeric(as.character(contig_df$fst_a_vs_c))
contig_df$fst_b_vs_c <- as.numeric(as.character(contig_df$fst_b_vs_c))
contig_df$total_sites <- as.numeric(as.character(contig_df$total_sites))
contig_df$biallelic_sites <- as.numeric(as.character(contig_df$biallelic_sites))
write.table(contig_df, file = "total_Fst_per_contig_within_cactorum.tsv", sep = "\t", quote = FALSE,
col.names = TRUE)
```

```
contig      fst_a_vs_b        fst_a_vs_c        fst_b_vs_c     
contig_1  :  1   Min.   :-0.1200   Min.   :-0.5200   Min.   :-0.6000  
contig_10 :  1   1st Qu.: 0.8505   1st Qu.: 0.9800   1st Qu.: 0.8380  
contig_100:  1   Median : 0.8905   Median : 0.9880   Median : 0.8850  
contig_101:  1   Mean   : 0.8755   Mean   : 0.9514   Mean   : 0.8479  
contig_102:  1   3rd Qu.: 0.9197   3rd Qu.: 0.9940   3rd Qu.: 0.9175  
contig_103:  1   Max.   : 1.0000   Max.   : 1.0000   Max.   : 1.0000  
(Other)   :137   NA's   :1                                            
total_sites      biallelic_sites
Min.   :  31782   Min.   :   2.0  
1st Qu.: 162975   1st Qu.:  63.5  
Median : 363204   Median : 203.0  
Mean   : 451514   Mean   : 247.7  
3rd Qu.: 638644   3rd Qu.: 365.5  
Max.   :1749730   Max.   :1115.0
```

The output fst file was downloaded to my local machine

```bash
cd /Users/armita/Downloads/Pc
mkdir -p Pc_popgen/within_cactorum
cd Pc_popgen/within_cactorum
scp cluster:/data/scratch/armita/idris/analysis/popgen/SNP_calling/summary_stats/within_cactorum/fst/total_Fst_per_contig_within_cactorum.tsv .
```

This file was opened in R studio for plotting

```R
library("ggplot2")
library("reshape")
df4 <- read.delim("~/Downloads/Pc/Pc_popgen/within_cactorum/total_Fst_per_contig_within_cactorum.tsv")
# colnames(df2) <- c("contig", "fst_i_vs_ii", "fst_i_vs_iii", "fst_ii_vs_iii", "biallelic_sites", "total_sites")
colnames(df4) <- c("contig", "a vs b", "a vs c", "b vs c", "biallelic_sites", "total_sites")
df4$a_vs_b[df4$a_vs_b < 0] <- 0
df4$a_vs_c[df4$a_vs_c < 0] <- 0
df4$b_vs_c[df4$b_vs_c < 0] <- 0
df5 <- melt(df4, id=c("contig","biallelic_sites", "total_sites"))
colnames(df5) <- c("contig", "biallelic_sites", "total_sites", "comparison", "fst")
p <- ggplot(df5, aes(x=comparison, y=fst)) +
geom_boxplot()
p <- p + ylab(expression(italic("F")["ST"]))
p <- p + theme(
axis.text.x = element_text(size=12, angle=-45),
axis.text.y = element_text(size=12),
axis.title.x = element_text(size=FALSE),
axis.title.y = element_text(size=14),
      )
p <- p + scale_y_continuous(minor_breaks = seq(0 , 1, 0.05), breaks = seq(0, 1, 0.25), limits=c(0, 1))
ggsave("total_Fst_per_contig_within_cactorum.pdf", p, width = 5, height = 3)
```


## Comparison within Pcac crown rot population

Slim down the vcf to the two populations (11-40 included in strawberry lineage):

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="17-21 10300 LV007 371 SCRP370 SCRP376 62471 P295 R36_14"
Prefix=within_strawberry
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 1.0 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```

```
After filtering, kept 13 out of 13 Individuals
Outputting VCF file...
After filtering, kept 1535 out of a possible 329672 Sites
Run Time = 6.00 seconds
```

## create the directory structure

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

mkdir -p $WorkDir/all
mkdir -p $WorkDir/gff
mkdir -p $WorkDir/contigs
```

###copy the "gff" folder containing gff files

```bash
CurDir=/data/scratch/armita/idris
Prefix=within_strawberry
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst
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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst


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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
cp $Reference $WorkDir/.
```

###Next step: in the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

cd $WorkDir/contigs
for File in $(ls *.fasta); do
   Prefix=${File%.fasta}
   mkdir $Prefix
   mv $File $Prefix/.
done
cd $CurDir
```

Open an R session

```bash
cd $WorkDir
~/prog/R/R-3.2.2/bin/R
```

```R
options(download.file.method = "wget")
install.packages("PopGenome")
install.packages("ggplot2")
library("PopGenome")
library("ggplot2")

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: Fst all/pairwise per gene or sliding window or per SNP
#More than one population needs to be defined, of course!
#Output files either per contig or the entire genome (prefix genome_)
Pc_Fxa_i <- c("415_1","15_7_1","15_13_1","404_1","12420_1","11-40_1", "415_2","15_7_2","15_13_2","404_2","12420_2","11-40_2")
Pc_Fxa_ii <- c("4040_1","416_1","4040_2","416_2")
Pc_Fxa_iii <- c("2003_3_1","4032_1","PC13_15_1","414_1","P421_1", "2003_3_2","4032_2","PC13_15_2","414_2","P421_2")

#In the output for pairwise FST, pop1, pop2 etc. refer to the order in which the populations have been listed here:
populations <- list(Pc_Fxa_i, Pc_Fxa_ii, Pc_Fxa_iii)
#Number of populations assigned above.
population_no <- length(populations)
pairs <- choose(population_no, 2)
population_names <- c("Pc_Fxa_i", "Pc_Fxa_ii", "Pc_Fxa_iii")
#Given in the same order, as above.
#Interval and jump size used in the sliding window analysis
# The sliding window size was set to 5000bp as there were 26963
# Snps in a 66Mb genome meaning ~ 1 snp per 2.5 Kb.
interval <-  10000
jump_size <-  interval / 10
#########################################################################

#Folder containing FASTA alignments in the current dir
gff <- "gff"
all_folders <- list.dirs("contigs", full.names = FALSE)
#Remove the gff folder from PopGenome contig analysis
contig_folders <- all_folders[all_folders != "gff"]

###Loop through each contig containing folder to calculate stats on each contig separately.
contig_df <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("contig", "fst_i_vs_ii", "fst_i_vs_iii", "fst_ii_vs_iii", "total_sites", "biallelic_sites"))))
for (dir in contig_folders[contig_folders != ""]){
 contig_folder <- paste("contigs/", dir, sep = "")
 GENOME.class <- readData(contig_folder, gffpath = FALSE, include.unknown = TRUE)
 # get.individuals(GENOME.class)
 GENOME.class <- set.populations(GENOME.class, populations, diploid = TRUE)

 # contig based analysis
 GENOME.class <- F_ST.stats(GENOME.class)
 GENOME.class@nuc.F_ST.pairwise
 #           contig_1.fasta
 # pop1/pop2      0.9000798
 contig_stats <- (c(dir, round(GENOME.class@nuc.F_ST.pairwise[1], 3), round(GENOME.class@nuc.F_ST.pairwise[2], 3), round(GENOME.class@nuc.F_ST.pairwise[3], 3), as.integer(GENOME.class@n.sites), as.integer(GENOME.class@n.biallelic.sites)))
 contig_stats_df <- as.data.frame(t(contig_stats))
 colnames(contig_stats_df) <- c("contig", "fst_i_vs_ii", "fst_i_vs_iii", "fst_ii_vs_iii", "total_sites", "biallelic_sites")
 contig_df <- rbind(contig_df, contig_stats_df)
}
contig_df$fst_i_vs_ii <- as.numeric(as.character(contig_df$fst_i_vs_ii))
contig_df$fst_i_vs_iii <- as.numeric(as.character(contig_df$fst_i_vs_iii))
contig_df$fst_ii_vs_iii <- as.numeric(as.character(contig_df$fst_ii_vs_iii))
contig_df$total_sites <- as.numeric(as.character(contig_df$total_sites))
contig_df$biallelic_sites <- as.numeric(as.character(contig_df$biallelic_sites))
write.table(contig_df, file = "total_Fst_per_contig_within_Fxa.tsv", sep = "\t", quote = FALSE,
col.names = TRUE)
```

```
contig               fst           total_sites      biallelic_sites
Length:142         Min.   :-0.1200   Min.   :  31782   Min.   :  2.00  
Class :character   1st Qu.: 0.8505   1st Qu.: 163568   1st Qu.: 50.25  
Mode  :character   Median : 0.8905   Median : 365662   Median :155.50  
               Mean   : 0.8755   Mean   : 454327   Mean   :189.86  
               3rd Qu.: 0.9197   3rd Qu.: 639063   3rd Qu.:288.00  
               Max.   : 1.0000   Max.   :1749730   Max.   :822.00  
```

The output fst file was downloaded to my local machine

```bash
cd /Users/armita/Downloads/Pc
mkdir -p Pc_popgen/within_strawberry
cd Pc_popgen/within_strawberry
scp cluster:/data/scratch/armita/idris/analysis/popgen/SNP_calling/summary_stats/within_strawberry/fst/total_Fst_per_contig_within_Fxa.tsv .
```

This file was opened in R studio for plotting

```R
library("ggplot2")
library("reshape")
df2 <- read.delim("~/Downloads/Pc/Pc_popgen/within_strawberry/total_Fst_per_contig_within_Fxa.tsv")
# colnames(df2) <- c("contig", "fst_i_vs_ii", "fst_i_vs_iii", "fst_ii_vs_iii", "biallelic_sites", "total_sites")
colnames(df2) <- c("contig", "i vs ii", "i vs iii", "ii vs iii", "biallelic_sites", "total_sites")
df2$i_vs_ii[df2$i_vs_ii < 0] <- 0
df2$i_vs_iii[df2$i_vs_iii < 0] <- 0
df2$ii_vs_iii[df2$ii_vs_iii < 0] <- 0
df3 <- melt(df2, id=c("contig","biallelic_sites", "total_sites"))
colnames(df3) <- c("contig", "biallelic_sites", "total_sites", "comparison", "fst")
# df2$comparison <- as.factor("apple vs strawberry isolates")
# df2$comparison <- as.factor("Pop. a vs b")
p <- ggplot(df3, aes(x=comparison, y=fst)) +
 geom_boxplot()
p <- p + ylab(expression(italic("F")["ST"]))
p <- p + theme(
 axis.text.x = element_text(size=12, angle=-45),
 axis.text.y = element_text(size=12),
 axis.title.x = element_text(size=FALSE),
 axis.title.y = element_text(size=14),
       )
p <- p + scale_y_continuous(minor_breaks = seq(0 , 1, 0.05), breaks = seq(0, 1, 0.25), limits=c(0, 1))
ggsave("total_Fst_per_contig_within_Fxa.pdf", p, width = 5, height = 3)
```

```
contig       i vs ii          i vs iii        ii vs iii     
contig_1  :  1   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
contig_10 :  1   1st Qu.:0.9600   1st Qu.:0.5990   1st Qu.:0.9475  
contig_100:  1   Median :0.9900   Median :0.9390   Median :0.9820  
contig_101:  1   Mean   :0.9001   Mean   :0.7244   Mean   :0.8997  
contig_102:  1   3rd Qu.:1.0000   3rd Qu.:0.9820   3rd Qu.:1.0000  
contig_103:  1   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
(Other)   :114   NA's   :13       NA's   :11       NA's   :6       
biallelic_sites    total_sites   
Min.   :  39583   Min.   : 1.00  
1st Qu.: 214998   1st Qu.: 4.00  
Median : 433546   Median : 9.00  
Mean   : 509578   Mean   :12.79  
3rd Qu.: 673452   3rd Qu.:19.00  
Max.   :1749730   Max.   :72.00  
```

Statistics calculated over the whole genome:

A fasta file was created of all the SNP sites across the genome:

```bash
# FastaList=$(ls contigs/*/*.fasta)
FastaList=$(ls contigs/*/*.fasta | sort -n -k3 -t'_' | head -n10)
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
$ProgDir/concat_fasta_alignments.py --inp_fasta $FastaList --out_fasta all/concat_contigs.fasta
```

```bash
CurDir=/data/scratch/armita/idris
Prefix=within_cactorum
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

Vcf=$(ls $WorkDir/../within_cactorum_filtered.recode.vcf)
Ploidy=2

ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen
$ProgDir/summary_stats/vcf_variants_to_fasta.py --inp_vcf $Vcf --out_fasta $WorkDir/all/concat_variants.fa
```


<!--
```bash
bgzip all/within_cactorum_filtered.recode.vcf
tabix -p vcf all/within_cactorum_filtered.recode.vcf.gz
``` -->

Open an R session

```bash
WorkDir=/data/scratch/armita/idris/analysis/popgen/SNP_calling/summary_stats/within_cactorum/fst
cd $WorkDir
~/prog/R/R-3.2.2/bin/R
```

```R
options(download.file.method = "wget")
install.packages("PopGenome")
install.packages("ggplot2")
library("PopGenome")
library("ggplot2")

# GENOME.class <- readVCF("all/within_cactorum_filtered.recode.vcf.gz", from=1, to=10000000,gffpath=FALSE, tid="contig_1", numcols = 1000000, approx = FALSE)

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: Fst all/pairwise per gene or sliding window or per SNP
#More than one population needs to be defined, of course!
#Output files either per contig or the entire genome (prefix genome_)
# Pc_Fxa <- c("415_1","15_7_1","15_13_1","404_1","12420_1","2003_3_1","4032_1","PC13_15_1","414_1","4040_1","416_1", "P421_1", "11-40_1", "415_2","15_7_2","15_13_2","404_2","12420_2","2003_3_2","4032_2","PC13_15_2","414_2","4040_2","416_2", "P421_2", "11-40_2")
# Pc_Mxd <- c("62471_1","P295_1","R36_14_1", "62471_2","P295_2","R36_14_2")
# Pc_LR <- c("17-21_1", "17-21_2")
# Pc_Fxa <- c("415.1","15_7.1","15_13.1","404.1","12420.1","2003_3.1","4032.1","PC13_15.1","414.1","4040.1","416.1", "P421.1", "11-40.1", "415.2","15_7.2","15_13.2","404.2","12420.2","2003_3.2","4032.2","PC13_15.2","414.2","4040.2","416.2", "P421.2", "11-40.2")
# Pc_Mxd <- c("62471.1","P295.1","R36_14.1", "62471.2","P295.2","R36_14.2")
# Pc_LR <- c("17-21.1", "17-21.2")
Pc_Fxa <- c("415","15_7","15_13","404","12420","2003_3","4032","PC13_15","414","4040","416", "P421", "11-40")
Pc_Mxd <- c("62471","P295","R36_14")
Pc_LR <- c("17-21")
#In the output for pairwise FST, pop1, pop2 etc. refer to the order in which the populations have been listed here:
populations <- list(Pc_Fxa, Pc_Mxd, Pc_LR)
#Number of populations assigned above.
population_no <- length(populations)
pairs <- choose(population_no, 2)
population_names <- c("Pc_Fxa", "Pc_Mxd", "Pc_LR")
#####################################

###Loop through each contig containing folder to calculate stats on each contig separately.
# contig_df <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("contig", "fst_a_vs_b", "fst_a_vs_c", "fst_b_vs_c", "total_sites", "biallelic_sites"))))

GENOME.class <- readData("all", gffpath = FALSE, include.unknown = TRUE)
# get.individuals(GENOME.class)
GENOME.class <- set.populations(GENOME.class, populations, diploid = TRUE)

# contig based analysis
GENOME.class <- F_ST.stats(GENOME.class)
GENOME.class@nuc.F_ST.pairwise
# concat_variants.fa
# pop1/pop2          0.8850618
# pop1/pop3          0.9830222
# pop2/pop3          0.8783476
get.F_ST(GENOME.class)
# haplotype.F_ST nucleotide.F_ST  Nei.G_ST Hudson.G_ST
# concat_variants.fa              0       0.9153229 0.1331719  0.04665282
# Hudson.H_ST Hudson.K_ST
# concat_variants.fa           0   0.9040782
GENOME.class <- diversity.stats(GENOME.class)
GENOME.class@nuc.diversity.within
#                       pop 1 pop 2 pop 3
# concat_variants.fa 530.1631  4494   156
# diversities need to be divided by the number of sites to give diversity per site.
GENOME.class@Pi
#                    pop 1 pop 2 pop 3
# concat_variants.fa     0     0     0
GENOME.class@hap.diversity.within
#                    pop 1 pop 2 pop 3
# concat_variants.fa     1     1     1
GENOME.class@hap.diversity.between
#           concat_variants.fa
# pop1/pop2                  1
# pop1/pop3                  1
# pop2/pop3                  1
get.sum.data(GENOME.class)
#                    n.sites n.biallelic.sites n.gaps n.unknowns n.valid.sites
# concat_variants.fa   35416             35416      0          0         35416
#                    n.polyallelic.sites trans.transv.ratio
# concat_variants.fa                   0           2.317969
```


<!--
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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

mkdir -p $WorkDir/all
mkdir -p $WorkDir/gff
mkdir -p $WorkDir/contigs
```

###copy the "gff" folder containing gff files

```bash
CurDir=/data/scratch/armita/idris
Prefix=core_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst
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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst


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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
cp $Reference $WorkDir/.
```

###Next step: in the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

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
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/fst

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


```R

library("PopGenome")
library("ggplot2")

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: Fst all/pairwise per gene or sliding window or per SNP
#More than one population needs to be defined, of course!
#Output files either per contig or the entire genome (prefix genome_)
Pc_Fxa <- c("415","15_7","15_13","404","12420","2003_3","4032","PC13_15","414","4040","416")
Pc_Mxd <- c("62471","P295","R36_14")
Pi_Fxa <- c("371","SCRP370","SCRP376")
#In the output for pairwise FST, pop1, pop2 etc. refer to the order in which the populations have been listed here:
populations <- list(Pc_Fxa, Pc_Mxd, Pi_Fxa)
#Number of populations assigned above.
population_no <- length(populations)
pairs <- choose(population_no, 2)
population_names <- c("Pc_Fxa", "Pc_Mxd", "Pi_Fxa")
#Given in the same order, as above.
#Interval and jump size used in the sliding window analysis
interval <-  1000
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

#### Gene based analysis
GENOME.class.split <- splitting.data(GENOME.class, subsites = "gene")
GENOME.class.split <- F_ST.stats(GENOME.class.split)
get.F_ST(GENOME.class.split)

FST_all <- GENOME.class.split@nuc.F_ST.vs.all
FST_all_d <- as.data.frame(FST_all)
FST_pairwise <- GENOME.class.split@nuc.F_ST.pairwise
Hudson_KST <- GENOME.class.split@Hudson.K_ST

for (i in seq_along(population_names)){
  file_hist <- paste(dir, "_", population_names[i], "_total_FST_per_gene.pdf",
  sep = "")
  fst_plot <- ggplot(FST_all_d, aes(x = FST_all_d[, i])) +
  geom_histogram(colour = "black", fill = "darkseagreen") + ggtitle(dir) +
  xlab(expression(paste("Total FST per gene"))) + ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(FST_all_d[, i], n = 10))
  ggsave(file_hist, fst_plot)
  file_table <- paste(dir, "_", population_names[i], "_total_FST_per_gene.txt",
  sep = "")
  file_table2 <- paste("genome_", population_names[i],
  "_total_FST_per_gene_all.txt", sep = "")
  current_gff <- paste(gff, "/", dir, ".gff", sep = "")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
      feature = FALSE, extract.gene.names = TRUE)
  fst_table <- cbind(gene_ids, FST_all[, i])
  write.table(fst_table, file = file_table, sep = "\t", quote = FALSE,
  col.names = FALSE)
  write.table(fst_table, file = file_table2, sep = "\t", quote = FALSE,
  col.names = FALSE, append = TRUE)
}

for (i in seq(pairs)){
  FST_pairwise_d <- as.data.frame(as.vector(FST_pairwise[i, ]))
  labelling <- gsub("/", "_vs_", row.names(FST_pairwise)[i])
  file_hist <- paste(dir, "_pairwise_FST_per_gene_", labelling, ".pdf",
  sep = "")
  fst_plot <- ggplot(FST_pairwise_d, aes(x = FST_pairwise_d[, 1])) +
  geom_histogram(colour = "black", fill = "cadetblue") + ggtitle(dir) +
  xlab(expression(paste("Pairwise FST per gene"))) + ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(FST_pairwise_d[, 1], n = 10))
  ggsave(file_hist, fst_plot)
  file_table <- paste(dir, "_pairwise_FST_per_gene_", labelling, ".txt",
  sep = "")
  file_table2 <- paste("genome_pairwise_FST_per_gene_", labelling, ".txt",
  sep = "")
  current_gff <- paste(gff, "/", dir, ".gff", sep = "")
  gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
      feature = FALSE, extract.gene.names = TRUE)
  fst_table <- cbind(gene_ids, FST_pairwise_d[, 1])
  write.table(fst_table, file = file_table, sep = "\t", quote = FALSE,
  col.names = FALSE, row.names = FALSE)
  write.table(fst_table, file = file_table2, sep = "\t", quote = FALSE,
  col.names = FALSE, row.names = FALSE, append = TRUE)
}

### Plot Hudson KST
Hudson_KST_d <- as.data.frame(as.vector(Hudson_KST[1, ]))
file_hist <- paste(dir, "_Hudson_KST_per_gene", ".pdf", sep = "")
fst_plot <- ggplot(Hudson_KST_d, aes(x = Hudson_KST_d[, 1])) +
geom_histogram(colour = "black", fill = "springgreen") + ggtitle(dir) +
xlab(expression(paste("Hudson KST per gene"))) + ylab("Number of genes") +
scale_x_continuous(breaks = pretty(Hudson_KST_d[, 1], n = 10))
ggsave(file_hist, fst_plot)
file_table <- paste(dir, "_Hudson_KST_per_gene", ".txt", sep = "")
file_table2 <- paste("genome_Hudson_KST_per_gene_all", ".txt", sep = "")
current_gff <- paste(gff, "/", dir, ".gff", sep = "")
gene_ids <- get_gff_info(GENOME.class.split, current_gff, chr = dir,
    feature = FALSE, extract.gene.names = TRUE)
fst_table <- cbind(gene_ids, Hudson_KST_d[, 1])
write.table(fst_table, file = file_table, sep = "\t", quote = FALSE,
col.names = FALSE, row.names = FALSE)
write.table(fst_table, file = file_table2, sep = "\t", quote = FALSE,
col.names = FALSE, row.names = FALSE, append = TRUE)

#### Sliding window analysis (interval)
GENOME.class.slide <- sliding.window.transform(GENOME.class, width = interval,
    jump = jump_size, type = 2, whole.data = TRUE)
GENOME.class.slide <- F_ST.stats(GENOME.class.slide, mode = "nucleotide")
FST_all_slide <- GENOME.class.slide@nuc.F_ST.vs.all
FST_pairwise_slide <- GENOME.class.slide@nuc.F_ST.pairwise
FST_all_slide_d <- as.data.frame(FST_all_slide)
#x axis
ids <- length(GENOME.class.slide@region.names)
xaxis <- seq(from = 1, to = ids, by = 1)

#Plot individual populations
for (i in seq_along(population_names)){
  file_slide <- paste(dir, "_", population_names[i],
  "_total_FST_per_sliding_window.pdf", sep = "")
  slide_plot <- ggplot(FST_all_slide_d, aes(x = xaxis,
      y = FST_all_slide_d[, i])) +
      geom_smooth(colour = "black", fill = "plum") + ggtitle(dir) +
      xlab("Contig coordinate (kbp)") + ylab("Total FST per interval") +
      scale_x_continuous(breaks = pretty(xaxis, n = 10))
  ggsave(file_slide, slide_plot)
  #write table with raw data
  slide_table <- paste(dir, "_", population_names[i],
  "_total_FST_per_sliding_window.txt", sep = "")
  write.table(FST_all_slide[, i], file = slide_table, sep = "\t",
  quote = FALSE, col.names = FALSE)
}

#Plot pairwise FST
for (i in seq(pairs)){
  FST_pairwise_slide_d <- as.data.frame(as.vector(FST_pairwise_slide[i, ]))
  labelling <- gsub("/", "_vs_", row.names(FST_pairwise_slide)[i])
  file_hist <- paste(dir, "_pairwise_FST_per_sliding_window_", labelling,
  ".pdf", sep = "")
  slide_plot <- ggplot(FST_pairwise_slide_d, aes(x = xaxis,
      y = FST_pairwise_slide_d[, 1])) +
      geom_smooth(colour = "black", fill = "slateblue") +
      ggtitle(dir) + xlab("Contig coordinate (kbp)") +
      ylab("Pairwise FST per interval") +
      scale_x_continuous(breaks = pretty(xaxis, n = 10))
  ggsave(file_slide, slide_plot)
  #write table with raw data
  slide_table <- paste(dir, "_pairwise_FST_per_sliding_window_", labelling,
  ".txt", sep = "")
  fst_table <- cbind(GENOME.class.slide@region.names, FST_pairwise_slide_d[, 1])
  write.table(fst_table, file = slide_table, sep = "\t", quote = FALSE,
  col.names = FALSE, row.names = FALSE)
}

}

##Genomewide results
for (i in seq_along(population_names)){
#Total FST
file_table2 <- paste("genome_", population_names[i],
"_total_FST_per_gene_all.txt", sep = "")
x <- as.data.frame(read.delim(file_table2))
file_hist <- paste("genome_", population_names[i],
"_total_FST_per_gene_all.pdf", sep = "")
fst_plot <- ggplot(x, aes(x = x[, 3])) +
geom_histogram(colour = "black", fill = "darkseagreen") +
xlab(expression(paste("Total FST per gene"))) + ylab("Number of genes") +
scale_x_continuous(breaks = pretty(x[, 3], n = 10))
ggsave(file_hist, fst_plot)
}

#Hudson KST
file_table2 <- paste("genome_Hudson_KST_per_gene_all", ".txt", sep = "")
x <- as.data.frame(read.delim(file_table2))
file_hist <- paste("genome_Hudson_KST_per_gene_all", ".pdf", sep = "")
fst_plot <- ggplot(x, aes(x = x[, 2])) +
geom_histogram(colour = "black", fill = "springgreen") +
xlab(expression(paste("Hudson KST per gene"))) + ylab("Number of genes") +
scale_x_continuous(breaks = pretty(x[, 2], n = 10))
ggsave(file_hist, fst_plot)


for (i in seq(pairs)){
  #Pairwise FST
  labelling <- gsub("/", "_vs_", row.names(FST_pairwise)[i])
  file_table2 <- paste("genome_pairwise_FST_per_gene_", labelling, ".txt",
  sep = "")
  x <- as.data.frame(read.delim(file_table2))
  file_hist <- paste("genome_pairwise_FST_per_gene_", labelling, ".pdf",
  sep = "")
  fst_plot <- ggplot(x, aes(x = x[, 2])) +
  geom_histogram(colour = "black", fill = "cadetblue") +
  xlab(expression(paste("Pairwise FST per gene"))) + ylab("Number of genes") +
  scale_x_continuous(breaks = pretty(x[, 2], n = 10))
  ggsave(file_hist, fst_plot)
}

```

-->
