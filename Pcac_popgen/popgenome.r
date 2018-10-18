

library(PopGenome)

GENOME.class <- readData("analysis/popgen/SNP_calling/popgenome/vcf", format = "VCF", gffpath="analysis/popgen/SNP_calling/popgenome/gff")
get.sum.data(GENOME.class)
GENOME.class@region.data

list1 <- GENOME.class@region.data@GeneSNPS
table(list1)['TRUE']
table(list1)['FALSE']

list2 <- GENOME.class@region.data@synonymous
#table(list2)['NaN']
table(list2)['TRUE']
list3 <- GENOME.class@region.data@CodingSNPS
table(list3)['TRUE'][1]
table(list3)['FALSE']

GENOME.class@region.data@populations2
GENOME.class <- set.populations(GENOME.class,
list(c("415","15_7","15_13","404","12420","2003_3","4032","PC13_15","414","4040","416"),c("62471","P295","R36_14"),c("371","SCRP370","SCRP376")), diploid=TRUE)
GENOME.class@region.data@populations

GENOME.class <- neutrality.stats(GENOME.class, FAST=TRUE)
get.neutrality(GENOME.class)

GENOME.class@Tajima.D
GENOME.class <- F_ST.stats(GENOME.class)
GENOME.class <- F_ST.stats.2(GENOME.class)

get.F_ST(GENOME.class)
GENOME.class@nucleotide.F_ST
get.F_ST(GENOME.class, pairwise=TRUE)
get.F_ST(GENOME.class, pairwise=TRUE)[1:3]

GENOME.class <- diversity.stats(GENOME.class)
get.diversity(GENOME.class)
GENOME.class@nuc.diversity.within

GENOME.class <- linkage.stats(GENOME.class)
get.linkage(GENOME.class)
get.linkage(GENOME.class)[1:3]
GENOME.class@Wall.B

GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="analysis/popgen/SNP_calling/popgenome/fasta/popgenome_input")
#list2 <- GENOME.class@region.data@synonymous
#table(list2)['NaN']
#table(list2)['FALSE']
GENOME.class <- MKT(GENOME.class, do.fisher.test=TRUE)
get.MKT(GENOME.class)
get.MKT(GENOME.class)[1]


GENOME.class <- detail.stats(GENOME.class)
get.detail(GENOME.class)


GENOME.class <- set.populations(GENOME.class, ref.chr="analysis/popgen/SNP_calling/popgenome/popgenome_input")



















# reading chromosome 1

vcf_handle <- vcf_open("analysis/popgen/SNP_calling/popgenome/input/414_v2_contigs_unmasked_filtered_no_errors.vcf")

GENOME.class <- readData("analysis/popgen/SNP_calling/popgenome/input/fasta", format = "fasta")

GENOME.class <- readData("analysis/popgen/SNP_calling/popgenome/input", format = "VCF")

GENOME.class <- readSNP("analysis/popgen/SNP_calling/popgenome/input")
# scan the data with consecutive windows
# window size: 1000 nucleotides (type=2)
# jump size: 1000 nucleotides (type=2)
> GENOME.class.slide <- sliding.window.transform(GENOME.class,1000,1000,type=2)
# calculate diversity statistics for all individuals
> GENOME.class.slide <- diversity.stats(GENOME.class.slide)
# Get the results ([[1]], because only one pop is defined)
> get.diversity(GENOME.class.slide)[[1]]
# alternative: directly access the nucleotide diversity
> plot(GENOME.class.slide@nuc.diversity.within)




# read chromosome 1 with the corresponding GFF-file
> GENOME.class <- readSNP("Arabidopsis", CHR=1, gffpath="Ara.gff")
# verify the nonsyn/syn SNPs (we need the reference sequence as a FASTA file!)
> GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="chr1.fas")
# split the data into exon regions
> GENOME.class.exons <- splitting.data(GENOME.class,subsites="exon")
# calculate the nonsynonymous diversities
> GENOME.class.exons <- diversity.stats(GENOME.class.exons, subsites="nonsyn")




appended_df <- read_delim("analysis/popgen/SNP_calling/414_v2_contigs_unmasked_filtered_no_errors.vcf",
     "\t", escape_double = FALSE, col_names = FALSE,
     trim_ws = TRUE)
colnames(appended_df) <- c("contig","position", "depth", "strain")
appended_df$depth_mod <- ifelse(appended_df$depth > 150, 150, appended_df$depth)
# appended_df$contig <- paste("Chr", appended_df$contig, sep = "")
# appended_df$strain <- factor(appended_df$strain, levels = c("12008", "12251", "12253", "12158", "12161"))


# install.packages("ggplot2")
library(ggplot2)
require(scales)
