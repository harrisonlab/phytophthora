#!/usr/bin/Rscript

##!/home/armita/prog/R/R-3.2.2/bin/Rscript

# R script for analysing a file of abundance of kmers.
# The input file should be a single column of abundance values.
# This script will produce a histogram of kmer abundance, a distribution plot
# and generate summary statistics.

#get config options
library(optparse)
library(scales)
library(RColorBrewer)
library(ggplot2)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)


opt_list <- list(
    make_option("--inp_all", type="character", help="tab seperated file containing 5' and 3' intergenic lengths"),
    make_option("--rxlr_subset", type="character", help="tab seperated file containing a subset of 5' and 3' intergenic lengths"),
    make_option("--crn_subset", type="character", help="tab seperated file containing a subset of 5' and 3' intergenic lengths"),
    make_option("--cazy_sec_subset", type="character", help="tab seperated file containing a subset of 5' and 3' intergenic lengths"),
    make_option("--PI_subset", type="character", help="tab seperated file containing a subset of 5' and 3' intergenic lengths"),
    make_option("--busco_subset", type="character", help="tab seperated file containing a subset of 5' and 3' intergenic lengths"),
    make_option("--cazy_nonsec_subset", type="character", help="tab seperated file containing a subset of 5' and 3' intergenic lengths"),
    make_option("--elicitin_subset", type="character", help="tab seperated file containing a subset of 5' and 3' intergenic lengths"),
    make_option("--prefix", type="character", help="output gene density plot in pdf format")
)
opt = parse_args(OptionParser(option_list=opt_list))

o = opt$prefix

o1 = paste(o, "_background.pdf", sep = "")
o2 = paste(o, "_subset.pdf", sep = "")

library(ggplot2)

df_all <- read.delim(file=opt$inp_all, header=F, sep="\t")
df_rxlr <- read.delim(file=opt$rxlr_subset, header=F, sep="\t")
df_crn <- read.delim(file=opt$crn_subset, header=F, sep="\t")
df_cazy_sec <- read.delim(file=opt$cazy_sec_subset, header=F, sep="\t")
df_pi <- read.delim(file=opt$PI_subset, header=F, sep="\t")
df_nlp <- read.delim(file=opt$busco_subset, header=F, sep="\t")
df_busco <- read.delim(file=opt$cazy_nonsec_subset, header=F, sep="\t")
df_cazy_nonsec <- read.delim(file=opt$elicitin_subset, header=F, sep="\t")
df_elicitin <- read.delim(file=opt$inp_all, header=F, sep="\t")

colnames(df_all) <- c("ID", "five_IG", "three_IG")
colnames(df_rxlr) <- c("ID", "five_IG", "three_IG")
colnames(df_crn) <- c("ID", "five_IG", "three_IG")
colnames(df_cazy_sec) <- c("ID", "five_IG", "three_IG")
colnames(df_pi) <- c("ID", "five_IG", "three_IG")
colnames(df_nlp) <- c("ID", "five_IG", "three_IG")
colnames(df_busco) <- c("ID", "five_IG", "three_IG")
colnames(df_cazy_nonsec) <- c("ID", "five_IG", "three_IG")
colnames(df_elicitin) <- c("ID", "five_IG", "three_IG")

density_plot <- ggplot(df_all, aes(df_all$five_IG, df_all$three_IG)) +
    stat_bin2d(bins = 100) +
    scale_y_continuous(trans='log2', limits=c(10,40826)) +
    scale_x_continuous(trans='log2', limits=c(10,40826)) +
    xlab("5' IG length") +
    ylab("3' IG length") +
    scale_fill_gradientn(colours=r)
ggsave(o1, density_plot, dpi=300, height=10, width=12)
#df2 <- read.delim(file=f_s, header=F, sep="\t")
#colnames(df2) <- c("ID", "five_IG", "three_IG")
#density_plot2 <-  density_plot + scale_fill_gradientn(colours=c("grey", "grey")) +
#    geom_point(data = df2, aes(df2$five_IG, df2$three_IG), colour = "red")
#ggsave(o2, density_plot2, dpi=300, height=10, width=12)



q()
