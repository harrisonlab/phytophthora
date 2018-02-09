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
    make_option("--prefix", type="character", help="output gene density plot in pdf format")
)
opt = parse_args(OptionParser(option_list=opt_list))
o = opt$prefix

o1 = paste(o, "_background.pdf", sep = "")
o2 = paste(o, "_subset.pdf", sep = "")

library(ggplot2)

df_all <- read.delim(file=opt$inp_all, header=F, sep="\t")


colnames(df_all) <- c("ID", "five_IG", "three_IG", "effector")
summary(df_all$effector)
df_all$effector<-as.factor(df_all$effector)
df_all$effector<-factor(df_all$effector, levels = c("Total", "RxLR", "CRN", "NLP", "Secreted CAZyme", "PI", "BUSCO", "Non-secreted CAZyme", "Elicitin"))
summary(df_all$effector)

summary(df_all$effector)
#df_all[1,1]

#df_total <- data.frame(df_all[df_all$effector=='Total' ])
#df_total <- subset(df_all, effector=='Total')

#density_plot <- ggplot(df_total, aes(df_total$five_IG, df1$three_IG)) +
#    stat_bin2d(bins = 100) +
#    scale_y_continuous(trans='log2', limits=c(10,40826)) +
#    scale_x_continuous(trans='log2', limits=c(10,40826)) +
#    xlab("5' IG length") +
#    ylab("3' IG length") +
#    scale_fill_gradientn(colours=r)
#
#    + scale_fill_gradientn(colours=c("grey", "grey")) +
#        geom_point(data = df2, aes(df2$five_IG, df2$three_IG), colour = "red")
#
#density_plot2 <- density_plot + scale_fill_gradientn(colours=c("grey", "grey")) +
#    geom_point(data = df_all, aes(df_all$five_IG, df_all$three_IG), colour = "red") +
#    facet_wrap( ~ df_all$effector, nrow = 3, ncol = 3)

density_plot <- ggplot(df_all, aes(df_all$five_IG, df_all$three_IG)) +
    stat_bin2d(bins = 100) +
    scale_y_continuous(trans='log2', limits=c(10,40826)) +
    scale_x_continuous(trans='log2', limits=c(10,40826)) +
    xlab("5' IG length") +
    ylab("3' IG length") +
    facet_wrap( ~ df_all$effector, nrow = 3, ncol = 3) +
    scale_fill_gradientn(colours=r)
ggsave(o1, density_plot, dpi=300, height=10, width=12)


q()
