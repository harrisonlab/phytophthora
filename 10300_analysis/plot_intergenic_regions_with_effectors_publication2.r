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
    make_option("--prefix", type="character", help="output gene density plot in pdf format"),
    make_option("--background_table", type="character", help="tab seperated file containing 5' and 3' intergenic lengths")
)
opt = parse_args(OptionParser(option_list=opt_list))


library(ggplot2)

df_all <- read.delim(file=opt$inp_all, header=F, sep="\t")


colnames(df_all) <- c("ID", "five_IG", "three_IG", "effector")
summary(df_all$effector)
df_all$effector<-as.factor(df_all$effector)
df_all$effector<-factor(df_all$effector, levels = c("Total", "RxLR", "CRN", "NLP", "Secreted CAZyme", "PI", "BUSCO", "Non-secreted CAZyme", "Elicitin"))
summary(df_all$effector)

#subset(df_all, effector=='RxLR')

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
#        geom_point(data = df_all, aes(df_all$five_IG, df_all$three_IG), colour = "red")
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



scatter_plot <- ggplot(df_all, aes(df_all$five_IG, df_all$three_IG)) +
    geom_point(size=0.5, colour = "red") +
    scale_y_continuous(trans='log2', limits=c(10,40826)) +
    scale_x_continuous(trans='log2', limits=c(10,40826)) +
    xlab("5' IG length") +
    ylab("3' IG length") +
    facet_wrap( ~ df_all$effector, nrow = 3, ncol = 3, strip.position = "top")

#density_plot <- ggplot(df_all[df_all$effector=='Total',], aes(df_all$five_IG, df_all$three_IG)) +
#    stat_bin2d(bins = 100) +
#    scale_y_continuous(trans='log2', limits=c(10,40826)) +
#    scale_x_continuous(trans='log2', limits=c(10,40826)) +
#    xlab("5' IG length") +
#    ylab("3' IG length")
#subset_df <- subset(df_all, effector!='Total')
#density_plot2 <-  density_plot + scale_fill_gradientn(colours=c("grey", "grey")) +
#    geom_point(data = subset_df, aes(five_IG, three_IG), colour = "red") +
#    facet_wrap( ~ df_all$effector, nrow = 3, ncol = 3)

density_plot1 <- ggplot(df_all, aes(df_all$five_IG, df_all$three_IG)) +
    stat_bin2d(bins = 100) +
    scale_y_continuous(trans='log2', limits=c(10,40826)) +
    scale_x_continuous(trans='log2', limits=c(10,40826)) +
    xlab("5' IG length") +
    ylab("3' IG length") +
    scale_fill_gradientn(colours=r) +
    facet_wrap( ~ df_all$effector, nrow = 3, ncol = 3)

df_total <- read.delim(file=opt$background_table, header=F, sep="\t")

colnames(df_total) <- c("ID", "five_IG", "three_IG", "effector")
df_total$effector<-as.factor(df_total$effector)
df_total$effector<-factor(df_total$effector, levels = c("Total", "RxLR", "CRN", "NLP", "Secreted CAZyme", "PI", "BUSCO", "Non-secreted CAZyme", "Elicitin"))

density_plot2 <- ggplot(df_total, aes(df_total$five_IG, df_total$three_IG)) +
    stat_bin2d(bins = 100) +
    scale_y_continuous(trans='log2', limits=c(10,40826)) +
    scale_x_continuous(trans='log2', limits=c(10,40826)) +
    xlab("5' IG length") +
    ylab("3' IG length") +
    scale_fill_gradientn(colours=c("grey", "grey")) +
    facet_wrap( ~ df_total$effector, nrow = 3, ncol = 3)


density_plot3 <- ggplot(df_all, aes(df_all$five_IG, df_all$three_IG)) +
    stat_bin2d(bins = 100) +
    scale_y_continuous(trans='log2', limits=c(10,40826)) +
    scale_x_continuous(trans='log2', limits=c(10,40826)) +
    xlab("5' IG length") +
    ylab("3' IG length") +
    scale_fill_gradientn(colours=c("white", "white")) +
    geom_point(data = df_all, aes(five_IG, three_IG), colour = "red", size = 0.5, shape = 15) +
    scale_y_continuous(trans='log2', limits=c(10,40826)) +
    scale_x_continuous(trans='log2', limits=c(10,40826)) +
    xlab("5' IG length") +
    ylab("3' IG length") +
    facet_wrap( ~ df_all$effector, nrow = 3, ncol = 3) +
    theme(panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank()
    )


o = opt$prefix

o1 = paste(o, "_plot_a.tiff", sep = "")
ggsave(o1, density_plot1, dpi=300, height=10, width=10)

o2 = paste(o, "_plot_b.tiff", sep = "")
ggsave(o2, density_plot2, dpi=300, height=10, width=10)

o3 = paste(o, "_plot_c.tiff", sep = "")
ggsave(o3, density_plot3, dpi=300, height=10, width=10)



q()
