#!/usr/bin/Rscript

# Plot a 5-way Venn diagram from a tab delimited file containing a matrix showing
 # presence /absence of orthogroups within 5 isolates.

 # This is intended to be used on the output of the orthoMCL pipeline following
 # building of the matrix using:
 # ~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL/orthoMCLgroups2tab.py

# This script requires the optparse R package. This can be downloaded by opening
# R and running the following command:
# install.packages("optparse",repos="http://cran.uk.r-project.org")
# When given the option, install this package to a local library.

# The script also requires the colorspace package. This can be downloaded by
# opening R and running the following command:
# install.packages(colorspace)

#get config options
library(optparse)
library(colorspace)
library(VennDiagram, lib.loc="/home/armita/R-packages/")
opt_list = list(
    make_option("--inp", type="character", help="tab seperated file containing matrix of presence of orthogroups"),
    make_option("--out", type="character", help="output venn diagram in pdf format")
)
opt = parse_args(OptionParser(option_list=opt_list))
f = opt$inp
o = opt$out

DEGs <- data.frame()
DEGs <- read.table(f, header = TRUE, sep='\t')

df1 <- DEGs
summary(df1)


area1=sum(df1[, "Mycelium_Mycelium"])
area2=sum(df1[, "Emily_12_hours"])
area3=sum(df1[, "Emily_48_hours"])
area4=sum(df1[, "Fenella_12_hours"])
area5=sum(df1[, "Fenella_48_hours"])

#print(area1, area2, area3, area4, area5)

colname1 <- paste("Mycelium_Mycelium")
colname2 <- paste("Emily_12_hours")
colname3 <- paste("Emily_48_hours")
colname4 <- paste("Fenella_12_hours")
colname5 <- paste("Fenella_48_hours")

# Species abreviation labels
label1 <- paste(colname1, sep="" )
label2 <- paste(colname2, sep="" )
label3 <- paste(colname3, sep="" )
label4 <- paste(colname4, sep="" )
label5 <- paste(colname5, sep="" )

n12=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_12_hours"] == 1))
n13=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_48_hours"] == 1))
n14=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Fenella_12_hours"] == 1))
n15=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Fenella_48_hours"] == 1))
n23=nrow(subset(df1, df1[,"Emily_12_hours"] == 1 & df1[,"Emily_48_hours"] == 1))
n24=nrow(subset(df1, df1[,"Emily_12_hours"] == 1 & df1[,"Fenella_12_hours"] == 1))
n25=nrow(subset(df1, df1[,"Emily_12_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n34=nrow(subset(df1, df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_12_hours"] == 1))
n35=nrow(subset(df1, df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n45=nrow(subset(df1, df1[,"Fenella_12_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n123=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_12_hours"] == 1 & df1[,"Emily_48_hours"] == 1))
n124=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_12_hours"] == 1 & df1[,"Fenella_12_hours"] == 1))
n125=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_12_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n134=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_12_hours"] == 1))
n135=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n145=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Fenella_12_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n234=nrow(subset(df1, df1[,"Emily_12_hours"] == 1 & df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_12_hours"] == 1))
n235=nrow(subset(df1, df1[,"Emily_12_hours"] == 1 & df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n245=nrow(subset(df1, df1[,"Emily_12_hours"] == 1 & df1[,"Fenella_12_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n345=nrow(subset(df1, df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_12_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n1234=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_12_hours"] == 1 & df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_12_hours"] == 1))
n1235=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_12_hours"] == 1 & df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n1245=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_12_hours"] == 1 & df1[,"Fenella_12_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n1345=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_12_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n2345=nrow(subset(df1, df1[,"Emily_12_hours"] == 1 & df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_12_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))
n12345=nrow(subset(df1, df1[,"Mycelium_Mycelium"] == 1 & df1[,"Emily_12_hours"] == 1 & df1[,"Emily_48_hours"] == 1 & df1[,"Fenella_12_hours"] == 1 & df1[,"Fenella_48_hours"] == 1))

summary(n12)
summary(n123)
summary(n1234)
summary(n12345)

pdf(o)
draw.quintuple.venn(
  area1, area2, area3, area4, area5,
  n12, n13, n14, n15, n23, n24, n25, n34, n35, n45,
  n123, n124, n125, n134, n135, n145, n234, n235, n245, n345,
  n1234, n1235, n1245, n1345, n2345,
  n12345,
  category = c(label1, label2, label3, label4, label5),
  lwd = rep(2, 5),
	lty = rep("solid", 5),
  col = rep("black", 5),
  fill = c(rainbow_hcl(5)),
  alpha = rep(0.5, 5),
  label.col = rep("black", 31),
  cex = rep(1, 31),
  fontface = rep("plain", 31),
  fontfamily = rep("serif", 31),
  cat.pos = c(0, 310, 215, 145, 60),
  cat.dist = rep(0.25, 5),
  cat.col = rep("black", 5),
  cat.cex = rep(1, 5),
  cat.fontface = rep("italic", 5),
  cat.fontfamily = rep("serif", 5),
  cat.just = rep(list(c(0.5, 0.5)), 5),
  rotation.degree = 0,
  rotation.centre = c(0.5, 0.5),
  ind = TRUE,
  margin = 0.15
)

dev.off()

warnings()
q()
