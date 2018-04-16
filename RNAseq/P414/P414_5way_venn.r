#!/usr/bin/Rscript

# Plot a 5-way Venn diagram from a tab delimited file containing a matrix showing
# presence of DEGs between 5 conditions of growth.

# This is intended to be used on the output of deseq2venn.py, after the creation of DEG lists by DeSeq

# The script also requires the colorspace package. This can be downloaded by
# opening R and running the following command:

#get config options
library("optparse")
library("colorspace")
library("VennDiagram")
opt_list <- list(
    make_option("--inp", type = "character",
    help = "tab seperated file containing matrix of DEGs"),
    make_option("--out", type = "character",
    help = "output venn diagram in pdf format")
)
opt <- parse_args(OptionParser(option_list = opt_list))
f <- opt$inp
o <- opt$out

DEGs <- data.frame()
DEGs <- read.table(f, header = TRUE)

df1 <- DEGs

df1$Myc <- 0
df1$E12 <- 0
df1$E48 <- 0
df1$F12 <- 0
df1$F48 <- 0

df1$Myc[DEGs[, "Mycelium_vs_Emily12h"] == 1 |
DEGs[, "Mycelium_vs_Emily48h"] == 1 |
DEGs[, "Mycelium_vs_Fenella12h"] == 1 |
DEGs[, "Mycelium_vs_Fenella48h"] == 1] <- 1


df1$E12[DEGs[, "Mycelium_vs_Emily12h"] == 1 |
DEGs[, "Emily12h_vs_Emily48h"] == 1 |
DEGs[, "Emily12h_vs_Fenella12h"] == 1 |
DEGs[, "Emily12h_vs_Fenella48h"] == 1] <- 1
df1$E48[DEGs[, "Mycelium_vs_Emily48h"] == 1 |
DEGs[, "Emily12h_vs_Emily48h"] == 1 |
DEGs[, "Fenella12h_vs_Emily48h"] == 1 |
DEGs[, "Emily48h_vs_Fenella48h"] == 1] <- 1
df1$F12[DEGs[, "Mycelium_vs_Fenella12h"] == 1 |
DEGs[, "Fenella12h_vs_Fenella48h"] == 1 |
DEGs[, "Emily12h_vs_Fenella12h"] == 1 |
DEGs[, "Fenella12h_vs_Emily48h"] == 1] <- 1
df1$F48[DEGs[, "Mycelium_vs_Fenella48h"] == 1 |
DEGs[, "Fenella12h_vs_Fenella48h"] == 1 |
DEGs[, "Emily12h_vs_Fenella48h"] == 1 |
DEGs[, "Emily48h_vs_Fenella48h"] == 1] <- 1

area1=sum(df1[, "Myc"])
area2=sum(df1[, "E12"])
area3=sum(df1[, "E48"])
area4=sum(df1[, "F48"])
area5=sum(df1[, "F12"])

# Species abreviation labels
label1 <- paste("Myc", sep="" )
label2 <- paste("E12", sep="" )
label3 <- paste("E48", sep="" )
label4 <- paste("F48", sep="" )
label5 <- paste("F12", sep="" )

n12=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E12"] == 1))
n13=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E48"] == 1))
n14=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"F48"] == 1))
n15=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"F12"] == 1))
n23=nrow(subset(df1, df1[,"E12"] == 1 & df1[,"E48"] == 1))
n24=nrow(subset(df1, df1[,"E12"] == 1 & df1[,"F48"] == 1))
n25=nrow(subset(df1, df1[,"E12"] == 1 & df1[,"F12"] == 1))
n34=nrow(subset(df1, df1[,"E48"] == 1 & df1[,"F48"] == 1))
n35=nrow(subset(df1, df1[,"E48"] == 1 & df1[,"F12"] == 1))
n45=nrow(subset(df1, df1[,"F48"] == 1 & df1[,"F12"] == 1))
n123=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E12"] == 1 & df1[,"E48"] == 1))
n124=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E12"] == 1 & df1[,"F48"] == 1))
n125=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E12"] == 1 & df1[,"F12"] == 1))
n134=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E48"] == 1 & df1[,"F48"] == 1))
n135=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E48"] == 1 & df1[,"F12"] == 1))
n145=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"F48"] == 1 & df1[,"F12"] == 1))
n234=nrow(subset(df1, df1[,"E12"] == 1 & df1[,"E48"] == 1 & df1[,"F48"] == 1))
n235=nrow(subset(df1, df1[,"E12"] == 1 & df1[,"E48"] == 1 & df1[,"F12"] == 1))
n245=nrow(subset(df1, df1[,"E12"] == 1 & df1[,"F48"] == 1 & df1[,"F12"] == 1))
n345=nrow(subset(df1, df1[,"E48"] == 1 & df1[,"F48"] == 1 & df1[,"F12"] == 1))
n1234=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E12"] == 1 & df1[,"E48"] == 1 & df1[,"F48"] == 1))
n1235=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E12"] == 1 & df1[,"E48"] == 1 & df1[,"F12"] == 1))
n1245=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E12"] == 1 & df1[,"F48"] == 1 & df1[,"F12"] == 1))
n1345=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E48"] == 1 & df1[,"F48"] == 1 & df1[,"F12"] == 1))
n2345=nrow(subset(df1, df1[,"E12"] == 1 & df1[,"E48"] == 1 & df1[,"F48"] == 1 & df1[,"F12"] == 1))
n12345=nrow(subset(df1, df1[,"Myc"] == 1 & df1[,"E12"] == 1 & df1[,"E48"] == 1 & df1[,"F48"] == 1 & df1[,"F12"] == 1))


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
  cat.dist = rep(0.2, 5),
  cat.col = rep("black", 5),
  cat.cex = rep(1, 5),
  cat.fontface = rep("plain", 5),
  cat.fontfamily = rep("serif", 5),
  cat.just = rep(list(c(0.5, 0.5)), 5),
  rotation.degree = 0,
  rotation.centre = c(0.5, 0.5),
  ind = TRUE,
  margin = 0.15
)

dev.off()
