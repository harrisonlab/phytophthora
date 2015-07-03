# Building a Venn Diagram

# Install the R package fro Venn Diagrams
# mkdir -p /home/armita/R-packages
# R
# > install.packages("VennDiagram", lib="/home/armita/R-packages/")
# > library(VennDiagram, lib.loc="/home/armita/R-packages/")

cat analysis/inparanoid/summary_tables/new_final_tab.csv | sed 's/\+/1/g' | sed 's/\-/0/g' > analysis/inparanoid/summary_tables/new_final_tab2.csv

R
library(VennDiagram, lib.loc="/home/armita/R-packages/")
data1 <- read.delim(file="analysis/inparanoid/summary_tables/new_final_tab2.csv", header=T, sep="\t")
attach(data1)

X10300_404 <- subset(data1$X404, X10300==1 & X404==1)

 
area1=sum(data1$X10300)
area2=sum(data1$X404)
area3=sum(data1$X414)
n12=nrow(subset(data1, X10300==1 & X404==1))
n23=nrow(subset(data1, X404==1 & X414==1))
n13=nrow(subset(data1, X10300==1 & X414==1))
n123=nrow(subset(data1, X10300==1 & X404==1 & X414==1))
pdf('analysis/inparanoid/summary_tables/venn_10300_404_414.pdf')
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
    category = c("10300", "404", "414"), rotation = 1, reverse = FALSE, euler.d = TRUE,
    scaled = TRUE, lwd = rep(3, 3), lty = rep("solid", 3),
    col = rep("black", 3), fill = c("yellow", "darkgreen", "darkgoldenrod"), alpha = rep(0.6, 3),
    label.col = rep("black", 7), cex = rep(3, 7), fontface = rep("bold", 7),
    fontfamily = rep("serif", 7), cat.pos = c(-40, 40, 180),
    cat.dist = c(0.05, 0.05, 0.025), cat.col = rep("black", 3),
    cat.cex = rep(2, 3), cat.fontface = rep("bold", 3),
    cat.fontfamily = rep("serif", 3),
    cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0)), cat.default.pos = "outer",
    cat.prompts = FALSE, rotation.degree = 0, rotation.centre = c(0.5, 0.5),
    ind = TRUE, sep.dist = 0.05, offset = 0)
dev.off()

# > library(VennDiagram, lib.loc="/home/armita/R-packages/")
# > draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
#     category = rep("", 3), rotation = 1, reverse = FALSE, euler.d = TRUE,
#     scaled = TRUE, lwd = rep(2, 3), lty = rep("solid", 3),
#     col = rep("black", 3), fill = NULL, alpha = rep(0.5, 3),
#     label.col = rep("black", 7), cex = rep(1, 7), fontface = rep("plain", 7),
#     fontfamily = rep("serif", 7), cat.pos = c(-40, 40, 180),
#     cat.dist = c(0.05, 0.05, 0.025), cat.col = rep("black", 3),
#     cat.cex = rep(1, 3), cat.fontface = rep("plain", 3),
#     cat.fontfamily = rep("serif", 3),
#     cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0)), cat.default.pos = "outer",
#     cat.prompts = FALSE, rotation.degree = 0, rotation.centre = c(0.5, 0.5),
#     ind = TRUE, sep.dist = 0.05, offset = 0, ...)