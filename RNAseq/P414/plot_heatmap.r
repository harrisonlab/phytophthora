#!/usr/bin/Rscript

annotTab_df <- read.delim("~/Downloads/heatmap/414_annotation_ncbi2.tsv")
View(`annotTab_df`)

install.packages("gplots")
library("gplots")
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

#-----
# Heatmap plotting function
#-----


# heatmap.2(
#   rxlr_mat,
#   trace="none",
#   col=r       # use on color palette defined earlier
# )
# # This plots the heamap but we dont have control of the order of the columns

plot_heatmap <- function(matrix, prefix) {
  nrows <- dim(matrix)[1]
  png(paste("../", prefix, "_heatmap.png", sep=""),    # create PNG for the heat map
   width = 5*300,        # 5 x 300 pixels
   # height = (nrows*32)+400+200,
   height = 10*300,
   res = 300,            # 300 pixels per inch
   pointsize = 12)       # smaller font size

   heatmap.2(
     matrix,
     trace="none",
     Colv=FALSE, # no col dendrogram and columns ordered by input
     dendrogram="row", # no col dendrogram and columns ordered by input
  #   lmat=rbind( c(0, 3, 4), c(2,1,1 ) ), lwid=c(1.5, 4, 2 ), # Change position and size of legend
     cexRow=0.75,
     cexCol=1,
    lhei = c(1,5), # Adjust relative height of key to heatmap
    density.info="none", # remove density info from key
      margins = c(6, 6),
     col=r       # use on color palette defined earlier
   )

  dev.off()               # close the PNG device
}
#---
# End plot_heatmap function
#---


#---
# Plot RxLR expression
#---

rxlr_df <- subset(annotTab_df, putative_rxlr=="RxLR" & DEG=="DEG", select=c(1,40:45,47))
rxlr_mat <- as.matrix(rxlr_df[c(2,5,3,6,4,7)])
rownames(rxlr_mat) <- rxlr_df[,1]
colnames(rxlr_mat) <- c('E12 vs M', 'F12 vs M', 'E48 vs M', 'F48 vs M', 'E48 vs E12', 'F48 vs F12')

plot_heatmap(rxlr_mat, "RxLR")

#---
# Plot CRN expression
#---

crn_df <- subset(annotTab_df, putative_crn=="CRN" & DEG=="DEG", select=c(1,40:45,47))
crn_mat <- as.matrix(crn_df[c(2,5,3,6,4,7)])
rownames(crn_mat) <- crn_df[,1]
colnames(crn_mat) <- c('E12 vs M', 'F12 vs M', 'E48 vs M', 'F48 vs M', 'E48 vs E12', 'F48 vs F12')

plot_heatmap(crn_mat, "CRN")


#---
# Plot CAZY expression
#---

cazy_df <- subset(annotTab_df, cazy!="" & DEG=="DEG", select=c(1,40:45,47))
cazy_mat <- as.matrix(cazy_df[c(2,5,3,6,4,7)])
rownames(cazy_mat) <- cazy_df[,1]
colnames(cazy_mat) <- c('E12 vs M', 'F12 vs M', 'E48 vs M', 'F48 vs M', 'E48 vs E12', 'F48 vs F12')

plot_heatmap(cazy_mat, "CAZY")


#---
# Plot Transcription factor expression
#---

tf_df <- subset(annotTab_df, transcriptionfactor!="" & DEG=="DEG", select=c(1,40:45,47))
tf_mat <- as.matrix(tf_df[c(2,5,3,6,4,7)])
rownames(tf_mat) <- tf_df[,1]
colnames(tf_mat) <- c('E12 vs M', 'F12 vs M', 'E48 vs M', 'F48 vs M', 'E48 vs E12', 'F48 vs F12')

plot_heatmap(tf_mat, "transcription_factor")


#---
# Plot expanded/contracted family expression
#---

expansion_df <- subset(annotTab_df, expansion_status!="" & DEG=="DEG", select=c(1,40:45,47))
expansion_mat <- as.matrix(expansion_df[c(2,5,3,6,4,7)])
rownames(expansion_mat) <- expansion_df[,1]
colnames(expansion_mat) <- c('E12 vs M', 'F12 vs M', 'E48 vs M', 'F48 vs M', 'E48 vs E12', 'F48 vs F12')

plot_heatmap(expansion_mat, "expanded-contracted_family")


#---
# Plot iterate through CAZY families plotting expression
#---
plot_level <- function(element){
  func_df <- subset(annotTab_df, cazy==element & DEG=="DEG", select=c(1,40:45,47))
  func_mat <- as.matrix(func_df[c(2,5,3,6,4,7)])
  rownames(func_mat) <- func_df[,1]
  colnames(func_mat) <- c('E12 vs M', 'F12 vs M', 'E48 vs M', 'F48 vs M', 'E48 vs E12', 'F48 vs F12')
  if (dim(func_mat)[1] > 2){
    element <- gsub(":", "-", element)
    plot_heatmap(func_mat, element)
  }
}

x <- levels(`414_annotation_ncbi2`$cazy)

lapply(x[-1], plot_level)

#---
# Plot iterate through apoplastic effector families plotting expression
#---
plot_apo_level <- function(element){
  func_df <- subset(annotTab_df, ipr_effectors==element & DEG=="DEG", select=c(1,40:45,47))
  func_mat <- as.matrix(func_df[c(2,5,3,6,4,7)])
  rownames(func_mat) <- func_df[,1]
  colnames(func_mat) <- c('E12 vs M', 'F12 vs M', 'E48 vs M', 'F48 vs M', 'E48 vs E12', 'F48 vs F12')
  if (dim(func_mat)[1] > 2){
    element <- gsub(":", "-", element)
    plot_heatmap(func_mat, element)
  }
}

x <- levels(`414_annotation_ncbi2`$ipr_effectors)

lapply(x[-1], plot_apo_level)


#---
# Iterate through expanded/contracted families plotting expression
#---
plot_exp_level <- function(element){
  func_df <- subset(annotTab_df, clade_expansion==element & DEG=="DEG", select=c(1,40:45,47))
  func_mat <- as.matrix(func_df[c(2,5,3,6,4,7)])
  rownames(func_mat) <- func_df[,1]
  colnames(func_mat) <- c('E12 vs M', 'F12 vs M', 'E48 vs M', 'F48 vs M', 'E48 vs E12', 'F48 vs F12')
  if (dim(func_mat)[1] > 2){
    element <- gsub(" ", "_", element)
    plot_heatmap(func_mat, element)
  }
}

x <- levels(`414_annotation_ncbi2`$clade_expansion)

lapply(x[-1], plot_exp_level)
