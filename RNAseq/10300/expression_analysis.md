

# DeSeq analysis of RNAseq data generated at NIBIO


The following commands were used to analyse the infection timecourse of
P. cactorum infecting F. vesca as well as controls including P. cactorum grown on
V8 juice media. Infected samples did not show elevated levels of read alignment
in comparison to mock infections. Those reads that did align were determined,
from the analysis below, to be a result of low-levels of contaminant reads from
the V8 media.

Evidence:
 * Reads in control samples
 * Samples with high expression at V8 showed high expression in controls and treatments


A second analysis is also presented below to normalise the expression values of
V8 samples on their own. This analysis was used for publication.

## Analysis using all RNAseq samples

```R

#===============================================================================
#       Load libraries
#===============================================================================

# Install packages if requied:
install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages("ggrepel", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu")
install.packages("gplots)
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")


# Load libraries
require("pheatmap")             
require(data.table)
require(DESeq2)
library("RColorBrewer")
#library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")
library("ggrepel")
library(Biostrings)
library(naturalsort)

#===============================================================================
#       Import data
#===============================================================================

#load tables into a "list of lists"
qq <- lapply(list.files("alignment/star/P.cactorum/10300/DeSeq","Sample_.*_featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))

# ensure the samples column is the same name as the treament you want to use:
qq[7]

#mm <- qq%>%Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by=c("Geneid","Chr","Start","End","Strand","Length")), .)

#merge the "list of lists" into a single table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

#convert data.table to data.frame for use with DESeq2
countData <- data.frame(m[,c(1,7:(ncol(m))),with=F])
rownames(countData) <- countData[,1]
countData <- countData[,-1]

#===============================================================================
#       Select a subset of samples (if required)
#===============================================================================


# indexes <- unique(gsub("(.*)_L00.*", "\\1", colnames(countData)))
#indexes <- c("Sample_13","Sample_14","Sample_15","Sample_16", "Sample_17", "Sample_18")

#countDataSubset <- sapply(indexes, function(xx) rowSums(countData[,grep(paste(xx,'_', sep = ""), names(countData)), drop=FALSE]))


#===============================================================================
#       Write a table of feature counts
#===============================================================================

#output countData
write.table(countData,"alignment/star/P.cactorum/10300/DeSeq/countData.txt",sep="\t",na="",quote=F)
# output countData with technical reps combined
#write.table(countDataSubset,"alignment/star/P.cactorum/10300/DeSeq/countDataCombined.txt",sep="\t",na="",quote=F)

#output gene details
write.table(m[,1:6,with=F],"alignment/star/P.cactorum/10300/DeSeq/genes.txt",sep="\t",quote=F,row.names=F)
# colnames(countData) <- sub("X","",colnames(countData)) countData <- countData[,colData$Sample]



#===============================================================================
#       Running analysis in DeSeq2 - importing data
#===============================================================================
# Import feature counts and information on gene models
# The columns of treatments in the feature counts has to be in the same order
# as the rows of treatments in the experimental design table
# As such, datasets are re-ordered before being combined.

unorderedColData <- read.table("alignment/star/P.cactorum/10300/DeSeq/P.cac_10300_RNAseq_design.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])
unorderedData <- read.table("alignment/star/P.cactorum/10300/DeSeq/countData.txt",header=T,sep="\t")
#unorderedData <- read.table("alignment/star/P.cactorum/10300/DeSeq/countDataCombined.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])


#===============================================================================
#       Setting the experimental design
#===============================================================================

colData$Group <- paste0(colData$Isolate,'_', colData$Timepoint)
design <- ~Group



#===============================================================================
#       Setting "size factors" for normalisation between treatments
#===============================================================================

dds <-     DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds, type = c("ratio")))
dds <- DESeq(dds, fitType="local")




#===============================================================================
#       Transforming data
#===============================================================================
# two approaches to normalisation were attempted:

vst<-varianceStabilizingTransformation(dds)
rld <- rlog( dds )

#===============================================================================
#       Investigating similarity between samples
#===============================================================================

#=
# Through heatmaps:
#=

pdf("alignment/star/P.cactorum/10300/DeSeq/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix,
  trace="none",  # turns off trace lines inside the heat map
  col=colours, # use on color palette defined earlier
  margins=c(12,12), # widens margins around plot
  srtCol=45,
  srtCol=45)
dev.off()

pdf("alignment/star/P.cactorum/10300/DeSeq/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

#=
# Through PCA plots:
#=


pdf("alignment/star/P.cactorum/10300/DeSeq/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Isolate", "Timepoint"))
dev.off()

pdf("alignment/star/P.cactorum/10300/DeSeq/PCA_rld.pdf")
plotPCA(rld,intgroup=c("Isolate","Timepoint"))
dev.off()


#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/star/P.cactorum/10300/DeSeq/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)


#===============================================================================
#       Differential gene expression
#===============================================================================

#=
# "Pcac_V8_media","Pcac_24_hpi":
#=

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Pcac_V8_media","Pcac_24_hpi"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.cactorum/10300/DeSeq/Pcac_V8_media_vs_Pcac_24_hpi.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/10300/DeSeq/Pcac_V8_media_vs_Pcac_24_hpi_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/10300/DeSeq/Pcac_V8_media_vs_Pcac_24_hpi_down.txt",sep="\t",na="",quote=F)

#=
# "Pcac_V8_media","Pcac_72_hpi":
#=

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Pcac_V8_media","Pcac_72_hpi"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.cactorum/10300/DeSeq/Pcac_V8_media_vs_Pcac_72_hpi.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/10300/DeSeq/Pcac_V8_media_vs_Pcac_72_hpi_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/10300/DeSeq/Pcac_V8_media_vs_Pcac_72_hpi_down.txt",sep="\t",na="",quote=F)

#=
# "Pcac_24_hpi","Pcac_72_hpi":
#=

alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","Pcac_24_hpi","Pcac_72_hpi"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.cactorum/10300/DeSeq/Pcac_24_hpi_vs_Pcac_72_hpi.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/10300/DeSeq/Pcac_24_hpi_vs_Pcac_72_hpi_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/10300/DeSeq/Pcac_24_hpi_vs_Pcac_72_hpi_down.txt",sep="\t",na="",quote=F)



#===============================================================================
#       Make a table of raw counts, normalised counts and fpkm values:
#===============================================================================

#=
# Raw counts:
#=

raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/star/P.cactorum/10300/DeSeq/raw_counts.txt",sep="\t",na="",quote=F)

#=
# Normalised counts:
#=
# Normalised counts may not always be most appropriate, particularly if a lot of 0s in
# background mapping reads from the host in samples with few reads from the pathogen lead to
# odd "size factors"

norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/star/P.cactorum/10300/DeSeq/normalised_counts.txt",sep="\t",na="",quote=F)

#=
# FPKM:
#=
# fpkm can either be counted by the total number of reads in the sample of esitumated
# from the "library normalisation factors" of the library ("robust" option). This will
# be influenced by your confidence in the normalisation accross your samples.

mygenes <- readDNAStringSet("gene_pred/annotation/P.cactorum/10300/10300_genes_incl_ORFeffectors.cdna.fasta")
t1 <- counts(dds)



#rownames(t1) <- gsub(".t.*", "", rownames(t1))

t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)

# Robust option
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.cactorum/10300/DeSeq/fpkm_norm_counts.txt",sep="\t",na="",quote=F)

# Total counts
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.cactorum/10300/DeSeq/fpkm_counts.txt",sep="\t",na="",quote=F)

```


## Analysis of samples from V8 juice



```R

#===============================================================================
#       Load libraries
#===============================================================================

# Install packages if requied:
install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )
install.packages("ggrepel", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu")
install.packages("gplots)
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
install.packages(naturalsort)


# Load libraries
require("pheatmap")             
require(data.table)
require(DESeq2)
library("RColorBrewer")
#library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")
library("ggrepel")
library(Biostrings)
library(naturalsort)

#===============================================================================
#       Import data
#===============================================================================

#load tables into a "list of lists"
qq <- lapply(list.files("alignment/star/P.cactorum/10300/DeSeq","Sample_.*_featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))

# ensure the samples column is the same name as the treament you want to use:
qq[7]

#merge the "list of lists" into a single table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

#convert data.table to data.frame for use with DESeq2
countData <- data.frame(m[,c(1,7:(ncol(m))),with=F])
rownames(countData) <- countData[,1]
countData <- countData[,-1]

#===============================================================================
#       Select a subset of samples (if required)
#===============================================================================

indexes <- c("Sample_16", "Sample_17", "Sample_18")
countDataSubset <- data.frame(countData[indexes])

#===============================================================================
#       Write a table of feature counts
#===============================================================================

#output countData
write.table(countDataSubset,"alignment/star/P.cactorum/10300/DeSeq/V8_countData.txt",sep="\t",na="",quote=F)

#output gene details
write.table(m[,1:6,with=F],"alignment/star/P.cactorum/10300/DeSeq/V8_genes.txt",sep="\t",quote=F,row.names=F)




#===============================================================================
#       Running analysis in DeSeq2 - importing data
#===============================================================================
# Import feature counts and information on gene models
# The columns of treatments in the feature counts has to be in the same order
# as the rows of treatments in the experimental design table
# As such, datasets are re-ordered before being combined.

unorderedColData <- read.table("alignment/star/P.cactorum/10300/DeSeq/P.cac_10300_RNAseq_design.txt",header=T,sep="\t")
rownames(unorderedColData) <- unorderedColData$Sample.name
unorderedColDataSubset <- unorderedColData[indexes,]

colData <- data.frame(unorderedColDataSubset[ order(unorderedColDataSubset$Sample.name),])
unorderedData <- read.table("alignment/star/P.cactorum/10300/DeSeq/V8_countData.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])


#===============================================================================
#       Setting the experimental design
#===============================================================================

#colData$Group <- paste0(colData$Isolate,'_', colData$Timepoint)
#design <- ~Group
design <- ~ 1

#===============================================================================
#       Setting "size factors" for normalisation between treatments
#===============================================================================

dds <-     DESeqDataSetFromMatrix(countData,colData,design)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds, type = c("ratio")))
dds <- DESeq(dds, fitType="local")




#===============================================================================
#       Transforming data
#===============================================================================
# two approaches to normalisation were attempted:

vst<-varianceStabilizingTransformation(dds)
rld <- rlog( dds )

#===============================================================================
#       Investigating similarity between samples
#===============================================================================

#=
# Through heatmaps:
#=

pdf("alignment/star/P.cactorum/10300/DeSeq/V8_heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix,
  trace="none",  # turns off trace lines inside the heat map
  col=colours, # use on color palette defined earlier
  margins=c(12,12), # widens margins around plot
  srtCol=45,
  srtCol=45)
dev.off()

pdf("alignment/star/P.cactorum/10300/DeSeq/V8_heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

#=
# Through PCA plots:
#=


pdf("alignment/star/P.cactorum/10300/DeSeq/V8_PCA_vst.pdf")
plotPCA(vst,intgroup=c("Isolate", "Timepoint"))
dev.off()

pdf("alignment/star/P.cactorum/10300/DeSeq/V8_PCA_rld.pdf")
plotPCA(rld,intgroup=c("Isolate","Timepoint"))
dev.off()


#Plot using rlog transformation, showing sample names:

#data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
#percentVar <- round(100 * attr(data, "percentVar"))

#pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
# geom_point(size=3) +
# xlab(paste0("PC1: ",percentVar[1],"% variance")) +
# ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
# coord_fixed()

#ggsave("alignment/star/P.cactorum/10300/DeSeq/V8_PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)


#===============================================================================
#       Make a table of raw counts, normalised counts and fpkm values:
#===============================================================================

#=
# Raw counts:
#=

raw_counts <- data.frame(counts(dds, normalized=FALSE))
# colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/star/P.cactorum/10300/DeSeq/V8_raw_counts.txt",sep="\t",na="",quote=F)

#=
# Normalised counts:
#=
# Normalised counts may not always be most appropriate, particularly if a lot of 0s in
# background mapping reads from the host in samples with few reads from the pathogen lead to
# odd "size factors"

norm_counts <- data.frame(counts(dds, normalized=TRUE))
#colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/star/P.cactorum/10300/DeSeq/V8_normalised_counts.txt",sep="\t",na="",quote=F)

#=
# FPKM:
#=
# fpkm can either be counted by the total number of reads in the sample of esitumated
# from the "library normalisation factors" of the library ("robust" option). This will
# be influenced by your confidence in the normalisation accross your samples.

mygenes <- readDNAStringSet("gene_pred/annotation/P.cactorum/10300/10300_genes_incl_ORFeffectors.cdna.fasta")
t1 <- counts(dds)


t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)

# Robust option
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
#colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.cactorum/10300/DeSeq/V8_fpkm_norm_counts.txt",sep="\t",na="",quote=F)

# Total counts
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
#colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.cactorum/10300/DeSeq/V8_fpkm_counts.txt",sep="\t",na="",quote=F)



```
