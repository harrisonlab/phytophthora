
```R
#install.packages("pheatmap", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

#install and load libraries
require("pheatmap")             
require(data.table)

#load tables into a "list of lists"
qq <- lapply(list.files("alignment/star/P.cactorum/414_v2/DeSeq","PRO.*_featurecounts.txt$",full.names=T,recursive=T),function(x) fread(x))

# ensure the samples column is the same name as the treament you want to use:
qq[7]

#mm <- qq%>%Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by=c("Geneid","Chr","Start","End","Strand","Length")), .)

#merge the "list of lists" into a single table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

#convert data.table to data.frame for use with DESeq2
countData <- data.frame(m[,c(1,7:(ncol(m))),with=F])
rownames(countData) <- countData[,1]
countData <- countData[,-1]

#indexes <- unique(gsub("(.*)_L00.*", "\\1", colnames(countData)))
indexes <- c("PRO1467_S4", "PRO1467_S5", "PRO1467_S6", "PRO1467_S14", "PRO1467_S15", "PRO1467_S22", "PRO1467_S7", "PRO1467_S20", "PRO1467_S9", "PRO1467_S17", "PRO1467_S18", "PRO1467_S28")

#countDataSubset <- sapply(indexes, function(xx) rowSums(countData[,grep(paste(xx,'_', sep = ""), names(countData)), drop=FALSE]))
# If featurecounts has multimapping reads divided by the number of mapping sites, these will be rounded:
#countDataSubset <- round(countDataSubset,0)
countData <- round(countData,0)


#output countData
write.table(countData,"alignment/star/P.cactorum/414_v2/DeSeq/countData.txt",sep="\t",na="",quote=F)
# output countData with technical reps combined
write.table(countDataSubset,"alignment/star/P.cactorum/414_v2/DeSeq/countDataCombined.txt",sep="\t",na="",quote=F)

#output gene details
write.table(m[,1:6,with=F],"alignment/star/P.cactorum/414_v2/DeSeq/genes.txt",sep="\t",quote=F,row.names=F)
# colnames(countData) <- sub("X","",colnames(countData)) countData <- countData[,colData$Sample]
```

Running DeSeq2

```R
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
require(DESeq2)

unorderedColData <- read.table("alignment/star/P.cactorum/414_v2/DeSeq/P.cactorum_RNAseq_design_parsed.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])
# unorderedData <- read.table("alignment/star/P.cactorum/414_v2/DeSeq/countData.txt",header=T,sep="\t")
unorderedData <- read.table("alignment/star/P.cactorum/414_v2/DeSeq/countDataCombined.txt",header=T,sep="\t")
countData <- data.frame(unorderedData[ , order(colnames(unorderedData))])
# colData$Group <- paste0(colData$Isolate,colData$Plant.Line,colData$Rep,colData$Flowcell,colData$Timepoint)
colData$Group <- paste0(colData$Isolate,colData$Plant.Line,colData$Timepoint)

design <- ~Group
#design <- colData$Group

#dds <-     DESeqDataSetFromMatrix(countDataSubset,colData,design)
dds <-     DESeqDataSetFromMatrix(countData,colData,design)*/
#sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds, type = c("iterate")))
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds, type = c("ratio")))
dds <- DESeq(dds, fitType="local")


#  #Run DESeq2 removing an outlier
#
#  library(DESeq2)
#  colData <- read.table("colData",header=T,sep="\t")
#  countData <- read.table("countData2",header=T,sep="\t")
#
#  colData$Group <- paste0(colData$Strain,colData$Light,colData$Time)
#  #Eliminate Frq08_DD24_rep3 sample from colData and countData
#  colData <- colData[!(colData$Sample=="Frq08_DD24_rep3"),]      
#  countData <- subset(countData, select=-Frq08_DD24_rep3)
#
#  design <- ~Group
#  dds <-  DESeqDataSetFromMatrix(countData,colData,design)
#  sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds))
#  dds <- DESeq(dds, fitType="local")
#
```

Sample Distances

```R
library("RColorBrewer")
install.packages("gplots)
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")
# install.packages("ggrepel", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu")
library("ggrepel")

vst<-varianceStabilizingTransformation(dds)

pdf("alignment/star/P.cactorum/414_v2/DeSeq/heatmap_vst.pdf", width=12,height=12)
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

# Sample distances measured with rlog transformation:

rld <- rlog( dds )

pdf("alignment/star/P.cactorum/414_v2/DeSeq/heatmap_rld.pdf")
sampleDists <- dist( t( assay(rld) ) )
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- paste(rld$Group)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()
```

PCA plots

```R
#vst<-varianceStabilizingTransformation(dds)
pdf("alignment/star/P.cactorum/414_v2/DeSeq/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Isolate", "Plant.Line", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/star/P.cactorum/414_v2/DeSeq/PCA_rld.pdf")
#plotPCA(rld,intgroup=c("Isolate", "Plant.Line", "Timepoint"))
plotPCA(rld,intgroup=indexes)
dev.off()

pdf("alignment/star/P.cactorum/414_v2/DeSeq/PCA_additional.pdf")

dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/star/P.cactorum/414_v2/DeSeq/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)
```

Analysis of gene expression
<!--
--Emily0 hours","P414Emily48 hours

```R
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","--Emily0 hours","P414Emily48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.cactorum/414_v2/DeSeq/--Emily0h_vs_P414Emily48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/414_v2/DeSeq/--Emily0h_vs_P414Emily48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/414_v2/DeSeq/--Emily0h_vs_P414Emily48_down.txt",sep="\t",na="",quote=F)
```

"--Fenella0 hours","P414Fenella48 hours"

```R
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","--Fenella0 hours","P414Fenella48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.upregulated <- sig.res[order(sig.res$log2FoldChange, decreasing = TRUE),]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.cactorum/414_v2/DeSeq/--Fenella0h_vs_P414Fenella48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/414_v2/DeSeq/--Fenella0h_vs_P414Fenella48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/414_v2/DeSeq/--Fenella0h_vs_P414Fenella48_down.txt",sep="\t",na="",quote=F)
```

--Emily0 hours","P414Emily12 hours

```R
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","--Emily0 hours","P414Emily12 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.cactorum/414_v2/DeSeq/--Emily0h_vs_P414Emily12.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/414_v2/DeSeq/--Emily0h_vs_P414Emily12_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/414_v2/DeSeq/--Emily0h_vs_P414Emily12_down.txt",sep="\t",na="",quote=F)
```

"--Fenella0 hours","P414Fenella12 hours"

```R
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","--Fenella0 hours","P414Fenella12 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.upregulated <- sig.res[order(sig.res$log2FoldChange, decreasing = TRUE),]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.cactorum/414_v2/DeSeq/--Fenella0h_vs_P414Fenella12.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/414_v2/DeSeq/--Fenella0h_vs_P414Fenella12_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/414_v2/DeSeq/--Fenella0h_vs_P414Fenella12_down.txt",sep="\t",na="",quote=F)
```
-->


"P414Emily12 hours","P414Emily48 hours"

```R
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","P414Emily12 hours","P414Emily48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.cactorum/414_v2/DeSeq/P414Emily12h_vs_P414Emily48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/414_v2/DeSeq/P414Emily12h_vs_P414Emily48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/414_v2/DeSeq/P414Emily12h_vs_P414Emily48_down.txt",sep="\t",na="",quote=F)
```

"P414Fenella12 hours,"P414Fenella48 hours"

```R
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","P414Fenella12 hours","P414Fenella48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.upregulated <- sig.res[order(sig.res$log2FoldChange, decreasing = TRUE),]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res.upregulated,"alignment/star/P.cactorum/414_v2/DeSeq/P414Fenella12h_vs_P414Fenella48.txt",sep="\t",na="",quote=F)
write.table(sig.res,"alignment/star/P.cactorum/414_v2/DeSeq/P414Fenella12h_vs_P414Fenella48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/414_v2/DeSeq/P414Fenella12h_vs_P414Fenella48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/414_v2/DeSeq/P414Fenella12h_vs_P414Fenella48_down.txt",sep="\t",na="",quote=F)
```

"P414Fenella48 hours","P414Emily48 hours"

```R
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","P414Fenella48 hours","P414Emily48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.upregulated <- sig.res[order(sig.res$log2FoldChange, decreasing = TRUE),]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.cactorum/414_v2/DeSeq/P414Fenella48h_vs_P414Emily48.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/414_v2/DeSeq/P414Fenella48h_vs_P414Emily48_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/414_v2/DeSeq/P414Fenella48h_vs_P414Emily48_down.txt",sep="\t",na="",quote=F)
```

"P414Fenella12 hours","P414Emily12 hours"

```R
alpha <- 0.05
res= results(dds, alpha=alpha,contrast=c("Group","P414Fenella12 hours","P414Emily12 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.upregulated <- sig.res[order(sig.res$log2FoldChange, decreasing = TRUE),]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/star/P.cactorum/414_v2/DeSeq/P414Fenella12h_vs_P414Emily12.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/star/P.cactorum/414_v2/DeSeq/P414Fenella12h_vs_P414Emily12_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/star/P.cactorum/414_v2/DeSeq/P414Fenella12h_vs_P414Emily12_down.txt",sep="\t",na="",quote=F)
```

Make a table of raw counts, normalised counts and fpkm values:

```R
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/star/P.cactorum/414_v2/DeSeq/raw_counts2.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/star/P.cactorum/414_v2/DeSeq/normalised_counts.txt",sep="\t",na="",quote=F)

library(Biostrings)
library(naturalsort)
mygenes <- readDNAStringSet("gene_pred/final_ncbi/P.cactorum/414_v2/final_ncbi/414_v2_genes_incl_ORFeffectors_renamed.cdna.fasta")
t1 <- counts(dds)
t1 <- mygenes[rownames(t1)]
rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)


# robust may be better set at fasle to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.cactorum/414_v2/DeSeq/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/star/P.cactorum/414_v2/DeSeq/fpkm_counts.txt",sep="\t",na="",quote=F)

```
