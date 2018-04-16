

```bash
export R_LIBS=/home/armita/R/x86_64-pc-linux-gnu-library/3.4
/home/deakig/R3.4/bin/R

```

```R
setwd("/data/scratch/armita/idris")
#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
load_all("~/pipelines/RNA-seq/scripts/myfunctions")
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)
```

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

```R
library(tximport)
library(rjson)
library(readr)

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files	    
txi.reps <- tximport(paste(list.dirs("alignment/salmon/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)	   

# get the sample names from the folders	    
mysamples <- list.dirs("alignment/salmon/DeSeq2",full.names=F,recursive=F)

# summarise to gene level (this can be done in the tximport step, but is easier to understand in two steps)
txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))
```


#==========================================================================================
#       Read sample metadata and annotations
#=========================================================================================
```R
# Read sample metadata	    
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts
unorderedColData <- read.table("alignment/salmon/DeSeq2/P.cactorum_RNAseq_design_parsed.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample.name),])

colData$Group <- paste0(colData$Plant.Line, " - ", colData$Timepoint)

# reorder colData for salmon 		 
#colData <- colData[mysamples,]

# reorder colData for featureCounts		 
#colData <- colData[colnames(countData),]

# get annotations     
#annotations <- fread("ip_annotations.txt", sep="\t",header=T)

```

Library normalisation

```R
# create dds object from Salmon counts and sample metadata (library size normalisation is taken from the length estimates)		 
dds <- DESeqDataSetFromTximport(txi.genes,colData,~1)
dds <- estimateSizeFactors(dds)
Group <- as.factor(dds$Group)
```


Sample Distances

```R
#install and load libraries
require("pheatmap")
require(data.table)
library("RColorBrewer")
#install.packages("gplots")
library("gplots", Sys.getenv("R_LIBS_USER"))
library("ggplot2")
# install.packages("ggrepel", Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu")
library("ggrepel")

vst<-varianceStabilizingTransformation(dds)

pdf("alignment/salmon/DeSeq2/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Group)
colnames(sampleDistMatrix) <- paste(vst$Group)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
#heatmap( sampleDistMatrix,
#  trace="none",  # turns off trace lines inside the heat map
#  col=colours, # use on color palette defined earlier
#  margins=c(12,12), # widens margins around plot
#  srtCol=45,
#  srtCol=45)
dev.off()

# Sample distances measured with rlog transformation:

rld <- rlog( dds )

pdf("alignment/salmon/DeSeq2/heatmap_rld.pdf")
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
pdf("alignment/salmon/DeSeq2/PCA_vst.pdf")
plotPCA(vst,intgroup=c("Isolate", "Plant.Line", "Timepoint"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/DeSeq2/PCA_rld.pdf")
plotPCA(vst,intgroup=c("Isolate", "Plant.Line", "Timepoint"))
dev.off()


#Plot using rlog transformation, showing sample names:

data <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
 geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
 coord_fixed()

ggsave("alignment/salmon/DeSeq2/PCA_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)
```


#       DESeq2 analysis

# Timepoint analysis

```bash
# set the significance level for BH adjustment	    
alpha <- 0.05

# Timepoint
#dds$Timepoint[is.na(dds$Timepoint)] <- 0
dds$Timepoint <- as.factor(dds$Timepoint)
design=~Timepoint
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
res_Timepoint <-  results(dds,alpha=alpha,parallel=T,contrast=c("Timepoint","12 hours","48 hours"))
res.merged <- data.frame()
res.merged <- rownames_to_column(as.data.frame(res_Timepoint))
colnames(res.merged)[1] <- "GENE"
write.table(res.merged,"alignment/salmon/DeSeq2/res_12_vs_48.tlb",sep="\t",quote=F,na="",row.names=F)

#res.merged <- left_join(rownames_to_column(as.data.frame(res_Timepoint)),annotations,by=c("rowname"="GENE"))
```

"Mycelium - Mycelium","Emily - 12 hours"

```R
alpha <- 0.05
design=~Group
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
res= results(dds, alpha=alpha,contrast=c("Group","Mycelium - Mycelium", "Emily - 12 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/DeSeq2/Mycelium_vs_Emily12h.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/DeSeq2/Mycelium_vs_Emily12h_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/DeSeq2/Mycelium_vs_Emily12h_down.txt",sep="\t",na="",quote=F)
```

"Mycelium - Mycelium","Emily - 48 hours"

```R
alpha <- 0.05
design=~Group
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
res= results(dds, alpha=alpha,contrast=c("Group","Mycelium - Mycelium", "Emily - 48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/DeSeq2/Mycelium_vs_Emily48h.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/DeSeq2/Mycelium_vs_Emily48h_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/DeSeq2/Mycelium_vs_Emily48h_down.txt",sep="\t",na="",quote=F)
```


"Emily - 12 hours","Emily - 48 hours"

```R
alpha <- 0.05
design=~Group
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
#dds$Timepoint[is.na(dds$Timepoint)] <- 0
res= results(dds, alpha=alpha,contrast=c("Group","Emily - 12 hours","Emily - 48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
#sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
#sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/salmon/DeSeq2/Emily12h_vs_Emily48h.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/DeSeq2/Emily12h_vs_Emily48h_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/DeSeq2/Emily12h_vs_Emily48h_down.txt",sep="\t",na="",quote=F)
```


"Mycelium - Mycelium","Fenella - 12 hours"

```R
alpha <- 0.05
design=~Group
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
res= results(dds, alpha=alpha,contrast=c("Group","Mycelium - Mycelium", "Fenella - 12 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/DeSeq2/Mycelium_vs_Fenella12h.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/DeSeq2/Mycelium_vs_Fenella12h_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/DeSeq2/Mycelium_vs_Fenella12h_down.txt",sep="\t",na="",quote=F)
```

"Mycelium - Mycelium","Fenella - 48 hours"

```R
alpha <- 0.05
design=~Group
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
res= results(dds, alpha=alpha,contrast=c("Group","Mycelium - Mycelium", "Fenella - 48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]

write.table(sig.res,"alignment/salmon/DeSeq2/Mycelium_vs_Fenella48h.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/DeSeq2/Mycelium_vs_Fenella48h_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/DeSeq2/Mycelium_vs_Fenella48h_down.txt",sep="\t",na="",quote=F)
```


"Fenella - 12 hours","Fenella - 48 hours"

```R
alpha <- 0.05
design=~Group
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
#dds$Timepoint[is.na(dds$Timepoint)] <- 0
res= results(dds, alpha=alpha,contrast=c("Group","Fenella - 12 hours","Fenella - 48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
#sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
#sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/salmon/DeSeq2/Fenella12h_vs_Fenella48h.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/DeSeq2/Fenella12h_vs_Fenella48h_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/DeSeq2/Fenella12h_vs_Fenella48h_down.txt",sep="\t",na="",quote=F)
```

# Cultivar Analysis

"Emily - 12 hours","Fenella - 12 hours"

```R
alpha <- 0.05
design=~Group
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
#dds$Timepoint[is.na(dds$Timepoint)] <- 0
res= results(dds, alpha=alpha,contrast=c("Group","Emily - 12 hours","Fenella - 12 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
#sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
#sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/salmon/DeSeq2/Emily12h_vs_Fenella12h.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/DeSeq2/Emily12h_vs_Fenella12h_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/DeSeq2/Emily12h_vs_Fenella12h_down.txt",sep="\t",na="",quote=F)
```

"Emily - 12 hours","Fenella - 48 hours"

```R
alpha <- 0.05
design=~Group
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
#dds$Timepoint[is.na(dds$Timepoint)] <- 0
res= results(dds, alpha=alpha,contrast=c("Group","Emily - 12 hours","Fenella - 48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
#sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
#sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/salmon/DeSeq2/Emily12h_vs_Fenella48h.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/DeSeq2/Emily12h_vs_Fenella48h_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/DeSeq2/Emily12h_vs_Fenella48h_down.txt",sep="\t",na="",quote=F)
```


"Fenella - 12 hours","Emily - 48 hours"

```R
alpha <- 0.05
design=~Group
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
#dds$Timepoint[is.na(dds$Timepoint)] <- 0
res= results(dds, alpha=alpha,contrast=c("Group","Fenella - 12 hours","Emily - 48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
#sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
#sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/salmon/DeSeq2/Fenella12h_vs_Emily48h.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/DeSeq2/Fenella12h_vs_Emily48h_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/DeSeq2/Fenella12h_vs_Emily48h_down.txt",sep="\t",na="",quote=F)
```


"Emily - 48 hours","Fenella - 48 hours"

```R
alpha <- 0.05
design=~Group
design(dds) <- design
dds <- DESeq(dds,parallel=T)
disp <- dispersions(dds)
dds$Group <- as.factor(dds$Group)
#dds$Timepoint[is.na(dds$Timepoint)] <- 0
res= results(dds, alpha=alpha,contrast=c("Group","Emily - 48 hours","Fenella - 48 hours"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
#Settings used: upregulated: min. 2x fold change, ie. log2foldchange min 1.
#               downregulated: min. 0.5x fold change, ie. log2foldchange max -1.
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# No threshold
#sig.res.upregulated2 <- sig.res[sig.res$log2FoldChange >0, ]
#sig.res.downregulated2 <- sig.res[sig.res$log2FoldChange <0, ]

write.table(sig.res,"alignment/salmon/DeSeq2/Emily48h_vs_Fenella48h.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/DeSeq2/Emily48h_vs_Fenella48h_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/DeSeq2/Emily48h_vs_Fenella48h_down.txt",sep="\t",na="",quote=F)
```


Make a table of raw counts, normalised counts and fpkm values:

```R
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Group)
write.table(raw_counts,"alignment/salmon/DeSeq2/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Group)
write.table(norm_counts,"alignment/salmon/DeSeq2/normalised_counts.txt",sep="\t",na="",quote=F)

#library(Biostrings)
#library(naturalsort)
#mygenes <- #readDNAStringSet("gene_pred/final_ncbi/P.cactorum/414_v2/final_ncbi/414_v2_genes_incl_ORFeffectors_renamed.cdna.fasta")
#t1 <- counts(dds)
#t1 <- mygenes[rownames(t1)]
#rowRanges(dds) <- GRanges(t1@ranges@NAMES,t1@ranges)


# robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/salmon/DeSeq2/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Group)
write.table(fpkm_counts,"alignment/salmon/DeSeq2/fpkm_counts.txt",sep="\t",na="",quote=F)
```
