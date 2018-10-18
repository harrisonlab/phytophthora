
```bash
echo "((SCRP370:1.0,(371:1.0,SCRP376:1.0)0.9:0.08245622453296164)1:1.6434064099824397,((17-21:0,(R36_14:1.0,(P295:1.0,62471:1.0)1:0.29696569109508353)1:0.32373562653903853)1:0.15949253028150512,(10300:0,(((2003_3:1.0,4032:1.0)0.54:0.03507127382577924,(P421:1.0,(414:1.0,PC13_15:1.0)0.48:0.027700345429534146)0.48:0.019840795963525615)0.61:0.03554532977083369,((4040:1.0,416:1.0)0.77:0.05715039539753075,(11-40:1.0,((12420:1.0,15_13:1.0)0.69:0.04700123987863236,(404:1.0,(15_7:1.0,415:1.0)0.45:0.015690963218303544)0.47:0.021387852567743337)0.52:0.02567263702771161)0.37:0.0046454517207186186)0.39:0.007662872745568983)0.88:0.07713662906987118)1:0.5668858014697449):1.6434064099824397);" > tmp.tre
```


```r
setwd("/Users/armita/Downloads/Pc/alignments/ASTRAL")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)

# tree <- read.tree("Pcac_phylogeny.consensus.scored_geneious2.tre")
tree <- read.tree("tmp.tre")



mydata <- read.csv("/Users/armita/Downloads/Pc/alignments/ASTRAL/traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$label
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]


t <- ggtree(tree) # Core tree
# t <- t + geom_treescale(offset=-0.5, fontsize = 3) # Add scalebar

# labels by values in another df
# t <- t %<+% mydata
# tips <- data.frame(t$data)
# tips$label <- tips$Isolate
# t <- t + geom_tiplab(data=tips, size=3, hjust=0, offset = +0.5)
# tips2 <- data.frame(t$data)
# tips2$label <- tips2$Host
# t <- t + geom_tiplab(data=tips2, size=3, hjust=0, offset = +0, fontface = "italic")

# Format nodes by values
# nodes <- data.frame(t$data)
#nodes <- nodes[!nodes$isTip,]
# nodes$label <- as.numeric(nodes[nodes$label,])
# as.numeric(nodes$label)
#nodes$label[nodes$label < 0.80] <- ''
# nodes$support[nodes$isTip] <- 'supported'
# nodes$support[(!nodes$isTip) & (nodes$label > 0.80)] <- 'supported'
# nodes$support[(!nodes$isTip) & (nodes$label < 0.80)] <- 'unsupported'
# nodes$support[(!nodes$isTip) & (nodes$label == '')] <- 'supported'
# t <- t + aes(linetype=nodes$support)
# nodes$label[nodes$label > 0.80] <- ''
# t <- t + geom_nodelab(data=nodes, size=2, hjust=-0.05) # colours as defined by col2rgb

# Add in a further set of labels
# tree_mod <- tree
# tree_mod$tip.label <- mydata$Source
# t <- t + geom_tiplab(data=tree_mod, aes(label=label), size=2, offset = +1)

# Annotate a clade with a bar line

t <- t + geom_cladelabel(node=23, label='Pi', align=T, colour='black', offset=+0.75, fontface = "italic")
t <- t + geom_cladelabel(node=23, label='R x i\n(3 isolates)', align=T, colour='black', offset=+0.125, fontface = "italic")
t <- t + geom_cladelabel(node=25, label='Pc', align=T, colour='black', offset=+0.75, fontface = "italic")
t <- t + geom_cladelabel(node=30, label='F x a\n(12 isolates)', align=T, colour='black', offset=+0.125, fontface = "italic", fontface = "italic")
t <- t + geom_cladelabel(node=27, label='M x d\n(3 isolates)', align=T, colour='black', offset=+0.125, fontface = "italic", fontface = "italic")
t <- t + geom_cladelabel(node=4, label='F x a*\n17-21', align=T, colour='black', offset=+0.125, fontface = "italic")
t <- t + geom_cladelabel(node=8, label='F x a\n10300', align=T, colour='black', offset=+0.125, fontface = "italic")


# t + geom_taxalink(25, 23, curvature = -0.5, arrow = arrow(),
#        arrow.fill = NULL, size = 0.5)


# t + edgelabels()

edge=data.frame(tree$edge, edge_num=1:nrow(tree$edge))
colnames(edge)=c("parent", "node", "edge_num")
t2 <- t %<+% edge + geom_label(aes(x=branch, label=node))
# t2
edge$gainloss <- n/a
# edge$gainloss[edge$node == 22] <- '1690 Pc gain or Pi loss\n1124 Pi gain or Pc loss'
edge$gainloss[edge$node == 25] <- '1690 Pc gain or Pi loss'
edge$gainloss[edge$node == 23] <- '1124 Pi gain or Pc loss'
edge$gainloss[edge$node == 29] <- '+45\n-21'
edge$gainloss[edge$node == 26] <- '+72\n-9'
edge$gainloss[edge$node == 27] <- '+47\n-19'
edge$gainloss[edge$node == 4] <- '+111\n-54'
edge$gainloss[edge$node == 8] <- '+491\n-264'
edge$gainloss[edge$node == 30] <- '+10\n-1'

t <- t %<+% edge + geom_label(aes(x=branch, label=gainloss))
t3 <- t %>% collapse(node=30)
t3 <- t3 %>% collapse(node=23)
t3 <- t3 %>% collapse(node=27)
t3

# Save as PDF and force a 'huge' size plot
ggsave("Pcac_gain-loss.pdf", width =20, height = 15, units = "cm", limitsize = FALSE)
````
