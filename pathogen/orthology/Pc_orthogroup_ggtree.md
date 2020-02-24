
<!-- ```bash
echo "((SCRP370:1.0,(371:1.0,SCRP376:1.0)0.9:0.08245622453296164)1:1.6434064099824397,((17-21:0,(R36_14:1.0,(P295:1.0,62471:1.0)1:0.29696569109508353)1:0.32373562653903853)1:0.15949253028150512,(10300:0,(((2003_3:1.0,4032:1.0)0.54:0.03507127382577924,(P421:1.0,(414:1.0,PC13_15:1.0)0.48:0.027700345429534146)0.48:0.019840795963525615)0.61:0.03554532977083369,((4040:1.0,416:1.0)0.77:0.05715039539753075,(11-40:1.0,((12420:1.0,15_13:1.0)0.69:0.04700123987863236,(404:1.0,(15_7:1.0,415:1.0)0.45:0.015690963218303544)0.47:0.021387852567743337)0.52:0.02567263702771161)0.37:0.0046454517207186186)0.39:0.007662872745568983)0.88:0.07713662906987118)1:0.5668858014697449):1.6434064099824397);" > tmp.tre
``` -->

```bash
echo "((((((4040:1.0,416:1.0)0.7:0.06531593380614709,10300:1.0)0.41:0.016352210995901117,11-40:1.0)0.39:0.008636177951949264,(((P421_v2:1.0,PC13_15:1.0)0.55:0.03690376847131205,(414:1.0,(404:1.0,(4032:1.0,2003_3:1.0)0.69:0.06080125345162024)0.4:0.0109659470879242)0.39:0.008636177951949264)0.46:0.021839095463741742,(((12420:1.0,15_7:1.0)0.42:0.01620261565102954,15_13:1.0)0.61:0.048706275070770744,415:1.0)0.29:0.0)0.5:0.030531643166989753)1:0.6671152295660177,(((P295:1.0,62471:1.0)1:0.2085908267031984,R36_14:1.0)1:0.30308722399627275,17-21:1.0)0.62:0.05430028347563276):2.068355017454495,((371:1.0,SCRP376:1.0)0.82:0.09167849201605893,SCRP370:1.0)1:2.068355017454495);" > /Users/armita/OneDrive\ -\ University\ of\ Greenwich/Armitage/Projects/Phytophthora\ EMR/Pcac_phylogeny.consensus.scored.tre.geneious_mod_gain-loss.newick

```


```r
setwd("/Users/armita/OneDrive\ -\ University\ of\ Greenwich/Armitage/Projects/Phytophthora\ EMR")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)


tree <- read.tree("Pcac_phylogeny.consensus.scored.tre.geneious_mod_gain-loss.newick")
tree$edge.length[tree$edge.length == 1] <- 0
mydata <- read.csv("/Users/armita/Downloads/Pc/alignments/ASTRAL/traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$label
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]


t <- ggtree(tree) # Core tree


edge=data.frame(tree$edge, edge_num=1:nrow(tree$edge))
colnames(edge)=c("parent", "node", "edge_num")

# Identify node numbers by making a separate plot.
t2 <- t %<+% edge + geom_label(aes(x=branch, label=node))
t2

edge$gainloss <- n/a
# edge$gainloss[edge$node == 22] <- '1690 Pc gain or Pi loss\n1124 Pi gain or Pc loss'
edge$gainloss[edge$node == 23] <- '1124 Pc gain or Pi loss'
edge$gainloss[edge$node == 40] <- '1690 Pi gain or Pc loss'
edge$gainloss[edge$node == 24] <- '+45\n-21'
edge$gainloss[edge$node == 37] <- '+72\n-9'
edge$gainloss[edge$node == 38] <- '+47\n-19'
edge$gainloss[edge$node == 18] <- '+111\n-54'


t <- t %<+% edge + geom_label(aes(x=branch, label=gainloss))
t3 <- t %>% collapse(node=40)
t3 <- t3 %>% collapse(node=24)
t3 <- t3 %>% collapse(node=38)
t3

# Save as PDF and force a 'huge' size plot
ggsave("Pcac_gain-loss.pdf", plot = t3, width =20, height = 15, units = "cm", limitsize = FALSE)
````
