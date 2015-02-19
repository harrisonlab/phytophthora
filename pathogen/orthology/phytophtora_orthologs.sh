# Phytophthora orthologous genes
ls -d gene_pred/augustus/P*/*
# gene_pred/augustus/P.cactorum/10300  gene_pred/augustus/P.cactorum/414       gene_pred/augustus/P.ideai/371
# gene_pred/augustus/P.cactorum/404    gene_pred/augustus/P.fragariae/JHVZ02
# gene_pred/augustus/P.cactorum/411    gene_pred/augustus/P.fragariae/SCRP245

set -- 10300 404 414 JHVZ02 371
	for a; do 
	shift
	for b; do
	printf "%s - %s\n" "$a" "$b"
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/sub_inparanoid.sh gene_pred/augustus/*/$a/"$a"_augustus_preds.aa gene_pred/augustus/*/$b/"$b"_augustus_preds.aa gene_pred/augustus/*/$a/"$a"_augustus_preds.gtf gene_pred/augustus/*/$b/"$b"_augustus_preds.gtf
	done 
done

cat analysis/inparanoid/*/*_seqs.txt | cut -f1 | sort | uniq > analysis/inparanoid/summary_tables/all_genes.txt
for FILEZ in $(ls analysis/inparanoid/*/sqltable.*); do 
	FILE_LIST="$FILE_LIST $FILEZ"
done


mkdir analysis/inparanoid/summary_tables
cat analysis/inparanoid/*/*_seqs.txt | cut -f1 | sort | uniq > analysis/inparanoid/summary_tables/all_genes.txt
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_orthology_tab.pl analysis/inparanoid/summary_tables/all_genes.txt analysis/inparanoid/*/sqltable.* > analysis/inparanoid/summary_tables/phytophthora_orthology_tab.csv
mkdir -p analysis/inparanoid/summary_tables/gene_orthologs
cd analysis/inparanoid/summary_tables/gene_orthologs
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/analyse_orthology_tab.pl ../../../../analysis/inparanoid/summary_tables/phytophthora_orthology_tab.csv -file ../../../../analysis/inparanoid/summary_tables/all_genes.txt -deep -print_once
cd ../../../../
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_final_ortho_tab.pl analysis/inparanoid/summary_tables/gene_orthologs/ 10300 404 414 371 JHVZ02 > analysis/inparanoid/summary_tables/phytophthora_orthology.csv
cat analysis/inparanoid/summary_tables/phytophthora_orthology.csv | sed 's/\+/1/g' | sed 's/\-/0/g' > analysis/inparanoid/summary_tables/phytophthora_orthology2.csv


# Making a 5-way Venn Diagram:
# Open R
R
# These commabds were executed in R

library(VennDiagram, lib.loc="/home/armita/R-packages/")
data1 <- read.delim(file="analysis/inparanoid/summary_tables/phytophthora_orthology2.csv", header=T, sep="\t")
attach(data1)

colnames(data1) <- c("ortho_group", "gene_name", 'X10300', 'X404', 'X414', 'X371', 'JHVZ02')

area1=sum(data1$X10300)
area2=sum(data1$X404)
area3=sum(data1$X414)
area4=sum(data1$X371)
area5=sum(data1$JHVZ02)
label1 <- paste('10300 (', area1, ')', sep ="")
label2 <- paste('404 (', area2, ')', sep ="")
label3 <- paste('414 (', area3, ')', sep ="")
label4 <- paste('371 (', area4, ')', sep ="")
label5 <- paste('JHVZ02 (', area5, ')', sep ="")

n12=nrow(subset(data1, X10300==1 & X404==1))
n13=nrow(subset(data1, X10300==1 & X414==1))
n14=nrow(subset(data1, X10300==1 & X371==1))
n15=nrow(subset(data1, X10300==1 & JHVZ02==1))
n23=nrow(subset(data1, X404==1 & X414==1))
n24=nrow(subset(data1, X404==1 & X371==1))
n25=nrow(subset(data1, X404==1 & JHVZ02==1))
n34=nrow(subset(data1, X414==1 & X371==1))
n35=nrow(subset(data1, X414==1 & JHVZ02==1))
n45=nrow(subset(data1, X371==1 & JHVZ02==1))
n123=nrow(subset(data1, X10300==1 & X404==1 & X414==1))
n124=nrow(subset(data1, X10300==1 & X404==1 & X371==1))
n125=nrow(subset(data1, X10300==1 & X404==1 & JHVZ02==1))
n134=nrow(subset(data1, X10300==1 & X414==1 & X371==1))
n135=nrow(subset(data1, X10300==1 & X414==1 & JHVZ02==1))
n145=nrow(subset(data1, X10300==1 & X371==1 & JHVZ02==1))
n234=nrow(subset(data1, X404==1 & X414==1 & X371==1))
n235=nrow(subset(data1, X404==1 & X414==1 & JHVZ02==1))
n245=nrow(subset(data1, X404==1 & X371==1 & JHVZ02==1))
n345=nrow(subset(data1, X414==1 & X371==1 & JHVZ02==1))
n1234=nrow(subset(data1, X10300==1 & X404==1 & X414==1 & X371==1))
n1235=nrow(subset(data1, X10300==1 & X404==1 & X414==1 & JHVZ02==1))
n1245=nrow(subset(data1, X10300==1 & X404==1 & X371==1 & JHVZ02==1))
n1345=nrow(subset(data1, X10300==1 & X414==1 & X371==1 & JHVZ02==1))
n2345=nrow(subset(data1, X404==1 & X414==1 & X371==1 & JHVZ02==1))
n12345=nrow(subset(data1, X10300==1 & X404==1 & X414==1 & X371==1 & JHVZ02==1))
pdf('analysis/inparanoid/summary_tables/phytoph_ortholog_venn.pdf')
draw.quintuple.venn(area1, area2, area3, area4, area5, n12, n13, n14, n15, n23, n24, n25, n34, n35, n45, n123, n124, n125, n134, n135, n145, n234, n235, n245, n345, n1234,
n1235, n1245, n1345, n2345, n12345, category = c(label1, label2, label3, label4, label5), lwd = rep(2, 5),
	lty = rep("solid", 5), col = rep("black", 5), fill = NULL, alpha = rep(0.5, 5), label.col = rep("black", 31), cex = rep(1, 31), fontface = rep("plain", 31),
    fontfamily = rep("serif", 31), cat.pos = c(0, 287.5, 215, 145, 70),
    cat.dist = rep(0.2, 5), cat.col = rep("black", 5), cat.cex = rep(1, 5),
    cat.fontface = rep("plain", 5), cat.fontfamily = rep("serif", 5),
    cat.just = rep(list(c(0.5, 0.5)), 5), rotation.degree = 0,
    rotation.centre = c(0.5, 0.5), ind = TRUE, margin = 0.15)
dev.off()
q()

printf "" > analysis/inparanoid/summary_tables/inparalog_groups.txt
for FILEZ in $(ls analysis/inparanoid/summary_tables/gene_orthologs/*); do
THISFILE=$(basename $FILEZ)
printf "$THISFILE\n" >> analysis/inparanoid/summary_tables/inparalog_groups.txt
cat $FILEZ | cut -f1 -d"|" | sort | uniq -c >> analysis/inparanoid/summary_tables/inparalog_groups.txt
done
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep '10300' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep '404' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep '414' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep '371' | wc -l
cat analysis/inparanoid/summary_tables/inparalog_groups.txt | grep -v '.txt' | grep -v '1 ' | grep 'JHVZ02' | wc -l



