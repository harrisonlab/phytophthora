
# 1. Alignment of Pcac raw reads vs the 36-14 genome

Alignment of reads from a single run:

```bash
  # Reference=$(ls ../../../../home/groups/harrisonlab/project_files/idris/repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w '62471')
  Reference=$(ls ../../../../home/groups/harrisonlab/project_files/idris/repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'R36_14')
  for StrainPath in $(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_dna/paired/P.*/* | grep -w -v -e '10300' -e '415' -e '416' -e 'PC13_15' -e '2003_3' | grep -e 'P.cactorum' -e 'P.idaei'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
    R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
    echo $F_Read
    echo $R_Read
    # OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_62471
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_R36_14
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
  done
```

Alignment of reads from multiple sequencing runs:

For isolates with two runs of data:

```bash
  # Reference=$(ls ../../../../home/groups/harrisonlab/project_files/idris/repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w '62471')
  Reference=$(ls ../../../../home/groups/harrisonlab/project_files/idris/repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'R36_14')
  for StrainPath in $(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_dna/paired/P.*/* | grep -e 'P.cactorum' -e 'P.idaei' | grep -w -e '2003_3' -e '415' -e '416' -e 'PC13_15'); do
    echo $StrainPath
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    F1_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1);
    R1_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1);
    F2_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n2 | tail -n1);
    R2_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n2 | tail -n1);
    echo $F1_Read
    echo $R1_Read
    echo $F2_Read
    echo $R2_Read
    # OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_62471
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_R36_14
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie_2lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $OutDir
  done
```

for isolates with three runs of data:

```bash
  # Reference=$(ls ../../../../home/groups/harrisonlab/project_files/idris/repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w '62471')
  Reference=$(ls ../../../../home/groups/harrisonlab/project_files/idris/repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'R36_14')
  # for StrainPath in $(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_dna/paired/P.*/* | grep -e 'P.cactorum' -e 'P.idaei' | grep -w -e '10300'); do
  for StrainPath in $(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_dna/paired/P.*/* | grep -e 'P.cactorum' -e 'P.idaei' | grep -w -e '404' -e '414'); do
    echo $StrainPath
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    F1_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1);
    R1_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1);
    F2_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n2 | tail -n1);
    R2_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n2 | tail -n1);
    F3_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n3 | tail -n1);
    R3_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n3 | tail -n1);
    echo $F1_Read
    echo $R1_Read
    echo $F2_Read
    echo $R2_Read
    echo $F3_Read
    echo $R3_Read
    # OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_62471
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_R36_14
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie_3lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $F3_Read $R3_Read $OutDir
  done
```


# 2. Pre SNP calling cleanup


## 2.1 Rename input mapping files in each folder by prefixing with the strain ID

```bash
  cd /data/scratch/armita/idris
  # for File in $(ls analysis/genome_alignment/bowtie/*/*/vs_62471/62471_contigs*sorted.bam | grep -e '404' -e '414'); do
  for File in $(ls analysis/genome_alignment/bowtie/*/*/vs_R36_14/R36_14_contigs*sorted.bam); do
    Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/popgen/$Organism/$Strain
    CurDir=$PWD
    mkdir -p $OutDir
    cd $OutDir
    # cp -s $CurDir/$File "$Strain"_vs_62471_aligned_sorted.bam
    cp -s $CurDir/$File "$Strain"_vs_R36_14_aligned_sorted.bam
    cd $CurDir
  done
```

## 2.2 Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used:
qsub $ProgDir/sub_pre_snp_calling.sh <INPUT SAM FILE> <SAMPLE_ID>

```bash
for Sam in $(ls analysis/popgen/*/*/*_vs_R36_14_aligned_sorted.bam); do
Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Sam | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
qsub $ProgDir/sub_pre_snp_calling.sh $Sam $Strain
done
```

# 3. Run SNP calling

#Runs a SNP calling script from Maria in order to be able to draw up a phylogeny
To change in each analysis:

<!-- ```bash
input=/home/groups/harrisonlab/project_files/phytophthora_fragariae/analysis/genome_alignment/bowtie
reference=repeat_masked/P.fragariae/Bc16/filtered_contigs_repmask/95m_contigs_unmasked.fa

filename=$(basename "$reference")
output="${filename%.*}.dict"
``` -->

##Prepare genome reference indexes required by GATK

```bash
# Reference=$(ls ../../../../home/groups/harrisonlab/project_files/idris/repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'R36_14')
Reference=$(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'R36_14')
OutName=$(echo $Reference | sed 's/.fa/.dict/g')
OutDir=$(dirname $Reference)
mkdir -p $OutDir
ProgDir=/home/sobczm/bin/picard-tools-2.5.0
java -jar $ProgDir/picard.jar CreateSequenceDictionary R=$Reference O=$OutName
samtools faidx $Reference
```

###Copy index file to same folder as BAM alignments
<!--
```bash
Reference=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa)
for AlignDir in $(ls -d analysis/popgen/P.*/*/); do
    Index="$Reference".dict
    Directory=analysis/genome_alignment/bowtie/*/$Strain/vs_62471/
    cp $Index $AlignDir/.
done
``` -->

Move to the directory where the output of SNP calling should be placed. Then
Start SNP calling with GATK.
The submission script required need to be custom-prepared for each analysis,
depending on what samples are being analysed. See inside the submission script
below:

```bash
# Isolate=62471
Isolate='R36_14'
Reference=$(ls /home/groups/harrisonlab/project_files/idris/repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w "${Isolate}")
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling
mkdir -p $OutDir
cd $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
qsub $ProgDir/sub_SNP_calling_multithreaded2.sh $Reference $Isolate
cd $CurDir
```

## Filter SNPs based on this region being present in all isolates

Only retain biallelic high-quality SNPS with no missing data (for any individual) for genetic analyses below (in some cases, may allow some missing data in order to retain more SNPs, or first remove poorly sequenced individuals with too much missing data and then filter the SNPs).

```bash
Vcf=$(ls analysis/popgen/SNP_calling/vs_R36_14/R36_14_contigs_softmasked_repeatmasker_TPSI_appended.vcf)
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
# mq=40
# qual=30
# dp=10
# gq=30
# na=0.95
# removeindel=Y
# $VcfLib/vcffilter -f "QUAL > $qual & MQ > $mq"
# $vcftools/vcftools --vcf temp.vcf --max-missing $na --remove-indels --recode --out ${filename%.vcf}_filtered
qsub $ProgDir/sub_vcf_parser.sh $Vcf 40 30 10 30 1 Y
```

```bash
# mv 62471_contigs_softmasked_repeatmasker_TPSI_appended_filtered.vcf analysis/popgen/SNP_calling/vs_62471/62471_contigs_softmasked_repeatmasker_TPSI_appended_filtered.vcf
mv R36_14_contigs_softmasked_repeatmasker_TPSI_appended_filtered.vcf analysis/popgen/SNP_calling/vs_R36_14/R36_14_contigs_softmasked_repeatmasker_TPSI_appended_filtered.vcf
```

```
After filtering, kept 20 out of 20 Individuals
Outputting VCF file...
After filtering, kept 65666 out of a possible 143123 Sites
Run Time = 9.00 seconds
```

## Remove sequencing errors from vcf files:

```bash
Isoalte="R36_14"
Vcf=$(ls analysis/popgen/SNP_calling/vs_${Isoalte}/${Isoalte}_contigs_softmasked_repeatmasker_TPSI_appended_filtered.vcf)
OutDir=$(dirname $Vcf)
Errors=$OutDir/${Isoalte}_error_SNPs.tsv
FilteredVcf=$OutDir/${Isoalte}_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
$ProgDir/flag_error_SNPs.py --ploidy 'diploid' --inp_vcf $Vcf --ref_isolate ${Isoalte} --errors $Errors --filtered $FilteredVcf
# echo "The number of probable errors from homozygous SNPs being called from reference illumina reads vs the reference assembly is:"
# cat $Errors
# cat $Errors | wc -l
echo "the number of SNPs before filtering errors:"
cat $Vcf | grep -v '#' | wc -l
echo "the number of SNPs after filtering errors:"
cat $FilteredVcf | grep -v '#' | wc -l
echo "These have been removed from the vcf file"
```

<!--
```
This is the output errors file but it can't be right.
```
```
contig_3	283491
contig_4	77401
contig_13	51801
contig_27	55400
contig_51	29459
contig_53	2727
contig_56	29940
contig_72	9040
contig_77	35680
contig_80	82020
``` -->

<!--
In some organisms, may want to thin (subsample) SNPs in high linkage diseqilibrium down to
1 SNP  per e.g. 10 kbp just for the population structure analyses.
```bash
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $input_vcf --thin 10000 --recode --out ${input_vcf%.vcf}_thinned
```
-->

## Collect VCF stats

General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  Isolate="R36_14"
  VcfTools=/home/sobczm/bin/vcftools/bin
  export PERL5LIB="$VcfTools:$PERL5LIB"
  # Vcf=$(ls analysis/popgen/SNP_calling/*_filtered_no_errors.vcf)
  # Stats=$(echo $Vcf | sed 's/.vcf/.stat/g')
  # perl $VcfTools/vcf-stats $Vcf > $Stats
  VcfFiltered=$(ls analysis/popgen/SNP_calling/vs_${Isolate}/${Isolate}_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```

Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
  Isolate="R36_14"
  for Vcf in $(ls analysis/popgen/SNP_calling/vs_${Isolate}/*filtered_no_errors.vcf); do
      ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
      $ProgDir/similarity_percentage.py $Vcf
  done
```

# Visualise the output as heatmap and clustering dendrogram
```bash
for Log in $(ls analysis/popgen/SNP_calling/vs_${Isolate}/*distance.log); do
  ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
  Rscript --vanilla $ProgDir/distance_matrix.R $Log
  mv Rplots.pdf analysis/popgen/SNP_calling/.
done
```


## Carry out PCA and plot the results

This step could not be carried out due to problems installing dependancies
<!--
```bash
for Vcf in $(ls analysis/popgen/SNP_calling/vs_62471/*filtered_no_errors.vcf); do
    echo $Vcf
    ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
    # Out=$(basename $Vcf)
    Out=analysis/popgen/SNP_calling
    echo $Out
    /home/deakig/R3.4/bin/Rscript --vanilla $ProgDir/pca.R $Vcf $Out/PCA.pdf
done
``` -->



## Calculate a NJ tree

These commands didnt work as P. idaei is too distant for sufficient sites to be shared
between isolates
<!--
based on all the SNPs. Outputs a basic display of the tree, plus a Newick file to be used for displaying the tree in FigTree and beautifying it.

cut down the vcf to remove P. idaei information.


```bash
for Vcf in $(ls analysis/popgen/SNP_calling/*_filtered_no_errors.vcf); do
echo $Vcf
Out=$(echo $Vcf | sed 's/.vcf/_Pcac_only.vcf/g')
ExludeList="371 SCRP370 SCRP376"
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExludeList > $Out
done
```


Remove all missing data for nj tree construction

```bash
  for Vcf in $(ls analysis/popgen/SNP_calling/*_Pcac_only.vcf); do
    echo $Vcf
    Out=$(basename $Vcf .vcf)
    echo $Out
    VcfTools=/home/sobczm/bin/vcftools/bin
    $VcfTools/vcftools --vcf $Vcf --mac 1 --max-missing 1.0 --recode --out analysis/popgen/SNP_calling/"$Out"_no_missing
  done
```

```bash
for Vcf in $(ls analysis/popgen/SNP_calling/*_Pcac_only_no_missing.recode.vcf); do
    echo $Vcf
    Ploidy=2
    # ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
    scripts=/home/adamst/git_repos/scripts/popgen/snp
    $ProgDir/nj_tree.sh $Vcf $Ploidy
    mv Rplots.pdf analysis/popgen/SNP_calling/NJ_tree.pdf
done
ls $PWD/analysis/popgen/SNP_calling/*.nwk
```

The Newick file was downloaded to my local machine imported
into geneious to be re-rooted, before exporting for reading
into ggtree

```R

```r
setwd("/Users/armita/Downloads/Pc/alignments/SNP")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)

tree <- read.tree("414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors_Pcac_only_no_missing.recode_haploidised_nj_geneious.nwk")



mydata <- read.csv("/Users/armita/Downloads/Pc/alignments/SNP/traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$label
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]


t <- ggtree(tree) # Core tree
t <- t + geom_treescale(offset=-0.5, fontsize = 3) # Add scalebar

# labels by values in another df
t <- t %<+% mydata
tips <- data.frame(t$data)
tips$label <- tips$Isolate
t <- t + geom_tiplab(data=tips, size=3, hjust=0, offset = +0.09)
tips2 <- data.frame(t$data)
tips2$label <- tips2$Host
t <- t + geom_tiplab(data=tips2, size=3, hjust=0, offset = +0, fontface = "italic")

# Format nodes by values
nodes <- data.frame(t$data)
#nodes <- nodes[!nodes$isTip,]
nodes$label <- as.numeric(nodes[nodes$label,])
as.numeric(nodes$label)
#nodes$label[nodes$label < 0.80] <- ''
nodes$support[nodes$isTip] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label > 80)] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label < 80)] <- 'unsupported'
nodes$support[(!nodes$isTip) & (nodes$label == '')] <- 'supported'
t <- t + aes(linetype=nodes$support)
nodes$label[nodes$label > 80] <- ''
t <- t + geom_nodelab(data=nodes, size=2, hjust=+0.05) # colours as defined by col2rgb

# Add in a further set of labels
# tree_mod <- tree
# tree_mod$tip.label <- mydata$Source
# t <- t + geom_tiplab(data=tree_mod, aes(label=label), size=2, offset = +1)

# Annotate a clade with a bar line
#t <- t + geom_cladelabel(node=24, label='P. i', align=T, colour='black', offset=+1, fontface = "italic")
#t <- t + geom_cladelabel(node=26, label='P. c', align=T, colour='black', offset=+1, fontface = "italic", fontface = "italic")

# Save as PDF and force a 'huge' size plot
ggsave("Pcac_SNP_phylogeny.pdf", width =20, height = 25, units = "cm", limitsize = FALSE)
```
 -->

# Identify SNPs in gene models:

Create custom SnpEff genome database

```bash
SnpEff=/home/sobczm/bin/snpEff
nano $SnpEff/snpEff.config
```


Add the following lines to the section with databases:

```
#---
# EMR Databases
#----
# Fus2 genome
Fus2v1.0.genome : Fus2
# Bc16 genome
Bc16v1.0.genome: BC-16
# P414 genome
P414v1.0.genome: 414
# 62471 genome
62471v1.0.genome: 62471
# R36_14 genome
R36_14v1.0.genome: R36_14
# SCRP371 genome
SCRP371v1.0.genome: SCRP371
```

Collect input files

```bash
Organism="P.cactorum"
Strain="R36_14"
DbName="R36_14v1.0"
# Strain="62471"
# DbName="62471v1.0"
ProjDir=/home/groups/harrisonlab/project_files/idris
Reference=$(ls $ProjDir/repeat_masked/${Organism}/${Strain}/filtered_contigs_repmask/${Strain}_contigs_unmasked.fa)
Gff=$(ls $ProjDir/gene_pred/final_incl_ORF/${Organism}/${Strain}/final_genes_genes_incl_ORFeffectors_renamed.gff3)
SnpEff=/home/sobczm/bin/snpEff
mkdir $SnpEff/data/${DbName}
cp $Reference $SnpEff/data/${DbName}/sequences.fa
cp $Gff $SnpEff/data/${DbName}/genes.gff

#Build database using GFF3 annotation
java -jar $SnpEff/snpEff.jar build -gff3 -v ${DbName}
```


## Annotate VCF files
```bash
Organism="P.cactorum"
Strain="R36_14"
DbName="R36_14v1.0"
# Strain="62471"
# DbName="62471v1.0"
# CurDir=/data/scratch/armita/idris
CurDir=/home/groups/harrisonlab/project_files/idris
cd $CurDir
for a in $(ls analysis/popgen/SNP_calling/vs_${Strain}/*_filtered_no_errors.vcf); do
    echo $a
    filename=$(basename "$a")
    Prefix=${filename%.vcf}
    OutDir=$(ls -d analysis/popgen/SNP_calling)
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 ${DbName} $a > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
    # mv 414_v2_contigs_unmasked_filtered* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf
    AllSnps=$(cat $OutDir/"$Prefix"_annotated.vcf | grep -v '#' | wc -l)
    GeneSnps=$(cat $OutDir/"$Prefix"_gene.vcf | grep -v '#' | wc -l)
    CdsSnps=$(cat $OutDir/"$Prefix"_coding.vcf | grep -v '#' | wc -l)
    NonsynSnps=$(cat $OutDir/"$Prefix"_nonsyn.vcf | grep -v '#' | wc -l)
    SynSnps=$(cat $OutDir/"$Prefix"_syn.vcf | grep -v '#' | wc -l)
done
```

# 3.0 Comparisons of groups to reference 62471 genome

# 3.2 P. cactorum ex. apple vs 62471

```bash
  Organism="P.cactorum"
  Strain="R36_14"
  DbName="R36_14v1.0"
  # Strain="62471"
  # DbName="62471v1.0"

  Prefix=Pc_apple_vs_${Strain}
  OutDir=analysis/popgen/SNP_calling/$Prefix
  mkdir -p $OutDir

  Vcf=$(ls analysis/popgen/SNP_calling/vs_${Strain}/*_filtered_no_errors.vcf)
  ExcludeList="12420 15_13 15_7 2003_3 4032 4040 404 415 416 371 SCRP370 SCRP376 414 4040 11-40 17-21 P421"
  VcfLib=/home/sobczm/bin/vcflib/bin
  $VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf

  VcfTools=/home/sobczm/bin/vcftools/bin
  $VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --non-ref-ac-any 1 --max-missing 0.95 --remove-indels --recode --out $OutDir/"$Prefix"_filtered_no_indels

  for Vcf in $(ls $OutDir/"$Prefix"_filtered_no_indels.recode.vcf); do
      echo $Vcf
      ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/summary_stats
      $ProgDir/annotate_snps_genome.sh $Vcf ${DbName}

      filename=$(basename "$Vcf")
      Prefix=$(echo $filename | sed 's/.vcf//g')
      SnpEff=/home/sobczm/bin/snpEff
      java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 ${DbName} $Vcf > $OutDir/"$Prefix"_annotated.vcf
      mv snpEff_genes.txt $OutDir/snpEff_genes_$Prefix.txt
      mv snpEff_summary.html $OutDir/snpEff_summary_$Prefix.html

      #Create subsamples of SNPs containing those in a given category
      #genic (includes 5', 3' UTRs)
      java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
      #coding
      java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/${filename%.vcf}_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
      #non-synonymous
      java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
      #synonymous
      java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
      #Four-fold degenrate sites (output file suffix: 4fd)
      ProgDir=/home/sobczm/bin/popgen/summary_stats
      python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf
  done
```

```
WARNINGS: Some warning were detected
Warning type	Number of warnings
WARNING_TRANSCRIPT_INCOMPLETE	2
WARNING_TRANSCRIPT_NO_STOP_CODON	1
```
