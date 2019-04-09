# Calculating Neutrality stats

This includes calculating the number of segregating sites, Watterson's Theta, Tajima's D, Fu & Li's F* and Fu & Li'd D


# All genes

## Tabindex the vcf files:


Slim down the vcf to the three populations:

```bash
cd /data/scratch/armita/idris
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="11-40 17-21 10300 LV007"
Prefix=core_isolates
OutDir=analysis/popgen/SNP_calling/summary_stats2/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
mkdir $OutDir/vcf
cp $OutDir/"$Prefix"_filtered.recode.vcf $OutDir/vcf/P414.vcf
bgzip $OutDir/vcf/P414.vcf
tabix -p vcf $OutDir/vcf/P414.vcf.gz

mkdir $OutDir/gff
Gff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.gff3)
cp $Gff $OutDir/gff/P414.gff
```


## Comparison between Pcac crown rot, apple and P. idaei populations

Slim down the vcf to the three populations:

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="11-40 17-21 10300 LV007"
Prefix=core_isolates
OutDir=analysis/popgen/SNP_calling/summary_stats2/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```


## create the directory structure

```bash
Prefix=core_isolates
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats2/$Prefix/neutrality_genes

mkdir -p $WorkDir/all
mkdir -p $WorkDir/gff
mkdir -p $WorkDir/contigs
```

###copy the "gff" folder containing gff files

```bash
CurDir=/data/scratch/armita/idris
cd $CurDir
Prefix=core_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats2/$Prefix/neutrality_genes
Gff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.gff3)
cd $WorkDir/gff
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/summary_stats
$ProgDir/split_gff_contig.sh $CurDir/$Gff
cd $CurDir
```

### Convert vcf files into fasta alignments

The vcf file needs to be converted into fasta alignment files. One file is
produced per reference contig

For the ploidy set <1|2|3>:
* 1 - haploid input
* 2 - diploid input, output as two separate phased haplotypes for each ind.
* 3 - diploid input, output as one sequence with ambiguity codes for each ind.

```bash
CurDir=/data/scratch/armita/idris
cd $CurDir
Prefix=core_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats2/$Prefix/neutrality_genes


Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
Vcf=$(ls analysis/popgen/SNP_calling/summary_stats2/$Prefix/${Prefix}_filtered.recode.vcf)
Ploidy=2

cd $WorkDir/contigs
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/summary_stats
python $ProgDir/vcf_to_fasta.py $CurDir/$Vcf $CurDir/$Reference $Ploidy
cd $CurDir
```


### This folder contains only contig FASTA files
### So create a new "contigs" directory to hold those files:

```bash
CurDir=/data/scratch/armita/idris
Prefix=core_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats2/$Prefix/neutrality_genes

Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
cp $Reference $WorkDir/.
```

###Next step: in the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
CurDir=/data/scratch/armita/idris
Prefix=core_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats2/$Prefix/neutrality_genes
cd $WorkDir/contigs
for File in $(ls *.fasta); do
    Prefix=${File%.fasta}
    mkdir $Prefix
    mv $File $Prefix/.
done
cd $CurDir
```


### Test if all contigs have a matching gff and remove any which do not

```bash
CurDir=/data/scratch/armita/idris
Prefix=core_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats2/$Prefix/neutrality_genes

cd $WorkDir
for File in $(ls $PWD/contigs/*/*.fasta); do
    Prefix=$(basename "$File" .fasta)
    expected_gff="$PWD/gff/${Prefix}.gff"
    if [ ! -f "$expected_gff" ]; then
       rm -rf $(dirname $File)
    fi
done
cd $CurDir
```

```bash
CurDir=/data/scratch/armita/idris
Prefix=core_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats2/$Prefix/neutrality_genes
mkdir $WorkDir/fasta
cd $WorkDir/fasta
for File in $(ls ../contigs/*/*.fasta); do
    cp $File .
done
cd $CurDir
```

```bash
Prefix=core_isolates
cd analysis/popgen/SNP_calling/summary_stats2/core_isolates/neutrality_genes
# ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen/popgenome_scripts
# Rscript --vanilla $ProgDir/calculate_nucleotide_diversity.R
```

Commands documented in calculate_nucleotide_diversity.R were used.

```R
library("PopGenome")
library("ggplot2")

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: number of segregating sites, Watterson's Theta, Tajima's D, Fu & Li's F*, Fu & Li'd D*
Pcac_Fa <- c(
"414_1", "15_13_1", "PC13_15_1", "404_1", "4040_1", "15_7_1", "2003_3_1", "416_1", "415_1", "4032_1", "P421_1", "12420_1",
"414_2", "15_13_2", "PC13_15_2", "404_2", "4040_2", "15_7_2", "2003_3_2", "416_2", "415_2", "4032_2", "P421_2", "12420_2"
)
Pcac_Md <- c(
  "P295_1", "62471_1", "R36_14_1",
  "P295_2", "62471_2", "R36_14_2"
)
Pi_ri <- c(
"371_1", "SCRP376_1", "SCRP370_1",
"371_2", "SCRP376_2", "SCRP370_2"
)
#Need to set argument diploid=TRUE if using diploid genomes in the below command:
populations <- list(Pcac_Fa, Pcac_Md, Pi_ri)
#Number of populations assigned above.
population_no <- length(populations)
population_names <- c("Pcac_Fa", "Pcac_Md", "Pi_ri")

GENOME.class <- readData("fasta/", format="fasta", gffpath = "gff/")
GENOME.class <- set.populations(GENOME.class, populations)
# GENOME.class <- neutrality.stats(GENOME.class)
# GENOME.class@Tajima.D
# GENOME.class <- diversity.stats(GENOME.class)
# GENOME.class@Pi
# GENOME.class@nuc.diversity.within
# GENOME.class@hap.diversity.within
#
# WHOLE <- concatenate.regions(data)
# WHOLE <- neutrality.stats(WHOLE)
# WHOLE@Tajima.D
# WHOLE <- diversity.stats(WHOLE)
# WHOLE@Pi
# WHOLE@nuc.diversity.within
# WHOLE@hap.diversity.within

GENOME.class.split <- splitting.data(GENOME.class, subsites = "gene")
GENOME.class.split <- diversity.stats(GENOME.class.split, pi=TRUE)
#GENOME.class.split@Pi
max(GENOME.class.split@Pi[,1])
max(GENOME.class.split@Pi[,2])
max(GENOME.class.split@Pi[,3])

TajimaD <- GENOME.class.split@Tajima.D

GENOME.class.split.nonsyn <- diversity.stats(GENOME.class.split, pi = TRUE,
    subsites = "nonsyn")
GENOME.class.split.syn <- diversity.stats(GENOME.class.split, pi = TRUE,
    subsites = "syn")
#Interval based
# GENOME.class.slide.nonsyn <- diversity.stats(GENOME.class.slide, pi = TRUE,
#     subsites = "nonsyn")
# GENOME.class.slide.syn <- diversity.stats(GENOME.class.slide, pi = TRUE,
#     subsites = "syn")

#Print output, gene-based
#Plot individual populations, gene-based
#Divide Pi per number of sites in the gene to calculate value per site
Pi_ns <- GENOME.class.split.nonsyn@Pi /
GENOME.class.split.syn@Pi / GENOME.class.split@n.sites
Pi_ns_d <- as.data.frame(Pi_ns)
# Write data
```

```
  > WHOLE@Tajima.D
                pop 1    pop 2    pop 3
  Concatenate 1.37229 1.403231 1.584488

  > WHOLE@Pi
  pop 1 pop 2 pop 3
  Concatenate     0     0     0
  > WHOLE@nuc.diversity.within
     pop 1 pop 2  pop 3
  Concatenate 546.1884  4494 2792.2
  > WHOLE@hap.diversity.within
  pop 1 pop 2 pop 3
  Concatenate     1     1     1
```


# RxLRs

## Comparison between Pcac crown rot, apple and P. idaei populations

Slim down the vcf to the three populations:

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="11-40 17-21 10300 LV007"
Prefix=core_isolates
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```


## create the directory structure

```bash
Prefix=core_isolates
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality_RxLR

mkdir -p $WorkDir/all
mkdir -p $WorkDir/gff
mkdir -p $WorkDir/contigs
```

###copy the "gff" folder containing gff files

```bash
CurDir=/data/scratch/armita/idris
Prefix=core_isolates
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality_RxLR
AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
RxlrHeaders=$WorkDir/RxLR_headers.txt
# cat $AnnotTab | cut -f1,14 | grep 'RxLR' | cut -f1 | cut -f1 -d '.' > $RxlrHeaders
cat $AnnotTab | cut -f1,14 | grep 'RxLR' | cut -f1 > $RxlrHeaders
Gff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.gff3)
RxlrGff=$WorkDir/RxLR.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $RxlrHeaders $Gff effector ID > $RxlrGff

cd $WorkDir/gff
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen
$ProgDir/summary_stats/split_gff_contig.sh $RxlrGff
cd $CurDir
```

### Convert vcf files into fasta alignments

The vcf file needs to be converted into fasta alignment files. One file is
produced per reference contig

For the ploidy set <1|2|3>:
* 1 - haploid input
* 2 - diploid input, output as two separate phased haplotypes for each ind.
* 3 - diploid input, output as one sequence with ambiguity codes for each ind.

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality_RxLR


Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
Vcf=$(ls analysis/popgen/SNP_calling/summary_stats/$Prefix/${Prefix}_filtered.recode.vcf)
Ploidy=2

cd $WorkDir/contigs
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen
python $ProgDir/summary_stats/vcf_to_fasta.py $CurDir/$Vcf $CurDir/$Reference $Ploidy
cd $CurDir
```


### This folder contains only contig FASTA files
### So create a new "contigs" directory to hold those files:

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality_RxLR

Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
cp $Reference $WorkDir/.
```

###Next step: in the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality_RxLR
cd $WorkDir/contigs
for File in $(ls *.fasta); do
    Prefix=${File%.fasta}
    mkdir $Prefix
    mv $File $Prefix/.
done
cd $CurDir
```


### Test if all contigs have a matching gff and remove any which do not

```bash
Prefix=core_isolates
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality_RxLR

cd $WorkDir
for File in $(ls $PWD/contigs/*/*.fasta); do
    Prefix=$(basename "$File" .fasta)
    expected_gff="$PWD/gff/${Prefix}.gff"
    if [ ! -f "$expected_gff" ]; then
       rm -rf $(dirname $File)
    fi
done
cd $CurDir
```

```bash
Prefix=core_isolates
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality_RxLR
mkdir $WorkDir/fasta
cd $WorkDir/fasta
for File in $(ls ../contigs/*/*.fasta); do
    cp $File .
done
cd $CurDir
```

```bash
Prefix=core_isolates
cd analysis/popgen/SNP_calling/summary_stats/$Prefix/neutrality_RxLR
```

```R
library("PopGenome")
library("ggplot2")

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: number of segregating sites, Watterson's Theta, Tajima's D, Fu & Li's F*, Fu & Li'd D*
Pcac_Fa <- c(
"414_1", "15_13_1", "PC13_15_1", "404_1", "4040_1", "15_7_1", "2003_3_1", "416_1", "415_1", "4032_1", "P421_1", "12420_1",
"414_2", "15_13_2", "PC13_15_2", "404_2", "4040_2", "15_7_2", "2003_3_2", "416_2", "415_2", "4032_2", "P421_2", "12420_2"
)
Pcac_Md <- c(
  "P295_1", "62471_1", "R36_14_1",
  "P295_2", "62471_2", "R36_14_2"
)
Pi_ri <- c(
"371_1", "SCRP376_1", "SCRP370_1",
"371_2", "SCRP376_2", "SCRP370_2"
)
#Need to set argument diploid=TRUE if using diploid genomes in the below command:
populations <- list(Pcac_Fa, Pcac_Md, Pi_ri)
#Number of populations assigned above.
population_no <- length(populations)
population_names <- c("Pcac_Fa", "Pcac_Md", "Pi_ri")

GENOME.class <- readData("fasta/", format="fasta", gffpath = "gff/")
GENOME.class <- set.populations(GENOME.class, populations)
GENOME.class@genelength
# GENOME.class <- neutrality.stats(GENOME.class)
# GENOME.class@Tajima.D
# GENOME.class <- diversity.stats(GENOME.class)
# GENOME.class@Pi
# GENOME.class@nuc.diversity.within
# GENOME.class@hap.diversity.within

# WHOLE <- concatenate.regions(GENOME.class)
# WHOLE@genelength
# WHOLE <- neutrality.stats(WHOLE)
# WHOLE@Tajima.D
# WHOLE <- diversity.stats(WHOLE, pi=TRUE)
# WHOLE@Pi
# WHOLE@nuc.diversity.within
# WHOLE@hap.diversity.within

GENOME.class.split <- splitting.data(GENOME.class, subsites = "gene")
GENOME.class.split <- diversity.stats(GENOME.class.split, pi=TRUE)
GENOME.class.split@Pi
max(GENOME.class.split@Pi[,1])
max(GENOME.class.split@Pi[,2])
max(GENOME.class.split@Pi[,3])

TajimaD <- GENOME.class.split@Tajima.D
```

```
WHOLE@Tajima.D
              pop 1    pop 2    pop 3
Concatenate 1.37229 1.403231 1.584488
> WHOLE@Pi
            pop 1 pop 2 pop 3
Concatenate     0     0     0
> WHOLE@nuc.diversity.within
               pop 1 pop 2  pop 3
Concatenate 546.1884  4494 2792.2
> WHOLE@hap.diversity.within
            pop 1 pop 2 pop 3
Concatenate     1     1     1
```


# CRNs

## Comparison between Pcac crown rot, apple and P. idaei populations


## create the directory structure

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/core_isolates/neutrality_CRN
mkdir -p $WorkDir/all
mkdir -p $WorkDir/gff
mkdir -p $WorkDir/contigs
```

###copy the "gff" folder containing gff files

```bash
CurDir=/data/scratch/armita/idris
cd $CurDir
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/core_isolates/neutrality_CRN
AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
Headers=$WorkDir/headers.txt
cat $AnnotTab | cut -f1,17 | grep 'CRN' | cut -f1 > $Headers
Gff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.gff3)
FeatGff=$WorkDir/features.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $Headers $Gff effector ID > $FeatGff

cd $WorkDir/gff
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen
$ProgDir/summary_stats/split_gff_contig.sh $FeatGff
cd $CurDir
```

### Convert vcf files into fasta alignments

The vcf file needs to be converted into fasta alignment files. One file is
produced per reference contig

For the ploidy set <1|2|3>:
* 1 - haploid input
* 2 - diploid input, output as two separate phased haplotypes for each ind.
* 3 - diploid input, output as one sequence with ambiguity codes for each ind.


Files were copied from the total genes analysis

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/core_isolates/neutrality_CRN
cp analysis/popgen/SNP_calling/summary_stats/core_isolates/neutrality_genes/contigs/contig_*.fasta $WorkDir/contigs
```


### This folder contains only contig FASTA files
### So create a new "contigs" directory to hold those files:

###Next step: in the folder "contigs" create subfolders, each to hold one contig FASTA file

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/core_isolates/neutrality_CRN
cd $WorkDir/contigs
for File in $(ls *.fasta); do
    Prefix=${File%.fasta}
    mkdir $Prefix
    mv $File $Prefix/.
done
cd $CurDir
```


### Test if all contigs have a matching gff and remove any which do not

```bash

CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/core_isolates/neutrality_CRN

cd $WorkDir
for File in $(ls $PWD/contigs/*/*.fasta); do
    Prefix=$(basename "$File" .fasta)
    expected_gff="$PWD/gff/${Prefix}.gff"
    if [ ! -f "$expected_gff" ]; then
       rm -rf $(dirname $File)
    fi
done
cd $CurDir
```

```bash
CurDir=/data/scratch/armita/idris
WorkDir=$CurDir/analysis/popgen/SNP_calling/summary_stats/core_isolates/neutrality_CRN
mkdir $WorkDir/fasta
cd $WorkDir/fasta
for File in $(ls ../contigs/*/*.fasta); do
    cp $File .
done
cd $CurDir
```

```bash
cd analysis/popgen/SNP_calling/summary_stats/core_isolates/neutrality_CRN
```

```R
library("PopGenome")
library("ggplot2")

######################## BEFORE RUNNING ################################
#Assign individuals to appropriate populations (or just 1!)
#This script calculates: number of segregating sites, Watterson's Theta, Tajima's D, Fu & Li's F*, Fu & Li'd D*
Pcac_Fa <- c(
"414_1", "15_13_1", "PC13_15_1", "404_1", "4040_1", "15_7_1", "2003_3_1", "416_1", "415_1", "4032_1", "P421_1", "12420_1",
"414_2", "15_13_2", "PC13_15_2", "404_2", "4040_2", "15_7_2", "2003_3_2", "416_2", "415_2", "4032_2", "P421_2", "12420_2"
)
Pcac_Md <- c(
  "P295_1", "62471_1", "R36_14_1",
  "P295_2", "62471_2", "R36_14_2"
)
Pi_ri <- c(
"371_1", "SCRP376_1", "SCRP370_1",
"371_2", "SCRP376_2", "SCRP370_2"
)
#Need to set argument diploid=TRUE if using diploid genomes in the below command:
populations <- list(Pcac_Fa, Pcac_Md, Pi_ri)
#Number of populations assigned above.
population_no <- length(populations)
population_names <- c("Pcac_Fa", "Pcac_Md", "Pi_ri")

Genome.Class <- readData("fasta/", format="fasta", gffpath = "gff/")
GENOME.class <- set.populations(GENOME.class, populations)
# GENOME.class <- neutrality.stats(GENOME.class)
# GENOME.class@Tajima.D
# GENOME.class <- diversity.stats(GENOME.class)
# GENOME.class@Pi
# GENOME.class@nuc.diversity.within
# GENOME.class@hap.diversity.within

WHOLE <- concatenate.regions(GENOME.class)
WHOLE <- neutrality.stats(WHOLE)
WHOLE@Tajima.D
WHOLE <- diversity.stats(WHOLE, pi=TRUE)
WHOLE@Pi
WHOLE@nuc.diversity.within
WHOLE@hap.diversity.within
```

```
WHOLE@Tajima.D
              pop 1    pop 2    pop 3
Concatenate 1.37229 1.403231 1.584488
> WHOLE@Pi
            pop 1 pop 2 pop 3
Concatenate     0     0     0
> WHOLE@nuc.diversity.within
               pop 1 pop 2  pop 3
Concatenate 546.1884  4494 2792.2
> WHOLE@hap.diversity.within
            pop 1 pop 2 pop 3
Concatenate     1     1     1
```
