# Population statistics calulated using VCFtools

This includes calculating the number of segregating sites, Watterson's Theta, Tajima's D, Fu & Li's F* and Fu & Li'd D

```bash
cd /data/scratch/armita/idris
```

## Comparison within all Pcac populations

Slim down the vcf to the three populations (11-40 included in strawberry lineage):

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="371 SCRP370 SCRP376"
Prefix=within_cactorum
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 1.0 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```

```
After filtering, kept 17 out of 17 Individuals
Outputting VCF file...
After filtering, kept 35428 out of a possible 329672 Sites
Run Time = 9.00 seconds
```

Set populations

```bash
Prefix=within_cactorum
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir $OutDir/vcftools
Pc_Fxa="415 15_7 15_13 404 12420 2003_3 4032 PC13_15 414 4040 416 P421 11-40"
for Strain in $(echo $Pc_Fxa); do
echo $Strain
done > $OutDir/vcftools/population_Pc_Fxa.txt
Pc_Mxd="62471 P295 R36_14"
for Strain in $(echo $Pc_Mxd); do
echo $Strain
done > $OutDir/vcftools/population_Pc_Mxd.txt
Pc_LR="17-21"
for Strain in $(echo $Pc_LR); do
echo $Strain
done > $OutDir/vcftools/population_Pc_LR.txt
```

```bash
Prefix=within_cactorum
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
Vcf=$(ls analysis/popgen/SNP_calling/summary_stats/within_cactorum/within_cactorum_filtered.recode.vcf)

# Fst
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $Vcf --weir-fst-pop $OutDir/vcftools/population_Pc_Fxa.txt --weir-fst-pop $OutDir/vcftools/population_Pc_Mxd.txt --weir-fst-pop $OutDir/vcftools/population_Pc_LR.txt --out $OutDir/vcftools/$Prefix

$VcfTools/vcftools --vcf $Vcf --weir-fst-pop $OutDir/vcftools/population_Pc_Fxa.txt --out $OutDir/vcftools/${Prefix}_Pc_Fxa --weir-fst-pop $OutDir/vcftools/population_Pc_Mxd.txt --out $OutDir/vcftools/${Prefix}_Pc_Fxa_vs_Mxd

# OUTPUT NUCLEOTIDE DIVERGENCE STATISTICS Pi
# Calculate a measure of heterozygosity on a per-individual basis
# Tajima’s D statistic in bins with size of the specified number
# p-value for each site from a Hardy-Weinberg Equilibrium test
# calculate and output a relatedness statistic based on Yang et al, Nature Genetics 2010
# Calculate the number and density of SNPs

VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $Vcf --site-pi $OutDir/vcftools/$Prefix
$VcfTools/vcftools --vcf $Vcf --het --out $OutDir/vcftools/$Prefix
$VcfTools/vcftools --vcf $Vcf --hardy --out $OutDir/vcftools/$Prefix
$VcfTools/vcftools --vcf $Vcf --TajimaD 100000 --out $OutDir/vcftools/$Prefix
$VcfTools/vcftools --vcf $Vcf --relatedness --out $OutDir/vcftools/$Prefix
$VcfTools/vcftools --vcf $Vcf --SNPdensity 100000 --out $OutDir/vcftools/$Prefix
```

heterozygocity stats - this may be more appropriately applied using SNPs just called within a single population.
```
INDV    O(HOM)  E(HOM)  N_SITES F
11-40   35351   26628.9 35416   0.99260
12420   35346   26628.9 35416   0.99203
15_13   35351   26628.9 35416   0.99260
15_7    35335   26628.9 35416   0.99078
17-21   35260   26628.9 35416   0.98225
2003_3  35350   26628.9 35416   0.99249
4032    35353   26628.9 35416   0.99283
404     35348   26628.9 35416   0.99226
4040    35340   26628.9 35416   0.99135
414     35350   26628.9 35416   0.99249
415     35351   26628.9 35416   0.99260
416     35340   26628.9 35416   0.99135
62471   35290   26628.9 35416   0.98566
P295    35292   26628.9 35416   0.98589
P421    35351   26628.9 35416   0.99260
PC13_15 35349   26628.9 35416   0.99238
R36_14  35270   26628.9 35416   0.98338
```

```
Keeping individuals in 'keep' list
After filtering, kept 17 out of 17 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.86067
Weir and Cockerham weighted Fst estimate: 0.93803
After filtering, kept 35428 out of a possible 35428 Sites
Run Time = 1.00 seconds


Keeping individuals in 'keep' list
After filtering, kept 16 out of 17 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.84511
Weir and Cockerham weighted Fst estimate: 0.93879
After filtering, kept 35428 out of a possible 35428 Sites
Run Time = 1.00 seconds
```


## Comparison within Pcac crown rot population

Slim down the vcf to the two populations (11-40 included in strawberry lineage):

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="17-21 10300 LV007 371 SCRP370 SCRP376 62471 P295 R36_14"
Prefix=within_strawberry
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 1.0 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```

```
After filtering, kept 13 out of 13 Individuals
Outputting VCF file...
After filtering, kept 1535 out of a possible 329672 Sites
Run Time = 6.00 seconds
```


Set populations

```bash
Prefix=within_strawberry
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir $OutDir/vcftools
Pc_Fxa_i="415 15_7 15_13 404 12420 11-40"
for Strain in $(echo $Pc_Fxa_i); do
echo $Strain
done > $OutDir/vcftools/population_Pc_Fxa_i.txt
Pc_Fxa_ii="4040 416"
for Strain in $(echo $Pc_Fxa_ii); do
echo $Strain
done > $OutDir/vcftools/population_Pc_Fxa_ii.txt
Pc_Fxa_iii="2003_3 4032 PC13_15 414 P421"
for Strain in $(echo $Pc_Fxa_iii); do
echo $Strain
done > $OutDir/vcftools/population_Pc_Fxa_iii.txt
```


```bash
Prefix=within_strawberry
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
Vcf=$(ls $OutDir/"$Prefix"_filtered.recode.vcf)

# Fst
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $Vcf --weir-fst-pop $OutDir/vcftools/population_Pc_Fxa_i.txt --weir-fst-pop $OutDir/vcftools/population_Pc_Fxa_ii.txt --weir-fst-pop $OutDir/vcftools/population_Pc_Fxa_iii.txt --out $OutDir/vcftools/$Prefix

# OUTPUT NUCLEOTIDE DIVERGENCE STATISTICS Pi
# Calculate a measure of heterozygosity on a per-individual basis
# Tajima’s D statistic in bins with size of the specified number
# p-value for each site from a Hardy-Weinberg Equilibrium test
# calculate and output a relatedness statistic based on Yang et al, Nature Genetics 2010
# Calculate the number and density of SNPs

VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $Vcf --site-pi $OutDir/vcftools/$Prefix
$VcfTools/vcftools --vcf $Vcf --het --out $OutDir/vcftools/$Prefix
$VcfTools/vcftools --vcf $Vcf --hardy --out $OutDir/vcftools/$Prefix
$VcfTools/vcftools --vcf $Vcf --TajimaD 100000 --out $OutDir/vcftools/$Prefix
$VcfTools/vcftools --vcf $Vcf --relatedness --out $OutDir/vcftools/$Prefix
$VcfTools/vcftools --vcf $Vcf --SNPdensity 100000 --out $OutDir/vcftools/$Prefix
```

```
After filtering, kept 13 out of 13 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.82335
Weir and Cockerham weighted Fst estimate: 0.92021
After filtering, kept 1535 out of a possible 1535 Sites
Run Time = 0.00 seconds
```

```
INDV    O(HOM)  E(HOM)  N_SITES F
11-40   1470    1004.8  1535    0.87740
12420   1465    1004.8  1535    0.86797
15_13   1470    1004.8  1535    0.87740
15_7    1454    1004.8  1535    0.84722
2003_3  1469    1004.8  1535    0.87551
4032    1472    1004.8  1535    0.88117
404     1467    1004.8  1535    0.87174
4040    1459    1004.8  1535    0.85665
414     1469    1004.8  1535    0.87551
415     1470    1004.8  1535    0.87740
416     1459    1004.8  1535    0.85665
P421    1470    1004.8  1535    0.87740
PC13_15 1468    1004.8  1535    0.87362
```
