# LD deay


## Within Crown rot isolates:

Make a slimmed vcf file of crown rot isolates (also done for 4GT)

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="11-40 17-21 P295 62471 R36_14 371 SCRP376 SCRP370 10300 LV007"
Prefix=crown-rot_isolates
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```

```bash
CurDir=$PWD
Vcf=$(ls analysis/popgen/SNP_calling/summary_stats/crown-rot_isolates/crown-rot_isolates_filtered.recode.vcf)
OutDir=$(dirname $Vcf)
cd $OutDir
ProgDir=/home/adamst/git_repos/scripts/popgen
qsub $ProgDir/snp/sub_beagle.sh $(basename $Vcf)
cd $CurDir
```
```
gunzip analysis/popgen/SNP_calling/summary_stats/crown-rot_isolates/crown-rot_isolates_filtered.recode_haplo.vcf.gz
```


```bash
  Vcf=$(ls analysis/popgen/SNP_calling/summary_stats/crown-rot_isolates/crown-rot_isolates_filtered.recode_haplo.vcf)
  Prefix=Pcac_Fa_vs_414
  OutDir=analysis/popgen/SNP_calling/LD/$Prefix
  mkdir -p $OutDir
  VcfTools=/home/sobczm/bin/vcftools/bin
  $VcfTools/vcftools --vcf $Vcf --hap-r2 --ld-window-bp 1000000 --out $OutDir/${Prefix}_ld_window_1mb
```


```bash
cat analysis/popgen/SNP_calling/LD/Pcac_Fa_vs_414/Pcac_Fa_vs_414_ld_window_100kb.hap.ld | cut -f2,3,5 > analysis/popgen/SNP_calling/LD/Pcac_Fa_vs_414/Pcac_Fa_vs_414_ld_window_100kb_slimmed.hap.ld

cat analysis/popgen/SNP_calling/LD/Pcac_Fa_vs_414/Pcac_Fa_vs_414_ld_window_1mb.hap.ld | cut -f2,3,5 > analysis/popgen/SNP_calling/LD/Pcac_Fa_vs_414/Pcac_Fa_vs_414_ld_window_100kb_slimmed.hap.ld
```

```R
library(ggplot2)
df1 <- read.delim("analysis/popgen/SNP_calling/LD/Pcac_Fa_vs_414/Pcac_Fa_vs_414_ld_window_1mb_slimmed.hap.ld", header=TRUE)
colnames(df1) <- c("PointA","PointB","R2")
df1$PointA <- as.numeric(df1$PointA)
df1$PointB <- as.numeric(df1$PointB)

df1$dist <- df1$PointB - df1$PointA

# Make Linkage ddecay model
distance<-df1$dist
LD.data<-df1$R2
n<-288
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))

# Calculating LD50 values (half maximum value)
df<-data.frame(distance,fpoints)
maxld<-max(LD.data)

#You could elucubrate if it's better to use the maximum ESTIMATED value of LD
#In that case you just set: maxld<-max(fpoints)
maxld<-max(fpoints)
h.decay<-maxld/2
half.decay.distance<-df$distance[which.min(abs(df$fpoints-h.decay))]

# Make plot
p <- ggplot(df1, aes(dist, R2))
p <- p + geom_line(stat = "summary_bin", binwidth = 1000)
#p <- p + geom_segment(aes(x = half.decay.distance, y = 0.1, xend = half.decay.distance, yend = 0), arrow = arrow(length = unit(0.5, "cm")), color = "red") */
ggsave('analysis/popgen/SNP_calling/LD/Pcac_Fa_vs_414/Pcac_Fa_vs_414_LD_1mb_plot.pdf', p)
b <- ggplot(df1, aes(x = 'Crown rot', y = R2)) +
 geom_boxplot()

b <- b + labs(x = '', y = "Nucleotide diversity (\u03C0)")
ggsave('Pcac_Pi_boxplot.jpg', p)
```

<!--
```R
library(ggplot2)
library("data.table")
df1 <- read.delim("../Pcac_Fa_vs_414_ld_window_100kb_slimmed.hap.ld", header=TRUE)
colnames(df1) <- c("PointA","PointB","R2")
df1$PointA <- as.numeric(df1$PointA)
df1$PointB <- as.numeric(df1$PointB)

df1$dist <- df1$PointB - df1$PointA
window_size <- 100000
bin_size <- 1000
Cstart <- 0.2
data <- c(tapply(df1$R2, cut(df1$dist, seq(0, window_size,
    by = bin_size)), mean))
data_df <- data.frame(data)
setDT(data_df, keep.rownames = TRUE)
colnames(data_df) <- c("Distance", "Rsqd")

brks <- seq(0, window_size, bin_size)
Max_val <- window_size / bin_size
ints <- seq(1, Max_val, by = 1)
data_df$midpoint <- (head(brks, -1) + diff(brks) / 2)[ints]
Rsqd <- data_df$Rsqd
midpoint <- data_df$midpoint
Cstart <- c(C = Cstart)
n<-288

# Fit binned data to Hills and Weir decay function (a non-linear model)
# Following code in script adapted from
# https://jujumaan.com/2017/07/15/linkage-disequilibrium-decay-plot/

modelC <- nls(Rsqd~ ( (10 + C * midpoint) / ( (2 + C * midpoint) * (11 + C *
    midpoint))) * (1 + ( (3 + C * midpoint) * (12 + 12 * C * midpoint + (C *
        midpoint) ^ 2)) / (n * (2 + C * midpoint) * (11 + C * midpoint))),
        start = Cstart, control = nls.control(maxiter = 100))

# extract rho, the recombination parameter, 4Nr
rho <- summary(modelC)$parameters[1]

# Use the new rho value to obtain LD values adjusted for their distances

newrsqd <- ( (10 + rho * data_df$midpoint) / ( (2 + rho * data_df$midpoint) *
(11 + rho * data_df$midpoint))) * (1 + ( (3 + rho * data_df$midpoint) * (12 +
    12 * rho * data_df$midpoint + (rho * data_df$midpoint) ^ 2)) / (n * (2 +
        rho * data_df$midpoint) * (11 + rho * data_df$midpoint)))


fitted_data <- data.frame(data_df$midpoint, newrsqd)
max_rsqd <- max(fitted_data$newrsqd)
half_decay <- max_rsqd * 0.5
half_decay_dist <- fitted_data$data_df.midpoint[which.min(abs(
    fitted_data$newrsqd - half_decay))]
fitted_data <- fitted_data[order(fitted_data$data_df.midpoint), ]

# Identify point where r^2 = 0.2

rsqd_pt2 <- fitted_data$data_df.midpoint[which.min(abs(fitted_data$newrsqd
    - 0.2))]

cat("Half decay distance of LD r^2:", half_decay_dist, units, "\n")
cat("Distance where r^2 = 0.2:", rsqd_pt2, units, "\n")
``` -->

```bash
# 6 chromsomes * 12 individuals = 72
# 10 chromsomes * 12 individuals = 120
# 24 contigs * 12 individuals = 288
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC154851/
for Size in 72 120 288; do
   LD_file=$(ls analysis/popgen/SNP_calling/LD/Pcac_Fa_vs_414/Pcac_Fa_vs_414_ld_window_100kb.hap.ld)
   # LD_file=$(ls analysis/popgen/SNP_calling/LD/Pcac_Fa_vs_414/Pcac_Fa_vs_414_ld_window_1mb.hap.ld)
   Out_file=$(dirname $LD_file)/r^2_decay_"$Size".pdf
   units=bp
   window_size=100000
   # window_size=1000000
   bin_size=1000
   Cstart=0.1
   ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
   echo "Sample size of $Size:"
   Rscript --vanilla $ProgDir/plot_LD_decay.R --out_file $Out_file --Chromosome_number $Size --LD_statistics $LD_file --units $units --window_size $window_size --bin_size $bin_size --Cstart $Cstart
   printf "\n"
done
```

```
Sample size of 72:
Half decay distance of LD r^2: 500 bp
Distance where r^2 = 0.2: 500 bp

Sample size of 120:
Half decay distance of LD r^2: 500 bp
Distance where r^2 = 0.2: 500 bp

Sample size of 288:
Half decay distance of LD r^2: 500 bp
Distance where r^2 = 0.2: 500 bp
```

## Within all P. cactorum isolates:

Make a slimmed vcf file of cactorum isolates

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="371 SCRP376 SCRP370"
Prefix=Pcactorum_isolates
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```

Phase vcf files

```bash
Prefix=Pcactorum_isolates
CurDir=$PWD
Vcf=$(ls analysis/popgen/SNP_calling/summary_stats/$Prefix/${Prefix}_filtered.recode.vcf)
OutDir=$(dirname $Vcf)
cd $OutDir
ProgDir=/home/adamst/git_repos/scripts/popgen
qsub $ProgDir/snp/sub_beagle.sh $(basename $Vcf)
cd $CurDir
```

unzip the output file
```bash
Prefix=Pcactorum_isolates
gunzip analysis/popgen/SNP_calling/summary_stats/$Prefix/${Prefix}_filtered.recode_haplo.vcf.gz
```

Calculate LD values
```bash
  Prefix=Pcactorum_isolates
  Vcf=$(ls analysis/popgen/SNP_calling/summary_stats/$Prefix/${Prefix}_filtered.recode_haplo.vcf)
  OutDir=analysis/popgen/SNP_calling/LD/$Prefix
  mkdir -p $OutDir
  VcfTools=/home/sobczm/bin/vcftools/bin
  $VcfTools/vcftools --vcf $Vcf --hap-r2 --ld-window-bp 100000 --out $OutDir/${Prefix}_ld_window_100kb
  $VcfTools/vcftools --vcf $Vcf --hap-r2 --ld-window-bp 1000000 --out $OutDir/${Prefix}_ld_window_1mb
```

```
After filtering, kept 35416 out of a possible 35428 Sites
After filtering, kept 35416 out of a possible 35428 Sites
```
```bash
Prefix=Pcactorum_isolates
# 6 chromsomes * 17 individuals = 102
# 10 chromsomes * 17 individuals = 170
# 24 contigs * 17 individuals = 408
for Size in 72 120 288; do
   LD_file=$(ls analysis/popgen/SNP_calling/LD/$Prefix/${Prefix}_ld_window_100kb.hap.ld)
   Out_file=$(dirname $LD_file)/r^2_decay_"$Size"_100kb.pdf
   units=bp
   window_size=100000
   bin_size=1000
   Cstart=0.1
   ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
   echo "Sample size of $Size:"
   Rscript --vanilla $ProgDir/plot_LD_decay.R --out_file $Out_file --Chromosome_number $Size --LD_statistics $LD_file --units $units --window_size $window_size --bin_size $bin_size --Cstart $Cstart
   printf "\n"

   Out_file=$(dirname $LD_file)/r^2_decay_"$Size"_1mb.pdf
   units=bp
   window_size=1000000
   bin_size=1000
   Cstart=0.1
   ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae/popgen_analysis
   echo "Sample size of $Size:"
   Rscript --vanilla $ProgDir/plot_LD_decay.R --out_file $Out_file --Chromosome_number $Size --LD_statistics $LD_file --units $units --window_size $window_size --bin_size $bin_size --Cstart $Cstart
   printf "\n"
done
```

```
Sample size of 72:
Half decay distance of LD r^2: 99500 bp
Distance where r^2 = 0.2: 99500 bp

Sample size of 72:
Half decay distance of LD r^2: 580500 bp
Distance where r^2 = 0.2: 777500 bp

Sample size of 120:
Half decay distance of LD r^2: 99500 bp
Distance where r^2 = 0.2: 99500 bp

Sample size of 120:
Half decay distance of LD r^2: 631500 bp
Distance where r^2 = 0.2: 828500 bp

Sample size of 288:
Half decay distance of LD r^2: 99500 bp
Distance where r^2 = 0.2: 99500 bp

Sample size of 288:
Half decay distance of LD r^2: 683500 bp
Distance where r^2 = 0.2: 882500 bp
```

Plotting the raw data:


```bash
Prefix=Pcactorum_isolates
cat analysis/popgen/SNP_calling/LD/$Prefix/${Prefix}_ld_window_100kb.hap.ld | cut -f2,3,5 > analysis/popgen/SNP_calling/LD/$Prefix/${Prefix}_ld_window_100kb_slimmed.hap.ld
cat analysis/popgen/SNP_calling/LD/$Prefix/${Prefix}_ld_window_1mb.hap.ld | cut -f2,3,5 > analysis/popgen/SNP_calling/LD/$Prefix/${Prefix}_ld_window_1mb_slimmed.hap.ld
```

```R
library(ggplot2)
df1 <- read.delim("analysis/popgen/SNP_calling/LD/Pcactorum_isolates/Pcactorum_isolates_ld_window_1mb_slimmed.hap.ld", header=TRUE)
colnames(df1) <- c("PointA","PointB","R2")
df1$PointA <- as.numeric(df1$PointA)
df1$PointB <- as.numeric(df1$PointB)

df1$dist <- df1$PointB - df1$PointA

# Make plot
p <- ggplot(df1, aes(dist, R2))
p <- p + geom_line(stat = "summary_bin", binwidth = 1000)
ggsave('analysis/popgen/SNP_calling/LD/Pcactorum_isolates/Pcactorum_isolates_LD_1mb_plot.pdf', p, width = 14, height = 7)

b <- ggplot(df1, aes(x = 'P.cactorum', y = R2)) +
 geom_boxplot()
b <- b + labs(x = '', y = "R^2")
ggsave('analysis/popgen/SNP_calling/LD/Pcactorum_isolates/Pcactorum_R2_boxplot.pdf', b)
```
