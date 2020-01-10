
# Identify the population structure within populations

This is done using faststructure.

Structure was investigated within all P. cactorum isolates and within
the crown rot P. cactorum population.

## Within all P. cactorum isolates:

Make a slimmed vcf file of P. cactorum isolates

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="371 SCRP376 SCRP370"
Prefix=Pcactorum_isolates
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/${Prefix}_filtered
```

## Convert from VCF to Plink's PED format

```bash
Prefix=Pcactorum_isolates
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix/fastStructure
mkdir -p $OutDir
cd $OutDir
cp ../${Prefix}_filtered.recode.vcf .
input_file=${Prefix}_filtered.recode.vcf
ProgDir=/home/vicker/programs/plink-1.90beta
$ProgDir/plink --allow-extra-chr --const-fid 0 --vcf $input_file --recode --make-bed --out ${input_file%.vcf} > ${input_file%.vcf}.log

# Run FastStructure with different k values
s=1
f=6
for k in $(seq $s $f); do
  ProgDir=/home/adamst/git_repos/scripts/popgen/snp
  qsub $ProgDir/sub_fast_structure.sh ${input_file%.vcf} $k
done
```
Choosing model complexity (K) among all the K values tested

```bash
structure=/home/sobczm/bin/fastStructure
python $structure/chooseK.py --input=${input_file%.vcf} > ${input_file%.vcf}_K_choice
```

```bash
#Generate sample lables
cut -f2 ${input_file%.vcf}.fam | cut -d " " -f2 > ${input_file%.vcf}.lab

#Draw output
# If errors here then close Xquartz on your local machine, open a new
# terminal session and try again.
for i in $(seq $s $f); do
  structure=/home/sobczm/bin/fastStructure
  python $structure/distruct_mod.py -K $i --input=${input_file%.vcf} --output=${input_file%.vcf}_${i}.svg --title K$i --popfile=${input_file%.vcf}.lab
done

rm sub_fast_structure.sh.*

cd $CurDir
```


## Within the P. cactorum cown rot isolates:

Make a slimmed vcf file of P. cactorum isolates




## Within Crown rot isolates:

Make a slimmed vcf file of crown rot isolates (also done for LD decay)

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

## Convert from VCF to Plink's PED format

```bash
Prefix=crown-rot_isolates
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix/fastStructure
mkdir -p $OutDir
cd $OutDir
cp ../${Prefix}_filtered.recode.vcf .
input_file=${Prefix}_filtered.recode.vcf
ProgDir=/home/vicker/programs/plink-1.90beta
$ProgDir/plink --allow-extra-chr --const-fid 0 --vcf $input_file --recode --make-bed --out ${input_file%.vcf} > ${input_file%.vcf}.log

# Run FastStructure with different k values
s=1
f=6
for k in $(seq $s $f); do
  ProgDir=/home/adamst/git_repos/scripts/popgen/snp
  qsub $ProgDir/sub_fast_structure.sh ${input_file%.vcf} $k
done
```
Choosing model complexity (K) among all the K values tested

```bash
structure=/home/sobczm/bin/fastStructure
python $structure/chooseK.py --input=${input_file%.vcf} > ${input_file%.vcf}_K_choice
```

```bash
#Generate sample lables
cut -f2 ${input_file%.vcf}.fam | cut -d " " -f2 > ${input_file%.vcf}.lab

#Draw output
# If errors here then close Xquartz on your local machine, open a new
# terminal session and try again.
for i in $(seq $s $f); do
  structure=/home/sobczm/bin/fastStructure
  python $structure/distruct_mod.py -K $i --input=${input_file%.vcf} --output=${input_file%.vcf}_${i}.svg --title K$i --popfile=${input_file%.vcf}.lab
done

rm sub_fast_structure.sh.*

cd $CurDir
```


## Within Strawberry clade isolates:

Make a slimmed vcf file of crown rot isolates (also done for LD decay)

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors.vcf)
ExcludeList="17-21 P295 62471 R36_14 371 SCRP376 SCRP370 10300 LV007"
Prefix=strawberry-lineage_isolates
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix
mkdir -p $OutDir
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered
```

## Convert from VCF to Plink's PED format

```bash
Prefix=strawberry-lineage_isolates
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling/summary_stats/$Prefix/fastStructure
mkdir -p $OutDir
cd $OutDir
cp ../${Prefix}_filtered.recode.vcf .
input_file=${Prefix}_filtered.recode.vcf
ProgDir=/home/vicker/programs/plink-1.90beta
$ProgDir/plink --allow-extra-chr --const-fid 0 --vcf $input_file --recode --make-bed --out ${input_file%.vcf} > ${input_file%.vcf}.log

# Run FastStructure with different k values
s=1
f=6
for k in $(seq $s $f); do
  ProgDir=/home/adamst/git_repos/scripts/popgen/snp
  qsub $ProgDir/sub_fast_structure.sh ${input_file%.vcf} $k
done
```
Choosing model complexity (K) among all the K values tested

```bash
structure=/home/sobczm/bin/fastStructure
python $structure/chooseK.py --input=${input_file%.vcf} > ${input_file%.vcf}_K_choice
```

```bash
#Generate sample lables
cut -f2 ${input_file%.vcf}.fam | cut -d " " -f2 > ${input_file%.vcf}.lab

#Draw output
# If errors here then close Xquartz on your local machine, open a new
# terminal session and try again.
for i in $(seq $s $f); do
  structure=/home/sobczm/bin/fastStructure
  python $structure/distruct_mod.py -K $i --input=${input_file%.vcf} --output=${input_file%.vcf}_${i}.svg --title K$i --popfile=${input_file%.vcf}.lab
done

rm sub_fast_structure.sh.*

cd $CurDir
```


The meanQ data for K=3 was downloaded to my local computer and plotted in R:

```R
library('ggplot2')
df1 <- read.table("~/Downloads/Pc/fastStructure_strawb-lineage/strawberry-lineage_isolates_filtered.recode.3.meanQ", quote="\"", comment.char="")
df_names <- read.table("~/Downloads/Pc/fastStructure_strawb-lineage/strawberry-lineage_isolates_filtered.recode.lab", quote="\"", comment.char="")
df_names$corrected_names <- df_names$V1
df_names$corrected_names <- gsub("414", "P414", df_names$corrected_names)
df_names$corrected_names <- gsub("15_7", "15-7", df_names$corrected_names)
df_names$corrected_names <- gsub("15_13", "15-13", df_names$corrected_names)
df_names$corrected_names <- gsub("PC13_15", "PC13/15", df_names$corrected_names)
df_names$corrected_names <- gsub("12420", "12-420", df_names$corrected_names)
df_names$corrected_names <- gsub("416", "P416", df_names$corrected_names)
df_names$corrected_names <- gsub("\\<404\\>", "P404", df_names$corrected_names)
df_names$corrected_names <- gsub("415", "P415", df_names$corrected_names)
df_names$corrected_names <- gsub("2003_3", "2003-3", df_names$corrected_names)
df_names$corrected_names <- factor(df_names$corrected_names,levels = c("P415", "15-7", "P404", "15-13", "12-420", "11-40", "P416", "4040", "PC13/15", "P414", "P421", "4032", "2003-3"))

df1$'Isolate' <- df_names$corrected_names
colnames(df1) <- c("i","ii","iii", "Isolate")
library(reshape2)
df2 <- melt(df1, id.vars=c("Isolate"))
colnames(df2) <- c("Isolate", "Population", "Proportion")


p <- ggplot() + geom_col(data = df2, aes(x = Isolate, y = Proportion, fill = Population))

p <- p + scale_y_continuous(limits=c(0, 1), expand = c(0.005, 0.005)) + geom_hline(yintercept=c(0,1))
p <- p + theme(
  axis.text.x = element_text(size=12, angle=-45),
  axis.text.y = element_text(size=12),
  axis.title.x = element_text(size=14),
  axis.title.y = element_text(size=14),
  legend.title = element_text(size=14),
  legend.text = element_text(size=12)
        )
ggsave("Pcactorum_strawberry_lineage_structure.pdf", p, width = 10, height = 5)

# Pcac data
df3 <- read.table("~/Downloads/Pc/fastStructure/Pcactorum_isolates_filtered.recode.3.meanQ", quote="\"", comment.char="")

df2_names <- read.table("~/Downloads/Pc/fastStructure/Pcactorum_isolates_filtered.recode.lab", quote="\"", comment.char="")
df2_names$corrected_names <- df2_names$V1
df2_names$corrected_names <- gsub("414", "P414", df2_names$corrected_names)
df2_names$corrected_names <- gsub("15_7", "15-7", df2_names$corrected_names)
df2_names$corrected_names <- gsub("15_13", "15-13", df2_names$corrected_names)
df2_names$corrected_names <- gsub("PC13_15", "PC13/15", df2_names$corrected_names)
df2_names$corrected_names <- gsub("12420", "12-420", df2_names$corrected_names)
df2_names$corrected_names <- gsub("416", "P416", df2_names$corrected_names)
df2_names$corrected_names <- gsub("\\<404\\>", "P404", df2_names$corrected_names)
df2_names$corrected_names <- gsub("415", "P415", df2_names$corrected_names)
df2_names$corrected_names <- gsub("2003_3", "2003-3", df2_names$corrected_names)
df2_names$corrected_names <- gsub("R36_14", "R36/14", df2_names$corrected_names)
df2_names$corrected_names <- factor(df2_names$corrected_names,levels = c("P415", "15-7", "P404", "15-13", "12-420", "11-40", "P416", "4040", "PC13/15", "P414", "P421", "4032", "2003-3", "62471", "P295", "R36/14", "17-21"))

df3$'Isolate' <- df2_names$corrected_names
colnames(df3) <- c("a","b","c", "Isolate")
library(reshape2)
df4 <- melt(df3, id.vars=c("Isolate"))
colnames(df4) <- c("Isolate", "Population", "Proportion")


p2 <- ggplot() + geom_col(data = df4, aes(x = Isolate, y = Proportion, fill = Population))

p2 <- p2 + scale_y_continuous(limits=c(0, 1), expand = c(0.005, 0.005)) + geom_hline(yintercept=c(0,1))
p2 <- p2 + theme(
  axis.text.x = element_text(size=12, angle=-45),
  axis.text.y = element_text(size=12),
  axis.title.x = element_text(size=14),
  axis.title.y = element_text(size=14),
  legend.title = element_text(size=14),
  legend.text = element_text(size=12)
        )

ggsave("Pcactorum_structure.pdf", p2, width = 10, height = 5)

layout(matrix(1:1, 1, 1))
p
p2
layout(matrix(1))
```
