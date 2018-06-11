
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
