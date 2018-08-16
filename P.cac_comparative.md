Phytophthora
============

Scripts used in the analysis of Phytophthora genomes

ands used during analysis of phytophthora genomes. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/idris

The following is a summary of the work presented in this Readme:
Data organisation:
  * Preparing data  
Draft Genome assembly
  * Data qc
  * Genome assembly
  * Repeatmasking
  * Gene prediction
  * Functional annotation
Genome analysis
  * Homology between predicted genes & published effectors


#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.


```bash
  cd /home/groups/harrisonlab/project_files/idris
  mkdir -p raw_dna/paired/P.cactorum/404/F
  mkdir -p raw_dna/paired/P.cactorum/404/R
  mkdir -p raw_dna/paired/P.cactorum/414/F
  mkdir -p raw_dna/paired/P.cactorum/414/R
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/cactorum/Cactorum/130614_P404
  cp $RawDat/cactp404_S3_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/404/F/.
  cp $RawDat/cactp404_S3_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/404/R/.
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/cactorum/Cactorum/130624_P404
  cp $RawDat/130624_cactp404_S3_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/404/F/.
  cp $RawDat/130624_cactp404_S3_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/404/R/.
  RawDat="/home/groups/harrisonlab/raw_data/raw_seq/cactorum/Cactorum/Cactorum\ 130120"
  cp /home/groups/harrisonlab/raw_data/raw_seq/cactorum/Cactorum/Cactorum\ 130120/cact414_S2_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/414/F/414_run1_F.fastq.gz
  cp /home/groups/harrisonlab/raw_data/raw_seq/cactorum/Cactorum/Cactorum\ 130120/cact414_S2_L001_R2_001.fastq.gz  raw_dna/paired/P.cactorum/414/R/414_run1_R.fastq.gz
  RawDat="/home/groups/harrisonlab/raw_data/raw_seq/cactorum/Cactorum/Cactorum\ 130517"
  cp /home/groups/harrisonlab/raw_data/raw_seq/cactorum/Cactorum/Cactorum\ 130517/cact414_S2_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/414/F/414_run2_F.fastq.gz
  cp /home/groups/harrisonlab/raw_data/raw_seq/cactorum/Cactorum/Cactorum\ 130517/cact414_S2_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/414/R/414_run2_R.fastq.gz
  mkdir -p raw_dna/paired/P.cactorum/415/F
  mkdir -p raw_dna/paired/P.cactorum/415/R
  mkdir -p raw_dna/paired/P.cactorum/416/F
  mkdir -p raw_dna/paired/P.cactorum/416/R
  mkdir -p raw_dna/paired/P.cactorum/62471/F/.
  mkdir -p raw_dna/paired/P.cactorum/62471/R/.
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/150716_M01678_0023_AB0YF
  cp $RawDat/Pcactorum415_S2_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/415/F/.
  cp $RawDat/Pcactorum415_S2_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/415/R/.
  cp $RawDat/Pcactorum416_S1_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/416/F/.
  cp $RawDat/Pcactorum416_S1_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/416/R/.
  RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160108_M01678_0039_AEMMF
  cp $RawDatDir/62471_S3_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/62471/F/.
  cp $RawDatDir/62471_S3_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/62471/R/.
  # Data run on 13/10/16 161010_M04465_0026_000000000-APP5B
  RawDatDir=/home/miseq_readonly/miseq_data/2016/RAW/161010_M04465_0026_000000000-APP5B/Data/Intensities/BaseCalls
  Species=P.cactorum
  Strain=15_7
  mkdir -p raw_dna/paired/$Species/$Strain/F
  mkdir -p raw_dna/paired/$Species/$Strain/R
  cp $RawDatDir/Pcact15-07_S2_L001_R1_001.fastq.gz raw_dna/paired/$Species/$Strain/F/.
  cp $RawDatDir/Pcact15-07_S2_L001_R2_001.fastq.gz raw_dna/paired/$Species/$Strain/R/.
  Species=P.cactorum
  Strain=15_13
  mkdir -p raw_dna/paired/$Species/$Strain/F
  mkdir -p raw_dna/paired/$Species/$Strain/R
  cp $RawDatDir/Pcact15-13_S3_L001_R1_001.fastq.gz raw_dna/paired/$Species/$Strain/F/.
  cp $RawDatDir/Pcact15-13_S3_L001_R2_001.fastq.gz raw_dna/paired/$Species/$Strain/R/.
  # Data run on 10/02/17 170210_M04465_0034_000000000-ATJ0E
  RawDatDir=/home/miseq_data/2017/RAW/170210_M04465_0034_000000000-ATJ0E/Data/Intensities/BaseCalls
  Date=170210
  Species=P.cactorum
  Strain=404
  mkdir -p raw_dna/paired/$Species/$Strain/F
  mkdir -p raw_dna/paired/$Species/$Strain/R
  cp $RawDatDir/P404_S3_L001_R1_001.fastq.gz raw_dna/paired/$Species/$Strain/F/"$Strain"_"$Date"_F.fastq.gz
  cp $RawDatDir/P404_S3_L001_R2_001.fastq.gz raw_dna/paired/$Species/$Strain/R/"$Strain"_"$Date"_R.fastq.gz
  Species=P.cactorum
  Strain=414
  mkdir -p raw_dna/paired/$Species/$Strain/F
  mkdir -p raw_dna/paired/$Species/$Strain/R
  cp $RawDatDir/P414_S4_L001_R1_001.fastq.gz raw_dna/paired/$Species/$Strain/F/"$Strain"_"$Date"_F.fastq.gz
  cp $RawDatDir/P414_S4_L001_R2_001.fastq.gz raw_dna/paired/$Species/$Strain/R/"$Strain"_"$Date"_R.fastq.gz
  Species=P.cactorum
  Strain=P295
  mkdir -p raw_dna/paired/$Species/$Strain/F
  mkdir -p raw_dna/paired/$Species/$Strain/R
  cp $RawDatDir/P295_S2_L001_R1_001.fastq.gz raw_dna/paired/$Species/$Strain/F/"$Strain"_"$Date"_F.fastq.gz
  cp $RawDatDir/P295_S2_L001_R2_001.fastq.gz raw_dna/paired/$Species/$Strain/R/"$Strain"_"$Date"_R.fastq.gz
  Species=P.cactorum
  Strain=12420
  mkdir -p raw_dna/paired/$Species/$Strain/F
  mkdir -p raw_dna/paired/$Species/$Strain/R
  cp $RawDatDir/12-420_S1_L001_R1_001.fastq.gz raw_dna/paired/$Species/$Strain/F/"$Strain"_"$Date"_F.fastq.gz
  cp $RawDatDir/12-420_S1_L001_R2_001.fastq.gz raw_dna/paired/$Species/$Strain/R/"$Strain"_"$Date"_R.fastq.gz

  # For missing data run for isolates 4032, 4040 and PC13/15
  #  RawDatDir=/home/miseq_data/2016/RAW/160311_M04465_0006_000000000-AKT5J/Data/Intensities/BaseCalls
  # Date=160311



  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160404_M004465_0008-ALVUT
  mkdir -p raw_dna/paired/P.cactorum/2003_3/F
  mkdir -p raw_dna/paired/P.cactorum/2003_3/R
  mkdir -p raw_dna/paired/P.idaei/SCRP370/F
  mkdir -p raw_dna/paired/P.idaei/SCRP370/R
  cp $RawDat/Pc20033_S2_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/2003_3/F/.
  cp $RawDat/Pc20033_S2_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/2003_3/R/.
  cp $RawDat/SCRP370_S3_L001_R1_001.fastq.gz  raw_dna/paired/P.idaei/SCRP370/F/.
  cp $RawDat/SCRP370_S3_L001_R2_001.fastq.gz  raw_dna/paired/P.idaei/SCRP370/R/.

  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160412_M04465_0010-AMLCU
  mkdir -p raw_dna/paired/P.rubi/SCRP249/F
  mkdir -p raw_dna/paired/P.rubi/SCRP249/R
  mkdir -p raw_dna/paired/P.rubi/SCRP324/F
  mkdir -p raw_dna/paired/P.rubi/SCRP324/R
  mkdir -p raw_dna/paired/P.rubi/SCRP333/F
  mkdir -p raw_dna/paired/P.rubi/SCRP333/R
  cp $RawDat/SCRP249_S1_L001_R1_001.fastq.gz raw_dna/paired/P.rubi/SCRP249/F/.
  cp $RawDat/SCRP249_S1_L001_R2_001.fastq.gz raw_dna/paired/P.rubi/SCRP249/R/.
  cp $RawDat/SCRP324_S2_L001_R1_001.fastq.gz raw_dna/paired/P.rubi/SCRP324/F/.
  cp $RawDat/SCRP324_S2_L001_R2_001.fastq.gz raw_dna/paired/P.rubi/SCRP324/R/.
  cp $RawDat/SCRP333_S3_L001_R1_001.fastq.gz raw_dna/paired/P.rubi/SCRP333/F/.
  cp $RawDat/SCRP333_S3_L001_R2_001.fastq.gz raw_dna/paired/P.rubi/SCRP333/R/.

  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160415_M004465_00011-AMLCL
  mkdir -p raw_dna/paired/P.cactorum/2003_3/F
  mkdir -p raw_dna/paired/P.cactorum/2003_3/R
  mkdir -p raw_dna/paired/P.cactorum/PC13_15/F
  mkdir -p raw_dna/paired/P.cactorum/PC13_15/R
  mkdir -p raw_dna/paired/P.idaei/SCRP376/F
  mkdir -p raw_dna/paired/P.idaei/SCRP376/R
  cp $RawDat/20033_S4_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/2003_3/F/.
  cp $RawDat/20033_S4_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/2003_3/R/.
  cp $RawDat/PC1315_S5_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/PC13_15/F/.
  cp $RawDat/PC1315_S5_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/PC13_15/R/.
  cp $RawDat/SCRP376_S1_L001_R1_001.fastq.gz raw_dna/paired/P.idaei/SCRP376/F/.
  cp $RawDat/SCRP376_S1_L001_R2_001.fastq.gz raw_dna/paired/P.idaei/SCRP376/R/.

  RawDat=/home/miseq_data/.tmp_nas_data/miseq_data/miseq_data/RAW/2016/160223_M04465_0003_000000000-AJ4TR/Data/Intensities/BaseCalls
  mkdir -p raw_dna/paired/P.cactorum/415/F
  mkdir -p raw_dna/paired/P.cactorum/415/R
  mkdir -p raw_dna/paired/P.cactorum/416/F
  mkdir -p raw_dna/paired/P.cactorum/416/R
  mkdir -p raw_dna/paired/P.cactorum/R36_14/F
  mkdir -p raw_dna/paired/P.cactorum/R36_14/R
  cp $RawDat/P415_S1_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/415/F/.
  cp $RawDat/P415_S1_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/415/R/.
  cp $RawDat/P416_S2_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/416/F/.
  cp $RawDat/P416_S2_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/416/R/.
  cp $RawDat/R3614_S3_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/R36_14/F/.
  cp $RawDat/R3614_S3_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/R36_14/R/.


  RawDat=/home/miseq_data/.tmp_nas_data/miseq_data/miseq_data/RAW/2016/160311_M04465_0006_000000000-AKT5J/Data/Intensities/BaseCalls
  mkdir -p raw_dna/paired/P.cactorum/4032/F
  mkdir -p raw_dna/paired/P.cactorum/4032/R
  mkdir -p raw_dna/paired/P.cactorum/4040/F
  mkdir -p raw_dna/paired/P.cactorum/4040/R
  mkdir -p raw_dna/paired/P.cactorum/PC13_15/F
  mkdir -p raw_dna/paired/P.cactorum/PC13_15/R
  cp $RawDat/Pcact-4032_S1_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/4032/F/.
  cp $RawDat/Pcact-4032_S1_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/4032/R/.
  cp $RawDat/Pcact-4040_S2_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/4040/F/.
  cp $RawDat/Pcact-4040_S2_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/4040/R/.
  cp $RawDat/Pcact-PC1315_S3_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/PC13_15/F/.
  cp $RawDat/Pcact-PC1315_S3_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/PC13_15/R/.

  RawDat=/data/seq_data/miseq/2017/RAW/171006_M04465_0050_000000000-B49T7/Data/Intensities/BaseCalls
  mkdir -p raw_dna/paired/P.cactorum/11-40/F
  mkdir -p raw_dna/paired/P.cactorum/11-40/R
  mkdir -p raw_dna/paired/P.cactorum/17-21/F
  mkdir -p raw_dna/paired/P.cactorum/17-21/R
  cp $RawDat/Pcact-CN11-40_S3_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/11-40/F/.
  cp $RawDat/Pcact-CN11-40_S3_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/11-40/R/.
  cp $RawDat/Pcact-CN17-21_S1_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/17-21/F/.
  cp $RawDat/Pcact-CN17-21_S1_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/17-21/R/.

  ProjectDir=/home/groups/harrisonlab/project_files/idris
  RawDat=/data/seq_data/miseq/2018/RAW/180213_M04465_0067_000000000-BJ4DR/Data/Intensities/BaseCalls
  OutDir=raw_dna/paired/P.cactorum/P421/F
  mkdir -p $OutDir
  cd $OutDir
  cp -s $RawDat/P421_S3_L001_R1_001.fastq.gz .
  cd $ProjectDir
  OutDir=raw_dna/paired/P.cactorum/P421/R
  mkdir -p $OutDir
  cd $OutDir
  cp $RawDat/P421_S3_L001_R2_001.fastq.gz .
  cd $ProjectDir
```

#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:


```bash
  # for RawData in $(ls raw_dna/paired/P.*/*/*/*.fastq.gz | grep '170210'); do
  for RawData in $(ls raw_dna/paired/P.cactorum/*/*/*.fastq.gz | grep -v '10300' | grep -e 'P421'); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```


Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf

Firstly, those strains with more than one run were identified:

```bash
for Strain in $(ls -d raw_dna/paired/P.*/* | grep -e 'P.cactorum' -e 'P.idaei'); do
NumReads=$(ls $Strain/F/*.gz | wc -l);
if [ $NumReads -gt 1 ]; then
echo "$Strain";
echo "$NumReads";
fi;
done
```

```bash
for StrainPath in $(ls -d raw_dna/paired/P.*/* | grep -v -w -e '10300' -e '404' -e '414' -e '415' -e '416' -e 'PC13_15' -e '2003_3' | grep -e 'P.cactorum' -e 'P.idaei' | grep -e 'P421'); do
# for StrainPath in $(ls -d raw_dna/paired/*/* | grep -e 'idaei'); do
Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
done		
echo $StrainPath
Read_F=$(ls $StrainPath/F/*.fastq.gz)
Read_R=$(ls $StrainPath/R/*.fastq.gz)
IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

Trimming was then performed for strains with two runs of data

```bash
for StrainPath in $(ls -d raw_dna/paired/P.cactorum/* | grep -w -e '415' -e '416' -e 'PC13_15' -e '2003_3'
| grep '2003_3'); do
Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
done		
echo $StrainPath
Read_F=$(ls $StrainPath/F/*.fastq.gz | head -n1 | tail -n1)
Read_R=$(ls $StrainPath/R/*.fastq.gz | head -n1 | tail -n1)
IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
Read_F=$(ls $StrainPath/F/*.fastq.gz | head -n2 | tail -n1)
Read_R=$(ls $StrainPath/R/*.fastq.gz | head -n2 | tail -n1)
IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

for three sets of data:

```bash
for StrainPath in $(ls -d raw_dna/paired/P.cactorum/* | grep -w -e '404' -e '414'); do
Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'rna_qc_' | grep 'qw' | wc -l)
done		
echo $StrainPath
Read_F=$(ls $StrainPath/F/*.fastq.gz | head -n1 | tail -n1)
Read_R=$(ls $StrainPath/R/*.fastq.gz | head -n1 | tail -n1)
IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
Read_F=$(ls $StrainPath/F/*.fastq.gz | head -n2 | tail -n1)
Read_R=$(ls $StrainPath/R/*.fastq.gz | head -n2 | tail -n1)
IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
Read_F=$(ls $StrainPath/F/*.fastq.gz | head -n3 | tail -n1)
Read_R=$(ls $StrainPath/R/*.fastq.gz | head -n3 | tail -n1)
IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

Data quality was visualised once again following trimming:

```bash
for TrimData in $(ls qc_dna/paired/P.cactorum/*/*/*.fq.gz | grep -e 'P421'); do
echo $TrimData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $TrimData;
done
```

Sequencing coveraqge was estimated:

```bash
for RawData in $(ls qc_dna/paired/P.*/*/*/*q.gz | grep -e 'P421'); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
# GenomeSz=65
GenomeSz=60
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

Find predicted coverage for these isolates:

```bash
for StrainDir in $(ls -d qc_dna/paired/P.*/* | grep -v -w -e '411' -e '10300_old' | grep -e 'cactorum' -e 'idaei' | grep -e '11-40' -e '17-21' -e 'P421'); do
Strain=$(basename $StrainDir)
printf "$Strain\t"
for File in $(ls qc_dna/paired/P.*/"$Strain"/*/*.txt); do
echo $(basename $File);
cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
done
```

```
  10300	206.15
  12420	25.6
  15_13	55.63
  15_7	59.71
  2003_3 69.29
  4032	53.14
  404	69.52
  4040	67.21
  414	76.49
  415	86.64
  416	106.94
  62471	72.21
  P295	39.29
  PC13_15	77.83
  R36_14	56.09
  371	71.89
  SCRP370	40.91
  SCRP376	49.19
  11-40	68.9
  17-21	78.86
  P421	48.22
```

<!--
kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

```bash
  for TrimPath in $(ls -d qc_dna/paired/P.*/* | grep -v -e 'cactorum' | grep 'idaei'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    TrimF=$(ls $TrimPath/F/*.fq.gz)
    TrimR=$(ls $TrimPath/R/*.fq.gz)
    echo $TrimF
    echo $TrimR
    qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimR
  done
```

mode kmer abundance prior to error correction was reported using the following
commands:

```bash
  for File in $(ls qc_dna/kmc/P.cactorum/*/*_true_kmer_summary.txt); do
    basename $File;
    tail -n3 $File | head -n1 ;
  done
```

```
  10300_true_kmer_summary.txt
  The mode kmer abundance is:  255
  404_true_kmer_summary.txt
  The mode kmer abundance is:  14
  414_true_kmer_summary.txt
  The mode kmer abundance is:  14
  415_true_kmer_summary.txt
  The mode kmer abundance is:  18
  416_true_kmer_summary.txt
  The mode kmer abundance is:  25
  62471_true_kmer_summary.txt
  The mode kmer abundance is:  5 <- incorrect thresholding - aprox. 40x
```

-->

#Assembly

Assembly was performed with:
* Spades

## Spades Assembly

Assembly was submitted for genomes with a single run of data

```bash
for StrainPath in $(ls -d qc_dna/paired/P.*/* | grep -v -w -e '10300' -e '404' -e '414' -e '415' -e '416' -e 'PC13_15' -e '2003_3'| grep -e 'P.cactorum' -e 'P.idaei' | grep -e '11-40' -e '17-21' -e 'P421'); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
OutDir=assembly/spades/$Organism/$Strain
Jobs=$(qstat | grep  -e 'subSpades' -e 'submit_SPA'  | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 5m
printf "."
Jobs=$(qstat | grep  -e 'subSpades' -e 'submit_SPA'  | grep 'qw' | wc -l)
done		
printf "\n"
echo $F_Read
echo $R_Read
qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct
done
```

For isolates with two runs of data:

```bash
for StrainPath in $(ls -d qc_dna/paired/P.cactorum/* | grep -w -e '415' -e '416' -e 'PC13_15' -e '2003_3' | grep -e '2003_3'); do
    echo $StrainPath
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1);
    TrimR1_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1);
    TrimF2_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n2 | tail -n1);
    TrimR2_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n2 | tail -n1);
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    OutDir=assembly/spades/$Organism/$Strain
    Jobs=$(qstat | grep -e 'subSpades' -e 'submit_SPA' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
    sleep 5m
    printf "."
    Jobs=$(qstat | grep -e 'subSpades' -e 'submit_SPA' | grep 'qw' | wc -l)
    done		
    printf "\n"
    qsub $ProgDir/subSpades_2lib_HiMem.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir correct
  done
```



for isolates with three runs of data:

```bash
for StrainPath in $(ls -d qc_dna/paired/P.cactorum/* | grep -w -e '404' -e '414'); do
    echo $StrainPath
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1);
    TrimR1_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1);
    TrimF2_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n2 | tail -n1);
    TrimR2_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n2 | tail -n1);
    TrimF3_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n3 | tail -n1);
    TrimR3_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n3 | tail -n1);
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    echo $TrimF3_Read
    echo $TrimR3_Read
    OutDir=assembly/spades/$Organism/$Strain
    Jobs=$(qstat | grep -e 'subSpades' -e 'submit_SPA' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
    sleep 5m
    printf "."
    Jobs=$(qstat | grep -e 'subSpades' -e 'submit_SPA' | grep 'qw' | wc -l)
    done		
    printf "\n"
    qsub $ProgDir/subSpades_3lib_HiMem.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir correct
  done
```


Assembly results were summarised using Quast:

```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs*/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -e '11-40' -e '17-21' -e 'P421'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    OutDir=$(dirname $Assembly)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir;
  done
```

Contigs were identified that had blast hits to non-phytophthora genomes
<!--
```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    # Exclude_db="bact,virus,hsref"
    # Exclude_db="paenibacillus"
    # Exclude_db="common_contaminants"
    Exclude_db="paenibacillus,bacillus,delftia"
    Good_db="phytoph"
    AssemblyDir=$(dirname $Assembly)
    # OutDir=$AssemblyDir/../deconseq_Paen
    # OutDir=$AssemblyDir/../deconseq_common
    OutDir=$AssemblyDir/../deconseq_common_seperated
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
  done
``` -->

```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -e '11-40' -e '17-21' -e 'P421' | grep 'P421_v2'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for Exclude_db in "bacillus" "delftia" "paenibacillus"; do
      Good_db="phytoph"
      AssemblyDir=$(dirname $Assembly)
      OutDir=$AssemblyDir/../deconseq_$Exclude_db
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
      qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
    done
  done
```

Results were summarised using the commands:

```bash
# for File in $(ls assembly/spades/P.*/*/deconseq/log.txt); do
# for File in $(ls assembly/spades/P.*/*/deconseq_common/log.txt); do
for Exclude_db in "bacillus" "delftia" "paenibacillus"; do
echo $Exclude_db
for File in $(ls assembly/spades/P.*/*/*/log.txt | grep "$Exclude_db" | grep -e '11-40' -e '17-21' -e 'P421'); do
Name=$(echo $File | rev | cut -f3 -d '/' | rev);
Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
Both=$(cat $File |cut -f2 | head -n2 | tail -n1);
Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
printf "$Name\t$Good\t$Both\t$Bad\n";
done
done
```
<!-- ```
  12420	5429	83	29
  15_13	5332	74	16
  15_7	5363	82	20
  2003_3	5558	68	31
  2003_4	5993	101	34
  4032	8260	91	51
  404	20003	70	121
  4040	5797	97	48
  414	5921	69	43
  415	5595	83	38
  416	5376	93	41
  62471	6978	100	48
  P295	5093	57	39
  PC13_15	5200	84	33
  R36_14	5854	102	42
  371	4679	27	5
  SCRP370	5303	38	17
  SCRP376	5228	70	17
``` -->
Results from Paen
```
12420	5508	3	30
15_13	5419	2	1
15_7	5457	1	7
2003_3	5473	5	111
4032	6705	28	1669
404	20143	2	49
4040	5324	3	615
414	5107	2	924
415	5707	2	7
416	5501	2	7
62471	6156	24	946
P295	5088	2	99
PC13_15	5275	3	39
R36_14	5957	3	38
371	4709	1	1
SCRP370	5355	2	1
SCRP376	5313	1	1
11-40	4573	1	1
17-21	4451	2	1
P421	8078	2	3
```
Results from delftia
```
12420	5534	1	6
15_13	5420	1	1
15_7	5462	0	3
2003_3	5583	0	6
4032	8396	2	4
4040	5938	2	2
404	20184	0	10
414	6029	0	4
415	5713	0	3
416	5507	0	3
62471	7123	1	2
P295	5179	1	9
PC13_15	5309	0	8
R36_14	5990	0	8
371	4709	1	1
SCRP370	5355	2	1
SCRP376	5313	1	1
11-40	4573	1	1
17-21	4452	1	1
P421	8065	2	16
```
Results from bacillus
```
12420	5529	3	9
15_13	5419	2	1
15_7	5461	1	3
2003_3	5577	1	11
4032	8330	2	70
4040	5890	5	47
404	20175	2	17
414	5995	1	37
415	5693	3	20
416	5486	2	22
62471	7056	5	65
P295	5162	3	24
PC13_15	5302	1	14
R36_14	5975	1	22
371	4710	0	1
SCRP370	5356	1	1
SCRP376	5314	0	1
bacillus
11-40	4573	1	1
17-21	4450	3	1
17-21	4451	2	1
P421	8077	3	3
P421	8078	2	3
```
<!--
Results from common conatmainants
```
12420	5513	6	22
15_13	5416	4	2
15_7	5457	4	4
2003_3	5546	5	38
4032	7969	10	423
4040	5646	15	281
404	20145	5	44
414	5812	7	214
415	5691	5	20
416	5483	5	22
62471	6788	5	333
P295	5121	6	62
PC13_15	5283	3	31
R36_14	5961	2	35
371	4709	1	1
SCRP370	5355	2	1
SCRP376	5312	1	2
``` -->

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -e '11-40' -e '17-21' -e 'P421' | grep 'P421_v2'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
for Exclude_db in "stenotrophomonas"; do
Good_db="phytoph"
AssemblyDir=$(dirname $Assembly)
OutDir=$AssemblyDir/../deconseq_$Exclude_db
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
done
done
```

Those contigs that were not identified as contaminants in each of the filtering
steps were retained.

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei'| grep -e '11-40' -e '17-21' -e 'P421' | grep -v -w 'P421'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
StrainDir=$(ls -d assembly/spades/$Organism/$Strain)
mkdir -p $StrainDir/deconseq_appended
for File in $(ls $StrainDir/deconseq_*/*cont.fa | grep -e "bacillus" -e "delftia" -e "Paen"); do
cat $File | grep '>'
done | sort | uniq | tr -d '>' > $StrainDir/deconseq_appended/exclude_list.txt
Instructions=$StrainDir/deconseq_appended/exclude_instructions.txt
printf "Exclude:\nSequence name, length, apparent source\n" > $Instructions
cat $StrainDir/deconseq_appended/exclude_list.txt | sed -r 's/$/\t.\tcontaminant/g' >> $Instructions
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $StrainDir/deconseq_appended/contigs_min_500bp_renamed.fasta --coord_file $Instructions
done
```

Strain P421 was searched against an additional database.

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei' |  grep 'P421'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
StrainDir=$(ls -d assembly/spades/$Organism/$Strain)
mkdir -p $StrainDir/deconseq_appended
for File in $(ls $StrainDir/deconseq_*/*cont.fa | grep -e "bacillus" -e "delftia" -e "Paen" -e "stenotrophomonas"); do
cat $File | grep '>'
done | sort | uniq | tr -d '>' > $StrainDir/deconseq_appended/exclude_list.txt
Instructions=$StrainDir/deconseq_appended/exclude_instructions.txt
printf "Exclude:\nSequence name, length, apparent source\n" > $Instructions
cat $StrainDir/deconseq_appended/exclude_list.txt | sed -r 's/$/\t.\tcontaminant/g' >> $Instructions
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $StrainDir/deconseq_appended/contigs_min_500bp_renamed.fasta --coord_file $Instructions
done
```

The number of reads aligning from P421 to potential contaminants were identified.


```bash
Fusarium=$(ls ../fusarium/assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa)
Stenotrophomonas=$(ls ../../../../../home/armita/prog/deconseq-standalone-0.4.3/database/stenotrophomonas/all_complete_GbStenotrophomonas.fasta)
for Reference in $(ls $Fusarium $Stenotrophomonas); do
for StrainPath in $(ls -d qc_dna/paired/*/* | grep '421'); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
Prefix=$(basename $Reference | sed 's/.fasta//g' | sed 's/.fa//g' | sed 's/.fna//g')
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_${Prefix}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie_unaligned.sh $Reference $F_Read $R_Read $OutDir
done
done
```

P241 was reassembled:

```bash
for StrainPath in $(ls -d analysis/genome_alignment/bowtie/P.*/* | grep 'P421'); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/vs_all_complete_GbStenotrophomonas/*.fq.1.gz)
R_Read=$(ls $StrainPath/vs_all_complete_GbStenotrophomonas/*.fq.2.gz)
OutDir=assembly/spades/$Organism/${Strain}_v2
echo $F_Read
echo $R_Read
qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct
done
```

Assembly results were summarised using Quast:

```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs*/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep 'P421_v2'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    OutDir=$(dirname $Assembly)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir;
  done
```


```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep 'P421_v2'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# for Exclude_db in "bacillus" "delftia" "paenibacillus" "stenotrophomonas" "FoL"; do
for Exclude_db in "FoL"; do
Good_db="phytoph"
AssemblyDir=$(dirname $Assembly)
OutDir=$AssemblyDir/../deconseq_$Exclude_db
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
done
done
```

Results were summarised using the commands:

```bash
for Exclude_db in "bacillus" "delftia" "paenibacillus" "stenotrophomonas" "FoL"; do
echo $Exclude_db
for File in $(ls assembly/spades/P.*/*/*/log.txt | grep "_${Exclude_db}" | grep 'P421_v2'); do
Name=$(echo $File | rev | cut -f3 -d '/' | rev);
Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
Both=$(cat $File |cut -f2 | head -n2 | tail -n1);
Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
printf "$Name\t$Good\t$Both\t$Bad\n";
done
done
```

```
bacillus
P421_v2	7706	2	8
delftia
P421_v2	7536	2	178
paenibacillus
P421_v2	7699	2	15
stenotrophomonas
P421_v2	7058	5	653
FoL
P421_v2	7189	45	482
```

Those contigs that were not identified as contaminants in each of the filtering
steps were retained.

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -w 'P421_v2'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
StrainDir=$(ls -d assembly/spades/$Organism/$Strain)
mkdir -p $StrainDir/deconseq_appended
for File in $(ls $StrainDir/deconseq_*/*cont.fa | grep -e "bacillus" -e "delftia" -e "Paen" -e "stenotrophomonas" -e "FoL"); do
cat $File | grep '>'
done | sort | uniq | tr -d '>' > $StrainDir/deconseq_appended/exclude_list.txt
Instructions=$StrainDir/deconseq_appended/exclude_instructions.txt
printf "Exclude:\nSequence name, length, apparent source\n" > $Instructions
cat $StrainDir/deconseq_appended/exclude_list.txt | sed -r 's/$/\t.\tcontaminant/g' >> $Instructions
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $StrainDir/deconseq_appended/contigs_min_500bp_renamed.fasta --coord_file $Instructions
done
```

Assembly stats were collected on filtered assemblies:

```bash
for Assembly in $(ls assembly/spades/P.*/*/deconseq_appended/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei'| grep 'P421_v2'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

Assembly stats were summarised and compared to previous assembly results using:

```bash
# for Assembly in $(ls assembly/spades/P.*/*/deconseq_Paen/report.tsv); do  
for Assembly in $(ls assembly/spades/P.*/*/deconseq_appended/report.tsv | grep -w '404'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/'| rev);
Size=$(cat $Assembly | grep 'Total length' | head -n1 | cut -f2);
OldAssembly=$(ls assembly/spades/P.*/$Strain/filtered_contigs*/report.tsv)
OldSize=$(cat $OldAssembly | grep 'Total length' | head -n1 | cut -f2);
printf "$Strain\t$Size\t$OldSize\n";
done
```

```
  11-40	59726247	59731760
  12420	60315503	66488881
  15_13	60703463	60708976
  15_7	60797060	60806527
  17-21	59761502	59767015
  2003_3	60800281	67797455
  4032	61561706	70105079
  4040	60657148	65177765
  404	75508939	82472226
  414	61116855	62710907
  415	60677782	64884339
  416	60502157	64802670
  62471	61265225	62823426
  P295	60707114	67293243
  P421	64527586	65473778
  P421_v2	60335997	64856938
  PC13_15	61130682	67857120
  R36_14	61656467	68350251
  371	60474213	60483665
  SCRP370	60582725	60561664
  SCRP376	60610468	60615981
```

Numbers of busco genes in each assembly were identified:

```bash
for Assembly in $(ls assembly/spades/P.*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
# for Assembly in $(ls assembly/spades/P.*/*/deconseq_appended/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei'| grep 'P421_v2'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB="Eukaryotic"
# OutDir=gene_pred/busco/$Organism/$Strain/assembly
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
<!--
```bash
  # for File in $(ls assembly/spades/*/*/deconseq_Paen/run_contigs_min_500bp_filtered_renamed/short_summary_contigs_min_500bp_filtered_renamed.txt  | grep -e 'P.cactorum' -e 'P.idaei'); do
  for File in $(ls assembly/spades/*/*/deconseq_common/run_contigs_min_500bp_filtered_renamed/short_summary_contigs_min_500bp_filtered_renamed.txt  | grep -e 'P.cactorum' -e 'P.idaei'); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
``` -->



A Bioproject and Biosample was made with NCBI genbank for submission of genomes.
Following the creation of these submissions, the .fasta assembly was uploaded
through the submission portal. A note was provided requesting that the assembly
be run through the contamination screen to aid a more detailed resubmission in
future. The returned FCSreport.txt was downloaded from the NCBI webportal and
used to correct the assembly to NCBI standards.

NCBI reports (FCSreport.txt) were manually downloaded to the following loactions:

```bash
for Assembly in $(ls assembly/spades/P.*/*/deconseq_appended/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep 'P421_v2'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
echo "$PWD/$NCBI_report_dir"
mkdir -p $NCBI_report_dir
done
```

The 12-420 genome was further edited by removing contig1
```bash
Strain='12420'
Organism='P.cactorum'
NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_submission/Contamination*.txt)
sed -i 's/Trim:/Exclude:/g' $NCBI_report
```


These downloaded files were used to correct assemblies:

```bash
for Assembly in $(ls assembly/spades/P.*/*/deconseq_appended/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -e 'P421_v2'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_submission/Contamination*.txt)
if [[ $NCBI_report ]]; then
echo "Contamination report found"
else
NCBI_report=genome_submission/$Organism/$Strain/initial_submission/no_edits.txt
printf "Exclude:\nSequence name, length, apparent source\n" > $NCBI_report
fi
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
mkdir -p $OutDir
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
# $ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
done
```

```bash
echo "Strain P421 has further edits to be made"
NCBI_report_dir=genome_submission/P.cactorum/P421_v2/second_submission
echo "$PWD/$NCBI_report_dir"
mkdir -p $NCBI_report_dir
# Move downloaded contamination report into this directory
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep 'P421_v2'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Round1Name=$(echo $Assembly | sed 's/contigs_min_500bp_renamed.fasta/round1.fasta/g')
mv $Assembly $Round1Name
NCBI_report=genome_submission/$Organism/$Strain/second_submission/Contamination*.txt
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
mkdir -p $OutDir
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --inp $Round1Name --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
done
#
# NCBI_report_dir=genome_submission/P.cactorum/P421_v2/third_submission
# echo "$PWD/$NCBI_report_dir"
# mkdir -p $NCBI_report_dir
# # Move downloaded contamination report into this directory
# for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep 'P421_v2'); do
# Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
# Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
# echo "$Organism - $Strain"
# Round2Name=$(echo $Assembly | sed 's/contigs_min_500bp_renamed.fasta/round2.fasta/g')
# mv $Assembly $Round2Name
# NCBI_report=genome_submission/$Organism/$Strain/third_submission/Contamination*.txt
# OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
# mkdir -p $OutDir
# ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
# $ProgDir/remove_contaminants.py --inp $Round2Name --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
# done
```



```bash
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep 'P421_v2'); do
    Kmer=$(echo $Assembly | rev | cut -f2 -d '/' | rev);
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    OutDir=$(dirname $Assembly)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir;
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB="Eukaryotic"
    OutDir=gene_pred/busco/$Organism/$Strain/assembly
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

Assembly stats were summarised using:
```bash
  for File in $(ls assembly/spades/P.cactorum/*/filtered_contigs/report.tsv); do
    echo "$File" | rev | cut -f3 -d '/' | rev ;
    cat $File | cut -f2;
  done
```

The following assemblies were accepted by ncbi following runs through the Paenobacillus deconseq database:

```
Pcac_12-420
Pcac_15-13
Pcac_4032
Pcac_4040
```
The following were run through the Paenobacillus, Bacillus, Delftia database:
```
R36-14
416
415
2003_3
```

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/spades/P.*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -e 'P.idaei' -e 'P.cactorum' | grep 'P421'); do
Strain=$(echo $BestAss | rev | cut -f3 -d '/' | rev)
Organism=$(echo $BestAss | rev | cut -f4 -d '/' | rev)
OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs_repmask
qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
for RepDir in $(ls -d repeat_masked/P.*/*/filtered_contigs_repmask | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
# printf "The number of bases masked by RepeatMasker:\t"
Repmasker=$(sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
# printf "The number of bases masked by TransposonPSI:\t"
TPSI=$(sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
# printf "The total number of masked bases are:\t"
Total=$(cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
printf "$Organism\t$Strain\t$Repmasker\t$TPSI\t$Total\n"
done
```

```
P.cactorum	11-40	13618425	3677249	14636093
P.cactorum	12420	13952648	3830927	14889294
P.cactorum	15_13	14871198	4311252	15878513
P.cactorum	15_7	14842223	4444413	15807291
P.cactorum	17-21	14153794	3637012	15026432
P.cactorum	2003_3	14869062	4312476	15913172
P.cactorum	4032	15067487	4189364	16046834
P.cactorum	4040	14914633	4118840	15860862
P.cactorum	404	12646327	4410591	13909457
P.cactorum	415	15088727	4308075	16093167
P.cactorum	416	14782441	4199154	15719961
P.cactorum	62471	14741873	4113154	15684662
P.cactorum	P295	14635224	4005338	15561536
P.cactorum	P421_v2	14071987	3526248	14909132
P.cactorum	PC13_15	14931997	4176759	15852338
P.cactorum	R36_14	15206511	4318141	16169407
P.idaei	371	15109356	4134856	16074815
P.idaei	SCRP370	15417500	4122740	16384703
P.idaei	SCRP376	15516157	4111645	16516659
```

# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

# Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

Busco has replaced CEGMA and was run to check gene space in assemblies

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
# Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
# Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB="Eukaryotic"
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/run_*_contigs_unmasked/short_summary_*.txt | grep -e 'P.cactorum' -e 'P.idaei'); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

```
P.cactorum	10300	283	3	17	303
P.cactorum	11-40	285	3	15	303
P.cactorum	12420	285	2	16	303
P.cactorum	15_13	286	2	15	303
P.cactorum	15_7	286	2	15	303
P.cactorum	17-21	286	2	15	303
P.cactorum	2003_3	286	2	15	303
P.cactorum	4032	285	3	15	303
P.cactorum	4040	286	2	15	303
P.cactorum	404	287	3	13	303
P.cactorum	414	285	2	16	303
P.cactorum	415	285	2	16	303
P.cactorum	416	286	2	15	303
P.cactorum	62471	285	3	15	303
P.cactorum	P295	285	3	15	303
P.cactorum	P421	285	3	15	303
P.cactorum	PC13_15	286	2	15	303
P.cactorum	R36_14	285	3	15	303
P.idaei	371	285	3	15	303
P.idaei	SCRP370	285	3	15	303
P.idaei	SCRP376	285	3	15	303
```


```bash
  for File in $(ls repeat_masked/*/*/filtered_contigs_repmask/report.tsv | grep -e 'P.cactorum' -e 'P.idaei'); do
  Strain=$(echo $File| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
  Contigs=$(cat $File | grep "contigs (>= 0 bp)" | cut -f2)
  Length=$(cat $File | grep "Total length (>= 0 bp)" | cut -f2)
  Largest=$(cat $File | grep "Largest contig" | cut -f2)
  N50=$(cat $File | grep "N50" | cut -f2)
  echo -e "$Organism\t$Strain\t$Contigs\t$Length\t$Largest\t$N50"
  done
```

```
P.cactorum	11-40	4571	59724358	288536	43828
P.cactorum	12420	5509	59887557	272591	37498
P.cactorum	15_13	5421	60703463	358283	41248
P.cactorum	15_7	5457	60796470	300516	43376
P.cactorum	17-21	4452	59759428	347645	49722
P.cactorum	2003_3	5467	60795145	353051	40273
P.cactorum	4032	6726	61561706	299428	38202
P.cactorum	4040	5326	60657148	238735	39295
P.cactorum	404	20136	75501778	346862	36598
P.cactorum	415	5695	60677265	313236	40742
P.cactorum	416	5485	60476717	251913	40465
P.cactorum	62471	6173	61265185	285881	36408
P.cactorum	P295	5088	60707037	357485	48651
P.cactorum	P421_v2	6581	59952001	239266	34103
P.cactorum	PC13_15	5278	61128622	346832	39779
P.cactorum	R36_14	5952	61655841	345260	36465
P.idaei	371	4720	60469986	295942	46381
P.idaei	SCRP370	5356	60582035	247223	39576
P.idaei	SCRP376	5313	60609233	230545	39406
```

#Gene prediction

Gene prediction was performed for P. cactorum genomes. Two gene prediction
approaches were used:

Gene prediction using Braker1
Prediction of all putative ORFs in the genome using the ORF finder (atg.pl)
approach.


## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to Fusarium genomes.
* qc of RNA seq data was performed as part of sequencing the 10300 genome:


### Aligning in house RNAseq data

#### Timecourse of infection

Maria performed alignment of RNAseq data vs the F. vesca and those reads that did not align were
used for alignment vs the P.cactorum genome.

make symbolic links to timecourse data

```bash
  for File in $(ls /home/sobczm/popgen/rnaseq/vesca_*/*.mate1.fq.gz); do
    Sample=$(echo $File | rev | cut -d '/' -f2 | rev | sed 's/vesca_//g')
    echo "$Sample"
    OutDir=qc_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07_no_vesca/$Sample/F
    mkdir -p "$OutDir"
    cp -s $File $OutDir/.
  done
  for File in $(ls /home/sobczm/popgen/rnaseq/vesca_*/*.mate2.fq.gz); do
    Sample=$(echo $File | rev | cut -d '/' -f2 | rev | sed 's/vesca_//g')
    echo "$Sample"
    OutDir=qc_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07_no_vesca/$Sample/R
    mkdir -p "$OutDir"
    cp -s $File $OutDir/.
  done
```

#### Mycelium

Qc of RNAseq data from mycelium and in-plant infection is detailed in commands for gene prediction of the PacBio assmbly under:
assembly/pacbio_assembly.md

Performed alignment of RNAseq data vs the F. vesca genome and those reads that did not align were
used for alignment vs the P.cactorum genome.

#### Aligning

Aligning mycelium data

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d ../../../../../data/scratch/armita/idris/alignment/star/fvesca/v1.1/*/* | grep 'mycelium'); do
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
done
printf "\n"
FileF=$(ls $RNADir/*.mate1.fq.gz)
FileR=$(ls $RNADir/*.mate2.fq.gz)
echo $FileF
echo $FileR
Prefix=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
Timepoint=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
echo "$Timepoint"
OutDir=../../../../../data/scratch/armita/idris/alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
done
```

Aligning timecourse data

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d qc_rna/paired/*/* | grep -e '_no_vesca'); do
FileNum=$(ls $RNADir/F/*.mate1.fq.gz | wc -l)
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
for num in $(seq 1 $FileNum); do
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
done
printf "\n"
FileF=$(ls $RNADir/F/*.mate1.fq.gz | head -n $num | tail -n1)
FileR=$(ls $RNADir/R/*.mate2.fq.gz | head -n $num | tail -n1)
echo $FileF
echo $FileR
Prefix=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint"
OutDir=../../../../../data/scratch/armita/idris/alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
done
done
```
<!--
Alignment stats were collected:

```bash
for File in $(ls alignment/star/P.cactorum/414/*/*/star_aligmentLog.final.out); do
Sample=$(echo $File | rev | cut -f2 -d '/' | rev);
ReadNumU=$(cat $File | grep 'Uniquely' | grep 'number' | cut -f2);
ReadPercU=$(cat $File | grep 'Uniquely' | grep '%' | cut -f2);
ReadNumM=$(cat $File | grep 'multiple' | grep 'Number' | cut -f2);
ReadPercM=$(cat $File | grep 'multiple' | grep '%' | cut -f2);
echo -e "$Sample""\t""$ReadNumU""\t""$ReadPercU""\t""$ReadNumM""\t""$ReadPercM";  
done
```

```

``` -->

Alignments were concatenated prior to gene prediction

```bash
qlogin -pe smp 8
cd /home/groups/harrisonlab/project_files/idris
for StrainDir in $(ls -d ../../../../../data/scratch/armita/idris/alignment/star/P.*/*); do
echo $StrainDir
BamFiles=$(ls $StrainDir/*/*/star_aligmentAligned.sortedByCoord.out.bam | grep -v -e 'PRO1467_S1_' -e 'PRO1467_S2_' -e 'PRO1467_S3_' -e 'PRO1467_S10_' -e 'PRO1467_S11_' -e 'PRO1467_S12_' | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=$StrainDir/concatenated
mkdir -p $OutDir
samtools merge -@ 8 -f $OutDir/concatenated.bam $BamFiles
done
logout
```


#### Braker prediction

Before braker predictiction was performed, I double checked that I had the
genemark key in my user area and copied it over from the genemark install
directory:

```bash
	ls ~/.gm_key
	cp /home/armita/prog/genemark/gm_key_64 ~/.gm_key
```

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    	echo "$Organism - $Strain"
    	mkdir -p alignment/$Organism/$Strain/concatenated
    	OutDir=gene_pred/braker/$Organism/"$Strain"_braker
    	AcceptedHits=$(ls ../../../../../data/scratch/armita/idris/alignment/$Organism/$Strain/concatenated/concatenated.bam)
    	GeneModelName="$Organism"_"$Strain"_braker
    	rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
    	qsub $ProgDir/sub_braker.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  	done
```


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    Jobs=$(qstat | grep 'sub_cuff' | grep 'qw'| wc -l)
    while [ $Jobs -gt 1 ]; do
    sleep 1m
    printf "."
    Jobs=$(qstat | grep 'sub_cuff' | grep 'qw'| wc -l)
    done
    printf "\n"
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
    mkdir -p $OutDir
    AcceptedHits=$(ls ../../../../../data/scratch/armita/idris/alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
  done
```

Secondly, genes were predicted using CodingQuary:

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | grep 'idaei'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Jobs=$(qstat | grep 'sub_Cod' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_Cod' | grep 'qw'| wc -l)
done
printf "\n"
OutDir=gene_pred/codingquary/$Organism/$Strain
CufflinksGTF=$(ls gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
  for BrakerGff in $(ls gene_pred/braker/P.*/*_braker/*/augustus.gff3 | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g' | sed 's/_braker//g')
    Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/"$Strain"_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
    PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
    AddDir=gene_pred/codingquary/$Organism/$Strain/additional
    FinalDir=gene_pred/codingquary/$Organism/$Strain/final
    AddGenesList=$AddDir/additional_genes.txt
    AddGenesGff=$AddDir/additional_genes.gff
    FinalGff=$AddDir/combined_genes.gff
    mkdir -p $AddDir
    mkdir -p $FinalDir

    bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
    bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
    $ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    # -
    # This section is edited
    $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuary_unspliced.gff3
    $ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuary_unspliced.gff3 > $FinalDir/final_genes_CodingQuary.gff3
    # -
    $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
    cp $BrakerGff $FinalDir/final_genes_Braker.gff3
    $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
    cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
    cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
    cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
    cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

    GffBraker=$FinalDir/final_genes_CodingQuary.gff3
    GffQuary=$FinalDir/final_genes_Braker.gff3
    GffAppended=$FinalDir/final_genes_appended.gff3
    cat $GffBraker $GffQuary > $GffAppended
  done
```

The final number of genes per isolate was observed using:
```bash
  for DirPath in $(ls -d gene_pred/codingquary/P.*/*/final | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    Strain=$(echo $DirPath| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $DirPath | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
    echo "";
  done
```

The number of genes predicted by Braker, supplimented by CodingQuary and in the
final combined dataset was shown:

```
P.cactorum - 11-40
18162
10665
28827

P.cactorum - 12420
18079
8152
26231

P.cactorum - 15_13
18008
8736
26744

P.cactorum - 15_7
17719
8885
26604

P.cactorum - 17-21
17859
10718
28577

P.cactorum - 2003_3
18241
11545
29786

P.cactorum - 4032
17790
9029
26819

P.cactorum - 4040
18098
11432
29530

P.cactorum - 404
23148
12618
35766

P.cactorum - 415
17932
8864
26796

P.cactorum - 416
17883
11508
29391

P.cactorum - 62471
18214
11382
29596

P.cactorum - P295
17776
8551
26327

P.cactorum - P421_v2
17819
8305
26124

P.cactorum - PC13_15
18147
8460
26607

P.cactorum - R36_14
18432
11497
29929

P.idaei - 371
17190
8197
25387

P.idaei - SCRP370
17291
10591
27882

P.idaei - SCRP376
17429
8061
25490
```


## Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing
Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	for Genome in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    echo "$Genome"
  	qsub $ProgDir/run_ORF_finder.sh $Genome
  done
```

The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
	ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
  for OrfGff in $(ls gene_pred/ORF_finder/P.*/*/*_ORF.gff | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    echo "$OrfGff"
  	OrfGffMod=$(echo $OrfGff | sed 's/.gff/.gff3/g')
  	$ProgDir/gff_corrector.pl $OrfGff > $OrfGffMod
  done
```





#Genomic analysis

## RxLR genes

Putative RxLR genes were identified within Augustus gene models using a number
of approaches:

 * A) From Augustus gene models - Signal peptide & RxLR motif  
 * B) From Augustus gene models - Hmm evidence of WY domains  
 * C) From Augustus gene models - Hmm evidence of RxLR effectors
 * D) From Augustus gene models - Hmm evidence of CRN effectors  
 <!-- * E) From ORF fragments - Signal peptide & RxLR motif  
 * F) From ORF fragments - Hmm evidence of WY domains  
 * G) From ORF fragments - Hmm evidence of RxLR effectors -->


 ### A) From Augustus gene models - Signal peptide & RxLR motif

 Required programs:
  * SigP
  * biopython

#### A.1) Signal peptide prediction using SignalP 2.0

 Proteins that were predicted to contain signal peptides were identified using
 the following commands:

 ```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | grep -v 'PC13_15'); do
    SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SplitDir=gene_pred/final_split/$Organism/$Strain
    mkdir -p $SplitDir
    BaseName="$Organism""_$Strain"_final
    $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
    for File in $(ls $SplitDir/*_final_*); do
      Jobs=$(qstat | grep 'pred_sigP' | grep -w 'qw'| wc -l)
      while [ $Jobs -gt 4 ]; do
        sleep 1
        printf "."
        Jobs=$(qstat | grep 'pred_sigP' | grep -w 'qw'| wc -l)
      done
      printf "\n"
      echo $File
      # qsub $ProgDir/pred_sigP.sh $File
      qsub $ProgDir/pred_sigP.sh $File signalp-3.0
      # qsub $ProgDir/pred_sigP.sh $File signalp-4.1
    done
  done
 ```

Problems with SignalP3 prediction of isolate PC13_15 were found to be associate
with the length of the isolate name. The name was shortened for prediction and
will be renamed again when results are concatenated.

 ```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | grep 'PC13'); do
    SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev | sed 's/_15//g')
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SplitDir=gene_pred/final_split/$Organism/$Strain
    mkdir -p $SplitDir
    BaseName="$Organism""_$Strain"_final
    $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
    for File in $(ls $SplitDir/*_final_*); do
      Jobs=$(qstat | grep 'pred_sigP' | grep -w 'qw'| wc -l)
      while [ $Jobs -gt 4 ]; do
        sleep 1
        printf "."
        Jobs=$(qstat | grep 'pred_sigP' | grep -w 'qw'| wc -l)
      done
      printf "\n"
      echo $File
      qsub $ProgDir/pred_sigP.sh $File
      qsub $ProgDir/pred_sigP.sh $File signalp-3.0
      qsub $ProgDir/pred_sigP.sh $File signalp-4.1
    done
  done
 ```


 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

 ```bash
  for SplitDir in $(ls -d gene_pred/final_split/P.*/* | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | grep -v PC13); do
    Strain=$(echo $SplitDir | cut -d '/' -f4)
    Organism=$(echo $SplitDir | cut -d '/' -f3)
    echo "$Organism - $Strain"
    InStringAA=''
    InStringNeg=''
    InStringTab=''
    InStringTxt=''
    for SigpDir in $(ls -d gene_pred/final_sig* | cut -f2 -d'/'); do
      for GRP in $(ls -l $SplitDir/*_final_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
        InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.aa";
        InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp_neg.aa";  
        InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.tab";
        InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.txt";
      done
      cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
      cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
      tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
      cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
    done
  done
```


 ```bash
for SplitDir in $(ls -d gene_pred/final_split/P.*/* | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | grep -w 'PC13'); do
Strain=$(echo $SplitDir | cut -d '/' -f4)
Organism=$(echo $SplitDir | cut -d '/' -f3)
echo "$Organism - $Strain"
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
for SigpDir in $(ls -d gene_pred/final_sig* | cut -f2 -d'/'); do
for GRP in $(ls -l $SplitDir/*_final_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.aa";
InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp_neg.aa";  
InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.tab";
InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.txt";
done
cat $InStringAA > gene_pred/$SigpDir/$Organism/${Strain}_15/${Strain}_15_aug_sp.aa
cat $InStringNeg > gene_pred/$SigpDir/$Organism/${Strain}_15/${Strain}_15_aug_neg_sp.aa
tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/${Strain}_15/${Strain}_15_aug_sp.tab
cat $InStringTxt > gene_pred/$SigpDir/$Organism/${Strain}_15/${Strain}_15aug_sp.txt
done
done
```

#### B.2) Prediction using Phobius

 Secreted proteins were also predicted using Phobius

 ```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/phobius/$Organism/$Strain
    mkdir -p $OutDir
    phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    $ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
    cat $OutDir/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d '>' > $OutDir/"$Strain"_phobius_headers.txt
  done
 ```


Secreted proteins from different sources were combined into a single file:

```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/combined_sigP/$Organism/$Strain
    mkdir -p $OutDir
    echo "The following number of sequences were predicted as secreted:"
    cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
    cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | wc -l
    echo "This represented the following number of unique genes:"
    cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
    cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
  done
```

```
P.cactorum - 11-40
The following number of sequences were predicted as secreted:
16631
This represented the following number of unique genes:
3273
P.cactorum - 12420
The following number of sequences were predicted as secreted:
16145
This represented the following number of unique genes:
3139
P.cactorum - 15_13
The following number of sequences were predicted as secreted:
16187
This represented the following number of unique genes:
3142
P.cactorum - 15_7
The following number of sequences were predicted as secreted:
16028
This represented the following number of unique genes:
3110
P.cactorum - 17-21
The following number of sequences were predicted as secreted:
17022
This represented the following number of unique genes:
3320
P.cactorum - 2003_3
The following number of sequences were predicted as secreted:
16784
This represented the following number of unique genes:
3289
P.cactorum - 4032
The following number of sequences were predicted as secreted:
16284
This represented the following number of unique genes:
3172
P.cactorum - 4040
The following number of sequences were predicted as secreted:
16816
This represented the following number of unique genes:
3310
P.cactorum - 404
The following number of sequences were predicted as secreted:
19599
This represented the following number of unique genes:
3924
P.cactorum - 415
The following number of sequences were predicted as secreted:
17325
This represented the following number of unique genes:
3166
P.cactorum - 416
The following number of sequences were predicted as secreted:
17468
This represented the following number of unique genes:
3283
P.cactorum - 62471
The following number of sequences were predicted as secreted:
17048
This represented the following number of unique genes:
3354
P.cactorum - P295
The following number of sequences were predicted as secreted:
16146
This represented the following number of unique genes:
3117
P.cactorum - P421_v2
The following number of sequences were predicted as secreted:
7065
This represented the following number of unique genes:
2972
P.cactorum - PC13_15
The following number of sequences were predicted as secreted:
16328
This represented the following number of unique genes:
3168
P.cactorum - R36_14
The following number of sequences were predicted as secreted:
17028
This represented the following number of unique genes:
3379
P.idaei - 371
The following number of sequences were predicted as secreted:
15536
This represented the following number of unique genes:
2918
P.idaei - SCRP370
The following number of sequences were predicted as secreted:
15635
This represented the following number of unique genes:
3104
P.idaei - SCRP376
The following number of sequences were predicted as secreted:
14942
This represented the following number of unique genes:
2936
```


Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
    qsub $ProgDir/submit_TMHMM.sh $Proteome
  done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt | grep -v -e '414' -e '10300'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
NonTmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $NonTmHeaders
SigP=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_secreted.fa)
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $NonTmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
echo "Number of SigP proteins:"
cat $SigP | grep '>' | wc -l
echo "Number without transmembrane domains:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
echo "Number of gene models:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l

# A text file was also made containing headers of proteins testing +ve
PosFile=$(ls gene_pred/trans_mem/$Organism/$Strain/"$Strain"_TM_genes_pos.txt)
TmHeaders=$(echo $PosFile | sed 's/.txt/_headers.txt/g')
cat $PosFile | cut -f1 > $TmHeaders

done
```

<!--

### C) From Augustus gene models - Effector-like structure identification using EffectorP

Required programs:
 * EffectorP.py

```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    BaseName="$Organism"_"$Strain"_EffectorP
    OutDir=analysis/effectorP/$Organism/$Strain
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
    qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
  done
```

Those genes that were predicted as secreted and tested positive by effectorP
were identified:

```bash
  for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
    cat $File | grep 'Effector' | cut -f1 > $Headers
    printf "EffectorP headers:\t"
    cat $Headers | wc -l
    Secretome=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_secreted.fa)
    {Strain}=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > ${Strain}
    {Strain}Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat ${Strain} | grep '>' | tr -d '>' > ${Strain}Headers
    printf "Secreted effectorP headers:\t"
    cat ${Strain}Headers | wc -l
    Gff=$(ls gene_pred/codingquary/$Organism/$Strain/*/final_genes_appended.gff3)
    EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl ${Strain}Headers $Gff effectorP ID > $EffectorP_Gff
  done
```
```
P.cactorum - 404
EffectorP headers:	17770
Secreted effectorP headers:	850
P.cactorum - 414
EffectorP headers:	18440
Secreted effectorP headers:	969
P.cactorum - 415
EffectorP headers:	20870
Secreted effectorP headers:	945
P.cactorum - 416
EffectorP headers:	21304
Secreted effectorP headers:	927
P.cactorum - 62471
EffectorP headers:	17022
Secreted effectorP headers:	836
P.idaei - 371
EffectorP headers:	17156
Secreted effectorP headers:	799
P.idaei - SCRP370
EffectorP headers:	16918
Secreted effectorP headers:	797
```

## D) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/CAZY/$Organism/$Strain
    mkdir -p $OutDir
    Prefix="$Strain"_CAZY
    CazyHmm=../../dbCAN/dbCAN-fam-HMMs.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
    qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
  done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $File)
echo "$Organism - $Strain"
ProgDir=/home/groups/harrisonlab/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
printf "number of CAZY genes identified:\t"
cat $CazyHeaders | wc -l
Gff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
CazyGff=$OutDir/"$Strain"_CAZY.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

SecretedProts=$(ls gene_pred/combined_sigP/$Organism/$Strain/*_all_secreted.fa)
SecretedHeaders=$(echo $SecretedProts | sed 's/.fa/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
printf "number of Secreted CAZY genes identified:\t"
cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
done
```
```
P.cactorum - 404
number of CAZY genes identified:	638
number of Secreted CAZY genes identified:	328
P.cactorum - 414
number of CAZY genes identified:	1049
number of Secreted CAZY genes identified:	456
P.cactorum - 415
number of CAZY genes identified:	764
number of Secreted CAZY genes identified:	353
P.cactorum - 416
number of CAZY genes identified:	758
number of Secreted CAZY genes identified:	348
P.cactorum - 62471
number of CAZY genes identified:	632
number of Secreted CAZY genes identified:	329
P.idaei - 371
number of CAZY genes identified:	593
number of Secreted CAZY genes identified:	284
P.idaei - SCRP370
number of CAZY genes identified:	602
number of Secreted CAZY genes identified:	288
```

Note - the CAZY genes identified may need further filtering based on e value and
cuttoff length - see below:

Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage

* For fungi, use E-value < 1e-17 and coverage > 0.45

* The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)
 -->

# D) Prediction of RxLRs from - Augustus/CodingQuary gene models

 The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

 The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
Proteome=$(ls gene_pred/codingquary/$Organism/$Strain/*/final_genes_combined.pep.fasta)
Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
mkdir -p $OutDir;
printf "\nstrain: $Strain\tspecies: $Organism\n";
printf "the total number of SigP gene is:\t";
cat $Secretome | grep '>' | wc -l;
printf "the number of unique SigP gene is:\t";
cat $Secretome | grep '>' | cut -f1 | tr -d ' '| sort | uniq | wc -l;

printf "the number of SigP-RxLR genes are:\t";
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_all_secreted_RxLR_regex.fa;
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $OutDir/"$Strain"_RxLR_regex.txt
cat $OutDir/"$Strain"_RxLR_regex.txt | wc -l

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa


printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $OutDir/"$Strain"_RxLR_EER_regex.txt
cat $OutDir/"$Strain"_RxLR_EER_regex.txt | wc -l
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa

printf "\n"

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_regex.txt
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_EER_regex.txt

cat $Gff | grep -w -f $OutDir/"$Strain"_RxLR_regex.txt> $OutDir/"$Strain"_RxLR_regex.gff3
cat $Gff | grep -w -f $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.gff3
done
```

```
strain: 11-40	species: P.cactorum
the total number of SigP gene is:	16631
the number of unique SigP gene is:	3273
the number of SigP-RxLR genes are:	279
the number of SigP-RxLR-EER genes are:	114


strain: 12420	species: P.cactorum
the total number of SigP gene is:	16145
the number of unique SigP gene is:	3139
the number of SigP-RxLR genes are:	273
the number of SigP-RxLR-EER genes are:	114


strain: 15_13	species: P.cactorum
the total number of SigP gene is:	16187
the number of unique SigP gene is:	3142
the number of SigP-RxLR genes are:	269
the number of SigP-RxLR-EER genes are:	112


strain: 15_7	species: P.cactorum
the total number of SigP gene is:	16028
the number of unique SigP gene is:	3110
the number of SigP-RxLR genes are:	266
the number of SigP-RxLR-EER genes are:	105


strain: 17-21	species: P.cactorum
the total number of SigP gene is:	17022
the number of unique SigP gene is:	3320
the number of SigP-RxLR genes are:	295
the number of SigP-RxLR-EER genes are:	131


strain: 2003_3	species: P.cactorum
the total number of SigP gene is:	16784
the number of unique SigP gene is:	3289
the number of SigP-RxLR genes are:	289
the number of SigP-RxLR-EER genes are:	119


strain: 4032	species: P.cactorum
the total number of SigP gene is:	16284
the number of unique SigP gene is:	3172
the number of SigP-RxLR genes are:	276
the number of SigP-RxLR-EER genes are:	112


strain: 4040	species: P.cactorum
the total number of SigP gene is:	16816
the number of unique SigP gene is:	3310
the number of SigP-RxLR genes are:	284
the number of SigP-RxLR-EER genes are:	119


strain: 404	species: P.cactorum
the total number of SigP gene is:	19599
the number of unique SigP gene is:	3924
the number of SigP-RxLR genes are:	296
the number of SigP-RxLR-EER genes are:	118


strain: 415	species: P.cactorum
the total number of SigP gene is:	17325
the number of unique SigP gene is:	3397
the number of SigP-RxLR genes are:	299
the number of SigP-RxLR-EER genes are:	120


strain: 416	species: P.cactorum
the total number of SigP gene is:	17468
the number of unique SigP gene is:	3469
the number of SigP-RxLR genes are:	306
the number of SigP-RxLR-EER genes are:	122


strain: 62471	species: P.cactorum
the total number of SigP gene is:	17048
the number of unique SigP gene is:	3354
the number of SigP-RxLR genes are:	292
the number of SigP-RxLR-EER genes are:	133


strain: P295	species: P.cactorum
the total number of SigP gene is:	16146
the number of unique SigP gene is:	3117
the number of SigP-RxLR genes are:	270
the number of SigP-RxLR-EER genes are:	122


strain: P421_v2	species: P.cactorum
the total number of SigP gene is:	7065
the number of unique SigP gene is:	2972
the number of SigP-RxLR genes are:	267
the number of SigP-RxLR-EER genes are:	115


strain: PC13_15	species: P.cactorum
the total number of SigP gene is:	16328
the number of unique SigP gene is:	3168
the number of SigP-RxLR genes are:	281
the number of SigP-RxLR-EER genes are:	121


strain: R36_14	species: P.cactorum
the total number of SigP gene is:	17028
the number of unique SigP gene is:	3379
the number of SigP-RxLR genes are:	288
the number of SigP-RxLR-EER genes are:	130


strain: 371	species: P.idaei
the total number of SigP gene is:	15536
the number of unique SigP gene is:	3049
the number of SigP-RxLR genes are:	255
the number of SigP-RxLR-EER genes are:	121


strain: SCRP370	species: P.idaei
the total number of SigP gene is:	15635
the number of unique SigP gene is:	3104
the number of SigP-RxLR genes are:	258
the number of SigP-RxLR-EER genes are:	116


strain: SCRP376	species: P.idaei
the total number of SigP gene is:	14942
the number of unique SigP gene is:	2936
the number of SigP-RxLR genes are:	248
the number of SigP-RxLR-EER genes are:	112
```


### G) From Secreted gene models - Hmm evidence of RxLR effectors

```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
    HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
    mkdir -p $OutDir
    HmmResults="$Strain"_RxLR_hmmer.txt
    hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
    echo "$Organism $Strain"
    cat $OutDir/$HmmResults | grep 'Initial search space'
    cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
    HmmFasta="$Strain"_RxLR_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
    Headers="$Strain"_RxLR_hmmer_headers.txt
    cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' | sort | uniq > $OutDir/$Headers
    Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
    cat $Gff | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
  done
```

```
  P.cactorum 404
  Initial search space (Z):              27775  [actual number of targets]
  Domain search space  (domZ):             127  [number of targets reported over threshold]
  P.cactorum 414
  Initial search space (Z):              32832  [actual number of targets]
  Domain search space  (domZ):             146  [number of targets reported over threshold]
  P.cactorum 415
  Initial search space (Z):              31944  [actual number of targets]
  Domain search space  (domZ):             132  [number of targets reported over threshold]
  P.cactorum 416
  Initial search space (Z):              32780  [actual number of targets]
  Domain search space  (domZ):             129  [number of targets reported over threshold]
  P.cactorum 62471
  Initial search space (Z):              27191  [actual number of targets]
  Domain search space  (domZ):             142  [number of targets reported over threshold]
  P.idaei 371
  Initial search space (Z):              27253  [actual number of targets]
  Domain search space  (domZ):             107  [number of targets reported over threshold]
  P.idaei SCRP370
  Initial search space (Z):              26983  [actual number of targets]
  Domain search space  (domZ):             105  [number of targets reported over threshold]
```


### F) Combining RxLRs from Regex and hmm searches


The total RxLRs are

```bash
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_RxLR_EER_regex.txt | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
Proteome=$(ls gene_pred/codingquary/$Organism/$Strain/*/final_genes_combined.pep.fasta)
HmmRxLR=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer_headers.txt
echo "$Organism - $Strain"
echo "Number of RxLRs identified by Regex:"
cat $RegexRxLR | sort | uniq | wc -l
echo "Number of RxLRs identified by Hmm:"
cat $HmmRxLR | sort | uniq | wc -l
echo "Number of RxLRs in combined dataset:"
cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
# echo "Number of RxLRs in both datasets:"
# cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
echo ""
# echo "Extracting RxLRs from datasets"
OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
mkdir -p $OutDir
cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_RxLR_headers.txt
Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
cat $Gff | grep -w -f $OutDir/"$Strain"_total_RxLR_headers.txt > $OutDir/"$Strain"_total_RxLR.gff
echo "Number of genes in the extracted gff file:"
cat $OutDir/"$Strain"_total_RxLR.gff | grep -w 'gene' | wc -l
done
```

```
P.cactorum - 11-40
Number of RxLRs identified by Regex:
114
Number of RxLRs identified by Hmm:
130
Number of RxLRs in combined dataset:
154
Number of genes in the extracted gff file:
154

P.cactorum - 12420
Number of RxLRs identified by Regex:
113
Number of RxLRs identified by Hmm:
124
Number of RxLRs in combined dataset:
145
Number of genes in the extracted gff file:
145

P.cactorum - 15_13
Number of RxLRs identified by Regex:
111
Number of RxLRs identified by Hmm:
119
Number of RxLRs in combined dataset:
144
Number of genes in the extracted gff file:
144

P.cactorum - 15_7
Number of RxLRs identified by Regex:
105
Number of RxLRs identified by Hmm:
115
Number of RxLRs in combined dataset:
136
Number of genes in the extracted gff file:
136

P.cactorum - 17-21
Number of RxLRs identified by Regex:
131
Number of RxLRs identified by Hmm:
134
Number of RxLRs in combined dataset:
165
Number of genes in the extracted gff file:
165

P.cactorum - 2003_3
Number of RxLRs identified by Regex:
119
Number of RxLRs identified by Hmm:
128
Number of RxLRs in combined dataset:
152
Number of genes in the extracted gff file:
152

P.cactorum - 4032
Number of RxLRs identified by Regex:
111
Number of RxLRs identified by Hmm:
120
Number of RxLRs in combined dataset:
143
Number of genes in the extracted gff file:
143

P.cactorum - 4040
Number of RxLRs identified by Regex:
118
Number of RxLRs identified by Hmm:
126
Number of RxLRs in combined dataset:
150
Number of genes in the extracted gff file:
150

P.cactorum - 404
Number of RxLRs identified by Regex:
117
Number of RxLRs identified by Hmm:
130
Number of RxLRs in combined dataset:
154
Number of genes in the extracted gff file:
154

P.cactorum - 415
Number of RxLRs identified by Regex:
120
Number of RxLRs identified by Hmm:
127
Number of RxLRs in combined dataset:
155
Number of genes in the extracted gff file:
152

P.cactorum - 416
Number of RxLRs identified by Regex:
121
Number of RxLRs identified by Hmm:
128
Number of RxLRs in combined dataset:
156
Number of genes in the extracted gff file:
152

P.cactorum - 62471
Number of RxLRs identified by Regex:
133
Number of RxLRs identified by Hmm:
139
Number of RxLRs in combined dataset:
167
Number of genes in the extracted gff file:
167

P.cactorum - P295
Number of RxLRs identified by Regex:
122
Number of RxLRs identified by Hmm:
130
Number of RxLRs in combined dataset:
155
Number of genes in the extracted gff file:
155

P.cactorum - P421_v2
Number of RxLRs identified by Regex:
114
Number of RxLRs identified by Hmm:
128
Number of RxLRs in combined dataset:
151
Number of genes in the extracted gff file:
151

P.cactorum - PC13_15
Number of RxLRs identified by Regex:
120
Number of RxLRs identified by Hmm:
131
Number of RxLRs in combined dataset:
154
Number of genes in the extracted gff file:
154

P.cactorum - R36_14
Number of RxLRs identified by Regex:
130
Number of RxLRs identified by Hmm:
141
Number of RxLRs in combined dataset:
168
Number of genes in the extracted gff file:
168

P.idaei - 371
Number of RxLRs identified by Regex:
120
Number of RxLRs identified by Hmm:
100
Number of RxLRs in combined dataset:
143
Number of genes in the extracted gff file:
133

P.idaei - SCRP370
Number of RxLRs identified by Regex:
115
Number of RxLRs identified by Hmm:
108
Number of RxLRs in combined dataset:
144
Number of genes in the extracted gff file:
144

P.idaei - SCRP376
Number of RxLRs identified by Regex:
112
Number of RxLRs identified by Hmm:
106
Number of RxLRs in combined dataset:
137
Number of genes in the extracted gff file:
137
```



```bash
for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_Aug_WY_hmmer.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_Aug_WY_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_Aug_WY_hmmer_headers.txt
# cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' | sort | uniq > $OutDir/$Headers
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | cut -f1 | sort | uniq > $OutDir/$Headers
# Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
# cat $Gff | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_WY_hmmer.gff
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/$Headers $Gff $HmmModel Name Augustus > $OutDir/"$Strain"_Aug_WY_hmmer.gff
done
```

<!-- 72 of 137 RxLRs were in the predicted 797 effectorP genes.
```bash
  cat analysis/RxLR_effectors/combined_evidence/P.idaei/SCRP370/SCRP370_total_RxLR_headers.txt ${Strain}Headers | sort | cut -f1 -d '.' | uniq -d | wc -l
``` -->

### D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:


```bash
  HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
  LFLAK_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm)
  DWL_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm)
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
    mkdir -p $OutDir
    echo "$Organism - $Strain"
    # Run hmm searches LFLAK domains
    CrinklerProts_LFLAK=$OutDir/"$Strain"_pub_CRN_LFLAK_hmm.txt
    hmmsearch -T0 $LFLAK_hmm $Proteome > $CrinklerProts_LFLAK
    cat $CrinklerProts_LFLAK | grep 'Initial search space'
    cat $CrinklerProts_LFLAK | grep 'number of targets reported over threshold'
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
    $ProgDir/hmmer2fasta.pl $CrinklerProts_LFLAK $Proteome > $OutDir/"$Strain"_pub_CRN_LFLAK_hmm.fa
    # Run hmm searches DWL domains
    CrinklerProts_DWL=$OutDir/"$Strain"_pub_CRN_DWL_hmm.txt
    hmmsearch -T0 $DWL_hmm $Proteome > $CrinklerProts_DWL
    cat $CrinklerProts_DWL | grep 'Initial search space'
    cat $CrinklerProts_DWL | grep 'number of targets reported over threshold'
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
    $ProgDir/hmmer2fasta.pl $CrinklerProts_DWL $Proteome > $OutDir/"$Strain"_pub_CRN_DWL_hmm.fa
    # Identify the genes detected in both models
    cat $OutDir/"$Strain"_pub_CRN_LFLAK_hmm.fa $OutDir/"$Strain"_pub_CRN_DWL_hmm.fa | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $OutDir/"$Strain"_pub_CRN_LFLAK_DWL.txt
    cat $OutDir/"$Strain"_pub_CRN_LFLAK_DWL.txt | wc -l
  done
```

```
P.cactorum - 11-40
Initial search space (Z):              28827  [actual number of targets]
Domain search space  (domZ):             109  [number of targets reported over threshold]
Initial search space (Z):              28827  [actual number of targets]
Domain search space  (domZ):              88  [number of targets reported over threshold]
59
P.cactorum - 12420
Initial search space (Z):              26231  [actual number of targets]
Domain search space  (domZ):             108  [number of targets reported over threshold]
Initial search space (Z):              26231  [actual number of targets]
Domain search space  (domZ):              86  [number of targets reported over threshold]
66
P.cactorum - 15_13
Initial search space (Z):              26744  [actual number of targets]
Domain search space  (domZ):             114  [number of targets reported over threshold]
Initial search space (Z):              26744  [actual number of targets]
Domain search space  (domZ):              99  [number of targets reported over threshold]
72
P.cactorum - 15_7
Initial search space (Z):              26604  [actual number of targets]
Domain search space  (domZ):             110  [number of targets reported over threshold]
Initial search space (Z):              26604  [actual number of targets]
Domain search space  (domZ):              88  [number of targets reported over threshold]
71
P.cactorum - 17-21
Initial search space (Z):              28577  [actual number of targets]
Domain search space  (domZ):             103  [number of targets reported over threshold]
Initial search space (Z):              28577  [actual number of targets]
Domain search space  (domZ):              90  [number of targets reported over threshold]
59
P.cactorum - 2003_3
Initial search space (Z):              29786  [actual number of targets]
Domain search space  (domZ):             111  [number of targets reported over threshold]
Initial search space (Z):              29786  [actual number of targets]
Domain search space  (domZ):              97  [number of targets reported over threshold]
61
P.cactorum - 4032
Initial search space (Z):              26819  [actual number of targets]
Domain search space  (domZ):             121  [number of targets reported over threshold]
Initial search space (Z):              26819  [actual number of targets]
Domain search space  (domZ):              96  [number of targets reported over threshold]
77
P.cactorum - 4040
Initial search space (Z):              29530  [actual number of targets]
Domain search space  (domZ):             122  [number of targets reported over threshold]
Initial search space (Z):              29530  [actual number of targets]
Domain search space  (domZ):              93  [number of targets reported over threshold]
73
P.cactorum - 404
Initial search space (Z):              35766  [actual number of targets]
Domain search space  (domZ):              96  [number of targets reported over threshold]
Initial search space (Z):              35766  [actual number of targets]
Domain search space  (domZ):              84  [number of targets reported over threshold]
51
P.cactorum - 415
Initial search space (Z):              26796  [actual number of targets]
Domain search space  (domZ):             123  [number of targets reported over threshold]
Initial search space (Z):              26796  [actual number of targets]
Domain search space  (domZ):              94  [number of targets reported over threshold]
75
P.cactorum - 416
Initial search space (Z):              29391  [actual number of targets]
Domain search space  (domZ):             130  [number of targets reported over threshold]
Initial search space (Z):              29391  [actual number of targets]
Domain search space  (domZ):             107  [number of targets reported over threshold]
78
P.cactorum - 62471
Initial search space (Z):              29596  [actual number of targets]
Domain search space  (domZ):             109  [number of targets reported over threshold]
Initial search space (Z):              29596  [actual number of targets]
Domain search space  (domZ):              90  [number of targets reported over threshold]
61
P.cactorum - P295
Initial search space (Z):              26327  [actual number of targets]
Domain search space  (domZ):             111  [number of targets reported over threshold]
Initial search space (Z):              26327  [actual number of targets]
Domain search space  (domZ):              99  [number of targets reported over threshold]
73
P.cactorum - P421_v2
Initial search space (Z):              26124  [actual number of targets]
Domain search space  (domZ):              97  [number of targets reported over threshold]
Initial search space (Z):              26124  [actual number of targets]
Domain search space  (domZ):              79  [number of targets reported over threshold]
58
P.cactorum - PC13_15
Initial search space (Z):              26607  [actual number of targets]
Domain search space  (domZ):             111  [number of targets reported over threshold]
Initial search space (Z):              26607  [actual number of targets]
Domain search space  (domZ):              99  [number of targets reported over threshold]
74
P.cactorum - R36_14
Initial search space (Z):              29929  [actual number of targets]
Domain search space  (domZ):             126  [number of targets reported over threshold]
Initial search space (Z):              29929  [actual number of targets]
Domain search space  (domZ):             108  [number of targets reported over threshold]
81
P.idaei - 371
Initial search space (Z):              25387  [actual number of targets]
Domain search space  (domZ):              75  [number of targets reported over threshold]
Initial search space (Z):              25387  [actual number of targets]
Domain search space  (domZ):              66  [number of targets reported over threshold]
41
P.idaei - SCRP370
Initial search space (Z):              27882  [actual number of targets]
Domain search space  (domZ):              61  [number of targets reported over threshold]
Initial search space (Z):              27882  [actual number of targets]
Domain search space  (domZ):              47  [number of targets reported over threshold]
29
P.idaei - SCRP376
Initial search space (Z):              25490  [actual number of targets]
Domain search space  (domZ):              57  [number of targets reported over threshold]
Initial search space (Z):              25490  [actual number of targets]
Domain search space  (domZ):              48  [number of targets reported over threshold]
26
```
<!--
23 of 92 P.idaei CRNs were in the effectorP dataset.
```bash
cat analysis/CRN_effectors/hmmer_CRN/P.idaei/SCRP370/SCRP370_pub_CRN_LFLAK_DWL.txt ${Strain}Headers | sort | cut -f1 -d '.' | uniq -d | wc -l
``` -->

Extract gff annotations for Crinklers:

```bash
  for CRNlist in $(ls analysis/CRN_effectors/hmmer_CRN/*/*/*_pub_CRN_LFLAK_DWL.txt | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    Strain=$(echo $CRNlist | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $CRNlist | rev | cut -f3 -d '/' | rev)
    OutName=$(echo $CRNlist | sed 's/.txt/.gff/g')
    echo "$Organism - $Strain"
    Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
    cat $CRNlist | sed -r 's/\.t.$//g' > tmp.txt
    cat $Gff | grep -w -f tmp.txt > $OutName
    rm tmp.txt
  done
```



### E) From ORF gene models - Signal peptide & RxLR motif

Required programs:
 * SigP
 * Phobius
 * biopython


#### E.1) Prediction using SignalP
Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | grep -e '415' -e '416'); do
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
SplitDir=gene_pred/ORF_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_ORF_preds
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_ORF_preds_*); do
Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
while [ $Jobs -gt 6 ]; do
sleep 1
printf "."
Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
done		
printf "\n"
echo $File
qsub $ProgDir/pred_sigP.sh $File
# qsub $ProgDir/pred_sigP.sh $File signalp-3.0
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```

 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

 ```bash
 for SplitDir in $(ls -d gene_pred/ORF_split/*/* | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | grep -e '415' -e '416'); do
 Strain=$(echo $SplitDir | cut -d '/' -f4)
 Organism=$(echo $SplitDir | cut -d '/' -f3)
 echo "$Organism - $Strain"
 for SigpDir in $(ls -d gene_pred/ORF_sig* | grep -v '_signalp-3.0' | cut -f2 -d'/'); do
 InStringAA=''
 InStringNeg=''
 InStringTab=''
 InStringTxt=''
 for GRP in $(ls -l $SplitDir/*_ORF_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
 InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.aa";
 InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp_neg.aa";
 InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.tab";
 InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.txt";
 done
 cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_ORF_sp.aa
 cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_ORF_neg_sp.aa
 tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_ORF_sp.tab
 cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_ORF_sp.txt
 done
 done
 ```
<!--
 ```bash
  qlogin
  WorkDir=/tmp/$USER/sigp/ORF
  mkdir -p $WorkDir

  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | head -n1); do
  SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)

  cp $Proteome $WorkDir/proteins.fa
  cd $WorkDir

  echo "Running using SignalP version: $SigP_Version"
  echo "Installing SignalP"
  Tarball=$(ls /home/armita/prog/signalP/signalp-3.0/signalp-3.0.Linux.tar.Z)
  cp $Tarball .
  tar -zxf *.tar.Z
  sed -i "s&SIGNALP=.*&SIGNALP=$PWD/signalp-3.0&g" signalp-3.0/signalp
  sed -i 's&AWK=nawk&AWK=/usr/bin/awk&g' signalp-3.0/signalp
  sed -i 's&AWK=/usr/bin/gawk&AWK=/usr/bin/awk&g' signalp-3.0/signalp
  sed -i 's&MAXSEQ=4000&MAXSEQ=1000000&g' signalp-3.0/signalp
  echo "SignalP installed"
  ls -lh signalp-3.0/syn/
  echo "Running test"
  signalp-3.0/signalp -t euk signalp-3.0/test/test.seq > signalp-3.0/tmp.txt
  echo "test run"
  signalp-3.0/signalp -t euk -f short -trunc 70 "proteins.fa" > "${Strain}"_sp.txt
  rm -rf signalp-3.0
  # $SigP_Version -t euk -f short -trunc 70 "proteins.fa" > "${Strain}"_sp.txt
  echo '----------------------------------------------------------------------' >> "${Strain}"_sp.txt
  PathToAnnotateSigP=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  $PathToAnnotateSigP/sigP_3.0_parser.py --inp_sigP "${Strain}"_sp.txt --out_tab "${Strain}"_sp.tab --out_fasta "${Strain}"_sp.aa --out_neg "${Strain}"_sp_neg.aa --inp_fasta "proteins.fa"
  OutDir=$ProjDir/gene_pred/ORF_signalp-3.0/$Organism/$Strain/.

  rm proteins.fa
  mkdir -p $OutDir
  cp * $OutDir/.
  cd $ProjDir
  done

 ``` -->


#### E.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
	for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/phobius/$Organism/$Strain
    mkdir -p $OutDir
  	phobius.pl $Proteome > $OutDir/"$Strain"_phobius_ORF.txt
  	cat $OutDir/"$Strain"_phobius_ORF.txt | grep -B1 'SIGNAL' | grep 'ID' | sed s'/ID.*g/g/g' > $OutDir/"$Strain"_phobius_headers_ORF.txt
  done
```

Secreted proteins from different sources were combined into a single file:

```bash
 for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
   Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   echo "$Organism - $Strain"
   OutDir=gene_pred/combined_sigP_ORF/$Organism/$Strain
   mkdir -p $OutDir
   echo "The following number of sequences were predicted as secreted:"
   cat gene_pred/ORF_sig*/$Organism/$Strain/*_ORF_sp.aa > $OutDir/"$Strain"_all_secreted.fa
   cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | tr -d '>' | tr -d ' ' | sed "s/HMM_score\t/HMM_score=\t/g" > $OutDir/"$Strain"_all_secreted_headers.txt
   cat $OutDir/"$Strain"_all_secreted_headers.txt | wc -l
   echo "This represented the following number of unique genes:"
   cat gene_pred/ORF_sig*/$Organism/$Strain/*_ORF_sp.aa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
   cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
 done
```


```
P.cactorum - 11-40
The following number of sequences were predicted as secreted:
17518
This represented the following number of unique genes:
17491
P.cactorum - 12420
The following number of sequences were predicted as secreted:
17579
This represented the following number of unique genes:
17525
P.cactorum - 15_13
The following number of sequences were predicted as secreted:
17687
This represented the following number of unique genes:
17660
P.cactorum - 15_7
The following number of sequences were predicted as secreted:
21929
This represented the following number of unique genes:
19425
P.cactorum - 17-21
The following number of sequences were predicted as secreted:
17481
This represented the following number of unique genes:
17469
P.cactorum - 2003_3
The following number of sequences were predicted as secreted:
17630
This represented the following number of unique genes:
17630
P.cactorum - 4032
The following number of sequences were predicted as secreted:
22067
This represented the following number of unique genes:
19641
P.cactorum - 4040
The following number of sequences were predicted as secreted:
21875
This represented the following number of unique genes:
19427
P.cactorum - 404
The following number of sequences were predicted as secreted:
43020
This represented the following number of unique genes:
29642
P.cactorum - 415
The following number of sequences were predicted as secreted:
41447
This represented the following number of unique genes:
27389
P.cactorum - 416
The following number of sequences were predicted as secreted:
42503
This represented the following number of unique genes:
27683
P.cactorum - 62471
The following number of sequences were predicted as secreted:
17830
This represented the following number of unique genes:
17816
P.cactorum - P295
The following number of sequences were predicted as secreted:
21889
This represented the following number of unique genes:
19461
P.cactorum - P421_v2
The following number of sequences were predicted as secreted:
17434
This represented the following number of unique genes:
17434
P.cactorum - PC13_15
The following number of sequences were predicted as secreted:
17837
This represented the following number of unique genes:
17837
P.cactorum - R36_14
The following number of sequences were predicted as secreted:
17959
This represented the following number of unique genes:
17959
P.idaei - 371
The following number of sequences were predicted as secreted:
37207
This represented the following number of unique genes:
25622
P.idaei - SCRP370
The following number of sequences were predicted as secreted:
17254
This represented the following number of unique genes:
17254
P.idaei - SCRP376
The following number of sequences were predicted as secreted:
17208
This represented the following number of unique genes:
17208
```

#### E.3) Prediction of RxLRs


Names of ORFs containing signal peptides were extracted from fasta files. This
included information on the position and hmm score of RxLRs.

```bash
	FastaFile=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp.aa
	SigP_headers=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_names.txt
	cat $FastaFile | grep '>' | sed -r 's/>//g' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | sed 's/--//g' > $SigP_headers
```

Due to the nature of predicting ORFs, some features overlapped with one another.
A single ORF was selected from each set of overlapped ORFs. This was was
selected on the basis of its SignalP Hmm score. Biopython was used to identify
overlaps and identify the ORF with the best signalP score.

Due to the nature of predicting ORFs, some features overlapped with one another.
A single ORF was selected from each set of overlapped ORFs. This was was
selected on the basis of its SignalP Hmm score. Biopython was used to identify
overlaps and identify the ORF with the best signalP score.

```bash
 for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff3 | grep -v -e 'atg' | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
   Organism=$(echo $ORF_Gff | rev |  cut -d '/' -f3 | rev) ;
   Strain=$(echo $ORF_Gff | rev | cut -d '/' -f2 | rev);
   OutDir=$(ls -d gene_pred/combined_sigP_ORF/$Organism/$Strain)
   echo "$Organism - $Strain"
   # SigP_fasta=$(ls $OutDir/"$Strain"_all_secreted.fa)
   SigP_headers=$(ls $OutDir/"$Strain"_all_secreted_headers.txt)
   ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)

   SigP_Gff=$OutDir/"$Strain"_all_secreted_unmerged.gff
   SigP_Merged_Gff=$OutDir/"$Strain"_all_secreted_merged.gff
   SigP_Merged_txt=$OutDir/"$Strain"_all_secreted_merged.txt
   SigP_Merged_AA=$OutDir/"$Strain"_all_secreted_merged.aa

   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/extract_gff_for_sigP_hits.pl $SigP_headers $ORF_Gff SigP Name > $SigP_Gff
   ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
   $ProgDir/make_gff_database.py --inp $SigP_Gff --db sigP_ORF.db
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF.db --id sigP_ORF --out sigP_ORF_merged.db --gff > $SigP_Merged_Gff
   cat $SigP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d'=' | rev > $SigP_Merged_txt
   # $ProgDir/extract_from_fasta.py --fasta $SigP_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
   $ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
 done
```

The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep -v -e 'atg' | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
mkdir -p $OutDir;
printf "\nstrain: $Strain\tspecies: $Organism\n";
printf "the number of SigP gene is:\t";
cat $Secretome | grep '>' | cut -f1 | tr -d ' ' | sort | uniq | wc -l
printf "the number of SigP-RxLR genes are:\t";
$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa;
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt
cat $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l
printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' '> $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt  $SigP_Merged_Gff 	RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
SigP_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt $SigP_Gff	RxLR_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt  $SigP_Gff	RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.gff
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt $SigP_Merged_Gff RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex.gff
RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_regex_merged.gff
RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_regex_merged.txt
RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_regex_merged.aa
RxLR_EER_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.gff
RxLR_EER_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.txt
RxLR_EER_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff --db sigP_ORF_RxLR.db
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.gff --db sigP_ORF_RxLR_EER.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR.db --id sigP_ORF_RxLR --out sigP_ORF_RxLR_merged.db --gff > $RxLR_Merged_Gff
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR_EER.db --id sigP_ORF_RxLR_EER --out sigP_ORF_RxLR_EER_merged.db --gff > $RxLR_EER_Merged_Gff
cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
cat $RxLR_EER_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_EER_Merged_txt
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_EER_Merged_txt > $RxLR_EER_Merged_AA
printf "Merged RxLR regex proteins:\t"
cat $RxLR_Merged_AA | grep '>' | wc -l
printf "Merged RxLR-EER regex proteins:\t"
cat $RxLR_EER_Merged_AA | grep '>' | wc -l
printf "\n"
done
```

```
strain: 11-40	species: P.cactorum
the number of SigP gene is:	17491
the number of SigP-RxLR genes are:	1007
the number of SigP-RxLR-EER genes are:	178
Merged RxLR regex proteins:	869
Merged RxLR-EER regex proteins:	160


strain: 12420	species: P.cactorum
the number of SigP gene is:	17525
the number of SigP-RxLR genes are:	1034
the number of SigP-RxLR-EER genes are:	175
Merged RxLR regex proteins:	889
Merged RxLR-EER regex proteins:	158


strain: 15_13	species: P.cactorum
the number of SigP gene is:	17660
the number of SigP-RxLR genes are:	1055
the number of SigP-RxLR-EER genes are:	177
Merged RxLR regex proteins:	912
Merged RxLR-EER regex proteins:	159


strain: 15_7	species: P.cactorum
the number of SigP gene is:	19425
the number of SigP-RxLR genes are:	1151
the number of SigP-RxLR-EER genes are:	175
Merged RxLR regex proteins:	987
Merged RxLR-EER regex proteins:	158


strain: 17-21	species: P.cactorum
the number of SigP gene is:	17469
the number of SigP-RxLR genes are:	1030
the number of SigP-RxLR-EER genes are:	176
Merged RxLR regex proteins:	890
Merged RxLR-EER regex proteins:	162


strain: 2003_3	species: P.cactorum
the number of SigP gene is:	17630
the number of SigP-RxLR genes are:	1048
the number of SigP-RxLR-EER genes are:	174
Merged RxLR regex proteins:	908
Merged RxLR-EER regex proteins:	157


strain: 4032	species: P.cactorum
the number of SigP gene is:	19641
the number of SigP-RxLR genes are:	1152
the number of SigP-RxLR-EER genes are:	182
Merged RxLR regex proteins:	989
Merged RxLR-EER regex proteins:	163


strain: 4040	species: P.cactorum
the number of SigP gene is:	19427
the number of SigP-RxLR genes are:	1136
the number of SigP-RxLR-EER genes are:	177
Merged RxLR regex proteins:	969
Merged RxLR-EER regex proteins:	158


strain: 404	species: P.cactorum
the number of SigP gene is:	29642
the number of SigP-RxLR genes are:	1633
the number of SigP-RxLR-EER genes are:	192
Merged RxLR regex proteins:	1371
Merged RxLR-EER regex proteins:	168


strain: 415	species: P.cactorum
the number of SigP gene is:	28547
the number of SigP-RxLR genes are:	1692
the number of SigP-RxLR-EER genes are:	214
Merged RxLR regex proteins:	1372
Merged RxLR-EER regex proteins:	184


strain: 416	species: P.cactorum
the number of SigP gene is:	29221
the number of SigP-RxLR genes are:	1700
the number of SigP-RxLR-EER genes are:	211
Merged RxLR regex proteins:	1343
Merged RxLR-EER regex proteins:	179


strain: 62471	species: P.cactorum
the number of SigP gene is:	17816
the number of SigP-RxLR genes are:	1069
the number of SigP-RxLR-EER genes are:	192
Merged RxLR regex proteins:	927
Merged RxLR-EER regex proteins:	175


strain: P295	species: P.cactorum
the number of SigP gene is:	19461
the number of SigP-RxLR genes are:	1144
the number of SigP-RxLR-EER genes are:	193
Merged RxLR regex proteins:	982
Merged RxLR-EER regex proteins:	174


strain: P421_v2	species: P.cactorum
the number of SigP gene is:	17434
the number of SigP-RxLR genes are:	1014
the number of SigP-RxLR-EER genes are:	180
Merged RxLR regex proteins:	872
Merged RxLR-EER regex proteins:	163


strain: PC13_15	species: P.cactorum
the number of SigP gene is:	17837
the number of SigP-RxLR genes are:	1056
the number of SigP-RxLR-EER genes are:	180
Merged RxLR regex proteins:	916
Merged RxLR-EER regex proteins:	162


strain: R36_14	species: P.cactorum
the number of SigP gene is:	17959
the number of SigP-RxLR genes are:	1079
the number of SigP-RxLR-EER genes are:	185
Merged RxLR regex proteins:	934
Merged RxLR-EER regex proteins:	169


strain: 371	species: P.idaei
the number of SigP gene is:	25622
the number of SigP-RxLR genes are:	1529
the number of SigP-RxLR-EER genes are:	222
Merged RxLR regex proteins:	1254
Merged RxLR-EER regex proteins:	189


strain: SCRP370	species: P.idaei
the number of SigP gene is:	17254
the number of SigP-RxLR genes are:	1044
the number of SigP-RxLR-EER genes are:	194
Merged RxLR regex proteins:	887
Merged RxLR-EER regex proteins:	172


strain: SCRP376	species: P.idaei
the number of SigP gene is:	17208
the number of SigP-RxLR genes are:	1039
the number of SigP-RxLR-EER genes are:	193
Merged RxLR regex proteins:	880
Merged RxLR-EER regex proteins:	171
```

### E5) From ORF gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search ORFs predicted with atg.pl. These were run with the following commands:


```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep -v -e 'atg' | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_ORF_WY_hmmer.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_ORF_WY_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_ORF_WY_hmmer_headers.txt
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
SigP_Merged_Gff=$(ls gene_pred/combined_sigP_ORF/$Organism/$Strain/${Strain}_all_secreted_merged.gff)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_WY_hmmer.gff
done
```

```
P.cactorum 11-40
Initial search space (Z):              17518  [actual number of targets]
Domain search space  (domZ):             101  [number of targets reported over threshold]
P.cactorum 12420
Initial search space (Z):              17579  [actual number of targets]
Domain search space  (domZ):             101  [number of targets reported over threshold]
P.cactorum 15_13
Initial search space (Z):              17687  [actual number of targets]
Domain search space  (domZ):             100  [number of targets reported over threshold]
P.cactorum 15_7
Initial search space (Z):              21929  [actual number of targets]
Domain search space  (domZ):             112  [number of targets reported over threshold]
P.cactorum 17-21
Initial search space (Z):              17481  [actual number of targets]
Domain search space  (domZ):             100  [number of targets reported over threshold]
P.cactorum 2003_3
Initial search space (Z):              17630  [actual number of targets]
Domain search space  (domZ):             102  [number of targets reported over threshold]
P.cactorum 4032
Initial search space (Z):              22067  [actual number of targets]
Domain search space  (domZ):             112  [number of targets reported over threshold]
P.cactorum 4040
Initial search space (Z):              21875  [actual number of targets]
Domain search space  (domZ):             109  [number of targets reported over threshold]
P.cactorum 404
Initial search space (Z):              43020  [actual number of targets]
Domain search space  (domZ):             208  [number of targets reported over threshold]
P.cactorum 415
Initial search space (Z):              41447  [actual number of targets]
Domain search space  (domZ):             234  [number of targets reported over threshold]
P.cactorum 416
Initial search space (Z):              42503  [actual number of targets]
Domain search space  (domZ):             236  [number of targets reported over threshold]
P.cactorum 62471
Initial search space (Z):              17830  [actual number of targets]
Domain search space  (domZ):             111  [number of targets reported over threshold]
P.cactorum P295
Initial search space (Z):              21889  [actual number of targets]
Domain search space  (domZ):             124  [number of targets reported over threshold]
P.cactorum P421_v2
Initial search space (Z):              17434  [actual number of targets]
Domain search space  (domZ):             103  [number of targets reported over threshold]
P.cactorum PC13_15
Initial search space (Z):              17837  [actual number of targets]
Domain search space  (domZ):             104  [number of targets reported over threshold]
P.cactorum R36_14
Initial search space (Z):              17959  [actual number of targets]
Domain search space  (domZ):             107  [number of targets reported over threshold]
P.idaei 371
Initial search space (Z):              37207  [actual number of targets]
Domain search space  (domZ):             209  [number of targets reported over threshold]
P.idaei SCRP370
Initial search space (Z):              17254  [actual number of targets]
Domain search space  (domZ):             106  [number of targets reported over threshold]
P.idaei SCRP376
Initial search space (Z):              17208  [actual number of targets]
Domain search space  (domZ):             102  [number of targets reported over threshold]
```


### E6) From ORF gene models - Hmm evidence of RxLR effectors

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_ORF_RxLR_hmmer_unmerged.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_ORF_RxLR_hmmer_headers_unmerged.txt
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
SigP_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_RxLR_hmmer_unmerged.gff3
RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.gff
RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.txt
RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_hmmer_unmerged.gff3 --db sigP_ORF_RxLR_hmm.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR_hmm.db --id sigP_ORF_RxLR_hmm --out sigP_ORF_RxLR_hmm_merged.db --gff > $RxLR_Merged_Gff
cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
printf "Merged RxLR-EER Hmm proteins:\t"
cat $RxLR_Merged_AA | grep '>' | wc -l
done
```

```
P.cactorum 11-40
Initial search space (Z):              17518  [actual number of targets]
Domain search space  (domZ):             137  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	124
P.cactorum 12420
Initial search space (Z):              17579  [actual number of targets]
Domain search space  (domZ):             137  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	124
P.cactorum 15_13
Initial search space (Z):              17687  [actual number of targets]
Domain search space  (domZ):             135  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	122
P.cactorum 15_7
Initial search space (Z):              21929  [actual number of targets]
Domain search space  (domZ):             154  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	124
P.cactorum 17-21
Initial search space (Z):              17481  [actual number of targets]
Domain search space  (domZ):             134  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	124
P.cactorum 2003_3
Initial search space (Z):              17630  [actual number of targets]
Domain search space  (domZ):             135  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	122
P.cactorum 4032
Initial search space (Z):              22067  [actual number of targets]
Domain search space  (domZ):             153  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	128
P.cactorum 4040
Initial search space (Z):              21875  [actual number of targets]
Domain search space  (domZ):             153  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	124
P.cactorum 404
Initial search space (Z):              43020  [actual number of targets]
Domain search space  (domZ):             275  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	127
P.cactorum 415
Initial search space (Z):              41447  [actual number of targets]
Domain search space  (domZ):             306  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	141
P.cactorum 416
Initial search space (Z):              42503  [actual number of targets]
Domain search space  (domZ):             310  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	143
P.cactorum 62471
Initial search space (Z):              17830  [actual number of targets]
Domain search space  (domZ):             146  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	134
P.cactorum P295
Initial search space (Z):              21889  [actual number of targets]
Domain search space  (domZ):             167  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	134
P.cactorum P421_v2
Initial search space (Z):              17434  [actual number of targets]
Domain search space  (domZ):             139  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	126
P.cactorum PC13_15
Initial search space (Z):              17837  [actual number of targets]
Domain search space  (domZ):             141  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	127
P.cactorum R36_14
Initial search space (Z):              17959  [actual number of targets]
Domain search space  (domZ):             144  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	133
P.idaei 371
Initial search space (Z):              37207  [actual number of targets]
Domain search space  (domZ):             279  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	132
P.idaei SCRP370
Initial search space (Z):              17254  [actual number of targets]
Domain search space  (domZ):             138  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	128
P.idaei SCRP376
Initial search space (Z):              17208  [actual number of targets]
Domain search space  (domZ):             136  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	126
```

### E7) Combining RxLRs from Regex and hmm searches


The total ORF RxLRs are

```bash
for RegexRxLREER in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_ORF_RxLR_EER_regex_merged.txt | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Organism=$(echo $RegexRxLREER | rev |  cut -d '/' -f3 | rev)
Strain=$(echo $RegexRxLREER | rev | cut -d '/' -f2 | rev)
Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF.gff3)
Proteome=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
# RegexRxLRfpkm=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers_renamed.txt)
HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_ORF_RxLR_hmm_merged.txt)
# echo "$Organism - $Strain"
# echo "Number of RxLR EERs identified by Regex:"
Regex=$(cat $RegexRxLREER | sort | uniq | wc -l)
# echo "Number of RxLRs identified by Regex with fpkm > 5"
# cat $RegexRxLRfpkm | sort | uniq | wc -l
# echo "Number of RxLRs identified by Hmm:"
Hmm=$(cat $HmmRxLR | sort | uniq | wc -l)
# echo "Number of RxLRs in combined dataset:"
# cat $RegexRxLREER $HmmRxLR $RegexRxLRfpkm | sort | uniq | wc -l
Combined=$(cat $RegexRxLREER $HmmRxLR | sort | uniq | wc -l)
# echo ""
# echo "Extracting RxLRs from datasets"
OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
mkdir -p $OutDir
# cat $RegexRxLREER $RegexRxLRfpkm $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
cat $RegexRxLREER $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_total_ORF_RxLR_headers.txt $Gff ORF_RxLR Name Augustus > $OutDir/"$Strain"_total_ORF_RxLR.gff
# echo "Number of genes in the extracted gff file:"
Extracted=$(cat $OutDir/"$Strain"_total_ORF_RxLR.gff | grep -w 'gene' | wc -l)
printf "$Organism\t$Strain\t$Regex\t$Hmm\t$Combined\t$Extracted\n"
done
```

```
P.cactorum	11-40	160	124	176	176
P.cactorum	12420	158	124	174	174
P.cactorum	15_13	159	122	174	174
P.cactorum	15_7	158	124	174	174
P.cactorum	17-21	162	124	178	178
P.cactorum	2003_3	157	122	172	172
P.cactorum	4032	163	128	179	179
P.cactorum	4040	158	124	174	174
P.cactorum	404	168	127	185	185
P.cactorum	415	184	141	203	203
P.cactorum	416	179	143	200	200
P.cactorum	62471	175	134	189	189
P.cactorum	P295	174	134	189	189
P.cactorum	P421_v2	163	126	179	179
P.cactorum	PC13_15	162	127	178	178
P.cactorum	R36_14	169	133	185	185
P.idaei	371	189	132	206	206
P.idaei	SCRP370	172	128	187	187
P.idaei	SCRP376	171	126	186	186
```

## 4.2.c Analysis of RxLR effectors - merger of Augustus / published genes with ORFs

Intersection between the coodinates of putative RxLRs from gene models and ORFs
were identified to determine the total number of RxLRs predicted in these
genomes.

The RxLR effectors from both Gene models and ORF finding approaches were
combined into a single file.

This step was complicated by the inconsistency in downloaded gff files for gene
models.


```bash
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/*/* | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$MergeDir/"$Strain"_total_RxLR.gff
AugTxt=$MergeDir/"$Strain"_total_RxLR_headers.txt
AugFa=$(ls gene_pred/codingquary/"$Species"/"$Strain"/final/final_genes_combined.pep.fasta)

ORFGff=$(ls $MergeDir/"$Strain"_total_ORF_RxLR.gff)
ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
ORFsTxt=$(ls $MergeDir/"$Strain"_total_ORF_RxLR_headers.txt)

ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_RxLR_EER_motif_hmm.gff
AugInORFs=$MergeDir/"$Strain"_AugInORFs_RxLR_EER_motif_hmm.gff
ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_RxLR_EER_motif_hmm.gff
AugUniq=$MergeDir/"$Strain"_Aug_Uniq_RxLR_EER_motif_hmm.gff
TotalRxLRsTxt=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.txt
TotalRxLRsGff=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.gff

bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq

# echo "$Species - $Strain"
# echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
CountORFsInAug=$(cat $ORFsInAug | grep -w 'gene' | wc -l)
# echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
CountAugInORFs=$(cat $AugInORFs | grep -w 'gene' | wc -l)
# echo "The number of RxLRs unique to ORF models:"
CountORFsUniq=$(cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l)
# cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
# echo "The number of RxLRs unique to Augustus models:"
CountAugUniq=$(cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l)
# echo "The total number of putative RxLRs are:"
cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalRxLRsTxt
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
cat $ORFsUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f3 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
CountTotalRxLR=$(cat $TotalRxLRsTxt | wc -l)
cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalRxLRsTxt > $TotalRxLRsGff

RxLRsFa=$MergeDir/"$Strain"_final_RxLR_EER.fa
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/unwrap_fasta.py --inp_fasta $AugFa | grep -A1 -w -f $AugTxt | grep -v -E '^--$' > $RxLRsFa
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalRxLRsTxt > $RxLRsFa
# echo "$Strain"
$ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalRxLRsTxt >> $RxLRsFa
# echo "$Strain"
# echo "The number of sequences extracted is"
CountExtractedRxLRs=$(cat $RxLRsFa | grep '>' | wc -l)
printf "$Species\t$Strain\t$CountORFsInAug\t$CountAugInORFs\t$CountORFsUniq\t$CountAugUniq\t$CountTotalRxLR\t$CountExtractedRxLRs\n"
done
```

```
P.cactorum	11-40	129	129	47	26	202	202
P.cactorum	12420	124	124	50	22	197	197
P.cactorum	15_13	120	120	54	24	199	199
P.cactorum	15_7	115	115	59	23	197	197
P.cactorum	17-21	136	136	42	29	207	207
P.cactorum	2003_3	128	128	44	25	198	198
P.cactorum	4032	122	122	57	22	203	203
P.cactorum	404	129	129	56	25	213	213
P.cactorum	4040	127	127	47	24	200	200
P.cactorum	415	126	126	77	26	229	229
P.cactorum	416	129	129	71	24	226	226
P.cactorum	62471	142	142	47	26	215	215
P.cactorum	P295	132	132	57	23	213	213
P.cactorum	P421_v2	126	126	53	25	205	205
P.cactorum	PC13_15	130	130	48	25	205	205
P.cactorum	R36_14	138	138	47	30	216	216
P.idaei	371	115	115	91	18	225	225
P.idaei	SCRP370	119	119	68	26	214	214
P.idaei	SCRP376	118	118	68	19	206	206
```




### H) From ORF gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in ORF gene models. This was done with the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/*/*/*.aa_cat.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
# Setting variables
Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
mkdir -p $OutDir
# Hmmer variables
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
# Searches for LFLAK domain
LFLAK_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm
HmmResultsLFLAK="$Strain"_ORF_CRN_LFLAK_unmerged_hmmer.txt
hmmsearch -T 0 $LFLAK_hmm $Proteome > $OutDir/$HmmResultsLFLAK
echo "Searching for LFLAK domains in: $Organism $Strain"
cat $OutDir/$HmmResultsLFLAK | grep 'Initial search space'
cat $OutDir/$HmmResultsLFLAK | grep 'number of targets reported over threshold'
HmmFastaLFLAK="$Strain"_ORF_CRN_LFLAK_unmerged_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResultsLFLAK $Proteome > $OutDir/$HmmFastaLFLAK
# Searches for DWL domain
DWL_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm
HmmResultsDWL="$Strain"_ORF_CRN_DWL_unmerged_hmmer.txt
hmmsearch -T 0 $DWL_hmm $Proteome > $OutDir/$HmmResultsDWL
echo "Searching for DWL domains in: $Organism $Strain"
cat $OutDir/$HmmResultsDWL | grep 'Initial search space'
cat $OutDir/$HmmResultsDWL | grep 'number of targets reported over threshold'
HmmFastaDWL="$Strain"_ORF_CRN_DWL_unmerged_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResultsDWL $Proteome > $OutDir/$HmmFastaDWL
# Identify ORFs found by both models
CommonHeaders=$OutDir/"$Strain"_ORF_CRN_DWL_LFLAK_unmerged_headers.txt
cat $OutDir/$HmmFastaLFLAK $OutDir/$HmmFastaDWL | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $CommonHeaders
echo "The number of CRNs common to both models are:"
cat $CommonHeaders | wc -l
# The sequences will be merged based upon the strength of their DWL domain score
# For this reason headers as they appear in the DWL fasta file were extracted
Headers=$OutDir/"$Strain"_CRN_hmmer_unmerged_headers.txt
cat $OutDir/$HmmFastaDWL | grep '>' | grep -w -f $CommonHeaders | tr -d '>' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | tr -d '-' | sed 's/hmm_score/HMM_score/g' > $Headers
# As we are dealing with JGI and Broad sequences, some features need formatting:
ORF_Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/*_ORF.gff3 | grep -v '_atg_')
# Gff features were extracted for each header
CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
# cat $Headers | cut -f1 > tmp.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $Headers $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
# $ProgDir/extract_gff_for_sigP_hits.pl tmp.txt $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
# Gff features were merged based upon the DWL hmm score
DbDir=analysis/databases/$Organism/$Strain
mkdir -p $DbDir
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $CRN_unmerged_Gff --db $DbDir/CRN_ORF.db
CRN_Merged_Gff=$OutDir/"$Strain"_CRN_merged_hmmer.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp $DbDir/CRN_ORF.db --id LFLAK_DWL_CRN --out $DbDir/CRN_ORF_merged.db --gff > $CRN_Merged_Gff
# Final results are reported:
echo "Number of CRN ORFs after merging:"
cat $CRN_Merged_Gff | grep 'gene' | wc -l
done
```

```
Searching for LFLAK domains in: P.cactorum 11-40
Initial search space (Z):             480478  [actual number of targets]
Domain search space  (domZ):             237  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 11-40
Initial search space (Z):             480478  [actual number of targets]
Domain search space  (domZ):             288  [number of targets reported over threshold]
The number of CRNs common to both models are:
139
Number of CRN ORFs after merging:
94
Searching for LFLAK domains in: P.cactorum 12420
Initial search space (Z):             481792  [actual number of targets]
Domain search space  (domZ):             252  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 12420
Initial search space (Z):             481792  [actual number of targets]
Domain search space  (domZ):             308  [number of targets reported over threshold]
The number of CRNs common to both models are:
148
Number of CRN ORFs after merging:
104
Searching for LFLAK domains in: P.cactorum 15_13
Initial search space (Z):             492365  [actual number of targets]
Domain search space  (domZ):             256  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 15_13
Initial search space (Z):             492365  [actual number of targets]
Domain search space  (domZ):             322  [number of targets reported over threshold]
The number of CRNs common to both models are:
159
Number of CRN ORFs after merging:
106
Searching for LFLAK domains in: P.cactorum 15_7
Initial search space (Z):             493630  [actual number of targets]
Domain search space  (domZ):             252  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 15_7
Initial search space (Z):             493630  [actual number of targets]
Domain search space  (domZ):             314  [number of targets reported over threshold]
The number of CRNs common to both models are:
150
Number of CRN ORFs after merging:
105
Searching for LFLAK domains in: P.cactorum 17-21
Initial search space (Z):             480195  [actual number of targets]
Domain search space  (domZ):             226  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 17-21
Initial search space (Z):             480195  [actual number of targets]
Domain search space  (domZ):             283  [number of targets reported over threshold]
The number of CRNs common to both models are:
135
Number of CRN ORFs after merging:
91
Searching for LFLAK domains in: P.cactorum 2003_3
Initial search space (Z):             492609  [actual number of targets]
Domain search space  (domZ):             252  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 2003_3
Initial search space (Z):             492609  [actual number of targets]
Domain search space  (domZ):             309  [number of targets reported over threshold]
The number of CRNs common to both models are:
150
Number of CRN ORFs after merging:
103
Searching for LFLAK domains in: P.cactorum 4032
Initial search space (Z):             497333  [actual number of targets]
Domain search space  (domZ):             253  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 4032
Initial search space (Z):             497333  [actual number of targets]
Domain search space  (domZ):             330  [number of targets reported over threshold]
The number of CRNs common to both models are:
157
Number of CRN ORFs after merging:
110
Searching for LFLAK domains in: P.cactorum 4040
Initial search space (Z):             490721  [actual number of targets]
Domain search space  (domZ):             257  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 4040
Initial search space (Z):             490721  [actual number of targets]
Domain search space  (domZ):             320  [number of targets reported over threshold]
The number of CRNs common to both models are:
156
Number of CRN ORFs after merging:
107
Searching for LFLAK domains in: P.cactorum 404
Initial search space (Z):             560045  [actual number of targets]
Domain search space  (domZ):             248  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 404
Initial search space (Z):             560045  [actual number of targets]
Domain search space  (domZ):             303  [number of targets reported over threshold]
The number of CRNs common to both models are:
147
Number of CRN ORFs after merging:
103
Searching for LFLAK domains in: P.cactorum 415
Initial search space (Z):             491086  [actual number of targets]
Domain search space  (domZ):             254  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 415
Initial search space (Z):             491086  [actual number of targets]
Domain search space  (domZ):             317  [number of targets reported over threshold]
The number of CRNs common to both models are:
154
Number of CRN ORFs after merging:
109
Searching for LFLAK domains in: P.cactorum 416
Initial search space (Z):             488742  [actual number of targets]
Domain search space  (domZ):             260  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 416
Initial search space (Z):             488742  [actual number of targets]
Domain search space  (domZ):             319  [number of targets reported over threshold]
The number of CRNs common to both models are:
157
Number of CRN ORFs after merging:
107
Searching for LFLAK domains in: P.cactorum 62471
Initial search space (Z):             493951  [actual number of targets]
Domain search space  (domZ):             249  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 62471
Initial search space (Z):             493951  [actual number of targets]
Domain search space  (domZ):             308  [number of targets reported over threshold]
The number of CRNs common to both models are:
153
Number of CRN ORFs after merging:
100
Searching for LFLAK domains in: P.cactorum P295
Initial search space (Z):             488637  [actual number of targets]
Domain search space  (domZ):             254  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum P295
Initial search space (Z):             488637  [actual number of targets]
Domain search space  (domZ):             327  [number of targets reported over threshold]
The number of CRNs common to both models are:
163
Number of CRN ORFs after merging:
109
Searching for LFLAK domains in: P.cactorum P421_v2
Initial search space (Z):             477118  [actual number of targets]
Domain search space  (domZ):             238  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum P421_v2
Initial search space (Z):             477118  [actual number of targets]
Domain search space  (domZ):             284  [number of targets reported over threshold]
The number of CRNs common to both models are:
142
Number of CRN ORFs after merging:
96
Searching for LFLAK domains in: P.cactorum PC13_15
Initial search space (Z):             496343  [actual number of targets]
Domain search space  (domZ):             251  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum PC13_15
Initial search space (Z):             496343  [actual number of targets]
Domain search space  (domZ):             321  [number of targets reported over threshold]
The number of CRNs common to both models are:
156
Number of CRN ORFs after merging:
109
Searching for LFLAK domains in: P.cactorum R36_14
Initial search space (Z):             499436  [actual number of targets]
Domain search space  (domZ):             260  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum R36_14
Initial search space (Z):             499436  [actual number of targets]
Domain search space  (domZ):             316  [number of targets reported over threshold]
The number of CRNs common to both models are:
159
Number of CRN ORFs after merging:
107
Searching for LFLAK domains in: P.idaei 371
Initial search space (Z):             487572  [actual number of targets]
Domain search space  (domZ):             217  [number of targets reported over threshold]
Searching for DWL domains in: P.idaei 371
Initial search space (Z):             487572  [actual number of targets]
Domain search space  (domZ):             254  [number of targets reported over threshold]
The number of CRNs common to both models are:
122
Number of CRN ORFs after merging:
86
Searching for LFLAK domains in: P.idaei SCRP370
Initial search space (Z):             487017  [actual number of targets]
Domain search space  (domZ):             220  [number of targets reported over threshold]
Searching for DWL domains in: P.idaei SCRP370
Initial search space (Z):             487017  [actual number of targets]
Domain search space  (domZ):             256  [number of targets reported over threshold]
The number of CRNs common to both models are:
125
Number of CRN ORFs after merging:
89
Searching for LFLAK domains in: P.idaei SCRP376
Initial search space (Z):             487331  [actual number of targets]
Domain search space  (domZ):             217  [number of targets reported over threshold]
Searching for DWL domains in: P.idaei SCRP376
Initial search space (Z):             487331  [actual number of targets]
Domain search space  (domZ):             260  [number of targets reported over threshold]
The number of CRNs common to both models are:
124
Number of CRN ORFs after merging:
89
```

Extract crinklers from ORFs and Braker/codingquary gene models


```bash
for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/*/* | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$(ls $MergeDir/"$Strain"_pub_CRN_LFLAK_DWL.gff)
AugFa=$(ls gene_pred/*/"$Species"/"$Strain"/final/final_genes_combined.pep.fasta)
ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer.gff3
ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_CRN_hmmer.bed
AugInORFs=$MergeDir/"$Strain"_AugInORFs_CRN_hmmer.bed
ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_CRN_hmmer.bed
AugUniq=$MergeDir/"$Strain"_Aug_Uniq_CRN_hmmer.bed
TotalCRNsTxt=$MergeDir/"$Strain"_final_CRN.txt
TotalCRNsGff=$MergeDir/"$Strain"_final_CRN.gff
TotalCRNsHeaders=$MergeDir/"$Strain"_Total_CRN_headers.txt
bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq
# echo "$Species - $Strain"

# echo "The number of ORF CRNs overlapping Augustus CRNs:"
CountORFsInAug=$(cat $ORFsInAug | grep -w -e 'transcript' -e 'mRNA' | wc -l)
# echo "The number of Augustus CRNs overlapping ORF CRNs:"
CountAugInORFs=$(cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA' | wc -l)
cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalCRNsTxt
# echo "The number of CRNs unique to ORF models:"
CountORFsUniq=$(cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' | wc -l)
cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt
# echo "The number of CRNs unique to Augustus models:"
CountAugUniq=$(cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l)
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt

cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalCRNsTxt > $TotalCRNsGff

CRNsFa=$MergeDir/"$Strain"_final_CRN.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalCRNsTxt > $CRNsFa
$ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalCRNsTxt >> $CRNsFa
# echo "The number of sequences extracted is"
CountNumExtracted=$(cat $CRNsFa | grep '>' | wc -l)
printf "$Species\t$Strain\t$CountORFsUniq\t$CountAugInORFs\t$CountORFsUniq\t$CountAugUniq\t$CountNumExtracted\n"
done
```

```
P.cactorum	11-40	40	54	40	5	99
P.cactorum	12420	42	62	42	4	108
P.cactorum	15_13	38	68	38	4	110
P.cactorum	15_7	40	67	40	4	111
P.cactorum	17-21	42	51	42	8	101
P.cactorum	2003_3	48	55	48	6	109
P.cactorum	4032	37	74	37	3	114
P.cactorum	404	57	49	57	4	110
P.cactorum	4040	41	66	41	7	114
P.cactorum	415	39	71	39	4	114
P.cactorum	416	41	70	41	8	119
P.cactorum	62471	48	52	48	9	109
P.cactorum	P295	44	65	44	8	117
P.cactorum	P421_v2	45	53	45	6	104
P.cactorum	PC13_15	40	70	40	4	114
P.cactorum	R36_14	37	70	37	11	118
P.idaei	371	46	41	46	1	88
P.idaei	SCRP370	60	30	60	0	90
P.idaei	SCRP376	63	26	63	0	89
```


# Making a combined file of Braker, Coding quary genes as well as additional ORF effector candidates


A gff file containing the combined Braker and CodingQuary genes as well as the
additional CRN and RxLR genes predicted by ORF analysis was made.

```bash
for GeneGff in $(ls gene_pred/*/P*/*/final/final_genes_appended.gff3 | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
GffOrfRxLR=$(ls analysis/RxLR_effectors/combined_evidence/$Organism/$Strain/*_ORFsUniq_RxLR_EER_motif_hmm.gff)
GffOrfCRN=$(ls analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/*_ORFsUniq_CRN_hmmer.bed)
Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
OutDir=gene_pred/final_incl_ORF/$Organism/$Strain
mkdir -p $OutDir
cat $GeneGff > $OutDir/${Strain}_genes_incl_ORFeffectors.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/add_ORF_features.pl $GffOrfRxLR $Assembly >> $OutDir/${Strain}_genes_incl_ORFeffectors.gff3
$ProgDir/add_ORF_features.pl $GffOrfCRN $Assembly >> $OutDir/${Strain}_genes_incl_ORFeffectors.gff3
# Make gene models from gff files.
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
# Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked_wrapped.fa)
$ProgDir/gff2fasta.pl $Assembly $OutDir/${Strain}_genes_incl_ORFeffectors.gff3 $OutDir/${Strain}_genes_incl_ORFeffectors
done
```

<!-- High confidence genes overlapped by ORFs to be inserted were manually inspected
to see whether the original gene model or the ORF should be kept. Preference was
given to the ORF effector candidate. Before this could be done bedtools was
used to identify which genes were intersected: -->

High confidence genes overlapped by ORFs to be inserted were manually inspected
to see whether the original gene model or the ORF should be kept. Manual inspection.
indicated that predicted gene models were, in general, represented the true gene.
As this approach was to be used for 20 genomes, not all conflicts were assessed
manually. Bedtools was used to identify which genes were intersected:

```bash
for GeneGff in $(ls gene_pred/*/P*/*/final/final_genes_appended.gff3 | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
  Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=gene_pred/final_incl_ORF/$Organism/$Strain
  GffOrfRxLR=$(ls analysis/RxLR_effectors/combined_evidence/$Organism/$Strain/${Strain}_ORFsUniq_RxLR_EER_motif_hmm.gff)
  GffOrfCRN=$(ls analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_ORFsUniq_CRN_hmmer.bed)
  Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  # Identify RxLR ORFs intersecting non-RxLR gene models
  bedtools intersect -wo -a $GeneGff -b $GffOrfRxLR | grep -e "AUGUSTUS.gene" -e "CodingQuarry_v2.0.gene" -e 'PGNCodingQuarry_v2.0' | grep "ORF_RxLR.gene"| cut -f1,9,18 | sed 's/ID=//g' | tr -d ';' > $OutDir/RxLR_ORFs_intersecting_non-RxLR_genes.txt
  # Identify CRN ORFs intersecting non-CRN gene models
  bedtools intersect -wo -a $GeneGff -b $GffOrfCRN | grep -e "AUGUSTUS.gene" -e "CodingQuarry_v2.0.gene" -e 'PGNCodingQuarry_v2.0' | grep "CRN_HMM.gene"| cut -f1,9,18 | sed 's/ID=//g' | tr -d ';' > $OutDir/CRN_ORFs_intersecting_non-CRN_genes.txt
done
```


Some ORF effectors were noted to overlap AugustusCodingQuarry gene models. These
were manually inspected in geneious and the following genes removed:

```bash
for CombinedGff in $(ls gene_pred/final_incl_ORF/*/*/*_genes_incl_ORFeffectors.gff3 | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | head -n1); do
Strain=$(echo $CombinedGff | rev | cut -f2 -d '/' | rev)
Organism=$(echo $CombinedGff | rev | cut -f3 -d '/' | rev)
# echo "$Organism - $Strain"
OutDir=$(dirname $CombinedGff)
cat $OutDir/RxLR_ORFs_intersecting_non-RxLR_genes.txt $OutDir/CRN_ORFs_intersecting_non-CRN_genes.txt | cut -f3 > $OutDir/exclude_list.txt
# echo "genes before filtering"
PreFilter=$(cat $CombinedGff | grep -w 'gene' | wc -l)
# echo "genes in list"
FilterList=$(cat $OutDir/exclude_list.txt | wc -l)
# echo "Unique genes in list"
UniqueFilterList=$(cat $OutDir/exclude_list.txt | sort | uniq | wc -l)
# echo "Number of genes removed"
GenesRemoved=$(cat $CombinedGff | grep -w -f $OutDir/exclude_list.txt | grep -w 'gene' | wc -l)
cat $CombinedGff | grep -w -v -f $OutDir/exclude_list.txt > $OutDir/${Strain}_genes_incl_ORFeffectors_filtered.gff3
# echo "final number of genes"
FinalGenes=$(cat $OutDir/${Strain}_genes_incl_ORFeffectors_filtered.gff3 | grep -w 'gene' | wc -l)
printf "$Organism\t$Strain\t$PreFilter\t$FilterList\t$UniqueFilterList\t$GenesRemoved\t$FinalGenes\n"
done
```

```
P.cactorum	11-40	28107	62	61	61	28046
P.cactorum	12420	25498	55	54	54	25444
P.cactorum	15_13	25996	54	54	54	25942
P.cactorum	15_7	25909	55	54	54	25855
P.cactorum	17-21	27926	58	57	57	27869
P.cactorum	2003_3	29014	66	64	64	28950
P.cactorum	4032	26113	52	51	51	26062
P.cactorum	4040	28727	58	57	57	28670
P.cactorum	404	35041	68	63	63	34978
P.cactorum	415	26102	61	60	60	26042
P.cactorum	416	28697	76	74	74	28623
P.cactorum	62471	28834	68	67	67	28767
P.cactorum	P295	25684	57	57	57	25627
P.cactorum	P421_v2	25629	63	63	63	25566
P.cactorum	PC13_15	25821	50	50	50	25771
P.cactorum	R36_14	29178	55	54	54	29124
P.idaei	371	24928	73	72	72	24856
P.idaei	SCRP370	27411	57	57	57	27354
P.idaei	SCRP376	24986	62	62	62	24924
```


```bash
 for GffAppended in $(ls gene_pred/final_incl_ORF/*/*/*_genes_incl_ORFeffectors_filtered.gff3 | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | head -n1); do
   Strain=$(echo $GffAppended | rev | cut -d '/' -f2 | rev)
   Organism=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
   echo "$Organism - $Strain"
   FinalDir=$(dirname $GffAppended)
   GffFiltered=$FinalDir/filtered_duplicates.gff
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
   $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
   GffRenamed=$FinalDir/final_genes_genes_incl_ORFeffectors_renamed.gff3
   LogFile=$FinalDir/final_genes_appended_renamed.log
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
   $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
   rm $GffFiltered
  Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
 $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_genes_incl_ORFeffectors_renamed
   # The proteins fasta file contains * instead of Xs for stop codons, these should
   # be changed
   sed -i 's/\*/X/g' $FinalDir/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta
 done
```
```bash
for Transcriptome in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Transcriptome | rev | cut -d '/' -f3 | rev)
GffFile=$(echo $Transcriptome | sed 's/.pep.fasta/.gff3/g')
Genes=$(cat $GffFile | grep -w 'gene' | wc -l)
Proteins=$(cat $Transcriptome | grep '>' | wc -l)
printf "$Organism\t$Strain\t$Genes\t$Proteins\n"
done
```


# Assessing gene space in predicted transcriptomes

```bash
for Transcriptome in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.gene.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Transcriptome | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/genes
qsub $ProgDir/sub_busco3.sh $Transcriptome $BuscoDB $OutDir
done
```

```
INFO    281 Complete BUSCOs (C)
INFO    272 Complete and single-copy BUSCOs (S)
INFO    9 Complete and duplicated BUSCOs (D)
INFO    7 Fragmented BUSCOs (F)
INFO    15 Missing BUSCOs (M)
INFO    303 Total BUSCO groups searched
```

```bash
 for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt); do
 Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
 Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
 Complete=$(cat $File | grep "(C)" | cut -f2)
 Single=$(cat $File | grep "(S)" | cut -f2)
 Fragmented=$(cat $File | grep "(F)" | cut -f2)
 Missing=$(cat $File | grep "(M)" | cut -f2)
 Total=$(cat $File | grep "Total" | cut -f2)
 echo -e "$Organism\t$Strain\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
 done
```


# Re-running Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	# cd /home/groups/harrisonlab/project_files/idris
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
  FileName=$(basename $Genes)
  SymLinkDir=$(dirname $Genes)/symlink_directory
  mkdir -p $SymLinkDir
  cp -s $PWD/$Genes $SymLinkDir/$FileName
	echo $Genes
	$ProgDir/sub_interproscan.sh $SymLinkDir/$FileName
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteome in $(ls gene_pred/final_incl_ORF/*/*/symlink_directory/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | grep -v '11-40'); do
Strain=$(echo $Proteome | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteome $InterProRaw
done
```


## B) SwissProt

```bash
 for Proteome in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
   Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   # Jobs=$(qstat | grep 'sub_swiss' | grep 'qw' | wc -l)
   # while [ $Jobs -gt 1 ]; do
   # sleep 1m
   # printf "."
   # Jobs=$(qstat | grep 'sub_swiss' | grep 'qw' | wc -l)
   # done		
   OutDir=gene_pred/swissprot/$Organism/$Strain
   SwissDbDir=../../uniprot/swissprot
   SwissDbName=uniprot_sprot
   ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
   qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
 done
```


# Genes with transcription factor annotations:


A list of PFAM domains, superfamily annotations used as part of the DBD database
and a further set of interproscan annotations listed by Shelest et al 2017 were made
http://www.transcriptionfactor.org/index.cgi?Domain+domain:all
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5415576/

```bash
for Interpro in $(ls gene_pred/interproscan/*/*/*_interproscan.tsv | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | grep -w '404'); do
Organism=$(echo $Interpro | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Interpro | rev | cut -f2 -d '/' | rev)
# echo "$Organism - $Strain"
OutDir=analysis/transcription_factors/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transcription_factors
$ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
# echo "total number of transcription factors"
cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
NumTF=$(cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l)
printf "$Organism\t$Strain\t$NumTF\n"
done
```

```
P.cactorum	11-40	968
P.cactorum	12420	904
P.cactorum	15_13	913
P.cactorum	15_7	921
P.cactorum	17-21	952
P.cactorum	2003_3	1011
P.cactorum	4032	924
P.cactorum	4040	978
P.cactorum	404	1064
P.cactorum	415	946
P.cactorum	416	1007
P.cactorum	62471	1003
P.cactorum	P295	928
P.cactorum	P421_v2	861
P.cactorum	PC13_15	932
P.cactorum	R36_14	1038
P.idaei	371	884
P.idaei	SCRP370	905
P.idaei	SCRP376	862
```



### Effector-like structure identification using EffectorP

Required programs:
* EffectorP.py

```bash
for Proteome in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
BaseName="$Organism"_"$Strain"_EffectorP
OutDir=analysis/effectorP/$Organism/$Strain
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
done
```


```bash
 for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
   Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
   Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
   # echo "$Organism - $Strain"
   Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
   cat $File | grep 'Effector' | cut -f1 > $Headers
   # printf "EffectorP headers:\t"
   NumEffP=$(cat $Headers | wc -l)
   printf "$Organism\t$Strain\t$NumEffP\n"
done
```

```
P.cactorum	11-40	19366
P.cactorum	12420	16432
P.cactorum	15_13	16682
P.cactorum	15_7	16710
P.cactorum	17-21	18894
P.cactorum	2003_3	20214
P.cactorum	4032	16916
P.cactorum	4040	20022
P.cactorum	404	26056
P.cactorum	415	17082
P.cactorum	416	19876
P.cactorum	62471	20048
P.cactorum	P295	16578
P.cactorum	P421_v2	16386
P.cactorum	PC13_15	16352
P.cactorum	R36_14	20220
P.idaei	371	15628
P.idaei	SCRP370	19020
P.idaei	SCRP376	16182
```

## CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/CAZY/$Organism/$Strain
mkdir -p $OutDir
Prefix="$Strain"_CAZY
CazyHmm=../../dbCAN/dbCAN-fam-HMMs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $File)
# echo "$Organism - $Strain"
ProgDir=/home/groups/harrisonlab/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
# printf "number of CAZY genes identified:\t"
CazymeNum=$(cat $CazyHeaders | wc -l)
printf "$Organism\t$Strain\t$CazymeNum\n"
Gff=$(ls gene_pred/final_incl_ORF/$Organism/$Strain/final_genes_genes_incl_ORFeffectors_renamed.gff3)
CazyGff=$OutDir/"$Strain"_CAZY.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

# SecretedProts=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_secreted.fa)
# # SecretedProts=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem_no_GPI.aa)
# SecretedHeaders=$(echo $SecretedProts | sed 's/.fa/_headers.txt/g' | sed 's/.aa/_headers.txt/g')
# cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
# CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
# $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
# printf "number of Secreted CAZY genes identified:\t"
# cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
# cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
done
```

```
P.cactorum	11-40	641
P.cactorum	12420	644
P.cactorum	15_13	640
P.cactorum	15_7	607
P.cactorum	17-21	646
P.cactorum	2003_3	623
P.cactorum	4032	667
P.cactorum	4040	646
P.cactorum	404	800
P.cactorum	415	624
P.cactorum	416	628
P.cactorum	62471	652
P.cactorum	P295	622
P.cactorum	P421_v2	647
P.cactorum	PC13_15	662
P.cactorum	R36_14	652
P.idaei	371	598
P.idaei	SCRP370	565
P.idaei	SCRP376	588
```
<!--
or, if those secreted proteins with TM domains or GPI anchors are removed:

```
P.cactorum - 10300
number of CAZY genes identified:	697
number of Secreted CAZY genes identified:	236
``` -->


Note - the CAZY genes identified may need further filtering based on e value and
cuttoff length - see below:

Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage

* For fungi, use E-value < 1e-17 and coverage > 0.45

* The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)

### Summary of CAZY families by organism


The following commands were used to identify b-glucan binding proteins as described in:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3627985/
<!--
```bash
 for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
   cat $File | grep -w -e 'GH5' -e 'GH55' -e 'GH81' -e 'GH16' -e 'GH3' -e 'GH17' -e 'GH72' | sed -r "s/\s+/\t/g" | cut -f1,4 |  sort | uniq > tmp1.txt
   cat tmp1.txt | cut -f2 > tmp2.txt
   cat tmp1.txt | grep -f tmp2.txt | cut -f1 | sort | uniq -c
 done
```

```
  28 GH16.hmm
  30 GH17.hmm
  22 GH3.hmm
  55 GH5.hmm
  13 GH72.hmm
  21 GH81.hmm
``` -->

## PhiBase genes:


Identifying PHIbase homologs

```bash
qlogin -pe smp 6
cd /home/groups/harrisonlab/project_files/idris
for QueryFasta in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
Strain=$(echo $QueryFasta | rev | cut -f2 -d '/' | rev)
Organism=$(echo $QueryFasta | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/blast_homology/$Organism/$Strain
mkdir -p $OutDir
dbFasta=$(ls /home/groups/harrisonlab/phibase/v4.4/phi_accessions.fa)
dbType="prot"
Prefix="${Strain}_phi_accessions"
Eval="1e-30"
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
blastx -num_threads 6 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
cat $OutDir/${Prefix}_hits.txt | grep 'effector' | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
done > log.txt
```


# Building Annotation Tables

```bash
OrthoIDs=analysis/orthology/orthomcl/Pcac_Pinf_publication/Pcac_Pinf_publication_isolateIDs.txt
printf \
"
Pc_CR1 414
Pc_CR2 12420
Pc_CR3 15_13
Pc_CR4 15_7
Pc_CR5 2003_3
Pc_CR6 4032
Pc_CR7 4040
Pc_CR8 404
Pc_CR9 415
Pc_CR10 416
Pc_CR11 P421_v2
Pc_CR12 PC13_15
Pc_CR13 10300
Pc_LR1 11-40
Pc_LR2 17-21
Pc_MD1 62471
Pc_MD2 P295
Pc_MD3 R36_14
Pi_RI1 371
Pi_RI2 SCRP370
Pi_RI3 SCRP376
" \
> $OrthoIDs

```


```bash
for GeneGff in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.gff3 | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300' | grep -v -e '415' -e '416' | grep '62471'); do
Strain=$(echo $GeneGff | rev | cut -f2 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
OutDir=gene_pred/annotation/$Organism/$Strain
mkdir -p $OutDir

# combine sigP2 headers
cat gene_pred/final_sigP/$Organism/$Strain/${Strain}_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' > $OutDir/sigP2_headers_combined.txt
cat gene_pred/ORF_sigP/$Organism/$Strain/${Strain}_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' >> $OutDir/sigP2_headers_combined.txt
# cat gene_pred/ORF_sigP/$Organism/$Strain/${Strain}_ORF_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' >> $OutDir/sigP2_headers_combined.txt
# cat gene_pred/ORF_sigP/$Organism/${Strain}_old/${Strain}_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' >> $OutDir/sigP2_headers_combined.txt

# combine sigP3 headers
cat gene_pred/final_signalp-3.0/$Organism/$Strain/${Strain}_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' '  > $OutDir/sigP3_headers_combined.txt
# cat gene_pred/ORF_signalp-3.0/$Organism/$Strain/${Strain}_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' >> $OutDir/sigP3_headers_combined.txt

# combine sigP4 headers
cat gene_pred/final_signalp-4.1/$Organism/$Strain/${Strain}_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' '  > $OutDir/sigP4_headers_combined.txt
cat gene_pred/ORF_signalp-4.1/$Organism/$Strain/${Strain}_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' >> $OutDir/sigP4_headers_combined.txt

cat analysis/phobius/$Organism/$Strain/${Strain}_phobius_headers.txt analysis/phobius/$Organism/$Strain/${Strain}_phobius_headers_ORF.txt | sed "s/^g_/contig_/g" > $OutDir/phobius_headers_combined.txt

cat analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/${Strain}_RxLR_EER_regex.fa analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/${Strain}_ORF_RxLR_EER_regex_merged.aa > $OutDir/${Strain}_combined_RxLR_EER_regex.fa

cat analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/${Strain}_RxLR_hmmer.fa analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/${Strain}_ORF_RxLR_hmmer.fa > $OutDir/${Strain}_combined_RxLR_hmmer.fa

cat analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain/${Strain}_Aug_WY_hmmer.fa analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain/${Strain}_ORF_WY_hmmer.fa > $OutDir/${Strain}_combined_WY_hmmer.fa

cat analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_ORFsUniq_CRN_hmmer.bed \
| grep -w 'transcript' | cut -f5 -d '=' | cut -f1 -d ';' \
> analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_ORFsUniq_CRN_hmmer.txt
cat analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_ORF_CRN_LFLAK_unmerged_hmmer.fa \
| grep -A1 -f analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_ORFsUniq_CRN_hmmer.txt \
| grep -v '\-\-' \
> analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_ORFsUniq_CRN_hmmer.fa

cat analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_pub_CRN_LFLAK_hmm.fa analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_ORFsUniq_CRN_hmmer.fa > $OutDir/${Strain}_combined_CRN_LFLAK_hmm.fa

cat analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_pub_CRN_DWL_hmm.fa analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_ORF_CRN_DWL_unmerged_hmmer.fa > $OutDir/${Strain}_combined_CRN_DWL_hmm.fa

Fasta=$(ls gene_pred/final_incl_ORF/$Organism/$Strain/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta)
OrfGff=$(ls gene_pred/ORF_finder/$Organism/$Strain/${Strain}_ORF.gff3)
GeneConversion=$(ls gene_pred/final_incl_ORF/$Organism/$Strain/final_genes_appended_renamed.log)
TFs=$(ls analysis/transcription_factors/$Organism/$Strain/"$Strain"_TF_domains.tsv)
SigP2=$(ls $OutDir/sigP2_headers_combined.txt)
SigP3=$(ls $OutDir/sigP3_headers_combined.txt)
SigP4=$(ls $OutDir/sigP4_headers_combined.txt)
Phobius=$(ls $OutDir/phobius_headers_combined.txt)
TM_out=$(ls gene_pred/trans_mem/$Organism/$Strain/"$Strain"_TM_genes_pos.txt)
# GPI_out=$(ls gene_pred/trans_mem/$Organism/$Strain/GPIsom/GPI_pos.fa)
RxLR_Motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/${Strain}_all_secreted_RxLR_regex.fa)
RxLR_Hmm=$(ls $OutDir/${Strain}_combined_RxLR_hmmer.fa)
RxLR_WY=$(ls $OutDir/${Strain}_combined_WY_hmmer.fa)
CRN_LFLAK=$(ls $OutDir/${Strain}_combined_CRN_LFLAK_hmm.fa)
CRN_DWL=$(ls $OutDir/${Strain}_combined_CRN_DWL_hmm.fa)
EffP_list=$(ls analysis/effectorP/$Organism/$Strain/"$Organism"_"$Strain"_EffectorP_headers.txt)
CAZY_list=$(ls gene_pred/CAZY/$Organism/$Strain/"$Strain"_CAZY.out.dm.ps)
PhiHits=$(ls analysis/blast_homology/$Organism/$Strain/"$Strain"_phi_accessions_hits_headers.txt)
InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vMar2018_tophit_parsed.tbl)
Orthology=$(ls /home/groups/harrisonlab/project_files/idris/analysis/orthology/orthomcl/Pcac_Pinf_publication/Pcac_Pinf_publication_orthogroups.txt)
OrthoIDs=$(ls analysis/orthology/orthomcl/Pcac_Pinf_publication/Pcac_Pinf_publication_isolateIDs.txt)
OrthoStrainID=$(cat $OrthoIDs | grep -w "$Strain" | cut -f1 -d ' ')
echo "$Strain - $OrthoStrainID"
OrthoStrainAll='Pc_CR1 Pc_CR2 Pc_CR3 Pc_CR4 Pc_CR5 Pc_CR6 Pc_CR7 Pc_CR8 Pc_CR9 Pc_CR10 Pc_CR11 Pc_CR12 Pc_CR13 Pc_LR1 Pc_LR2 Pc_MD1 Pc_MD2 Pc_MD3 Pi_RI1 Pi_RI2 Pi_RI3'

ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/gene_annotation
$ProgDir/illumina_annotation_table.py \
--protein_fasta $Fasta \
--conversion_log $GeneConversion \
--ORF_gff $OrfGff \
--genes_gff $GeneGff \
--SigP2 $SigP2 \
--SigP3 $SigP3 \
--SigP4 $SigP4 \
--phobius $Phobius \
--TM_list $TM_out \
--RxLR_regex $RxLR_Motif \
--RxLR_hmm $RxLR_Hmm \
--CRN_dwl $CRN_DWL \
--CRN_lflak $CRN_LFLAK \
--EffP_list $EffP_list \
--PhiHits $PhiHits \
--CAZY $CAZY_list \
--TFs $TFs \
--InterPro $InterPro \
--Swissprot $SwissProt \
--orthogroups $Orthology \
--strain_id $OrthoStrainID  \
--OrthoMCL_all $OrthoStrainAll \
> $OutDir/"$Strain"_annotation_ncbi2.tsv
done
```

Analysing annotation tables

```bash
cat $OutDir/"$Strain"_annotation_ncbi.tsv | cut -f25 | grep -v "^$" | grep -v 'ipr_effectors' | sort | uniq -c | sort -nr
```


```bash
for File in $(ls gene_pred/annotation/*/*/*_annotation_ncbi.tsv | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  GeneNum=$(cat $File | cut -f1 | tail -n+2 | cut -f1 -d '.' | uniq | wc -l)
  ProtNum=$(cat $File | cut -f1 | tail -n+2 | uniq | wc -l)
  Secreted=$(cat $File | cut -f11 | tail -n+2 | grep 'Yes' | wc -l)
  EffP=$(cat $File | cut -f11,18 | tail -n+2 | grep "Yes.Yes" | wc -l)
  Cazy=$(cat $File | cut -f11,19 | tail -n+2 | grep "Yes.CAZY"| wc -l)
  RxLR=$(cat $File | cut -f14 | tail -n+2 | grep "RxLR" |wc -l)
  CRN=$(cat $File | cut -f17 | tail -n+2 | grep "CRN" | wc -l)
  TFs=$(cat $File | cut -f20 | tail -n+2 | grep -v "^$" | wc -l)
  Elicitin=$(cat $File | cut -f11,25 | tail -n+2 | grep 'Yes' | grep "MAMP:Elicitin" | wc -l)
  Transglutaminase=$(cat $File | cut -f11,25 | tail -n+2 | grep 'Yes' | grep "MAMP:Transglutaminase" | wc -l)
  NLP=$(cat $File | cut -f11,25 | tail -n+2 | grep 'Yes' | grep "Apoplastic:NLP" | wc -l)
  Kazal=$(cat $File | cut -f11,25 | tail -n+2 | grep 'Yes' | grep "Apoplastic:Protease inhibitor (Kazal-type)" | wc -l)
  Cathepsin=$(cat $File | cut -f11,25 | tail -n+2 | grep 'Yes' | grep "Apoplastic:Protease inhibitor (cathepsin)" | wc -l)
  Cystatin=$(cat $File | cut -f11,25 | tail -n+2 | grep 'Yes' | grep "Apoplastic:Protease inhibitor (cystatin-like)" | wc -l)
  GlucanaseInhibitor=$(cat $File | cut -f11,25 | grep 'Yes' | tail -n+2 | grep "Apoplastic:Glucanase inhibitor" | wc -l)
  Phytotoxin=$(cat $File | cut -f11,25 | tail -n+2 | grep 'Yes' | grep "Apoplastic:Phytotoxin" | wc -l)
  Cutinase=$(cat $File | cut -f11,25 | tail -n+2 | grep 'Yes' | grep "Apoplastic:Cutinase" | wc -l)
  printf "$Organism\t$Strain\t$GeneNum\t$ProtNum\t$Secreted\t$EffP\t$Cazy\t$RxLR\t$CRN"
  printf "\t$TFs\t$Elicitin\t$Transglutaminase\t$NLP\t$Kazal\t$Cathepsin\t$Cystatin\t$GlucanaseInhibitor\t$Phytotoxin\t$Cutinase\n"
done
```


Apple expansion

```bash
AnnotTab=$(ls gene_pred/annotation/P.cactorum/62471/62471_annotation_ncbi2.tsv)
OutDir=gene_pred/annotation/P.cactorum/62471
Strain=62471
cat $AnnotTab | grep -w -e 'apple expanded' -e 'crown rot loss'  | cut -f21 | sort | uniq | wc -l
cat $AnnotTab | grep -w -e 'apple expanded' -e 'crown rot loss' | wc -l
cat $AnnotTab | grep -w -e 'apple expanded'| cut -f21 | sort | uniq | wc -l
cat $AnnotTab | grep -w -e 'apple expanded' | wc -l
cat $AnnotTab | grep -w -e 'crown rot loss' | cut -f21 | sort | uniq | wc -l
cat $AnnotTab | grep -w -e 'crown rot loss' | wc -l

cat $AnnotTab | grep -w -e 'apple expanded' > $OutDir/"$Strain"_annotation_ncbi_md_expanded.tsv
cat $OutDir/"$Strain"_annotation_ncbi_md_expanded.tsv | grep 'RxLR' | wc -l
cat $OutDir/"$Strain"_annotation_ncbi_md_expanded.tsv | cut -f28 | wc -l

cat $AnnotTab | grep -w -e 'crown rot loss' > $OutDir/"$Strain"_annotation_ncbi_cr_loss.tsv
cat $OutDir/"$Strain"_annotation_ncbi_cr_loss.tsv | grep 'RxLR' | wc -l
```
