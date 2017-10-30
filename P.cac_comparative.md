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
```

#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:


```bash
  # for RawData in $(ls raw_dna/paired/P.*/*/*/*.fastq.gz | grep '170210'); do
  for RawData in $(ls raw_dna/paired/P.cactorum/*/*/*.fastq.gz | grep -v '10300' | grep -e '11-40' -e '17-21'); do
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
for StrainPath in $(ls -d raw_dna/paired/P.*/* | grep -v -w -e '10300' -e '404' -e '414' -e '415' -e '416' -e 'PC13_15' -e '2003_3' | grep -e 'P.cactorum' -e 'P.idaei' | grep -e '11-40' -e '17-21'); do
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
for TrimData in $(ls qc_dna/paired/P.cactorum/*/*/*.fq.gz | grep -e '11-40' -e '17-21'); do
echo $TrimData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $TrimData;
done
```

Sequencing coveraqge was estimated:

```bash
for RawData in $(ls qc_dna/paired/P.*/*/*/*q.gz | grep -e '11-40' -e '17-21'); do
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
for StrainDir in $(ls -d qc_dna/paired/P.*/* | grep -v -w -e '411' -e '10300_old' | grep -e 'cactorum' -e 'idaei' | grep -e '11-40' -e '17-21'); do
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
for StrainPath in $(ls -d qc_dna/paired/P.*/* | grep -v -w -e '10300' -e '404' -e '414' -e '415' -e '416' -e 'PC13_15' -e '2003_3'| grep -e 'P.cactorum' -e 'P.idaei'); do
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
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs*/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep 'SCRP376'); do
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
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for Exclude_db in "bacillus" "delftia" "Paen"; do
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
for Exclude_db in "bacillus" "delftia" "Paen"; do
  echo $Exclude_db
  for File in $(ls assembly/spades/P.*/*/*/log.txt | grep "$Exclude_db"); do
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

Those contigs that were not identified as contaminants in each of the filtering
steps were retained.

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
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

Assembly stats were collected on filtered assemblies:

```bash
  # for Assembly in $(ls assembly/spades/P.*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
  for Assembly in $(ls assembly/spades/P.*/*/deconseq_appended/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
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
for Assembly in $(ls assembly/spades/P.*/*/deconseq_appended/report.tsv); do  
Strain=$(echo $Assembly | rev | cut -f3 -d '/'| rev);
Size=$(cat $Assembly | grep 'Total length' | head -n1 | cut -f2);
OldAssembly=$(ls assembly/spades/P.*/$Strain/filtered_contigs*/report.tsv)
OldSize=$(cat $OldAssembly | grep 'Total length' | head -n1 | cut -f2);
printf "$Strain\t$Size\t$OldSize\n";
done
```

```
  12420	60315503	66488881
  15_13	60703463	60708976
  15_7	60797060	60806527
  2003_3	60800281	67797455
  4032	61561706	70105079
  4040	60657148	65177765
  404	75508939	82472226
  414	61116855	62710907
  415	60677782	64884339
  416	60502157	64802670
  62471	61265225	62823426
  P295	60707114	67293243
  PC13_15	61130682	67857120
  R36_14	61656467	68350251
  371	60474213	60483665
  SCRP370	60582725	60561664
  SCRP376	60610468	60615981
```

Numbers of busco genes in each assembly were identified:

```bash
# for Assembly in $(ls assembly/spades/P.*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
for Assembly in $(ls assembly/spades/P.*/*/deconseq_appended/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB="Eukaryotic"
# OutDir=gene_pred/busco/$Organism/$Strain/assembly
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
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
for Assembly in $(ls assembly/spades/P.*/*/deconseq_appended/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
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
for Assembly in $(ls assembly/spades/P.*/*/deconseq_appended/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep '12420'); do
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
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
    Kmer=$(echo $Assembly | rev | cut -f2 -d '/' | rev);
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    # OutDir=assembly/spades/$Organism/$Strain/filtered_contigs;
    OutDir=$(dirname $Assembly)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir;
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
for BestAss in $(ls assembly/spades/P.*/*/*/contigs_min_500bp_renamed.fasta | grep -e 'P.idaei' -e 'P.cactorum' | grep -e '2003_3'); do
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
for RepDir in $(ls -d repeat_masked/P.*/*/filtered_contigs_repmask | grep -e 'P.cactorum' -e 'P.idaei' | grep -v '12420'); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
printf "$Organism\t$Strain\n"
printf "The number of bases masked by RepeatMasker:\t"
sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The number of bases masked by TransposonPSI:\t"
sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The total number of masked bases are:\t"
cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
echo
done
```

```
  P.cactorum	404
  The number of bases masked by RepeatMasker:	14847140
  The number of bases masked by TransposonPSI:	4124170
  The total number of masked bases are:	15712478

  P.cactorum	414
  The number of bases masked by RepeatMasker:	16199492
  The number of bases masked by TransposonPSI:	4951212
  The total number of masked bases are:	17347836

  P.cactorum	415
  The number of bases masked by RepeatMasker:	13217011
  The number of bases masked by TransposonPSI:	4691282
  The total number of masked bases are:	14475772

  P.cactorum	416
  The number of bases masked by RepeatMasker:	15190662
  The number of bases masked by TransposonPSI:	5235395
  The total number of masked bases are:	16481354

  P.cactorum	62471
  The number of bases masked by RepeatMasker:	14879505
  The number of bases masked by TransposonPSI:	4104109
  The total number of masked bases are:	15851505

  P.idaei	371
  The number of bases masked by RepeatMasker:	15158762
  The number of bases masked by TransposonPSI:	4133282
  The total number of masked bases are:	16190027

  P.idaei	SCRP370
  The number of bases masked by RepeatMasker:	15585647
  The number of bases masked by TransposonPSI:	4125579
  The total number of masked bases are:	16524671
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
<!--
```bash
  # for Genome in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v 'P.cactorum'); do
  for Genome in $(ls assembly/spades/P.*/*/filtered_contigs_repmask/contigs_min_500bp_renamed.fasta | grep -e 'P.cactorum'); do
    Strain=$(echo $Genome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Genome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/cegma/$Organism/$Strain
    Prefix="$Strain"_dna
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
    qsub $ProgDir/sub_cegma.sh $Genome $Prefix $OutDir
  done
```

```bash
  for Genome in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Genome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Genome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/cegma/$Organism/$Strain
    Prefix="$Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
    qsub $ProgDir/sub_cegma.sh $Genome $Prefix $OutDir
  done
```

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/P*/*/*_dna_cegma.completeness_report); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
``` -->

Busco has replaced CEGMA and was run to check gene space in assemblies

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
# Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
# Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB="Eukaryotic"
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt | grep -e 'P.cactorum' -e 'P.idaei'); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
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


#### Aligning

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e 'P.fragariae' | grep 'A4'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNA in $(ls qc_rna/raw_rna/genbank/*/*/*_trim.fq.gz); do
      Timepoint=$(echo $RNA | rev | cut -f1 -d '/' | rev | sed 's/_trim.*//g')
      echo "$Timepoint"
      OutDir=alignment/$Organism/$Strain/$Timepoint
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
      qsub $ProgDir/tophat_alignment_unpaired.sh $Assembly $RNA $OutDir
    done
  done
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
    for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep 'A4'); do
    	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    	echo "$Organism - $Strain"
    	mkdir -p alignment/$Organism/$Strain/concatenated
    	samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
      alignment/$Organism/$Strain/SRR1206032/accepted_hits.bam \
    	alignment/$Organism/$Strain/SRR1206033/accepted_hits.bam
    	# OutDir=gene_pred/braker/$Organism/"$Strain"_braker
    	# AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
    	# GeneModelName="$Organism"_"$Strain"_braker
    	# rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    	# ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
    	# qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  	done
```


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa | grep 'A4'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```

Secondly, genes were predicted using CodingQuary:

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa | grep 'A4'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/codingquary/$Organism/$Strain
CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
  for BrakerGff in $(ls gene_pred/braker/P.*/*_braker/*/augustus.gff3 | grep -v -e '414' -e 'fragariae' | grep 'P.idaei'); do
    Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g' | sed 's/_braker//g')
    Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/"$Strain"_contigs_softmasked.fa)
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
  for DirPath in $(ls -d gene_pred/codingquary/P.*/*/final | grep -v -w -e '414'); do
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
  P.cactorum - 404
  24296
  3479
  27775

  P.cactorum - 415
  28984
  2960
  31944

  P.cactorum - 416
  29746
  3034
  32780

  P.cactorum - 62471
  24682
  2509
  27191

  P.idaei - 371
  24667
  2586
  27253

  P.idaei - SCRP370
  24421
  2562
  26983
```


## Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing
Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	for Genome in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -e 'cactorum' -e 'idaei'); do
    echo "$Genome"
  	qsub $ProgDir/run_ORF_finder.sh $Genome
  done
```

The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
	ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
  for OrfGff in $(ls gene_pred/ORF_finder/P.*/*/*_ORF.gff | grep -v 'atg' | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '10300'); do
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
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -v -w -e '414' | grep 'P.idaei'); do
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
      qsub $ProgDir/pred_sigP.sh $File
      qsub $ProgDir/pred_sigP.sh $File signalp-4.1
    done
  done
 ```

 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

 ```bash
  for SplitDir in $(ls -d gene_pred/final_split/P.*/* | grep -v '414'); do
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

#### B.2) Prediction using Phobius

 Secreted proteins were also predicted using Phobius

 ```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -v -w -e '414' | grep 'P.idaei'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/phobius/$Organism/$Strain
    mkdir -p $OutDir
    phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    $ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
  done
 ```


Secreted proteins from different sources were combined into a single file:

```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta); do
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
  P.cactorum - 404
  The following number of sequences were predicted as secreted:
  9463
  This represented the following number of unique genes:
  3177
  P.cactorum - 414
  The following number of sequences were predicted as secreted:
  6891
  This represented the following number of unique genes:
  4069
  P.cactorum - 415
  The following number of sequences were predicted as secreted:
  10500
  This represented the following number of unique genes:
  3557
  P.cactorum - 416
  The following number of sequences were predicted as secreted:
  10342
  This represented the following number of unique genes:
  3543
  P.cactorum - 62471
  The following number of sequences were predicted as secreted:
  9456
  This represented the following number of unique genes:
  3150
  P.idaei - 371
  The following number of sequences were predicted as secreted:
  8813
  This represented the following number of unique genes:
  3002
  P.idaei - SCRP370
  The following number of sequences were predicted as secreted:
  8792
  This represented the following number of unique genes:
  3006
```


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
    OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
    OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
    printf "Secreted effectorP headers:\t"
    cat $OutFileHeaders | wc -l
    Gff=$(ls gene_pred/codingquary/$Organism/$Strain/*/final_genes_appended.gff3)
    EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
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


# D) Prediction of RxLRs from - Augustus/CodingQuary gene models

 The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

 The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash

for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep -w -e '414' -e 'P.cactorum' -e 'P.idaei' | grep -v '10300'); do
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

strain: 404	species: P.cactorum
the total number of SigP gene is:	9463
the number of unique SigP gene is:	3177
the number of SigP-RxLR genes are:	272
the number of SigP-RxLR-EER genes are:	118


strain: 414	species: P.cactorum
the total number of SigP gene is:	6891
the number of unique SigP gene is:	4069
the number of SigP-RxLR genes are:	325
the number of SigP-RxLR-EER genes are:	147


strain: 415	species: P.cactorum
the total number of SigP gene is:	10500
the number of unique SigP gene is:	3557
the number of SigP-RxLR genes are:	281
the number of SigP-RxLR-EER genes are:	125


strain: 416	species: P.cactorum
the total number of SigP gene is:	10342
the number of unique SigP gene is:	3543
the number of SigP-RxLR genes are:	278
the number of SigP-RxLR-EER genes are:	122


strain: 62471	species: P.cactorum
the total number of SigP gene is:	9456
the number of unique SigP gene is:	3150
the number of SigP-RxLR genes are:	287
the number of SigP-RxLR-EER genes are:	134


strain: 371	species: P.idaei
the total number of SigP gene is:	8813
the number of unique SigP gene is:	3002
the number of SigP-RxLR genes are:	247
the number of SigP-RxLR-EER genes are:	112


strain: SCRP370	species: P.idaei
the total number of SigP gene is:	8792
the number of unique SigP gene is:	3006
the number of SigP-RxLR genes are:	238
the number of SigP-RxLR-EER genes are:	111
```


### G) From Secreted gene models - Hmm evidence of RxLR effectors

```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414' -e 'idaei'); do
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
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_RxLR_EER_regex.txt | grep -v -e 'Aug' -e '10300' | grep -e 'P.idaei' -e 'P.cactorum' | grep -v '414'); do
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
  P.cactorum - 404
  Number of RxLRs identified by Regex:
  118
  Number of RxLRs identified by Hmm:
  127
  Number of RxLRs in combined dataset:
  149
  P.cactorum - 414
  Number of RxLRs identified by Regex:
  146
  Number of RxLRs identified by Hmm:
  145
  Number of RxLRs in combined dataset:
  173
  P.cactorum - 415
  Number of RxLRs identified by Regex:
  125
  Number of RxLRs identified by Hmm:
  132
  Number of RxLRs in combined dataset:
  159
  P.cactorum - 416
  Number of RxLRs identified by Regex:
  122
  Number of RxLRs identified by Hmm:
  129
  Number of RxLRs in combined dataset:
  155
  P.cactorum - 62471
  Number of RxLRs identified by Regex:
  134
  Number of RxLRs identified by Hmm:
  142
  Number of RxLRs in combined dataset:
  169
  P.idaei - 371
  Number of RxLRs identified by Regex:
  112
  Number of RxLRs identified by Hmm:
  107
  Number of RxLRs in combined dataset:
  136
  P.idaei - SCRP370
  Number of RxLRs identified by Regex:
  111
  Number of RxLRs identified by Hmm:
  105
  Number of RxLRs in combined dataset:
  137

```
72 of 137 RxLRs were in the predicted 797 effectorP genes.
```bash
  cat analysis/RxLR_effectors/combined_evidence/P.idaei/SCRP370/SCRP370_total_RxLR_headers.txt $OutFileHeaders | sort | cut -f1 -d '.' | uniq -d | wc -l
```

### D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:


```bash
  HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
  LFLAK_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm)
  DWL_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm)
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414' -e 'P.idaei'); do
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
  P.cactorum - 404
  Initial search space (Z):              27775  [actual number of targets]
  Domain search space  (domZ):             156  [number of targets reported over threshold]
  Initial search space (Z):              27775  [actual number of targets]
  Domain search space  (domZ):             129  [number of targets reported over threshold]
  113
  P.cactorum - 414
  Initial search space (Z):              32832  [actual number of targets]
  Domain search space  (domZ):             218  [number of targets reported over threshold]
  Initial search space (Z):              32832  [actual number of targets]
  Domain search space  (domZ):             189  [number of targets reported over threshold]
  170
  P.cactorum - 415
  Initial search space (Z):              31944  [actual number of targets]
  Domain search space  (domZ):             201  [number of targets reported over threshold]
  Initial search space (Z):              31944  [actual number of targets]
  Domain search space  (domZ):             169  [number of targets reported over threshold]
  147
  P.cactorum - 416
  Initial search space (Z):              32780  [actual number of targets]
  Domain search space  (domZ):             206  [number of targets reported over threshold]
  Initial search space (Z):              32780  [actual number of targets]
  Domain search space  (domZ):             182  [number of targets reported over threshold]
  155
  P.cactorum - 62471
  Initial search space (Z):              27191  [actual number of targets]
  Domain search space  (domZ):             154  [number of targets reported over threshold]
  Initial search space (Z):              27191  [actual number of targets]
  Domain search space  (domZ):             132  [number of targets reported over threshold]
  112
  P.idaei - 371
  Initial search space (Z):              27253  [actual number of targets]
  Domain search space  (domZ):             127  [number of targets reported over threshold]
  Initial search space (Z):              27253  [actual number of targets]
  Domain search space  (domZ):              98  [number of targets reported over threshold]
  87
  P.idaei - SCRP370
  Initial search space (Z):              26983  [actual number of targets]
  Domain search space  (domZ):             129  [number of targets reported over threshold]
  Initial search space (Z):              26983  [actual number of targets]
  Domain search space  (domZ):             103  [number of targets reported over threshold]
  92
```

23 of 92 P.idaei CRNs were in the effectorP dataset.
```bash
cat analysis/CRN_effectors/hmmer_CRN/P.idaei/SCRP370/SCRP370_pub_CRN_LFLAK_DWL.txt $OutFileHeaders | sort | cut -f1 -d '.' | uniq -d | wc -l
```

Extract gff annotations for Crinklers:

```bash
  for CRNlist in $(ls analysis/CRN_effectors/hmmer_CRN/*/*/*_pub_CRN_LFLAK_DWL.txt | grep -e 'P.idaei' -e 'P.cactorum' | grep -v -e '10300'); do
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
	for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '10300'); do
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
  		qsub $ProgDir/pred_sigP.sh $File signalp-4.1
  	done
  done
```


 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

```bash
for SplitDir in $(ls -d gene_pred/ORF_split/*/* | grep -w -e 'P.cactorum' -e 'P.idaei' | grep -v -e '10300' -e '414_v2'); do
Strain=$(echo $SplitDir | cut -d '/' -f4)
Organism=$(echo $SplitDir | cut -d '/' -f3)
echo "$Organism - $Strain"
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
for SigpDir in $(ls -d gene_pred/ORF_sig* | cut -f2 -d'/'); do
for GRP in $(ls -l $SplitDir/*_ORF_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.aa";
InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp_neg.aa";
InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.tab";
InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.txt";
done
cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
done
done
```


#### E.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
	for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v '10300' | grep -w -e '404' -e '414' -e '415' -e '416'); do
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
  for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '10300' -e '414_v2' | grep -w -e '404' -e '414' -e '415' -e '416'); do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/combined_sigP_ORF/$Organism/$Strain
    mkdir -p $OutDir
    echo "The following number of sequences were predicted as secreted:"
    # cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
    cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa > $OutDir/"$Strain"_all_secreted.fa
    cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | tr -d '>' | tr -d ' ' | sed "s/HMM_score\t/HMM_score=\t/g" > $OutDir/"$Strain"_all_secreted_headers.txt
    cat $OutDir/"$Strain"_all_secreted_headers.txt | wc -l
    echo "This represented the following number of unique genes:"
    # cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    # $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
    cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
  done
```

```
  P.cactorum - 404
  The following number of sequences were predicted as secreted:
  55547
  This represented the following number of unique genes:
  26231
  P.cactorum - 414
  The following number of sequences were predicted as secreted:
  70912
  This represented the following number of unique genes:
  33601
  P.cactorum - 415
  The following number of sequences were predicted as secreted:
  60266
  This represented the following number of unique genes:
  28714
  P.cactorum - 416
  The following number of sequences were predicted as secreted:
  61712
  This represented the following number of unique genes:
  29425
  P.cactorum - 62471
  The following number of sequences were predicted as secreted:
  55686
  This represented the following number of unique genes:
  26324
  P.idaei - 371
  The following number of sequences were predicted as secreted:
  54493
  This represented the following number of unique genes:
  25750
  P.idaei - SCRP370
  The following number of sequences were predicted as secreted:
  54601
  This represented the following number of unique genes:
  25802
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

```bash
  for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff3 | grep -v -e 'atg' | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '10300' -e '414_v2' | grep -v -w -e '404' -e '414' -e '415' -e '416'); do
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

#<--- Progress up to here
The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '10300' -e '414_v2' | grep -v -w -e '404' -e '414' -e '415' -e '416'); do
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
mkdir -p $OutDir;
printf "\nstrain: $Strain\tspecies: $Organism\n";
printf "the number of SigP gene is:\t";
cat $Secretome | grep '>' | wc -l;
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
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt  $SigP_Gff	RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt $SigP_Merged_Gff RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex.gff
RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.gff
RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.txt
RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff --db sigP_ORF_RxLR.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR.db --id sigP_ORF_RxLR --out sigP_ORF_RxLR_merged.db --gff > $RxLR_Merged_Gff
cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
printf "Merged RxLR-EER regex proteins:\t"
cat $RxLR_Merged_AA | grep '>' | wc -l
printf "\n"
done
```

```
strain: 414	species: P.cactorum
the number of SigP gene is:	70912
the number of SigP-RxLR genes are:	1876
the number of SigP-RxLR-EER genes are:	220
Merged RxLR-EER regex proteins:197

strain: 404	species: P.cactorum
the number of SigP gene is:	55547
the number of SigP-RxLR genes are:	1519
the number of SigP-RxLR-EER genes are:	194
Merged RxLR-EER regex proteins:	170

strain: 415	species: P.cactorum
the number of SigP gene is:	60266
the number of SigP-RxLR genes are:	1591
the number of SigP-RxLR-EER genes are:	192
Merged RxLR-EER regex proteins:	168

strain: 416	species: P.cactorum
the number of SigP gene is:	61712
the number of SigP-RxLR genes are:	1621
the number of SigP-RxLR-EER genes are:	192
Merged RxLR-EER regex proteins:	168
```


### E5) From ORF gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search ORFs predicted with atg.pl. These were run with the following commands:


```bash
for Secretome in $(ls gene_pred/ORF_sigP/P.cactorum/10300/10300_ORF_sp_merged.aa); do
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
SigP_Merged_Gff=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_WY_hmmer.gff
done
```

<!-- P.cactorum 10300
Initial search space (Z):              15271  [actual number of targets]
Domain search space  (domZ):             113  [number of targets reported over threshold] -->

* P.cactorum 10300
* Initial search space (Z):              14767  [actual number of targets]
* Domain search space  (domZ):             113  [number of targets reported over threshold]


### E6) From ORF gene models - Hmm evidence of RxLR effectors

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep -w -e '414' | grep -v -w -e '404' -e '415' -e '416'); do
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
  P.cactorum 404
  Initial search space (Z):              55547  [actual number of targets]
  Domain search space  (domZ):             417  [number of targets reported over threshold]
  Merged RxLR-EER Hmm proteins:	129
  P.cactorum 414
  Initial search space (Z):              70912  [actual number of targets]
  Domain search space  (domZ):             495  [number of targets reported over threshold]
  Merged RxLR-EER Hmm proteins:	154
  P.cactorum 415
  Initial search space (Z):              60266  [actual number of targets]
  Domain search space  (domZ):             421  [number of targets reported over threshold]
  Merged RxLR-EER Hmm proteins:	129
  P.cactorum 416
  Initial search space (Z):              61712  [actual number of targets]
  Domain search space  (domZ):             412  [number of targets reported over threshold]
  Merged RxLR-EER Hmm proteins:	128
```

* P.cactorum 10300
* Initial search space (Z):              14767  [actual number of targets]
* Domain search space  (domZ):             144  [number of targets reported over threshold]

### E7) Combining RxLRs from Regex and hmm searches


The total RxLRs are

```bash
  for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_ORF_RxLR_EER_regex_merged.txt | grep -v -e 'Aug' -e '10300' | grep -e 'P.idaei' -e 'P.cactorum' | grep -e '414' -e '404' -e '415' -e '416'); do
    Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
    Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF.gff3)
    Proteome=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
    HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_ORF_RxLR_hmm_merged.txt)
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
    cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
    # cat $Gff | grep -w -f $OutDir/"$Strain"_total_ORF_RxLR_headers.txt > $OutDir/"$Strain"_total_ORF_RxLR.gff
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_total_ORF_RxLR_headers.txt $Gff ORF_RxLR Name Augustus > $OutDir/"$Strain"_total_ORF_RxLR.gff
    echo "Number of genes in the extracted gff file:"
    cat $OutDir/"$Strain"_total_ORF_RxLR.gff | grep -w 'gene' | wc -l
  done
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
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/*/* | grep -e 'P.idaei' -e 'P.cactorum' | grep -w -e '414' -e '404' -e '415' -e '416'); do
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

echo "$Species - $Strain"
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l
# cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
echo "The total number of putative RxLRs are:"
cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalRxLRsTxt
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
cat $ORFsUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f3 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
cat $TotalRxLRsTxt | wc -l
cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalRxLRsTxt > $TotalRxLRsGff

RxLRsFa=$MergeDir/"$Strain"_final_RxLR_EER.fa
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/unwrap_fasta.py --inp_fasta $AugFa | grep -A1 -w -f $AugTxt | grep -v -E '^--$' > $RxLRsFa
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalRxLRsTxt > $RxLRsFa
# echo "$Strain"
$ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalRxLRsTxt >> $RxLRsFa
# echo "$Strain"
echo "The number of sequences extracted is"
cat $RxLRsFa | grep '>' | wc -l
done
```

```
  P.cactorum - 404
  The number of ORF RxLRs overlapping Augustus RxLRs:
  129
  The number of Augustus RxLRs overlapping ORF RxLRs:
  129
  The number of RxLRs unique to ORF models:
  58
  The number of RxLRs unique to Augustus models:
  20
  The total number of putative RxLRs are:
  207
  P.cactorum - 414
  The number of ORF RxLRs overlapping Augustus RxLRs:
  153
  The number of Augustus RxLRs overlapping ORF RxLRs:
  153
  The number of RxLRs unique to ORF models:
  62
  The number of RxLRs unique to Augustus models:
  20
  The total number of putative RxLRs are:
  235
  P.cactorum - 415
  The number of ORF RxLRs overlapping Augustus RxLRs:
  134
  The number of Augustus RxLRs overlapping ORF RxLRs:
  134
  The number of RxLRs unique to ORF models:
  52
  The number of RxLRs unique to Augustus models:
  25
  The total number of putative RxLRs are:
  211
  P.cactorum - 416
  The number of ORF RxLRs overlapping Augustus RxLRs:
  130
  The number of Augustus RxLRs overlapping ORF RxLRs:
  130
  The number of RxLRs unique to ORF models:
  55
  The number of RxLRs unique to Augustus models:
  25
  The total number of putative RxLRs are:
  210
```



### H) From ORF gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in ORF gene models. This was done with the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/*/*/*.aa_cat.fa | grep -w -e 'P.cactorum' -e 'P.idaei' | grep -v -e 'atg' -e '10300' -e '414_v2' | grep -v -w -e '404' -e '414' -e '415' -e '416'); do
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
ORF_Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/*_ORF.gff3)
# Gff features were extracted for each header
CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $Headers $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
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
Searching for LFLAK domains in: P.cactorum 404
Initial search space (Z):             487049  [actual number of targets]
Domain search space  (domZ):             238  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 404
Initial search space (Z):             487049  [actual number of targets]
Domain search space  (domZ):             293  [number of targets reported over threshold]
The number of CRNs common to both models are:
139
Number of CRN ORFs after merging:
97
Searching for LFLAK domains in: P.cactorum 414
Initial search space (Z):             631759  [actual number of targets]
Domain search space  (domZ):             342  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 414
Initial search space (Z):             631759  [actual number of targets]
Domain search space  (domZ):             455  [number of targets reported over threshold]
The number of CRNs common to both models are:
223
Number of CRN ORFs after merging:
155
Searching for LFLAK domains in: P.cactorum 415
Initial search space (Z):             542852  [actual number of targets]
Domain search space  (domZ):             310  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 415
Initial search space (Z):             542852  [actual number of targets]
Domain search space  (domZ):             375  [number of targets reported over threshold]
The number of CRNs common to both models are:
178
Number of CRN ORFs after merging:
130
Searching for LFLAK domains in: P.cactorum 416
Initial search space (Z):             557317  [actual number of targets]
Domain search space  (domZ):             319  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 416
Initial search space (Z):             557317  [actual number of targets]
Domain search space  (domZ):             404  [number of targets reported over threshold]
The number of CRNs common to both models are:
200
Number of CRN ORFs after merging:
141
Searching for LFLAK domains in: P.cactorum 62471
Initial search space (Z):             488325  [actual number of targets]
Domain search space  (domZ):             245  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 62471
Initial search space (Z):             488325  [actual number of targets]
Domain search space  (domZ):             308  [number of targets reported over threshold]
The number of CRNs common to both models are:
153
Number of CRN ORFs after merging:
100
Searching for LFLAK domains in: P.idaei 371
Initial search space (Z):             487686  [actual number of targets]
Domain search space  (domZ):             217  [number of targets reported over threshold]
Searching for DWL domains in: P.idaei 371
Initial search space (Z):             487686  [actual number of targets]
Domain search space  (domZ):             254  [number of targets reported over threshold]
The number of CRNs common to both models are:
122
Number of CRN ORFs after merging:
86
Searching for LFLAK domains in: P.idaei SCRP370
Initial search space (Z):             486809  [actual number of targets]
Domain search space  (domZ):             220  [number of targets reported over threshold]
Searching for DWL domains in: P.idaei SCRP370
Initial search space (Z):             486809  [actual number of targets]
Domain search space  (domZ):             256  [number of targets reported over threshold]
The number of CRNs common to both models are:
125
Number of CRN ORFs after merging:
89
```



Extract crinklers from published gene models


```bash
  for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/*/* | grep -w -e 'P.cactorum' -e 'P.idaei' | grep -v -e 'atg' -e '10300' -e '414_v2' | grep -w -e '404' -e '414' -e '415' -e '416'); do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    AugGff=$(ls $MergeDir/"$Strain"_pub_CRN_LFLAK_DWL.gff)
    AugFa=$(ls gene_pred/codingquary/"$Species"/"$Strain"/final/final_genes_combined.pep.fasta)
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
    echo "$Species - $Strain"

    echo "The number of ORF CRNs overlapping Augustus CRNs:"
    cat $ORFsInAug | grep -w -e 'transcript' -e 'mRNA' | wc -l
    echo "The number of Augustus CRNs overlapping ORF CRNs:"
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA' | wc -l
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalCRNsTxt
    echo "The number of CRNs unique to ORF models:"
    cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' | wc -l
    cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt
    echo "The number of CRNs unique to Augustus models:"
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt

    cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalCRNsTxt > $TotalCRNsGff

    CRNsFa=$MergeDir/"$Strain"_final_CRN.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalCRNsTxt > $CRNsFa
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalCRNsTxt >> $CRNsFa
    echo "The number of sequences extracted is"
    cat $CRNsFa | grep '>' | wc -l

  done
```

```
P.cactorum - 404
The number of ORF CRNs overlapping Augustus CRNs:
93
The number of Augustus CRNs overlapping ORF CRNs:
93
The number of CRNs unique to ORF models:
4
The number of CRNs unique to Augustus models:
19
The number of sequences extracted is
113
P.cactorum - 414
The number of ORF CRNs overlapping Augustus CRNs:
151
The number of Augustus CRNs overlapping ORF CRNs:
151
The number of CRNs unique to ORF models:
4
The number of CRNs unique to Augustus models:
17
The number of sequences extracted is
173
P.cactorum - 415
The number of ORF CRNs overlapping Augustus CRNs:
126
The number of Augustus CRNs overlapping ORF CRNs:
126
The number of CRNs unique to ORF models:
4
The number of CRNs unique to Augustus models:
20
The number of sequences extracted is
147
P.cactorum - 416
The number of ORF CRNs overlapping Augustus CRNs:
137
The number of Augustus CRNs overlapping ORF CRNs:
137
The number of CRNs unique to ORF models:
4
The number of CRNs unique to Augustus models:
17
The number of sequences extracted is
155

```


## 4. 2 Ananlysis of RxLR effectors

Due to RxLR effectors being predicted from a number of sources the number of
unique RxLRs were identified from motif and Hmm searches within gene models.

Details on the commands run to identify this can be found within this repository
in 10300_analysis/effector_charactisation.md
 -->
