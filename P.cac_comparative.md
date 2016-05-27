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
```

#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:


```bash
  for RawData in $(ls raw_dna/paired/P.cactorum/*/*/*.fastq.gz | grep -v '10300'); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```


Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
  for Strain in 404 414 415 416 62471; do
    echo $Strain
    Read_F=$(ls raw_dna/paired/P.*/$Strain/F/*.fastq.gz)
    Read_R=$(ls raw_dna/paired/P.*/$Strain/R/*.fastq.gz)
    IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
    qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
  done
```

Trimming was then performed for strains with multiple runs of data

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
	IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
	echo "414"
	StrainPath=raw_dna/paired/P.cactorum/414
	ReadsF=$(ls $StrainPath/F/414_run1_F.fastq.gz)
	ReadsR=$(ls $StrainPath/R/414_run1_R.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	StrainPath=raw_dna/paired/P.cactorum/414
	ReadsF=$(ls $StrainPath/F/414_run2_F.fastq.gz)
	ReadsR=$(ls $StrainPath/R/414_run2_R.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
```

Data quality was visualised once again following trimming:

```bash
  for RawData in $(ls qc_dna/paired/P.cactorum/*/*/*q.gz | grep -w -e '414'); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

```bash
  for TrimPath in $(ls -d qc_dna/paired/P.cactorum/*); do
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


#Assembly

Assembly was performed with:
* Spades

## Spades Assembly

```bash
for StrainPath in $(ls -d qc_dna/paired/P.cactorum/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*.fq.gz)
R_Read=$(ls $StrainPath/R/*.fq.gz)
OutDir=assembly/spades/$Organism/$Strain
Jobs=$(qstat | grep 'submit_SPA' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
	sleep 5m
	printf "."
	Jobs=$(qstat | grep 'submit_SPA' | grep 'qw' | wc -l)
done		
printf "\n"
echo $F_Read
echo $R_Read
qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct 10
done
```


```bash
  for StrainPath in $(ls -d qc_dna/paired/P.*/414); do  
    echo $StrainPath
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $StrainPath/F/414_run1_F_trim.fq.gz);
    TrimR1_Read=$(ls $StrainPath/R/414_run1_R_trim.fq.gz);
    TrimF2_Read=$(ls $StrainPath/F/414_run2_F_trim.fq.gz);
    TrimR2_Read=$(ls $StrainPath/R/414_run2_R_trim.fq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    OutDir=assembly/spades/$Organism/$Strain
    qsub $ProgDir/subSpades_2lib_HiMem.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir correct 10
  done
```

<!-- ```bash
for Assembly in $(ls assembly/spades/P.cactorum/*/*/*500bp.fasta); do
Kmer=$(echo $Assembly | rev | cut -f2 -d '/' | rev);
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
OutDir=assembly/spades/$Organism/$Strain/filtered_contigs; qsub $ProgDir/sub_quast.sh $Assembly $OutDir;
done
``` -->
