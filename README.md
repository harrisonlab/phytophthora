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
  mkdir -p raw_dna/paired/P.cactorum/415/F
  mkdir -p raw_dna/paired/P.cactorum/415/R
  mkdir -p raw_dna/paired/P.cactorum/416/F
  mkdir -p raw_dna/paired/P.cactorum/416/R
  mkdir -p raw_dna/paired/P.fragariae/A4/F
  mkdir -p raw_dna/paired/P.fragariae/A4/R
  mkdir -p raw_dna/paired/P.fragariae/SCRP245_v2/F
  mkdir -p raw_dna/paired/P.fragariae/SCRP245_v2/R
  mkdir -p raw_dna/paired/P.fragariae/Bc23/F
  mkdir -p raw_dna/paired/P.fragariae/Bc23/R
  mkdir -p raw_dna/paired/P.fragariae/Nov5/F
  mkdir -p raw_dna/paired/P.fragariae/Nov5/R
  mkdir -p raw_dna/paired/P.fragariae/Nov77/F
  mkdir -p raw_dna/paired/P.fragariae/Nov77/R
  mkdir -p raw_dna/paired/P.fragariae/ONT3/F
  mkdir -p raw_dna/paired/P.fragariae/ONT3/R
  mkdir -p raw_dna/paired/P.cactorum/2003_3/F
  mkdir -p raw_dna/paired/P.cactorum/2003_3/R
  mkdir -p raw_dna/paired/P.idaei/SCRP370/F
  mkdir -p raw_dna/paired/P.idaei/SCRP370/R
  mkdir -p raw_dna/paired/P.rubi/SCRP249/F
  mkdir -p raw_dna/paired/P.rubi/SCRP249/R
  mkdir -p raw_dna/paired/P.rubi/SCRP324/F
  mkdir -p raw_dna/paired/P.rubi/SCRP324/R
  mkdir -p raw_dna/paired/P.rubi/SCRP333/F
  mkdir -p raw_dna/paired/P.rubi/SCRP333/R
  mkdir -p raw_dna/paired/P.cactorum/2003_4/F
  mkdir -p raw_dna/paired/P.cactorum/2003_4/R
  mkdir -p raw_dna/paired/P.cactorum/PC13_15/F
  mkdir -p raw_dna/paired/P.cactorum/PC13_15/R


  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/150716_M01678_0023_AB0YF
  cp $RawDat/Pcactorum415_S2_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/415/F/.
  cp $RawDat/Pcactorum415_S2_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/415/R/.
  cp $RawDat/Pcactorum416_S1_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/416/F/.
  cp $RawDat/Pcactorum416_S1_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/416/R/.
  cp $RawDat/PfragariaeA4_S3_L001_R1_001.fastq.gz raw_dna/paired/P.fragariae/A4/F/.
  cp $RawDat/PfragariaeA4_S3_L001_R2_001.fastq.gz raw_dna/paired/P.fragariae/A4/R/.
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/150925_M01678_0029_AC669
  cp $RawDat/Pfrag-SCRP245_S3_L001_R1_001.fastq.gz raw_dna/paired/P.fragariae/SCRP245_v2/F/.
  cp $RawDat/Pfrag-SCRP245_S3_L001_R2_001.fastq.gz raw_dna/paired/P.fragariae/SCRP245_v2/R/.
  cp $RawDat/Pfrag-Bc23_S2_L001_R1_001.fastq.gz raw_dna/paired/P.fragariae/Bc23/F/.
  cp $RawDat/Pfrag-Bc23_S2_L001_R2_001.fastq.gz raw_dna/paired/P.fragariae/Bc23/R/.
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/150918_M01678_0028_AC60K
  cp $RawDat/Pfrag-Nov-5_S2_L001_R1_001.fastq.gz raw_dna/paired/P.fragariae/Nov5/F/.
  cp $RawDat/Pfrag-Nov-5_S2_L001_R2_001.fastq.gz raw_dna/paired/P.fragariae/Nov5/R/.
  cp $RawDat/Pfrag-Nov-77_S3_L001_R1_001.fastq.gz raw_dna/paired/P.fragariae/Nov77/F/.
  cp $RawDat/Pfrag-Nov-77_S3_L001_R2_001.fastq.gz raw_dna/paired/P.fragariae/Nov77/R/.
  cp $RawDat/Pfrag-ONT-3_S1_L001_R1_001.fastq.gz raw_dna/paired/P.fragariae/ONT3/F/.
  cp $RawDat/Pfrag-ONT-3_S1_L001_R2_001.fastq.gz raw_dna/paired/P.fragariae/ONT3/R/.
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160404_M004465_0008-ALVUT
  cp $RawDat/Pc20033_S2_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/2003_3/F/.
  cp $RawDat/Pc20033_S2_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/2003_3/R/.
  cp $RawDat/SCRP370_S3_L001_R1_001.fastq.gz  raw_dna/paired/P.idaei/SCRP370/F/.
  cp $RawDat/SCRP370_S3_L001_R2_001.fastq.gz  raw_dna/paired/P.idaei/SCRP370/R/.
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160412_M04465_0010-AMLCU
  cp $RawDat/SCRP249_S1_L001_R1_001.fastq.gz raw_dna/paired/P.rubi/SCRP249/F/.
  cp $RawDat/SCRP249_S1_L001_R2_001.fastq.gz raw_dna/paired/P.rubi/SCRP249/R/.
  cp $RawDat/SCRP324_S2_L001_R1_001.fastq.gz raw_dna/paired/P.rubi/SCRP324/F/.
  cp $RawDat/SCRP324_S2_L001_R2_001.fastq.gz raw_dna/paired/P.rubi/SCRP324/R/.
  cp $RawDat/SCRP333_S3_L001_R1_001.fastq.gz raw_dna/paired/P.rubi/SCRP333/F/.
  cp $RawDat/SCRP333_S3_L001_R2_001.fastq.gz raw_dna/paired/P.rubi/SCRP333/R/.
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160415_M004465_00011-AMLCL
  cp $RawDat/20033_S4_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/2003_4/F/.
  cp $RawDat/20033_S4_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/2003_4/R/.
  cp $RawDat/PC1315_S5_L001_R1_001.fastq.gz raw_dna/paired/P.cactorum/PC13_15/F/.
  cp $RawDat/PC1315_S5_L001_R2_001.fastq.gz raw_dna/paired/P.cactorum/PC13_15/R/.
```


#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
# 415
for RawData in $(ls raw_dna/paired/P.cactorum/415/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
# 416
for RawData in $(ls raw_dna/paired/P.cactorum/416/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
# A4
for RawData in $(ls raw_dna/paired/P.fragariae/A4/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
# SCRP245_v2
for RawData in $(ls raw_dna/paired/P.fragariae/SCRP245_v2/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
# Bc23
for RawData in $(ls raw_dna/paired/P.fragariae/Bc23/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
# Nov5
for RawData in $(ls raw_dna/paired/P.fragariae/Nov5/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
# Nov77
for RawData in $(ls raw_dna/paired/P.fragariae/Nov77/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
# ONT3
for RawData in $(ls raw_dna/paired/P.fragariae/ONT3/*/*.fastq.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
```



Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf

```bash
for Strain in 414; do
echo $Strain
Read_F=$(ls raw_dna/paired/P.*/$Strain/F/cact414_130517_S2_L001_R1_001.fastq.gz)
Read_R=$(ls raw_dna/paired/P.*/$Strain/R/cact414_130517_S2_L001_R2_001.fastq.gz)
IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

```bash
  for Strain in 415 416 A4; do
    echo $Strain
    Read_F=$(ls raw_dna/paired/P.*/$Strain/F/*.fastq.gz)
    Read_R=$(ls raw_dna/paired/P.*/$Strain/R/*.fastq.gz)
    IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
    qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
  done
```

```bash
  for Strain in SCRP245_v2 Bc23; do
    echo $Strain
    Read_F=$(ls raw_dna/paired/P.*/$Strain/F/*.fastq.gz)
    Read_R=$(ls raw_dna/paired/P.*/$Strain/R/*.fastq.gz)
    IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
    qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
  done
```

```bash
  for Strain in Nov5 Nov77 ONT3; do
    echo $Strain
    Read_F=$(ls raw_dna/paired/P.*/$Strain/F/*.fastq.gz)
    Read_R=$(ls raw_dna/paired/P.*/$Strain/R/*.fastq.gz)
    IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
    qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
  done
```

Data quality was visualised once again following trimming:

```bash
  for RawData in $(ls qc_dna/paired/P.cactorum/415/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
  # 416
  for RawData in $(ls qc_dna/paired/P.cactorum/416/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
  # A4
  for RawData in $(ls qc_dna/paired/P.fragariae/A4/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
  # SCRP245_v2
  for RawData in $(ls qc_dna/paired/P.fragariae/SCRP245_v2/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
  for RawData in $(ls qc_dna/paired/P.fragariae/SCRP245_v2/*/SCRP245_v2_no_adapt.fq); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
  # Bc23
  for RawData in $(ls qc_dna/paired/P.fragariae/Bc23/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
  # Nov5
  for RawData in $(ls qc_dna/paired/P.fragariae/Nov5/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
  # Nov77
  for RawData in $(ls qc_dna/paired/P.fragariae/Nov77/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
  # ONT3
  for RawData in $(ls qc_dna/paired/P.fragariae/ONT3/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

``bash
  qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/qc_trimmomatic.sh raw_dna/paired/P.fragariae/SCRP245_v2/F/Pfrag-SCRP245_S3_L001_R1_001.fastq.gz raw_dna/paired/P.fragariae/SCRP245_v2/R/Pfrag-SCRP245_S3_L001_R2_001.fastq.gz /home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa qc_dna/paired/P.fragariae/SCRP245_v2
```
```bash
  for RawData in $(ls qc_dna/paired/P.fragariae/SCRP245_v2/*/*.fq); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```

<!-- ## Error correction

First run error correction. (This job is CPU intensive rather than RAM intensive
and will run on any node of the cluster).

```bash
  for Strain in 415 416 A4 SCRP245_v2 Bc23; do
    echo $Strain
    Trim_F=$(ls qc_dna/paired/P.*/$Strain/F/*.fq.gz)
    Trim_R=$(ls qc_dna/paired/P.*/$Strain/R/*.fq.gz)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    OutDir=$(dirname $Trim_F | sed 's/F/corrected/')
    echo $OutDir
    qsub $ProgDir/sub_spades_correction.sh $F_Read $R_Read $OutDir
  done
``` -->

kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
  # for Strain in 415 416 A4 SCRP245_v2 Bc23 Nov5 Nov77; do
  for Strain in ONT3; do
    echo $Strain
    Trim_F=$(ls qc_dna/paired/P.*/$Strain/F/*.fq.gz)
    Trim_R=$(ls qc_dna/paired/P.*/$Strain/R/*.fq.gz)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
  done
```

** Estimated Genome Size is:

** Esimated Coverage is:


#Assembly
Assembly was performed using: Velvet / Abyss / Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis


```bash
  for Strain in 415 416 A4 SCRP245_v2 Bc23 Nov5 Nov77 ONT3; do
    F_Read=$(ls qc_dna/paired/P.*/$Strain/F/*.fq.gz)
    R_Read=$(ls qc_dna/paired/P.*/$Strain/R/*.fq.gz)
    CovCutoff='10'
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Species=$(echo $F_Read | rev | cut -f4 -d '/' | rev)
    OutDir=assembly/spades/$Species/$Strain
    echo $Species
    echo $Strain
    qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct $CovCutoff
    # qsub $ProgDir/submit_dipSPAdes.sh $F_Read $R_Read $OutDir correct $CovCutoff
  done
```

Quast

```bash
  for Strain in $(ls -d assembly/spades/*/* | rev | cut -f1 -d'/' | rev); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    Assembly=$(ls -d assembly/spades/*/$Strain/filtered_contigs/contigs_min_500bp.fasta)
    Species=$(echo $Assembly | rev | cut -f4 -d'/' | rev)
    OutDir=$(ls -d assembly/spades/*/$Strain/filtered_contigs)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

Assemblies were summarised to allow the best assembly to be determined by eye.

** Assembly stats are:
  * Assembly size:
  * N50:153669
  * N80:
  * N20:
  * Longest contig:687738
  **

As SPADes was run with the option to autodetect a minimum coverage the assembly was assessed to identify the coverage of assembled contigs. This was done using the following command:
```bash
	BestAss=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp.fasta
	cat $BestAss | grep '>' | cut -f6 -d'_' | sort -n | cut -f1 -d '.' | sort -n | uniq -c | less
```
From this it was determined that SPades could not be trusted to set its own minimum threshold for coverage.
In future an option will be be used to set a coverage for spades.
In the meantime contigs with a coverage lower than 10 were filtered out using the following commands:

```bash
	Headers=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.txt
	cat $BestAss | grep '>' | grep -E -v 'cov_.\..*_' > $Headers
	FastaMinCov=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.fasta
	cat $BestAss | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -A1 -f $Headers | grep -v -E '^\-\-' > $FastaMinCov
```

```bash
	~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants/remove_contaminants.py --inp ../neonectria_ditissima/assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.fasta  --out assembly/spades/N.galligena/R0905_v2/filtered_contigs/contigs_min_500bp_10x_filtered_renamed.fasta  --coord_file editfile.tab
```
We run Quast again.
```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
	Assembly=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.fasta
	OutDir=assembly/spades/N.ditissima/R0905_v2/contigs_min_500bp_10x_headers
	qsub $ProgDir/sub_quast.sh $Assembly $OutDir
```



<!--
## 5. BLAST Searches

## 5.1 Identifying SIX genes

Protein sequence of previously characterised SIX genes used to BLAST against
assemlies.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Query=analysis/blast_homology/oomycete_avr_genes/appended_oomycete_avr_cds.fasta
	for Assembly in $(ls assembly/spades/*/*/*/contigs_min_500bp.fasta); do
		echo $Assembly
		qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
	done
```

## 5.2 Looking for consensus Phytophthora sequencing consortium genes.

The phytophthora sequencing consortium has predicted genes for P. cactorum using
maker. The sequencing consortium gene models were blasted against the 10300
genome to identify consensus between these gene models and Braker1 gene models.

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  Query=analysis/blast_homology/seq_consortium_genes/Pcact_combine.fasta
  for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
    echo $Assembly
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```

-->
