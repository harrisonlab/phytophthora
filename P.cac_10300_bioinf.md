#Phytophthora cactorum isolate 10300
==========

Scripts used for the analysis of the P. cactorum isolate 10300.
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/idris
. As such this script relies upon adhering to a 
specific directory structure.

The following is a summary of the work presented in this Readme.

The following processes were applied to Alternaria genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation

Analyses performed on these genomes involved BLAST searching for:
Known Avr genes
Ab initio prediction of putative RxLR, CRN and SSC effectors.
Upregulated genes during an RNAseq timecourse.


#Making local repositories for data

```shell
	Organism=P.cactorum
	Strain=10300
	ProjDir=/home/groups/harrisonlab/project_files/idris
	cd $ProjDir
	mkdir -p raw_dna/paired/$Organism/$Strain/F
	mkdir -p raw_dna/paired/$Organism/$Strain/R
	mkdir -p raw_dna/mate_paired/$Organism/$Strain/F
	mkdir -p raw_dna/mate_paired/$Organism/$Strain/R
```
Move raw data into local repositories:
```shell
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads
	cp $RawDatDir/060612_ID141_0001_650MVAAXX/Pcactorum_ID141_lane3_1Kb_R1.fastq raw_dna/paired/P.cactorum/10300/F/.
	cp $RawDatDir/060612_ID141_0001_650MVAAXX/Pcactorum_ID141_lane3_1Kb_R2.fastq raw_dna/paired/P.cactorum/10300/R/.
	cp $RawDatDir/060612_ID141_0001_650MVAAXX/Pcactorum_ID141_lane4_300bp_R1.fastq raw_dna/paired/P.cactorum/10300/F/.
	cp $RawDatDir/060612_ID141_0001_650MVAAXX/Pcactorum_ID141_lane4_300bp_R2.fastq raw_dna/paired/P.cactorum/10300/R/.
	cp $RawDatDir/170212_ID136_0001_64H5HAAXX/Pcactorum_ID136_lane4_300bp_R1.fastq raw_dna/paired/P.cactorum/10300/F/.
	cp $RawDatDir/170212_ID136_0001_64H5HAAXX/Pcactorum_ID136_lane4_300bp_R2.fastq raw_dna/paired/P.cactorum/10300/R/.
	cp $RawDatDir/170212_ID136_0001_64H5HAAXX/Pcactorum_ID136_lane5_1Kb_R1.fastq raw_dna/paired/P.cactorum/10300/F/.
	cp $RawDatDir/170212_ID136_0001_64H5HAAXX/Pcactorum_ID136_lane5_1Kb_R2.fastq raw_dna/paired/P.cactorum/10300/R/.
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/cactorum/Cactorum
	cp $RawDatDir/131008_matepair_Pc/Pcact10300_S2_L001_R1_001.fastq.gz raw_dna/mate-paired/P.cactorum/10300/F/.
	cp $RawDatDir/131008_matepair_Pc/Pcact10300_S2_L001_R2_001.fastq.gz raw_dna/mate-paired/P.cactorum/10300/R/.
```



#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```shell
	for RawData in $(ls raw_dna/paired/P.cactorum/10300/*/*.fastq*); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $RawData
	done
	for RawData in $(ls raw_dna/mate-paired/P.cactorum/10300/*/*001.fastq.gz); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```


Trimming was performed on data to trim adapters from 
sequences and remove poor quality data. This was done with fastq-mcf

```shell
	for ReadsF in $(ls raw_dna/paired/P.cactorum/10300/F/*.fastq); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
		ReadsR=$(echo $ReadsF | sed -E s%/F/%/R/%g | sed s/_R1/_R2/g)
		ls $ReadsF
		ls $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
	for StrainPath in $(ls -d raw_dna/mate-paired/P.cactorum/10300); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_TruSeq_mate_adapters.fa
		ReadsF=$(ls $StrainPath/F/*001.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*001.fastq.gz)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
```

Data quality was visualised once again following trimming:
```shell
	for RawData in $(ls qc_dna/*/P.cactorum/10300/*/*.fq.gz); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```


Long mate pair sequences were still in the reverse complement.
These were corrected using fastx toolkit's reverse complement
function.

The fastx toolkit is installed on the cluster, but it is the 0.0.12 version.
There are two problems with this install, firstly it is outdated and secondly
it is not installed on the /home/ directory so is not copied onto worked nodes
when using the SGE. The 2014 release 0.0.14 was installed locally using the
following commands:
```shell
	cd ~/prog
	wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2
	tar -jxvf  fastx_toolkit-0.0.14.tar.bz2
	cd fastx_toolkit-0.0.14
	make --prefix /home/armita/prog/fastx_toolkit-0.0.14 && make && 
	# Then the the bin directory was added to my profile
```


The following commands were used to reverse complement the mate-pair reads,
and then submit the reversed reads to fastqc for visualisation:
```shell
	cd /home/groups/harrisonlab/project_files/idris
	cat qc_dna/mate-paired/P.cactorum/10300/F/Pcact10300_S2_L001_R1_001_trim.fq.gz | gunzip -cf | fastx_reverse_complement -Q33 -o qc_dna/mate-paired/P.cactorum/10300/F/Pcact10300_S2_L001_R1_001_trim_rev.fq.gz
	cat qc_dna/mate-paired/P.cactorum/10300/R/Pcact10300_S2_L001_R2_001_trim.fq.gz | gunzip -cf | fastx_reverse_complement -Q33 -o qc_dna/mate-paired/P.cactorum/10300/R/Pcact10300_S2_L001_R2_001_trim_rev.fq.gz
	for RevRreads in $(ls qc_dna/mate-paired/P.cactorum/10300/*/*_trim_rev.fq.gz); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RevRreads; 
		qsub $ProgDir/run_fastqc.sh $RevRreads
	done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

```shell
	for TrimPath in $(ls -d qc_dna/paired/P.cactorum/10300); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF=$(ls $TrimPath/F/*.fq.gz)
		TrimR=$(ls $TrimPath/R/*.fq.gz)
		TrimMatePath=$(echo $TrimPath | sed 's/paired/mate-paired/g')
		TrimMateF=$(ls $TrimMatePath/F/*.fq.gz)
		TrimMateR=$(ls $TrimMatePath/R/*.fq.gz)
		echo $TrimF $TrimMateF
		echo $TrimR $TrimMateR
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimMateF $TrimR $TrimMateR
	done
```
<!-- 

#Assembly

Assembly was performed using Velvet

A range of hash lengths were used and the best assembly selected for subsequent analysis

Velvet had been previously compiled on a user-specific install with the follwing command:
```
	cd ~/prog/velvet_1.2.08
	make CATEGORIES=5 MAXKMERLENGTH=201 OPENMP=1
```

```shell
	for TrimPath in qc_dna/paired/P.cactorum/10300; do
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/velvet
	Strain=$(printf $TrimPath | rev | cut -f1 -d '/' | rev)

	MinHash=41
	MaxHash=81
	HashStep=2
	TrimF=$(ls $TrimPath/F/*.fastq*)
	TrimR=$(ls $TrimPath/R/*.fastq*)
	GenomeSz=70
	ExpCov=8
	MinCov=2
	InsLgth=600
	qsub $ProgDir/submit_velvet_range.sh $MinHash $MaxHash $HashStep \
	$TrimF $TrimR $GenomeSz $ExpCov $MinCov $InsLgth
	done
	
```


Assemblies were summarised to allow the best assembly to be determined by eye

```shell
	for StrainPath in $(ls -d assembly/velvet/P.cactorum/10300); do
		printf "N50\tMax_contig_size\tNumber of bases in contigs\tNumber of contigs\tNumber of contigs >=1kb\tNumber of contigs in N50\tNumber of bases in contigs >=1kb\tGC Content of contigs\n" > $StrainPath/P.cactorum_10300_assembly_stats.csv
		for StatsFile in $(ls $StrainPath/*/stats.txt); do 
			cat $StatsFile | rev | cut -f1 -d ' ' | rev | paste -d '\t' -s >> $StrainPath/P.cactorum_10300_assembly_stats.csv
		done
	done
	tail -n+1 assembly/velvet/P.cactorum/10300/P.cactorum_10300_assembly_stats.csv > assembly/velvet/P.cactorum_10300_assembly_stats.csv
```

 -->
