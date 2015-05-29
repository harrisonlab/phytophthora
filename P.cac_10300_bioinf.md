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
	qlogin
	cd /home/groups/harrisonlab/project_files/idris
	cat qc_dna/mate-paired/P.cactorum/10300/F/Pcact10300_S2_L001_R1_001_trim.fq.gz | gunzip -cf | fastx_reverse_complement -z -Q33 -o qc_dna/mate-paired/P.cactorum/10300/F/Pcact10300_S2_L001_R1_001_trim_rev.fq.gz
	cat qc_dna/mate-paired/P.cactorum/10300/R/Pcact10300_S2_L001_R2_001_trim.fq.gz | gunzip -cf | fastx_reverse_complement -z -Q33 -o qc_dna/mate-paired/P.cactorum/10300/R/Pcact10300_S2_L001_R2_001_trim_rev.fq.gz
	logout
	cd /home/groups/harrisonlab/project_files/idris
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

 

#Assembly


## Velvet assembly


Assembly was performed using Velvet

A range of hash lengths were used and the best assembly selected for subsequent analysis

Velvet is installed in the master files directory of the cluster. There it was configured
to accept sequence data from two genomic libraries and to run with a max kmer length of 151.
A local install of velvet in my user profile was recompiled to accept up to 5 genomic libraries
as inputs and to accept longer maximum kmer lengths. The commands used to recompile velvet were 
as follows:
```shell
	cd ~/prog/velvet_1.2.08
	make CATEGORIES=5 MAXKMERLENGTH=201 OPENMP=1 OMP_NUM_THREADS=15
```

This local install of velvet was used to assemble the 10300 genome.
The assembly job was submitted to the SGE via a dedicated script for this assembly.

The command used to submit this assembly to the SGE was:
```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/assembly
	qsub $ProgDir/velvet_10300_assembly.sh
```

The script can be viewed in the repository harrison_lab/phytophthora/assembly .
However, for reference, the variables set in this script were:

```shell
	CurPath=$PWD
	ProgDir=/home/armita/prog/velvet_1.2.08
	TrimPath=qc_dna/paired/P.cactorum/10300
	MatePath=qc_dna/mate-paired/P.cactorum/10300
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/velvet

	MinHash=41
	MaxHash=81
	HashStep=2
	GenomeSz=70
	ExpCov=70
	MinCov=20
	Lib1InsLgth=300
	Lib1F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane4_300bp_R1_trim.fq.gz
	Lib1R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane4_300bp_R2_trim.fq.gz
	Lib2InsLgth=1000
	Lib2F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane5_1Kb_R1_trim.fq.gz
	Lib2R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane5_1Kb_R2_trim.fq.gz
	Lib3InsLgth=1000
	Lib3F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane3_1Kb_R1_trim.fq.gz  
	Lib3R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane3_1Kb_R2_trim.fq.gz  
	Lib4InsLgth=300
	Lib4F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane4_300bp_R1_trim.fq.gz
	Lib4R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane4_300bp_R2_trim.fq.gz
	Lib5InsLgth=5000
	Lib5F=$CurPath/$MatePath/F/Pcact10300_S2_L001_R1_001_trim_rev.fq.gz
	Lib5R=$CurPath/$MatePath/R/Pcact10300_S2_L001_R2_001_trim_rev.fq.gz
```

Unfortunately this assembly required more RAm than available on the 96Gb worker node.

As such Velvet can not be used until a node with more RAM become available. 

## Abyss assembly

Assembly was performed using Abyss.

A SGE script was written to perform assembly using Abyss on the cluster.

```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/assembly
	qsub $ProgDir/abyss_10300_assembly.sh
```

However due to openmpi not being installed on the cluster this scrip didn't work.


Abyss did work when using qlogin to work on the worker nodes of the cluster.
Therefore the following commands were used to perform assembly while using 
qlogin within a 'screen' session.

These were the commands used:

```shell
	screen -a
	qlogin -pe smp 16 -l virtual_free=4G

	#---	Step 1		---
	# 		Set Variables
	#----------------------

	Strain=10300
	Organism=P.cactorum
	KmerSz=31

	CurPath=/home/groups/harrisonlab/project_files/idris
	TrimPath=qc_dna/paired/P.cactorum/10300
	MatePath=qc_dna/mate-paired/P.cactorum/10300

	AssemblyName="$Strain"_abyss_"$KmerSz"
	WorkDir=/tmp/"$Strain"_assembly
	OutDir=$CurPath/assembly/abyss/$Organism/$Strain/$AssemblyName


	<!-- 
	CurPath=/home/groups/harrisonlab/project_files/alternaria
	Lib1F=$CurPath/tmp_F_1.fq.gz
	Lib1R=$CurPath/tmp_R_1.fq.gz

	Lib2F=$CurPath/tmp_F_2.fq.gz
	Lib2R=$CurPath/tmp_R_2.fq.gz

	Lib3F=$CurPath/tmp_F_3.fq.gz  
	Lib3R=$CurPath/tmp_R_3.fq.gz

	Lib4F=$CurPath/tmp_F_4.fq.gz
	Lib4R=$CurPath/tmp_R_4.fq.gz

	Lib5F=$CurPath/tmp_F_5.fq.gz
	Lib5R=$CurPath/tmp_R_5.fq.gz
	 -->

	Lib1F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane4_300bp_R1_trim.fq.gz
	Lib1R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane4_300bp_R2_trim.fq.gz

	Lib2F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane5_1Kb_R1_trim.fq.gz
	Lib2R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane5_1Kb_R2_trim.fq.gz

	Lib3F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane3_1Kb_R1_trim.fq.gz  
	Lib3R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane3_1Kb_R2_trim.fq.gz  

	Lib4F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane4_300bp_R1_trim.fq.gz
	Lib4R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane4_300bp_R2_trim.fq.gz

	Lib5F=$CurPath/$MatePath/F/Pcact10300_S2_L001_R1_001_trim_rev.fq.gz
	Lib5R=$CurPath/$MatePath/R/Pcact10300_S2_L001_R2_001_trim_rev.fq.gz


	#---	Step 2		---
	# 		Copy data onto
	#		Worker node
	#----------------------

	mkdir -p $WorkDir
	cd $WorkDir

	cp $Lib1F Lib1_1.fq.gz
	cp $Lib1R Lib1_2.fq.gz

	cp $Lib2F Lib2_1.fq.gz
	cp $Lib2R Lib2_2.fq.gz

	cp $Lib3F Lib3_1.fq.gz
	cp $Lib3R Lib3_2.fq.gz

	cp $Lib4F Lib4_1.fq.gz
	cp $Lib4R Lib4_2.fq.gz

	cp $Lib5F Lib5_1.fq.gz
	cp $Lib5R Lib5_2.fq.gz


	#---	Step 3		---
	# 		Assemble
	#----------------------


	abyss-pe k=$KmerSz np=16 j=16 name=$AssemblyName lib='pe1 pe2 pe3 pe4 mp5' pe1='Lib1_1.fq.gz Lib1_2.fq.gz' pe2='Lib2_1.fq.gz Lib2_2.fq.gz' pe3='Lib3_1.fq.gz Lib3_2.fq.gz' pe4='Lib4_1.fq.gz Lib4_2.fq.gz' mp5='Lib5_1.fq.gz Lib5_2.fq.gz'


	#---	Step 4		---
	# 		Cleanup
	#----------------------

	rm Lib1_1.fq.gz
	rm Lib1_2.fq.gz

	rm Lib2_1.fq.gz
	rm Lib2_2.fq.gz

	rm Lib3_1.fq.gz
	rm Lib3_2.fq.gz

	rm Lib4_1.fq.gz
	rm Lib4_2.fq.gz

	rm Lib5_1.fq.gz
	rm Lib5_2.fq.gz

	mkdir -p $OutDir
	cp -r $WorkDir/* $OutDir/.
	rm -r $WorkDir



	#---	Step 5		---
	# 		Exit
	#----------------------

	logout

	# Exit screen using crt+a ctrl+d

``` 

<!--

The script 

```shell
	
	ProgArgs="$MinHash $MaxHash $HashStep $GenomeSz $ExpCov $MinCov " \
	. "$Lib1InsLgth $Lib1F $Lib1R " . "$Lib2InsLgth $Lib2F $Lib2R " \
	. "$Lib3InsLgth $Lib3F $Lib3R " . "$Lib4InsLgth $Lib4F $Lib4R " \
	. "$Lib5InsLgth $Lib5F $Lib5R"	
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
