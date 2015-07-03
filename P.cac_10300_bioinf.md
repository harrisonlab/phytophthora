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

Unfortunately this assembly required more RAM than available on the 96Gb worker node.

As such Velvet can not be used until a node with more RAM become available. 

## Abyss assembly

###Trimming matepair junction adapters

Abyss is very sensitive to the presence of junction adapter sequences in
matepair datasets. To remove the occurence of any junction adapters the program
nextclip was used.
```shell
	qlogin
	cd /home/groups/harrisonlab/project_files/idris
	InF=qc_dna/mate-paired/P.cactorum/10300/F/Pcact10300_S2_L001_R1_001_trim.fq.gz
	InR=qc_dna/mate-paired/P.cactorum/10300/R/Pcact10300_S2_L001_R2_001_trim.fq.gz
	Organism=P.cactorum
	Strain=10300
	OutName="$Organism"_"$Strain"_trim_clip
	CurDir=$PWD
	WorkDir=/tmp/nextclip_no_rev
	mkdir -p $WorkDir
	cd $WorkDir
	cp $CurDir/$InF F_Read.fq.gz
	cp $CurDir/$InR R_Read.fq.gz
	gunzip *.gz
	nextclip -i F_Read.fq -j R_Read.fq -o "$OutName" > logfile.txt
	rm F_Read.fq
	rm R_Read.fq
	for File in (ls *.fq); do
		gzip *.fq
	done
	cp -r $WorkDir $CurDir/qc_dna/mate-paired/P.cactorum/10300/.
	rm -r $WorkDir
```


###Abyss assembly

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
	qlogin -pe smp 16 -l virtual_free=0.9G

	#---	Step 1		---
	# 		Set Variables
	#----------------------

	Strain=10300
	Organism=P.cactorum
	# KmerSz=31

	CurPath=/home/groups/harrisonlab/project_files/idris
	TrimPath=qc_dna/paired/P.cactorum/10300
	MatePath=qc_dna/mate-paired/P.cactorum/10300

	AssemblyName="$Strain"_abyss
	WorkDir=/tmp/"$Strain"_assembly
	OutDir=$CurPath/assembly/abyss/$Organism/$Strain

	Lib1F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane4_300bp_R1_trim.fq.gz
	Lib1R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane4_300bp_R2_trim.fq.gz

	Lib2F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane5_1Kb_R1_trim.fq.gz
	Lib2R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane5_1Kb_R2_trim.fq.gz

	Lib3F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane3_1Kb_R1_trim.fq.gz  
	Lib3R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane3_1Kb_R2_trim.fq.gz  

	Lib4F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane4_300bp_R1_trim.fq.gz
	Lib4R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane4_300bp_R2_trim.fq.gz

	Lib5F=$CurPath/$MatePath/nextclip/P.cactorum_10300_trim_rev_clip_D_R1.fastq.gz
	Lib5R=$CurPath/$MatePath/nextclip/P.cactorum_10300_trim_rev_clip_D_R2.fastq.gz


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

	# for KmerSz in 31 35 41 45 51 55 61 64; do
	for KmerSz in 60 61 62 63; do
		AssemblyDir="$AssemblyName"_"$KmerSz"
		mkdir -p $AssemblyDir
		echo "Running Abyss with kmer size:\t $KmerSz\n" 2>&1 |  tee -a $WorkDir/"$Organism"_"$Strain"_Abyss.log
		abyss-pe -C $AssemblyDir k=$KmerSz np=16 j=16 name=$AssemblyName lib='pe1 pe2 pe3 pe4' mp='mp5' pe1='../Lib1_1.fq.gz ../Lib1_2.fq.gz' pe2='../Lib2_1.fq.gz ../Lib2_2.fq.gz' pe3='../Lib3_1.fq.gz ../Lib3_2.fq.gz' pe4='../Lib4_1.fq.gz ../Lib4_2.fq.gz' mp5='../Lib5_1.fq.gz ../Lib5_2.fq.gz' 2>&1 |  tee -a $WorkDir/"$Organism"_"$Strain"_Abyss.log
	done

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

The results of abyss assemblies at different kmers were summarised using the commands:
```shell
	for AbyssDir in $(ls -d assembly/abyss/P.cactorum/10300/10300_abyss_*); do
		ls $AbyssDir/10300_abyss-stats.csv
		cat $AbyssDir/10300_abyss-stats.csv | tail -n +2 | sed 's/,/\t/g'
	done > assembly/abyss/P.cactorum/10300/assembly_stats.csv
```

The output files from abyss did not summarise the number of contaigs
that were assembled longer than 1000bp. These were determined using 
the program 
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
	for Unitigs in $(ls assembly/abyss/P.cactorum/10300/10300_abyss_*/10300_abyss-scaffolds.fa); do
		echo $Unitigs
		Contigs1Kb=$(echo $Unitigs | sed 's/.fa/_1000bp.fa/g')
		$ProgDir/filter_abyss_contigs.py $Unitigs 1000 > $Contigs1Kb
		printf "the number of assembled contigs are:\t"
		cat $Contigs1Kb | grep '>' | wc -l
		printf "The number of bases in these contigs (including N/n characters) are:\t"
		cat $Contigs1Kb | grep -v '>' | tr -d ' \t\n\r\f' | wc -c
		printf "The number of bases in these contigs (excluding N/n characters) are:\t"
		cat $Contigs1Kb | grep -v '>' | tr -d ' \t\n\r\f' | tr -d 'nN' | wc -c
	done
```

From this the kmer length of 51 was identified as the best assembly due to 
having a high number of bases in the assembly, with low contig numbers and 
showing the greatest N20-N80 values of the assemblies.

##Dip-spades assembly

 Mulitple assembly programs were trialed to identify which resulted in the best assembly syayistics. 
 The program dip-spades was  run. Spades is designed for bacterial and small genomes. Dipspades has 
 been designed to work on small diploid genomes.
 
 The commands used to run dip-spades were as follows:
 
 
```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/assembly
	qsub $ProgDir/dip-spades_10300_assembly.sh
```
 
<!-- 

 ```shell
	 screen -a
	 qlogin -pe smp 8 -l virtual_free=4G
 
 	#---	Step 1		---
	# 		Set Variables
	#----------------------

	Strain=10300
	Organism=P.cactorum
	# KmerSz=31

	CurPath=/home/groups/harrisonlab/project_files/idris
	TrimPath=qc_dna/paired/P.cactorum/10300
	MatePath=qc_dna/mate-paired/P.cactorum/10300

	AssemblyName="$Strain"_dip-spades
	WorkDir=/tmp/"$Strain"_dip-spades
	OutDir=$CurPath/assembly/dip-spades/$Organism/$Strain

	Lib1F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane4_300bp_R1_trim.fq.gz
	Lib1R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane4_300bp_R2_trim.fq.gz

	Lib2F=$CurPath/$TrimPath/F/Pcactorum_ID136_lane5_1Kb_R1_trim.fq.gz
	Lib2R=$CurPath/$TrimPath/R/Pcactorum_ID136_lane5_1Kb_R2_trim.fq.gz

	Lib3F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane3_1Kb_R1_trim.fq.gz  
	Lib3R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane3_1Kb_R2_trim.fq.gz  

	Lib4F=$CurPath/$TrimPath/F/Pcactorum_ID141_lane4_300bp_R1_trim.fq.gz
	Lib4R=$CurPath/$TrimPath/R/Pcactorum_ID141_lane4_300bp_R2_trim.fq.gz

	Lib5F=$CurPath/$MatePath/nextclip/P.cactorum_10300_trim_rev_clip_D_R1.fastq.gz
	Lib5R=$CurPath/$MatePath/nextclip/P.cactorum_10300_trim_rev_clip_D_R2.fastq.gz


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

	echo "Running dip-spades" 2>&1 |  tee -a $WorkDir/"$Organism"_"$Strain"_dipspades.log
	dipspades.py --pe1-1 Lib1_1.fq.gz --pe1-2 Lib1_2.fq.gz \
	--pe2-1 Lib2_1.fq.gz --pe2-2 Lib2_2.fq.gz \
	--pe3-1 Lib3_1.fq.gz --pe3-2 Lib3_2.fq.gz \
	--pe4-1 Lib4_1.fq.gz --pe4-2 Lib4_2.fq.gz \
	--mp1-1 Lib5_2.fq.gz --mp1-2 Lib5_1.fq.gz \
	-t 16 -m 96 -k 21,33,55,77 --careful -o $AssemblyName 2>&1 |  tee -a $WorkDir/"$Organism"_"$Strain"_dipspades.log

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
	cp -r $WorkDir/$AssemblyName/* $OutDir/.
	rm -r $WorkDir



	#---	Step 5		---
	# 		Exit
	#----------------------

	logout

	# Exit screen using crt+a ctrl+d

 ```

 -->

#Repeat masking

repeat masking was performed on the abyss (kmer51) assembly of 
the P. cactorum 10300 genome. The commands used were as follows:


Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking
	
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
	BestAss10300=assembly/abyss/P.cactorum/10300/10300_abyss_51/10300_abyss-scaffolds_1000bp.fa

	qsub $ProgDir/rep_modeling.sh $BestAss10300
	qsub $ProgDir/transposonPSI.sh $BestAss10300
```

##Parsing Assemblies
```
for Genome in $(ls repeat_masked/P.cactorum/10300/*/*_contigs_*.fa | grep -v 'parsed'); do 
ls $Genome; 
OutName=$(echo $Genome | sed 's/.fa/_parsed.fa/g')
cat $Genome | sed 's/>/>NODE_/g' | sed 's/ .*//' > $OutName
done
```


#Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified. 
	Gene model training
		- Gene models were trained for the 10300 repeatmasked geneome using assembled RNAseq data and predicted CEGMA genes.
	Gene prediction
		- Gene models were used to predict genes in the 103033 genome. This used RNAseq data as hints for gene models.

#Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

This was first performed on the 10300 unmasked assembly:

```shell
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
cd /home/groups/harrisonlab/project_files/idris/
for Genome in $(ls repeat_masked/P.cactorum/10300/*/*_contigs_unmasked_parsed.fa); do 
echo $Genome; 
qsub $ProgDir/sub_cegma.sh $Genome dna;
done
```

These results were then moved to a directory for unmasked genomes
```shell
	mv gene_pred/cegma/P.cactorum/10300 gene_pred/cegma/P.cactorum/10300_unmasked
```

The analysis was then repeated for the 10300 repeatmasked genome:
```shell
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
cd /home/groups/harrisonlab/project_files/idris/
for Genome in $(ls repeat_masked/P.cactorum/10300/*/*_contigs_hardmasked_parsed.fa); do 
echo $Genome; 
qsub $ProgDir/sub_cegma.sh $Genome dna;
done
```

These results were moved to a directory for repeatmasked data.
```shell
	mv gene_pred/cegma/P.cactorum/10300 gene_pred/cegma/P.cactorum/10300_hardmasked
```

Outputs were summarised using the commands:
```shell
	for File in $(ls gene_pred/cegma/P.cactorum/10300*/*_dna_cegma.completeness_report); do 
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev); 
		Species=$(echo $File | rev | cut -f3 -d '/' | rev); 
		printf "$Species\t$Strain\n"; 
		cat $File | head -n18 | tail -n+4;printf "\n"; 
	done > gene_pred/cegma/P.cactorum/10300_cegma_results_dna_summary.txt
	
	less gene_pred/cegma/P.cactorum/10300_cegma_results_dna_summary.txt
```








<!-- 
#Gene model training

Data quality was visualised using fastqc:
```shell
	for RawData in raw_rna/paired/*/*/*/*.fastq.gz; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from 
sequences and remove poor quality data. This was done with fastq-mcf

```shell
	for StrainPath in raw_rna/paired/*/*; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*.fastq.gz)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters RNA
	done
```

Data quality was visualised once again following trimming:
```shell
	for TrimData in qc_rna/paired/*/*/*/*.fastq.gz; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $TrimData
	done
```	

RNAseq data was assembled into transcriptomes using Trinity
```shell
	for StrainPath in qc_rna/paired/*/*; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/transcriptome_assembly
		ReadsF=$(ls $StrainPath/F/*.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*.fastq.gz)	
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/transcriptome_assembly_trinity.sh $ReadsF $ReadsR
	done
```	
Gene training was performed using RNAseq data. The cluster can not run this script using qlogin. As such it was run on the head node (-naughty) using screen.
Training for 650 and 1166 was performed in two instances of screen and occassionally viewed to check progress over time.
(screen is detached after opening using ctrl+a then ctrl+d. - if just ctrl+d is pressed the instance of screen is deleted. - be careful)
```shell
	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	Assembly650=assembly/trinity/A.alternata_ssp._gaisen/650/650_rna_contigs/Trinity.fasta 
	Genome650=repeat_masked/A.alternata_ssp._gaisen/650/A.alternata_ssp._gaisen_650_67_repmask/650_contigs_unmasked.fa
	$ProgDir/training_by_transcriptome.sh $Assembly650 $Genome650

	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	Assembly1166=assembly/trinity/A.alternata_ssp._tenuissima/1166/1166_rna_contigs/Trinity.fasta 
	Genome1166=repeat_masked/A.alternata_ssp._tenuissima/1166/A.alternata_ssp._tenuissima_1166_43_repmask/1166_contigs_unmasked.fa
	$ProgDir/training_by_transcriptome.sh $Assembly1166 $Genome1166
```

Quality of Trinity assemblies were assessed using Cegma to assess gene-space within the transcriptome
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	for Transcriptome in $(ls assembly/trinity/A.*/*/*_rna_contigs/Trinity.fasta); do  
		echo $Transcriptome;  
		qsub $ProgDir/sub_cegma.sh $Transcriptome rna; 
	done
```
Outputs were summarised using the commands:
```shell
	for File in $(ls gene_pred/cegma/A.alternata_ssp._*/*/*_rna_cegma.completeness_report); do 
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev); 
		Species=$(echo $File | rev | cut -f3 -d '/' | rev); 
		printf "$Species\t$Strain\n"; 
		cat $File | head -n18 | tail -n+4;printf "\n"; 
	done > gene_pred/cegma/cegma_results_rna_summary.txt
	
	less gene_pred/cegma/cegma_results_rna_summary.txt
```
	
 -->
 
 
 <!-- 
 
#Gene prediction

Gene prediction was performed for A. alternata isolates. 
RNAseq reads were used as Hints for the location of CDS. 
A concatenated dataset of both ssp. tenuissima and ssp. gaisen RNAseq reads were used as hints for all strains.
Genes were predicted for ssp. tenuissima using the gene model trained to ssp. tenunissima. 
Genes were predicted for ssp. gaisen using the gene model trained to ssp. gaisen. 
Genes were predicted for ssp. arborescens using the gene model trained to ssp. tenuisima.


```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	mkdir -p qc_rna/concatenated
	RnaFiles=$(ls qc_rna/paired/A.alternata_ssp._*/*/*/*.fastq.gz | paste -s -d ' ')
	ConcatRna=qc_rna/concatenated/A.alternata_RNA_650_1166.fa.gz
	cat $RnaFiles > $ConcatRna
	
	for Genome in repeat_masked/A.alternata_ssp._*/*/*/*_contigs_unmasked.fa; do
		Species=$(printf $Genome | rev | cut -f4 -d '/' | rev)
		if [ "$Species" == 'A.alternata_ssp._tenuissima' ]; then 
			GeneModel=A.alternata_ssp._tenuissima_1166
		elif [ "$Species" == 'A.alternata_ssp._gaisen' ]; then 
			GeneModel=A.alternata_ssp._gaisen_650
		elif [ "$Species" == 'A.alternata_ssp._arborescens' ]; then 
			GeneModel=A.alternata_ssp._tenuissima_1166
		fi
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRna $GeneModel
	done
```
 -->


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
