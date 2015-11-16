#Phytophthora cactorum isolate 10300
==========

Scripts used for the analysis of the P. cactorum isolate 10300.
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/idris
. As such this script relies upon adhering to a
specific directory structure.

The following is a summary of the work presented in this Readme.

The following processes were applied to the P. cactorum 10300 genome prior to analysis:
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
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
		ReadsR=$(echo $ReadsF | sed -E s%/F/%/R/%g | sed s/_R1/_R2/g)
		ls $ReadsF
		ls $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
	for StrainPath in $(ls -d raw_dna/mate-paired/P.cactorum/10300); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
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
<!--
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
-->

## Abyss assembly

###Trimming matepair junction adapters

Abyss is very sensitive to the presence of junction adapter sequences in
matepair datasets. To remove the occurence of any junction adapters the program
nextclip was used.

Cateory A reads contain adapters in both reads of a pair. Category B and C reads
contain an adapter in a single read of the pair. Category D reads don't contain
adapters in the either read of the pair.
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
	for File in $(ls *.fastq); do
		gzip $File
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
	qlogin -l h=blacklace11 -l virtual_free=12G -pe smp 8

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

	Lib5F=$CurPath/$MatePath/nextclip_no_rev/P.cactorum_10300_trim_clip_D_R1.fastq.gz
	Lib5R=$CurPath/$MatePath/nextclip_no_rev/P.cactorum_10300_trim_clip_D_R2.fastq.gz


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

	for KmerSz in 41 47 49 51 53 55 61; do
	# for KmerSz in 60 61 62 63; do
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

The output files from abyss did not summarise the number of contigs
that were assembled longer than 500bp. These were determined using
the program
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
		for Unitigs in $(ls assembly/abyss/P.cactorum/10300/10300_abyss_*/10300_abyss-scaffolds.fa); do
		echo $Unitigs
		Contigs500=$(echo $Unitigs | sed 's/.fa/_500bp.fa/g')
		$ProgDir/filter_abyss_contigs.py $Unitigs 500 > $Contigs500
		printf "the number of assembled contigs are:\t"
		cat $Contigs500 | grep '>' | wc -l
		printf "The number of bases in these contigs (including N/n characters) are:\t"
		cat $Contigs500 | grep -v '>' | tr -d ' \t\n\r\f' | wc -c
		printf "The number of bases in these contigs (excluding N/n characters) are:\t"
		cat $Contigs500 | grep -v '>' | tr -d ' \t\n\r\f' | tr -d 'nN' | wc -c
	done
```

From this the kmer length of 51 was identified as the best assembly due to
having a high number of bases in the assembly, with low contig numbers and
showing the greatest N20-N80 values of the assemblies.


### Run Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    for Assembly in $(ls assembly/abyss/P.cactorum/10300/10300_abyss_*/10300_abyss-scaffolds_500bp.fa); do
		Kmer=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/abyss/$Organism/$Strain/$Kmer
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

### Renaming contigs

Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/abyss/P.cactorum/10300/10300_abyss_53/10300_abyss-scaffolds_500bp.fa); do
		Kmer=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/abyss/$Organism/$Strain/$Kmer
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/10300_abyss-scaffolds_500bp_renamed.fa --coord_file tmp.csv
  done
  rm tmp.csv
```

#Repeat masking

repeat masking was performed on the abyss (kmer51) assembly of
the P. cactorum 10300 genome. The commands used were as follows:


Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assembly was used to perform repeatmasking

```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
	BestAss10300=assembly/abyss/P.cactorum/10300/10300_abyss_53/10300_abyss-scaffolds_500bp_renamed.fa

	qsub $ProgDir/rep_modeling.sh $BestAss10300
	qsub $ProgDir/transposonPSI.sh $BestAss10300
```


# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained for the 10300 repeatmasked geneome using assembled RNAseq data and predicted CEGMA genes.
	Gene prediction
		- Gene models were used to predict genes in the 103033 genome. This used RNAseq data as hints for gene models.

## Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

This was first performed on the 10300 unmasked assembly:

```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	Genome=repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
	qsub $ProgDir/sub_cegma.sh $Genome dna
```

## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.
 * The commands used to do this can be found in /gene_prediction/10300_braker1_prediction.md

## Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing
Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	Genome=repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
	qsub $ProgDir/run_ORF_finder.sh $Genome
	# ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen
	# qsub $ProgDir/path_pipe.sh $Genome
```

The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
	ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
	ORF_Gff=gene_pred/ORF_finder/P.cactorum/10300/10300_ORF.gff
	ORF_Gff_mod=gene_pred/ORF_finder/P.cactorum/10300/10300_ORF_corrected.gff3
	$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
```

<!-- ```bash
	for Genome in $(ls repeat_masked/P.cactorum/10300/10300_abyss_51_repmask/10300_contigs_unmasked_parsed.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
		ConcatRna=qc_rna/paired/genbank/P.cactorum/10300_genbank_appended.fastq
		Organism=P.cactorum
		Strain=10300
		echo $Genome
		GeneModel=P.cactorum_10300
		OutDir=gene_pred/augustus/Model-"$GeneModel"/$Organism/$Strain
		qsub $ProgDir/submit_augustus.sh $GeneModel $Genome false $OutDir
		GeneModel=fusarium
		OutDir=gene_pred/augustus/Model-"$GeneModel"/$Organism/$Strain
		qsub $ProgDir/submit_augustus.sh $GeneModel $Genome false $OutDir
		GeneModel=P.cactorum_10300
		OutDir=gene_pred/augustus/Model-"$GeneModel"_hints/$Organism/$Strain
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRna $GeneModel $OutDir
	done
``` -->


#Functional annotation

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	cd /home/groups/harrisonlab/project_files/idris
	getAnnoFasta.pl gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.gff
	Proteome=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	ProgDir/sub_interproscan.sh $Proteome
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	PredGenes=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
	InterProRaw=gene_pred/interproscan/10300/P.cactorum/raw
	OutDir=gene_pred/interproscan/P.cactorum/10300
	Organism=P.cactorum
	Strain=10300
	mkdir -p $OutDir
	printf "" > $OutDir/"$Strain"_interproscan.tsv
	printf "" > $OutDir/"$Strain"_interproscan.xml
	printf "" > $OutDir/interpro_features.gff
	printf "" > $OutDir/"$Strain"_interpro.gff3
	for File in $(ls -v $InterProRaw/*_split_*.tsv); do
		cat $File >> $OutDir/"$Strain"_interproscan.tsv
	done
	for File in $(ls -v $InterProRaw/*_split_*.xml); do
		cat $File >> $OutDir/"$Strain"_interproscan.xml
	done
	for File in $(ls -v $InterProRaw/*_split_*.gff3); do
		FastaStart=$(cat $File | grep -E "^##FASTA" -n | cut -d ':' -f1)
		cat $File | head -n "$FastaStart" | grep -v -E "^#" >> $OutDir/interpro_features.gff
	done
	cat $OutDir/interpro_features.gff $PredGenes >> $OutDir/"$Strain"_interpro.gff3
	rm $OutDir/interpro_features.gff
```
#Genomic analysis

## RxLR genes

Putative RxLR genes were identified within Augustus gene models using a number
of approaches:

 * A) From Augustus gene models - Signal peptide & RxLR motif  
 * B) From Augustus gene models - Hmm evidence of WY domains  
 * C) From Augustus gene models - Hmm evidence of RxLR effectors
 * D) From Augustus gene models - Hmm evidence of CRN effectors  
 * E) From ORF fragments - Signal peptide & RxLR motif  
 * F) From ORF fragments - Hmm evidence of WY domains  
 * G) From ORF fragments - Hmm evidence of RxLR effectors


### A) From Augustus gene models - Signal peptide & RxLR motif

Required programs:
 * SigP
 * biopython


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
	# getAnnoFasta.pl gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.gff
	Proteome=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
	Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
	SplitDir=gene_pred/braker_split/$Organism/$Strain
	mkdir -p $SplitDir
	BaseName="$Organism""_$Strain"_braker_preds
	$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
	for File in $(ls $SplitDir/*_braker_preds_*); do
		Jobs=$(qstat | grep 'pred_sigP' | wc -l)
		while [ $Jobs -ge 32 ]; do
			sleep 10
			printf "."
			Jobs=$(qstat | grep 'pred_sigP' | wc -l)
		done
		printf "\n"
		echo $File
		qsub $ProgDir/pred_sigP.sh $File
	done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
	SplitDir=gene_pred/braker_split/P.cactorum/10300
	Strain=$(echo $SplitDir | cut -d '/' -f4)
	Organism=$(echo $SplitDir | cut -d '/' -f3)
	InStringAA=''
	InStringNeg=''
	InStringTab=''
	InStringTxt=''
	for GRP in $(ls -l $SplitDir/*_braker_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
	InStringAA="$InStringAA gene_pred/braker_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.aa";  
	InStringNeg="$InStringNeg gene_pred/braker_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp_neg.aa";  
	InStringTab="$InStringTab gene_pred/braker_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.tab";
	InStringTxt="$InStringTxt gene_pred/braker_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.txt";  
	done
	cat $InStringAA > gene_pred/braker_sigP/$Organism/$Strain/"$Strain"_aug_sp.aa
	cat $InStringNeg > gene_pred/braker_sigP/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
	tail -n +2 -q $InStringTab > gene_pred/braker_sigP/$Organism/$Strain/"$Strain"_aug_sp.tab
	cat $InStringTxt > gene_pred/braker_sigP/$Organism/$Strain/"$Strain"_aug_sp.txt
	Headers=gene_pred/braker_sigP/$Organism/$Strain/"$Strain"_aug_sp_headers.txt
	# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	# BrakerGff=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.gff
	# ExtractedGff=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus_extracted.gff
	# cat $BrakerGff | grep -v '#' > $ExtractedGff
	# SigPGff=gene_pred/braker_sigP/$Organism/$Strain/"$Strain"_aug_sp.gff
	# cat gene_pred/braker_sigP/$Organism/$Strain/"$Strain"_aug_sp.aa | grep '>' | tr -d '>' | cut -f1 -d ' ' > $Headers
	# $ProgDir/gene_list_to_gff.pl $Headers $ExtractedGff SigP Name Augustus > $SigPGff
```
The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
	ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors;
	for Secretome in $(ls gene_pred/braker_sigP/P.cactorum/10300/10300_aug_sp.aa); do
		Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
		Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
		OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
		mkdir -p $OutDir;
		printf "\nstrain: $Strain\tspecies: $Organism\n";
		printf "the number of SigP gene is:\t";
		cat $Secretome | grep '>' | wc -l;
		printf "the number of SigP-RxLR genes are:\t";
		$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa;
		cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' > $OutDir/"$Strain"_Aug_RxLR_regex.txt
		cat $OutDir/"$Strain"_Aug_RxLR_regex.txt | wc -l
		printf "the number of SigP-RxLR-EER genes are:\t";
		cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' > $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt
		cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt | wc -l
		printf "\n"
		GeneModels=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.gff
		cat $GeneModels | grep -w -f $OutDir/"$Strain"_Aug_RxLR_regex.txt > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
		cat $GeneModels | grep -w -f $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt > $OutDir/"$Strain"_Aug_RxLR_EER_regex.gff3
	done
```

strain: 10300	species: P.cactorum
the number of SigP gene is:	2204
the number of SigP-RxLR genes are:	207
the number of SigP-RxLR-EER genes are:	109

### B) From Augustus gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search gene models predicted with Braker1. These were run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_Aug_WY_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_Aug_WY_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
		Headers="$Strain"_Aug_WY_hmmer_headers.txt
		cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' > $OutDir/$Headers
		GeneModels=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.gff
		cat $GeneModels | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_WY_hmmer.gff3
	done
```

P.cactorum 10300
Initial search space (Z):              20689  [actual number of targets]
Domain search space  (domZ):             171  [number of targets reported over threshold]

### C) From Augustus gene models - Hmm evidence of RxLR effectors
```bash
	for Proteome in $(ls gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa); do
		ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
		HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_Aug_RxLR_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_Aug_RxLR_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
		Headers="$Strain"_Aug_RxLR_hmmer_headers.txt
		cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' > $OutDir/$Headers
		# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		# Col2=cropped.hmm
		GeneModels=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.gff
		# $ProgDir/gene_list_to_gff.pl $OutDir/$Headers $GeneModels $Col2 Name > $OutDir/"$Strain"_Aug_RxLR_hmmer.gff3
		cat $GeneModels | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_RxLR_hmmer.gff3
	done
```

P.cactorum 10300
Initial search space (Z):              20689  [actual number of targets]
Domain search space  (domZ):             124  [number of targets reported over threshold]

### D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
for Proteome in $(ls gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_Aug_CRN_hmmer.txt
hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_aug_CRN_hmmer_out.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
Headers="$Strain"_Aug_RxLR_hmmer_headers.txt
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' > $OutDir/$Headers
GeneModels=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.gff
cat $GeneModels | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_CRN_hmmer.gff3
done
```

P.cactorum 10300
Initial search space (Z):              20689  [actual number of targets]
Domain search space  (domZ):              92  [number of targets reported over threshold]


### E) From ORF gene models - Signal peptide & RxLR motif

Required programs:
 * SigP
 * biopython


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
	Proteome=gene_pred/ORF_finder/P.cactorum/10300/10300.aa_cat.fa
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
	Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
	SplitDir=gene_pred/ORF_split/$Organism/$Strain
	mkdir -p $SplitDir
	BaseName="$Organism""_$Strain"_ORF_preds
	$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
	for File in $(ls $SplitDir/*_ORF_preds_*); do
		Jobs=$(qstat | grep 'pred_sigP' | wc -l)
		while [ $Jobs -ge 32 ]; do
			sleep 10
			printf "."
			Jobs=$(qstat | grep 'pred_sigP' | wc -l)
		done
		printf "\n"
		echo $File
		qsub $ProgDir/pred_sigP.sh $File
	done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
	SplitDir=gene_pred/ORF_split/P.cactorum/10300
	Strain=$(echo $SplitDir | cut -d '/' -f4)
	Organism=$(echo $SplitDir | cut -d '/' -f3)
	InStringAA=''
	InStringNeg=''
	InStringTab=''
	InStringTxt=''
	for GRP in $(ls -l $SplitDir/*_ORF_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
		InStringAA="$InStringAA gene_pred/ORF_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.aa";  
		InStringNeg="$InStringNeg gene_pred/ORF_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp_neg.aa";  
		InStringTab="$InStringTab gene_pred/ORF_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.tab";
		InStringTxt="$InStringTxt gene_pred/ORF_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.txt";  
	done
	cat $InStringAA > gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp.aa
	cat $InStringNeg > gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_neg_sp.aa
	tail -n +2 -q $InStringTab > gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp.tab
	cat $InStringTxt > gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp.txt
```

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
overlapps and identify the ORF with the best signalP score.

```bash
	Organism=P.cactorum
	Strain=10300
	ORF_Gff=gene_pred/ORF_finder/P.cactorum/10300/10300_ORF_corrected.gff3

	SigP_fasta=gene_pred/ORF_sigP/P.cactorum/10300/10300_ORF_sp.aa
	SigP_headers=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_names.txt

	SigP_Gff=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_unmerged.gff
	SigP_Merged_Gff=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.gff
	SigP_Merged_txt=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.txt
	SigP_Merged_AA=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.aa

	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	$ProgDir/extract_gff_for_sigP_hits.pl $SigP_headers $ORF_Gff SigP Name >$SigP_Gff
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	$ProgDir/make_gff_database.py --inp $SigP_Gff --db sigP_ORF.db
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF.db --id sigP_ORF --out sigP_ORF_merged.db --gff > $SigP_Merged_Gff
	cat $SigP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d'=' | rev > $SigP_Merged_txt
	$ProgDir/extract_from_fasta.py --fasta $SigP_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
```

The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
	for Secretome in $(ls gene_pred/ORF_sigP/P.cactorum/10300/10300_ORF_sp_merged.aa); do
		ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
		Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
		Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
		OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
		mkdir -p $OutDir;
		printf "\nstrain: $Strain\tspecies: $Organism\n";
		printf "the number of SigP gene is:\t";
		cat $Secretome | grep '>' | wc -l;
		printf "the number of SigP-RxLR genes are:\t";
		$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa;
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/"$Strain"_ORF_RxLR_regex.txt
		cat $OutDir/"$Strain"_ORF_RxLR_regex.txt | wc -l
		printf "the number of SigP-RxLR-EER genes are:\t";
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' '> $OutDir/"$Strain"_ORF_RxLR_EER_regex.txt
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.txt | wc -l
		printf "\n"
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_regex.txt  $SigP_Merged_Gff RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex.gff
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex.txt $SigP_Merged_Gff RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex.gff
	done
```

* strain: 10300	species: P.cactorum
* the number of SigP gene is: 15271
* the number of SigP-RxLR genes are: 935
* the number of SigP-RxLR-EER genes are: 170

### F) From ORF gene models - Hmm evidence of WY domains
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

P.cactorum 10300
Initial search space (Z):              15271  [actual number of targets]
Domain search space  (domZ):             113  [number of targets reported over threshold]

<!--
```bash
	for Proteome in $(ls gene_pred/ORF_finder/P.cactorum/10300/10300.aa_cat.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
		HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_WY_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_WY_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
		Headers="$Strain"_ORF_WY_hmmer_headers.txt
		cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
		SigP_Merged_Gff=gene_pred/ORF_finder/$Organism/$Strain/10300_ORF_corrected.gff3
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_WY_hmmer.gff
	done
```

P.cactorum 10300
Initial search space (Z):             459307  [actual number of targets]
Domain search space  (domZ):            1124  [number of targets reported over threshold] -->

### G) From ORF gene models - Hmm evidence of RxLR effectors

```bash
	for Secretome in $(ls gene_pred/ORF_sigP/P.cactorum/10300/10300_ORF_sp_merged.aa); do
		ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
		HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
		Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_RxLR_hmmer.txt
		hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
		Headers="$Strain"_ORF_RxLR_hmmer_headers.txt
		cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
		SigP_Merged_Gff=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.gff
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_RxLR_hmmer.gff3
	done
```

* P.cactorum 10300
* Initial search space (Z):              15271  [actual number of targets]
* Domain search space  (domZ):             145  [number of targets reported over threshold]

<!--
```bash
	for Proteome in $(ls gene_pred/ORF_finder/P.cactorum/10300/10300.aa_cat.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
		HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_RxLR_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
		Headers="$Strain"_ORF_RxLR_hmmer_headers.txt
		cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
		SigP_Merged_Gff=gene_pred/ORF_finder/$Organism/$Strain/10300_ORF_corrected.gff3
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_RxLR_hmmer.gff3
	done
```

P.cactorum 10300
Initial search space (Z):             459307  [actual number of targets]
Domain search space  (domZ):             365  [number of targets reported over threshold] -->


### H) From ORF gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in ORF gene models. This was done with the following commands:

```bash
	for Secretome in $(ls gene_pred/ORF_sigP/P.cactorum/10300/10300_ORF_sp_merged.aa); do
		ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
		HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
		Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_CRN_hmmer.txt
		hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_CRN_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
		Headers="$Strain"_CRN_hmmer_headers.txt
		cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
		SigP_Merged_Gff=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.gff
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_CRN_hmmer.gff3
	done
```
P.cactorum 10300
Initial search space (Z):              15271  [actual number of targets]
Domain search space  (domZ):              18  [number of targets reported over threshold]

<!--
```bash
	for Proteome in $(ls gene_pred/ORF_finder/P.cactorum/10300/10300.aa_cat.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
		HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_CRN_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_CRN_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
		Headers="$Strain"_CRN_hmmer_headers.txt
		cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
		SigP_Merged_Gff=gene_pred/ORF_finder/$Organism/$Strain/10300_ORF_corrected.gff3
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_CRN_hmmer.gff3
	done
```

P.cactorum 10300
Initial search space (Z):             459307  [actual number of targets]
Domain search space  (domZ):             225  [number of targets reported over threshold] -->

## 5. BLAST Searches

## 5.1 Identifying RxLR homolgs

nucleotide sequence of previously characterised RxLR genes were used to perform
BLAST searches against assemlies.

A query file was built from the RxLRs identified from transcriptome sequencing
in Chen et al 2014.
```bash
	cat analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.csv | grep -v 'Additional' | grep -v 'Index' | grep -v -e '^,' | cut -f 2,5 -d ',' | sed -r 's/^/>/g' | sed 's/,/\n/g' > analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.fa
```

### 5.1.a Blasting for Chen et al P. cactorum RxLRs:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Assembly=repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
	Query=analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.fa
	qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
	Query=analysis/blast_homology/oomycete_avr_genes/appended_oomycete_avr_cds.fasta
	qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	BlastHits=analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_avr_cds.fasta_homologs.csv
	HitsGff=analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_avr_cds.fasta_homologs.gff
	Column2=Avr_homolog
	NumHits=1
	$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
	BlastHits=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.csv
	HitsGff=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.gff
	Column2=Chen_RxLR
	NumHits=1
	$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
```

### 5.1.a Blasting for Hmmer predicted RxLRs from publihsed genomes.

These RxLRs were predicted from RxLR Hmm motifs.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Assembly=repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
	PinfRxLR=analysis/RxLR_effectors/hmmer_RxLR/P.infestans/T30-4/T30-4_pub_RxLR_hmmer.fa
	qsub $ProgDir/blast_pipe.sh $PinfRxLR dna $Assembly
	PparRxLR=analysis/RxLR_effectors/hmmer_RxLR/P.parisitica/310/310_pub_RxLR_hmmer.fa
	qsub $ProgDir/blast_pipe.sh $PinfRxLR dna $Assembly
	PcapRxLR=analysis/RxLR_effectors/hmmer_RxLR/P.capsici/LT1534/LT1534_pub_RxLR_hmmer.fa
	qsub $ProgDir/blast_pipe.sh $PinfRxLR dna $Assembly
	PsojRxLR=analysis/RxLR_effectors/hmmer_RxLR/P.sojae/67593/67593_pub_RxLR_hmmer.fa
	qsub $ProgDir/blast_pipe.sh $PinfRxLR dna $Assembly
```


Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	BlastHits=analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_avr_cds.fasta_homologs.csv
	HitsGff=analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_avr_cds.fasta_homologs.gff
	Column2=Avr_homolog
	NumHits=1
	$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
	BlastHits=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.csv
	HitsGff=analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.gff
	Column2=Chen_RxLR
	NumHits=1
	$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
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


# 6 Consolidating data

key files were grouped together into a single location. This directory was can
be shared with collaborators for further analysis.

```bash
	ProjDir=/home/groups/harrisonlab/project_files/idris

	# P. cactorum - Assembly
	Assembly_fa_Pcac=$ProjDir/repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
	Cegma_Gff=gene_pred/cegma/P.cactorum/10300/10300_dna_cegma.cegma.gff
	Cegma_report_Pcac=gene_pred/cegma/P.cactorum/10300/10300_dna_cegma.completeness_report
	Repeatmasked_txt_Pcac=$ProjDir/repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_hardmasked.fa.tbl
	Repeatmasked_Gff_Pcac=$ProjDir/repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_hardmasked.gff
	AssemblyDir=$ProjDir/collaboration/P.cactorum/10300/assembly
	mkdir -p $AssemblyDir
	cp $Assembly_fa_Pcac $AssemblyDir/.
	cp $Cegma_Gff $AssemblyDir/.
	cp $Cegma_report_Pcac $AssemblyDir/.
	cp $Repeatmasked_txt_Pcac $AssemblyDir/.
	cp $Repeatmasked_Gff_Pcac $AssemblyDir/.

	# P. cactorum - Braker1 genes
	Genes_aa_Pcac=$ProjDir/gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
	Genes_Gff_Pcac=$ProjDir/gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.gff
	InterPro_tsv_Pcac=$ProjDir/gene_pred/interproscan/P.cactorum/10300/10300_interproscan.tsv
	BrakerDir=$ProjDir/collaboration/P.cactorum/10300/gene_pred/braker
	mkdir -p $BrakerDir
	cp $Genes_aa_Pcac $BrakerDir/.
	cp $Genes_Gff_Pcac $BrakerDir/.
	cp $InterPro_tsv_Pcac $BrakerDir/.
	# P. cactorum - ORF genes
	ORF_aa_Pcac=$ProjDir/gene_pred/ORF_finder/P.cactorum/10300/10300.aa_cat.fa
	ORF_Gff_Pcac=$ProjDir/gene_pred/ORF_finder/P.cactorum/10300/10300_ORF_corrected.gff3
	ORFDir=$ProjDir/collaboration/P.cactorum/10300/gene_pred/ORFs
	mkdir -p $ORFDir
	cp $ORF_aa_Pcac $ORFDir/.
	cp $ORF_Gff_Pcac $ORFDir/.
	# P.cactorum - Secreted proteins
	SigP_fa_Braker_Pcac=$ProjDir/gene_pred/braker_sigP/P.cactorum/10300/10300_aug_sp.aa
	SigP_Gff_Braker_Pcac=
	SigP_fa_ORF_Pcac=$ProjDir/gene_pred/ORF_sigP/P.cactorum/10300/10300_ORF_sp.aa
	SigP_Gff_ORF_Pcac=$ProjDir/gene_pred/ORF_sigP/P.cactorum/10300/10300_ORF_sp_merged.gff
	SecretedDir=$ProjDir/collaboration/P.cactorum/10300/effectors/secreted_proteins
	mkdir -p $SecretedDir
	cp $SigP_fa_Braker_Pcac $SecretedDir/.
	cp $SigP_Gff_Braker_Pcac $SecretedDir/.
	cp $SigP_fa_ORF_Pcac $SecretedDir/.
	cp $SigP_Gff_ORF_Pcac $SecretedDir/.
	# P. cactorum - RxLR effectors
	RxLR_EER_fa_Braker_Pcac=$ProjDir/analysis/RxLR_effectors/RxLR_EER_regex_finder/P.cactorum/10300/10300_Aug_RxLR_EER_regex.fa
	RxLR_EER_Gff_Braker_Pcac=$ProjDir/analysis/RxLR_effectors/RxLR_EER_regex_finder/P.cactorum/10300/10300_Aug_RxLR_EER_regex.gff3
	RxLR_hmm_fa_Braker_Pcac=$ProjDir/analysis/RxLR_effectors/hmmer_RxLR/P.cactorum/10300/10300_Aug_RxLR_hmmer.fa
	RxLR_hmm_Gff_Braker_Pcac=$ProjDir/analysis/RxLR_effectors/hmmer_RxLR/P.cactorum/10300/10300_Aug_RxLR_hmmer.gff3
	RxLR_WY_fa_Braker_Pcac=$ProjDir/analysis/RxLR_effectors/hmmer_WY/P.cactorum/10300/10300_Aug_WY_hmmer.fa
	RxLR_WY_Gff_Braker_Pcac=$ProjDir/analysis/RxLR_effectors/hmmer_WY/P.cactorum/10300/10300_Aug_WY_hmmer.gff3
	RxLR_EER_fa_ORF_Pcac=$ProjDir/analysis/RxLR_effectors/RxLR_EER_regex_finder/P.cactorum/10300/10300_ORF_RxLR_EER_regex.fa
	RxLR_EER_Gff_ORF_Pcac=$ProjDir/analysis/RxLR_effectors/RxLR_EER_regex_finder/P.cactorum/10300/10300_ORF_RxLR_EER_regex.gff
	RxLR_hmm_fa_ORF_Pcac=$ProjDir/analysis/RxLR_effectors/hmmer_RxLR/P.cactorum/10300/10300_ORF_RxLR_hmmer.fa
	RxLR_hmm_Gff_ORF_Pcac=$ProjDir/analysis/RxLR_effectors/hmmer_RxLR/P.cactorum/10300/10300_ORF_RxLR_hmmer.gff3
	RxLR_WY_fa_ORF_Pcac=$ProjDir/analysis/RxLR_effectors/hmmer_WY/P.cactorum/10300/10300_ORF_WY_hmmer.fa
	RxLR_WY_Gff_ORF_Pcac=$ProjDir/analysis/RxLR_effectors/hmmer_WY/P.cactorum/10300/10300_ORF_WY_hmmer.gff
	RxLRDir=$ProjDir/collaboration/P.cactorum/10300/effectors/RxLRs
	mkdir -p $RxLRDir
	cp $RxLR_EER_fa_Braker_Pcac $RxLRDir/.
	cp $RxLR_EER_Gff_Braker_Pcac $RxLRDir/.
	cp $RxLR_hmm_fa_Braker_Pcac $RxLRDir/.
	cp $RxLR_hmm_Gff_Braker_Pcac $RxLRDir/.
	cp $RxLR_WY_fa_Braker_Pcac $RxLRDir/.
	cp $RxLR_WY_Gff_Braker_Pcac $RxLRDir/.
	cp $RxLR_EER_fa_ORF_Pcac $RxLRDir/.
	cp $RxLR_EER_Gff_ORF_Pcac $RxLRDir/.
	cp $RxLR_hmm_fa_ORF_Pcac $RxLRDir/.
	cp $RxLR_hmm_Gff_ORF_Pcac $RxLRDir/.
	cp $RxLR_WY_fa_ORF_Pcac $RxLRDir/.
	cp $RxLR_WY_Gff_ORF_Pcac $RxLRDir/.
	# P. cactorum - Crinkler effectors
	RxLR_CRN_fa_Braker_Pcac=$ProjDir/analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_aug_CRN_hmmer_out.fa
	RxLR_CRN_Gff_Braker_Pcac=$ProjDir/analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_Aug_CRN_hmmer.gff3
	RxLR_CRN_fa_ORF_Pcac=$ProjDir/analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORF_CRN_hmmer_out.fa
	RxLR_CRN_Gff_ORF_Pcac=$ProjDir/analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORF_CRN_hmmer.gff3
	CrinklerDir=$ProjDir/collaboration/P.cactorum/10300/effectors/crinklers
	mkdir -p $CrinklerDir
	cp $RxLR_CRN_fa_Braker_Pcac $CrinklerDir/.
	cp $RxLR_CRN_Gff_Braker_Pcac $CrinklerDir/.
	cp $RxLR_CRN_fa_ORF_Pcac $CrinklerDir/.
	cp $RxLR_CRN_Gff_ORF_Pcac $CrinklerDir/.
	# P. cactorum - Blast of known pathogenicity genes
	Chen_BlastHits_csv_Pcac=$ProjDir/analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.csv
	Chen_BlastHits_gff_Pcac=$ProjDir/analysis/blast_homology/P.cactorum/10300/10300_chen_et_al_2014_RxLR.fa_homologs.gff
	BlastDir=$ProjDir/collaboration/P.cactorum/10300/effectors/blast_hits
	mkdir -p $BlastDir
	cp $Chen_BlastHits_csv_Pcac $BlastDir/.
	cp $Chen_BlastHits_gff_Pcac $BlastDir/.

	# P. infestans
	# P. infestans - Assembly
	Assembly_fa_Pinf=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
	AssemblyDir=$ProjDir/collaboration/P.infestans/T30-4/assembly
	mkdir -p $AssemblyDir
	cp $Assembly_fa_Pcac $AssemblyDir/.
	# P. infestans - Published genes
	Genes_aa_Pinf=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
	GenesDir=$ProjDir/collaboration/P.infestans/T30-4/gene_pred/published
	mkdir -p $GenesDir
	cp $Genes_aa_Pcac $GenesDir/.
	# P. infestans - ORF genes
	ORF_aa_Pinf=gene_pred/ORF_finder/P.infestans/T30-4/T30-4.aa_cat.fa
	ORF_Gff_Pinf=gene_pred/ORF_finder/P.infestans/T30-4/T30-4_ORF_merged_corrected.gff
	ORFDir=$ProjDir/collaboration/P.infestans/T30-4/gene_pred/ORFs
	mkdir -p $ORFDir
	cp $ORF_aa_Pinf $ORFDir/.
	cp $ORF_Gff_Pinf $ORFDir/.
	# P.infestans - Secreted proteins
	SigP_fa_Pub_Pinf=gene_pred/published_sigP/P.infestans/T30-4/T30-4_pub_sp.aa
	SigP_fa_ORF_Pinf=gene_pred/ORF_sigP/P.infestans/T30-4/T30-4_ORF_sp.aa
	SigP_Gff_ORF_Pinf=gene_pred/ORF_sigP/P.infestans/T30-4/T30-4_ORF_sp_merged.gff
	SecretedDir=$ProjDir/collaboration/P.infestans/T30-4/effectors/secreted_proteins
	mkdir -p $SecretedDir
	cp $SigP_fa_Pub_Pinf $SecretedDir/.
	cp $SigP_fa_ORF_Pinf $SecretedDir/.
	cp $SigP_Gff_ORF_Pinf $SecretedDir/.
	# P. infestans - RxLR effectors
	RxLR_EER_fa_Pub_Pinf=$ProjDir/analysis/RxLR_effectors/RxLR_EER_regex_finder/P.infestans/T30-4/T30-4_pub_RxLR_EER_regex.fa
	RxLR_EER_Gff_Pub_Pinf=$ProjDir/analysis/RxLR_effectors/RxLR_EER_regex_finder/P.infestans/T30-4/T30-4_pub_RxLR_EER_regex.gff3
	RxLR_hmm_fa_Pub_Pinf=$ProjDir/analysis/RxLR_effectors/hmmer_RxLR/P.infestans/T30-4/T30-4_pub_RxLR_hmmer.fa
	RxLR_hmm_Gff_Pub_Pinf=$ProjDir/analysis/RxLR_effectors/hmmer_RxLR/P.infestans/T30-4/T30-4_pub_RxLR_hmmer.gff3
	RxLR_WY_fa_Pub_Pinf=$ProjDir/analysis/RxLR_effectors/hmmer_WY/P.infestans/T30-4/T30-4_pub_WY_hmmer.fa
	RxLR_WY_Gff_Pub_Pinf=$ProjDir/analysis/RxLR_effectors/hmmer_WY/P.infestans/T30-4/T30-4_pub_WY_hmmer.gff3
	RxLR_EER_fa_ORF_Pinf=$ProjDir/analysis/RxLR_effectors/RxLR_EER_regex_finder/P.infestans/T30-4/T30-4_ORF_RxLR_EER_regex.fa
	RxLR_EER_Gff_ORF_Pinf=$ProjDir/analysis/RxLR_effectors/RxLR_EER_regex_finder/P.infestans/T30-4/T30-4_ORF_RxLR_EER_regex.gff
	RxLR_hmm_fa_ORF_Pinf=$ProjDir/analysis/RxLR_effectors/hmmer_RxLR/P.infestans/T30-4/T30-4_ORF_RxLR_hmmer.fa
	RxLR_hmm_Gff_ORF_Pinf=$ProjDir/analysis/RxLR_effectors/hmmer_RxLR/P.infestans/T30-4/T30-4_ORF_RxLR_hmmer.gff3
	RxLR_WY_fa_ORF_Pinf=$ProjDir/analysis/RxLR_effectors/hmmer_WY/P.infestans/T30-4/T30-4_ORF_WY_hmmer.fa
	RxLR_WY_Gff_ORF_Pinf=$ProjDir/analysis/RxLR_effectors/hmmer_WY/P.infestans/T30-4/T30-4_ORF_WY_hmmer.gff
	RxLRDir=$ProjDir/collaboration/P.infestans/T30-4/effectors/RxLRs
	mkdir -p $RxLRDir
	cp $RxLR_EER_fa_Pub_Pinf $RxLRDir/.
	cp $RxLR_EER_Gff_Pub_Pinf $RxLRDir/.
	cp $RxLR_hmm_fa_Pub_Pinf $RxLRDir/.
	cp $RxLR_hmm_Gff_Pub_Pinf $RxLRDir/.
	cp $RxLR_WY_fa_Pub_Pinf $RxLRDir/.
	cp $RxLR_WY_Gff_Pub_Pinf $RxLRDir/.
	cp $RxLR_EER_fa_ORF_Pinf $RxLRDir/.
	cp $RxLR_EER_Gff_ORF_Pinf $RxLRDir/.
	cp $RxLR_hmm_fa_ORF_Pinf $RxLRDir/.
	cp $RxLR_hmm_Gff_ORF_Pinf $RxLRDir/.
	cp $RxLR_WY_fa_ORF_Pinf $RxLRDir/.
	cp $RxLR_WY_Gff_ORF_Pinf $RxLRDir/.
	# P. infestans - Crinkler effectors
	RxLR_CRN_fa_Pub_Pinf=$ProjDir/analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_pub_CRN_hmmer_out.fa
	RxLR_CRN_Gff_Pub_Pinf=$ProjDir/analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_pub_CRN_hmmer.gff3
	RxLR_CRN_fa_ORF_Pinf=$ProjDir/analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_ORF_CRN_hmmer_out.fa
	RxLR_CRN_Gff_ORF_Pinf=$ProjDir/analysis/CRN_effectors/hmmer_CRN/P.infestans/T30-4/T30-4_ORF_CRN_hmmer.gff3
	CrinklerDir=$ProjDir/collaboration/P.cactorum/10300/effectors/crinklers
	mkdir -p $CrinklerDir
	cp $RxLR_CRN_fa_Pub_Pinf $CrinklerDir/.
	cp $RxLR_CRN_Gff_Pub_Pinf $CrinklerDir/.
	cp $RxLR_CRN_fa_ORF_Pinf $CrinklerDir/.
	cp $RxLR_CRN_Gff_ORF_Pinf $CrinklerDir/.

	# P. parasitica
	# P. capsica
	# P. sojae

	sscp -r collaboration armita@149.155.32.11:/home/sftp_chroot/Pcactorum/.
```
