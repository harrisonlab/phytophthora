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
  for RawData in $(ls raw_dna/paired/P./*/*/*.fastq.gz); do
  # for RawData in $(ls raw_dna/paired/P.cactorum/*/*/*.fastq.gz | grep -v '10300'); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```


Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
  for StrainPath in $(ls -d raw_dna/paired/*/* | grep -e '415' -e '416' -e '62471'); do
  # for StrainPath in $(ls -d raw_dna/paired/*/* | grep -e 'idaei'); do
    echo $StrainPath
    Read_F=$(ls $StrainPath/F/*.fastq.gz)
    Read_R=$(ls $StrainPath/R/*.fastq.gz)
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
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
  IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa

  echo "404"
  StrainPath=raw_dna/paired/P.cactorum/404
  ReadsF=$(ls $StrainPath/F/130624_cactp404_S3_L001_R1_001.fastq.gz)
  ReadsR=$(ls $StrainPath/R/130624_cactp404_S3_L001_R2_001.fastq.gz)
  qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
  StrainPath=raw_dna/paired/P.cactorum/404
  ReadsF=$(ls $StrainPath/F/cactp404_S3_L001_R1_001.fastq.gz)
  ReadsR=$(ls $StrainPath/R/cactp404_S3_L001_R2_001.fastq.gz)
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


#Assembly

Assembly was performed with:
* Spades

## Spades Assembly

```bash
  for StrainPath in $(ls -d qc_dna/paired/P.cactorum/* | grep -e '62471'); do
  # for StrainPath in $(ls -d qc_dna/paired/P.idaei/*); do
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
    qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct 10
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
  for StrainPath in $(ls -d qc_dna/paired/P.*/404); do  
    echo $StrainPath
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $StrainPath/F/130624_cactp404_S3_L001_R1_001_trim.fq.gz);
    TrimR1_Read=$(ls $StrainPath/R/130624_cactp404_S3_L001_R2_001_trim.fq.gz);
    TrimF2_Read=$(ls $StrainPath/F/cactp404_S3_L001_R1_001_trim.fq.gz);
    TrimR2_Read=$(ls $StrainPath/R/cactp404_S3_L001_R2_001_trim.fq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    OutDir=assembly/spades/$Organism/$Strain
    qsub $ProgDir/subSpades_2lib_HiMem.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir correct 10
  done
```


Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep 'P.idaei'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```


```bash
  for Assembly in $(ls assembly/spades/P.*/*/*/contigs_min_500bp_renamed.fasta | grep -e 'P.idaei' -e 'P.cactorum'); do
    Kmer=$(echo $Assembly | rev | cut -f2 -d '/' | rev);
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir;
  done
```



# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
	for BestAss in $(ls assembly/spades/P.*/*/*/contigs_min_500bp_renamed.fasta | grep -e 'P.idaei' -e 'P.cactorum'); do
		qsub $ProgDir/rep_modeling.sh $BestAss
		qsub $ProgDir/transposonPSI.sh $BestAss
	done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
for RepDir in $(ls -d repeat_masked/P.*/*/filtered_contigs_repmask | grep -e 'P.cactorum'); do
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
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  cd /home/groups/harrisonlab/project_files/idris
  for Genome in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -e 'P.cactorum'); do
    echo $Genome;
    qsub $ProgDir/sub_cegma.sh $Genome dna;
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
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -v -w -e '414' -e 'P.fragariae'); do
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

```bash
    for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -v -w -e '414' -e 'P.fragariae'); do
    	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    	echo "$Organism - $Strain"
    	mkdir -p alignment/$Organism/$Strain/concatenated
    	samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
      alignment/$Organism/$Strain/SRR1206032/accepted_hits.bam \
    	alignment/$Organism/$Strain/SRR1206033/accepted_hits.bam
    	OutDir=gene_pred/braker/$Organism/"$Strain"_braker
    	AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
    	GeneModelName="$Organism"_"$Strain"_braker
    	rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
    	qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  	done
```


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
	for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -v -w -e '414' -e 'P.fragariae'); do
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
	for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa | grep -w -e '414'); do
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
	# for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff3 | grep -w -e 'Fus2'); do
	for BrakerGff in $(ls gene_pred/braker/P.*/*_braker_pacbio/*/augustus.gff3 | grep -e '414'); do
		Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g')
		Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		# BrakerGff=gene_pred/braker/$Organism/$Strain/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
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
		# GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/additional/additional_genes.gff
		# GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/out/PredictedPass.gff3

		$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
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

		# cat $BrakerGff $AddDir/additional_gene_parsed.gff3 | bedtools sort > $FinalGff
	done
```

The final number of genes per isolate was observed using:
```bash
  for DirPath in $(ls -d gene_pred/codingquary/P.*/*/final | grep -w -e '414'); do
    echo $DirPath;
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
    echo "";
  done
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
	ORF_Gff=gene_pred/ORF_finder/P.cactorum/10300/10300_ORF.gff
	ORF_Gff_mod=gene_pred/ORF_finder/P.cactorum/10300/10300_ORF_corrected.gff3
	$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
```
