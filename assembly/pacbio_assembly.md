
## Data extraction


for P.cactorum data:
```bash
  cd /home/groups/harrisonlab/project_files/idris
  RawDatDir=/home/harrir/projects/pacbio_test/p_cact
  mkdir -p raw_dna/pacbio/P.cactorum/414
  cp -r $RawDatDir/A07_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/B07_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/G06_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/H06_1 raw_dna/pacbio/P.cactorum/414/.
  OutDir=raw_dna/pacbio/P.cactorum/414/extracted
  mkdir -p $OutDir
  cat raw_dna/pacbio/P.cactorum/414/*/Analysis_Results/*.subreads.fastq | gzip -cf > $OutDir/concatenated_pacbio.fastq.gz
  #
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/pacbio/Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage.tar.gz
  mkdir -p raw_dna/pacbio/P.cactorum/414
  cp -r $RawDat raw_dna/pacbio/P.cactorum/414/.
  cd raw_dna/pacbio/P.cactorum/414
  tar -zxvf Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage.tar.gz
  # Data for 414 was contained in E02_1, F02_1 and G02_1
  cat \
  Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage/E02_1/Analysis_Results/*.subreads.fastq \
  Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage/F02_1/Analysis_Results/*.subreads.fastq \
  Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage/G02_1/Analysis_Results/*.subreads.fastq \
  | gzip -cf > extracted/concatenated_pacbio_extra_coverage.fastq.gz
  rm Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage.tar.gz
  rm -r Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage
```
<!--
for P. fragariae data (commands for tom to run)

```bash
  cd /home/groups/harrisonlab/project_files/phytophthora_fragariae
  RawDatDir=/home/harrir/projects/pacbio_test/p_frag
  mkdir -p raw_dna/pacbio/P.fragariae/Bc16
  cp -r $RawDatDir/C07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  cp -r $RawDatDir/D07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  cp -r $RawDatDir/E07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  cp -r $RawDatDir/F07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  OutDir=raw_dna/pacbio/P.fragariae/Bc16/extracted
  mkdir -p $OutDir
  cat raw_dna/pacbio/P.fragariae/Bc16/*/Analysis_Results/*.subreads.fastq > $OutDir/concatenated_pacbio.fastq
``` -->

<!-- Sequencing depth and possible contamination was identified using kmer counting
kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

```bash
  for Reads in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
    echo $Reads
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/kmc_kmer_counting.sh $Reads
  done
``` -->

## Assembly


### Canu assembly

<!--
```bash
  Reads=$(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq.gz)
  Run1=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq.gz)
  Run2=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz)
  Reads=raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_both.fastq.gz
  cat $Run1 $Run2 > $Reads
  GenomeSz="65m"
  Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir=assembly/canu/$Organism/"$Strain"_3
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
```

```bash
  Run2=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz)
  GenomeSz="65m"
  Strain=$(echo $Run2 | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Run2 | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir=assembly/canu/$Organism/"$Strain"_4
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu.sh $Run2 $GenomeSz $Prefix $OutDir
```
-->

```bash
  Run1=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq.gz)
  Run2=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz)
  GenomeSz="65m"
  Strain=$(echo $Run1 | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Run1 | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu2
  OutDir=assembly/canu/$Organism"_2"/$Strain
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu_2lib.sh $Run1 $Run2 $GenomeSz $Prefix $OutDir
```


```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/canu/*/*/*.contigs.fasta | grep -e '_2' -e '_3'); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

Assemblies were polished using Pilon

```bash
  for Assembly in $(ls assembly/canu/*/*/*_canu.contigs.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
    OutDir=assembly/canu/$Organism/$Strain/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/canu/*/*/polished/pilon.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/canu/$Organism/$Strain/polished
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

### Spades Assembly

For P. cactorum

```bash
  # for PacBioDat in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
for PacBioDat in $(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz); do
echo $StrainPath
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
TrimF2_Read=$(ls $IlluminaDir/F/414_run2_F_trim.fq.gz);
TrimR2_Read=$(ls $IlluminaDir/R/414_run2_R_trim.fq.gz);
TrimF3_Read=$(ls $IlluminaDir/F/414_170210_F_trim.fq.gz);
TrimR3_Read=$(ls $IlluminaDir/R/414_170210_R_trim.fq.gz);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
# OutDir=assembly/spades_pacbio/$Organism/$Strain
OutDir=assembly/spades_pacbio/$Organism/"$Strain"_4
qsub $ProgDir/subSpades_3lib_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir 28
done
```

Contigs shorter than 500bp were removed from the assembly

```bash
for Contigs in $(ls assembly/spades_pacbio/*/*/contigs.fasta); do
AssemblyDir=$(dirname $Contigs)
mkdir $AssemblyDir/filtered_contigs
FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
$FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs/contigs_min_500bp.fasta
done
```

Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades_pacbio/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```



## Merging pacbio and hybrid assemblies

```bash
  # for PacBioAssembly in $(ls assembly/canu/*/*/polished/pilon.fasta); do
  for PacBioAssembly in $(ls assembly/canu/P.cactorum/414*/414_canu.contigs.fasta | grep '_4'); do
    Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
    # HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
    HybridAssembly=$(ls assembly/spades_pacbio/P.cactorum/414_4/contigs.fasta)
    OutDir=assembly/merged_canu_spades/$Organism/$Strain
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
    qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir
  done
```

This merged assembly was polished using Pilon

```bash
  for Assembly in $(ls assembly/merged_canu_spades/*/*/merged.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir
  done
```

Contigs were renamed in accordance with ncbi recomendations.

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/merged_canu_spades/*/*/polished/pilon.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/filtered_contigs
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/merged_canu_spades/*/*/filtered_contigs/contigs_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

# Preliminary analysis

## Checking PacBio coverage against P414 contigs

The accuracy of PacBio assembly pipelines is currently unknown. To help identify
regions that may have been missassembled the pacbio reads were aligned back to
the assembled genome. Coverage was determined using bedtools genomecov and
regions with low coverage flagged using a python script flag_low_coverage.py.
These low coverage regions were visually inspected using IGV.

```bash
  Assembly=$(ls assembly/merged_canu_spades/P.cactorum/414/filtered_contigs/contigs_min_500bp_renamed.fasta)
  Reads=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq)
  OutDir=analysis/genome_alignment/bwa/P.cactorum/414/vs_414
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
  #
  # AlignedBam=$OutDir/414_contigs_renamed.fasta_aligned_sorted.bam.gz
  # CoverageTxt=$OutDir/414_bp_genome_cov.txt
  # bedtools genomecov -max 5 -bga -d -ibam $AlignedBam -g $Assembly > $CoverageTxt
  #
  # Threshold=5
  # FlaggedRegions=$OutDir/414_flagged_regions.txt
  # $ProgDir/flag_low_coverage.py --genomecov $CoverageTxt --min $Threshold > $FlaggedRegions
```

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  for BestAss in $(ls assembly/merged_canu_spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -w -e '414_v2'); do
    qsub $ProgDir/rep_modeling.sh $BestAss
    qsub $ProgDir/transposonPSI.sh $BestAss
  done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
  for RepDir in $(ls -d repeat_masked/P.*/*/filtered_contigs_repmask | grep -w -e '414_v2'); do
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
  P.cactorum	414
  The number of bases masked by RepeatMasker:	16199492
  The number of bases masked by TransposonPSI:	4951212
  The total number of masked bases are:	17347836
```


# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

## Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  for Genome in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
    echo $Genome;
    qsub $ProgDir/sub_cegma.sh $Genome dna;
  done
```

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/*/*/*_dna_cegma.completeness_report | grep -w -e '414_v2'); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
```

#Gene prediction

Gene prediction was performed for the P. cactorum genome. Two gene prediction
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
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
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

<!--
#### Braker prediction

```bash
    for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414'); do
    	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    	echo "$Organism - $Strain"
    	mkdir -p alignment/$Organism/$Strain/concatenated
    	samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
      alignment/$Organism/$Strain/SRR1206032/accepted_hits.bam \
    	alignment/$Organism/$Strain/SRR1206033/accepted_hits.bam
    	OutDir=gene_pred/braker/$Organism/"$Strain"_braker_pacbio
    	AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
    	GeneModelName="$Organism"_"$Strain"_braker_pacbio
    	rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_pacbio
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
	for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414'); do
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
	for Genome in $(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa); do
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

#Functional annotation

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
	cd /home/groups/harrisonlab/project_files/idris
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Proteins in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e 'F414'); do
		Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		echo $Strain
		InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
		$ProgDir/append_interpro.sh $Proteins $InterProRaw
	done
```

## B) SwissProt

```bash
  screen -a
  qlogin
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  OutDir=$ProjDir/gene_pred/swissprot/$Species/$Strain/
  mkdir -p $OutDir
  blastp \
  -db /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot \
  -query $ProjDir/$Proteome \
  -out $OutDir/swissprot_v2015_10_hits.tbl \
  -evalue 1e-100 \
  -outfmt 6 \
  -num_threads 16 \
  -num_alignments 10
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
# * E) From ORF fragments - Signal peptide & RxLR motif  
# * F) From ORF fragments - Hmm evidence of WY domains  
# * G) From ORF fragments - Hmm evidence of RxLR effectors


 ### A) From Augustus gene models - Signal peptide & RxLR motif

 Required programs:
  * SigP
  * biopython

 #### A.1) Signal peptide prediction using SignalP 2.0

 Proteins that were predicted to contain signal peptides were identified using
 the following commands:

 ```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
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
      while [ $Jobs -gt 1 ]; do
        sleep 3
        printf "."
        Jobs=$(qstat | grep 'pred_sigP' | grep -w 'qw'| wc -l)
      done
      printf "\n"
      echo $File
      # qsub $ProgDir/pred_sigP.sh $File
      qsub $ProgDir/pred_sigP.sh $File signalp-4.1
    done
  done
 ```

 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:
```bash
  SplitDir=gene_pred/final_split/P.cactorum/414
  Strain=$(echo $SplitDir | cut -d '/' -f4)
  Organism=$(echo $SplitDir | cut -d '/' -f3)
  InStringAA=''
  InStringNeg=''
  InStringTab=''
  InStringTxt=''
  for SigpDir in $(ls -d gene_pred/final_sigP* | cut -f2 -d'/'); do
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
```

#### B.2) Prediction using Phobius

 Secreted proteins were also predicted using Phobius

 ```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/phobius/$Organism/$Strain
    mkdir -p $OutDir
    phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
    # cat $OutDir/"$Strain"_phobius.txt | grep -B1 'SIGNAL' | grep 'ID' | sed s'/ID.*g/g/g' > $OutDir/"$Strain"_phobius_headers.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    $ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
  done
 ```


Secreted proteins from different sources were combined into a single file:

```bash
for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  OutDir=gene_pred/combined_sigP/$Organism/$Strain
  cat gene_pred/*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | tr -d ' ' > $OutDir/"$Strain"_all_secreted.fa
  cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | sort -g | rev | uniq -s42 | rev | tr -d '>' > $OutDir/"$Strain"_secreted_unique_headers.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $OutDir/"$Strain"_all_secreted.fa --headers $OutDir/"$Strain"_secreted_unique_headers.txt > $OutDir/"$Strain"_secreted.fa
done
```


```bash
for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  OutDir=gene_pred/combined_sigP/$Organism/$Strain
  mkdir -p $OutDir
  cat gene_pred/*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
  cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
done
```

 The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

 The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash

for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep -w -e '414'); do
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
cat $Secretome | grep '>' | cut -f1 | sort | uniq | wc -l;

printf "the number of SigP-RxLR genes are:\t";
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_all_secreted_RxLR_regex.fa;
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | cut -f1 | tr -d '>' | sort -g | uniq > $OutDir/"$Strain"_RxLR_regex.txt
cat $OutDir/"$Strain"_RxLR_regex.txt | wc -l

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa


printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sort -g | uniq > $OutDir/"$Strain"_RxLR_EER_regex.txt
cat $OutDir/"$Strain"_RxLR_EER_regex.txt | wc -l
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa

printf "\n"

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_regex.txt
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_EER_regex.txt
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_Aug_RxLR_regex.txt $Gff RxLR_regex ID AUGUSTUS > $OutDir/"$Strain"_Aug_RxLR_regex.gff
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt $GeneModels RxLR_regex ID AUGUSTUS > $OutDir/"$Strain"_Aug_RxLR_EER_regex.gff

cat $Gff | grep -w -f $OutDir/"$Strain"_Aug_RxLR_regex.txt > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
cat $Gff | grep -w -f $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt > $OutDir/"$Strain"_Aug_RxLR_EER_regex.gff3
done
```
```
  strain: 414	species: P.cactorum
  the total number of SigP gene is:	10267
  the number of unique SigP gene is:	5555
  the number of SigP-RxLR genes are:	449
  the number of SigP-RxLR-EER genes are:	204
```


### G) From Secreted gene models - Hmm evidence of RxLR effectors

```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
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
  P.cactorum 414
  Initial search space (Z):              32832  [actual number of targets]
  Domain search space  (domZ):             146  [number of targets reported over threshold]
```

### F) Combining RxLRs from Regex and hmm searches


The total RxLRs are

```bash
  for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_RxLR_EER_regex.txt | grep -v 'Aug' | grep '414'); do
    Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
    Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
    Proteome=$(ls gene_pred/codingquary/$Organism/$Strain/*/final_genes_combined.pep.fasta)
    HmmRxLR=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer_headers.txt
    echo "Number of RxLRs identified by Regex:"
    cat $RegexRxLR | wc -l
    echo "Number of RxLRs identified by Hmm:"
    cat $HmmRxLR | wc -l
    echo "Number of RxLRs in combined dataset:"
    cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
    echo "Number of RxLRs in both datasets:"
    cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
    echo "Extracting RxLRs from datasets"
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
  Number of RxLRs identified by Regex:
  204
  Number of RxLRs identified by Hmm:
  145
  Number of RxLRs in combined dataset:
  230
  Number of RxLRs in both datasets:
  118
  Extracting RxLRs from datasets
  Number of genes in the extracted gff file:
  230
```

### D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:


```bash
HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
LFLAK_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm)
DWL_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm)
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
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
  P.cactorum - 414
  Initial search space (Z):              32832  [actual number of targets]
  Domain search space  (domZ):             218  [number of targets reported over threshold]
  Initial search space (Z):              32832  [actual number of targets]
  Domain search space  (domZ):             189  [number of targets reported over threshold]
  170
```

Extract gff annotations for Crinklers:

```bash
  for CRNlist in $(ls analysis/CRN_effectors/hmmer_CRN/*/*/*_pub_CRN_LFLAK_DWL.txt | grep '414'); do
    Strain=$(echo $CRNlist | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $CRNlist | rev | cut -f3 -d '/' | rev)
    OutName=$(echo $CRNlist | sed 's/.txt/.gff/g')
    echo "$Organism - $Strain"
    Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
    cat $CRNlist | sed -r 's/\.t.$//g' > tmp.txt
    cat $Gff | grep -w -f tmp.txt >  $OutName
    rm tmp.txt
  done
```
-->



Contigs 1, 94, 113 and 131 were also identified as contaminants from a
Paenibacilluss sp.

Fistly the number of genes and RxLRs on these contigs were identified from results so far:



```bash
echo "The total number of transcripts was:"
cat gene_pred/codingquary/P.cactorum/414/final/final_genes_appended.gff3 | grep -w 'mRNA' | wc -l
echo "The number of transcripts on contigs 1,94, 113 and 131 were:"
cat gene_pred/codingquary/P.cactorum/414/final/final_genes_appended.gff3 | grep -w 'mRNA' | grep -v -w -e 'contig_1' -e 'contig_94' -e 'contig_113' -e 'contig_131' | wc -l

echo "The number of secreted genes excluding contaminant contigs:"

cat gene_pred/combined_sigP/P.cactorum/414/414_all_secreted.fa | grep '>' | tr -d '>' > gene_pred/combined_sigP/P.cactorum/414/414_all_secreted_headers.txt

SigP_headers=gene_pred/combined_sigP/P.cactorum/414/414_all_secreted_headers.txt
Gff=gene_pred/codingquary/P.cactorum/414/final/final_genes_appended.gff3
SigP_Gff=gene_pred/combined_sigP/P.cactorum/414/414_secreted_unique.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $SigP_headers $Gff SigP ID > $SigP_Gff

cat $SigP_Gff | grep 'mRNA' | grep -v -w -e 'contig_1' -e 'contig_94' -e 'contig_113' -e 'contig_131' |  wc -l

echo "The number of ORFs were predicted, excluding contmainant contigs:"
cat gene_pred/ORF_finder/P.cactorum/414/414_ORF.gff3 | grep -w 'transcript' | grep -v -w -e 'contig_1' -e 'contig_94' -e 'contig_113' -e 'contig_131' | wc -l

echo "The number of RxLRs predicted, excluding contmainant contigs:"
cat analysis/RxLR_effectors/combined_evidence/P.cactorum/414/414_total_RxLR.gff | grep -w 'gene' | grep -v -w -e 'contig_1' -e 'contig_94' -e 'contig_113' -e 'contig_131' | wc -l

echo "The number of CRNs predicted, on LS contigs was observed using:"
cat analysis/CRN_effectors/hmmer_CRN/P.cactorum/414/414_pub_CRN_LFLAK_DWL.gff | grep -w -e 'contig_1' -e 'contig_94' -e 'contig_113' -e 'contig_131' | wc -l
```

Contigs were renamed in accordance with ncbi recomendations.
Then contaminant contigs were removed from the assembly and analysis re-run:

```bash
  touch tmp.csv
  printf "contig_1\t.\texclude\t.\t.\tConaminant contig\n" >> tmp.csv
  printf "contig_94\t.\texclude\t.\t.\tConaminant contig\n" >> tmp.csv
  printf "contig_113\t.\texclude\t.\t.\tConaminant contig\n" >> tmp.csv
  printf "contig_131\t.\texclude\t.\t.\tConaminant contig\n" >> tmp.csv

  for Assembly in $(ls assembly/merged_canu_spades/*/*/filtered_contigs/contigs_min_500bp_renamed.fasta | grep -v -e '_v2'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/merged_canu_spades/$Organism/"$Strain"_v2/filtered_contigs
    mkdir -p $OutDir
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/merged_canu_spades/*/*_v2/filtered_contigs/contigs_min_500bp_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=assembly/merged_canu_spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  for BestAss in $(ls assembly/merged_canu_spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -w -e '414_v2'); do
    qsub $ProgDir/rep_modeling.sh $BestAss
    qsub $ProgDir/transposonPSI.sh $BestAss
  done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
  for RepDir in $(ls -d repeat_masked/P.*/*/filtered_contigs_repmask | grep -w -e '414'); do
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

```


# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

## Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  # for Genome in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
  for Genome in $(ls assembly/merged_canu_spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -w -e '414_v2'); do
    echo $Genome;
    qsub $ProgDir/sub_cegma.sh $Genome dna;
  done
```

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/*/*/*_dna_cegma.completeness_report | grep -w -e '414_v2'); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
```

#Gene prediction

Gene prediction was performed for the P. cactorum genome. Two gene prediction
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
  # for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
  for Assembly in $(ls assembly/merged_canu_spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -w -e '414_v2'); do
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
    for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414'); do
    	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    	echo "$Organism - $Strain"
    	mkdir -p alignment/$Organism/$Strain/concatenated
    	samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
      alignment/$Organism/$Strain/SRR1206032/accepted_hits.bam \
    	alignment/$Organism/$Strain/SRR1206033/accepted_hits.bam
    	OutDir=gene_pred/braker/$Organism/"$Strain"_braker_pacbio
    	AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
    	GeneModelName="$Organism"_"$Strain"_braker_pacbio
    	rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_pacbio
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
	for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414'); do
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
	# for Genome in $(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa); do
  for Assembly in $(ls assembly/merged_canu_spades/*/*/filtered_contigs/contigs_min_500bp_renamed.fasta | grep -v -e '_v2'); do
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
#Functional annotation

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
	cd /home/groups/harrisonlab/project_files/idris
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Proteins in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e 'F414'); do
		Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		echo $Strain
		InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
		$ProgDir/append_interpro.sh $Proteins $InterProRaw
	done
```

## B) SwissProt

```bash
  screen -a
  qlogin
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  OutDir=$ProjDir/gene_pred/swissprot/$Species/$Strain/
  mkdir -p $OutDir
  blastp \
  -db /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot \
  -query $ProjDir/$Proteome \
  -out $OutDir/swissprot_v2015_10_hits.tbl \
  -evalue 1e-100 \
  -outfmt 6 \
  -num_threads 16 \
  -num_alignments 10
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
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
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
      while [ $Jobs -gt 1 ]; do
        sleep 3
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
  SplitDir=gene_pred/final_split/P.cactorum/414
  Strain=$(echo $SplitDir | cut -d '/' -f4)
  Organism=$(echo $SplitDir | cut -d '/' -f3)
  InStringAA=''
  InStringNeg=''
  InStringTab=''
  InStringTxt=''
  for SigpDir in $(ls -d gene_pred/final_sigP* | cut -f2 -d'/'); do
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
```

#### B.2) Prediction using Phobius

 Secreted proteins were also predicted using Phobius

 ```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/phobius/$Organism/$Strain
    mkdir -p $OutDir
    phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
    # cat $OutDir/"$Strain"_phobius.txt | grep -B1 'SIGNAL' | grep 'ID' | sed s'/ID.*g/g/g' > $OutDir/"$Strain"_phobius_headers.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    $ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
  done
 ```


Secreted proteins from different sources were combined into a single file:

```bash
for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  OutDir=gene_pred/combined_sigP/$Organism/$Strain
  cat gene_pred/*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | tr -d ' ' > $OutDir/"$Strain"_all_secreted.fa
  cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | sort -g | rev | uniq -s42 | rev | tr -d '>' > $OutDir/"$Strain"_secreted_unique_headers.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $OutDir/"$Strain"_all_secreted.fa --headers $OutDir/"$Strain"_secreted_unique_headers.txt > $OutDir/"$Strain"_secreted.fa
done
```


```bash
for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
  Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
  OutDir=gene_pred/combined_sigP/$Organism/$Strain
  mkdir -p $OutDir
  cat gene_pred/*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
  cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
done
```

 The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

 The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash

for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep -w -e '414'); do
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
cat $Secretome | grep '>' | cut -f1 | sort | uniq | wc -l;

printf "the number of SigP-RxLR genes are:\t";
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_all_secreted_RxLR_regex.fa;
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | cut -f1 | tr -d '>' | sort -g | uniq > $OutDir/"$Strain"_RxLR_regex.txt
cat $OutDir/"$Strain"_RxLR_regex.txt | wc -l

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa


printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sort -g | uniq > $OutDir/"$Strain"_RxLR_EER_regex.txt
cat $OutDir/"$Strain"_RxLR_EER_regex.txt | wc -l
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa

printf "\n"

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_regex.txt
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_EER_regex.txt
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_Aug_RxLR_regex.txt $Gff RxLR_regex ID AUGUSTUS > $OutDir/"$Strain"_Aug_RxLR_regex.gff
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt $GeneModels RxLR_regex ID AUGUSTUS > $OutDir/"$Strain"_Aug_RxLR_EER_regex.gff

cat $Gff | grep -w -f $OutDir/"$Strain"_Aug_RxLR_regex.txt > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
cat $Gff | grep -w -f $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt > $OutDir/"$Strain"_Aug_RxLR_EER_regex.gff3
done
```
```
  strain: 414	species: P.cactorum
  the total number of SigP gene is:	10267
  the number of unique SigP gene is:	5555
  the number of SigP-RxLR genes are:	449
  the number of SigP-RxLR-EER genes are:	204
```


### G) From Secreted gene models - Hmm evidence of RxLR effectors

```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
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
  P.cactorum 414
  Initial search space (Z):              32832  [actual number of targets]
  Domain search space  (domZ):             146  [number of targets reported over threshold]
```

### F) Combining RxLRs from Regex and hmm searches


The total RxLRs are

```bash
  for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_RxLR_EER_regex.txt | grep -v 'Aug' | grep '414'); do
    Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
    Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
    Proteome=$(ls gene_pred/codingquary/$Organism/$Strain/*/final_genes_combined.pep.fasta)
    HmmRxLR=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer_headers.txt
    echo "Number of RxLRs identified by Regex:"
    cat $RegexRxLR | wc -l
    echo "Number of RxLRs identified by Hmm:"
    cat $HmmRxLR | wc -l
    echo "Number of RxLRs in combined dataset:"
    cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
    echo "Number of RxLRs in both datasets:"
    cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
    echo "Extracting RxLRs from datasets"
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
  Number of RxLRs identified by Regex:
  204
  Number of RxLRs identified by Hmm:
  145
  Number of RxLRs in combined dataset:
  230
  Number of RxLRs in both datasets:
  118
  Extracting RxLRs from datasets
  Number of genes in the extracted gff file:
  230
```

### D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:


```bash
HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
LFLAK_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm)
DWL_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm)
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414'); do
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
  P.cactorum - 414
  Initial search space (Z):              32832  [actual number of targets]
  Domain search space  (domZ):             218  [number of targets reported over threshold]
  Initial search space (Z):              32832  [actual number of targets]
  Domain search space  (domZ):             189  [number of targets reported over threshold]
  170
```

Extract gff annotations for Crinklers:

```bash
  for CRNlist in $(ls analysis/CRN_effectors/hmmer_CRN/*/*/*_pub_CRN_LFLAK_DWL.txt | grep '414'); do
    Strain=$(echo $CRNlist | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $CRNlist | rev | cut -f3 -d '/' | rev)
    OutName=$(echo $CRNlist | sed 's/.txt/.gff/g')
    echo "$Organism - $Strain"
    Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
    cat $CRNlist | sed -r 's/\.t.$//g' > tmp.txt
    cat $Gff | grep -w -f tmp.txt >  $OutName
    rm tmp.txt
  done
```
