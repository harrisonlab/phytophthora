#Published Phythopthora isolates
==========

Scripts used for the analysis of the P. cactorum isolate 10300.
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/idris
. As such this script relies upon adhering to a
specific directory structure.

The following is a summary of the work presented in this Readme.

The following processes were applied to the P. cactorum 10300 genome prior to analysis:
<!-- Data qc -->
<!-- Genome assembly -->
Repeatmasking
Gene prediction
Functional annotation

Analyses performed on these genomes involved BLAST searching for:
Known Avr genes
Ab initio prediction of putative RxLR, CRN and SSC effectors.



#Making local repositories for data

Data was downloaded for P. infestans isoalte T30-4 from:
http://protists.ensembl.org/Phytophthora_infestans/Info/Index
ftp://ftp.ensemblgenomes.org/pub/protists/release-29/fasta/phytophthora_infestans/dna/
ftp://ftp.ensemblgenomes.org/pub/protists/release-29/fasta/phytophthora_infestans
This represented 4921 contigs and 17787 proteins
<!-- The previous release is available at:
http://www.broadinstitute.org/annotation/genome/phytophthora_infestans/MultiDownloads.html
The current release represents 4,921 contigs and 18140 proteins. -->
This data was moved into the following directories:

```bash
  Organism=P.infestans
	Strain=T30-4
	ProjDir=/home/groups/harrisonlab/project_files/idris
	cd $ProjDir
  mkdir -p assembly/external_group/$Organism/$Strain/cdna
  mkdir -p assembly/external_group/$Organism/$Strain/cds
  mkdir -p assembly/external_group/$Organism/$Strain/dna
  mkdir -p assembly/external_group/$Organism/$Strain/ncrna
  mkdir -p assembly/external_group/$Organism/$Strain/pep
```

Data was downloaded for P. parasitica isolate 310 from
http://www.broadinstitute.org/annotation/genome/Phytophthora_parasitica/MultiDownloads.html
This represented 708 contigs and 20822 proteins
<!-- The current release is available at:
http://www.ncbi.nlm.nih.gov/Traces/wgs/?val=AGFV02
The current release is 82 Mb, represents 708 contigs and 27,942 genes (Unpublished) -->
This data was moved into the following directories:

```bash
  Organism=P.parasitica
	Strain=310
	ProjDir=/home/groups/harrisonlab/project_files/idris
	cd $ProjDir
  mkdir -p assembly/external_group/$Organism/$Strain/cdna
  mkdir -p assembly/external_group/$Organism/$Strain/cds
  mkdir -p assembly/external_group/$Organism/$Strain/dna
  mkdir -p assembly/external_group/$Organism/$Strain/ncrna
  mkdir -p assembly/external_group/$Organism/$Strain/pep
```

Data was downloaded for P. capsici isolate LT1534 from:
http://genome.jgi.doe.gov/Phyca11/Phyca11.download.ftp.html
<!-- This represented 917 contigs and 19805 genes
This was the current release. -->
The current release is 56 Mb represents 917 contigs and 19805 genes.
This data was moved into the following directories:

```bash
  Organism=P.capsici
	Strain=LT1534
	ProjDir=/home/groups/harrisonlab/project_files/idris
	cd $ProjDir
	mkdir -p assembly/external_group/$Organism/$Strain/dna
  mkdir -p assembly/external_group/$Organism/$Strain/pep
  mkdir -p assembly/external_group/$Organism/$Strain/annot
```

Data was downloaded for P. sojae isolate  from:
<!-- http://protists.ensembl.org/Phytophthora_sojae/Info/Index
ftp://ftp.ensemblgenomes.org/pub/protists/release-29/fasta/phytophthora_sojae/dna/
ftp://ftp.ensemblgenomes.org/pub/protists/release-29/fasta/phytophthora_sojae
This represented 1810 contigs and 18969 genes
The current version of the genome was downloaded from: -->
http://genome.jgi.doe.gov/Physo3/Physo3.info.html
The current release is 82Mb represents 83 contigs and 26584 genes.
This data was moved into the following directories:

```bash
  Organism=P.sojae
	Strain=P6497
	ProjDir=/home/groups/harrisonlab/project_files/idris
	cd $ProjDir
  mkdir -p assembly/external_group/$Organism/$Strain/cdna
  mkdir -p assembly/external_group/$Organism/$Strain/cds
  mkdir -p assembly/external_group/$Organism/$Strain/dna
  mkdir -p assembly/external_group/$Organism/$Strain/ncrna
  mkdir -p assembly/external_group/$Organism/$Strain/pep
```

## Parsing data

To run the gene prediction and functional annotation scripts all spaces and pipe
symbols had to be removed from the headers of assembly fasta files. This was
performed using the following commands:

```bash
  for File in $(ls assembly/external_group/P.*/*/dna/*.genome.fa); do
    OutFile=$(echo $File | sed 's/.fa/.parsed.fa/g');
    echo $OutFile; cat $File | sed 's/ /_/g' | sed 's/|/_/g' > $OutFile;
  done
```

# Assessing published genome assemblies

```bash
  Pinf_ass=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
  Ppar_ass=assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra-310_2_supercontigs.fasta
  Pcap_ass=assembly/external_group/P.capsici/LT1534/dna/Phyca11_unmasked_genomic_scaffolds.fasta
  Psoj_ass=assembly/external_group/P.sojae/P6497/dna/Physo3_AssemblyScaffolds.fasta
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $Pinf_ass $Ppar_ass $Pcap_ass $Psoj_ass; do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/external_group/$Organism/$Strain/dna/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```



# Repeat masking

Repeat masking was performed on published genomes to allow comparison of
repetative between species, including P. cactorum.

```bash
  # Pinf_ass=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
  # Ppar_ass=phytophthora_parasitica_inra-310_2_supercontigs_parsed.fasta
  Pcap_ass=assembly/external_group/P.capsici/LT1534/dna/Phyca11_unmasked_genomic_scaffolds.fasta
  Psoj_ass=assembly/external_group/P.sojae/P6497/dna/Physo3_AssemblyScaffolds.fasta
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  # for Genome in $Pinf_ass $Ppar_ass $Pcap_ass $Psoj_ass; do
  for Genome in $Ppar_ass $Pcap_ass $Psoj_ass; do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Genome
    qsub $ProgDir/transposonPSI.sh $Genome
  done
```

```bash
  Pinf_ass=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
  Pinf_ass_mod=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed_for_repmask.fa
  Ppar_ass=assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra-310_2_supercontigs.fasta
  Ppar_ass_mod=assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra-310_2_supercontigs_parsed.fasta
  cat $Pinf_ass | cut -f1 -d ':' | sed 's/\./_/g'> $Pinf_ass_mod
  cat $Ppar_ass | cut -f1 -d ' ' > $Ppar_ass_mod
  for Genome in $Pinf_ass_mod $Ppar_ass_mod; do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Genome
    qsub $ProgDir/transposonPSI.sh $Genome
  done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
  RepPcac=repeat_masked/P.cactorum/10300/10300_abyss_53_repmask
  RepPinf=repeat_masked/P.infestans/T30-4/dna_repmask
  RepPpar=repeat_masked/P.parisitica/310/dna_repmask
  RepPcap=repeat_masked/P.capsici/LT1534/dna_repmask    
  RepPsoj=repeat_masked/P.sojae/P6497/dna_repmask
  for RepDir in $RepPcac $RepPinf $RepPpar $RepPcap $RepPsoj; do
  # for RepDir in $RepPsoj; do
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
  P.cactorum	10300
  The number of bases masked by RepeatMasker:	9932231
  The number of bases masked by TransposonPSI:	2318917
  The total number of masked bases are:	10813641
  P.infestans	T30-4
  The number of bases masked by RepeatMasker:	123431647
  The number of bases masked by TransposonPSI:	28631065
  The total number of masked bases are:	152062712
  P.parisitica	310
  The number of bases masked by RepeatMasker:	5712393
  The number of bases masked by TransposonPSI:	2413506
  The total number of masked bases are:	6999606
  P.capsici	LT1534
  The number of bases masked by RepeatMasker:	12603095
  The number of bases masked by TransposonPSI:	3902251
  The total number of masked bases are:	13568722
  P.sojae	67593
  The number of bases masked by RepeatMasker:	22097217
  The number of bases masked by TransposonPSI:	5746805
  The total number of masked bases are:	23698005
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

```bash
  Pinf_ass=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
  Ppar_ass=assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra-310_2_supercontigs_parsed.fasta
  Pcap_ass=assembly/external_group/P.capsici/LT1534/dna/Phyca11_unmasked_genomic_scaffolds.fasta
  Psoj_ass=assembly/external_group/P.sojae/P6497/dna/Physo3_AssemblyScaffolds.fasta
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  for Genome in $Pinf_ass $Ppar_ass $Pcap_ass $Psoj_ass; do
  	qsub $ProgDir/sub_cegma.sh $Genome dna
  done
```

<!-- ## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.
 * The commands used to do this can be found in /gene_prediction/10300_braker1_prediction.md -->

## Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing
Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```bash
  Pinf_ass=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
  Ppar_ass=assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra-310_2_supercontigs_parsed.fasta
  Pcap_ass=assembly/external_group/P.capsici/LT1534/dna/Phyca11_unmasked_genomic_scaffolds.fasta
  Psoj_ass=assembly/external_group/P.sojae/P6497/dna/Physo3_AssemblyScaffolds.fasta

  for Genome in $Pinf_ass $Ppar_ass $Pcap_ass $Psoj_ass; do
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    qsub $ProgDir/run_ORF_finder.sh $Genome
  done
```

The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
for ORF_Gff in $(ls gene_pred/ORF_finder/P.*/*/*_ORF.gff | grep -v '_atg_'); do
Strain=$(echo $ORF_Gff | rev | cut -f2 -d '/' | rev)
Organism=$(echo $ORF_Gff | rev | cut -f3 -d '/' | rev)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
ORF_Gff_mod=gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF_corrected.gff3
$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
done
```

<!-- ```bash
  for ORF_Gff in $(ls analysis/rxlr_atg/P.infestans/T30-4/T30-4_ORF.gff); do
    Strain=$(echo $ORF_Gff | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $ORF_Gff | rev | cut -f3 -d '/' | rev)
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
    mkdir -p gene_pred/ORF_finder/$Organism/$Strain
    ORF_Gff_mod=gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF_corrected.gff3
    $ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
  done
``` -->



#Functional annotation

#Genomic analysis


## BLAST searches for previously identified effectors


BLast searches were used to identify the presence of previously identified
effectors in Phytophthora genomes.



```bash
  Pcac_ass=repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
  Pinf_ass=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
  Ppar_ass=assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra-310_2_supercontigs_parsed.fasta
  Pcap_ass=assembly/external_group/P.capsici/LT1534/dna/Phyca11_unmasked_genomic_scaffolds.fasta
  Psoj_ass=assembly/external_group/P.sojae/P6497/dna/Physo3_AssemblyScaffolds.fasta
  for Assembly in $(ls $Pcac_ass $Pinf_ass $Ppar_ass $Pcap_ass $Psoj_ass); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  	Effector=analysis/blast_homology/oomycete_avr_genes/appended_oomycete_avr_cds.fasta
    OutDir=analysis/blast_homology/$Organism/$Strain
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    qsub $ProgDir/run_blast2csv.sh $Effector dna $Assembly $OutDir
  done
```


Convert top blast hits into gff annotations

```bash
  for BlastHitsCsv in $(ls analysis/blast_homology/*/*/*appended_oomycete_avr_cds.fasta_hits.csv); do
    Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
    Column2=Avr_gene_homolog
    NumHits=1
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
  done
```

The blast hits were summarised in a single table for all the genomes. The top
identity of the top blast hit in relation to the enquire query sequence was
presented for each blast hit.

```bash
  OutFile=analysis/blast_homology/oomycete_avr_genes/oomycete_avr_summary.tab
  echo "Organism" > tmp2.tab
  cat analysis/blast_homology/P.sojae/P6497/P6497_appended_oomycete_avr_cds.fasta_hits.csv | cut -f1 >> tmp2.tab
  for BlastHits in $(ls analysis/blast_homology/*/*/*appended_oomycete_avr_cds.fasta_hits.csv); do
    Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
    echo "$Organism" > tmp.tab
    echo "$Strain" >> tmp.tab
    cat $BlastHits | cut -f9  | tail -n +2 >> tmp.tab
    paste tmp2.tab tmp.tab > $OutFile
    cp $OutFile tmp2.tab
  done
  rm tmp.tab
  rm tmp2.tab
```
<!--
```bash
	for HitsGff in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.gff | grep -v 'trinity' | grep -w 'Fus2'); do
		Strain=$(echo $HitsGff | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		# GffBraker=gene_pred/codingquary/$Organism/$Strain/final/final_genes_Braker.gff3
		# GffQuary=gene_pred/codingquary/$Organism/$Strain/final/final_genes_CodingQuary.gff3
		GffAppended=gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3
		OutDir=$(dirname $HitsGff)
		SixIntersect=$OutDir/"$Strain"_Fo_path_genes_CRX.fa_hit_genes.bed
		# bedtools intersect -wo -a $HitsGff -b $GffBraker > $SixIntersect
		# bedtools intersect -wo -a $HitsGff -b $GffQuary >> $SixIntersect
		bedtools intersect -wao -a $HitsGff -b $GffAppended > $SixIntersect
		bedtools intersect -wao -a $HitsGff -b $GffAppended | cut -f9,18 | grep -v 'Parent'
		echo ""
	done > analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX_hit_genes_summary.tab
``` -->




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


### A) From Published gene models - Signal peptide & RxLR motif

Required programs:
 * SigP
 * biopython


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
Pinf_pep=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
Ppar_pep=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins.pep.all.fa
Pcap_pep=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins.fasta
Psoj_pep=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401.aa.fasta

  for Proteome in $Pinf_pep $Ppar_pep $Pcap_pep $Psoj_pep; do
    echo "$Proteome"
    SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SplitDir=gene_pred/published_split/$Organism/$Strain
    mkdir -p $SplitDir
    BaseName="$Organism""_$Strain"_published_preds
    $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
    for File in $(ls $SplitDir/*_published_preds_*); do
      Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
      while [ $Jobs -gt 1 ]; do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
      done
      printf "\n"
      echo $File
      qsub $ProgDir/pred_sigP.sh $File
      # qsub $ProgDir/pred_sigP.sh $File signalp-4.1
    done
  done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
  for SplitDir in $(ls -d gene_pred/published_split/P.*/*); do
    Strain=$(echo $SplitDir | rev | cut -d '/' -f1 | rev)
    Organism=$(echo $SplitDir | rev | cut -d '/' -f2 | rev)
    InStringAA=''
    InStringNeg=''
    InStringTab=''
    InStringTxt=''
    for GRP in $(ls -l $SplitDir/*_published_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
      InStringAA="$InStringAA gene_pred/published_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_published_preds_$GRP""_sp.aa";  
      InStringNeg="$InStringNeg gene_pred/published_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_published_preds_$GRP""_sp_neg.aa";  
      InStringTab="$InStringTab gene_pred/published_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_published_preds_$GRP""_sp.tab";
      InStringTxt="$InStringTxt gene_pred/published_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_published_preds_$GRP""_sp.txt";
    done
    cat $InStringAA > gene_pred/published_sigP/$Organism/$Strain/"$Strain"_pub_sp.aa
    cat $InStringNeg > gene_pred/published_sigP/$Organism/$Strain/"$Strain"_pub_neg_sp.aa
    tail -n +2 -q $InStringTab > gene_pred/published_sigP/$Organism/$Strain/"$Strain"_pub_sp.tab
    cat $InStringTxt > gene_pred/published_sigP/$Organism/$Strain/"$Strain"_pub_sp.txt
  done
```

The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash

  for Secretome in $(ls gene_pred/published_sigP/*/P6497/*pub_sp.aa); do
    ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors;
    Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
    Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
    OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
    mkdir -p $OutDir;
    printf "\nstrain: $Strain\tspecies: $Organism\n";
    printf "the number of SigP gene is:\t";
    cat $Secretome | grep '>' | wc -l;
    printf "the number of SigP-RxLR genes are:\t";
    $ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_pub_RxLR_regex.fa;
    cat $OutDir/"$Strain"_pub_RxLR_regex.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' > $OutDir/"$Strain"_pub_RxLR_regex.txt
    cat $OutDir/"$Strain"_pub_RxLR_regex.txt | wc -l
    printf "the number of SigP-RxLR-EER genes are:\t";
    cat $OutDir/"$Strain"_pub_RxLR_regex.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' > $OutDir/"$Strain"_pub_RxLR_EER_regex.txt
    cat $OutDir/"$Strain"_pub_RxLR_EER_regex.txt | wc -l
    printf "\n"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    # $ProgDir/extract_from_fasta.py --fasta $OutDir/"$Strain"_pub_RxLR_regex.fa --headers $OutDir/"$Strain"_pub_RxLR_EER_regex.txt > $OutDir/"$Strain"_pub_RxLR_EER_regex.fa
    # GeneModels=$(ls assembly/external_group/P.*/$Strain/pep/*.gff*)
    # cat $GeneModels | grep -w -f $OutDir/"$Strain"_pub_RxLR_regex.txt > $OutDir/"$Strain"_pub_RxLR_regex.gff3
    # cat $GeneModels | grep -w -f $OutDir/"$Strain"_pub_RxLR_EER_regex.txt > $OutDir/"$Strain"_pub_RxLR_EER_regex.gff3
  done
```
```
  strain: LT1534	species: P.capsici
  the number of SigP gene is:	1650
  the number of SigP-RxLR genes are:	161
  the number of SigP-RxLR-EER genes are:	108


  strain: T30-4	species: P.infestans
  the number of SigP gene is:	2187
  the number of SigP-RxLR genes are:	510
  the number of SigP-RxLR-EER genes are:	357


  strain: 310	species: P.parisitica
  the number of SigP gene is:	2316
  the number of SigP-RxLR genes are:	312
  the number of SigP-RxLR-EER genes are:	206


  strain: P6497	species: P.sojae
  the number of SigP gene is:	2816
  the number of SigP-RxLR genes are:	416
  the number of SigP-RxLR-EER genes are:	240
```


### B) From published gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search gene models predicted with Braker1. These were run with the following commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
  HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
  for Proteome in $Pinf_pep $Ppar_pep $Pcap_pep $Psoj_pep; do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
    mkdir -p $OutDir
    HmmResults="$Strain"_pub_WY_hmmer.txt
    hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
    echo "$Organism $Strain"
    cat $OutDir/$HmmResults | grep 'Initial search space'
    cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
    HmmFasta="$Strain"_pub_WY_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
    Headers="$Strain"_pub_WY_hmmer_headers.txt
    cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' > $OutDir/$Headers
  # GeneModels=$(ls assembly/external_group/P.*/$Strain/pep/*.gff*)
  # cat $GeneModels | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_pub_WY_hmmer.gff3
  done
```
```
  P.infestans T30-4
  Initial search space (Z):              17787  [actual number of targets]
  Domain search space  (domZ):             267  [number of targets reported over threshold]
  P.parisitica 310
  Initial search space (Z):              20822  [actual number of targets]
  Domain search space  (domZ):             257  [number of targets reported over threshold]
  P.capsici LT1534
  Initial search space (Z):              19805  [actual number of targets]
  Domain search space  (domZ):             106  [number of targets reported over threshold]
  P.sojae P6497
  Initial search space (Z):              26584  [actual number of targets]
  Domain search space  (domZ):             285  [number of targets reported over threshold]
```

#### B.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
  Pinf_pep=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
  Ppar_pep=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins.pep.all.fa
  Pcap_pep=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins.fasta
  Psoj_pep=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401.aa.fasta

  cat $Pinf_pep | cut -f1 -d ' ' > assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26_parsed.pep.all.fa
  cat $Ppar_pep | cut -f1 -d ' ' > assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins_parsed.pep.all.fa
  cat $Pcap_pep | sed "s/jgi|Phyca11|.*|//g" | cut -f4 -d '|' > assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins_parsed.fasta
  cat $Psoj_pep | sed "s/jgi|Physo3||.*|//g" > assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401_parsed.aa.fasta

  Pinf_pep=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26_parsed.pep.all.fa
  Ppar_pep=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins_parsed.pep.all.fa
  Pcap_pep=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins_parsed.fasta
  Psoj_pep=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401_parsed.aa.fasta

  for Proteome in $Pinf_pep $Ppar_pep $Pcap_pep $Psoj_pep; do
    echo "$Proteome"
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/phobius/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    cat $Proteome | sed -e "s/*$//g" > $OutDir/tmp.fa
    phobius.pl $OutDir/tmp.fa > $OutDir/"$Strain"_phobius.txt
    cat $OutDir/"$Strain"_phobius.txt | grep -B1 'SIGNAL' | grep 'ID' | sed 's/ID   //g' > $OutDir/"$Strain"_phobius_headers.txt
    rm $OutDir/tmp.fa
    # qsub $ProgDir/sub_phobius.sh $Proteome $OutDir
  done
  cat analysis/phobius/P.infestans/T30-4/T30-4_phobius.txt | grep -B1 'SIGNAL' | grep 'ID' | sed 's/ID   //g' > analysis/phobius/P.infestans/T30-4/T30-4_phobius_headers.txt
  cat analysis/phobius/P.parisitica/310/310_phobius.txt | grep -B1 'SIGNAL' | grep 'ID' | sed 's/ID   //g' > analysis/phobius/P.parisitica/310/310_phobius_headers.txt
  cat analysis/phobius/P.capsici/LT1534/LT1534_phobius.txt | grep -B1 'SIGNAL' | grep 'ID' | sed 's/ID   //g' > analysis/phobius/P.capsici/LT1534/LT1534_phobius_headers.txt
  cat analysis/phobius/P.sojae/P6497/P6497_phobius.txt | grep -B1 'SIGNAL' | grep 'ID' | sed 's/ID   //g' | cut -f3 -d '|' > analysis/phobius/P.sojae/P6497/P6497_phobius_headers.txt
```

### C) From Augustus gene models - Hmm evidence of RxLR effectors
```bash
  for Proteome in $Pinf_pep $Ppar_pep $Pcap_pep $Psoj_pep; do
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
    HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
    mkdir -p $OutDir
    HmmResults="$Strain"_pub_RxLR_hmmer.txt
    hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
    echo "$Organism $Strain"
    cat $OutDir/$HmmResults | grep 'Initial search space'
    cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
    HmmFasta="$Strain"_pub_RxLR_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
    Headers="$Strain"_pub_RxLR_hmmer_headers.txt
    cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' > $OutDir/$Headers
		# # ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		# Col2=cropped.hmm
		# GeneModels=$(ls assembly/external_group/P.*/$Strain/pep/*.gff*)
		# # $ProgDir/gene_list_to_gff.pl $OutDir/$Headers $GeneModels $Col2 Name > $OutDir/"$Strain"_pub_RxLR_hmmer.gff3
		# cat $GeneModels | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_pub_RxLR_hmmer.gff3
	done
```
```
  P.infestans T30-4
  Initial search space (Z):              17787  [actual number of targets]
  Domain search space  (domZ):             290  [number of targets reported over threshold]
  P.parisitica 310
  Initial search space (Z):              20822  [actual number of targets]
  Domain search space  (domZ):             217  [number of targets reported over threshold]
  P.capsici LT1534
  Initial search space (Z):              19805  [actual number of targets]
  Domain search space  (domZ):              84  [number of targets reported over threshold]
  P.sojae P6497
  Initial search space (Z):              26584  [actual number of targets]
  Domain search space  (domZ):             280  [number of targets reported over threshold]
```

### D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:

<!-- ```bash
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
  HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
  for Proteome in $Pinf_pep $Ppar_pep $Pcap_pep $Psoj_pep; do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
    mkdir -p $OutDir
    HmmResults="$Strain"_pub_CRN_hmmer.txt
    hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
    echo "$Organism $Strain"
    cat $OutDir/$HmmResults | grep 'Initial search space'
    cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
    HmmFasta="$Strain"_pub_CRN_hmmer_out.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
    # Headers="$Strain"_pub_RxLR_hmmer_headers.txt
    # cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' > $OutDir/$Headers
    # GeneModels=$(ls assembly/external_group/P.*/$Strain/pep/*.gff*)
    # cat $GeneModels | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_pub_CRN_hmmer.gff3
  done
```
```
  P.infestans T30-4
  Initial search space (Z):              17787  [actual number of targets]
  Domain search space  (domZ):             187  [number of targets reported over threshold]
  P.parisitica 310
  Initial search space (Z):              20822  [actual number of targets]
  Domain search space  (domZ):              42  [number of targets reported over threshold]
  P.capsici LT1534
  Initial search space (Z):              19805  [actual number of targets]
  Domain search space  (domZ):             106  [number of targets reported over threshold]
  P.sojae P6497
  Initial search space (Z):              26584  [actual number of targets]
  Domain search space  (domZ):             152  [number of targets reported over threshold]
``` -->

```bash
  Pcac_pep=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
  Pinf_pep=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
  Ppar_pep=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins.pep.all.fa
  Pcap_pep=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins.fasta
  Psoj_pep=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401.aa.fasta
  for Proteome in $Pcac_pep $Pinf_pep $Ppar_pep $Pcap_pep $Psoj_pep; do
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
  P.cactorum - 10300
  Initial search space (Z):              20689  [actual number of targets]
  Domain search space  (domZ):              95  [number of targets reported over threshold]
  Initial search space (Z):              20689  [actual number of targets]
  Domain search space  (domZ):              82  [number of targets reported over threshold]
  67
  P.infestans - T30-4
  Initial search space (Z):              17787  [actual number of targets]
  Domain search space  (domZ):             193  [number of targets reported over threshold]
  Initial search space (Z):              17787  [actual number of targets]
  Domain search space  (domZ):             168  [number of targets reported over threshold]
  167
  P.parisitica - 310
  Initial search space (Z):              20822  [actual number of targets]
  Domain search space  (domZ):              44  [number of targets reported over threshold]
  Initial search space (Z):              20822  [actual number of targets]
  Domain search space  (domZ):              39  [number of targets reported over threshold]
  31
  P.capsici - LT1534
  Initial search space (Z):              19805  [actual number of targets]
  Domain search space  (domZ):             104  [number of targets reported over threshold]
  Initial search space (Z):              19805  [actual number of targets]
  Domain search space  (domZ):              88  [number of targets reported over threshold]
  82
  P.sojae - P6497
  Initial search space (Z):              26584  [actual number of targets]
  Domain search space  (domZ):             170  [number of targets reported over threshold]
  Initial search space (Z):              26584  [actual number of targets]
  Domain search space  (domZ):             133  [number of targets reported over threshold]
  108
```

Gff annotations were extracted for these putative crinklers:

```bash
  PinfPubGff=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts.gff3
  for GeneGff in $PinfPubGff; do
    echo "$GeneGff"
    Strain=$(echo "$GeneGff" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$GeneGff" | rev | cut -f4 -d '/' | rev)
    CrnDir=$(ls -d analysis/CRN_effectors/hmmer_CRN/$Species/$Strain)
    Source="pub"
    CRN_hmm_txt=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.txt
    CRN_hmm_gff=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.gff
    echo "$Species - $Strain"
    cat $GeneGff | grep -w -f $CRN_hmm_txt > $CRN_hmm_gff
    cat $CRN_hmm_gff | grep 'exon' | cut -f9 | cut -f2 -d ':' | sort | uniq | wc -l
  done
```

Gff annotations were extracted for these putative crinklers:

```bash
  PinfPubGff=assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts.gff3
  for GeneGff in $PinfPubGff; do
    echo "$GeneGff"
    Strain=$(echo "$GeneGff" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$GeneGff" | rev | cut -f4 -d '/' | rev)
    CrnDir=$(ls -d analysis/CRN_effectors/hmmer_CRN/$Species/$Strain)
    Source="pub"
    CRN_hmm_txt=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.txt
    CRN_hmm_gff=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.gff
    echo "$Species - $Strain"
    cat $GeneGff | grep -w -f $CRN_hmm_txt > $CRN_hmm_gff
    cat $CRN_hmm_gff | grep 'exon' | cut -f9 | cut -f2 -d ':' | sort | uniq | wc -l
  done
```

Fasta sequences for CRNs were extracted

```bash
  Pcac_pep=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
  Pinf_pep=assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
  Ppar_pep=assembly/external_group/P.parisitica/310/pep/phytophthora_parasitica_inra-310_2_proteins.pep.all.fa
  Pcap_pep=assembly/external_group/P.capsici/LT1534/pep/Phyca11_filtered_proteins.fasta
  Psoj_pep=assembly/external_group/P.sojae/P6497/pep/Physo3_GeneCatalog_proteins_20110401.aa.fasta
  for AugFa in $Pcac_pep $Pinf_pep $Ppar_pep $Pcap_pep $Psoj_pep; do
    Strain=$(echo "$AugFa" | rev | cut -f3 -d '/' | rev)
    Species=$(echo "$AugFa" | rev | cut -f4 -d '/' | rev)
    ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
    MergeDir=analysis/CRN_effectors/hmmer_CRN/$Species/$Strain
    TotalCRNsHeaders=$MergeDir/"$Strain"_Total_CRN_headers.txt
    CRNsFa=$MergeDir/"$Strain"_Total_CRN.fa
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
    $ProgDir/unwrap_fasta.py --inp_fasta $AugFa | grep -A1 -w -f $TotalCRNsHeaders | grep -v -E '^--$' > $CRNsFa
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalCRNsHeaders >> $CRNsFa
    echo "$Strain"
    echo "The number of sequences extracted is"
    cat $CRNsFa | grep '>' | wc -l
  done
```

Signal peptides were identified in the predicted Crinklers using SignalP2.0:

```bash
  for Crinklers in $(ls analysis/CRN_effectors/hmmer_CRN/*/*/*_Total_CRN.fa); do
    echo "$Crinklers"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    Strain=$(echo $Crinklers | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Crinklers | rev | cut -f3 -d '/' | rev)
    qsub $ProgDir/pred_sigP.sh $Crinklers
  done
```

```bash
  for File in $(ls gene_pred/hmmer_CRN_sigP/P*/*/*/*_Total_CRN.fa_sp.aa); do
    echo $(basename $File)
    cat $File | grep '>' | wc -l
  done
```

### E) From ORF gene models - Signal peptide & RxLR motif

Required programs:
 * SigP
 * biopython


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
  Pinf_ORF=gene_pred/ORF_finder/P.infestans/T30-4/T30-4.aa_cat.fa
  Ppar_ORF=gene_pred/ORF_finder/P.parisitica/310/310.aa_cat.fa
  Pcap_ORF=gene_pred/ORF_finder/P.capsici/LT1534/LT1534.aa_cat.fa
  Psoj_ORF=gene_pred/ORF_finder/P.sojae/P6497/P6497.aa_cat.fa
  for Proteome in $Pinf_ORF $Ppar_ORF $Pcap_ORF $Psoj_ORF; do
  # for Proteome in $Psoj_ORF; do
    echo "$Proteome"
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
      while [ $Jobs -gt 1 ]; do
        sleep 5
        printf "."
        Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
      done
      printf "\n"
      echo $File
      qsub $ProgDir/pred_sigP.sh $File
      # qsub $ProgDir/pred_sigP.sh $File signalp-4.1
    done
  done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
  for SplitDir in $(ls -d gene_pred/ORF_split/P.*/*); do
    # for SplitDir in $(ls -d gene_pred/ORF_split/P.*/P6497); do
  	Strain=$(echo $SplitDir | rev | cut -d '/' -f1 | rev)
  	Organism=$(echo $SplitDir | rev | cut -d '/' -f2 | rev)
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
  done
```

Names of ORFs containing signal peptides were extracted from fasta files. This
included information on the position and hmm score of RxLRs.

```bash
  for FastaFile in $(ls gene_pred/ORF_sigP/*/*/*_ORF_sp.aa); do
    Strain=$(echo $FastaFile | rev | cut -d '/' -f2 | rev)
    Organism=$(echo $FastaFile | rev | cut -d '/' -f3 | rev)
    echo "$Strain"
    SigP_headers=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_names.txt
    cat $FastaFile | grep '>' | sed -r 's/>//g' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | sed 's/--//g' > $SigP_headers
  done
```
```bash
  for FastaFile in $(ls gene_pred/ORF_sigP/*/*/*_ORF_sp.aa | grep -e 'infestans' -e 'parisitica'); do
    Strain=$(echo $FastaFile | rev | cut -d '/' -f2 | rev)
    Organism=$(echo $FastaFile | rev | cut -d '/' -f3 | rev)
    echo "$Strain"
    SigP_headers=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_names.txt
    cat $FastaFile | grep '>' | sed -r 's/>//g' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | sed 's/--//g' | sed 's/:/_a_/g' | sed 's/supercont1./supercont1_b_/g' | sed 's/Supercontig_2./Supercontig_2_c_/g' > $SigP_headers
  done
```

Due to the nature of predicting ORFs, some features overlapped with one another.
A single ORF was selected from each set of overlapped ORFs. This was was
selected on the basis of its SignalP Hmm score. Biopython was used to identify
overlapps and identify the ORF with the best signalP score.

```bash
  for SigP_fasta in $(ls gene_pred/ORF_sigP/P.*/*/*_ORF_sp.aa | grep -v -e 'infestans' -e 'parisitica'); do
    # for SigP_fasta in $(ls gene_pred/ORF_sigP/P.*/P6497/*_ORF_sp.aa); do
    Strain=$(echo $SigP_fasta | rev | cut -d '/' -f2 | rev)
    Organism=$(echo $SigP_fasta | rev | cut -d '/' -f3 | rev)
    echo "$Strain"
    ORF_Gff=gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF_corrected.gff3
    SigP_fasta=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp.aa
    SigP_headers=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_names.txt
    SigP_Gff=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_unmerged.gff
    SigP_Merged_Gff=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.gff
    SigP_Merged_txt=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.txt
    SigP_Merged_AA=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.aa

    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $SigP_headers $ORF_Gff SigP Name > $SigP_Gff
    ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $SigP_Gff --db sigP_ORF.db
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF.db --id sigP_ORF --out sigP_ORF_merged.db --gff > $SigP_Merged_Gff
    cat $SigP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d'=' | rev > $SigP_Merged_txt
    $ProgDir/extract_from_fasta.py --fasta $SigP_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
  done
```

```bash
  for SigP_fasta in $(ls gene_pred/ORF_sigP/P.*/*/*_ORF_sp.aa | grep -e 'infestans' -e 'parisitica'); do
    Strain=$(echo $SigP_fasta | rev | cut -d '/' -f2 | rev)
    Organism=$(echo $SigP_fasta | rev | cut -d '/' -f3 | rev)
    echo "$Strain"
    ORF_Gff=gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF_corrected.gff3
    SigP_fasta=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp.aa
    SigP_headers=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_names.txt
    SigP_Gff=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_unmerged.gff
    SigP_Merged_Gff=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.gff
    SigP_Merged_txt=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.txt
    SigP_Merged_AA=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_merged.aa

    cat $ORF_Gff | sed 's/:/_a_/g' | sed 's/supercont1./supercont1_b_/g' | sed 's/Supercontig_2./Supercontig_2_c_/g' > tmp.gff
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $SigP_headers tmp.gff SigP Name > $SigP_Gff
    ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $SigP_Gff --db sigP_ORF.db
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF.db --id sigP_ORF --out sigP_ORF_merged.db --gff | sed 's/_a_/:/g' | sed 's/supercont1_b_/supercont1./g' | sed 's/Supercontig_2_c_/Supercontig_2./g' > $SigP_Merged_Gff
    cat $SigP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d'=' | rev > $SigP_Merged_txt
    $ProgDir/extract_from_fasta.py --fasta $SigP_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
  done
```



The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
  for Secretome in $(ls gene_pred/ORF_sigP/P.*/*/*_ORF_sp_merged.aa); do
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

```
strain: 10300	species: P.cactorum
the number of SigP gene is:	14767
the number of SigP-RxLR genes are:	916
the number of SigP-RxLR-EER genes are:	169


strain: LT1534	species: P.capsici
the number of SigP gene is:	14007
the number of SigP-RxLR genes are:	964
the number of SigP-RxLR-EER genes are:	244


strain: T30-4	species: P.infestans
the number of SigP gene is:	36713
the number of SigP-RxLR genes are:	2786
the number of SigP-RxLR-EER genes are:	408


strain: 310	species: P.parisitica
the number of SigP gene is:	13652
the number of SigP-RxLR genes are:	915
the number of SigP-RxLR-EER genes are:	288


strain: P6497	species: P.sojae
the number of SigP gene is:	24495
the number of SigP-RxLR genes are:	1742
the number of SigP-RxLR-EER genes are:	277
```

#### E.2) SigP ORF Prediction using Phobius

Secreted proteins were also predicted using Phobius. This was done in a screen
session on the head node of the cluster.

```bash
for FastaFile in $(ls gene_pred/ORF_sigP/*/*/*_ORF_sp.aa | grep -e 'P.infestans' -e 'P.parisitica' -e 'P.capsici' -e 'P.sojae'); do
Strain=$(echo $FastaFile | rev | cut -d '/' -f2 | rev)
Organism=$(echo $FastaFile | rev | cut -d '/' -f3 | rev)
echo "$Strain"
OutDir=analysis/phobius/$Organism/$Strain
mkdir -p $OutDir
phobius.pl $FastaFile > $OutDir/"$Strain"_phobius_ORF.txt
cat $OutDir/"$Strain"_phobius_ORF.txt | grep -B1 'SIGNAL' | grep 'ID' | sed s'/ID   //g' > $OutDir/"$Strain"_phobius_headers_ORF.txt
done
```

### F) From ORF gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search ORFs predicted with atg.pl. These were run with the following commands:


```bash
  for Secretome in $(ls gene_pred/ORF_sigP/P.*/*/*_ORF_sp_merged.aa); do
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

```
  P.cactorum 10300
  Initial search space (Z):              14767  [actual number of targets]
  Domain search space  (domZ):             113  [number of targets reported over threshold]
  P.capsici LT1534
  Initial search space (Z):              14007  [actual number of targets]
  Domain search space  (domZ):             109  [number of targets reported over threshold]
  P.infestans T30-4
  Initial search space (Z):              36713  [actual number of targets]
  Domain search space  (domZ):             239  [number of targets reported over threshold]
  P.parisitica 310
  Initial search space (Z):              13652  [actual number of targets]
  Domain search space  (domZ):             163  [number of targets reported over threshold]
  P.sojae P6497
  Initial search space (Z):              24495  [actual number of targets]
  Domain search space  (domZ):             158  [number of targets reported over threshold]
```


### G) From ORF gene models - Hmm evidence of RxLR effectors

```bash
  for Secretome in $(ls gene_pred/ORF_sigP/P.*/*/*_ORF_sp_merged.aa); do
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
```
  P.cactorum 10300
  Initial search space (Z):              14767  [actual number of targets]
  Domain search space  (domZ):             144  [number of targets reported over threshold]
  P.capsici LT1534
  Initial search space (Z):              14007  [actual number of targets]
  Domain search space  (domZ):             158  [number of targets reported over threshold]
  P.infestans T30-4
  Initial search space (Z):              36713  [actual number of targets]
  Domain search space  (domZ):             280  [number of targets reported over threshold]
  P.parisitica 310
  Initial search space (Z):              13652  [actual number of targets]
  Domain search space  (domZ):             233  [number of targets reported over threshold]
  P.sojae P6497
  Initial search space (Z):              24495  [actual number of targets]
  Domain search space  (domZ):             227  [number of targets reported over threshold]
```


### H) From ORF gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in ORF gene models. This was done with the following commands:

```bash
  for Proteome in $(ls gene_pred/ORF_finder/*/*/*.aa_cat.fa); do
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
    Headers="$Strain"_CRN_hmmer_unmerged_headers.txt
    cat $OutDir/$HmmFastaDWL | grep '>' | grep -w -f $CommonHeaders | tr -d '>' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | tr -d '-' | sed 's/hmm_score/HMM_score/g' > $OutDir/$Headers
    # As we are dealing with JGI and Broad sequences, some headers need formatting:
    cat $OutDir/$Headers | sed 's/:/_a_/g' | sed 's/supercont1./supercont1_b_/g' | sed 's/Supercontig_2./Supercontig_c_/g' > tmp.txt
    # As we are dealing with JGI and Broad sequences, some features need formatting:
    ORF_Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/*_ORF_corrected.gff3)
    cat $ORF_Gff | sed 's/:/_a_/g' | sed 's/supercont1./supercont1_b_/g' | sed 's/Supercontig_2./Supercontig_c_/g' > tmp.gff
    # Gff features were extracted for each header
    CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl tmp.txt tmp.gff CRN_HMM Name > $CRN_unmerged_Gff
    # Gff features were merged based upon the DWL hmm score
    DbDir=analysis/databases/$Organism/$Strain
    mkdir -p $DbDir
    ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $CRN_unmerged_Gff --db $DbDir/CRN_ORF.db
    CRN_Merged_Gff=$OutDir/"$Strain"_CRN_merged_hmmer.gff3
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp $DbDir/CRN_ORF.db --id LFLAK_DWL_CRN --out $DbDir/CRN_ORF_merged.db --gff > $CRN_Merged_Gff
    # As we are dealing with JGI and Broad sequences, some features need formatting:
    sed -i 's/_a_/:/g' $CRN_Merged_Gff
    sed -i 's/supercont1_b_/supercont1./g' $CRN_Merged_Gff
    sed -i 's/Supercontig_c_/Supercontig_2./g' $CRN_Merged_Gff
    # Final results are reported:
    echo "Number of CRN ORFs after merging:"
    cat $CRN_Merged_Gff | grep 'gene' | wc -l
    # Temporary files containing the reformatted headers and features were deleted
    rm tmp.txt
    rm tmp.gff
  done
```
<!--
```
  P.cactorum 10300
  Initial search space (Z):             443642  [actual number of targets]
  Domain search space  (domZ):             206  [number of targets reported over threshold]
  Number of CRN ORFs after merging:
  101
  P.capsici LT1534
  Initial search space (Z):             437865  [actual number of targets]
  Domain search space  (domZ):             307  [number of targets reported over threshold]
  Number of CRN ORFs after merging:
  159
  P.infestans T30-4
  Initial search space (Z):            1363381  [actual number of targets]
  Domain search space  (domZ):             586  [number of targets reported over threshold]
  Number of CRN ORFs after merging:
  367
  P.parisitica 310
  Initial search space (Z):             410940  [actual number of targets]
  Domain search space  (domZ):             139  [number of targets reported over threshold]
  Number of CRN ORFs after merging:
  53
  P.sojae P6497
  Initial search space (Z):             696564  [actual number of targets]
  Domain search space  (domZ):             398  [number of targets reported over threshold]
  Number of CRN ORFs after merging:
  212
``` -->
```
Searching for LFLAK domains in: P.cactorum 10300
Initial search space (Z):             443642  [actual number of targets]
Domain search space  (domZ):             183  [number of targets reported over threshold]
Searching for DWL domains in: P.cactorum 10300
Initial search space (Z):             443642  [actual number of targets]
Domain search space  (domZ):             205  [number of targets reported over threshold]
The number of CRNs common to both models are:
85
Number of CRN ORFs after merging:
58
Searching for LFLAK domains in: P.capsici LT1534
Initial search space (Z):             437865  [actual number of targets]
Domain search space  (domZ):             296  [number of targets reported over threshold]
Searching for DWL domains in: P.capsici LT1534
Initial search space (Z):             437865  [actual number of targets]
Domain search space  (domZ):             361  [number of targets reported over threshold]
The number of CRNs common to both models are:
174
Number of CRN ORFs after merging:
104
Searching for LFLAK domains in: P.infestans T30-4
Initial search space (Z):            1363381  [actual number of targets]
Domain search space  (domZ):             568  [number of targets reported over threshold]
Searching for DWL domains in: P.infestans T30-4
Initial search space (Z):            1363381  [actual number of targets]
Domain search space  (domZ):             760  [number of targets reported over threshold]
The number of CRNs common to both models are:
373
Number of CRN ORFs after merging:
255
Searching for LFLAK domains in: P.parisitica 310
Initial search space (Z):             410940  [actual number of targets]
Domain search space  (domZ):             104  [number of targets reported over threshold]
Searching for DWL domains in: P.parisitica 310
Initial search space (Z):             410940  [actual number of targets]
Domain search space  (domZ):             105  [number of targets reported over threshold]
The number of CRNs common to both models are:
47
Number of CRN ORFs after merging:
26
Searching for LFLAK domains in: P.sojae P6497
Initial search space (Z):             696564  [actual number of targets]
Domain search space  (domZ):             414  [number of targets reported over threshold]
Searching for DWL domains in: P.sojae P6497
Initial search space (Z):             696564  [actual number of targets]
Domain search space  (domZ):             506  [number of targets reported over threshold]
The number of CRNs common to both models are:
230
Number of CRN ORFs after merging:
147
```
