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

```shell
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

```shell
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

```shell
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

```shell
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

```bash
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
  # for Proteome in $Pinf_ORF $Ppar_ORF $Pcap_ORF $Psoj_ORF; do
  for Proteome in $Psoj_ORF; do
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
<-- P. soj
The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
  # for SplitDir in $(ls -d gene_pred/ORF_split/P.*/*); do
    for SplitDir in $(ls -d gene_pred/ORF_split/P.*/P6497); do
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
  for FastaFile in $(ls gene_pred/ORF_sigP/*/*/*_ORF_sp.aa | grep -e 'sojae' -e 'infestans' -e 'parisitica'); do
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
for SigP_fasta in $(ls gene_pred/ORF_sigP/P.*/*/*_ORF_sp.aa | grep -v -e 'sojae' -e 'infestans' -e 'parisitica'); do
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
  for SigP_fasta in $(ls gene_pred/ORF_sigP/P.*/310/*_ORF_sp.aa | grep -e 'sojae' -e 'infestans' -e 'parisitica'); do
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
  the number of SigP gene is:	15271
  the number of SigP-RxLR genes are:	935
  the number of SigP-RxLR-EER genes are:	170


  strain: LT1534	species: P.capsici
  the number of SigP gene is:	14865
  the number of SigP-RxLR genes are:	1004
  the number of SigP-RxLR-EER genes are:	254


  strain: T30-4	species: P.infestans
  the number of SigP gene is:	37637
  the number of SigP-RxLR genes are:	2825
  the number of SigP-RxLR-EER genes are:	417


  strain: 310	species: P.parisitica
  the number of SigP gene is:	14526
  the number of SigP-RxLR genes are:	949
  the number of SigP-RxLR-EER genes are:	309


  strain: 67593	species: P.sojae
  the number of SigP gene is:	24222
  the number of SigP-RxLR genes are:	1727
  the number of SigP-RxLR-EER genes are:	282
```


 ### F) From ORF gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search ORFs predicted with atg.pl. These were run with the following commands:


```bash
  for Secretome in $(ls gene_pred/ORF_sigP/P.*/310/*_ORF_sp_merged.aa); do
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
  Initial search space (Z):              15271  [actual number of targets]
  Domain search space  (domZ):             113  [number of targets reported over threshold]
  P.capsici LT1534
  Initial search space (Z):              14865  [actual number of targets]
  Domain search space  (domZ):             117  [number of targets reported over threshold]
  P.infestans T30-4
  Initial search space (Z):              37637  [actual number of targets]
  Domain search space  (domZ):             244  [number of targets reported over threshold]
  P.parisitica 310
  Initial search space (Z):              14526  [actual number of targets]
  Domain search space  (domZ):             175  [number of targets reported over threshold]
  P.sojae 67593
  Initial search space (Z):              24222  [actual number of targets]
  Domain search space  (domZ):             159  [number of targets reported over threshold]
```


### G) From ORF gene models - Hmm evidence of RxLR effectors

```bash
  for Secretome in $(ls gene_pred/ORF_sigP/P.*/310/*_ORF_sp_merged.aa); do
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
  Initial search space (Z):              15271  [actual number of targets]
  Domain search space  (domZ):             145  [number of targets reported over threshold]
  P.capsici LT1534
  Initial search space (Z):              14865  [actual number of targets]
  Domain search space  (domZ):             165  [number of targets reported over threshold]
  P.infestans T30-4
  Initial search space (Z):              37637  [actual number of targets]
  Domain search space  (domZ):             284  [number of targets reported over threshold]
  P.parisitica 310
  Initial search space (Z):              14526  [actual number of targets]
  Domain search space  (domZ):             244  [number of targets reported over threshold]
  P.sojae 67593
  Initial search space (Z):              24222  [actual number of targets]
  Domain search space  (domZ):             228  [number of targets reported over threshold]
```

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

<!-- ```bash
  for Secretome in $(ls gene_pred/ORF_sigP/P.*/*/*_ORF_sp_merged.aa); do
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
```
  P.cactorum 10300
  Initial search space (Z):              15271  [actual number of targets]
  Domain search space  (domZ):              18  [number of targets reported over threshold]
  P.capsici LT1534
  Initial search space (Z):              14865  [actual number of targets]
  Domain search space  (domZ):              44  [number of targets reported over threshold]
  P.infestans T30-4
  Initial search space (Z):              26646  [actual number of targets]
  Domain search space  (domZ):               4  [number of targets reported over threshold]
  P.parisitica 310
  Initial search space (Z):              14526  [actual number of targets]
  Domain search space  (domZ):              14  [number of targets reported over threshold]
  P.sojae 67593
  Initial search space (Z):              24222  [actual number of targets]
  Domain search space  (domZ):              75  [number of targets reported over threshold]
``` -->

```bash
  for Proteome in $(ls gene_pred/ORF_finder/*/*/*.aa_cat.fa); do
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
    HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
    mkdir -p $OutDir
    HmmResults="$Strain"_ORF_CRN_unmerged_hmmer.txt
    hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
    echo "$Organism $Strain"
    cat $OutDir/$HmmResults | grep 'Initial search space'
    cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
    HmmFasta="$Strain"_ORF_CRN_hmmer_unmerged_out.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
    Headers="$Strain"_CRN_hmmer_unmerged_headers.txt
    cat $OutDir/$HmmFasta | grep '>' | tr -d '>' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | tr -d '-' | sed 's/hmm_score/HMM_score/g' > $OutDir/$Headers
    cat $OutDir/$Headers | sed 's/:/_a_/g' | sed 's/supercont1./supercont1_b_/g' | sed 's/Supercontig_2./Supercontig_c_/g' > tmp.txt
    ORF_Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/*_ORF_corrected.gff3)
    cat $ORF_Gff | sed 's/:/_a_/g' | sed 's/supercont1./supercont1_b_/g' | sed 's/Supercontig_2./Supercontig_c_/g' > tmp.gff
    CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl tmp.txt tmp.gff CRN_HMM Name > $CRN_unmerged_Gff
    DbDir=analysis/databases/$Organism/$Strain
    mkdir -p $DbDir
    ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $CRN_unmerged_Gff --db $DbDir/CRN_ORF.db
    CRN_Merged_Gff=$OutDir/"$Strain"_CRN_merged_hmmer.gff3
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp $DbDir/CRN_ORF.db --id $HmmModel --out $DbDir/CRN_ORF_merged.db --gff > $CRN_Merged_Gff
    sed -i 's/_a_/:/g' $CRN_Merged_Gff
    sed -i 's/supercont1_b_/supercont1./g' $CRN_Merged_Gff
    sed -i 's/Supercontig_c_/Supercontig_2./g' $CRN_Merged_Gff
    echo "Number of CRN ORFs after merging:"
    cat $CRN_Merged_Gff | grep 'gene' | wc -l
    rm tmp.txt
    rm tmp.gff
  done
```

```
  P.cactorum 10300
  Initial search space (Z):             459307  [actual number of targets]
  Domain search space  (domZ):             225  [number of targets reported over threshold]
  Number of CRN ORFs after merging:
  115
  P.capsici LT1534
  Initial search space (Z):             478905  [actual number of targets]
  Domain search space  (domZ):             356  [number of targets reported over threshold]
  Number of CRN ORFs after merging:
  181
  P.infestans T30-4
  Initial search space (Z):            1415148  [actual number of targets]
  Domain search space  (domZ):             618  [number of targets reported over threshold]
  Number of CRN ORFs after merging:
  384
  P.parisitica 310
  Initial search space (Z):             442892  [actual number of targets]
  Domain search space  (domZ):             151  [number of targets reported over threshold]
  Number of CRN ORFs after merging:
  59
  P.sojae 67593
  Initial search space (Z):             689510  [actual number of targets]
  Domain search space  (domZ):             384  [number of targets reported over threshold]
  Number of CRN ORFs after merging:
  203
```


<!--
## 4. 2 Ananlysis of RxLR effectors

Due to RxLR effectors being predicted from a number of sources the number of
unique RxLRs were identified from motif and Hmm searches within gene models.
221 RxLR effectors were predicted in total from the P.cactorum genome. Of these,
90 were shared between both datasets.

```bash
  for InDir in $(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/*/310); do
    Strain=$(echo "$InDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$InDir" | rev | cut -f2 -d '/' | rev)
    RxLR_motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Species/$Strain/"$Strain"*_RxLR_EER_regex.fa | grep -v 'ORF')
    RxLR_hmm=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Species/$Strain/"$Strain"*_RxLR_hmmer.fa | grep -v 'ORF')
    echo "$Species - $Strain"
    echo "Total number of RxLRs in predicted genes:"
    cat $RxLR_motif $RxLR_hmm | grep '>' | cut -f1 | sort | uniq | wc -l
    echo "Total number of RxLRs shared between prediciton sources:"
    cat $RxLR_motif $RxLR_hmm | grep '>' | cut -f1 | sort | uniq -d | wc -l
    echo ""
  done
```

```
  P.cactorum - 10300
  Total number of RxLRs in predicted genes:
  291
  Total number of RxLRs shared between prediciton sources:
  90

  P.capsici - LT1534
  Total number of RxLRs in predicted genes:
  245
  Total number of RxLRs shared between prediciton sources:
  0

  P.infestans - T30-4
  Total number of RxLRs in predicted genes:
  870
  Total number of RxLRs shared between prediciton sources:
  70

  P.parisitica - 310
  Total number of RxLRs in predicted genes:
  609
  Total number of RxLRs shared between prediciton sources:
  80

  P.sojae - 67593
  Total number of RxLRs in predicted genes:
  531
  Total number of RxLRs shared between prediciton sources:
  173
``` -->
