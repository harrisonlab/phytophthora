

```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/assembly
qsub $ProgDir/spades_10300_assembly.sh
```



Assembly results were summarised using Quast:

```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs*/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep '10300'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    OutDir=$(dirname $Assembly)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir;
  done
```
<!--
```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for Exclude_db in "bacillus" "delftia" "Paen"; do
      Good_db="phytoph"
      AssemblyDir=$(dirname $Assembly)
      OutDir=$AssemblyDir/../deconseq_$Exclude_db
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
      qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
    done
  done
```

Results were summarised using the commands:

```bash
# for File in $(ls assembly/spades/P.*/*/deconseq/log.txt); do
# for File in $(ls assembly/spades/P.*/*/deconseq_common/log.txt); do
for Exclude_db in "bacillus" "delftia" "Paen"; do
  echo $Exclude_db
  for File in $(ls assembly/spades/P.*/*/*/log.txt | grep "$Exclude_db"); do
    Name=$(echo $File | rev | cut -f3 -d '/' | rev);
    Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
    Both=$(cat $File |cut -f2 | head -n2 | tail -n1);
    Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
    printf "$Name\t$Good\t$Both\t$Bad\n";
  done
done
```
-->



Numbers of busco genes in each assembly were identified:

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs*/contigs_min_500bp.fasta | grep -e 'P.cactorum' -e 'P.idaei' | grep '10300'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB="Eukaryotic"
# OutDir=gene_pred/busco/$Organism/$Strain/assembly
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

The Spades genome assembly was considered to be worse than the abyss assembly




The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_hardmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```
Number of masked bases:
10813641
```

## Gene prediction

#### Aligning


```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for FileF in $(ls qc_rna/paired/genbank/*/F/*.fastq); do
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
done
printf "\n"
FileR=$(echo $FileF | sed 's&/F/&/R/&g' | sed 's/F.fastq/R.fastq/g')
echo $FileF
echo $FileR
Prefix=$(echo $FileF | rev | cut -f3 -d '/' | rev)
# Timepoint=$(echo $FileF | rev | cut -f2 -d '/' | rev)
Timepoint="genbank"
#echo "$Timepoint"
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
done
```

Accepted hits .bam file were concatenated and indexed for use for gene model training:


```bash
for OutDir in $(ls -d alignment/star/*/* | grep '1177'); do
  Strain=$(echo $OutDir | rev | cut -d '/' -f1 | rev)
  Organism=$(echo $OutDir | rev | cut -d '/' -f2 | rev)
  echo "$Organism - $Strain"
  # For all alignments
  BamFiles=$(ls $OutDir/treatment/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
  mkdir -p $OutDir/concatenated
  samtools merge -f $OutDir/concatenated/concatenated.bam $BamFiles
done
```


#### Braker prediction

```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
mkdir -p alignment/$Organism/$Strain/concatenated
OutDir=gene_pred/braker/$Organism/"$Strain"_braker
AcceptedHits=$(ls alignment/star/$Organism/$Strain/genbank/P.cactorum/star_aligmentAligned.sortedByCoord.out.bam)
GeneModelName="$Organism"_"$Strain"_braker
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

** Number of genes predicted:  **


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=$(ls alignment/star/$Organism/$Strain/genbank/P.cactorum/star_aligmentAligned.sortedByCoord.out.bam)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```


Secondly, genes were predicted using CodingQuary:

```bash
  for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain
    CufflinksGTF=$(ls gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf)
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
  done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
  for BrakerGff in $(ls gene_pred/braker/*/*_braker/*/augustus.gff3 | grep '10300'); do
    Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
    Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_softmasked_repeatmasker_TPSI_appended.fa  | grep '10300_abyss_53_repmask')
    CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
    PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
    AddDir=gene_pred/codingquary/$Organism/$Strain/additional
    FinalDir=gene_pred/final/$Organism/$Strain/final
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
    # -
    # This section is edited
    $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuary_unspliced.gff3
    $ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuary_unspliced.gff3 > $FinalDir/final_genes_CodingQuary.gff3
    # -
    $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
    cp $BrakerGff $FinalDir/final_genes_Braker.gff3
    $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
    cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
    cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
    cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
    cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta


    GffBraker=$FinalDir/final_genes_Braker.gff3
    GffQuary=$FinalDir/final_genes_CodingQuary.gff3
    GffAppended=$FinalDir/final_genes_appended.gff3
    cat $GffBraker $GffQuary > $GffAppended
  done
```

In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified


```bash
  for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep '10300'); do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    FinalDir=gene_pred/final/$Organism/$Strain/final
    GffFiltered=$FinalDir/filtered_duplicates.gff
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
    GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
    LogFile=$FinalDir/final_genes_appended_renamed.log
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered
    Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_softmasked_repeatmasker_TPSI_appended.fa  | grep '10300_abyss_53_repmask')
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed
    # The proteins fasta file contains * instead of Xs for stop codons, these should
    # be changed
    sed -i 's/\*/X/g' gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta
  done
```


## Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing
Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
    echo "$Assembly"
  	qsub $ProgDir/run_ORF_finder.sh $Assembly
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
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta | grep '10300'); do
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/final_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_final
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_final_*); do
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
while [ $Jobs -gt 10 ]; do
sleep 1
printf "."
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
done
printf "\n"
echo $File
qsub $ProgDir/pred_sigP.sh $File
qsub $ProgDir/pred_sigP.sh $File signalp-3.0
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
 ```

 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

 ```bash
  for SplitDir in $(ls -d gene_pred/final_split/P.*/* | grep '10300'); do
    Strain=$(echo $SplitDir | cut -d '/' -f4)
    Organism=$(echo $SplitDir | cut -d '/' -f3)
    echo "$Organism - $Strain"
    InStringAA=''
    InStringNeg=''
    InStringTab=''
    InStringTxt=''
    for SigpDir in $(ls -d gene_pred/final_sig* | cut -f2 -d'/'); do
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
  done
```

#### B.2) Prediction using Phobius

 Secreted proteins were also predicted using Phobius

 ```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta | grep '10300'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/phobius/$Organism/$Strain
    mkdir -p $OutDir
    phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    $ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
  done
 ```


Secreted proteins from different sources were combined into a single file:

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta | grep '10300'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/combined_sigP/$Organism/$Strain
    mkdir -p $OutDir
    echo "The following number of sequences were predicted as secreted:"
    cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
    cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | wc -l
    echo "This represented the following number of unique genes:"
    cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
    cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
  done
```

```

```


### C) From Augustus gene models - Effector-like structure identification using EffectorP

Required programs:
 * EffectorP.py

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta | grep '10300'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
BaseName="$Organism"_"$Strain"_EffectorP
OutDir=analysis/effectorP/$Organism/$Strain
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
done
```

Those genes that were predicted as secreted and tested positive by effectorP
were identified:

```bash
  for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
    cat $File | grep 'Effector' | cut -f1 > $Headers
    printf "EffectorP headers:\t"
    cat $Headers | wc -l
    Secretome=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_secreted.fa)
    OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
    OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
    printf "Secreted effectorP headers:\t"
    cat $OutFileHeaders | wc -l
    Gff=$(ls gene_pred/codingquary/$Organism/$Strain/*/final_genes_appended.gff3)
    EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
  done
```
```

```

## D) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta | grep '10300'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/CAZY/$Organism/$Strain
mkdir -p $OutDir
Prefix="$Strain"_CAZY
CazyHmm=../../dbCAN/dbCAN-fam-HMMs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $File)
echo "$Organism - $Strain"
ProgDir=/home/groups/harrisonlab/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
printf "number of CAZY genes identified:\t"
cat $CazyHeaders | wc -l
Gff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
CazyGff=$OutDir/"$Strain"_CAZY.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

SecretedProts=$(ls gene_pred/combined_sigP/$Organism/$Strain/*_all_secreted.fa)
SecretedHeaders=$(echo $SecretedProts | sed 's/.fa/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
printf "number of Secreted CAZY genes identified:\t"
cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
done
```

```

```

Note - the CAZY genes identified may need further filtering based on e value and
cuttoff length - see below:

Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage

* For fungi, use E-value < 1e-17 and coverage > 0.45

* The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)


# D) Prediction of RxLRs from - Augustus/CodingQuary gene models

 The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

 The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash

for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep -w -e '414' -e 'P.cactorum' -e 'P.idaei' | grep -v '10300'); do
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
cat $Secretome | grep '>' | cut -f1 | tr -d ' '| sort | uniq | wc -l;

printf "the number of SigP-RxLR genes are:\t";
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_all_secreted_RxLR_regex.fa;
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $OutDir/"$Strain"_RxLR_regex.txt
cat $OutDir/"$Strain"_RxLR_regex.txt | wc -l

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa


printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $OutDir/"$Strain"_RxLR_EER_regex.txt
cat $OutDir/"$Strain"_RxLR_EER_regex.txt | wc -l
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa

printf "\n"

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_regex.txt
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_EER_regex.txt

cat $Gff | grep -w -f $OutDir/"$Strain"_RxLR_regex.txt> $OutDir/"$Strain"_RxLR_regex.gff3
cat $Gff | grep -w -f $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.gff3
done
```

```

```


### G) From Secreted gene models - Hmm evidence of RxLR effectors

```bash
  for Proteome in $(ls gene_pred/codingquary/*/*/*/final_genes_combined.pep.fasta | grep -w -e '414' -e 'idaei'); do
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

```


### F) Combining RxLRs from Regex and hmm searches


The total RxLRs are

```bash
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_RxLR_EER_regex.txt | grep -v -e 'Aug' -e '10300' | grep -e 'P.idaei' -e 'P.cactorum' | grep -v '414'); do
Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
Proteome=$(ls gene_pred/codingquary/$Organism/$Strain/*/final_genes_combined.pep.fasta)
HmmRxLR=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer_headers.txt
echo "$Organism - $Strain"
echo "Number of RxLRs identified by Regex:"
cat $RegexRxLR | sort | uniq | wc -l
echo "Number of RxLRs identified by Hmm:"
cat $HmmRxLR | sort | uniq | wc -l
echo "Number of RxLRs in combined dataset:"
cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
# echo "Number of RxLRs in both datasets:"
# cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
echo ""
# echo "Extracting RxLRs from datasets"
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

```
72 of 137 RxLRs were in the predicted 797 effectorP genes.
```bash
  cat analysis/RxLR_effectors/combined_evidence/P.idaei/SCRP370/SCRP370_total_RxLR_headers.txt $OutFileHeaders | sort | cut -f1 -d '.' | uniq -d | wc -l
```

### D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:


```bash
  HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
  LFLAK_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm)
  DWL_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm)
  for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta | grep -w -e '10300'); do
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
Initial search space (Z):              24161  [actual number of targets]
Domain search space  (domZ):             117  [number of targets reported over threshold]
Initial search space (Z):              24161  [actual number of targets]
Domain search space  (domZ):              94  [number of targets reported over threshold]
74
```
<!--
23 of 92 P.idaei CRNs were in the effectorP dataset.
```bash
cat analysis/CRN_effectors/hmmer_CRN/P.idaei/SCRP370/SCRP370_pub_CRN_LFLAK_DWL.txt $OutFileHeaders | sort | cut -f1 -d '.' | uniq -d | wc -l
``` -->

Extract gff annotations for Crinklers:

```bash
  for CRNlist in $(ls analysis/CRN_effectors/hmmer_CRN/*/*/*_pub_CRN_LFLAK_DWL.txt | grep -e '10300'); do
    Strain=$(echo $CRNlist | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $CRNlist | rev | cut -f3 -d '/' | rev)
    OutName=$(echo $CRNlist | sed 's/.txt/.gff/g')
    echo "$Organism - $Strain"
    Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended.gff3)
    cat $CRNlist | sed -r 's/\.t.$//g' > tmp.txt
    cat $Gff | grep -w -f tmp.txt > $OutName
    rm tmp.txt
  done
```



### E) From ORF gene models - Signal peptide & RxLR motif

Required programs:
 * SigP
 * Phobius
 * biopython


#### E.1) Prediction using SignalP
Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep '10300'); do
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
while [ $Jobs -gt 15 ]; do
sleep 1
printf "."
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
done		
printf "\n"
echo $File
qsub $ProgDir/pred_sigP.sh $File
qsub $ProgDir/pred_sigP.sh $File signalp-3.0
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```


 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:

```bash
for SplitDir in $(ls -d gene_pred/ORF_split/*/* | grep -w -e 'P.cactorum' -e 'P.idaei' | grep -v -e '10300' -e '414_v2'); do
Strain=$(echo $SplitDir | cut -d '/' -f4)
Organism=$(echo $SplitDir | cut -d '/' -f3)
echo "$Organism - $Strain"
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
for SigpDir in $(ls -d gene_pred/ORF_sig* | cut -f2 -d'/'); do
for GRP in $(ls -l $SplitDir/*_ORF_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.aa";
InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp_neg.aa";
InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.tab";
InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.txt";
done
cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
done
done
```


#### E.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
	for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep '10300'); do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/phobius/$Organism/$Strain
    mkdir -p $OutDir
  	phobius.pl $Proteome > $OutDir/"$Strain"_phobius_ORF.txt
  	cat $OutDir/"$Strain"_phobius_ORF.txt | grep -B1 'SIGNAL' | grep 'ID' | sed s'/ID.*g/g/g' > $OutDir/"$Strain"_phobius_headers_ORF.txt
  done
```

Secreted proteins from different sources were combined into a single file:

```bash
  for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '10300' -e '414_v2' | grep -w -e '404' -e '414' -e '415' -e '416'); do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/combined_sigP_ORF/$Organism/$Strain
    mkdir -p $OutDir
    echo "The following number of sequences were predicted as secreted:"
    # cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
    cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa > $OutDir/"$Strain"_all_secreted.fa
    cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | tr -d '>' | tr -d ' ' | sed "s/HMM_score\t/HMM_score=\t/g" > $OutDir/"$Strain"_all_secreted_headers.txt
    cat $OutDir/"$Strain"_all_secreted_headers.txt | wc -l
    echo "This represented the following number of unique genes:"
    # cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    # $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
    cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
  done
```

```

```

#### E.3) Prediction of RxLRs


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
overlaps and identify the ORF with the best signalP score.

```bash
  for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff3 | grep -v -e 'atg' | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '10300' -e '414_v2' | grep -v -w -e '404' -e '414' -e '415' -e '416'); do
    Organism=$(echo $ORF_Gff | rev |  cut -d '/' -f3 | rev) ;
    Strain=$(echo $ORF_Gff | rev | cut -d '/' -f2 | rev);
    OutDir=$(ls -d gene_pred/combined_sigP_ORF/$Organism/$Strain)
    echo "$Organism - $Strain"
    # SigP_fasta=$(ls $OutDir/"$Strain"_all_secreted.fa)
    SigP_headers=$(ls $OutDir/"$Strain"_all_secreted_headers.txt)
    ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)

    SigP_Gff=$OutDir/"$Strain"_all_secreted_unmerged.gff
    SigP_Merged_Gff=$OutDir/"$Strain"_all_secreted_merged.gff
    SigP_Merged_txt=$OutDir/"$Strain"_all_secreted_merged.txt
    SigP_Merged_AA=$OutDir/"$Strain"_all_secreted_merged.aa

    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $SigP_headers $ORF_Gff SigP Name > $SigP_Gff
    ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $SigP_Gff --db sigP_ORF.db
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF.db --id sigP_ORF --out sigP_ORF_merged.db --gff > $SigP_Merged_Gff
    cat $SigP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d'=' | rev > $SigP_Merged_txt
    # $ProgDir/extract_from_fasta.py --fasta $SigP_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
    $ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
  done
```

#<--- Progress up to here
The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '10300' -e '414_v2' | grep -v -w -e '404' -e '414' -e '415' -e '416'); do
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
mkdir -p $OutDir;
printf "\nstrain: $Strain\tspecies: $Organism\n";
printf "the number of SigP gene is:\t";
cat $Secretome | grep '>' | wc -l;
printf "the number of SigP-RxLR genes are:\t";
$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa;
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt
cat $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l
printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' '> $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt  $SigP_Merged_Gff 	RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
SigP_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt  $SigP_Gff	RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt $SigP_Merged_Gff RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex.gff
RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.gff
RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.txt
RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff --db sigP_ORF_RxLR.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR.db --id sigP_ORF_RxLR --out sigP_ORF_RxLR_merged.db --gff > $RxLR_Merged_Gff
cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
printf "Merged RxLR-EER regex proteins:\t"
cat $RxLR_Merged_AA | grep '>' | wc -l
printf "\n"
done
```

```

```


### E5) From ORF gene models - Hmm evidence of WY domains
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

<!-- P.cactorum 10300
Initial search space (Z):              15271  [actual number of targets]
Domain search space  (domZ):             113  [number of targets reported over threshold] -->

* P.cactorum 10300
* Initial search space (Z):              14767  [actual number of targets]
* Domain search space  (domZ):             113  [number of targets reported over threshold]


### E6) From ORF gene models - Hmm evidence of RxLR effectors

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep -w -e '414' | grep -v -w -e '404' -e '415' -e '416'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_ORF_RxLR_hmmer_unmerged.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_ORF_RxLR_hmmer_headers_unmerged.txt
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
SigP_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_RxLR_hmmer_unmerged.gff3
RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.gff
RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.txt
RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_hmmer_unmerged.gff3 --db sigP_ORF_RxLR_hmm.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR_hmm.db --id sigP_ORF_RxLR_hmm --out sigP_ORF_RxLR_hmm_merged.db --gff > $RxLR_Merged_Gff
cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
printf "Merged RxLR-EER Hmm proteins:\t"
cat $RxLR_Merged_AA | grep '>' | wc -l
done
```
```

```

* P.cactorum 10300
* Initial search space (Z):              14767  [actual number of targets]
* Domain search space  (domZ):             144  [number of targets reported over threshold]

### E7) Combining RxLRs from Regex and hmm searches


The total RxLRs are

```bash
  for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_ORF_RxLR_EER_regex_merged.txt | grep -v -e 'Aug' -e '10300' | grep -e 'P.idaei' -e 'P.cactorum' | grep -e '414' -e '404' -e '415' -e '416'); do
    Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
    Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF.gff3)
    Proteome=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
    HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_ORF_RxLR_hmm_merged.txt)
    echo "$Organism - $Strain"
    echo "Number of RxLRs identified by Regex:"
    cat $RegexRxLR | sort | uniq | wc -l
    echo "Number of RxLRs identified by Hmm:"
    cat $HmmRxLR | sort | uniq | wc -l
    echo "Number of RxLRs in combined dataset:"
    cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
    # echo "Number of RxLRs in both datasets:"
    # cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
    echo ""
    # echo "Extracting RxLRs from datasets"
    OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
    mkdir -p $OutDir
    cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
    # cat $Gff | grep -w -f $OutDir/"$Strain"_total_ORF_RxLR_headers.txt > $OutDir/"$Strain"_total_ORF_RxLR.gff
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_total_ORF_RxLR_headers.txt $Gff ORF_RxLR Name Augustus > $OutDir/"$Strain"_total_ORF_RxLR.gff
    echo "Number of genes in the extracted gff file:"
    cat $OutDir/"$Strain"_total_ORF_RxLR.gff | grep -w 'gene' | wc -l
  done
```

## 4.2.c Analysis of RxLR effectors - merger of Augustus / published genes with ORFs

Intersection between the coodinates of putative RxLRs from gene models and ORFs
were identified to determine the total number of RxLRs predicted in these
genomes.

The RxLR effectors from both Gene models and ORF finding approaches were
combined into a single file.

This step was complicated by the inconsistency in downloaded gff files for gene
models.


```bash
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/*/* | grep -e 'P.idaei' -e 'P.cactorum' | grep -w -e '414' -e '404' -e '415' -e '416'); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$MergeDir/"$Strain"_total_RxLR.gff
AugTxt=$MergeDir/"$Strain"_total_RxLR_headers.txt
AugFa=$(ls gene_pred/codingquary/"$Species"/"$Strain"/final/final_genes_combined.pep.fasta)

ORFGff=$(ls $MergeDir/"$Strain"_total_ORF_RxLR.gff)
ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
ORFsTxt=$(ls $MergeDir/"$Strain"_total_ORF_RxLR_headers.txt)

ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_RxLR_EER_motif_hmm.gff
AugInORFs=$MergeDir/"$Strain"_AugInORFs_RxLR_EER_motif_hmm.gff
ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_RxLR_EER_motif_hmm.gff
AugUniq=$MergeDir/"$Strain"_Aug_Uniq_RxLR_EER_motif_hmm.gff
TotalRxLRsTxt=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.txt
TotalRxLRsGff=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.gff

bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq

echo "$Species - $Strain"
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l
# cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
echo "The total number of putative RxLRs are:"
cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalRxLRsTxt
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
cat $ORFsUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f3 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
cat $TotalRxLRsTxt | wc -l
cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalRxLRsTxt > $TotalRxLRsGff

RxLRsFa=$MergeDir/"$Strain"_final_RxLR_EER.fa
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/unwrap_fasta.py --inp_fasta $AugFa | grep -A1 -w -f $AugTxt | grep -v -E '^--$' > $RxLRsFa
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalRxLRsTxt > $RxLRsFa
# echo "$Strain"
$ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalRxLRsTxt >> $RxLRsFa
# echo "$Strain"
echo "The number of sequences extracted is"
cat $RxLRsFa | grep '>' | wc -l
done
```

```

```



### H) From ORF gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in ORF gene models. This was done with the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/*/*/*.aa_cat.fa | grep -w -e 'P.cactorum' -e 'P.idaei' | grep -v -e 'atg' -e '10300' -e '414_v2' | grep -v -w -e '404' -e '414' -e '415' -e '416'); do
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
Headers=$OutDir/"$Strain"_CRN_hmmer_unmerged_headers.txt
cat $OutDir/$HmmFastaDWL | grep '>' | grep -w -f $CommonHeaders | tr -d '>' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | tr -d '-' | sed 's/hmm_score/HMM_score/g' > $Headers
# As we are dealing with JGI and Broad sequences, some features need formatting:
ORF_Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/*_ORF.gff3)
# Gff features were extracted for each header
CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $Headers $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
# Gff features were merged based upon the DWL hmm score
DbDir=analysis/databases/$Organism/$Strain
mkdir -p $DbDir
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $CRN_unmerged_Gff --db $DbDir/CRN_ORF.db
CRN_Merged_Gff=$OutDir/"$Strain"_CRN_merged_hmmer.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp $DbDir/CRN_ORF.db --id LFLAK_DWL_CRN --out $DbDir/CRN_ORF_merged.db --gff > $CRN_Merged_Gff
# Final results are reported:
echo "Number of CRN ORFs after merging:"
cat $CRN_Merged_Gff | grep 'gene' | wc -l
done
```

```

```



Extract crinklers from published gene models


```bash
  for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/*/* | grep -w -e 'P.cactorum' -e 'P.idaei' | grep -v -e 'atg' -e '10300' -e '414_v2' | grep -w -e '404' -e '414' -e '415' -e '416'); do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    AugGff=$(ls $MergeDir/"$Strain"_pub_CRN_LFLAK_DWL.gff)
    AugFa=$(ls gene_pred/codingquary/"$Species"/"$Strain"/final/final_genes_combined.pep.fasta)
    ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
    ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer.gff3
    ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_CRN_hmmer.bed
    AugInORFs=$MergeDir/"$Strain"_AugInORFs_CRN_hmmer.bed
    ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_CRN_hmmer.bed
    AugUniq=$MergeDir/"$Strain"_Aug_Uniq_CRN_hmmer.bed
    TotalCRNsTxt=$MergeDir/"$Strain"_final_CRN.txt
    TotalCRNsGff=$MergeDir/"$Strain"_final_CRN.gff
    TotalCRNsHeaders=$MergeDir/"$Strain"_Total_CRN_headers.txt
    bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
    bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
    bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
    bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq
    echo "$Species - $Strain"

    echo "The number of ORF CRNs overlapping Augustus CRNs:"
    cat $ORFsInAug | grep -w -e 'transcript' -e 'mRNA' | wc -l
    echo "The number of Augustus CRNs overlapping ORF CRNs:"
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA' | wc -l
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalCRNsTxt
    echo "The number of CRNs unique to ORF models:"
    cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' | wc -l
    cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt
    echo "The number of CRNs unique to Augustus models:"
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt

    cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalCRNsTxt > $TotalCRNsGff

    CRNsFa=$MergeDir/"$Strain"_final_CRN.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalCRNsTxt > $CRNsFa
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalCRNsTxt >> $CRNsFa
    echo "The number of sequences extracted is"
    cat $CRNsFa | grep '>' | wc -l

  done
```

```
```