

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

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Single=$(cat $File | grep "(S)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
  done
```




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
<!--
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
``` -->


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
    cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_appended_renamed.pep.fasta
    cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
    cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
    cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta


    GffBraker=$FinalDir/final_genes_Braker.gff3
    GffQuary=$FinalDir/final_genes_CodingQuary.gff3
    GffAppended=$FinalDir/final_genes_appended_renamed.gff3
    cat $GffBraker $GffQuary > $GffAppended
  done
```

In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified


```bash
  for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3 | grep '10300'); do
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
  for OrfGff in $(ls gene_pred/ORF_finder/P.*/*/*_ORF.gff | grep -v 'atg' | grep -e 'P.cactorum' -e 'P.idaei' | grep -e '414' -e '404'); do
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
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '10300'); do
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
    for SigpDir in $(ls -d gene_pred/final_sig* | cut -f2 -d'/'); do
      InStringAA=''
      InStringNeg=''
      InStringTab=''
      InStringTxt=''
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
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '10300'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/phobius/$Organism/$Strain
    mkdir -p $OutDir
    phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    $ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
    cat $OutDir/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d '>' > $OutDir/"$Strain"_phobius_headers.txt
  done
 ```


Secreted proteins from different sources were combined into a single file:

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '10300'); do
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
P.cactorum - 10300
The following number of sequences were predicted as secreted:
17606
This represented the following number of unique genes:
3292
```

Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '10300'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
qsub $ProgDir/submit_TMHMM.sh $Proteome
done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt | grep '10300'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
NonTmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $NonTmHeaders
SigP=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_secreted.fa)
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $NonTmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
echo "Number of SigP proteins:"
cat $SigP | grep '>' | wc -l
echo "Number without transmembrane domains:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
echo "Number of gene models:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l

# A text file was also made containing headers of proteins testing +ve
PosFile=$(ls gene_pred/trans_mem/$Organism/$Strain/"$Strain"_TM_genes_pos.txt)
TmHeaders=$(echo $PosFile | sed 's/.txt/_headers.txt/g')
cat $PosFile | cut -f1 > $TmHeaders

done
```

```
P.cactorum - 10300
Number of SigP proteins:
3292
Number without transmembrane domains:
2244
Number of gene models:
2233
```

Proteins containing GPI anchors were also removed using GPIsom


These proteins were identified through submitting the combined protein file to
the webserver at: http://gpi.unibe.ch


An output directory was made to download the file to:
"GPI anchored (C&N-term signal) (SignalP):"

```bash
  for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '10300'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/trans_mem/$Organism/$Strain/GPIsom
    mkdir -p $OutDir
  done
```

The link to the proteins results was followed embedded in the text of
"GPI anchored (C&N-term signal) (SignalP):" thes +GPI anchor proteins
were copied and pasted into:

<!-- The link to the results was followed embedded in the text of
"Seqs with C-terminal signal (GPI-SOM):" thes +GPI anchor proteins
were copied and pasted into: -->

```bash
  nano gene_pred/trans_mem/P.cactorum/10300/GPIsom/GPI_pos.fa
```

Those proteins with GPI anchors were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/GPIsom/GPI_pos.fa | grep '10300'); do
Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
TmHeaders=$(echo "$File" | sed 's/.fa/.txt/g')
cat $File | grep '>' | cut -f1 -d ' ' | sed 's/>//g' > $TmHeaders
SigP=$(ls gene_pred/combined_sigP/$Organism/$Strain/*_sp_no_trans_mem.aa)
SigPHeaders=gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_sp_no_trans_mem_headers.txt
cat $SigP | grep '>' | cut -f1 | sed 's/>//g'> $SigPHeaders
GoodHeaders=$(echo "$File" | sed 's/_pos.fa/_neg.txt/g')
cat $SigPHeaders | grep -v -f $TmHeaders > $GoodHeaders
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# cat $SigP | grep -v -A1 -f $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $GoodHeaders  > $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa
echo "Number of SigP proteins:"
cat $SigP | grep '>' |wc -l
echo "Number with GPI anchors in entire proteome:"
cat $TmHeaders | wc -l
echo "Number without GPI anchors:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa | grep '>' | wc -l
echo "Number of gene models:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l
done
```

```
  P.cactorum - 10300
  Number of SigP proteins:
  2244
  Number with GPI anchors in entire proteome:
  483
  Number without GPI anchors:
  1984
  Number of gene models:
  1973
```



### C) From Augustus gene models - Effector-like structure identification using EffectorP

Required programs:
 * EffectorP.py

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '10300'); do
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

Note - this doesnt exclude proteins with TM domains or GPI anchors

```bash
  for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt | grep '10300'); do
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
    Gff=$(ls gene_pred/final/$Organism/$Strain/*/final_genes_appended_renamed.gff3)
    EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
  done
```
```
P.cactorum - 10300
EffectorP headers:	13994
Secreted effectorP headers:	845
```

## D) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '10300'); do
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
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep '10300'); do
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
Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
CazyGff=$OutDir/"$Strain"_CAZY.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

# SecretedProts=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_secreted.fa)
SecretedProts=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem_no_GPI.aa)
SecretedHeaders=$(echo $SecretedProts | sed 's/.fa/_headers.txt/g' | sed 's/.aa/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
printf "number of Secreted CAZY genes identified:\t"
cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
done
```

```
P.cactorum - 10300
number of CAZY genes identified:	697
number of Secreted CAZY genes identified:	372
```

or, if those secreted proteins with TM domains or GPI anchors are removed:

```
P.cactorum - 10300
number of CAZY genes identified:	697
number of Secreted CAZY genes identified:	236
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

### Summary of CAZY families by organism


```bash
for CAZY in $(ls gene_pred/CAZY/*/*/*_CAZY.out.dm.ps | grep '10300'); do
  Strain=$(echo $CAZY | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $CAZY | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $CAZY)
  echo "$Organism - $Strain"
  Secreted=$(ls gene_pred/combined_sigP/$Organism/$Strain/*_sp_no_trans_mem_no_GPI_headers.txt)
  Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/CAZY
  $ProgDir/summarise_CAZY.py --cazy $CAZY --inp_secreted $Secreted --inp_gff $Gff --summarise_family --trim_gene_id 2 --kubicek_2014
done
```

```
  P.cactorum - 10300
  Polygalacturonase - 17
  A-Arabinosidases - 5
  Xylanases - 12
  Polygalacturonate lyases - 34
  B-Glycosidases - 16
  Cellulases - 26
  other - 126
```

The following commands were used to identify b-glucan binding proteins as described in:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3627985/

```bash
  for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep '10300'); do
    cat $File | grep -w -e 'GH5' -e 'GH55' -e 'GH81' -e 'GH16' -e 'GH3' -e 'GH17' -e 'GH72' | sed -r "s/\s+/\t/g" | cut -f1,4 |  sort | uniq > tmp1.txt
    cat tmp1.txt | cut -f2 > tmp2.txt
    cat tmp1.txt | grep -f tmp3.txt | cut -f1 | sort | uniq -c
  done
```

```
  2 GH16.hmm
  18 GH17.hmm
  14 GH3.hmm
  14 GH5.hmm
  11 GH81.hmm
```

# D) Prediction of RxLRs from - Augustus/CodingQuary gene models

 The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

 The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep '10300'); do
Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
Proteome=$(ls gene_pred/final/$Organism/$Strain/*/final_genes_appended_renamed.pep.fasta)
Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain"
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
  strain: 10300	species: P.cactorum
  the total number of SigP gene is:	10309
  the number of unique SigP gene is:	3292
  the number of SigP-RxLR genes are:	294
  the number of SigP-RxLR-EER genes are:	130
```


### G) From Predicted gene models - Hmm evidence of RxLR effectors

```bash
  for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '10300'); do
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
    Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
    cat $Gff | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
  done
```

```
P.cactorum 10300
Initial search space (Z):              24161  [actual number of targets]
Domain search space  (domZ):             148  [number of targets reported over threshold]
```


### F) Combining RxLRs from Regex and hmm searches


The total RxLRs are

```bash
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_RxLR_EER_regex.txt | grep -v -e 'Aug' -e 'ORF' | grep '10300'); do
Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
Proteome=$(ls gene_pred/final/$Organism/$Strain/*/final_genes_appended_renamed.pep.fasta)
HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer_headers.txt | grep -v -e 'Aug' -e 'ORF')
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
Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
cat $Gff | grep -w -f $OutDir/"$Strain"_total_RxLR_headers.txt > $OutDir/"$Strain"_total_RxLR.gff
echo "Number of genes in the extracted gff file:"
cat $OutDir/"$Strain"_total_RxLR.gff | grep -w 'gene' | wc -l
done
```

```
P.cactorum - 10300
Number of RxLRs identified by Regex:
130
Number of RxLRs identified by Hmm:
148
Number of RxLRs in combined dataset:
176
Number of genes in the extracted gff file:
176
```



```bash
for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep '10300'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_Aug_WY_hmmer.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_Aug_WY_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_Aug_WY_hmmer_headers.txt
# cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' | sort | uniq > $OutDir/$Headers
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | cut -f1 | sort | uniq > $OutDir/$Headers
Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
cat $Gff | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_WY_hmmer.gff
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/$Headers $Gff $HmmModel Name Augustus > $OutDir/"$Strain"_Aug_WY_hmmer.gff
done
```


<!-- 72 of 137 RxLRs were in the predicted 797 effectorP genes.
```bash
  cat analysis/RxLR_effectors/combined_evidence/P.idaei/SCRP370/SCRP370_total_RxLR_headers.txt $OutFileHeaders | sort | cut -f1 -d '.' | uniq -d | wc -l
``` -->

### D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:


```bash
  HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
  LFLAK_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm)
  DWL_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm)
  for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep -w -e '10300'); do
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
    Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
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
for SplitDir in $(ls -d gene_pred/ORF_split/*/* | grep '10300'); do
Strain=$(echo $SplitDir | cut -d '/' -f4)
Organism=$(echo $SplitDir | cut -d '/' -f3)
echo "$Organism - $Strain"
for SigpDir in $(ls -d gene_pred/ORF_sig* | cut -f2 -d'/'); do
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
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
  for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep '10300'); do
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
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
    cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
  done
```

```
P.cactorum - 10300
The following number of sequences were predicted as secreted:
55959
This represented the following number of unique genes:
27593
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

Due to the nature of predicting ORFs, some features overlapped with one another.
A single ORF was selected from each set of overlapped ORFs. This was was
selected on the basis of its SignalP Hmm score. Biopython was used to identify
overlaps and identify the ORF with the best signalP score.

```bash
  for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff3 | grep -v -e 'atg' | grep '10300'); do
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


```
  8494
```

The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep '10300'); do
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
mkdir -p $OutDir;
printf "\nstrain: $Strain\tspecies: $Organism\n";
printf "the number of SigP gene is:\t";
cat $Secretome | grep '>' | cut -f1 | tr -d ' ' | sort | uniq | wc -l
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
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt $SigP_Gff	RxLR_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt  $SigP_Gff	RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.gff
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt $SigP_Merged_Gff RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex.gff
RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_regex_merged.gff
RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_regex_merged.txt
RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_regex_merged.aa
RxLR_EER_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.gff
RxLR_EER_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.txt
RxLR_EER_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff --db sigP_ORF_RxLR.db
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.gff --db sigP_ORF_RxLR_EER.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR.db --id sigP_ORF_RxLR --out sigP_ORF_RxLR_merged.db --gff > $RxLR_Merged_Gff
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR_EER.db --id sigP_ORF_RxLR_EER --out sigP_ORF_RxLR_EER_merged.db --gff > $RxLR_EER_Merged_Gff
cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
cat $RxLR_EER_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_EER_Merged_txt
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_EER_Merged_txt > $RxLR_EER_Merged_AA
printf "Merged RxLR regex proteins:\t"
cat $RxLR_Merged_AA | grep '>' | wc -l
printf "Merged RxLR-EER regex proteins:\t"
cat $RxLR_EER_Merged_AA | grep '>' | wc -l
printf "\n"
done
```

```
strain: 10300	species: P.cactorum
the number of SigP gene is:	27593
the number of SigP-RxLR genes are:	1589
the number of SigP-RxLR-EER genes are:	211
Merged RxLR regex proteins:	1314
Merged RxLR-EER regex proteins:	182
```

Quantification of ORF RxLRs was also performed.


#### Aligning


```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FileF=$(ls qc_rna/raw_rna/genbank/P.cactorum/F/SRR1206032_trim.fq.gz)
FileR=$(ls qc_rna/raw_rna/genbank/P.cactorum/R/SRR1206033_trim.fq.gz)
echo $FileF
echo $FileR
Prefix="genbank"
Timepoint="treatment"
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
```

<!--


```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
InGff=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain*/*_ORF_RxLR_regex_merged.gff)
InGff_features=$(echo $InGff | sed 's/.gff/_features.gff/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/add_ORF_features.pl $InGff $Assembly >> $InGff_features
for BamFile in $(ls alignment/star/P.cactorum/10300/treatment/genbank/star_aligmentAligned.sortedByCoord.out.bam); do
# OutDir=alignment/star/$Organism/$Strain/featureCounts_genbank
OutDir=$(dirname $BamFile)
mkdir -p $OutDir
# Prefix=$(echo $BamFile | rev | cut -f3 -d '/' | rev | sed 's/vesca_//g')
Prefix=$Strain
echo $Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
# echo "$BamFile $Gff $OutDir $Prefix"
qsub $ProgDir/sub_featureCounts.sh $BamFile $InGff_features $OutDir $Prefix
done
done
```

Those RxLRs with evidence of expression in planta (fpkm >5) were included in
ORF RxLR gene models:

```bash
  Organism="P.cactorum"
  Strain="10300"
  OutDir=$(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain)
  for FeatureCounts in $(ls alignment/star/$Organism/$Strain/treatment/genbank/*_featurecounts.txt); do
    OutFile=$(echo $FeatureCounts | sed 's/.txt/_fpkm_>5.txt/g')
  # cat $FeatureCounts | tail -n+3 | cut -f1,7 | grep -v -e '\s0$' -e '\s1$' -e '\s2$' -e '\s3$' -e '\s4$' > $OutFile
    cat $FeatureCounts | tail -n+3 | awk '$7>5 {print $0}' > $OutFile
  cat $OutFile
  done | cut -f1 | sort | uniq > $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers.txt
  cat $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers.txt | wc -l
  SigP_Gff=$(ls gene_pred/combined_sigP_ORF/P.cactorum/10300/10300_all_secreted_unmerged.gff)
  cat $SigP_Gff | grep -f $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers.txt | cut -f4 -d '=' | cut -f1 -d ';' > $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers_renamed.txt
  cat $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers_renamed.txt | wc -l
```

The number of Aditional RxLRs from those  secreted ORF RxLRs with a fpkm >5 was:

```
  904
``` -->



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
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep '10300'); do
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
P.cactorum 10300
Initial search space (Z):              55959  [actual number of targets]
Domain search space  (domZ):             471  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	145
```

<!--
* P.cactorum 10300
* Initial search space (Z):              14767  [actual number of targets]
* Domain search space  (domZ):             144  [number of targets reported over threshold]
 -->

### E7) Combining RxLRs from Regex and hmm searches


 The total ORF RxLRs are

 ```bash
for RegexRxLREER in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_ORF_RxLR_EER_regex_merged.txt | grep '10300'); do
 Organism=$(echo $RegexRxLREER | rev |  cut -d '/' -f3 | rev)
 Strain=$(echo $RegexRxLREER | rev | cut -d '/' -f2 | rev)
 Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF.gff3)
 Proteome=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
 # RegexRxLRfpkm=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers_renamed.txt)
 HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_ORF_RxLR_hmm_merged.txt)
 echo "$Organism - $Strain"
 echo "Number of RxLR EERs identified by Regex:"
 cat $RegexRxLREER | sort | uniq | wc -l
 # echo "Number of RxLRs identified by Regex with fpkm > 5"
 # cat $RegexRxLRfpkm | sort | uniq | wc -l
 echo "Number of RxLRs identified by Hmm:"
 cat $HmmRxLR | sort | uniq | wc -l
 echo "Number of RxLRs in combined dataset:"
 # cat $RegexRxLREER $HmmRxLR $RegexRxLRfpkm | sort | uniq | wc -l
 cat $RegexRxLREER $HmmRxLR | sort | uniq | wc -l
 echo ""
 # echo "Extracting RxLRs from datasets"
 OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
 mkdir -p $OutDir
 # cat $RegexRxLREER $RegexRxLRfpkm $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
 cat $RegexRxLREER $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
 ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
 $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_total_ORF_RxLR_headers.txt $Gff ORF_RxLR Name Augustus > $OutDir/"$Strain"_total_ORF_RxLR.gff
 echo "Number of genes in the extracted gff file:"
 cat $OutDir/"$Strain"_total_ORF_RxLR.gff | grep -w 'gene' | wc -l
done
```
<!--
When including SigP RxLRs with an fpkm >5

```
P.cactorum - 10300
The number of ORF RxLRs overlapping Augustus RxLRs:
149
The number of Augustus RxLRs overlapping ORF RxLRs:
149
The number of RxLRs unique to ORF models:
864
The number of RxLRs unique to Augustus models:
27
The total number of putative RxLRs are:
1040
The number of sequences extracted is
1040
```
-->

```
P.cactorum - 10300
Number of RxLR EERs identified by Regex:
182
Number of RxLRs identified by Hmm:
145
Number of RxLRs in combined dataset:
205

Number of genes in the extracted gff file:
205
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
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/*/* | grep '10300'); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$MergeDir/"$Strain"_total_RxLR.gff
AugTxt=$MergeDir/"$Strain"_total_RxLR_headers.txt
AugFa=$(ls gene_pred/final/"$Species"/"$Strain"/final/final_genes_appended_renamed.pep.fasta)

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
  P.cactorum - 10300
  The number of ORF RxLRs overlapping Augustus RxLRs:
  149
  The number of Augustus RxLRs overlapping ORF RxLRs:
  149
  The number of RxLRs unique to ORF models:
  56
  The number of RxLRs unique to Augustus models:
  27
  The total number of putative RxLRs are:
  232
  The number of sequences extracted is
  232
```



### H) From ORF gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in ORF gene models. This was done with the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/*/*/*.aa_cat.fa | grep '10300'); do
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
ORF_Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/*_ORF.gff3 | grep -v '_atg_')
# Gff features were extracted for each header
CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
# cat $Headers | cut -f1 > tmp.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $Headers $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
# $ProgDir/extract_gff_for_sigP_hits.pl tmp.txt $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
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
```

Extract crinklers from ORFs and Braker/codingquary gene models


```bash
  for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/*/* | grep '10300'); do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    AugGff=$(ls $MergeDir/"$Strain"_pub_CRN_LFLAK_DWL.gff)
    AugFa=$(ls gene_pred/final/"$Species"/"$Strain"/final/final_genes_appended_renamed.pep.fasta)
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
  P.cactorum - 10300
  The number of ORF CRNs overlapping Augustus CRNs:
  55
  The number of Augustus CRNs overlapping ORF CRNs:
  56
  The number of CRNs unique to ORF models:
  3
  The number of CRNs unique to Augustus models:
  18
  The number of sequences extracted is
  77
```



## C.ii) SSCP

Small secreted cysteine rich proteins were identified within secretomes. These
proteins may be identified by EffectorP, but this approach allows direct control
over what constitutes a SSCP.

```bash

for Secretome in $(ls gene_pred/combined_sigP/*/*/*_secreted.fa | grep -v 'all_secreted' |  grep '10300'); do
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/sscp/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_all_results.fa
cat $OutDir/"$Strain"_sscp_all_results.fa | grep 'Yes' > $OutDir/"$Strain"_sscp.fa
printf "number of SSC-rich genes:\t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' | cut -f1 -d '.' | sort | uniq | wc -l
printf "Number of effectors predicted by EffectorP:\t"
EffectorP=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP_secreted_headers.txt)
cat $EffectorP | wc -l
printf "Number of SSCPs predicted by both effectorP and this approach: \t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' > $OutDir/"$Strain"_sscp_headers.txt
cat $OutDir/"$Strain"_sscp_headers.txt $EffectorP | cut -f1 | sort | uniq -d | wc -l
echo ""
done
```

```
P.cactorum - 10300
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	388
number of SSC-rich genes:	386
Number of effectors predicted by EffectorP:	845
Number of SSCPs predicted by both effectorP and this approach: 	223
```
<!--
```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep '10300'); do
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"

OutDir=analysis/sscp_ORF/$Organism/$Strain
mkdir -p $OutDir
printf "the number of SigP gene is:\t";
cat $Secretome | grep '>' | cut -f1 | sort | uniq | wc -l
printf "the number of SSCP genes are:\t";
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_ORF.fa

cat $OutDir/"$Strain"_sscp_ORF.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/"$Strain"_sscp_ORF_headers_unmerged.txt
cat $OutDir/"$Strain"_sscp_ORF_headers_unmerged.txt | tr -d ' ' | sort | uniq | wc -l

SigP_Gff=$(ls gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_sscp_ORF_headers_unmerged.txt  $SigP_Gff	sscp_filter.py Name Augustus > $OutDir/"$Strain"_sscp_ORF_unmerged.gff

SSCP_Merged_Gff=$OutDir/"$Strain"_ORF_sscp_merged.gff
SSCP_Merged_txt=$OutDir/"$Strain"_ORF_sscp_merged.txt
SSCP_Merged_AA=$OutDir/"$Strain"_ORF_sscp_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_sscp_ORF_unmerged.gff --db sigP_ORF_sscp.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_sscp.db --id sigP_ORF_sscp --out sigP_ORF_sscp_merged.db --gff > $SSCP_Merged_Gff
cat $SSCP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $SSCP_Merged_txt

ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $SSCP_Merged_txt > $SSCP_Merged_AA
printf "Merged SSCP proteins:\t"
cat $SSCP_Merged_AA | grep '>' | wc -l
printf "\n"
done
```
-->



# Making a combined file of Braker, Coding quary genes as well as additional ORF effector candidates


A gff file containing the combined Braker and CodingQuary genes as well as the
additional CRN and RxLR genes predicted by ORF analysis was made.

```bash
GeneGff=$(ls gene_pred/final/P.cactorum/10300/final/final_genes_appended_renamed.gff3)
GffOrfRxLR=$(ls analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_ORFsUniq_RxLR_EER_motif_hmm.gff)
GffOrfCRN=$(ls analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORFsUniq_CRN_hmmer.bed)
Assembly=$(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_softmasked_repeatmasker_TPSI_appended.fa)
OutDir=gene_pred/final_incl_ORF/P.cactorum/10300
mkdir -p $OutDir
cat $GeneGff > $OutDir/10300_genes_incl_ORFeffectors.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/add_ORF_features.pl $GffOrfRxLR $Assembly >> $OutDir/10300_genes_incl_ORFeffectors.gff3
$ProgDir/add_ORF_features.pl $GffOrfCRN $Assembly >> $OutDir/10300_genes_incl_ORFeffectors.gff3
# Make gene models from gff files.
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
Assembly=$(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_softmasked_repeatmasker_TPSI_appended.fa)
$ProgDir/gff2fasta.pl $Assembly $OutDir/10300_genes_incl_ORFeffectors.gff3 $OutDir/10300_genes_incl_ORFeffectors
```

High confidence genes overlapped by ORFs to be inserted were manually inspected
to see whether the original gene model or the ORF should be kept. Preference was
given to the ORF effector candidate. BEfore this could be done bedtools was
used to identify which genes were intersected:

```bash
# GeneGff=$(ls gene_pred/final/P.cactorum/10300/final/final_genes_appended.gff3)
GeneGff=$(ls gene_pred/final/P.cactorum/10300/final/final_genes_appended_renamed.gff3)
GffOrfRxLR=$(ls analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_ORFsUniq_RxLR_EER_motif_hmm.gff)
GffOrfCRN=$(ls analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORFsUniq_CRN_hmmer.bed)
Assembly=$(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_softmasked_repeatmasker_TPSI_appended.fa)
bedtools intersect -wo -a $GeneGff -b $GffOrfRxLR | grep -e "AUGUSTUS.gene" | grep "ORF_RxLR.gene"| cut -f1,9,18
bedtools intersect -wo -a $GeneGff -b $GffOrfCRN | grep -e "AUGUSTUS.gene" | grep "CRN_HMM.gene"| cut -f1,9,18

```

```
contig_6	ID=g553;	ID=ORF_R_6618
contig_12	ID=g1094;	ID=ORF_F_12367
contig_12	ID=g1095;	ID=ORF_F_12367
contig_27	ID=g2191;	ID=ORF_F_23817
contig_28	ID=g2212;	ID=ORF_R_23034
contig_32	ID=g2443;	ID=ORF_R_25613
contig_37	ID=g2771;	ID=ORF_F_29498
contig_39	ID=g2890;	ID=ORF_F_30553
contig_50	ID=g3520;	ID=ORF_R_35309
contig_70	ID=g4524;	ID=ORF_R_46076
contig_105	ID=g6026;	ID=ORF_F_61898
contig_140	ID=g7303;	ID=ORF_F_74889
contig_167	ID=g8212;	ID=ORF_R_84315
contig_177	ID=g8502;	ID=ORF_R_87287
contig_184	ID=g8735;	ID=ORF_F_89337
contig_295	ID=g11604;	ID=ORF_R_116876
contig_338	ID=g12540;	ID=ORF_R_126392
contig_376	ID=g13291;	ID=ORF_R_133618
contig_446	ID=g14477;	ID=ORF_F_144470
contig_481	ID=g15058;	ID=ORF_F_149698
contig_505	ID=g15388;	ID=ORF_F_152796
contig_510	ID=g15458;	ID=ORF_F_153549
contig_511	ID=g15477;	ID=ORF_F_153700
contig_598	ID=g16576;	ID=ORF_R_164307
contig_621	ID=g16841;	ID=ORF_F_166062
contig_653	ID=g17169;	ID=ORF_F_169259
contig_757	ID=g18179;	ID=ORF_R_179006
contig_993	ID=g19772;	ID=ORF_F_193355
contig_1054     ID=g18249;      ID=ORF_R_196780
contig_1263     ID=g18988;      ID=ORF_F_203498
contig_1339     ID=g19176;      ID=ORF_R_205898
contig_1865     ID=g20046;      ID=ORF_F_213350
contig_1054	ID=g20072;	ID=ORF_R_196780
contig_1263	ID=g20935;	ID=ORF_F_203498
contig_1339	ID=g21173;	ID=ORF_R_205898
contig_1865	ID=g22216;	ID=ORF_F_213350
```

Some ORF effectors were noted to overlap AugustusCodingQuary gene models. These
were manually inspected in geneious and the following genes removed:

```bash

RemoveAugList="g553 g2191 g2212 g2443 g2890 g4524 g7303 g12540 g14477 g15058 g18179 g19772"
RemoveOrfList="ORF_F_12367 ORF_F_29498 ORF_R_35309 ORF_F_61898 ORF_R_84315 ORF_R_87287 ORF_F_89337 ORF_R_116876 ORF_R_133618 ORF_F_152796 ORF_F_153549 ORF_F_153700 ORF_R_164307 ORF_F_166062 ORF_F_169259 ORF_R_196780 ORF_F_203498 ORF_R_205898 ORF_F_213350"
echo "$RemoveAugList $RemoveOrfList" | sed 's/ /\n/g' > tmp.txt

OutDir=gene_pred/final_incl_ORF/P.cactorum/10300
cat $OutDir/10300_genes_incl_ORFeffectors.gff3  | grep -w 'gene' | wc -l
cat $OutDir/10300_genes_incl_ORFeffectors.gff3 | grep -w -v -f tmp.txt > $OutDir/10300_genes_incl_ORFeffectors_filtered.gff3
cat tmp.txt | wc -l
cat $OutDir/10300_genes_incl_ORFeffectors_filtered.gff3  | grep -w 'gene' | wc -l
```


```bash
  for GffAppended in $(ls gene_pred/final_incl_ORF/P.cactorum/10300/10300_genes_incl_ORFeffectors_filtered.gff3); do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    FinalDir=$(dirname $GffAppended)
    GffFiltered=$FinalDir/filtered_duplicates.gff
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
    GffRenamed=$FinalDir/final_genes_genes_incl_ORFeffectors_renamed.gff3
    LogFile=$FinalDir/final_genes_appended_renamed.log
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered
    Assembly=$(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_genes_incl_ORFeffectors_renamed
    # The proteins fasta file contains * instead of Xs for stop codons, these should
    # be changed
    sed -i 's/\*/X/g' $FinalDir/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta
  done
```

No duplicate genes were found.
```bash

```


# Assessing gene space in predicted transcriptomes

```bash
for Transcriptome in $(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.gene.fasta); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/genes
qsub $ProgDir/sub_busco3.sh $Transcriptome $BuscoDB $OutDir
done
```

```
INFO    281 Complete BUSCOs (C)
INFO    272 Complete and single-copy BUSCOs (S)
INFO    9 Complete and duplicated BUSCOs (D)
INFO    7 Fragmented BUSCOs (F)
INFO    15 Missing BUSCOs (M)
INFO    303 Total BUSCO groups searched
```

```bash
  for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Single=$(cat $File | grep "(S)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
  done
```



Make a list of BUSCO genes:

```bash
BuscoHits=$(ls gene_pred/busco/P.cactorum/414_v2/genes/run_final_genes_combined.gene/single_copy_busco_sequences/EOG0937*.fna)
OutDir=$(ls -d gene_pred/busco/P.cactorum/414_v2/genes/run_final_genes_combined.gene)
cat $BuscoHits | grep '>' | cut -f3 -d ':' > $OutDir/busco_single_copy_gene_headers.txt
```


# Re-running Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
  # mkdir -p gene_pred/final_genes/P.cactorum/414_v2/final_ncbi
  # gene_pred/final_ncbi/P.cactorum/414_v2/414_v2_genes_incl_ORFeffectors_* gene_pred/final_genes/P.cactorum/414_v2/final_ncbi/.
	screen -a
	cd /home/groups/harrisonlab/project_files/idris
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteome in $(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteome $InterProRaw
done
```


## B) SwissProt

```bash
  for Proteome in $(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    OutDir=gene_pred/swissprot/$Organism/$Strain
    SwissDbDir=../../uniprot/swissprot
    SwissDbName=uniprot_sprot
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
    qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
  done
```


# Additional effector candidates:

## Identifying homologs to previously characterised effectors

<!--
```bash
repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Query=$(ls analysis/blast_homology/oomycete_avr_genes/appended_oomycete_effectors_cds.fasta)
	qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
  for BlastHits in $(ls analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_avr_cds.fasta_hits.csv); do
    Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)  
    Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    HitsGff=analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_avr_cds.fasta_hits.gff
    Column2=Avr_homolog
    NumHits=10
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
  done
```

```bash
  GeneModels=$(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.gff3)
  BlastHits=$(ls analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_effector_cds.fasta_homologs.gff)
  bedtools intersect -wb -a $GeneModels -b $BlastHits | grep -w 'mRNA' | cut -f9,18 | tr -d ';' | tr -d '"' | sed 's/ID=//g' | sed -r "s/Parent=.*\t/\t/g" > analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_effector_cds.fasta_homologs.headers.txt

``` -->

```bash
qlogin -l h=blacklace11.blacklace
cd /home/groups/harrisonlab/project_files/idris
dbFasta=$(ls analysis/blast_homology/oomycete_avr_genes/appended_oomycete_effectors_cds.fasta)
dbType="nucl"
QueryFasta=$(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
Prefix="10300_appended_oomycete_effectors"
Eval="1e-30"
# WorkDir=$TMPDIR
OutDir=analysis/blast_homology/P.cactorum/10300
mkdir -p $OutDir
#-------------------------------------------------------
# 		Step 1.		Make database
#-------------------------------------------------------
makeblastdb -in $dbFasta -input_type fasta -dbtype nucl -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
#-------------------------------------------------------
# 		Step 2.		Blast Search
#-------------------------------------------------------
tblastx -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
#-------------------------------------------------------
# 		Step 3.		Summarise hits
#-------------------------------------------------------
cat $OutDir/${Prefix}_hits.txt | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
```

<!--
```bash
for Proteome in $(ls gene_pred/braker/P.cactorum/10300/P.cactorum/final_genes_Braker.cds.fasta ); do
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Query=$(ls analysis/blast_homology/oomycete_avr_genes/appended_oomycete_effectors_cds.fasta)
	qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
done
``` -->


## PhiBase genes:


Identifying PHIbase homologs
<!--
The PHIbase database was searched against the assembled genomes using tBLASTx. The default output location is:

```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
# Version 4.4 November 10th 2017
Version=v4.4
DbDir=$(ls -d ../../phibase/$Version)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/blast_pipe.sh $DbDir/phi_accessions.fa protein $Assembly
done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
  for BlastHits in $(ls analysis/blast_homology/P.cactorum/10300/10300_phi_accessions.fa_homologs.csv); do
    Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)  
    Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    HitsGff=analysis/blast_homology/P.cactorum/10300/10300_phi_accessions.fa_homologs_homologs.gff
    Column2=phibase_homolog
    NumHits=10
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
  done
```

```bash
  GeneModels=$(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.gff3)
  BlastHits=$(ls analysis/blast_homology/P.cactorum/10300/10300_phi_accessions.fa_homologs_homologs.gff)
  bedtools intersect -wb -a $GeneModels -b $BlastHits | grep -w 'mRNA' | cut -f9,18 | tr -d ';' | tr -d '"' | sed 's/ID=//g' | sed -r "s/Parent=.*\t/\t/g" > analysis/blast_homology/P.cactorum/10300/10300_phi_accessions.fa_homologs_homologs.headers.txt
```
 -->

```bash
qlogin -l h=blacklace11.blacklace
cd /home/groups/harrisonlab/project_files/idris
dbFasta=$(ls ../../phibase/v4.4/phi_accessions.fa)
dbType="prot"
QueryFasta=$(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
Prefix="10300_phi_accessions"
Eval="1e-30"
 # WorkDir=$TMPDIR
 OutDir=analysis/blast_homology/P.cactorum/10300
 mkdir -p $OutDir
 #-------------------------------------------------------
 # 		Step 1.		Make database
 #-------------------------------------------------------
 makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
 #-------------------------------------------------------
 # 		Step 2.		Blast Search
 #-------------------------------------------------------
 blastx -num_threads 20 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
 #-------------------------------------------------------
 # 		Step 3.		Summarise hits
 #-------------------------------------------------------
 cat $OutDir/${Prefix}_hits.txt | grep 'effector' | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
```

## Chen et al 2014 homologs:

All Chen 2014 genes

```bash
cd /home/groups/harrisonlab/project_files/idris
OutDir=analysis/blast_homology/chen_2014
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/GB/GX/GBGX01/GBGX01.1.fsa_nt.gz
gunzip GBGX01.1.fsa_nt.gz
mv GBGX01.1.fsa_nt 10300_transcripts_chen_2014.fasta
cd /home/groups/harrisonlab/project_files/idris

qlogin -l h=blacklace10.blacklace -pe smp 16
cd /home/groups/harrisonlab/project_files/idris
# dbFasta=$(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
dbFasta=$(ls analysis/blast_homology/chen_2014/10300_transcripts_chen_2014.fasta)
dbType="nucl"
# QueryFasta=$(ls analysis/blast_homology/chen_2014/10300_transcripts_chen_2014.fasta)
QueryFasta=$(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
Prefix="10300_chen_transcripts"
Eval="1e-30"
 # WorkDir=$TMPDIR
 OutDir=analysis/blast_homology/P.cactorum/10300
 mkdir -p $OutDir
 #-------------------------------------------------------
 # 		Step 1.		Make database
 #-------------------------------------------------------
 makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
 #-------------------------------------------------------
 # 		Step 2.		Blast Search
 #-------------------------------------------------------
 tblastx -num_threads 16 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
 #-------------------------------------------------------
 # 		Step 3.		Summarise hits
 #-------------------------------------------------------
 cat $OutDir/${Prefix}_hits.txt | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
```



Chen 2014 RxLRs
```bash
qlogin -pe smp 16
cd /home/groups/harrisonlab/project_files/idris
dbFasta=$(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
# dbFasta=$(ls analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.fa)
dbType="nucl"
QueryFasta=$(ls analysis/blast_homology/oomycete_avr_genes/chen_et_al_2014_RxLR.fa)
# QueryFasta=$(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
Prefix="10300_chen_RxLRs"
Eval="1e-30"
 # WorkDir=$TMPDIR
OutDir=analysis/blast_homology/P.cactorum/10300
mkdir -p $OutDir
#-------------------------------------------------------
# 		Step 1.		Make database
#-------------------------------------------------------
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
#-------------------------------------------------------
# 		Step 2.		Blast Search
#-------------------------------------------------------
tblastx -num_threads 16 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
#-------------------------------------------------------
# 		Step 3.		Summarise hits
#-------------------------------------------------------
cat $OutDir/${Prefix}_hits.txt | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
```


# Build annotation table


```bash
Organism=P.cactorum
Strain=10300
OutDir=analysis/gene_tables/$Organism/$Strain
mkdir -p $OutDir

# combine sigP2 headers
cat gene_pred/final_sigP/P.cactorum/10300/10300_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' > $OutDir/sigP2_headers_combined.txt
cat gene_pred/ORF_sigP/P.cactorum/10300/10300_ORF_sp_headers.txt | tr -d ' ' >> $OutDir/sigP2_headers_combined.txt

# combine sigP3 headers
cat gene_pred/final_signalp-3.0/P.cactorum/10300/10300_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' '  > $OutDir/sigP3_headers_combined.txt
cat gene_pred/ORF_signalp-3.0/P.cactorum/10300/10300_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' >> $OutDir/sigP3_headers_combined.txt

# combine sigP4 headers
cat gene_pred/final_signalp-4.1/P.cactorum/10300/10300_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' '  > $OutDir/sigP4_headers_combined.txt
cat gene_pred/ORF_signalp-4.1/P.cactorum/10300/10300_aug_sp.aa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' >> $OutDir/sigP4_headers_combined.txt

cat analysis/phobius/P.cactorum/10300/10300_phobius_headers.txt analysis/phobius/P.cactorum/10300/10300_phobius_headers_ORF.txt > $OutDir/10300_phobius_headers_combined.txt

# cat analysis/RxLR_effectors/RxLR_EER_regex_finder/P.cactorum/10300/10300_all_secreted_RxLR_regex.fa | grep '>' | cut -f1 | tf -d '>' > $OutDir/10300_RxLR_regex_non-ORF_headers.txt

cat analysis/RxLR_effectors/RxLR_EER_regex_finder/P.cactorum/10300/10300_RxLR_EER_regex.fa analysis/RxLR_effectors/RxLR_EER_regex_finder/P.cactorum/10300/10300_ORF_RxLR_EER_regex_merged.aa > $OutDir/10300_combined_RxLR_EER_regex.fa

cat analysis/RxLR_effectors/hmmer_RxLR/P.cactorum/10300/10300_RxLR_hmmer.fa analysis/RxLR_effectors/hmmer_RxLR/P.cactorum/10300/10300_ORF_RxLR_hmmer.fa > $OutDir/10300_combined_RxLR_hmmer.fa

cat analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain/10300_Aug_WY_hmmer.fa analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain/10300_ORF_WY_hmmer.fa > $OutDir/10300_combined_WY_hmmer.fa

cat analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORFsUniq_CRN_hmmer.bed \
	| grep -w 'transcript' | cut -f5 -d '=' | cut -f1 -d ';' \
	> analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORFsUniq_CRN_hmmer.txt
cat analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORF_CRN_LFLAK_unmerged_hmmer.fa \
	| grep -A1 -f analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORFsUniq_CRN_hmmer.txt \
	| grep -v '\-\-' \
	> analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORFsUniq_CRN_hmmer.fa

cat analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_pub_CRN_LFLAK_hmm.fa analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORFsUniq_CRN_hmmer.fa > $OutDir/10300_combined_CRN_LFLAK_hmm.fa

cat analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_pub_CRN_DWL_hmm.fa analysis/CRN_effectors/hmmer_CRN/P.cactorum/10300/10300_ORF_CRN_DWL_unmerged_hmmer.fa > $OutDir/10300_combined_CRN_DWL_hmm.fa

Organism=P.cactorum
Strain=10300
OutDir=analysis/gene_tables/$Organism/$Strain

GeneGff=gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.gff3
GeneFasta=gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta
OrfGff=$(ls gene_pred/ORF_finder/P.cactorum/10300/10300_ORF.gff3)
ConversionLog=$(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_appended_renamed.log)
SigP2=$(ls $OutDir/sigP2_headers_combined.txt)
SigP3=$(ls $OutDir/sigP3_headers_combined.txt)
SigP4=$(ls $OutDir/sigP4_headers_combined.txt)
PhobiusTxt=$(ls $OutDir/10300_phobius_headers_combined.txt)
TmhmmTxt=$(ls gene_pred/trans_mem/P.cactorum/10300/10300_TM_genes_pos_headers.txt)
GpiTxt=$(ls gene_pred/trans_mem/P.cactorum/10300/GPIsom/GPI_pos.txt)
RxLR_Motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/P.cactorum/10300/10300_all_secreted_RxLR_regex.fa)
RxLR_EER_Motif=$(ls $OutDir/10300_combined_RxLR_EER_regex.fa)
RxLR_Hmm=$(ls $OutDir/10300_combined_RxLR_hmmer.fa)
RxLR_WY=$(ls $OutDir/10300_combined_WY_hmmer.fa)
CRN_LFLAK=$(ls $OutDir/10300_combined_CRN_LFLAK_hmm.fa)
CRN_DWL=$(ls $OutDir/10300_combined_CRN_DWL_hmm.fa)
CAZY_Hmm=$(ls gene_pred/CAZY/*/*/*CAZY.out.dm.ps | grep '10300')
# RawCounts=alignment/star/P.cactorum/10300/DeSeq/V8_raw_counts.txt
# Fpkm=alignment/star/P.cactorum/10300/DeSeq/V8_fpkm_norm_counts.txt
PhiHits=$(ls analysis/blast_homology/P.cactorum/10300/10300_phi_accessions_hits_headers.txt)
AvrHits=$(ls analysis/blast_homology/P.cactorum/10300/10300_appended_oomycete_effectors_hits_headers.txt)
ChenHits=$(ls analysis/blast_homology/P.cactorum/10300/10300_chen_transcripts_hits_headers.txt)
InterPro=$(ls gene_pred/interproscan/final_incl_ORF/P.cactorum/P.cactorum_interproscan.tsv)
Swissprot=$(ls gene_pred/swissprot/P.cactorum/10300/swissprot_vJul2016_tophit_parsed.tbl)
OrthoName='Pcac'
OrthoAll='Pcac Pinf Ppar Pcap Psoj'
OrthoFile=$(ls analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_publication/formatted/Results_Jan09/Orthogroups.txt)
# OrthoFile=$(ls analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj_publication/*orthogroups.txt)
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/10300_analysis
$ProgDir/10300_gene_tables.py --gff_format gff3 --ortho_name $OrthoName --ortho_all $OrthoAll --ortho_file $OrthoFile --gene_gff $GeneGff --ORF_gff $OrfGff --gene_fasta $GeneFasta --conversion_list $ConversionLog --SigP2 $SigP2 --SigP3 $SigP3 --SigP4 $SigP4 --phobius $PhobiusTxt --TMHMM $TmhmmTxt --GPI $GpiTxt --RxLR_motif $RxLR_Motif --RxLR_EER_motif $RxLR_EER_Motif --RxLR_Hmm $RxLR_Hmm --RxLR_WY $RxLR_WY --CRN_LFLAK $CRN_LFLAK --CRN_DWL $CRN_DWL --CAZY $CAZY_Hmm --PhiHits $PhiHits --AvrHits $AvrHits --ChenHits $ChenHits --InterPro $InterPro --Swissprot $Swissprot > $OutDir/10300_gene_table_final.tsv
```

# Genomic analysis


Identification of NLP necrosis inducing proteins:

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | grep 'NPP1' | wc -l
cat $AnnotTab | grep 'NPP1' > $OutDir/10300_gene_table_NLP.tsv
cat $OutDir/10300_gene_table_NLP.tsv | grep 'secreted' | wc -l
cat $OutDir/10300_gene_table_NLP.tsv | grep '#NPP1#' | wc -l
cat $OutDir/10300_gene_table_NLP.tsv | grep '#NPP1#' | grep 'secreted' | wc -l
cat $OutDir/10300_gene_table_NLP.tsv | grep 'secreted' | grep "GBGX" | wc -l
cat $AnnotTab | grep -e 'AF356840' -e 'AF352031' > $OutDir/10300_gene_table_NLP_Avr_homologs.tsv

```

Proteases:

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
# Number of secreted proteases:
cat $AnnotTab | grep 'IPR001254' | wc -l
cat $AnnotTab | grep 'IPR001254' | grep 'secreted' | wc -l
cat $AnnotTab | grep 'IPR001254' | grep 'secreted' > $OutDir/10300_gene_table_protease.tsv
# Number with homologs in Chen 2014:
cat $OutDir/10300_gene_table_protease.tsv | grep 'GBGX' | wc -l
cat $OutDir/10300_gene_table_protease.tsv | grep 'GBGX' | cut -f8 | sort | uniq -c | wc -l
```

Protease inhibitors:

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | grep -i 'kazal' | wc -l
cat $AnnotTab | grep 'IPR002350' | wc -l
cat $AnnotTab | grep 'IPR002350' > $OutDir/10300_gene_table_EPI_kazal.tsv
cat $OutDir/10300_gene_table_EPI_kazal.tsv | grep 'secreted' | wc -l
cat $AnnotTab | grep -i 'cathepsin' | wc -l
cat $AnnotTab | grep 'IPR013201' | wc -l
cat $AnnotTab | grep 'IPR013201' > $OutDir/10300_gene_table_EPI_cathepsin.tsv
cat $OutDir/10300_gene_table_EPI_cathepsin.tsv | grep 'secreted' | wc -l
cat $AnnotTab | grep -e 'IPR000010' -e 'IPR027214' -e 'IPR018073' -e 'IPR020381' -e 'PF00031' -e '#EPIC' -e "EPIC._" | wc -l
cat $AnnotTab | grep -e 'IPR000010' -e 'IPR027214' -e 'IPR018073' -e 'IPR020381' -e 'PF00031' -e '#EPIC' -e "EPIC._"  > $OutDir/10300_gene_table_EPI_cysteine.tsv
cat $OutDir/10300_gene_table_EPI_cysteine.tsv | grep 'secreted' | wc -l
# Number identified by BLAST searches:
cat $AnnotTab | grep "GIP._P" | wc -l
cat $AnnotTab | grep 'IPR001314' | wc -l
cat $OutDir/10300_gene_table_protease.tsv | grep "GIP._P" | wc -l
cat $AnnotTab | grep "GIP._P" | grep 'secreted' | wc -l
cat $AnnotTab | grep -e 'PITG_13636' -e 'PITG_21456' > $OutDir/10300_gene_table_GIP.tsv

# Number with homologs in Chen 2014:
cat $OutDir/10300_gene_table_EPI_kazal.tsv | grep 'GBGX' | wc -l
cat $OutDir/10300_gene_table_EPI_kazal.tsv | grep 'GBGX' | cut -f8 | sort | uniq -c | wc -l
cat $OutDir/10300_gene_table_EPI_cathepsin.tsv | grep 'GBGX' | wc -l
cat $OutDir/10300_gene_table_EPI_cathepsin.tsv | grep 'GBGX' | cut -f8 | sort | uniq -c | wc -l
```


```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | grep 'PF09461' | wc -l
```

Elicitins:

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | grep -i 'elicitin' | wc -l
cat $AnnotTab | grep 'IPR002200' | wc -l
cat $AnnotTab | grep 'IPR002200' > $OutDir/10300_gene_table_elicitin.tsv
cat $AnnotTab | grep 'IPR002200' | grep 'secreted' | wc -l
cat $OutDir/10300_gene_table_elicitin.tsv | grep 'GBGX' | wc -l
cat $OutDir/10300_gene_table_elicitin.tsv | grep 'GBGX' | cut -f8 | sort | uniq -c | wc -l
```

TGases (R-glutaminyl-peptide:amine-γ-glutamyltransferase):

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | grep -i 'elicitor' | wc -l
cat $AnnotTab | grep 'IPR032048' | wc -l
cat $AnnotTab | grep 'IPR032048' > $OutDir/10300_gene_table_TGase.tsv
cat $AnnotTab | grep 'IPR032048' | grep 'secreted' | wc -l
cat $OutDir/10300_gene_table_TGase.tsv | grep 'GBGX' | wc -l
cat $OutDir/10300_gene_table_TGase.tsv | grep 'GBGX' | cut -f8 | sort | uniq -c | wc -l
```

RxLRs

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | grep -w 'RxLR' | tail -n+1 | wc -l
cat $AnnotTab | grep -w 'RxLR' > $OutDir/10300_gene_table_RxLR.tsv
cat $OutDir/10300_gene_table_RxLR.tsv | grep 'GBGX' | wc -l
cat $OutDir/10300_gene_table_RxLR.tsv | grep 'GBGX' | cut -f8 | sort | uniq -c | wc -l

cat $OutDir/10300_gene_table_RxLR.tsv | cut -f1,20 | grep 'Yes' | wc -l
cat $OutDir/10300_gene_table_RxLR.tsv | cut -f1,35 | sort | uniq -c | sort -nr
cat $OutDir/10300_gene_table_RxLR.tsv | cut -f1,25,26,27,31,35 | grep ";" > $OutDir/RxLR_notable_annotations.tsv
cat $OutDir/10300_gene_table_RxLR.tsv | grep -w 'Pcac(1)' | grep -w 'Pinf(1)' | wc -l
cat $OutDir/10300_gene_table_RxLR.tsv | grep -w 'Pcac(1)' | grep -w 'Psoj(1)' | wc -l
cat $OutDir/10300_gene_table_RxLR.tsv | grep -w 'Pcac(1)' | grep -w 'Ppar(1)' | wc -l
cat $OutDir/10300_gene_table_RxLR.tsv | grep -w 'Pcac(1)' | grep -w 'Pcap(1)' | wc -l

cat $OutDir/10300_gene_table_RxLR.tsv | cut -f25,26,27 | sort | uniq -c | sort -nr | less -S

cat $OutDir/10300_gene_table_RxLR.tsv | grep -i 'Avr' | cut -f1,25,26,27,30,31,35 | wc -l
cat $OutDir/10300_gene_table_RxLR.tsv | cut -f1,30 | grep -i 'Avr' | wc -l
```

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
Avr1='PITG_16663'
cat $AnnotTab | grep 'PITG_16663T' | sed "s/^/P. infestans Avr1 ortholog\t/g" > $OutDir/RxLR_orthologs.tsv
Avr2='PITG_05121 PITG_06077 PITG_07500 PITG_07499 PITG_08278 PITG_22870 PITG_08943 PITG_23009 PITG_13956 PITG_23008 PITG_13940 PITG_13930 PITG_20025 PITG_21645 PITG_21949'
cat $AnnotTab | grep -e 'PITG_05121T' -e 'PITG_06077T' -e 'PITG_07500T' -e 'PITG_07499T' -e 'PITG_08278T' -e 'PITG_22870T' -e 'PITG_08943T' -e 'PITG_23009T' -e 'PITG_13956T' -e 'PITG_23008T' -e 'PITG_13940T' -e 'PITG_13930T' -e 'PITG_20025T' -e 'PITG_21645T' -e 'PITG_21949' | sed "s/^/P. infestans Avr2 ortholog\t/g" >> $OutDir/RxLR_orthologs.tsv
Avr3a='PITG_14368 PITG_14371 PITG_14374'
cat $AnnotTab | grep -e 'PITG_14368T' -e 'PITG_14371T' -e 'PITG_14374 ' | sed "s/^/P. infestans Avr3a ortholog\t/g" >> $OutDir/RxLR_orthologs.tsv
Avr3bHomolog='PITG_15732'
cat $AnnotTab | grep 'PITG_15732T' | sed "s/^/P. infestans Avr3b ortholog\t/g" >> $OutDir/RxLR_orthologs.tsv
Avr4='PITG_07387'
cat $AnnotTab | grep 'PITG_07387T' | sed "s/^/P. infestans Avr4 ortholog\t/g" >> $OutDir/RxLR_orthologs.tsv
Avrblb1='PITG_21388'
cat $AnnotTab | grep 'PITG_21388T' | sed "s/^/P. infestans Avrblb1 ortholog\t/g" >> $OutDir/RxLR_orthologs.tsv
Avrblb2='PITG_04085 PITG_04086 PITG_20300 PITG_04090 PITG_20303 PITG_20301 PITG_18683'
cat $AnnotTab | grep -e 'PITG_04085T' -e 'PITG_04086T' -e 'PITG_20300T' -e 'PITG_04090T' -e 'PITG_20303T' -e 'PITG_20301T' -e 'PITG_18683T' | sed "s/^/P. infestans Avrblb2 ortholog\t/g" >> $OutDir/RxLR_orthologs.tsv
```


Crinklers:

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | grep -w 'CRN' | wc -l
cat $AnnotTab | grep -w 'CRN' > $OutDir/10300_gene_table_CRN.tsv

echo "The number of crinklers with SignalP predictions:"
cat $OutDir/10300_gene_table_CRN.tsv
echo "The number of crinklers with phobius predictions:"

cat $OutDir/10300_gene_table_CRN.tsv | grep 'GBGX' | wc -l
cat $OutDir/10300_gene_table_CRN.tsv | grep 'GBGX' | cut -f8 | sort | uniq -c | wc -l

```

CAZY proteins:

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | grep 'CAZY' | tail -n+2 | wc -l
# cat $AnnotTab | cut -f1,12,13,14,24,31 | grep 'CAZY' > $OutDir/10300_gene_table_CAZY.tsv
cat $AnnotTab | grep 'CAZY' > $OutDir/10300_gene_table_CAZY.tsv
cat $AnnotTab | grep 'CAZY' | cut -f12,13,14 | grep 'Yes' | tail -n+2 | wc -l
cat $AnnotTab | grep 'CAZY' | cut -f12,13,14,27 | grep 'Yes' | tail -n+2 > $OutDir/10300_gene_table_CAZY_secreted.tsv
cat $AnnotTab | grep 'CAZY' | grep 'secreted' > $OutDir/10300_gene_table_CAZY_secreted.tsv

cat $OutDir/10300_gene_table_CAZY_secreted.tsv | grep 'GBGX' | wc -l
cat $OutDir/10300_gene_table_CAZY_secreted.tsv | grep 'GBGX' | tail -n+2 | cut -f '34' | sort | uniq -c | wc -l

# cat $OutDir/10300_gene_table_CAZY.tsv | cut -f27 | sort | uniq -c | sort -nr > $OutDir/10300_gene_table_CAZY_hmm_models.txt
# cat $AnnotTab | grep 'CAZY' | cut -f12,13,14,27 | grep 'Yes' | tail -n+2 | cut -f4 | sort | uniq -c | sort -nr > $OutDir/10300_gene_table_CAZY_hmm_models.txt
cat $OutDir/10300_gene_table_CAZY_secreted.tsv | cut -f27 | sort | uniq -c | sort -nr > $OutDir/10300_gene_table_CAZY_hmm_models.txt
cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep 'GH'
cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep -v 'GH' | grep 'CBM'
cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep 'AA'
cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep 'CE'
cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep -v 'GH' | grep 'GT'
cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep 'PL'
```

Cutinases:

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | grep 'PF01083' | wc -l
cat $AnnotTab | grep 'PF01083' | grep 'secreted' | wc -l
```

Avr gene homologs:

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | cut -f1,30 | grep '_' | wc -l
cat $AnnotTab | cut -f1,19,23,25,26,30,35 | grep -P "_.*\t" > $OutDir/10300_gene_table_avr.tsv
```



# Determine lengths of intergenic regions

```bash
AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=analysis/intergenic_regions/P.cactorum/10300
mkdir -p $OutDir

GeneGff=$(ls gene_pred/final_incl_ORF/P.cactorum/10300/final_genes_genes_incl_ORFeffectors_renamed.gff3)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/intergenic_regions
$ProgDir/find_intergenic_regions.py --Gff $GeneGff > $OutDir/10300_intergenic_regions.txt

# RxLR IG distance plots and test
RxLR_list=$OutDir/RxLR_headers.txt
cat $AnnotTab | grep -w 'RxLR' | cut -f1 | cut -f1 -d '.' | sort | uniq > $RxLR_list
cat $OutDir/10300_intergenic_regions.txt | grep -w -f $RxLR_list > $OutDir/10300_RxLR_intergenic_regions.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/10300_analysis
$ProgDir/plot_intergenic_regions_with_effectors.r --inp_all $OutDir/10300_intergenic_regions.txt --inp_subset $OutDir/10300_RxLR_intergenic_regions.txt --prefix $OutDir/10300_RxLR

#[1] 0
#[1] 2
#[1] 0

# CRN IG distance plots and test
CRN_list=$OutDir/CRN_headers.txt
cat $AnnotTab | grep -w 'CRN' | cut -f1 | cut -f1 -d '.' | sort | uniq > $CRN_list
cat $OutDir/10300_intergenic_regions.txt | grep -w -f $CRN_list > $OutDir/10300_CRN_intergenic_regions.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/10300_analysis
# $ProgDir/plot_intergenic_regions_with_effectors.r --inp_all $OutDir/10300_intergenic_regions.txt --inp_subset $OutDir/10300_CRN_intergenic_regions.txt --prefix $OutDir/10300_CRN

#[1] 5635
#[1] 25
#[1] 326

# Elicitin IG distance plots and test
Effector="elicitin"
Headers=$OutDir/${Effector}_headers.txt
cat $AnnotTab | grep 'IPR002200' | cut -f1,12,13,14 | grep 'Yes' | cut -f1 -d '.' | sort | uniq > $Headers
cat $OutDir/10300_intergenic_regions.txt | grep -w -f $Headers > $OutDir/10300_${Effector}_intergenic_regions.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/10300_analysis
# $ProgDir/plot_intergenic_regions_with_effectors.r --inp_all $OutDir/10300_intergenic_regions.txt --inp_subset $OutDir/10300_${Effector}_intergenic_regions.txt --prefix $OutDir/10300_${Effector}

# [1] 9774
# [1] 9699
# [1] 9947


# Protease Inhibitor IG distance plots and test
Effector="EPI"
Headers=$OutDir/${Effector}_headers.txt
# cat $AnnotTab | grep -e 'IPR002350' -e 'IPR013201' | cut -f1,12,13,14 | grep 'Yes' | cut -f1 | cut -f1 -d '.' | sort | uniq > $Headers
cat $AnnotTab | grep -e 'IPR002350' -e 'IPR013201' -e 'IPR000010' -e 'IPR027214' -e 'IPR018073' -e 'IPR020381' -e 'PF00031' -e '#EPIC' -e'EPIC._' -e 'IPR001314' -e 'IPR032048' | cut -f1,12,13,14 | grep 'Yes' | cut -f1 | cut -f1 -d '.' | sort | uniq > $Headers
cat $OutDir/10300_intergenic_regions.txt | grep -w -f $Headers > $OutDir/10300_${Effector}_intergenic_regions.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/10300_analysis
# $ProgDir/plot_intergenic_regions_with_effectors.r --inp_all $OutDir/10300_intergenic_regions.txt --inp_subset $OutDir/10300_${Effector}_intergenic_regions.txt --prefix $OutDir/10300_${Effector}

#[1] 7649
#[1] 5089
#[1] 7085

# NLP IG distance plots and test
Effector="NLP"
Headers=$OutDir/${Effector}_headers.txt
cat $AnnotTab | grep -e 'NPP1' -e 'AF356840' -e 'AF352031' | cut -f1,12,13,14 | grep 'Yes' | cut -f1 | cut -f1 -d '.' | sort | uniq > $Headers
cat $OutDir/10300_intergenic_regions.txt | grep -w -f $Headers > $OutDir/10300_${Effector}_intergenic_regions.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/10300_analysis
# $ProgDir/plot_intergenic_regions_with_effectors.r --inp_all $OutDir/10300_intergenic_regions.txt --inp_subset $OutDir/10300_${Effector}_intergenic_regions.txt --prefix $OutDir/10300_${Effector}

# [1] 3689
# [1] 2267
# [1] 3019

# CAZY IG distance plots and test
Effector="CAZY"
Headers=$OutDir/${Effector}_headers.txt
cat $AnnotTab | grep -w -e 'CAZY' | cut -f1,12,13,14 | grep 'Yes' | cut -f1 | cut -f1 -d '.' | sort | uniq > $Headers
cat $OutDir/10300_intergenic_regions.txt | grep -w -f $Headers > $OutDir/10300_${Effector}_intergenic_regions.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/10300_analysis
# $ProgDir/plot_intergenic_regions_with_effectors.r --inp_all $OutDir/10300_intergenic_regions.txt --inp_subset $OutDir/10300_${Effector}_intergenic_regions.txt --prefix $OutDir/10300_${Effector}

#[1] 30
#[1] 2264
#[1] 81

Effector="CAZY_non_sec"
Headers=$OutDir/${Effector}_headers.txt
cat $AnnotTab | grep -w -e 'CAZY' | cut -f1,12,13,14 | grep -v 'Yes' | cut -f1 | cut -f1 -d '.' | sort | uniq > $Headers
cat $OutDir/10300_intergenic_regions.txt | grep -w -f $Headers > $OutDir/10300_${Effector}_intergenic_regions.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/10300_analysis
# $ProgDir/plot_intergenic_regions_with_effectors.r --inp_all $OutDir/10300_intergenic_regions.txt --inp_subset $OutDir/10300_${Effector}_intergenic_regions.txt --prefix $OutDir/10300_${Effector}

#[1] 9970
#[1] 9978
#[1] 9997


Effector="Busco"
Headers=$OutDir/${Effector}_headers.txt
cat gene_pred/busco/final_incl_ORF/P.cactorum/genes/run_final_genes_genes_incl_ORFeffectors_renamed.gene/single_copy_busco_sequences/*.fna | grep '>' | cut -f3 -d ':'| sort | uniq > $Headers
cat $OutDir/10300_intergenic_regions.txt | grep -w -f $Headers > $OutDir/10300_${Effector}_intergenic_regions.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/10300_analysis
# $ProgDir/plot_intergenic_regions_with_effectors.r --inp_all $OutDir/10300_intergenic_regions.txt --inp_subset $OutDir/10300_${Effector}_intergenic_regions.txt --prefix $OutDir/10300_${Effector}

#[1] 9999
#[1] 9999
#[1] 9999


# IG plot for publication

AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
OutDir=analysis/intergenic_regions/P.cactorum/10300

rxlr=$(ls $OutDir/10300_RxLR_intergenic_regions.txt)
crn=$(ls $OutDir/10300_CRN_intergenic_regions.txt)
cazy_sec=$(ls $OutDir/10300_CAZY_intergenic_regions.txt)
pi=$(ls $OutDir/10300_EPI_intergenic_regions.txt)
nlp=$(ls $OutDir/10300_NLP_intergenic_regions.txt)
cazy_nonsec=$(ls $OutDir/10300_CAZY_non_sec_intergenic_regions.txt)
busco=$(ls $OutDir/10300_Busco_intergenic_regions.txt)
elicitin=$(ls $OutDir/10300_elicitin_intergenic_regions.txt)

cat $OutDir/10300_intergenic_regions.txt | sed "s/$/\tTotal/g"  | cut -f1,2,3,5 > $OutDir/10300_IG_matrix_table.txt
for File in $(ls $OutDir/10300_*_intergenic_regions.txt); do
  Filename=$(basename $File)
  Effector=$(echo $Filename | sed 's/10300_//g' | sed 's/_intergenic_regions.txt//g' | sed 's/elicitin/Elicitin/g' | sed 's/CAZY_non_sec/Non-secreted CAZyme/g' | sed 's/CAZY/Secreted CAZyme/g' | sed 's/Busco/BUSCO/g' | sed 's/EPI/PI/g')
  echo $File
  echo $Effector
  cat $File | sed "s/$/\t$Effector/g"  | cut -f1,2,3,5 >> $OutDir/10300_IG_matrix_table.txt
done

# Creat a table of the 'Total' values repeated for each effector to act as a
# background for the IG density plots
cat $OutDir/10300_intergenic_regions.txt | sed "s/$/\tTotal/g"  | cut -f1,2,3,5 > $OutDir/10300_IG_matrix_table_for_plot_background.txt
for File in $(ls $OutDir/10300_*_intergenic_regions.txt); do
  Filename=$(basename $File)
  Effector=$(echo $Filename | sed 's/10300_//g' | sed 's/_intergenic_regions.txt//g' | sed 's/elicitin/Elicitin/g' | sed 's/CAZY_non_sec/Non-secreted CAZyme/g' | sed 's/CAZY/Secreted CAZyme/g' | sed 's/Busco/BUSCO/g' | sed 's/EPI/PI/g')
  echo $File
  echo $Effector
  cat $OutDir/10300_intergenic_regions.txt | sed "s/$/\t$Effector/g"  | cut -f1,2,3,5 >> $OutDir/10300_IG_matrix_table_for_plot_background.txt
done

# Scratch=/data/scratch/armita
# OutDir=$Scratch/analysis/intergenic_regions/P.cactorum/10300
# mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/10300_analysis
# $ProgDir/plot_intergenic_regions_with_effectors_publication.r --inp_all $OutDir/10300_intergenic_regions.txt --rxlr_subset $rxlr --crn_subset $crn --cazy_sec_subset $cazy_sec --PI_subset $pi --busco_subset $busco --cazy_nonsec_subset $cazy_nonsec --elicitin_subset $elicitin --prefix $OutDir/10300_IG_matrix_plot
$ProgDir/plot_intergenic_regions_with_effectors_publication2.r --inp_all $OutDir/10300_IG_matrix_table.txt --background_table $OutDir/10300_IG_matrix_table_for_plot_background.txt --prefix $OutDir/10300_IG_matrix_plot



AnnotTab=$(ls analysis/gene_tables/P.cactorum/10300/10300_gene_table_final.tsv)
IG_file=$(ls analysis/intergenic_regions/P.cactorum/10300/10300_intergenic_regions.txt)


cat $AnnotTab | grep -w -C1 "CRN" | 'CRN'



for Effector in 'CAZY' 'CRN'; do
EffectorList=$(ls analysis/intergenic_regions/P.cactorum/10300/10300_${Effector}_intergenic_regions.txt)
OutDir=analysis/intergenic_regions/P.cactorum/10300/${Effector}
mkdir -p $OutDir
cat analysis/intergenic_regions/P.cactorum/10300/10300_${Effector}_intergenic_regions.txt | grep '+' > $OutDir/10300_${Effector}_intergenic_regions_+.txt
cat $IG_file | grep -B1 -f $OutDir/10300_${Effector}_intergenic_regions_+.txt | awk 'NR%3==2' | cut -f1 > $OutDir/10300_${Effector}_intergenic_regions_+_fiveprime_neighbour.txt
cat analysis/intergenic_regions/P.cactorum/10300/10300_${Effector}_intergenic_regions.txt | grep '-' > $OutDir/10300_${Effector}_intergenic_regions_-.txt
cat $IG_file | grep -A1 -f $OutDir/10300_${Effector}_intergenic_regions_-.txt | awk 'NR%3==2' | cut -f1 > $OutDir/10300_${Effector}_intergenic_regions_-_fiveprime_neighbour.txt
cat $OutDir/10300_${Effector}_intergenic_regions_+_fiveprime_neighbour.txt $OutDir/10300_${Effector}_intergenic_regions_-_fiveprime_neighbour.txt > $OutDir/10300_${Effector}_intergenic_regions_fiveprime_neighbour.txt

cat $AnnotTab | grep -w -f $OutDir/10300_${Effector}_intergenic_regions_fiveprime_neighbour.txt > $OutDir/10300_${Effector}_intergenic_regions_fiveprime_neighbour_annot.txt



cat analysis/intergenic_regions/P.cactorum/10300/10300_${Effector}_intergenic_regions.txt | grep '-' > $OutDir/10300_${Effector}_intergenic_regions_-.txt
cat $IG_File | grep -B1 -f $OutDir/10300_${Effector}_intergenic_regions_-.txt | awk 'NR%3==1' | cut -f1 > $OutDir/10300_${Effector}_intergenic_regions_-_fiveprime_neighbour.txt
cat analysis/intergenic_regions/P.cactorum/10300/10300_${Effector}_intergenic_regions.txt | grep '+' > $OutDir/10300_${Effector}_intergenic_regions_+.txt
cat $IG_file | grep -A1 -f $OutDir/10300_${Effector}_intergenic_regions_+.txt | awk 'NR%3==2' | cut -f1 > $OutDir/10300_${Effector}_intergenic_regions_+_fiveprime_neighbour.txt
cat $OutDir/10300_${Effector}_intergenic_regions_-_fiveprime_neighbour.txt $OutDir/10300_${Effector}_intergenic_regions_+_fiveprime_neighbour.txt > $OutDir/10300_${Effector}_intergenic_regions_fiveprime_neighbour.txt

cat $AnnotTab | grep -w -f $OutDir/10300_${Effector}_intergenic_regions_fiveprime_neighbour.txt > $OutDir/10300_${Effector}_intergenic_regions_fiveprime_neighbour_annot.txt

done

cat $OutDir/10300_CRN_intergenic_regions.txt | grep '-' > $OutDir/10300_CRN_intergenic_regions_-.txt


```
