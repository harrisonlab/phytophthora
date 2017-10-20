

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
