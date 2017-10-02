
# Strucutral variant discovery:

Structural variant discovery in P.cactorum and P. idae genomes.

Work was performed in the following directory:

```bash
OutDir=analysis/popgen/indel_calling
mkdir -p $OutDir

```

# Extension of the GATK SNP discovery pipeline


Filtering was performed to remove all SNPs from the vcf file produced in 01 and only retain indels:

```bash
Vcf=$(ls analysis/popgen/SNP_calling/414_v2_contigs_unmasked.vcf)
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
qsub $ProgDir/sub_vcf_parser_only_indels.sh $Vcf 40 30 10 30 1
```

```bash
OutDir=analysis/popgen/indel_calling
mv 414_v2_contigs_unmasked_filtered_only_indels.vcf $OutDir/414_v2_contigs_unmasked_only_indels.vcf
```

## Remove sequencing errors from vcf files:

```bash
Vcf=$(ls analysis/popgen/indel_calling/414_v2_contigs_unmasked_only_indels.vcf)
OutDir=$(dirname $Vcf)
Errors=$OutDir/414_error_SNPs.tsv
FilteredVcf=$OutDir/414_v2_contigs_unmasked_filtered_no_indels_no_errors.vcf
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
$ProgDir/flag_error_SNPs.py --inp_vcf $Vcf --ref_isolate 414 --errors $Errors --filtered $FilteredVcf
echo "The number of probable errors from homozygous SNPs being called from reference illumina reads vs the reference assembly is:"
cat $Errors | wc -l
echo "These have been removed from the vcf file"
```

```
  355
```


## Collect VCF stats

General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  VcfTools=/home/sobczm/bin/vcftools/bin
  export PERL5LIB="$VcfTools:$PERL5LIB"
  VcfFiltered=$(ls analysis/popgen/indel_calling/414_v2_contigs_unmasked_filtered_no_indels_no_errors.vcf)
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```


# Identify SNPs in gene models:

Create custom SnpEff genome database

```bash
SnpEff=/home/sobczm/bin/snpEff
nano $SnpEff/snpEff.config
```


Add the following lines to the section with databases:

```
#---
# EMR Databases
#----
# Fus2 genome
Fus2v1.0.genome : Fus2
# Bc16 genome
Bc16v1.0.genome: BC-16
# P414 genome
P414v1.0.genome: 414
```

Collect input files

```bash
Reference=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa)
Gff=$(ls gene_pred/final_ncbi/P.cactorum/414_v2/final_ncbi/414_v2_genes_incl_ORFeffectors_renamed.gff3)
SnpEff=/home/sobczm/bin/snpEff
mkdir $SnpEff/data/P414v1.0
cp $Reference $SnpEff/data/P414v1.0/sequences.fa
cp $Gff $SnpEff/data/P414v1.0/genes.gff

#Build database using GFF3 annotation
java -jar $SnpEff/snpEff.jar build -gff3 -v P414v1.0
```


## Annotate VCF files
```bash
CurDir=/home/groups/harrisonlab/project_files/idris
cd $CurDir
for a in $(ls analysis/popgen/indel_calling/414_v2_contigs_unmasked_filtered_no_indels_no_errors.vcf); do
    echo $a
    filename=$(basename "$a")
    Prefix=${filename%.vcf}
    OutDir=$(ls -d analysis/popgen/indel_calling)
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $a > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
    # mv 414_v2_contigs_unmasked_filtered* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf
    AllSnps=$(cat $OutDir/"$Prefix"_annotated.vcf | grep -v '#' | wc -l)
    GeneSnps=$(cat $OutDir/"$Prefix"_gene.vcf | grep -v '#' | wc -l)
    CdsSnps=$(cat $OutDir/"$Prefix"_coding.vcf | grep -v '#' | wc -l)
    NonsynSnps=$(cat $OutDir/"$Prefix"_nonsyn.vcf | grep -v '#' | wc -l)
    SynSnps=$(cat $OutDir/"$Prefix"_syn.vcf | grep -v '#' | wc -l)
    #-
    # SNPs in effectors
    #-
    AnnotaTable=$(ls gene_pred/annotation/P.cactorum/414_v2/414_v2_gene_table_incl_exp.tsv)
    Busco=$(ls gene_pred/busco/P.cactorum/414_v2/genes/run_final_genes_combined.gene/busco_single_copy_gene_headers.txt)
    RxLR=$(ls gene_pred/annotation/P.cactorum/414_v2/renamed_RxLR.txt)
    CRN=$(ls analysis/CRN_effectors/hmmer_CRN/P.cactorum/414_v2/414_v2_final_CRN_ID.txt)
    cat $AnnotaTable | cut -f1,12 | tail -n+2 | grep 'Yes' | cut -f1 > $RxLR
    cat $AnnotaTable | cut -f1,13 | tail -n+2 | grep 'Yes' | cut -f1 > $CRN
    #-
    # syn SNPs in effectors:
    #-
    SynVcf=$OutDir/"$Prefix"_syn.vcf
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
    BuscoOut=$OutDir/"$Prefix"_syn_Busco.vcf
    $ProgDir/vcf_extract_genes.py --vcf $SynVcf --gene_list $Busco > $BuscoOut
    BuscoSynSnps=$(cat $BuscoOut | grep -v '#' | wc -l)
    RxlrOut=$OutDir/"$Prefix"_syn_RxLR.vcf
    $ProgDir/vcf_extract_genes.py --vcf $SynVcf --gene_list $RxLR > $RxlrOut
    RxlrSynSnps=$(cat $RxlrOut | grep -v '#' | wc -l)
    CrnOut=$OutDir/"$Prefix"_CRN.vcf
    $ProgDir/vcf_extract_genes.py --vcf $SynVcf --gene_list $CRN > $CrnOut
    CrnSynSnps=$(cat $CrnOut | grep -v '#' | wc -l)  
    #-
    # non-syn SNPs in effectors:
    #-
    NonSynVcf=$OutDir/"$Prefix"_nonsyn.vcf
    BuscoOut=$OutDir/"$Prefix"_nonsyn_Busco.vcf
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
    $ProgDir/vcf_extract_genes.py --vcf $NonSynVcf --gene_list $Busco > $BuscoOut
    BuscoNonSynSnps=$(cat $BuscoOut | grep -v '#' | wc -l)
    RxlrOut=$OutDir/"$Prefix"_nonsyn_RxLR.vcf
    $ProgDir/vcf_extract_genes.py --vcf $NonSynVcf --gene_list $RxLR > $RxlrOut
    RxlrNonSynSnps=$(cat $RxlrOut | grep -v '#' | wc -l)
    CrnOut=$OutDir/"$Prefix"_nonsyn_CRN.vcf
    $ProgDir/vcf_extract_genes.py --vcf $NonSynVcf --gene_list $CRN > $CrnOut
    CrnNonSynSnps=$(cat $CrnOut | grep -v '#' | wc -l)
    printf "Comparison\tAllSnps\tGeneSnps\tCdsSnps\tSynSnps\tNonsynSnps\tBuscoSynSnps\tBuscoNonSynSnps\tRxlrSynSnps\tRxlrNonSynSnps\tCrnSynSnps\tCrnNonSynSnps\n"
    printf "$Prefix\t$AllSnps\t$GeneSnps\t$CdsSnps\t$SynSnps\t$NonsynSnps\t$BuscoSynSnps\t$BuscoNonSynSnps\t$RxlrSynSnps\t$RxlrNonSynSnps\t$CrnSynSnps\t$CrnNonSynSnps\n"

    #-
    # Make venn diagrams
    # -
    # These relate to the number of SNPs differing from the P414
    # reference in eahc group
for Vcf in $(ls analysis/popgen/indel_calling/*_nonsyn*.vcf | grep -v -e 'recode' -e '.vcf_'); do
Prefix=$(echo $Vcf | sed 's/.vcf//g')
Group1="12420 15_13 15_7 2003_3 4032 404 4040 414 415 416 62471"
Group2="PC13_15 P295 R36_14"
Group3="371 SCRP370 SCRP376"
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
$ProgDir/vcf_2_venn.py --vcf $Vcf --g1_name Pc_Fxa --g1_isolates $Group1 --g2_name Pc_Mxd --g2_isolates $Group2 --g3_name Pi_Ri --g3_isolates $Group3 --prefix $Prefix
done
done
```

```
Comparison  AllSnps	GeneSnps	CdsSnps	SynSnps	NonsynSnps	BuscoSynSnps	BuscoNonSynSnps	RxlrSynSnps	RxlrNonSynSnps	CrnSynSnps	CrnNonSynSnps
414_v2_contigs_unmasked_filtered_no_indels_no_errors	26510	2723	4	0	4	0	0	0	0	0	0
```

The four genes were identified as:

```
cat gene_pred/annotation/P.cactorum/414_v2/414_v2_gene_table_incl_exp.tsv | grep -w -e 'g7540' -e 'g9518' -e 'g20514' -e 'g22621' | less -S
```

# Lumpy structural variant discovery:

 lumpy is specifically geared toward identifying indels and uses more lines of evidence. Also, GATK will detect only small indels, whereas lumpy can detect bigger structural variants.

Run BWA-mem


```bash
CurDir=$PWD
Reference=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa)
for StrainPath in $(ls -d qc_dna/paired/P.*/* | grep -e 'P.idaei' -e 'P.cactorum'|  grep -v '10300' | grep '12420'); do
  Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
  Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
  echo $Strain
  echo $Organism
    ReadsF=$(ls $StrainPath/F/*fq.gz)
    ReadsR=$(ls $StrainPath/R/*fq.gz)
    ConcatTmpDir=tmp_concat_dir
    mkdir -p $ConcatTmpDir
    ConcatF=$ConcatTmpDir/"$Strain"_F_reads.fq.gz
    ConcatR=$ConcatTmpDir/"$Strain"_R_reads.fq.gz
    cat $ReadsF > $ConcatF
    cat $ReadsR > $ConcatR
    OutDir=analysis/popgen/indel_calling/alignments2
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
    qsub $ProgDir/sub_prep_lumpy.sh $Strain $CurDir/$Reference $ConcatF $ConcatR $OutDir
done
```

```bash
  ConcatTmpDir=tmp_concat_dir
  rm -r $ConcatTmpDir
  OutDir=analysis/popgen/indel_calling/alignments
  mkdir -p $OutDir
  mv *_rg.bam* $OutDir/.
  mv *_rg_discordants.bam*  $OutDir/.
  mv *_rg_splitters.bam* $OutDir/.
```

```bash
  CurDir=$PWD
  AlignDir=analysis/popgen/indel_calling/alignments
  cd $AlignDir
  for Strain in $(ls *_sorted.bam | grep -v -e 'split' -e 'discordant' | sed 's/_sorted.bam//g'); do
    # Jobs=$(qstat | grep 'sub_lumpy' | wc -l)
    # while [ $Jobs -gt 0 ]; do
    # printf "."
    # sleep 10s
    # Jobs=$(qstat | grep 'sub_lumpy' | wc -l)
    # done		
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
  # qsub $ProgDir/sub_lumpy.sh "$Strain"_lumpy $Strain
  $ProgDir/sub_lumpy.sh "$Strain"_lumpy $Strain
  done
  cd $CurDir
```
