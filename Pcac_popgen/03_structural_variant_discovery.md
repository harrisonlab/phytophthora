
# Strucutral variant discovery:

Structural variant discovery in P.cactorum and P. idae genomes.

.vcf file formats are complicated. Documentations should be read to understand
outputs from SNP and SV calling:
https://samtools.github.io/hts-specs/VCFv4.2.pdf



Work was performed in the following directory:

```bash
OutDir=analysis/popgen/indel_calling
mkdir -p $OutDir

```
<!--
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
    # java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
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
``` -->

# structural variant discovery:

 lumpy and svaba are specifically geared toward identifying indels and structural variants and uses more lines of evidence that gatk. Also, GATK will detect only small indels, whereas lumpy and svaba can detect bigger structural variants. Lumpy is an alignment based SV caller and svaba is an assembly based SV caller.


<!--
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
  AlignDir=analysis/popgen/indel_calling/alignments2
  cd $AlignDir
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
  qsub $ProgDir/sub_lumpy.sh lumpy_out


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
``` -->






For this analysis svaba was used.
<!--
Run BWA-mem

```bash
CurDir=$PWD
Reference=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa)
for StrainPath in $(ls -d qc_dna/paired/P.*/* | grep -e 'P.idaei' -e 'P.cactorum' | grep -v '10300'); do
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
# cat $ReadsF > $ConcatF
# cat $ReadsR > $ConcatR
OutDir=analysis/popgen/indel_calling/alignments3
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
qsub $ProgDir/sub_bwa.sh $Strain $CurDir/$Reference $ConcatF $ConcatR $OutDir
done
``` -->

Alignments from SNP calling commands (detailed in 01_Pcac_alignments.md) were used.

Symbolic links were made, putting a link to each alignment in a single directory.

```bash
  OutDir=analysis/popgen/indel_calling/alignments
  mkdir -p $OutDir
  for File in $(ls analysis/genome_alignment/bowtie/*/*/vs_414/*sorted); do
    Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    cp -s $PWD/$File $OutDir/"$Strain"_vs_414_aligned_sorted.bam
    samtools index -@ 4 $OutDir/"$Strain"_vs_414_aligned_sorted.bam
  done
```

```bash
  Prefix=Pcac_svaba
  Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  AlignDir=analysis/popgen/indel_calling/alignments
  OutDir=analysis/popgen/indel_calling/svaba
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
  qsub $ProgDir/sub_svaba.sh $Prefix $Reference $AlignDir $OutDir
```

Total number of SVs:

```bash
for Vcf in $(ls analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.*.vcf | grep -v -e 'unfiltered' -e 'filtered' -e 'no_errors'); do
  echo "$Vcf"
cat $Vcf | grep -v '#' | grep -v 'FILTER' | wc -l
cat $Vcf | grep -v -e '#' -e 'FILTER' | cut -f7 | sort | uniq -c | sort -nr
done
```

```
analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.indel.vcf
74424
  74424 PASS
analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.sv.vcf
164
    164 PASS
```

Filter vcf files to remove low quality calls.

```bash
for Vcf in $(ls analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.*.vcf | grep -v -e 'unfiltered' -e 'filtered' -e 'no_errors'); do
  OutName=$(echo $Vcf | sed 's/.vcf/.filtered.vcf/g')
  echo "$OutName"
  cat $Vcf | grep -e '#' -e 'FILTER' -e 'PASS' > $OutName
  cat $OutName | grep -v -e '#' -e 'FILTER' | wc -l
done
```

```
analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.indel.filtered.vcf
74424
analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.sv.filtered.vcf
164
```

 ```bash
   Reference=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
   Gff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.gff3)
   SnpEff=/home/sobczm/bin/snpEff
   mkdir $SnpEff/data/P414v1.0
   cp $Reference $SnpEff/data/P414v1.0/sequences.fa
   cp $Gff $SnpEff/data/P414v1.0/genes.gff

   #Build database using GFF3 annotation
   java -jar $SnpEff/snpEff.jar build -gff3 -v P414v1.0
 ```


## Remove sequencing errors from vcf files:

```bash
for Vcf in $(ls analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.*.vcf | grep -v 'unfiltered' | grep 'filtered'); do
OutDir=$(dirname $Vcf)
Prefix=$(basename $Vcf .vcf)
Errors=$OutDir/${Prefix}_errors.tsv
FilteredVcf=$OutDir/${Prefix}_no_errors.vcf
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
$ProgDir/flag_error_SNPs.py --inp_vcf $Vcf --ref_isolate 414_vs_414_aligned_sorted.bam --errors $Errors --filtered $FilteredVcf
echo $Prefix
echo "The number of probable errors from homozygous variants being called from reference illumina reads vs the reference assembly is:"
cat $Errors | wc -l
echo "These have been removed from the vcf file"
done
```

```
Pcac_svaba_sv.svaba.indel.filtered
The number of probable errors from homozygous variants being called from reference illumina reads vs the reference assembly is:
276
These have been removed from the vcf file
Pcac_svaba_sv.svaba.sv.filtered
The number of probable errors from homozygous variants being called from reference illumina reads vs the reference assembly is:
1
These have been removed from the vcf file
```


## Collect VCF stats

General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
VcfTools=/home/sobczm/bin/vcftools/bin
export PERL5LIB="$VcfTools:$PERL5LIB"
for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_no_errors.vcf | grep -v 'unfiltered' | grep 'filtered' | grep 'no_errors'); do
OutDir=$(dirname $Vcf)
Prefix=$(basename $Vcf .vcf)
echo $Prefix
perl $VcfTools/vcf-stats $Vcf > $OutDir/$Prefix.stats
done
```

Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
  for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_no_errors.vcf | grep -v 'unfiltered' | grep 'filtered' | grep 'no_errors'); do
      ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
      $ProgDir/similarity_percentage.py $Vcf
  done
```

# Visualise the output as heatmap and clustering dendrogram
```bash
for Log in $(ls analysis/popgen/indel_calling/svaba/*distance.log | grep -v 'unfiltered' | grep 'filtered' | grep 'no_errors'); do
  OutDir=$(dirname $Log)
  Prefix=$(basename $Log .log)
  echo $Prefix
  ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
  Rscript --vanilla $ProgDir/distance_matrix.R $Log
  mv Rplots.pdf $OutDir/$Prefix.pdf
done
```


# Identify SVs in gene models:

A custom SNPeff database was created in 01_Pcac_alignments.md

## Annotate VCF files

http://snpeff.sourceforge.net/SnpEff_manual.html#run

```bash
CurDir=/data/scratch/armita/idris
cd $CurDir
for a in $(ls analysis/popgen/indel_calling/svaba/*_no_errors.vcf | grep -v 'unfiltered' | grep 'filtered' | grep 'no_errors'); do
    echo $a
    filename=$(basename "$a")
    Prefix=${filename%.vcf}
    OutDir=$(dirname $a)
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $a > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
    # mv 414_v2_contigs_unmasked_filtered* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has '3_prime_UTR_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has 'splice_region_variant') || (ANN[*].EFFECT has 'splice_donor_variant') || (ANN[*].EFFECT has 'splice_acceptor_variant') || (ANN[*].EFFECT has 'frameshift_variant') || (ANN[*].EFFECT has 'disruptive_inframe_deletion') || (ANN[*].EFFECT has 'inframe_deletion') || (ANN[*].EFFECT has 'disruptive_inframe_insertion') || (ANN[*].EFFECT has 'inframe_insertion') || (ANN[*].EFFECT has 'stop_lost') ||
    (ANN[*].EFFECT has 'stop_gained') || (ANN[*].EFFECT has 'exon_loss_variant') || (ANN[*].EFFECT has 'transcript_ablation')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf

    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'frameshift_variant') || (ANN[*].EFFECT has 'disruptive_inframe_deletion') || (ANN[*].EFFECT has 'inframe_deletion') || (ANN[*].EFFECT has 'disruptive_inframe_insertion') || (ANN[*].EFFECT has 'inframe_insertion') || (ANN[*].EFFECT has 'stop_lost') ||
    (ANN[*].EFFECT has 'stop_gained') || (ANN[*].EFFECT has 'exon_loss_variant') || (ANN[*].EFFECT has 'transcript_ablation')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'frameshift_variant') || (ANN[*].EFFECT has 'disruptive_inframe_deletion') || (ANN[*].EFFECT has 'inframe_deletion') || (ANN[*].EFFECT has 'disruptive_inframe_insertion') || (ANN[*].EFFECT has 'inframe_insertion') || (ANN[*].EFFECT has 'stop_lost') ||
    (ANN[*].EFFECT has 'stop_gained') || (ANN[*].EFFECT has 'exon_loss_variant') || (ANN[*].EFFECT has 'transcript_ablation')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
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
    # #-
    # # SNPs in BUSCO
    # #-
    # AnnotaTable=$(ls gene_pred/annotation/P.cactorum/414_v2/414_v2_gene_table_incl_exp.tsv)
    # AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi2.tsv)
    # Busco=$(ls gene_pred/busco/P.cactorum/414_v2/genes/run_final_genes_combined.gene/busco_single_copy_gene_headers.txt)
    # RxLR=$(ls gene_pred/annotation/P.cactorum/414_v2/renamed_RxLR.txt)
    # CRN=$(ls analysis/CRN_effectors/hmmer_CRN/P.cactorum/414_v2/414_v2_final_CRN_ID.txt)
    # cat $AnnotaTable | cut -f1,12 | tail -n+2 | grep 'Yes' | cut -f1 > $RxLR
    # cat $AnnotaTable | cut -f1,13 | tail -n+2 | grep 'Yes' | cut -f1 > $CRN
    # #-
    # # syn SNPs in effectors:
    # #-
    # SynVcf=$OutDir/"$Prefix"_syn.vcf
    # ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
    # BuscoOut=$OutDir/"$Prefix"_syn_Busco.vcf
    # $ProgDir/vcf_extract_genes.py --vcf $SynVcf --gene_list $Busco > $BuscoOut
    # BuscoSynSnps=$(cat $BuscoOut | grep -v '#' | wc -l)
    # RxlrOut=$OutDir/"$Prefix"_syn_RxLR.vcf
    # $ProgDir/vcf_extract_genes.py --vcf $SynVcf --gene_list $RxLR > $RxlrOut
    # RxlrSynSnps=$(cat $RxlrOut | grep -v '#' | wc -l)
    # CrnOut=$OutDir/"$Prefix"_CRN.vcf
    # $ProgDir/vcf_extract_genes.py --vcf $SynVcf --gene_list $CRN > $CrnOut
    # CrnSynSnps=$(cat $CrnOut | grep -v '#' | wc -l)  
    # #-
    # # non-syn SNPs in effectors:
    # #-
    # NonSynVcf=$OutDir/"$Prefix"_nonsyn.vcf
    # BuscoOut=$OutDir/"$Prefix"_nonsyn_Busco.vcf
    # ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
    # $ProgDir/vcf_extract_genes.py --vcf $NonSynVcf --gene_list $Busco > $BuscoOut
    # BuscoNonSynSnps=$(cat $BuscoOut | grep -v '#' | wc -l)
    # RxlrOut=$OutDir/"$Prefix"_nonsyn_RxLR.vcf
    # $ProgDir/vcf_extract_genes.py --vcf $NonSynVcf --gene_list $RxLR > $RxlrOut
    # RxlrNonSynSnps=$(cat $RxlrOut | grep -v '#' | wc -l)
    # CrnOut=$OutDir/"$Prefix"_nonsyn_CRN.vcf
    # $ProgDir/vcf_extract_genes.py --vcf $NonSynVcf --gene_list $CRN > $CrnOut
    # CrnNonSynSnps=$(cat $CrnOut | grep -v '#' | wc -l)
    # printf "Comparison\$AllSnps\tGeneSnps\tCdsSnps\tSynSnps\tNonsynSnps\tBuscoSynSnps\tBuscoNonSynSnps\tRxlrSynSnps\tRxlrNonSynSnps\tCrnSynSnps\tCrnNonSynSnps\n"
    # printf "$Prefix\t$AllSnps\t$GeneSnps\t$CdsSnps\t$SynSnps\t$NonsynSnps\t$BuscoSynSnps\t$BuscoNonSynSnps\t$RxlrSynSnps\t$RxlrNonSynSnps\t$CrnSynSnps\t$CrnNonSynSnps\n"
# done

cat $OutDir/"$Prefix"_annotated.vcf | grep -v '#' | cut -f8 | cut -f2 -d '|' | sed 's/&/\n/g' | sort | uniq -c | sort -nr

# #-
# # Make venn diagrams
# # -
# # These relate to the number of SNPs differing from the P414
# # reference in eahc group
# for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_nonsyn*.vcf | grep -v -e 'recode' -e '.vcf_'); do
# Prefix=$(echo $Vcf | sed 's/.vcf//g')
# echo $Prefix
# # Group1="12420 15_13 15_7 2003_3 4032 404 4040 414 415 416 62471"
# # Group2="PC13_15 P295 R36_14"
# # Group3="371 SCRP370 SCRP376"
# Group1="12420_sorted.bam 15_13_sorted.bam 15_7_sorted.bam 2003_3_sorted.bam 4032_sorted.bam 404_sorted.bam 4040_sorted.bam 414_sorted.bam 415_sorted.bam 416_sorted.bam 62471_sorted.bam"
# Group2="PC13_15_sorted.bam P295_sorted.bam R36_14_sorted.bam"
# Group3="371_sorted.bam SCRP370_sorted.bam SCRP376_sorted.bam"
# ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
# $ProgDir/vcf_2_venn.py --vcf $Vcf --g1_name Pc_Fxa --g1_isolates $Group1 --g2_name Pc_Mxd --g2_isolates $Group2 --g3_name Pi_Ri --g3_isolates $Group3 --prefix $Prefix
done
```

```
# InDel
58159 intergenic_region
 6282 frameshift_variant
 5322 intron_variant
 1844 disruptive_inframe_deletion
  928 inframe_deletion
  868 splice_region_variant
  854 disruptive_inframe_insertion
  799 inframe_insertion
  173 stop_lost
  142 start_lost
  121 stop_gained
   57 splice_acceptor_variant
   55 splice_donor_variant
    2 exon_loss_variant
    1 transcript_ablation

# SV
95 intergenic_region
22 inframe_insertion
22 disruptive_inframe_insertion
19 frameshift_variant
 6 intron_variant
 2 splice_region_variant
 1 stop_lost
```

genic variants:
* intron_variant
* 3_prime_UTR_variant
* 5_prime_UTR_variant
* splice_region_variant
* splice_donor_variant
* splice_acceptor_variant
* frameshift_variant
* disruptive_inframe_deletion
* inframe_deletion
* disruptive_inframe_insertion
* inframe_insertion
* stop_lost
* stop_gained
* exon_loss_variant
* transcript_ablation
coding
* frameshift_variant
* disruptive_inframe_deletion
* inframe_deletion
* disruptive_inframe_insertion
* inframe_insertion
* stop_lost
* stop_gained
* exon_loss_variant
* transcript_ablation
non-synonymous
* frameshift_variant
* disruptive_inframe_deletion
* inframe_deletion
* disruptive_inframe_insertion
* inframe_insertion
* stop_lost
* stop_gained
* exon_loss_variant
* transcript_ablation
synonymous
*


## Annotate VCF files

http://snpeff.sourceforge.net/SnpEff_manual.html#run

```bash
CurDir=/data/scratch/armita/idris
cd $CurDir
for a in $(ls analysis/popgen/indel_calling/svaba/*_no_errors.vcf | grep -v 'unfiltered' | grep 'filtered' | grep 'no_errors'); do
    echo $a
    filename=$(basename "$a")
    Prefix=${filename%.vcf}
    OutDir=$(dirname $a)
    SnpEff=/home/sobczm/bin/snpEff
    # java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $a > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
    # mv 414_v2_contigs_unmasked_filtered* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant' ) || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_region_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')
" $OutDir/${Prefix}_annotated.vcf > $OutDir/${Prefix}_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${Prefix}_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${Prefix}_annotated.vcf > $OutDir/${Prefix}_nonsyn.vcf
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

cat $OutDir/"$Prefix"_annotated.vcf | grep -v '#' | cut -f8 | cut -f2 -d '|' | sed 's/&/\n/g' | sort | uniq -c | sort -nr

done
```

```
analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.indel.filtered_no_errors.vcf
58159 intergenic_region
 6282 frameshift_variant
 5322 intron_variant
 1844 disruptive_inframe_deletion
  928 inframe_deletion
  868 splice_region_variant
  854 disruptive_inframe_insertion
  799 inframe_insertion
  173 stop_lost
  142 start_lost
  121 stop_gained
   57 splice_acceptor_variant
   55 splice_donor_variant
    2 exon_loss_variant
    1 transcript_ablation
analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.sv.filtered_no_errors.vcf
   95 intergenic_region
   22 inframe_insertion
   22 disruptive_inframe_insertion
   19 frameshift_variant
    6 intron_variant
    2 splice_region_variant
    1 stop_lost
```

Phase variants

```bash
for Vcf in $(ls analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.*.filtered_no_errors.vcf); do
  OutName=$(echo $Vcf | sed 's/.filtered_no_errors.vcf/_phased.tsv/g')
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen/popgenome_scripts
  $ProgDir/phase_Pc_vcf.py --inp_vcf $Vcf > $OutName
  echo $(basename $OutName)
  cat $OutName | cut -f1 | sort | uniq -c | sort -nr
done
```

For Non-Syn SNPs

```bash
for Vcf in $(ls analysis/popgen/indel_calling/svaba/Pcac_svaba_sv.svaba.*.filtered_no_errors_nonsyn.vcf | head -n1); do
  OutName=$(echo $Vcf | sed 's/.filtered_no_errors_nonsyn.vcf/_phased_nonsyn.tsv/g')
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen/popgenome_scripts
  $ProgDir/phase_Pc_vcf.py --inp_vcf $Vcf > $OutName
  OutTxt=$(echo $Vcf | sed 's/.filtered_no_errors_nonsyn.vcf/_phased_nonsyn.txt/g')
  echo $(basename ${OutTxt})
  cat $OutName | cut -f1,9 | cut -f1,7 -d'|' | sed "s/\t.*|/\t/g"  > $OutTxt
  cat $OutTxt | cut -f1 | sort | uniq -c | sort -nr
done
```

```
Pcac_svaba_sv.svaba.indel_phased_nonsyn.txt
   5335 species private fixed
   3351 P. idaei private unfixed
    638 apple private unfixed
    342 leather rot variant
    287 apple private fixed
    280 ancestral variation crown rot fixed
    126 crown rot private unfixed
    125 P. cactorum private unfixed
     96 crown rot private fixed
     80 P. cactorum private crown rot fixed
     53 ancestral variation differentially fixed
     47 ancestral variation unfixed
     16 ancestral variation apple fixed
      2 P. cactorum private apple fixed
Pcac_svaba_sv.svaba.sv_phased_nonsyn.txt
     20 leather rot variant
     17 P. idaei private unfixed
     11 species private fixed
      6 crown rot private unfixed
      6 apple private unfixed
      2 P. cactorum private unfixed
      1 ancestral variation crown rot fixed
```

```bash
Phased=$(ls analysis/popgen/SNP_calling/414_contigs_softmasked_repeatmasker_TPSI_appended_filtered_no_errors_nonsyn_phased.tsv)
cat $Phased | grep -w -f analysis/popgen/SNP_calling/RxLR_genes.txt | cut -f1 | sort | uniq -c | sort -nr
```

```
284 species private fixed
   7 crown rot private fixed
   6 P. idaei private unfixed
   5 apple private fixed
   4 leather rot variant
   1 triallelic SNP
   1 P. cactorum private crown rot fixed
   1 apple private unfixed
```

on local computer:
```bash
cat 414_annotation_ncbi3.tsv | grep -w 'RxLR' | grep 'crown rot private fixed' > strawberry_unique_RxLR.tsv
cat strawberry_unique_RxLR.tsv | wc -l
cat 414_annotation_ncbi3.tsv | grep -w 'RxLR' | grep 'apple private fixed' > apple_unique_RxLR.tsv
cat apple_unique_RxLR.tsv | wc -l
cat 414_annotation_ncbi3.tsv | grep -w 'CRN' | grep 'crown rot private fixed' > strawberry_unique_CRN.tsv
cat strawberry_unique_CRN.tsv | wc -l
cat 414_annotation_ncbi3.tsv | grep -w 'CRN' | grep 'apple private fixed' > apple_unique_CRN.tsv
cat apple_unique_CRN.tsv | wc -l

cat 414_annotation_ncbi3.tsv | grep -w 'RxLR' | grep 'crown rot private' > strawberry_unique_RxLR.tsv
cat strawberry_unique_RxLR.tsv | wc -l
cat 414_annotation_ncbi3.tsv | grep -w 'RxLR' | grep 'apple private' > apple_unique_RxLR.tsv
cat apple_unique_RxLR.tsv | wc -l
cat 414_annotation_ncbi3.tsv | grep -w 'CRN' | grep 'crown rot private' > strawberry_unique_CRN.tsv
cat strawberry_unique_CRN.tsv | wc -l
cat 414_annotation_ncbi3.tsv | grep -w 'CRN' | grep 'apple private' > apple_unique_CRN.tsv
cat apple_unique_CRN.tsv | wc -l
```

## Filter vcf file by group

Rename vcf entries

```bash
  for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_no_errors.vcf | grep -v 'unfiltered' | grep 'filtered' | grep 'no_errors'); do
  VcfNew="${Vcf%.vcf}_renamed.vcf"
  # echo $VcfNew
  cat $Vcf | sed 's/_vs_414_aligned_sorted.bam//g' > $VcfNew
  done
```




### P. idaei vs P414

```bash
Prefix=Pi_vs_P414
OutDir=analysis/popgen/indel_calling/svaba/$Prefix
mkdir -p $OutDir

for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_no_errors_renamed.vcf); do
ExcludeList="12420 15_13 15_7 2003_3 4032 404 415 416 62471 PC13_15 P295 R36_14 414 4040 11-40 17-21 P421"
VarType=$(basename $Vcf | cut -f3 -d '.')
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/${Prefix}_${VarType}.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/${Prefix}_${VarType}.vcf --non-ref-ac-any 1 --max-missing 0.95 --recode --out $OutDir/${Prefix}_${VarType}_filtered_no_errors
done

for Vcf in $(ls $OutDir/"$Prefix"_*_filtered_no_errors.recode.vcf); do
    echo $Vcf
    # ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/summary_stats
    # $ProgDir/annotate_snps_genome.sh $Vcf P414v1.0

    filename=$(basename "$Vcf")
    OutPrefix="${filename%.vcf}"
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $Vcf > $OutDir/"$OutPrefix"_annotated.vcf

    #Create subsamples of SNPs containing those in a given category

    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant' ) || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_region_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')
" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${filename%.vcf}_annotated.vcf > $OutDir/"$OutPrefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/${OutPrefix}_syn.vcf
done
```


### P. cactorum ex. apple vs P414

```bash
Prefix=Pc_apple_vs_P414
OutDir=analysis/popgen/indel_calling/svaba/$Prefix
mkdir -p $OutDir

for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_no_errors_renamed.vcf); do
ExcludeList="12420 15_13 15_7 2003_3 4032 4040 404 415 416 371 SCRP370 SCRP376 414 4040 11-40 17-21 P421"
VarType=$(basename $Vcf | cut -f3 -d '.')
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/${Prefix}_${VarType}.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/${Prefix}_${VarType}.vcf --non-ref-ac-any 1 --max-missing 0.95 --recode --out $OutDir/${Prefix}_${VarType}_filtered_no_errors
done

for Vcf in $(ls $OutDir/"$Prefix"_*_filtered_no_errors.recode.vcf); do
    echo $Vcf
    filename=$(basename "$Vcf")
    OutPrefix="${filename%.vcf}"
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $Vcf > $OutDir/"$OutPrefix"_annotated.vcf

    #Create subsamples of SNPs containing those in a given category

    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant' ) || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_region_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')
" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${filename%.vcf}_annotated.vcf > $OutDir/"$OutPrefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$OutPrefix"_annotated.vcf > $OutDir/"$OutPrefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$OutPrefix"_syn.vcf
done
```


### P. cactorum ex. strawberry vs P414

```bash
Prefix=Pc_strawberry_vs_P414
OutDir=analysis/popgen/indel_calling/svaba/$Prefix
mkdir -p $OutDir

for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_no_errors_renamed.vcf); do
ExcludeList="62471 P295 R36_14 371 SCRP370 SCRP376 11-40 17-21"
VarType=$(basename $Vcf | cut -f3 -d '.')
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/${Prefix}_${VarType}.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/${Prefix}_${VarType}.vcf --non-ref-ac-any 1 --max-missing 0.95 --recode --out $OutDir/${Prefix}_${VarType}_filtered_no_errors
done

for Vcf in $(ls $OutDir/"$Prefix"_*_filtered_no_errors.recode.vcf); do
    echo $Vcf
    filename=$(basename "$Vcf")
    OutPrefix="${filename%.vcf}"
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $Vcf > $OutDir/"$OutPrefix"_annotated.vcf

    #Create subsamples of SNPs containing those in a given category

    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant' ) || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_region_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')
" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${filename%.vcf}_annotated.vcf > $OutDir/"$OutPrefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$OutPrefix"_annotated.vcf > $OutDir/"$OutPrefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$OutPrefix"_syn.vcf
done
```



### P414 vs P414

```bash
Prefix=P414_vs_P414
OutDir=analysis/popgen/indel_calling/svaba/$Prefix
mkdir -p $OutDir

for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_no_errors_renamed.vcf); do
ExcludeList="12420 15_13 15_7 2003_3 4032 404 415 416 PC13_15 62471 P295 R36_14 371 SCRP370 SCRP376 4040 11-40 17-21 P421"
VarType=$(basename $Vcf | cut -f3 -d '.')
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/${Prefix}_${VarType}.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/${Prefix}_${VarType}.vcf --non-ref-ac-any 1 --max-missing 0.95 --recode --out $OutDir/${Prefix}_${VarType}_filtered_no_errors
done

for Vcf in $(ls $OutDir/"$Prefix"_*_filtered_no_errors.recode.vcf); do
    echo $Vcf
    filename=$(basename "$Vcf")
    OutPrefix="${filename%.vcf}"
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $Vcf > $OutDir/"$OutPrefix"_annotated.vcf

    #Create subsamples of SNPs containing those in a given category

    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant' ) || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_region_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')
" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${filename%.vcf}_annotated.vcf > $OutDir/"$OutPrefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$OutPrefix"_annotated.vcf > $OutDir/"$OutPrefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$OutPrefix"_syn.vcf
done
```

### Leather rot vs P414

```bash
Prefix=Pc_leather_rot_vs_P414
OutDir=analysis/popgen/indel_calling/svaba/$Prefix
mkdir -p $OutDir

for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_no_errors_renamed.vcf); do
ExcludeList="12420 15_13 15_7 2003_3 4032 404 415 416 PC13_15 62471 P295 R36_14 414 371 SCRP370 SCRP376 4040 P421"
VarType=$(basename $Vcf | cut -f3 -d '.')
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/${Prefix}_${VarType}.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/${Prefix}_${VarType}.vcf --non-ref-ac-any 1 --max-missing 0.95 --recode --out $OutDir/${Prefix}_${VarType}_filtered_no_errors
done

for Vcf in $(ls $OutDir/"$Prefix"_*_filtered_no_errors.recode.vcf); do
    echo $Vcf
    filename=$(basename "$Vcf")
    OutPrefix="${filename%.vcf}"
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $Vcf > $OutDir/"$OutPrefix"_annotated.vcf

    #Create subsamples of SNPs containing those in a given category

    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant' ) || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_region_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')
" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${filename%.vcf}_annotated.vcf > $OutDir/"$OutPrefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$OutPrefix"_annotated.vcf > $OutDir/"$OutPrefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$OutPrefix"_syn.vcf
done
```

### Leather rot 17-21 vs P414

```bash
Prefix=Pc_17-21_vs_P414
OutDir=analysis/popgen/indel_calling/svaba/$Prefix
mkdir -p $OutDir

for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_no_errors_renamed.vcf); do
ExcludeList="12420 15_13 15_7 2003_3 4032 404 415 416 PC13_15 62471 P295 R36_14 414 371 SCRP370 SCRP376 4040 P421 11-40"
VarType=$(basename $Vcf | cut -f3 -d '.')
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/${Prefix}_${VarType}.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/${Prefix}_${VarType}.vcf --non-ref-ac-any 1 --max-missing 0.95 --recode --out $OutDir/${Prefix}_${VarType}_filtered_no_errors
done

for Vcf in $(ls $OutDir/"$Prefix"_*_filtered_no_errors.recode.vcf); do
    echo $Vcf
    filename=$(basename "$Vcf")
    OutPrefix="${filename%.vcf}"
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $Vcf > $OutDir/"$OutPrefix"_annotated.vcf

    #Create subsamples of SNPs containing those in a given category

    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant' ) || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_region_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')
" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${filename%.vcf}_annotated.vcf > $OutDir/"$OutPrefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$OutPrefix"_annotated.vcf > $OutDir/"$OutPrefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$OutPrefix"_syn.vcf
done
```

### Leather rot 11-40 vs P414

```bash
Prefix=Pc_11-40_vs_P414
OutDir=analysis/popgen/indel_calling/svaba/$Prefix
mkdir -p $OutDir

for Vcf in $(ls analysis/popgen/indel_calling/svaba/*_no_errors_renamed.vcf); do
ExcludeList="12420 15_13 15_7 2003_3 4032 404 415 416 PC13_15 62471 P295 R36_14 414 371 SCRP370 SCRP376 4040 P421 17-21"
VarType=$(basename $Vcf | cut -f3 -d '.')
VcfLib=/home/sobczm/bin/vcflib/bin
$VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/${Prefix}_${VarType}.vcf
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $OutDir/${Prefix}_${VarType}.vcf --non-ref-ac-any 1 --max-missing 0.95 --recode --out $OutDir/${Prefix}_${VarType}_filtered_no_errors
done

for Vcf in $(ls $OutDir/"$Prefix"_*_filtered_no_errors.recode.vcf); do
    echo $Vcf
    filename=$(basename "$Vcf")
    OutPrefix="${filename%.vcf}"
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $Vcf > $OutDir/"$OutPrefix"_annotated.vcf

    #Create subsamples of SNPs containing those in a given category

    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant' ) || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_region_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')
" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${filename%.vcf}_annotated.vcf > $OutDir/"$OutPrefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'start_lost') || (ANN[0].EFFECT has 'stop_gained') || (ANN[0].EFFECT has 'stop_lost') || (ANN[0].EFFECT has 'inframe_deletion') || (ANN[0].EFFECT has 'inframe_insertion') || (ANN[0].EFFECT has 'transcript_ablation')" $OutDir/${OutPrefix}_annotated.vcf > $OutDir/${OutPrefix}_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$OutPrefix"_annotated.vcf > $OutDir/"$OutPrefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$OutPrefix"_syn.vcf
done
```


## Summarise SNP effects

```bash
cd /data/scratch/armita/idris/
AnnotaTable=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi2.tsv)
Busco=$(ls gene_pred/busco/P.cactorum/414/genes/run_final_genes_genes_incl_ORFeffectors_renamed.gene/busco_single_copy_gene_headers.txt)

OutDir=analysis/popgen/SNP_calling
mkdir $OutDir
# RxLR=$(ls gene_pred/annotation/P.cactorum/414/renamed_RxLR.txt)
# CRN=$(ls analysis/CRN_effectors/hmmer_CRN/P.cactorum/414_v2/414_v2_final_CRN_ID.txt)
cat $AnnotaTable | grep -w 'RxLR' | cut -f1 > $OutDir/RxLR_genes.txt
cat $AnnotaTable | grep -w 'CRN' | cut -f1 > $OutDir/CRN_genes.txt
  printf "Comparison\tvariant\tAllSnps\tGeneSnps\tCdsSnps\tSynSnps\tNonsynSnps\tBuscoSynSnps\tBuscoNonSynSnps\tRxlrSynSnps\tRxlrNonSynSnps\tCrnSynSnps\tCrnNonSynSnps\n"
for Type in "indel" "sv"; do
  for Folder in $(ls -d analysis/popgen/indel_calling/svaba/*_vs_P414*); do
  # for Folder in $(ls -d analysis/popgen/indel_calling/svaba/*_vs_P414*  | grep -e '17-21' -e '11-40'); do

    Comparison=$(echo $Folder | rev | cut -f1 -d '/' | rev)
    AllSnps=$(cat $Folder/*_${Type}_filtered_no_errors.recode_annotated.vcf | grep -v '#' | wc -l)
    GeneSnps=$(cat $Folder/*_${Type}_filtered_no_errors.recode_gene.vcf | grep -v '#' | wc -l)
    CdsSnps=$(cat $Folder/*_${Type}_filtered_no_errors.recode_coding.vcf | grep -v '#' | wc -l)
    NonsynSnps=$(cat $Folder/*_${Type}_filtered_no_errors.recode_nonsyn.vcf | grep -v '#' | wc -l)
    SynSnps=$(cat $Folder/*_${Type}_filtered_no_errors.recode_syn.vcf | grep -v '#' | wc -l)
    #syn SNPs in effectors:
    BuscoOut=$Folder/${Comparison}_${Type}_filtered_no_errors_syn_Busco.vcf
    cat $Folder/*_${Type}_filtered_no_errors.recode_syn.vcf | grep -w -f $Busco > $BuscoOut
    BuscoSynSnps=$(cat $BuscoOut | wc -l)
    RxlrOut=$Folder/${Comparison}_${Type}_filtered_no_errors_syn_RxLR.vcf
    cat $Folder/*_${Type}_filtered_no_errors.recode_syn.vcf | grep -f $OutDir/RxLR_genes.txt > $RxlrOut
    RxlrSynSnps=$(cat $RxlrOut | wc -l)
    CrnOut=$Folder/${Comparison}_${Type}_filtered_no_errors_syn_CRN.vcf
    cat $Folder/*_${Type}_filtered_no_errors.recode_syn.vcf | grep -f $OutDir/CRN_genes.txt > $CrnOut
    CrnSynSnps=$(cat $CrnOut | wc -l)  
    # non-syn SNPs in effectors:
    BuscoOut=$Folder/${Comparison}_${Type}_filtered_no_errors_nonsyn_Busco.vcf
    cat $Folder/*_${Type}_filtered_no_errors.recode_nonsyn.vcf | grep -w -f $Busco > $BuscoOut
    BuscoNonSynSnps=$(cat $BuscoOut | wc -l)
    RxlrOut=$Folder/${Comparison}_${Type}_filtered_no_errors_nonsyn_RxLR.vcf
    cat $Folder/*_${Type}_filtered_no_errors.recode_nonsyn.vcf | grep -f $OutDir/RxLR_genes.txt > $RxlrOut
    RxlrNonSynSnps=$(cat $RxlrOut | wc -l)
    CrnOut=$Folder/${Comparison}_${Type}_filtered_no_errors_nonsyn_CRN.vcf
    cat $Folder/*_${Type}_filtered_no_errors.recode_nonsyn.vcf | grep -f $OutDir/CRN_genes.txt > $CrnOut
    CrnNonSynSnps=$(cat $CrnOut | wc -l)
    printf "$Comparison\t${Type}\t$AllSnps\t$GeneSnps\t$CdsSnps\t$SynSnps\t$NonsynSnps\t$BuscoSynSnps\t$BuscoNonSynSnps\t$RxlrSynSnps\t$RxlrNonSynSnps\t$CrnSynSnps\t$CrnNonSynSnps\n"
  done
  # Comparison=$(echo $Folder | rev | cut -f1 -d '/' | rev)
  # AllSnps=$(cat $Folder/*_filtered_no_errors.recode_annotated.vcf | grep -v '#' | wc -l)
  # GeneSnps=$(cat $Folder/*_filtered_no_errors.recode_gene.vcf | grep -v '#' | wc -l)
  # CdsSnps=$(cat $Folder/*_filtered_no_errors.recode_coding.vcf | grep -v '#' | wc -l)
  # NonsynSnps=$(cat $Folder/*_filtered_no_errors.recode_nonsyn.vcf | grep -v '#' | wc -l)
  # SynSnps=$(cat $Folder/*_filtered_no_errors.recode_syn.vcf | grep -v '#' | wc -l)
  # #syn SNPs in effectors:
  # BuscoOut=$Folder/${Comparison}_filtered_no_errors_syn_Busco.vcf
  # cat $Folder/*_${Type}_filtered_no_errors.recode_syn.vcf | grep -w -f $Busco > $BuscoOut
  # BuscoSynSnps=$(cat $BuscoOut | wc -l)
  # RxlrOut=$Folder/${Comparison}_filtered_no_errors_syn_RxLR.vcf
  # cat $Folder/*_${Type}_filtered_no_errors.recode_syn.vcf | grep -f $OutDir/RxLR_genes.txt > $RxlrOut
  # RxlrSynSnps=$(cat $RxlrOut | wc -l)
  # CrnOut=$Folder/${Comparison}_filtered_no_errors_syn_CRN.vcf
  # cat $Folder/*_${Type}_filtered_no_errors.recode_syn.vcf | grep -f $OutDir/CRN_genes.txt > $CrnOut
  # CrnSynSnps=$(cat $CrnOut | wc -l)  
  # # non-syn SNPs in effectors:
  # BuscoOut=$Folder/${Comparison}_filtered_no_errors_nonsyn_Busco.vcf
  # cat $Folder/*_${Type}_filtered_no_errors.recode_nonsyn.vcf | grep -w -f $Busco > $BuscoOut
  # BuscoNonSynSnps=$(cat $BuscoOut | wc -l)
  # RxlrOut=$Folder/${Comparison}_filtered_no_errors_nonsyn_RxLR.vcf
  # cat $Folder/*_${Type}_filtered_no_errors.recode_nonsyn.vcf | grep -f $OutDir/RxLR_genes.txt > $RxlrOut
  # RxlrNonSynSnps=$(cat $RxlrOut | wc -l)
  # CrnOut=$Folder/${Comparison}_filtered_no_errors_nonsyn_CRN.vcf
  # cat $Folder/*_filtered_no_errors.recode_nonsyn.vcf | grep -f $OutDir/CRN_genes.txt > $CrnOut
  # CrnNonSynSnps=$(cat $CrnOut | wc -l)
  # printf "$Comparison\t${Type}\t$AllSnps\t$GeneSnps\t$CdsSnps\t$SynSnps\t$NonsynSnps\t$BuscoSynSnps\t$BuscoNonSynSnps\t$RxlrSynSnps\t$RxlrNonSynSnps\t$CrnSynSnps\t$CrnNonSynSnps\n"
done
```

```
Comparison      Variant All Gene        Cds Syn Nonsyn      BuscoSyn    BuscoNonSyn RxlrSyn RxlrNonSyn  CrnSyn      CrnNonSyn
P414_vs_P414    indel   712     162     131     0       131     0       0       0       0       0       2
Pc_11-40_vs_P414        indel   1057    241     191     0       191     0       0       0       1       0       2
Pc_17-21_vs_P414        indel   9200    1791    1299    0       1299    0       6       0       14      0       6
Pc_apple_vs_P414        indel   11383   2249    1659    0       1659    0       6       0       12      0       12
Pc_leather_rot_vs_P414  indel   9647    1883    1370    0       1370    0       6       0       14      0       7
Pc_strawberry_vs_P414   indel   1760    394     316     0       316     0       1       0       2       0       5
Pi_vs_P414      indel   63362   14005   9282    0       9282    0       64      0       81      0       32
P414_vs_P414    sv      0       0       0       0       0       0       0       0       0       0       0
Pc_11-40_vs_P414        sv      0       0       0       0       0       0       0       0       0       0       0
Pc_17-21_vs_P414        sv      16      6       6       0       6       0       0       0       0       0       0
Pc_apple_vs_P414        sv      16      9       9       0       9       0       0       0       0       0       0
Pc_leather_rot_vs_P414  sv      16      6       6       0       6       0       0       0       0       0       0
Pc_strawberry_vs_P414   sv      12      8       8       0       8       0       0       0       0       0       0
Pi_vs_P414      sv      92      33      29      0       29      0       0       0       0       0       0
```


# 3.0 Comparisons of groups to reference P414 genome
<!--
```bash
AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi2.tsv)

  for Vcf in $(ls analysis/popgen/indel_calling/svaba/*.filtered_no_errors_annotated.vcf | head -n1); do
    summarise variants
  done

```
-->

<!--
# 3.1 P. idaei vs P414

```bash
  Prefix=Pi_vs_P414
  OutDir=analysis/popgen/indel_calling/svaba/$Prefix
  mkdir -p $OutDir

  for Vcf in $(ls analysis/popgen/indel_calling/svaba/*.filtered_no_errors_annotated.vcf | head -n1); do
    Variant=$(echo $Vcf | rev |  cut -f3 -d '.' | rev)
    echo "$Variant"
    x="_vs_414_aligned_sorted.bam"
    ExcludeList="12420$x 15_13$x 15_7$x 2003_3$x 4032$x 404$x 415$x 416$x 62471$x PC13_15$x P295$x R36_14$x 414$x 4040$x 11-40$x 17-21$x P421$x"
    VcfLib=/home/sobczm/bin/vcflib/bin
    $VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/${Prefix}_${Variant}.vcf
    VcfTools=/home/sobczm/bin/vcftools/bin
    $VcfTools/vcftools --vcf $OutDir/${Prefix}_${Variant}.vcf --max-non-ref-ac 1 --out $OutDir/${Prefix}_${Variant}_filtered
    # $VcfTools/vcftools --vcf $OutDir/${Prefix}_${Variant}.vcf --maf 0.001 --recode --out $OutDir/${Prefix}_${Variant}_filtered
  done



  for Vcf in $(ls $OutDir/"$Prefix"_filtered_no_indels.recode.vcf); do
      echo $Vcf
      ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/summary_stats
      $ProgDir/annotate_snps_genome.sh $Vcf P414v1.0

      filename=$(basename "$Vcf")
      Prefix=$(echo $filename | sed 's/.vcf//g')
      SnpEff=/home/sobczm/bin/snpEff
      java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 P414v1.0 $Vcf > $OutDir/"$Prefix"_annotated.vcf
      mv snpEff_genes.txt $OutDir/snpEff_genes_$Prefix.txt
      mv snpEff_summary.html $OutDir/snpEff_summary_$Prefix.html

      #Create subsamples of SNPs containing those in a given category

      #genic (includes 5', 3' UTRs)
      java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
      #coding
      java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/${filename%.vcf}_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
      #non-synonymous
      java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
      #synonymous
      java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
      #Four-fold degenrate sites (output file suffix: 4fd)
      ProgDir=/home/sobczm/bin/popgen/summary_stats
      python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf
  done
``` -->
