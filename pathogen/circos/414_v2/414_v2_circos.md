Commands for generating a circos plot for P. cac isolate 414

```bash
  OutDir=analysis/circos/P.cactorum/414_final
  mkdir -p $OutDir
  # Convert the Fus2 genome into circos format
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Genome=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
  $ProgDir/fasta2circos.py --genome $Genome --contig_prefix "" > $OutDir/414_genome.txt

  # Make 100kb windows for plots
  $ProgDir/fasta2gff_windows.py --genome $Genome  > $OutDir/414_100kb_windows.gff
  # Identify GC content in 100kb windows
  $ProgDir/gc_content2circos.py --genome $Genome --gff $OutDir/414_100kb_windows.gff > $OutDir/414_GC_scatterplot.txt

  # Identify gene density in 100Kb windows
  GeneGff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/414_genes_incl_ORFeffectors.gff3)
  bedtools coverage -a $GeneGff -b $OutDir/414_100kb_windows.gff > $OutDir/414_gene_density.bed
  # Convert coverage bed files into circos format
  $ProgDir/coverage_bed2circos.py --bed $OutDir/414_gene_density.bed > $OutDir/414_gene_density_lineplot.txt

  # # Plot location of transposons as a scatterplot
  # Repeatmasker_gff=repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_transposonmasked.gff
  # TransposonPSI_gff=repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa.TPSI.allHits.chains.gff3
  # $ProgDir/gff2circos_scatterplot.py --gff $Repeatmasker_gff --feature similarity --value 0.5 > $OutDir/414_transposon_scatterplot.txt
  # $ProgDir/gff2circos_scatterplot.py --gff $TransposonPSI_gff --feature translated_nucleotide_match --value 0.5 >> $OutDir/414_transposon_scatterplot.txt

  # Plot location of transposons as a line plot
  Repeatmasker_gff=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_transposonmasked.gff)
  TransposonPSI_gff=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa.TPSI.allHits.chains.gff3)

  bedtools coverage -a $Repeatmasker_gff -b $OutDir/414_100kb_windows.gff > $OutDir/414_transposon_density.bed
    # bedtools coverage -a $TransposonPSI_gff -b $OutDir/414_100kb_windows.gff > $OutDir/414_transposon_density.bed
  # Convert coverage bed files into circos format
  $ProgDir/coverage_bed2circos.py --bed $OutDir/414_transposon_density.bed --per_X_bp "100000" > $OutDir/414_transposon_density_lineplot.txt


  # Plot location of RxLR genes as a scatterplot
  RxLR_gff=$(ls analysis/RxLR_effectors/combined_evidence/P.cactorum/414/414_total_RxLR.gff)
  $ProgDir/gff2circos_scatterplot.py --gff $RxLR_gff --feature gene --value 0.66 > $OutDir/414_RxLR_scatterplot.txt
  # Plot location of CRN genes as a scatterplot
  CRN_gff=$(ls analysis/CRN_effectors/hmmer_CRN/P.cactorum/414/414_pub_CRN_LFLAK_DWL.gff)
  $ProgDir/gff2circos_scatterplot.py --gff $CRN_gff --feature gene --value 0.33 > $OutDir/414_CRN_scatterplot.txt

 # Note - this next step needs to be ran in a qlogin session with high RAM

 # Location of DEG effectors:
 GeneGff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/414_genes_incl_ORFeffectors.gff3)
 AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
 cat $AnnotTab | grep -w 'DEG' | grep -w 'RxLR' | cut -f1 > $OutDir/DEG_RxLRs_headers.txt
 cat $GeneGff  | grep -f $OutDir/DEG_RxLRs_headers.txt > $OutDir/DEG_RxLRs.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $OutDir/DEG_RxLRs.gff --feature mRNA --value 0.66 > $OutDir/DEG_RxLRs_scatterplot.txt

cat $AnnotTab | grep -w 'DEG' | grep -w 'CRN' | cut -f1 > $OutDir/DEG_CRNs_headers.txt
cat $GeneGff  | grep -f $OutDir/DEG_CRNs_headers.txt > $OutDir/DEG_CRNs.gff
$ProgDir/gff2circos_scatterplot.py --gff $OutDir/DEG_CRNs.gff --feature mRNA --value 0.66 > $OutDir/DEG_CRNs_scatterplot.txt


#  Convert FoC MiSeq reads aligning in 100kb windows into coverage stats
OutDir=analysis/circos/P.cactorum/414_final
for ReadsBam in $(ls assembly/merged_SMARTdenovo_spades/P.cactorum/414/polished/aligned_MiSeq/contigs_min_500bp_renamed.fasta_aligned.bam); do
Strain=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Organism=$(echo $ReadsBam | rev | cut -f5 -d '/' | rev)
# AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -a $OutDir/414_100kb_windows.gff -b $ReadsBam > $OutDir/"$Strain"_coverage_vs_414.bed
# Convert coverage bed files into circos format
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/coverage_bed2circos.py --bed $OutDir/"$Strain"_coverage_vs_414.bed > $OutDir/"$Strain"_coverage_vs_414_scatterplot.txt
done

for Vcf in $(ls analysis/popgen/SNP_calling/*/*_filtered_no_indels.recode.vcf); do
  echo $Vcf
  Filename=$(basename $Vcf)
  OutFile=$(echo $Filename | sed 's/.vcf/_scatterplot.txt/g')
  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
  $ProgDir/vcf_2circosplot.py  --inp_vcf $Vcf > $OutDir/$OutFile
done

  circos -conf /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/circos/414_v2/414_v2_circos.conf -outputdir $OutDir
```

Note - Bedtools overestimates coverage from paired end data, treating an overlapping paired read as 2X coverage.
You should therefore be careful to use coverage plots as a representation of areas that show coverage rather
than a quantification of coverage in an area. Discussed at: https://www.biostars.org/p/172179/
