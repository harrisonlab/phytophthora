Commands for generating a circos plot for P. cac isolate 414


```bash
  OutDir=analysis/circos/P.cactorum/414
  mkdir -p $OutDir
  # Convert the Fus2 genome into circos format
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Genome=assembly/merged_canu_spades/P.cactorum/414/filtered_contigs/contigs_min_500bp_renamed.fasta
  $ProgDir/fasta2circos.py --genome $Genome --contig_prefix "" > $OutDir/414_genome.txt

  # Make 100kb windows for plots
  $ProgDir/fasta2gff_windows.py --genome $Genome  > $OutDir/414_100kb_windows.gff
  # Identify GC content in 100kb windows
  $ProgDir/gc_content2circos.py --genome $Genome --gff $OutDir/414_100kb_windows.gff > $OutDir/414_GC_scatterplot.txt

  # Identify gene density in 100Kb windows
  GeneGff=gene_pred/codingquary/P.cactorum/414/final/final_genes_appended.gff3
  bedtools coverage -a $GeneGff -b $OutDir/414_100kb_windows.gff > $OutDir/414_gene_density.bed
  # Convert coverage bed files into circos format
  $ProgDir/coverage_bed2circos.py --bed $OutDir/414_gene_density.bed > $OutDir/414_gene_density_lineplot.txt

  # Plot location of RxLR genes as a scatterplot
  RxLR_gff=analysis/RxLR_effectors/combined_evidence/P.cactorum/414/414_total_RxLR.gff
  $ProgDir/gff2circos_scatterplot.py --gff $RxLR_gff --feature gene --value 0.66 > $OutDir/414_RxLR_scatterplot.txt
  # Plot location of CRN genes as a scatterplot
  CRN_gff=analysis/CRN_effectors/hmmer_CRN/P.cactorum/414/414_pub_CRN_LFLAK_DWL.gff
  $ProgDir/gff2circos_scatterplot.py --gff $CRN_gff --feature gene --value 0.33 > $OutDir/414_CRN_scatterplot.txt

  # Convert FoC MiSeq reads aligning in 100kb windows into coverage stats
  for ReadsBam in $(ls analysis/genome_alignment/bowtie/P.*/*/vs_414/414_contigs_unmasked.fa_aligned_sorted.bam); do
    Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
    # AlignDir=$(dirname $ReadsBam)
    echo "$Organism - $Strain"
    bedtools coverage -abam $ReadsBam -b $OutDir/414_100kb_windows.gff > $OutDir/"$Strain"_coverage_vs_414.bed
    # Convert coverage bed files into circos format
    $ProgDir/coverage_bed2circos.py --bed $OutDir/"$Strain"_coverage_vs_414.bed > $OutDir/"$Strain"_coverage_vs_414_scatterplot.txt
  done

  circos -conf /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/circos/414_circos.conf -outputdir $OutDir
```
