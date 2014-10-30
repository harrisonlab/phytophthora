#!/usr/bin/bash

for infile in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa); do
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/oomycete_avr_genes/appended_oomycete_avr_cds.fasta  dna $infile
done

for infile in $(ls analysis/blast_homology/*/*/*_appended_oomycete_avr_cds.fasta_homologs.csv); do
OUTFILE=$(echo $infile | sed 's/fasta.*/gff/')
/home/armita/git_repos/emr_repos/tools/pathogen/blast/blast2gff.pl oomycete_avr_gene $infile > $OUTFILE
done

/home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_differentials.sh analysis/blast_homology/P.cactorum/404/404_appended_oomycete_avr_cds.fasta_homologs.csv 