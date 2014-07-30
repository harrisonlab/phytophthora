#!/usr/bin/bash

QUERY=analysis/blast_homology/P.inf_AVR2.fa

PROGRAM=/home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh 

qsub $PROGRAM $QUERY assembly/velvet/P.cactorum/404/404_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/P.cactorum/414/414_assembly.41/sorted_contigs.fa 

qsub $PROGRAM $QUERY assembly/velvet/P.idaei/371/371_assembly.41/sorted_contigs.fa 

qsub $PROGRAM $QUERY assembly/genbank/P.fragariae/JHVZ02/assembly_version2/JHVZ02.1.fsa_nt