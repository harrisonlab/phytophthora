#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 6
#$ -l h_vmem=4G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen/popgenome_scripts
Rscript --vanilla $ProgDir/calculate_nucleotide_diversity.R
