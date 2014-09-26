#!/usr/bin/bash

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/P.cactorum/404/404_assembly.51/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/P.cactorum/414/414_assembly.61/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/P.ideai/371/371_assembly.41/sorted_contigs.fa


qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/P.cactorum/404/404_assembly.51/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/P.cactorum/414/414_assembly.61/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/P.ideai/371/371_assembly.41/sorted_contigs.fa
