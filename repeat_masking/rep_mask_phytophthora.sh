#!/usr/bin/bash

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/P.cactorum/404/404_assembly.51/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/P.cactorum/414/414_assembly.61/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/P.ideai/371/371_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/kamoun_lab/P.cactorum/10300/version1/Pcactorum_contigs_singletons_14436.fa 
#qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/genbank/P.fragariae/JHVZ02/assembly_version2/JHVZ02.1.fsa_nt
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/genbank/P.fragariae/JHVZ02/assembly_version2_ed/JHVZ02_renamed.fa
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/P.fragariae/SCRP245/SCRP245_assembly.61/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/rep_modeling.sh assembly/velvet/P.cactorum/411/411_assembly.41/sorted_contigs.fa

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/P.cactorum/404/404_assembly.51/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/P.cactorum/414/414_assembly.61/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/velvet/P.ideai/371/371_assembly.41/sorted_contigs.fa

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking/transposonPSI.sh assembly/genbank/P.fragariae/JHVZ02/assembly_version2_ed/JHVZ02_renamed.fa