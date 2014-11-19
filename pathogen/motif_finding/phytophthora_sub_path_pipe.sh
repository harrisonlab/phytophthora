SCRIPT_DIR=/home/armita/git_repos/emr_repos/tools/pathogen/path_pipe.sh

qsub $SCRIPT_DIR assembly/velvet/P.cactorum/404/404_assembly.41/sorted_contigs.fa 
qsub $SCRIPT_DIR assembly/velvet/P.cactorum/414/414_assembly.41/sorted_contigs.fa  
qsub $SCRIPT_DIR assembly/kamoun_lab/P.cactorum/10300/version1/Pcactorum_contigs_singletons_14436.fa
qsub $SCRIPT_DIR assembly/velvet/P.ideai/371/371_assembly.41/sorted_contigs.fa 
qsub $SCRIPT_DIR assembly/genbank/P.fragariae/JHVZ02/assembly_version2/JHVZ02.1.fsa_nt