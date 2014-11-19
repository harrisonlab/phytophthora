#!/bin/bash
# Phytophthora submit path pipe.

for infile in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/path_pipe.sh $infile
done
