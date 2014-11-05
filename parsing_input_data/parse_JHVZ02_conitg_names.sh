#!/bin/bash

# commands used to parse JHVZ02 fasta headers

cat repeat_masked/P.fragariae/JHVZ02/assembly_version2_repmask/JHVZ02_contigs_hardmasked.fa | sed 's/.*|.*|.*|.*| />/' | sed 's/,.*//' | sed 's/ /_/g' |  sed 's/Phytophthora_/P./' | sed 's/strain_//' | sed 's/3runs_c/contig_/' > repeat_masked/P.fragariae/JHVZ02/assembly_version2_repmask/JHVZ02_contigs_hardmasked.fa