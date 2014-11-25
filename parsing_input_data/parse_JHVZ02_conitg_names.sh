#!/bin/bash

# commands used to parse JHVZ02 fasta headers

# commands used to previously rename masked contigs
#cat repeat_masked/P.fragariae/JHVZ02/assembly_version2_repmask/JHVZ02_contigs_hardmasked.fa | sed 's/.*|.*|.*|.*| />/' | sed 's/,.*//' | sed 's/ /_/g' |  sed 's/Phytophthora_/P./' | sed 's/strain_//' | sed 's/3runs_c/contig_/' > repeat_masked/P.fragariae/JHVZ02/assembly_version2_repmask/JHVZ02_contigs_hardmasked.fa

#commands used to rename downloaded contigs before masking.
cat assembly/genbank/P.fragariae/JHVZ02/assembly_version2/JHVZ02.1.fsa_nt | sed 's/gi|.*|gb|/P.frag_309.62_NODE_/' | sed 's/| Phytophthora fragariae strain 309.62 3runs//' | sed 's/, whole genome shotgun sequence//' > assembly/genbank/P.fragariae/JHVZ02/assembly_version2_ed/JHVZ02_renamed.fa