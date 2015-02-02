#!/bin/bash
# Training augustus for P. catcorum, P. idaei and P. fragariae

#qsub /home/armita/git_repos/emr_repos/scripts/alternaria/gene_pred/Gene_pred_pipe.sh raw_rna/genbank/mixed_species/mixed_strains/phytophthora_ESTs.fasta repeat_masked/P.cactorum/404/404_assembly.51_repmask/404_contigs_hardmasked.fa


fastq-dump raw_rna/genbank/P.cactorum/10300/SRR1206032.sra
mv SRR1206032.fastq raw_rna/genbank/P.cactorum/10300/.
fastq-dump raw_rna/genbank/P.cactorum/10300/SRR1206033.sra 
mv SRR1206033.fastq raw_rna/genbank/P.cactorum/10300/.
fastq-dump --split-3 raw_rna/genbank/P.cactorum/10300/SRR1206032.sra -O raw_rna/genbank/P.cactorum/10300/
fastq-dump --split-3 raw_rna/genbank/P.cactorum/10300/SRR1206033.sra -O raw_rna/genbank/P.cactorum/10300/
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh raw_rna/genbank/P.cactorum/10300/SRR1206032.fastq raw_rna/genbank/P.cactorum/10300/SRR1206033.fastq /home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa 

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/transcriptome_assembly/transcriptome_assembly_trinity.sh qc_/paired/genbank/P.cactorum/F/P.cactorum_qc_F.fastq qc_/paired/genbank/P.cactorum/R/P.cactorum_qc_R.fastq

/home/armita/git_repos/emr_repos/scripts/alternaria/gene_pred/Gene_pred_pipe.sh assembly/trinity/paired/P.cactorum/10300_rna_contigs/Trinity.fasta repeat_masked/P.cactorum/10300/version1_repmask/10300_contigs_unmasked.fa