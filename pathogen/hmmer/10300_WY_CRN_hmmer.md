#10300 Hmmer commands
==========
These commands were used to search Phytophthora proteins 
for characteristic motifs and domains using hmmer.


#Discovery WY motifs.

WY motifs are found in about half of P. infestans RxLR proteins.
A hmm model has been published in previous studies to detect WY motifs in proteins.
This was run on isolates sequenced at EMR.

The hmm model has previously been created using hmmbuild on a 
multiple sequence alignment file of RxLR proteins containing WY 
domains but with the signal peptide removed.
The multiple sequence alignment would have been in stockholm format.

The predicted phytophthora proteins are searched for homology
to the hmm model.

The previously devloped hmm model was used to identify genes 
with WY domains in:
a) Predicted Augustus proteins
b) Predicted Augustus proteins containing SigP and RxLR motifs
c) Both sets of results from (a) and (b)
d) Predicted ORF fragments
e) Predicted ORF fragments containing SigP and RxLR motifs


a) Predicted Augustus proteins
```shell
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	Proteome=../gene_pred/augustus/P.cactorum/10300/10300_augustus_preds.aa
	HmmResults=10300_Aug_WY_hmmer_out.txt
	hmmsearch $HmmModel $Proteome > $HmmResults
```

The number of proteins that passed the threshold for WY domains were:
```shell
	cat $HmmResults | grep -B500 'inclusion threshold' | tail -n +15 | head -n -1 | wc -l
```
The P. cactorum 10300 proteome contained 75 proteins that passes the WY-domain threshold. 

The number of genes that tested +ve but didn't pass the
inclusion threshold were counted using the following command:
```shell
	cat $HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
```
A further 18 genes tested +ve but didn't pass the
inclusion threshold.

b) Predicted Augustus proteins containing SigP and RxLR motifs

```shell
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	Proteome=../analysis/sigP_rxlr/P.cactorum/10300/10300_sigP_RxLR.fa
	HmmResults=10300_Aug_RxLR_WY_hmmer_out.txt
	hmmsearch $HmmModel $Proteome > $HmmResults
	cat $HmmResults | grep -B500 'inclusion threshold' | tail -n +15 | head -n -1 | wc -l
```
Of the 125 predicted genes containing signal peptides and RxLRs, 40 contained WY domains.

The number of genes that tested +ve but didn't pass the
inclusion threshold were counted using the following command:
```shell
	cat $HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
```
One additional gene tested +ve but didn't pass the
inclusion threshold.

c) Comparison of results from (a) and (b)
The number of proteins containing WY domains in the 10300
proteome was compared to the number detected in 10300 proteins
containing SigP & RxLRs.
The following command was used to identify genes that were present in 
both outputs.
```shell
	for File in $(ls *_Aug_*WY_hmmer_out.txt); do 
        cat $File | grep -B500 'inclusion threshold' | tail -n +15 | head -n -1 | cut -f1 | grep -E -w -o 'g.* ' | sed 's/ //g'; 
    done | sort | uniq -d | wc -l
```
This showed that 29 genes were present in both outputs, passing the threshold.
```shell
	for File in $(ls *_Aug_*WY_hmmer_out.txt); do 
		cat $HmmResults | tail -n +15 | grep -B500 'Domain annotation for each sequence' | grep -v 'inclusion threshold' | head -n -3 | cut -f1 | grep -E -w -o 'g.* ' | sed 's/ //g'; 
	done | sort | uniq -d | wc -l
```
This showed that 41 genes were present in both outputs, irrespective of threshold.


d) Predicted ORF fragments
ORF fragments predicted from the atg.pl script 
were searched for WY domains:
```shell
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	Proteome=../analysis/rxlr_atg/P.cactorum/10300/10300.aa_cat.fa
	HmmResults=10300_ORF_WY_hmmer_out.txt
	hmmsearch $HmmModel $Proteome > $HmmResults
	cat $HmmResults | grep -B500 'inclusion threshold' | tail -n +15 | head -n -1 | wc -l
```
Of the 456656 ORF fragments searched, 486 sequences
showed homology to WY models. 
Visual inspection of these hits showed many of these 
features had a neighbouring feature that also passed 
the threshold for similarity to the model.
The number of genes that tested +ve but didn't pass the
inclusion threshold were counted using the following command:
```shell
	cat $HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
```
229 ORF fragments didn't pass the inclusion threshold.

e) Predicted ORF fragments containing SigP and RxLR motifs
ORF fragments predicted from the atg.pl script that also
possessed SigP and RxLR motifs were searched for WY domains:
```shell
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	Proteome=../analysis/rxlr_atg/P.cactorum/10300/10300_sp_rxlr.fa
	HmmResults=10300_ORF_RxLR_WY_hmmer_out.txt
	hmmsearch $HmmModel $Proteome > $HmmResults
	cat $HmmResults | grep -B500 'inclusion threshold' | tail -n +15 | head -n -1 | wc -l
```
Of the 1252 ORF fragments that possess RxLR and SigP domains, 75 also contained WY domains.

The number of genes that tested +ve but didn't pass the
inclusion threshold were counted using the following command:
```shell
	cat $HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
```
9 ORF fragments didn't pass the inclusion threshold.


#Discovery of CRN motifs

A hmm model to detect crinklers has been developed in previous studies for P. infestans.
The CRN model is a composite of three models. The first model describes the 
recombination domain containing a LxLFLAK motif in the first 60 aa of the protein.
The second describes a 'DI' domain that follows the recombination domain. The third
model describes a HVLVVVP motif in the protein known as the DWL domain.

The CRN hmm model was used to screen:
a) Predicted proteins
b) Predicted crinklers (using motif searches)
c) Predicted ORF fragments.


a) Predicted proteins
The proteome as predicted by Augustus was 
searched for homology to Crinkler models.
The following commands were run:
```shell
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
	Proteome=../gene_pred/augustus/P.cactorum/10300/10300_augustus_preds.aa
	HmmResults=10300_CRN_hmmer_out.txt
	hmmsearch $HmmModel $Proteome > $HmmResults
	cat $HmmResults | grep -B500 'inclusion threshold' | tail -n +15 | head -n -1 | wc -l
```
Of 17801 proteins predicted this search returned
64 proteins with homology to Crinkler models.

The number of genes that tested +ve but didn't pass the
inclusion threshold were counted using the following command:
```shell
	cat $HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
```
5 genes didn't pass the inclusion threshold.

b) Predicted crinklers
Crinklers predicted using searches for HVLVVVP and 
LxFLAK motifs were searched using the crinkler model.
```shell
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
	Proteome=../analysis/motif_search/P.cactorum/10300/10300_LxLFLAK_HVLVVVP.fa
	HmmResults=10300_Aug_CRN_hmmer_out.txt
	hmmsearch $HmmModel $Proteome > $HmmResults
```
As these results didn't return any hits with scores below
the inclusion threshold the previous command to view results
didn't work. In this case the following command was used:
```shell 
	cat $HmmResults | grep -B500 'Domain annotation for each sequence' | tail -n +15 | head -n -3 | wc -l
```
Of the 16 putative crinklers, 16 of these had support from the hmm model.


c) Predicted ORF fragments
ORF fragments as predicted by atg.pl were searched using
the Crinkler model. These sequences had not been filtered 
by presence of a SigP or any motif.
```shell
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
	Proteome=../analysis/rxlr_atg/P.cactorum/10300/10300.aa_cat.fa
	HmmResults=10300_ORF_CRN_hmmer_out.txt
	hmmsearch $HmmModel $Proteome > $HmmResults
	cat $HmmResults | grep -B500 'inclusion threshold' | tail -n +15 | head -n -1 | wc -l
```
Of the 456656 ORF fragments searched, 120 sequences
showed homology to Crinkler models. 
Visual inspection of these hits showed many of these 
features had a neighbouring feature that also passed 
the threshold for similarity to the model.

There were a number of ORF fragments that didn't pass 
the inclusion threshold. These were counted using the following command:
```shell
	cat $HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
```
47 genes were +ve but did not pass the inclusion threshold. 

