#Identification of motifs & domains in Phytophthora proteomes using Hmmer
=====
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

```shell
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls analysis/rxlr_atg/P.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/hmmer/CRN/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_WY_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep -B500 'inclusion threshold' | tail -n +15 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
	done
```
<!-- 
	P.cactorum 10300
	486
	257
	P.cactorum 404
	486
	218
	P.cactorum 411
	52
	18
	P.cactorum 414
	486
	230
	P.fragariae JHVZ02
	486
	282
	P.fragariae SCRP245
	1
	2
	P.ideai 371
	486
	220
 -->


#Prediction of proteins with Crinkler domains in Phytophthora proteomes

A hmm model to detect crinklers has been developed in previous studies for P. infestans.
The CRN model is a composite of three models. The first model describes the 
recombination domain containing a LxLFLAK motif in the first 60 aa of the protein.
The second describes a 'DI' domain that follows the recombination domain. The third
model describes a HVLVVVP motif in the protein known as the DWL domain.
```shell
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
	for Proteome in $(ls analysis/rxlr_atg/P.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/hmmer/CRN/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_CRN_hmmer_out.txt
		hmmsearch $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep -B500 'inclusion threshold' | tail -n +15 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
	done
```

<!-- 
	P.cactorum 10300
	120
	47
	P.cactorum 404
	128
	23
	P.cactorum 411
	10
	8
	P.cactorum 414
	147
	44
	P.fragariae JHVZ02
	192
	46
	P.fragariae SCRP245
	0
	8
	P.ideai 371
	105
	69
 -->