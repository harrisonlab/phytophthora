# Assess publicly available sequences


The following genomes were publicly available

```shell
	ls assembly/external_group/*/*/dna/*.genome.fa
```
```
	assembly/external_group/P.fragariae/309-62/dna/P.fragariae_309.62.genome.fa
	assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.fa
	assembly/external_group/P.kernoviae/00238-432/dna/Phytophthora_kernoviae.GCA_000333075.1.26.dna.genome.fa
	assembly/external_group/P.lateralis/MPF4/dna/Phytophthora_lateralis.GCA_000318465.1.26.dna.genome.fa
	assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra_310.i2.scaffolds.genome.fa
	assembly/external_group/P.parisitica/329/dna/phytophthora_parasitica_329.assembly.genome.fa
	assembly/external_group/P.parisitica/chvinca01/dna/phyt_para_chvinca01.i1.scaffolds.genome.fa
	assembly/external_group/P.parisitica/cj01a1/dna/phyt_para_cj01a1.1.scaffolds.genome.fa
	assembly/external_group/P.parisitica/cj02b3/dna/phyt_para_cj02b3.i1.scaffolds.genome.fa
	assembly/external_group/P.parisitica/cj05e6/dna/phyt_para_cj05e6.i1.scaffolds.genome.fa
	assembly/external_group/P.parisitica/iac_01_95/dna/phyt_para_iac_01_95.i1.scaffolds.genome.fa
	assembly/external_group/P.parisitica/p10297/dna/phyt_para_p10297.1.scaffolds.genome.fa
	assembly/external_group/P.parisitica/p1569/dna/phyt_para_p1569.1.scaffolds.genome.fa
	assembly/external_group/P.parisitica/p1976/dna/phyt_para_p1976.1.scaffolds.genome.fa
	assembly/external_group/P.ramorum/164328/dna/Phytophthora_ramorum.ASM14973v1.26.dna.genome.fa
	assembly/external_group/P.sojae/67593/dna/Phytophthora_sojae.ASM14975v1.26.dna.genome.fa
```

# assess gene space in assemblies
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	for Genome in $(ls assembly/external_group/*/*/dna/*.genome.fa); do
		echo $Genome
		qsub $ProgDir/sub_cegma.sh $Genome dna; 
	done
```

To run the path pipe script all spaces and pipe symbols had to be removed
from the headers of fasta files. This was performed using the following commands: 

```shell
	for File in $(ls assembly/external_group/P.*/T30-4/dna/*.genome.fa); do 
		OutFile=$(echo $File | sed 's/.fa/.parsed.fa/g'); 
		echo $OutFile; cat $File | sed 's/ /_/g' | sed 's/|/_/g' > $OutFile; 
	done
```

#Augustus gene prediction

## Gene prediction

```shell
	for Genome in $(ls assembly/external_group/*/*/dna/*.genome.parsed.fa); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
		ConcatRNA=qc_rna/paired/genbank/P.cactorum/10300_genbank_appended.fastq
		GeneModel=P.cactorum_10300
		echo $Genome
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRNA $GeneModel
	done
```


##Predict secreted proteins

Proteins carrying secretion signals were predicted from Augustus gene models.
This approach used SignalP 3.0. The commands used are shown below:


```
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls gene_pred/augustus/P.*/*/*_augustus_preds.aa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		InName="$Organism""_$Strain""_proteins.fa"
		SplitDir=gene_pred/sigP_aug/$Organism/$Strain
		mkdir -p $SplitDir
		cp $Proteome $SplitDir/$InName
		cd $SplitDir
		$SplitfileDir/splitfile_500.pl $InName
		rm $InName
		cd $CurPath
		for file in $(ls $SplitDir/*_split*); do
			Jobs=$(qstat | grep 'pred_sigP' | wc -l)
			while [ $Jobs -ge 32 ]; do
				sleep 10
				printf "."
				Jobs=$(qstat | grep 'pred_sigP' | wc -l)
			done
			printf "\n"
			echo $file
			qsub $ProgDir/pred_sigP.sh $file
		done
	done
```


Note - these commands include a while loop that delays submission of pred_sigP.sh 
jobs to the SGE if there are 32 or more pred_sigP.sh jobs already running. The commands
telling it to do this are:

```shell
	Jobs=$(qstat | grep 'pred_sigP' | wc -l)
	while [ $Jobs -ge 32 ]; do
		sleep 10
		printf "."
		Jobs=$(qstat | grep 'pred_sigP' | wc -l)
	done
```


The batch files of predicted proteins needed to be combined into a single file for each strain.
This was done with the following commands:
```shell
	for SplitDir in $(ls -d gene_pred/sigP_aug/P.*/*); do
		Strain=$(echo $SplitDir | cut -d '/' -f4)
		Organism=$(echo $SplitDir | cut -d '/' -f3)
		InStringAA=''
		InStringNeg=''
		InStringTab=''
		InStringTxt=''
		for GRP in $(ls -l $SplitDir/*.fa_split_* | rev | cut -d '_' -f1 | rev | sort -n); do  
			InStringAA="$InStringAA gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_proteins.fa_split_$GRP""_sp.aa";  
			InStringNeg="$InStringNeg gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_proteins.fa_split_$GRP""_sp_neg.aa";  
			InStringTab="$InStringTab gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_proteins.fa_split_$GRP""_sp.tab"; 
			InStringTxt="$InStringTxt gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_proteins.fa_split_$GRP""_sp.txt";  
		done
		cat $InStringAA > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_sp.aa
		cat $InStringNeg > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
		tail -n +2 -q $InStringTab > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_sp.tab
		cat $InStringTxt > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_sp.txt
	done
	cp -r gene_pred/sigP/* gene_pred/sigP_aug
```


## RxLR prediction

###Motif searching

RxLRs were predicted using the program rxlr_finder.py. This program uses a regular
expression to identify proteins that contain an RxLR motif within the amino acids
between the signal peptide cleavage site and 100aa downstream:

```shell
	for Pathz in $(ls gene_pred/sigP_aug/P.*/*/*_aug_sp.aa); do 
		ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/rxlr; 
		Strain=$(echo $Pathz | cut -d '/' -f4); 
		Organism=$(echo $Pathz | cut -d '/' -f3) ; 
		OutDir=analysis/sigP_rxlr/"$Organism"/"$Strain"; 
		mkdir -p $OutDir; 
		printf "\nstrain: $Strain\tspecies: $Organism\n"; 
		printf "the number of SigP gene is:\t"; 
		cat $Pathz | grep '>' | wc -l; 
		printf "the number of SigP-RxLR genes are:\t"; 
		$ProgDir/rxlr_finder.py $Pathz > $OutDir/"$Strain"_aug_RxLR_finder.fa; 
		cat $OutDir/"$Strain"_aug_RxLR_finder.fa | grep '>' | wc -l; 
	done
```

###Domain searching

Hmm models for the WY domain contained in many RxLRs were used to search gene
models predicted with Augustus. These were run with the following commands:

```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls gene_pred/augustus/P.*/*/*_augustus_preds.aa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/hmmer/WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_aug_WY_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep -B500 'inclusion threshold' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
		HmmFasta="$Strain"_aug_WY_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta	
	done
```


##Crinkler Prediction

###Motif identification

Crinkler motifs in Augustus gene models were identified by searching for 
presence of two Crinkler motifs. This was performed using the following 
commands:

```shell
	for Proteome in $(ls gene_pred/augustus/P.*/*/*_augustus_preds.aa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/CRN/$Organism/$Strain
		mkdir -p $OutDir
		OutFile=$OutDir/"$Strain"_aug_LxLFLAK_HVLVVVP.fa
		cat $Proteome | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -B1 "L.LFLAK" | grep -B1 "HVLVVVP" | grep -v '^\-\-' | sed "s/>/>$Strain\_/g" > $OutFile
		printf "Number of Crinklers in $Organism $Strain\t"
		cat $OutFile | grep '>' | wc -l
	done
```

###Domain searching

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:


```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
	for Proteome in $(ls gene_pred/augustus/P.*/*/*_augustus_preds.aa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/hmmer/CRN/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_aug_CRN_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep -B500 'inclusion threshold' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
		HmmFasta="$Strain"_aug_CRN_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta	
	done
```

## On Repeatmasked Genomes

The number of genes predicted in T-30 was considerably more than in publihsed gene models.
To determine if this was a result of repetative sequences affecting the gene models, gene 
prediction was performed on the repeat-masked dataset for P. infestans.

Firstly, the genes predicted using the unmasked genome were moved to a new directory.
```shell
	mv gene_pred/augustus/P.infestans/T30-4 gene_pred/augustus/P.infestans/T30-4_unmasked
```

The repeatmasked genome was unziped and parsed using the following commands:

```shell
	cd /home/groups/harrisonlab/project_files/idris
	gunzip assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna_rm.genome.fa.gz
```
Gene prediction was then performed using Augustus 
(using a P.cactorum gene model):
```shell
	for Genome in $(ls assembly/external_group/*/T30-4/dna/*.genome.parsed.fa); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
		ConcatRNA=qc_rna/paired/genbank/P.cactorum/10300_genbank_appended.fastq
		GeneModel=P.cactorum_10300
		echo $Genome
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRNA $GeneModel
	done
```

Signal peptides were predicted within these gene models using:

###Predict secreted proteins

Proteins carrying secretion signals were predicted from Augustus gene models.
This approach used SignalP 3.0. Firstly the previous signal peptides predicted
from unmasked data needed to be moved.

```shell
 mv gene_pred/sigP/P.infestans/T30-4 gene_pred/sigP/P.infestans/T30-4_unmasked
```

The commands to perform SigP prediction are shown below:

```shell
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls gene_pred/augustus/P.*/T30-4/*_augustus_preds.aa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		InName="$Organism""_$Strain""_proteins.fa"
		SplitDir=gene_pred/sigP_aug/$Organism/$Strain
		mkdir -p $SplitDir
		cp $Proteome $SplitDir/$InName
		cd $SplitDir
		$SplitfileDir/splitfile_500.pl $InName
		rm $InName
		cd $CurPath
		for file in $(ls $SplitDir/*_split*); do
			Jobs=$(qstat | grep 'pred_sigP' | wc -l)
			while [ $Jobs -ge 32 ]; do
				sleep 10
				printf "."
				Jobs=$(qstat | grep 'pred_sigP' | wc -l)
			done
			printf "\n"
			echo $file
			qsub $ProgDir/pred_sigP.sh $file
		done
	done
```

###RXLR prediction (repeat masked genome)

Firstly, the RxLR predictions that were made from unmasked genomes were moved
```shell
	mv analysis/sigP_rxlr/P.infestans/T30-4 analysis/sigP_rxlr/P.infestans/T30-4_unmasked
	mv analysis/hmmer/WY/P.infestans/T30-4 analysis/hmmer/WY/P.infestans/T30-4_unmasked
```
RxLR prediction was performed using the commands:

####Motif searching

RxLRs were predicted using the program rxlr_finder.py. This program uses a regular
expression to identify proteins that contain an RxLR motif within the amino acids
between the signal peptide cleavage site and 100aa downstream:

```shell
	for Pathz in $(ls gene_pred/sigP_aug/P.*/T30-4/*_aug_sp.aa); do 
		ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/rxlr; 
		Strain=$(echo $Pathz | cut -d '/' -f4); 
		Organism=$(echo $Pathz | cut -d '/' -f3) ; 
		OutDir=analysis/sigP_rxlr/"$Organism"/"$Strain"; 
		mkdir -p $OutDir; 
		printf "\nstrain: $Strain\tspecies: $Organism\n"; 
		printf "the number of SigP gene is:\t"; 
		cat $Pathz | grep '>' | wc -l; 
		printf "the number of SigP-RxLR genes are:\t"; 
		$ProgDir/rxlr_finder.py $Pathz > $OutDir/"$Strain"_aug_RxLR_finder.fa; 
		cat $OutDir/"$Strain"_aug_RxLR_finder.fa | grep '>' | wc -l; 
	done
```

####Domain searching

Hmm models for the WY domain contained in many RxLRs were used to search gene
models predicted with Augustus. These were run with the following commands:

```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls gene_pred/augustus/P.*/T30-4/*_augustus_preds.aa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/hmmer/WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_aug_WY_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep -B500 'inclusion threshold' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
		HmmFasta="$Strain"_aug_WY_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta	
	done
```

###Crinkler prediction (repeat masked genome)

Firstly, the CRN predictions that were made from unmasked genomes were moved
```shell
	mv analysis/CRN/P.infestans/T30-4 analysis/CRN/P.infestans/T30-4_unmasked
	mv analysis/hmmer/CRN/P.infestans/T30-4 analysis/hmmer/CRN/P.infestans/T30-4_unmasked
```

Crinkler prediction was performed using the commands:

####Motif identification

Crinkler motifs in masked P. infestans Augustus gene models were identified by searching for 
presence of two Crinkler motifs. This was performed using the following 
commands:

```shell
	for Proteome in $(ls gene_pred/augustus/P.*/T30-4/*_augustus_preds.aa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/CRN/$Organism/$Strain
		mkdir -p $OutDir
		OutFile=$OutDir/"$Strain"_aug_LxLFLAK_HVLVVVP.fa
		cat $Proteome | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -B1 "L.LFLAK" | grep -B1 "HVLVVVP" | grep -v '^\-\-' | sed "s/>/>$Strain\_/g" > $OutFile
		printf "Number of Crinklers in $Organism $Strain\t"
		cat $OutFile | grep '>' | wc -l
	done
```

###Domain searching

A hmm model relating to crinkler domains was used to identify putative crinklers
in masked P. infestans Augustus gene models. This was done with the following commands:


```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
	for Proteome in $(ls gene_pred/augustus/P.*/T30-4/*_augustus_preds.aa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/hmmer/CRN/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_aug_CRN_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep -B500 'inclusion threshold' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
		HmmFasta="$Strain"_aug_CRN_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta	
	done
```


# atg.pl path pipe ORF Prediction


##path pipe & RxLR motif prediction

The RxLR path pipe was run, predicting open reading frames in the genome, 
identifying ORFs with signal peptides and identifying those proteins that
also carry an RxLR motif.


As described earlier, to run the path pipe script all spaces and pipe symbols had to be removed
from the headers of fasta files. This was performed using the following commands: 

```shell
	for File in $(ls assembly/external_group/P.*/*/dna/*.genome.fa); do 
		OutFile=$(echo $File | sed 's/.fa/.parsed.fa/g'); 
		echo $OutFile; cat $File | sed 's/ /_/g' | sed 's/|/_/g' > $OutFile; 
	done
```

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing
Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen
	for Genome in $(ls assembly/external_group/*/*/dna/*.genome.parsed.fa); do
		echo $Genome
		qsub $ProgDir/path_pipe.sh $Genome 
	done
```

##RxLR Prediction

The number of ORF fragments, ORFs containing SigPs and ORFs containing 
SigP & RxLR motifs were identified.

```shell
	for OrfDir in $(ls -d analysis/rxlr_atg/P.*/*); do
		OrfFrags=$(ls $OrfDir/*.aa_cat.fa)
		SigpFrags=$(ls $OrfDir/*.sp.pve)
		RxlrFrags=$(ls $OrfDir/*_sp_rxlr.fa)
		printf "Now in:\t"
		printf "$OrfDir\n"
		printf "The number of ORF fragments are:\t"
		cat $OrfFrags | grep '>' | wc -l
		printf "The number conatining SigPs are:\t"
		cat $SigpFrags | grep '>' | wc -l
		printf "The number containing SigPs & RxLRs are:\t"
		cat $RxlrFrags | grep '>' | wc -l
	done
```


###Motif searching

To ensure that RxLR motifs were predicted using the same method 
for ORF fragments and released gene models RxLRs were predicted
using the program rxlr_finder.py:

```shell
	for Pathz in $(ls analysis/rxlr_atg/P.*/*/*.sp.pve); do 
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/rxlr; 
	Strain=$(echo $Pathz | cut -d '/' -f4); 
	Organism=$(echo $Pathz | cut -d '/' -f3) ; 
	OutDir=analysis/rxlr_atg/"$Organism"/"$Strain"; 
	mkdir -p $OutDir; 
	printf "\nstrain: $Strain\tspecies: $Organism\n"; 
	printf "the number of SigP gene is:\t"; 
	cat $Pathz | grep '>' | wc -l; 
	printf "the number of SigP-RxLR genes are:\t"; 
	$ProgDir/rxlr_finder.py $Pathz > $OutDir/"$Strain"_RxLR_finder.fa; 
	cat $OutDir/"$Strain"_RxLR_finder.fa | grep '>' | wc -l; 
	done
```


###Domain searching


```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls analysis/rxlr_atg/P.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/hmmer/WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_WY_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep -B500 'inclusion threshold' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
		HmmFasta="$Strain"_ORF_WY_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta	
	done
```

##Crinkler Prediction

###Motif identification

```shell
	for Proteome in $(ls analysis/rxlr_atg/P.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/CRN/$Organism/$Strain
		mkdir -p $OutDir
		OutFile=$OutDir/"$Strain"_ORF_LxLFLAK_HVLVVVP.fa
		cat $Proteome | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -B1 "L.LFLAK" | grep -B1 "HVLVVVP" | grep -v '\-\-' | sed "s/>/>$Strain\_/g" > $OutFile
		printf "Number of Crinklers in $Organism $Strain\t"
		cat $OutFile | grep '>' | wc -l
	done
```

###Domain searching


```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
	for Proteome in $(ls analysis/rxlr_atg/P.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/hmmer/CRN/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_CRN_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep -B500 'inclusion threshold' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
		HmmFasta="$Strain"_ORF_CRN_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta	
	done
```





#Assess gene predicted gene models

##Predict secreted proteins

```
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls assembly/external_group/P.*/*/*/*.pep.all.fa); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		InName="$Organism""_$Strain""_proteins.fa"
		SplitDir=gene_pred/sigP/$Organism/$Strain
		mkdir -p $SplitDir
		cp $Proteome $SplitDir/$InName
		cd $SplitDir
		$SplitfileDir/splitfile_500.pl $InName
		rm $InName
		cd $CurPath
		for file in $(ls $SplitDir/*_split*); do 
			echo $file
			qsub $ProgDir/pred_sigP.sh $file
		done
	done
```

The batch files of predicted proteins needed to be combined into a single file for each strain.
This was done with the following commands:
```shell
	for Pathz in $(ls -d gene_pred/SigP/*/*); do
		Strain=$(echo $Pathz | cut -d '/' -f4)
		Organism=$(echo $Pathz | cut -d '/' -f3)
		InStringAA=''
		InStringNeg=''
		InStringTab=''
		InStringTxt=''
		for GRP in $(ls -l gene_pred/SigP/$Organism/$Strain/*.fa_split_* | cut -d '_' -f6 | sort -n); do  
			InStringAA="$InStringAA gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_proteins.fa_split_$GRP""_sp.aa";  
			InStringNeg="$InStringNeg gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_proteins.fa_split_$GRP""_sp_neg.aa";  
			InStringTab="$InStringTab gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_proteins.fa_split_$GRP""_sp.tab"; 
			InStringTxt="$InStringTxt gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_proteins.fa_split_$GRP""_sp.txt";  
		done
		cat $InStringAA > gene_pred/sigP/$Organism/$Strain/"$Strain"_sp.aa
		cat $InStringNeg > gene_pred/sigP/$Organism/$Strain/"$Strain"_neg_sp.aa
		tail -n +2 -q $InStringTab > gene_pred/sigP/$Organism/$Strain/"$Strain"_sp.tab
		cat $InStringTxt > gene_pred/sigP/$Organism/$Strain/"$Strain"_sp.txt
	done
```

# Predict pathogenicity genes in published gene models

##RxLRs
```shell
	for Pathz in $(ls -d gene_pred/sigP/*/*); do
		ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/rxlr
		Strain=$(echo $Pathz | cut -d '/' -f4)
		Organism=$(echo $Pathz | cut -d '/' -f3) 
		OutDir=analysis/sigP_rxlr/"$Organism"/"$Strain"
		mkdir -p $OutDir
		printf "\nstrain: $Strain\tspecies: $Organism\n"
		printf "the number of SigP gene is:\t"
		cat $Pathz/"$Strain"_sp.aa | grep '>' | wc -l
		printf "the number of SigP-RxLR genes are:\t"
		$ProgDir/rxlr_finder.py $Pathz/"$Strain"_sp.aa > $OutDir/"$Strain"_sp_RxLR.fa
		cat $OutDir/"$Strain"_sp_RxLR.fa | grep '>' | wc -l
	done
```

###Motif identification

RxLR motifs were predicted using the program rxlr_finder.py:

```shell
	for Pathz in $(ls -d gene_pred/sigP/*/*); do 
		ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/rxlr; 
		Strain=$(echo $Pathz | cut -d '/' -f4); 
		Organism=$(echo $Pathz | cut -d '/' -f3) ; 
		OutDir=analysis/sigP_rxlr/"$Organism"/"$Strain"; 
		mkdir -p $OutDir; 
		printf "\nstrain: $Strain\tspecies: $Organism\n"; 
		printf "the number of SigP gene is:\t"; 
		cat $Pathz/"$Strain"_sp.aa | grep '>' | wc -l; 
		printf "the number of SigP-RxLR genes are:\t"; 
		$ProgDir/rxlr_finder.py $Pathz/"$Strain"_sp.aa > $OutDir/"$Strain"_sp_RxLR.fa; 
		cat $OutDir/"$Strain"_sp_RxLR.fa | grep '>' | wc -l; 
	done
```


###Domain searching


```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls assembly/external_group/P.*/*/*/*.pep.all.fa); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=analysis/hmmer/WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_Published_WY_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep -B500 'inclusion threshold' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
		HmmFasta="$Strain"_Published_WY_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta	
	done
```


##Crinklers

###Motif identification

```shell
	for Proteome in $(ls assembly/external_group/P.*/*/*/*.pep.all.fa); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=analysis/CRN/$Organism/$Strain
		mkdir -p $OutDir
		OutFile=$OutDir/"$Strain"_LxLFLAK_HVLVVVP.fa
		cat $Proteome | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -B1 "L.LFLAK" | grep -B1 "HVLVVVP" | grep -v '\-\-' | sed "s/>/>$Strain\_/g" > $OutFile
		printf "Number of Crinklers in $Organism $Strain\t"
		cat $OutFile | grep '>' | wc -l
	done
```

###Domain searching


```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
	for Proteome in $(ls assembly/external_group/P.*/*/*/*.pep.all.fa); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=analysis/hmmer/CRN/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_Published_CRN_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep -B500 'inclusion threshold' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | grep -A500 'inclusion threshold' | grep -B500 'Domain annotation for each sequence' | tail -n +2 | head -n -3 | wc -l
		HmmFasta="$Strain"_Published_CRN_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta	
	done
```



# Annotating pathogenicity genes


##Extract features from Augustus predicitons

###RxLR motif predictions

Gff features for RxLRs were extracted from Augustus gff output
files using a list of names of putative pathogenicity genes. This list was built
using the following commands:

```shell
	for RxLR_File in $(ls analysis/sigP_rxlr/P*/*/*_aug_RxLR_finder.fa); do
		Organism=$(echo $RxLR_File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $RxLR_File | rev | cut -f2 -d '/' | rev)
		echo $Strain
		OutFile=analysis/sigP_rxlr/"$Organism"/"$Strain"/"$Strain"_aug_RxLR_finder_names.txt
		cat $RxLR_File | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutFile
	done
```

This list was then used to extract gff features using the program gene_list_to_gff.pl.
The commands used to run this were:

Note: it was realised that the Augustus gff features were in gff3 format rather 
than gtf format.

```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Col2=rxlr_finder.py
	for GeneNames in $(ls analysis/sigP_rxlr/P.*/*/*_aug_RxLR_finder_names.txt); do
		Organism=$(echo $GeneNames | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $GeneNames | rev | cut -f2 -d '/' | rev)
		echo "$Strain"
		GeneModels=gene_pred/augustus/"$Organism"/"$Strain"/*_augustus_preds.gtf
		OutFile=analysis/sigP_rxlr/"$Organism"/"$Strain"/"$Strain"_aug_RxLR_finder.gff
		$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 ID > $OutFile
	done
```

###RxLR WY domain predictions


Gff features for WY domains-containing proteins were extracted from Augustus gff output
files using a list of names of putative pathogenicity genes. This list was built
using the following commands:

```shell
	for WY_File in $(ls analysis/hmmer/WY/P*/*/*_aug_WY_hmmer_out.fa); do
		Organism=$(echo $WY_File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $WY_File | rev | cut -f2 -d '/' | rev)
		echo $Strain
		OutFile=analysis/hmmer/WY/"$Organism"/"$Strain"/"$Strain"_aug_WY_hmmer_names.txt
		cat $WY_File | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutFile
	done
```

This list was then used to extract gff features using the program gene_list_to_gff.pl.
The commands used to run this were:

Note: it was realised that the Augustus gff features were in gff3 format rather 
than gtf format.

```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Col2=WY_hmmer
	for GeneNames in $(ls analysis/hmmer/WY/P.*/*/*_aug_WY_hmmer_names.txt); do
		Organism=$(echo $GeneNames | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $GeneNames | rev | cut -f2 -d '/' | rev)
		echo "$Strain"
		GeneModels=gene_pred/augustus/"$Organism"/"$Strain"/*_augustus_preds.gtf
		OutFile=analysis/hmmer/WY/"$Organism"/"$Strain"/"$Strain"_aug_WY_hmmer.gff
		$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 ID > $OutFile
	done
```

###Crinkler motif predictions

Gff features for Crinkler motif-containing proteins were extracted from Augustus gff output
files using a list of names of putative pathogenicity genes. This list was built
using the following commands:

Note: Unlike the fasta files above, these accessions have been renamed with the 
name of the strain preceeding the gene ID. As such, addtional steps have been added
to filter the gene names.

```shell
	for MotifCRN_File in $(ls analysis/CRN/P*/*/*_aug_LxLFLAK_HVLVVVP.fa); do
		Organism=$(echo $MotifCRN_File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $MotifCRN_File | rev | cut -f2 -d '/' | rev)
		echo $Strain
		OutFile=analysis/CRN/"$Organism"/"$Strain"/"$Strain"_aug_LxLFLAK_HVLVVVP_names.txt
		cat $MotifCRN_File | grep '>' | cut -f1 | sed 's/>//g' | rev | cut -f1 -d '_' | rev | sed 's/ //g' > $OutFile
	done
```

This list was then used to extract gff features using the program gene_list_to_gff.pl.
The commands used to run this were:

Note: it was realised that the Augustus gff features were in gff3 format rather 
than gtf format.

```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Col2=CRN_motif
	for GeneNames in $(ls analysis/CRN/*/*/*_aug_LxLFLAK_HVLVVVP_names.txt); do
		Organism=$(echo $GeneNames | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $GeneNames | rev | cut -f2 -d '/' | rev)
		echo "$Strain"
		GeneModels=gene_pred/augustus/"$Organism"/"$Strain"/*_augustus_preds.gtf
		OutFile=analysis/CRN/"$Organism"/"$Strain"/"$Strain"_aug_LxLFLAK_HVLVVVP.gff
		$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 ID > $OutFile
	done
```

###Crinkler motif predictions

Gff features for Crinkler motif-containing proteins were extracted from Augustus gff output
files using a list of names of putative pathogenicity genes. This list was built
using the following commands:

```shell
	for HmmCRN_File in $(ls analysis/hmmer/CRN/P*/*/*_aug_CRN_hmmer_out.fa); do
		Organism=$(echo $HmmCRN_File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $HmmCRN_File | rev | cut -f2 -d '/' | rev)
		echo $Strain
		OutFile=analysis/hmmer/CRN/"$Organism"/"$Strain"/"$Strain"_aug_CRN_hmmer_names.txt
		cat $HmmCRN_File | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutFile
	done
```

This list was then used to extract gff features using the program gene_list_to_gff.pl.
The commands used to run this were:

Note: it was realised that the Augustus gff features were in gff3 format rather 
than gtf format.

```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Col2=CRN_hmm
	for GeneNames in $(ls analysis/hmmer/CRN/*/*/*_aug_CRN_hmmer_names.txt); do
		Organism=$(echo $GeneNames | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $GeneNames | rev | cut -f2 -d '/' | rev)
		echo "$Strain"
		GeneModels=gene_pred/augustus/"$Organism"/"$Strain"/*_augustus_preds.gtf
		OutFile=analysis/hmmer/CRN/"$Organism"/"$Strain"/"$Strain"_aug_CRN_hmmer.gff
		$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 ID > $OutFile
	done
```





##Extract features from atg.pl predictions

Gff features from atg.pl were corrected to .gff3 format. 
This was done using the following commands:

```
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	for StrainDir in $(ls -d analysis/rxlr_atg/P*/*); do
		Organism=$(echo $StrainDir | rev | cut -f2 -d '/' | rev)
		Strain=$(echo $StrainDir | rev | cut -f1 -d '/' | rev)
		GffFile=$(ls "$StrainDir"/*_ORF.gff | grep -v '_atg_')
		echo $Strain
		OutFile=analysis/rxlr_atg/"$Organism"/"$Strain"/"$Strain"_ORF.gff3
		$ProgDir/gff_corrector.pl $GffFile > $OutFile
	done
```

###RxLR motif predictions

for RxLRs were extracted from atg.pl gff output
files using a list of names of putative pathogenicity genes. This list was built
using the following commands:

```shell
	for RxLR_File in $(ls analysis/rxlr_atg/P*/*/*_sp_rxlr.fa); do
		Organism=$(echo $RxLR_File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $RxLR_File | rev | cut -f2 -d '/' | rev)
		echo $Strain
		OutFile=analysis/rxlr_atg/"$Organism"/"$Strain"/"$Strain"_sp_rxlr_names.txt
		cat $RxLR_File | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutFile
	done
```

This list was then used to extract gff features using the program gene_list_to_gff.pl.
The commands used to run this were:

Note: it was realised that the Augustus gff features were in gff3 format rather 
than gtf format.
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Col2=atg_RxLR
	for GeneNames in $(ls analysis/rxlr_atg/*/*/*_sp_rxlr_names.txt); do
		Organism=$(echo $GeneNames | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $GeneNames | rev | cut -f2 -d '/' | rev)
		echo "$Strain"
		GeneModels=analysis/rxlr_atg/"$Organism"/"$Strain"/"$Strain"_ORF.gff3
		OutFile=analysis/rxlr_atg/"$Organism"/"$Strain"/"$Strain"_ORF_sp_rxlr.gff3
		$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 Name > $OutFile
	done
```

###RxLR WY domain predictions

Gff features for proteins containing WY domains were extracted from atg.pl gff output
files using a list of names of putative pathogenicity genes. This list was built
using the following commands:

```shell
	for WY_File in $(ls analysis/hmmer/WY/P*/*/*_ORF_WY_hmmer_out.fa); do
		Organism=$(echo $WY_File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $WY_File | rev | cut -f2 -d '/' | rev)
		echo $Strain
		OutFile=analysis/hmmer/WY/"$Organism"/"$Strain"/"$Strain"_ORF_WY_hmmer_names.txt
		cat $WY_File | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutFile
	done
```

This list was then used to extract gff features using the program gene_list_to_gff.pl.
The commands used to run this were:

Note: it was realised that the Augustus gff features were in gff3 format rather 
than gtf format.
$ProgDir/gene_list_to_gff.pl out_names.txt out2.gff atg_RxLR Name > out3.gff
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Col2=atg_RxLR
	for GeneNames in $(ls analysis/hmmer/WY/P*/*/*_ORF_WY_hmmer_names.txt); do
		Organism=$(echo $GeneNames | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $GeneNames | rev | cut -f2 -d '/' | rev)
		echo "$Strain"
		GeneModels=analysis/rxlr_atg/"$Organism"/"$Strain"/"$Strain"_ORF.gff3
		OutFile=analysis/hmmer/WY/"$Organism"/"$Strain"/"$Strain"_ORF_WY_hmmer.gff3
		$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 Name > $OutFile
	done
```
###Crinkler motif predictions

Putative Crinkler motif containing ORFs were extracted from atg.pl gff output
files using a list of names of putative pathogenicity genes. This list was built
using the following commands:

```shell
	for CRN_File in $(ls analysis/CRN/P*/*/*_ORF_LxLFLAK_HVLVVVP.fa); do
		Organism=$(echo $CRN_File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $CRN_File | rev | cut -f2 -d '/' | rev)
		echo $Strain
		OutFile=analysis/CRN/"$Organism"/"$Strain"/"$Strain"_ORF_LxLFLAK_HVLVVVP_names.txt
		cat $CRN_File | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutFile
	done
```

This list was then used to extract gff features using the program gene_list_to_gff.pl.
The commands used to run this were:

Note: it was realised that the Augustus gff features were in gff3 format rather 
than gtf format.
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Col2=atg_CRN
	for GeneNames in $(ls analysis/CRN/*/*/*_ORF_LxLFLAK_HVLVVVP_names.txt); do
		Organism=$(echo $GeneNames | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $GeneNames | rev | cut -f2 -d '/' | rev)
		echo "$Strain"
		GeneModels=analysis/rxlr_atg/"$Organism"/"$Strain"/"$Strain"_ORF.gff3
		OutFile=analysis/CRN/"$Organism"/"$Strain"/"$Strain"_ORF_LxLFLAK_HVLVVVP.gff3
		$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 Name > $OutFile
	done
```

###RxLR CRN domain predictions

Gff features for proteins containing crinkler domains were extracted from atg.pl gff output
files using a list of names of putative pathogenicity genes. This list was built
using the following commands:

```shell
	for CRN_File in $(ls analysis/hmmer/CRN/P*/*/*_ORF_CRN_hmmer_out.fa); do
		Organism=$(echo $CRN_File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $CRN_File | rev | cut -f2 -d '/' | rev)
		echo $Strain
		OutFile=analysis/hmmer/CRN/"$Organism"/"$Strain"/"$Strain"_ORF_CRN_hmmer_names.txt
		cat $CRN_File | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutFile
	done
```

This list was then used to extract gff features using the program gene_list_to_gff.pl.
The commands used to run this were:

Note: it was realised that the Augustus gff features were in gff3 format rather 
than gtf format.
$ProgDir/gene_list_to_gff.pl out_names.txt out2.gff atg_RxLR Name > out3.gff
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Col2=hmm_CRN
	for GeneNames in $(ls analysis/hmmer/CRN/P*/*/*_ORF_CRN_hmmer_names.txt); do
		Organism=$(echo $GeneNames | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $GeneNames | rev | cut -f2 -d '/' | rev)
		echo "$Strain"
		GeneModels=analysis/rxlr_atg/"$Organism"/"$Strain"/"$Strain"_ORF.gff3
		OutFile=analysis/hmmer/CRN/"$Organism"/"$Strain"/"$Strain"_ORF_CRN_hmmer.gff3
		$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 Name > $OutFile
	done
```




# Comparing Augustus predictions and atg.pl predictions


## atg.pl RxLR overlap

As >3000 RxLR containing ORFs were predicted in P.infestans 
using the atg.pl script, the number of overlapping features 
were identified in the gff file.

First the gff features for RxLR containing ORFs were extracted
from the gene models predicted with atg.pl

```shell
	RxlrFasta=analysis/rxlr_atg/P.infestans/T30-4/T30-4_sp_rxlr.fa
	RxlrHeaders=analysis/rxlr_atg/P.infestans/T30-4/T30-4_sp_rxlr_headers.txt
	OrfGff=analysis/rxlr_atg/P.infestans/T30-4/T30-4_ORF.gff 
	RxlrGff=analysis/rxlr_atg/P.infestans/T30-4/T30-4_sp_rxlr.gff
	cat $RxlrFasta | grep '>' | cut -f1 | sed 's/>//g' > $RxlrHeaders
	cat $OrfGff | grep -w -f $RxlrHeaders > $RxlrGff
```

The bedtools program was used to do identify overlap between gff features.
```shell
```	


# Benchmarking

RxLRs predicted in published studies were compared to those 
predicted in house on published gene models.


This was first performed on P. infestans.

Files containing the headers of putative effectors were made using the following commands:

```shell
	mkdir -p analysis/benchmarking/P.infestans/T30-4/rxlr
	Pinf_pub_RxLR=analysis/benchmarking/P.infestans/T30-4/rxlr/P.inf_annotated_RxLR_headers.txt
	Pinf_pred_RxLR=analysis/benchmarking/P.infestans/T30-4/rxlr/P.inf_pipeline_RxLR_headers.txt
	Pinf_pred_WY=analysis/benchmarking/P.infestans/T30-4/rxlr/P.inf_WY_RxLR_headers.txt
	Pinf_pred_mixed=analysis/benchmarking/P.infestans/T30-4/rxlr/P.inf_mixed_headers.txt
	cat assembly/external_group/P.infestans/T30-4/rxlr/P.inf_RxLR_parsed.csv | cut -f2 | grep 'PITG' > $Pinf_pub_RxLR
	cat analysis/sigP_rxlr/P.infestans/T30-4/T30-4_sp_RxLR.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/T0 //g' > $Pinf_pred_RxLR
	cat analysis/hmmer/WY/P.infestans/T30-4/T30-4_Published_WY_hmmer_out.txt | grep 'pep:known' | sed -e 's/ \+/\t/g' | cut -f10 | sed 's/T0//g' > $Pinf_pred_WY
	cat $Pinf_pred_RxLR $Pinf_pred_WY | sort | uniq > $Pinf_pred_mixed
```

Duplicate header names were identified between RxLR files:

```shell
	OutFile=analysis/benchmarking/P.infestans/T30-4/rxlr/Pinf_shared_RxLR_stats.txt
	printf "The number of annoated RxLRs in P.inf are:\t" > $OutFile
	cat "$Pinf_pub_RxLR" | wc -l >> $OutFile
	printf "The number of RxLRs predicted using motif searches are:\t" >> $OutFile
	cat "$Pinf_pred_RxLR" | wc -l >> $OutFile
	printf "The number of RxLRs shared between public predictions and motif searches are: \t" >> $OutFile
	cat "$Pinf_pub_RxLR" "$Pinf_pred_RxLR" | sort | uniq -d | wc -l >> $OutFile
	printf "The number of proteins with hmm model hits to WY profiles are:\t" >> $OutFile
	cat "$Pinf_pred_WY" | wc -l >> $OutFile
	printf "The number of RxLRs shared between public predictions and WY hits are: \t" >> $OutFile
	cat "$Pinf_pub_RxLR" "$Pinf_pred_WY" | sort | uniq -d | wc -l >> $OutFile
	printf "The number of RxLRs shared between all predicted RxLRs from motif and hmm searches are:\t" >> $OutFile
	cat "$Pinf_pred_mixed" | wc -l  >> $OutFile
	printf "The number of RxLRs shared between public predictions and all predicted RxLRs from motif and hmm searches are:\t" >> $OutFile
	cat "$Pinf_pub_RxLR" "$Pinf_pred_mixed" | sort | uniq -d | wc -l >> $OutFile
```
This shows that of the 454 RxLRs predicted in P.inf our RxLR motif analysis and WY domain searches also identify 434 of these.


Features of the 20 unpredicted proteins were examined. These proteins were identified by 
saving the headers of matched proteins to a new file. The original list of published RxLRs
were filtered to remove those proteins with a match.

```shell
	Pinf_confirmed_RxLR=analysis/benchmarking/P.infestans/T30-4/rxlr/P.inf_confirmed_RxLR_headers.txt
	Pinf_remainder_RxLR=analysis/benchmarking/P.infestans/T30-4/rxlr/P.inf_remainder_RxLR_headers.txt
	cat "$Pinf_pub_RxLR" "$Pinf_pred_mixed" | sort | uniq -d > $Pinf_confirmed_RxLR
	cat assembly/external_group/P.infestans/T30-4/rxlr/P.inf_RxLR_parsed.csv | grep 'PITG' | grep -v -f "$Pinf_confirmed_RxLR" > $Pinf_remainder_RxLR
```

This revealed that the remaining 20 annotated RxLRs in the Pinfestans genome did not contain
an RxLR domain, or had an RxLR variant.
