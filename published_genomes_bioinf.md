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
		SplitDir=gene_pred/sigP/$Organism/$Strain/split_aug
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
<!-- 

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
 -->


## RxLR prediction






# atg.pl path pipe ORF Prediction


##path pipe & RxLR motif prediction

The RxLR path pipe was run, predicting open reading frames in the genome, 
identifying ORFs with signal peptides and identifying those proteins that
also carry an RxLR motif.


To run the path pipe script all spaces and pipe symbols had to be removed
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
	for Pathz in $(ls analysis/rxlr_atg/P.*/309-62/*.sp.pve); do 
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
