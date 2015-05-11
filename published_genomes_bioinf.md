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


#Predict gene modoels


ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
ConcatRNA=qc_rna/paired/genbank/P.cactorum/10300_genbank_appended.fastq
GeneModel=P.cactorum_10300
for Genome in $(ls assembly/external_group/*/*/dna/*.genome.fa); do echo $Genome; qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRNA $GeneModel; done





# atg.pl path pipe


##path pipe & RxLR motif prediction

```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen
	for Genome in $(ls assembly/external_group/*/*/dna/*.genome.fa); do
		echo $Genome
		qsub $ProgDir/path_pipe.sh $Genome; 
	done
```

##WY domain prediciton

##CRN prediction from motifs

##CRN prediction from domains





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
		SplitDir=gene_pred/SigP/$Organism/$Strain
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

# Predict pathogenicity genes in published gene models

##RxLRs

###Motif identification

###Domain searching


```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls assembly/external_group/P.*/*/*/*.pep.all.fa); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
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


##Crinklers

###Motif identification

###Domain searching


```shell
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/Phyt_annot_CRNs_D1.hmm
	for Proteome in $(ls assembly/external_group/P.*/*/*/*.pep.all.fa); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
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
