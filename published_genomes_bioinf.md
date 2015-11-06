# Assess publicly available sequences


The following genomes were publicly available

```bash
	ls assembly/external_group/*/*/dna/*.genome.fa
```
```bash
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
```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	for Genome in $(ls assembly/external_group/*/*/dna/*.genome.fa); do
		echo $Genome
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```

<!--
The 310 genome's scaffolds only contained numbers stopping cegma from running properly.
The script was as ajusted as so:

```bash
	cat assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra_310.i2.scaffolds.genome.parsed.fa | sed 's/>/>NODE_/g' > assembly/external_group/P.parisitica/310/dna/phytophthora_parasitica_inra_310.i2.scaffolds.genome.parsed2.fa
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	for Genome in $(ls assembly/external_group/P.parisitica/310/dna/*_310.i2.scaffolds.genome.parsed2.fa); do
		echo $Genome
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```
 -->

To run the path pipe script all spaces and pipe symbols had to be removed
from the headers of fasta files. This was performed using the following commands:

```bash
	for File in $(ls assembly/external_group/P.*/T30-4/dna/*.genome.fa); do
		OutFile=$(echo $File | sed 's/.fa/.parsed.fa/g');
		echo $OutFile; cat $File | sed 's/ /_/g' | sed 's/|/_/g' > $OutFile;
	done
```

#Repeat masking

Rpeatmasking was performed on assemblies:


```bash
	for Genome in $(ls assembly/external_group/P.*/*/dna/*.genome.parsed.fa | grep -v 'rm'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
		qsub $ProgDir/rep_modeling.sh $Genome
		qsub $ProgDir/transposonPSI.sh $Genome
	done
```
<!--
The P infestans genome contained contigs with headers
longer than 50 characters in length. This prevented the
rep_modeling script from working. To get around this the
Fasta headers had to be corrected and the repeat masking
resubmitted.

```bash
	Genome=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed.fa
	Genome_parsed=assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna.genome.parsed2.fa
	cat $Genome | sed 's/.*supercontig_supercontig:/>/g' > $Genome_parsed
	for Genome in $(ls $Genome_parsed); do
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
		qsub $ProgDir/rep_modeling.sh $Genome
		qsub $ProgDir/transposonPSI.sh $Genome
	done
```
 -->


#Augustus gene prediction
#Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained for the 10300 repeatmasked geneome using assembled RNAseq data and predicted CEGMA genes.
	Gene prediction
		- Gene models were used to predict genes in the 10300 genome. This used RNAseq data as hints for gene models.

<!--
##Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

This was first performed on the published unmasked assemblies:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
cd /home/groups/harrisonlab/project_files/idris/
for Genome in $(ls repeat_masked/P.parisitica/310/*/*_contigs_unmasked_parsed.fa); do
echo $Genome;
qsub $ProgDir/sub_cegma.sh $Genome dna;
done
```

These results were then moved to a directory for unmasked genomes
```bash
	mv gene_pred/cegma/P.cactorum/10300 gene_pred/cegma/P.cactorum/10300_unmasked
```

The analysis was then repeated for the 10300 repeatmasked genome:
```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
cd /home/groups/harrisonlab/project_files/idris/
for Genome in $(ls repeat_masked/P.cactorum/10300/*/*_contigs_hardmasked_parsed.fa); do
echo $Genome;
qsub $ProgDir/sub_cegma.sh $Genome dna;
done
```

These results were moved to a directory for repeatmasked data.
```bash
	mv gene_pred/cegma/P.cactorum/10300 gene_pred/cegma/P.cactorum/10300_hardmasked
```

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/P.cactorum/10300*/*_dna_cegma.completeness_report); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/P.cactorum/10300_cegma_results_dna_summary.txt

	less gene_pred/cegma/P.cactorum/10300_cegma_results_dna_summary.txt
```
 -->

## Gene prediction

```bash
	for Genome in $(ls assembly/external_group/*/*/dna/*.genome.parsed.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
		ConcatRNA=qc_rna/paired/genbank/P.cactorum/10300_genbank_appended.fastq
		GeneModel=P.cactorum_10300
		echo $Genome
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRNA $GeneModel
	done
```


Gene predicted was also performed on repeatmasked assemblies


```bash
	for File in $(ls assembly/external_group/P.*/*/dna/*.dna_rm.*.fa.gz); do
		OutFile=$(echo $File | sed 's/.fa.gz/.parsed.fa/g');
		echo $OutFile;
		cat $File | gunzip -fc | sed 's/ /_/g' | sed 's/|/_/g' > $OutFile;
	done
```

```bash
	for Genome in $(ls assembly/external_group/*/*/dna/*.dna_rm.*.parsed.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
		ConcatRNA=qc_rna/paired/genbank/P.cactorum/10300_genbank_appended.fastq
		GeneModel=P.cactorum_10300
		echo $Genome
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRNA $GeneModel
	done
```

The P. parisitica genome 310 could not be used for gene prediction as it's nodes
needed to be further parsed (as they only contained numbers. These were modified
and resubmitted as follows:

```bash
	for Genome in $(ls repeat_masked/P.parisitica/310/dna_repmask/310_contigs_hardmasked.fa); do
		ModGenome=$(echo $Genome | sed 's/hardmasked.fa/hardmasked_parsed.fa/g')
		cat $Genome | sed 's/>/>Node_/g' > $ModGenome
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
		ConcatRNA=qc_rna/paired/genbank/P.cactorum/10300_genbank_appended.fastq
		GeneModel=P.cactorum_10300
		echo $Genome
		qsub $ProgDir/augustus_pipe.sh $ModGenome $ConcatRNA $GeneModel
	done
```

When these predictions had finished the output augustus directory was renamed to reflect
that these genes had predicted from masked assemblies.

```bash
	mv gene_pred/augustus gene_pred/augustus_masked
```


<!--
Gene prediction was also performed on the P. infestans T30-4 genome following
Repeatmasking by our pipeline. This

```bash
	for Genome in $(ls repeat_masked/P.infestans/T30-4/dna_repmask/T30-4_contigs_hardmasked.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
		ConcatRNA=qc_rna/paired/genbank/P.cactorum/10300_genbank_appended.fastq
		GeneModel=P.cactorum_10300
		echo $Genome
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRNA $GeneModel
		mv gene_pred/augustus gene_pred/augustus_mask_in_house
	done
```

The 13465 gene models produced from repeatmasking were compared to the 17787 gene
models from the published repeatmasked genome.

```bash
	ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
	ProjDir=/home/groups/harrisonlab/project_files/idris
	WorkDir=$ProjDir/analysis/benchmarking/P.infestans/T30-4/repeatmasked
	mkdir -p $WorkDir
	cd $WorkDir

	Gff_Pinf_published=$ProjDir/assembly/external_group/P.infestans/T30-4/pep/phytophthora_infestans_t30-4_1_transcripts.gff3
	Gff_Pinf_in_house=$ProjDir/gene_pred/augustus_mask_in_house/P.infestans/T30-4/T30-4_augustus_preds.gtf
	Gff_Pinf_in_house_mod=$WorkDir/T30-4_augustus_preds_parsed.gtf
	DB_A=$WorkDir/T30-4_published.db
	DB_B=$WorkDir/T30-4_in_house.db
	A_ID=$WorkDir/T30-4_publihsed_IDs.txt
	B_ID=$WorkDir/T30-4_inhouse_IDs.txt
	DB_A_ID=$WorkDir/T30-4_publihsed_IDs.db
	DB_B_ID=$WorkDir/T30-4_inhouse_IDs.db
	DB_merge=$WorkDir/T30-4_merged.db

	cat $Gff_Pinf_in_house | sed 's/ASM14294v1:supercont/Supercontig_/g' | sed "s/:.*\tAUGUSTUS/\tAUGUSTUS/g" > $Gff_Pinf_in_house_mod

	$ProgDir/make_gff_database.py --inp $Gff_Pinf_published --db $DB_A
	$ProgDir/make_gff_database.py --inp $Gff_Pinf_in_house_mod --db $DB_B

	$ProgDir/get_db_id.py --db $DB_A --type gene --out $A_ID
	$ProgDir/note2db.py --in_db $DB_A --out_db $DB_A_ID --id_file $A_ID --str published_gene --attribute ID
	$ProgDir/get_db_id.py --db $DB_B --type gene --out $B_ID
	$ProgDir/note2db.py --in_db $DB_B --out_db $DB_B_ID --id_file $B_ID --str in_house_gene --attribute ID

	$ProgDir/merge_db.py --inp $DB_A_ID $DB_B_ID --db $DB_merge
	$ProgDir/report_overlaps.py --inp $DB_merge --A PI_T30-4_FINAL_CALLGENES_4 --B AUGUSTUS
			#	The total number of Augustus genes are:	18179
			#	The total number of atg genes are:	13443
			#	Into this many features:	11815
	$ProgDir/report_overlaps.py --inp $DB_merge --B PI_T30-4_FINAL_CALLGENES_4 --A AUGUSTUS
			#	The total number of Augustus genes are:	13443
			#	The total number of atg genes are:	18179
			#	Into this many features:	11600
	```


$ProgDir/effectors2genemodels.sh \
	$ProjDir/gene_pred/augustus_unmasked/P.infestans/T30-4_unmasked/T30-4_augustus_preds.gtf \
	$ProjDir/analysis/rxlr_atg_unmasked/P.infestans/T30-4_unmasked/T30-4_ORF_sp_rxlr.gff3 \
	$ProjDir/analysis/hmmer/WY/P.infestans/T30-4_unmasked/T30-4_unmasked_ORF_WY_hmmer.gff3 \
	$ProjDir/analysis/sigP_rxlr/P.infestans/T30-4_unmasked/T30-4_unmasked_aug_RxLR_finder_names.txt \
	$ProjDir/analysis/hmmer/WY/P.infestans/T30-4_unmasked/T30-4_unmasked_aug_WY_hmmer_names.txt \
	T30-4_Aug_ORF_full_rxlr.db \
	T30-4_effectors.gff \
	2>&1 | tee T30-4_logfile.txt
```
```

 -->


##Predict secreted proteins

Proteins carrying secretion signals were predicted from Augustus gene models.
This approach used SignalP 3.0. The commands used are shown below:


```bash
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

```bash
	Jobs=$(qstat | grep 'pred_sigP' | wc -l)
	while [ $Jobs -ge 32 ]; do
		sleep 10
		printf "."
		Jobs=$(qstat | grep 'pred_sigP' | wc -l)
	done
```


The batch files of predicted proteins needed to be combined into a single file for each strain.
This was done with the following commands:
```bash
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

```bash
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

The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identfy RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
	for Pathz in $(ls gene_pred/sigP_aug/P.*/*/*_aug_sp.aa); do
		ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors;
		Strain=$(echo $Pathz | cut -d '/' -f4);
		Organism=$(echo $Pathz | cut -d '/' -f3) ;
		OutDir=analysis/sigP_rxlr/"$Organism"/"$Strain";
		mkdir -p $OutDir;
		printf "\nstrain: $Strain\tspecies: $Organism\n";
		printf "the number of SigP gene is:\t";
		cat $Pathz | grep '>' | wc -l;
		printf "the number of SigP-RxLR genes are:\t";
		$ProgDir/RxLR_EER_regex_finder.py $Pathz > $OutDir/"$Strain"_aug_RxLR_EER_regex.fa;
		cat $OutDir/"$Strain"_aug_RxLR_EER_regex.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_aug_RxLR_regex.txt
		cat $OutDir/"$Strain"_aug_RxLR_regex.txt | wc -l
		printf "the number of SigP-RxLR-EER genes are:\t";
		cat $OutDir/"$Strain"_aug_RxLR_EER_regex.fa | grep '>' | grep 'EER_motif_start' |  cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_aug_RxLR_EER_regex.txt
		cat $OutDir/"$Strain"_aug_RxLR_EER_regex.txt | wc -l
		printf "\n"
	done
```
Results were as follows:

```
	strain: 10300	species: P.cactorum
	the number of SigP gene is:	1856
	the number of SigP-RxLR genes are:	133
	the number of SigP-RxLR-EER genes are:	61


	strain: 404	species: P.cactorum
	the number of SigP gene is:	1655
	the number of SigP-RxLR genes are:	119
	the number of SigP-RxLR-EER genes are:	56


	strain: 414	species: P.cactorum
	the number of SigP gene is:	1688
	the number of SigP-RxLR genes are:	125
	the number of SigP-RxLR-EER genes are:	58


	strain: JHVZ02	species: P.fragariae
	the number of SigP gene is:	2098
	the number of SigP-RxLR genes are:	209
	the number of SigP-RxLR-EER genes are:	110


	strain: 371	species: P.ideai
	the number of SigP gene is:	1440
	the number of SigP-RxLR genes are:	90
	the number of SigP-RxLR-EER genes are:	44


	strain: T30-4	species: P.infestans
	the number of SigP gene is:	1841
	the number of SigP-RxLR genes are:	142
	the number of SigP-RxLR-EER genes are:	66


	strain: 00238-432	species: P.kernoviae
	the number of SigP gene is:	1314
	the number of SigP-RxLR genes are:	100
	the number of SigP-RxLR-EER genes are:	61


	strain: MPF4	species: P.lateralis
	the number of SigP gene is:	1394
	the number of SigP-RxLR genes are:	102
	the number of SigP-RxLR-EER genes are:	59


	strain: 310	species: P.parisitica
	the number of SigP gene is:	1570
	the number of SigP-RxLR genes are:	130
	the number of SigP-RxLR-EER genes are:	71


	strain: chvinca01	species: P.parisitica
	the number of SigP gene is:	1495
	the number of SigP-RxLR genes are:	119
	the number of SigP-RxLR-EER genes are:	64


	strain: cj01a1	species: P.parisitica
	the number of SigP gene is:	1652
	the number of SigP-RxLR genes are:	144
	the number of SigP-RxLR-EER genes are:	85


	strain: cj02b3	species: P.parisitica
	the number of SigP gene is:	1499
	the number of SigP-RxLR genes are:	110
	the number of SigP-RxLR-EER genes are:	60


	strain: cj05e6	species: P.parisitica
	the number of SigP gene is:	1503
	the number of SigP-RxLR genes are:	114
	the number of SigP-RxLR-EER genes are:	63


	strain: iac_01_95	species: P.parisitica
	the number of SigP gene is:	1499
	the number of SigP-RxLR genes are:	119
	the number of SigP-RxLR-EER genes are:	69


	strain: p10297	species: P.parisitica
	the number of SigP gene is:	1620
	the number of SigP-RxLR genes are:	150
	the number of SigP-RxLR-EER genes are:	93


	strain: p1569	species: P.parisitica
	the number of SigP gene is:	1640
	the number of SigP-RxLR genes are:	140
	the number of SigP-RxLR-EER genes are:	80


	strain: p1976	species: P.parisitica
	the number of SigP gene is:	1615
	the number of SigP-RxLR genes are:	138
	the number of SigP-RxLR-EER genes are:	83


	strain: 164328	species: P.ramorum
	the number of SigP gene is:	2244
	the number of SigP-RxLR genes are:	221
	the number of SigP-RxLR-EER genes are:	147


	strain: 67593	species: P.sojae
	the number of SigP gene is:	2691
	the number of SigP-RxLR genes are:	278
	the number of SigP-RxLR-EER genes are:	152
```

###Domain searching

Hmm models for the WY domain contained in many RxLRs were used to search gene
models predicted with Augustus. These were run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls gene_pred/augustus_unmasked/P.*/*/*_augustus_preds.aa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/hmmer/WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_aug_WY_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_aug_WY_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```
Results were as follows:

```
	P.cactorum 10300
	Initial search space (Z):              17801  [actual number of targets]
	Domain search space  (domZ):              93  [number of targets reported over threshold]
	P.cactorum 404
	Initial search space (Z):              16704  [actual number of targets]
	Domain search space  (domZ):              80  [number of targets reported over threshold]
	P.cactorum 414
	Initial search space (Z):              17141  [actual number of targets]
	Domain search space  (domZ):              83  [number of targets reported over threshold]
	P.fragariae 309-62
	Initial search space (Z):              18903  [actual number of targets]
	Domain search space  (domZ):             103  [number of targets reported over threshold]
	P.fragariae JHVZ02
	Initial search space (Z):              26596  [actual number of targets]
	Domain search space  (domZ):             133  [number of targets reported over threshold]
	P.ideai 371
	Initial search space (Z):              15901  [actual number of targets]
	Domain search space  (domZ):              67  [number of targets reported over threshold]
	P.infestans T30-4
	Initial search space (Z):              18326  [actual number of targets]
	Domain search space  (domZ):             121  [number of targets reported over threshold]
	P.infestans T30-4_unmasked
	Initial search space (Z):              38832  [actual number of targets]
	Domain search space  (domZ):             124  [number of targets reported over threshold]
	P.kernoviae 00238-432
	Initial search space (Z):              12862  [actual number of targets]
	Domain search space  (domZ):              53  [number of targets reported over threshold]
	P.lateralis MPF4
	Initial search space (Z):              12923  [actual number of targets]
	Domain search space  (domZ):              68  [number of targets reported over threshold]
	P.parisitica 310
	Initial search space (Z):              14054  [actual number of targets]
	Domain search space  (domZ):             118  [number of targets reported over threshold]
	P.parisitica chvinca01
	Initial search space (Z):              13589  [actual number of targets]
	Domain search space  (domZ):             101  [number of targets reported over threshold]
	P.parisitica cj01a1
	Initial search space (Z):              14771  [actual number of targets]
	Domain search space  (domZ):             113  [number of targets reported over threshold]
	P.parisitica cj02b3
	Initial search space (Z):              13590  [actual number of targets]
	Domain search space  (domZ):              97  [number of targets reported over threshold]
	P.parisitica cj05e6
	Initial search space (Z):              13500  [actual number of targets]
	Domain search space  (domZ):              97  [number of targets reported over threshold]
	P.parisitica iac_01_95
	Initial search space (Z):              13551  [actual number of targets]
	Domain search space  (domZ):              98  [number of targets reported over threshold]
	P.parisitica p10297
	Initial search space (Z):              14779  [actual number of targets]
	Domain search space  (domZ):             127  [number of targets reported over threshold]
	P.parisitica p1569
	Initial search space (Z):              14934  [actual number of targets]
	Domain search space  (domZ):             119  [number of targets reported over threshold]
	P.parisitica p1976
	Initial search space (Z):              14833  [actual number of targets]
	Domain search space  (domZ):             121  [number of targets reported over threshold]
	P.ramorum 164328
	Initial search space (Z):              18326  [actual number of targets]
	Domain search space  (domZ):             194  [number of targets reported over threshold]
	P.sojae 67593
	Initial search space (Z):              25649  [actual number of targets]
	Domain search space  (domZ):             220  [number of targets reported over threshold]
```

### From Augustus gene models - Hmm evidence of RxLR effectors

```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
	for Proteome in $(ls gene_pred/augustus_unmasked/P.*/*/*_augustus_preds.aa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_Aug_RxLR_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"__Aug_RxLR_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

Results were as follows:

Note: The number of genes used in the P. infestans genome shown here was noted
to indicate that this gene prediction was performed on the published masked
and unmasked gene models.

```
	P.cactorum 10300
	Initial search space (Z):              17801  [actual number of targets]
	Domain search space  (domZ):              70  [number of targets reported over threshold]
	P.cactorum 404
	Initial search space (Z):              16704  [actual number of targets]
	Domain search space  (domZ):              64  [number of targets reported over threshold]
	P.cactorum 414
	Initial search space (Z):              17141  [actual number of targets]
	Domain search space  (domZ):              66  [number of targets reported over threshold]
	P.fragariae 309-62
	Initial search space (Z):              18903  [actual number of targets]
	Domain search space  (domZ):              91  [number of targets reported over threshold]
	P.fragariae JHVZ02
	Initial search space (Z):              26596  [actual number of targets]
	Domain search space  (domZ):             120  [number of targets reported over threshold]
	P.ideai 371
	Initial search space (Z):              15901  [actual number of targets]
	Domain search space  (domZ):              44  [number of targets reported over threshold]
	P.infestans T30-4
	Initial search space (Z):              18326  [actual number of targets]
	Domain search space  (domZ):              70  [number of targets reported over threshold]
	P.infestans T30-4_unmasked
	Initial search space (Z):              38832  [actual number of targets]
	Domain search space  (domZ):              77  [number of targets reported over threshold]
	P.kernoviae 00238-432
	Initial search space (Z):              12862  [actual number of targets]
	Domain search space  (domZ):              71  [number of targets reported over threshold]
	P.lateralis MPF4
	Initial search space (Z):              12923  [actual number of targets]
	Domain search space  (domZ):              68  [number of targets reported over threshold]
	P.parisitica 310
	Initial search space (Z):              14054  [actual number of targets]
	Domain search space  (domZ):              80  [number of targets reported over threshold]
	P.parisitica chvinca01
	Initial search space (Z):              13589  [actual number of targets]
	Domain search space  (domZ):              69  [number of targets reported over threshold]
	P.parisitica cj01a1
	Initial search space (Z):              14771  [actual number of targets]
	Domain search space  (domZ):              99  [number of targets reported over threshold]
	P.parisitica cj02b3
	Initial search space (Z):              13590  [actual number of targets]
	Domain search space  (domZ):              67  [number of targets reported over threshold]
	P.parisitica cj05e6
	Initial search space (Z):              13500  [actual number of targets]
	Domain search space  (domZ):              71  [number of targets reported over threshold]
	P.parisitica iac_01_95
	Initial search space (Z):              13551  [actual number of targets]
	Domain search space  (domZ):              79  [number of targets reported over threshold]
	P.parisitica p10297
	Initial search space (Z):              14779  [actual number of targets]
	Domain search space  (domZ):             108  [number of targets reported over threshold]
	P.parisitica p1569
	Initial search space (Z):              14934  [actual number of targets]
	Domain search space  (domZ):              92  [number of targets reported over threshold]
	P.parisitica p1976
	Initial search space (Z):              14833  [actual number of targets]
	Domain search space  (domZ):              97  [number of targets reported over threshold]
	P.ramorum 164328
	Initial search space (Z):              18326  [actual number of targets]
	Domain search space  (domZ):             158  [number of targets reported over threshold]
	P.sojae 67593
	Initial search space (Z):              25649  [actual number of targets]
	Domain search space  (domZ):             173  [number of targets reported over threshold]
```


##Crinkler Prediction

###Motif identification

Crinkler motifs in Augustus gene models were identified by searching for
presence of two Crinkler motifs. This was performed using the following
commands:

```bash
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


```bash
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
		cat $OutDir/$HmmResults | sed '/inclusion threshold/q' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | sed '1,/inclusion threshold/d' | sed '/Domain annotation for each sequence/q' | tail -n +2 | head -n -3 | wc -l		HmmFasta="$Strain"_aug_CRN_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

## On Repeatmasked Genomes

The number of genes predicted in T-30 was considerably more than in publihsed gene models.
To determine if this was a result of repetative sequences affecting the gene models, gene
prediction was performed on the repeat-masked dataset for P. infestans.

Firstly, the genes predicted using the unmasked genome were moved to a new directory.
```bash
	mv gene_pred/augustus/P.infestans/T30-4 gene_pred/augustus/P.infestans/T30-4_unmasked
```

The repeatmasked genome was unziped and parsed using the following commands:

```bash
	cd /home/groups/harrisonlab/project_files/idris
	gunzip assembly/external_group/P.infestans/T30-4/dna/Phytophthora_infestans.ASM14294v1.26.dna_rm.genome.fa.gz
```
Gene prediction was then performed using Augustus
(using a P.cactorum gene model):
```bash
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

```bash
 mv gene_pred/sigP/P.infestans/T30-4 gene_pred/sigP/P.infestans/T30-4_unmasked
```

The commands to perform SigP prediction are shown below:

```bash
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
```bash
	mv analysis/sigP_rxlr/P.infestans/T30-4 analysis/sigP_rxlr/P.infestans/T30-4_unmasked
	mv analysis/hmmer/WY/P.infestans/T30-4 analysis/hmmer/WY/P.infestans/T30-4_unmasked
```
RxLR prediction was performed using the commands:

####Motif searching

RxLRs were predicted using the program rxlr_finder.py. This program uses a regular
expression to identify proteins that contain an RxLR motif within the amino acids
between the signal peptide cleavage site and 100aa downstream:

```bash
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


The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identfy RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Pathz in $(ls gene_pred/sigP_aug/P.*/*/*_aug_sp.aa); do
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors;
Strain=$(echo $Pathz | cut -d '/' -f4);
Organism=$(echo $Pathz | cut -d '/' -f3) ;
OutDir=analysis/rxlr_atg/"$Organism"/"$Strain";
mkdir -p $OutDir;
printf "\nstrain: $Strain\tspecies: $Organism\n";
printf "the number of SigP gene is:\t";
cat $Pathz | grep '>' | wc -l;
printf "the number of SigP-RxLR genes are:\t";
$ProgDir/RxLR_EER_regex_finder.py $Pathz > $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa;
cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_Aug_RxLR_regex.txt
cat $OutDir/"$Strain"_Aug_RxLR_regex.txt | wc -l
printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa | grep '>' | grep 'EER_motif_start' |  cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt
cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt | wc -l
printf "\n"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
Col2=RxLR_EER_regex_finder.py
GeneNames=$OutDir/"$Strain"_Aug_RxLR_regex.txt
GeneModels=analysis/rxlr_atg/"$Organism"/"$Strain"/"$Strain"_Aug.gff
$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 Name > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
done
```


####Domain searching

Hmm models for the WY domain contained in many RxLRs were used to search gene
models predicted with Augustus. These were run with the following commands:

```bash
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
		cat $OutDir/$HmmResults | sed '/inclusion threshold/q' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | sed '1,/inclusion threshold/d' | sed '/Domain annotation for each sequence/q' | tail -n +2 | head -n -3 | wc -l		HmmFasta="$Strain"_aug_WY_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

###Crinkler prediction (repeat masked genome)

Firstly, the CRN predictions that were made from unmasked genomes were moved
```bash
	mv analysis/CRN/P.infestans/T30-4 analysis/CRN/P.infestans/T30-4_unmasked
	mv analysis/hmmer/CRN/P.infestans/T30-4 analysis/hmmer/CRN/P.infestans/T30-4_unmasked
```

Crinkler prediction was performed using the commands:

####Motif identification

Crinkler motifs in masked P. infestans Augustus gene models were identified by searching for
presence of two Crinkler motifs. This was performed using the following
commands:

```bash
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


```bash
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
		cat $OutDir/$HmmResults | sed '/inclusion threshold/q' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | sed '1,/inclusion threshold/d' | sed '/Domain annotation for each sequence/q' | tail -n +2 | head -n -3 | wc -l
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

```bash
	for File in $(ls assembly/external_group/P.*/*/dna/*.genome.fa); do
		OutFile=$(echo $File | sed 's/.fa/.parsed.fa/g');
		echo $OutFile; cat $File | sed 's/ /_/g' | sed 's/|/_/g' > $OutFile;
	done
```

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing
Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen
	for Genome in $(ls assembly/external_group/*/*/dna/*.genome.parsed.fa); do
		echo $Genome
		qsub $ProgDir/path_pipe.sh $Genome
	done
```

This pipeline was also run on repeatmasked genomes:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen
	for Genome in $(ls assembly/external_group/*/*/dna/*.dna_rm.*.parsed.fa); do
	echo $Genome
	qsub $ProgDir/path_pipe.sh $Genome
	done
```

##RxLR Prediction

The number of ORF fragments, ORFs containing SigPs and ORFs containing
SigP & RxLR motifs were identified.

```bash
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

```bash
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


The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identfy RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
	for Pathz in $(ls analysis/rxlr_atg_unmasked/P.*/*/*.sp.pve); do
		ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors;
		Strain=$(echo $Pathz | cut -d '/' -f4);
		Organism=$(echo $Pathz | cut -d '/' -f3) ;
		OutDir=analysis/rxlr_atg_unmasked/"$Organism"/"$Strain";
		mkdir -p $OutDir;
		printf "\nstrain: $Strain\tspecies: $Organism\n";
		printf "the number of SigP gene is:\t";
		cat $Pathz | grep '>' | wc -l;
		printf "the number of SigP-RxLR genes are:\t";
		$ProgDir/RxLR_EER_regex_finder.py $Pathz > $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa;
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_ORF_RxLR_regex.txt
		cat $OutDir/"$Strain"_ORF_RxLR_regex.txt | wc -l
		printf "the number of SigP-RxLR-EER genes are:\t";
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa | grep '>' | grep 'EER_motif_start' |  cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_ORF_RxLR_EER_regex.txt
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.txt | wc -l
		printf "\n"
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		Col2=RxLR_EER_regex_finder.py
		GeneNames=$OutDir/"$Strain"_ORF_RxLR_regex.txt
		GeneModels=analysis/rxlr_atg_unmasked/"$Organism"/"$Strain"/"$Strain"_ORF.gff3
		$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 Name > $OutDir/"$Strain"_ORF_RxLR_regex.gff3
	done
```

###Domain searching


```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
		HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
		for Proteome in $(ls analysis/rxlr_atg_unmasked/P.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/hmmer/WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_WY_hmmer_out.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

### From Augustus gene models - Hmm evidence of RxLR effectors

```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
	for Proteome in $(ls analysis/rxlr_atg_unmasked/P.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_RxLR_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

### F) From ORF fragments - Hmm evidence of RxLR effectors
```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
	for Proteome in $(ls analysis/rxlr_atg_unmasked/P.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_RxLR_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```
Results were as follows:
```
	P.cactorum 10300
	Initial search space (Z):             456656  [actual number of targets]
	Domain search space  (domZ):             366  [number of targets reported over threshold]
	P.cactorum 404
	Initial search space (Z):             432691  [actual number of targets]
	Domain search space  (domZ):             321  [number of targets reported over threshold]
	P.cactorum 414
	Initial search space (Z):             445395  [actual number of targets]
	Domain search space  (domZ):             326  [number of targets reported over threshold]
	P.fragariae 309-62
	Initial search space (Z):             571132  [actual number of targets]
	Domain search space  (domZ):             479  [number of targets reported over threshold]
	P.fragariae SCRP245
	Initial search space (Z):             259266  [actual number of targets]
	Domain search space  (domZ):              32  [number of targets reported over threshold]
	P.ideai 371
	Initial search space (Z):             416004  [actual number of targets]
	Domain search space  (domZ):             336  [number of targets reported over threshold]
	P.infestans T30-4_unmasked
	Initial search space (Z):            1415148  [actual number of targets]
	Domain search space  (domZ):             859  [number of targets reported over threshold]
	P.kernoviae 00238-432
	Initial search space (Z):             406661  [actual number of targets]
	Domain search space  (domZ):             369  [number of targets reported over threshold]
	P.lateralis MPF4
	Initial search space (Z):             386192  [actual number of targets]
	Domain search space  (domZ):             330  [number of targets reported over threshold]
	P.parisitica 310
	Initial search space (Z):             442892  [actual number of targets]
	Domain search space  (domZ):             606  [number of targets reported over threshold]
	P.parisitica chvinca01
	Initial search space (Z):             377541  [actual number of targets]
	Domain search space  (domZ):             565  [number of targets reported over threshold]
	P.parisitica cj01a1
	Initial search space (Z):             426915  [actual number of targets]
	Domain search space  (domZ):             648  [number of targets reported over threshold]
	P.parisitica cj02b3
	Initial search space (Z):             378215  [actual number of targets]
	Domain search space  (domZ):             529  [number of targets reported over threshold]
	P.parisitica cj05e6
	Initial search space (Z):             373670  [actual number of targets]
	Domain search space  (domZ):             548  [number of targets reported over threshold]
	P.parisitica iac_01_95
	Initial search space (Z):             376075  [actual number of targets]
	Domain search space  (domZ):             556  [number of targets reported over threshold]
	P.parisitica p10297
	Initial search space (Z):             429443  [actual number of targets]
	Domain search space  (domZ):             654  [number of targets reported over threshold]
	P.parisitica p1569
	Initial search space (Z):             436049  [actual number of targets]
	Domain search space  (domZ):             636  [number of targets reported over threshold]
	P.parisitica p1976
	Initial search space (Z):             431194  [actual number of targets]
	Domain search space  (domZ):             650  [number of targets reported over threshold]
	P.ramorum 164328
	Initial search space (Z):             521713  [actual number of targets]
	Domain search space  (domZ):             601  [number of targets reported over threshold]
	P.sojae 67593
	Initial search space (Z):             689510  [actual number of targets]
	Domain search space  (domZ):             738  [number of targets reported over threshold]
```

##Crinkler Prediction

###Motif identification

```bash
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


```bash
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
		cat $OutDir/$HmmResults | sed '/inclusion threshold/q' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | sed '1,/inclusion threshold/d' | sed '/Domain annotation for each sequence/q' | tail -n +2 | head -n -3 | wc -l
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
```bash
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
```bash
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

```bash
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


```bash
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
		cat $OutDir/$HmmResults | sed '/inclusion threshold/q' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | sed '1,/inclusion threshold/d' | sed '/Domain annotation for each sequence/q' | tail -n +2 | head -n -3 | wc -l
		HmmFasta="$Strain"_Published_WY_hmmer_out.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```


##Crinklers

###Motif identification

```bash
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


```bash
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
		cat $OutDir/$HmmResults | sed '/inclusion threshold/q' | tail -n +16 | head -n -1 | wc -l
		cat $OutDir/$HmmResults | sed '1,/inclusion threshold/d' | sed '/Domain annotation for each sequence/q' | tail -n +2 | head -n -3 | wc -l
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

```bash
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

```bash
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

```bash
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

```bash
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

```bash
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

```bash
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

```bash
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

```bash
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

```bash
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

```bash
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
```bash
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

```bash
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
```bash
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

```bash
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
```bash
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

```bash
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
```bash
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

#Combining gff files

##Adding notes to effector gff files
Before combining gff files from multiple sources, the program used to predict
the gff features were named as 'Notes' in the attributes column of the gff file.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Note=RxLR_motif
	for File in $(ls analysis/sigP_rxlr/P.*/*/*_aug_RxLR_finder.gff); do
		Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
		OutFile=analysis/sigP_rxlr/$Organism/$Strain/"$Strain"_aug_RxLR_finder_source.gff
		$ProgDir/add_note_to_gff.pl $File $Note > $OutFile
	done
```

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Note=WY_hmmer
	for File in $(ls analysis/hmmer/WY/P.*/*/*_aug_WY_hmmer.gff); do
		Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
		OutFile=analysis/hmmer/WY/$Organism/$Strain/"$Strain"_aug_WY_hmmer_source.gff
		$ProgDir/add_note_to_gff.pl $File $Note > $OutFile
	done
```

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
	Note=CRN_motif
	for File in $(ls analysis/CRN/P.*/*/*_aug_LxLFLAK_HVLVVVP.gff); do
		Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
		OutFile=analysis/CRN/$Organism/$Strain/"$Strain"_LxLFLAK_HVLVVVP_source.gff
		$ProgDir/add_note_to_gff.pl $File $Note > $OutFile
	done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
Note=CRN_hmm
for File in $(ls analysis/hmmer/CRN/P.*/*/*_aug_CRN_hmmer.gff); do
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
OutFile=analysis/CRN/$Organism/$Strain/"$Strain"_aug_CRN_hmmer_source.gff
$ProgDir/add_note_to_gff.pl $File $Note > $OutFile
done
```

# Comparing Augustus predictions and atg.pl predictions


## atg.pl RxLR overlap

As >3000 RxLR containing ORFs were predicted in P.infestans
using the atg.pl script, the number of overlapping features
were identified in the gff file.

First the gff features for RxLR containing ORFs were extracted
from the gene models predicted with atg.pl

```bash
	RxlrFasta=analysis/rxlr_atg/P.infestans/T30-4/T30-4_sp_rxlr.fa
	RxlrHeaders=analysis/rxlr_atg/P.infestans/T30-4/T30-4_sp_rxlr_headers.txt
	OrfGff=analysis/rxlr_atg/P.infestans/T30-4/T30-4_ORF.gff
	RxlrGff=analysis/rxlr_atg/P.infestans/T30-4/T30-4_sp_rxlr.gff
	cat $RxlrFasta | grep '>' | cut -f1 | sed 's/>//g' > $RxlrHeaders
	cat $OrfGff | grep -w -f $RxlrHeaders > $RxlrGff
```

## High confidence effectors

Gff files of predicted effectors from Augustus and ORF predictions were combined
using a set of tools written for gffutils. Gene models were annotated with
notes if the gene was a predicted effector. Gene models from ORF finder were
merged if the showed evidence of tiling. Gene models from ORF finder were
merged with Augustus gene models. Finally those features with notes identifying
them as effector candidates were extracted.

The commands to run this were:

```bash
	for Strain in $(ls -d analysis/atg_sp_rxlr_unmasked/P.*/* | rev | cut -f1 -d '/' | rev); do
		echo $Strain;
		# Orf_Gff=$(ls gene_pred/ORF_finder/*/$Strain/"$Strain"_ORF.gff);
		Orf_Gff=analysis/atg_sp_rxlr_unmasked/P.*/$Strain/"$Strain"_ORF.gff
		Aug_Gff=gene_pred/augustus_unmasked/P.*/*/"$Strain"_augustus_preds.gtf
		# Aug_Gff=$(ls gene_pred/augustus/*/"$Strain"/"$Strain"_augustus_preds.gtf | grep -v 'old');
		echo $Orf_Gff;
		echo $Aug_Gff;
		ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff;
		qsub $ProgDir/merge_phytophthora_effectors.sh $Orf_Gff $Aug_Gff;
	done
```  


# Benchmarking

RxLRs predicted in published studies were compared to those
predicted in house on published gene models.


This was first performed on P. infestans.

Files containing the headers of putative effectors were made using the following commands:

```bash
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

```bash
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

```bash
	Pinf_confirmed_RxLR=analysis/benchmarking/P.infestans/T30-4/rxlr/P.inf_confirmed_RxLR_headers.txt
	Pinf_remainder_RxLR=analysis/benchmarking/P.infestans/T30-4/rxlr/P.inf_remainder_RxLR_headers.txt
	cat "$Pinf_pub_RxLR" "$Pinf_pred_mixed" | sort | uniq -d > $Pinf_confirmed_RxLR
	cat assembly/external_group/P.infestans/T30-4/rxlr/P.inf_RxLR_parsed.csv | grep 'PITG' | grep -v -f "$Pinf_confirmed_RxLR" > $Pinf_remainder_RxLR
```

This revealed that the remaining 20 annotated RxLRs in the Pinfestans genome did not contain
an RxLR domain, or had an RxLR variant.



# Whisman Regex
###ORF
```bash
	for File in $(ls analysis/rxlr_atg_unmasked/*/*/*.sp.pve); do
	Strain=$(echo $File | rev | cut -f2 -d'/' | rev);
	printf "$Strain\n";
	OutDir=$(dirname $File);
	mkdir -p $OutDir
	/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/rxlr/whisson_rxlr_finder.py $File > $OutDir/"$Strain"_whisson_ORF_rxlr.fa;
	printf "Number of SigP RxLR ORFs: \t";
	cat $OutDir/"$Strain"_whisson_ORF_rxlr.fa | grep '>' | wc -l;
	printf "Number of SigP RxLR ORFs with EER motifs: \t";
	cat $OutDir/"$Strain"_whisson_ORF_rxlr.fa | grep '>' | grep 'EER_motif_start(' | wc -l;
	done
```
##Augustus
```bash
	for File in $(ls gene_pred/sigP_aug/P.*/*/*_aug_sp.aa); do
	Organism=$(echo $File | rev | cut -f3 -d'/' | rev);
	Strain=$(echo $File | rev | cut -f2 -d'/' | rev);
	printf "$Strain\n";
	OutDir=analysis/whisson/$Organism/$Strain;
	mkdir -p $OutDir
	/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/rxlr/whisson_rxlr_finder.py $File > $OutDir/"$Strain"_whisson_Aug_rxlr.fa;
	printf "Number of SigP RxLR ORFs: \t";
	cat $OutDir/"$Strain"_whisson_Aug_rxlr.fa | grep '>' | wc -l;
	printf "Number of SigP RxLR ORFs with EER motifs: \t";
	cat $OutDir/"$Strain"_whisson_Aug_rxlr.fa | grep '>' | grep 'EER_motif_start(' | wc -l;
	done
```

##Published models
```bash
	for File in $(ls assembly/external_group/P.*/*/pep/*.fa); do
	Organism=$(echo $File | rev | cut -f4 -d'/' | rev);
	Strain=$(echo $File | rev | cut -f3 -d'/' | rev);
	printf "$Strain\n";
	OutDir=analysis/whisson/$Organism/"$Strain"_published;
	mkdir -p $OutDir
	/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/rxlr/whisson_rxlr_finder.py $File > $OutDir/"$Strain"_whisson_Pub_rxlr.fa;
	printf "Number of SigP RxLR ORFs: \t";
	cat $OutDir/"$Strain"_whisson_Pub_rxlr.fa | grep '>' | wc -l;
	printf "Number of SigP RxLR ORFs with EER motifs: \t";
	cat $OutDir/"$Strain"_whisson_Pub_rxlr.fa | grep '>' | grep 'EER_motif_start(' | wc -l;
	done
```
