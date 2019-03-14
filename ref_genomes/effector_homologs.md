# Effector homologs

Identification of known effector gene homologs in reference genomes:


# download ref genomes


A tsv file was made containing Species, strain IDs and links to fasta files
for downloading genomes:

```bash
ls
mkdir assembly/genbank_03-19
# nano assembly/genbank_03-19/ref_genomes.tsv
Lines=$(cat assembly/genbank_03-19/ref_genomes.tsv | wc -l)
CurDir=/home/groups/harrisonlab/project_files/idris
for Num in $(seq 2 $Lines); do
  # echo $Num
  Species=$(cat assembly/genbank_03-19/ref_genomes.tsv | tail -n +$Num | head -n1 | cut -f1 | sed 's/ //g')
  Strain=$(cat assembly/genbank_03-19/ref_genomes.tsv | tail -n +$Num | head -n1 | cut -f2 | sed 's/ /_/g')
  Link=$(cat assembly/genbank_03-19/ref_genomes.tsv | tail -n +$Num | head -n1 | cut -f5)
  OutDir=assembly/genbank_03-19/$Species/$Strain
  mkdir -p $OutDir
  cd $OutDir
  wget $Link
  gunzip *.gz
  cd $CurDir
done
$(cat assembly/genbank_03-19/ref_genomes.tsv | cut -f4 | tail -n +2); do
  wget $Link
cat assembly/genbank_03-19/ref_genomes.tsv | cut -f1 | less
```

## Gene prediction - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	for Assembly in $(ls assembly/genbank_03-19/P.*/*/*.fsa_nt); do
    echo "$Assembly"
    OutDir=$(dirname $Assembly | sed 's/assembly/gene_pred/g')
  	qsub $ProgDir/run_ORF_finder.sh $Assembly $OutDir
  done
```



## Identifying homologs to previously characterised effectors

```bash
qlogin -pe smp=4
cd /home/groups/harrisonlab/project_files/idris
for Assembly in $(ls assembly/genbank_03-19/P.*/*/*.fsa_nt); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=analysis/blast_homology/genbank_03-19/$Organism/$Strain/Avr_homologs
  mkdir -p $OutDir
  AvrFasta=$(ls analysis/blast_homology/oomycete_avr_genes/appended_oomycete_effectors_cds.fasta)
  dbType="nucl"
  # CdsFasta=$(ls gene_pred/final_incl_ORF/$Organism/$Strain/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta)
  Eval="1e-10"
  #-------------------------------------------------------
  # 		Step 1.		Blast cds vs avr gene database
  #-------------------------------------------------------
  # Prefix="${Strain}_vs_appended_oomycete_effectors"
  # makeblastdb -in $AvrFasta -input_type fasta -dbtype nucl -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
  # tblastx -db $OutDir/$Prefix.db -query $CdsFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
  # cat $OutDir/${Prefix}_hits.txt | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt

  #-------------------------------------------------------
  # 		Step 2.		Blast avr gene database vs genome
  #-------------------------------------------------------
  Prefix="appended_oomycete_effectors_vs_${Strain}"
  makeblastdb -in $Assembly -input_type fasta -dbtype nucl -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
  tblastx -db $OutDir/$Prefix.db -query $AvrFasta -outfmt 6 -num_alignments 5 -out $OutDir/${Prefix}_hits.txt -evalue $Eval -num_threads 16
  cat $OutDir/${Prefix}_hits.txt | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
  #
  # #-------------------------------------------------------
  # # 		Step 3.		reciprocal hits
  # #-------------------------------------------------------
  # cat $OutDir/${Strain}_vs_appended_oomycete_effectors_hits_headers.txt > $OutDir/all_hits.txt
  # cat $OutDir/appended_oomycete_effectors_vs_${Strain}_hits_headers.txt | awk ' { t = $1; $1 = $2; $2 = t; print; } ' OFS=$'\t' >> $OutDir/all_hits.txt
  #
  # cat $OutDir/all_hits.txt | sort | uniq -d > $OutDir/reciprocal_best_hits.txt
  # printf "Number of reciprocal blast pairs: "
  # cat $OutDir/reciprocal_best_hits.txt | wc -l
  # rm $OutDir/all_hits.txt
done
```

```bash
rm analysis/blast_homology/genbank_03-19/Avr3a_hits_headers.tsv
rm analysis/blast_homology/genbank_03-19/Avr3b_hits_headers.tsv
for File in $(ls analysis/blast_homology/genbank_03-19/*/*/Avr_homologs/appended_oomycete_effectors_vs_*_hits.txt); do
  Organism=$(echo $File | rev |cut -f 4 -d '/' | rev)
  Strain=$(echo $File | rev |cut -f 3 -d '/' | rev)
  echo $Organism - $Strain
  cat $File | grep 'Avr3a' | wc -l
  cat $File | grep 'Avr3a' | sed "s/^/$Organism\t$Strain\t/g" >> analysis/blast_homology/genbank_03-19/Avr3a_hits_headers.tsv
  cat $File | grep 'Avr3b' | wc -l
  cat $File | grep 'Avr3b' | sed "s/^/$Organism\t$Strain\t/g" >> analysis/blast_homology/genbank_03-19/Avr3b_hits_headers.tsv
done
```

### Searches using reciprocal blasts vs ORFs


```bash
qlogin -pe smp 16 -l virtual_free=1G
cd /home/groups/harrisonlab/project_files/idris
for Assembly in $(ls assembly/genbank_03-19/P.*/*/*.fsa_nt | head -n1); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=analysis/blast_homology/genbank_03-19/$Organism/$Strain/Avr_homologs2
  mkdir -p $OutDir
  AvrFasta=$(ls analysis/blast_homology/oomycete_avr_genes/appended_oomycete_effectors_cds_extended.fasta)
  dbType="nucl"
  CdsFasta=$(ls gene_pred/genbank_03-19/$Organism/$Strain./${Organism}_nuc.fa)
  Eval="1e-10"
  #-------------------------------------------------------
  # 		Step 1.		Blast cds vs avr gene database
  #-------------------------------------------------------
  Prefix="${Strain}_vs_appended_oomycete_effectors"
  makeblastdb -in $AvrFasta -input_type fasta -dbtype nucl -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
  tblastx -db $OutDir/$Prefix.db -query $CdsFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval  -num_threads 16
  cat $OutDir/${Prefix}_hits.txt | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt

  #-------------------------------------------------------
  # 		Step 2.		Blast avr gene database vs genome
  #-------------------------------------------------------
  Prefix="appended_oomycete_effectors_vs_${Strain}"
  makeblastdb -in $CdsFasta -input_type fasta -dbtype nucl -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
  tblastx -db $OutDir/$Prefix.db -query $AvrFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval -num_threads 16
  cat $OutDir/${Prefix}_hits.txt | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt

  #-------------------------------------------------------
  # 		Step 3.		reciprocal hits
  #-------------------------------------------------------
  cat $OutDir/${Strain}_vs_appended_oomycete_effectors_hits_headers.txt > $OutDir/all_hits.txt
  cat $OutDir/appended_oomycete_effectors_vs_${Strain}_hits_headers.txt | awk ' { t = $1; $1 = $2; $2 = t; print; } ' OFS=$'\t' >> $OutDir/all_hits.txt

  cat $OutDir/all_hits.txt | sort | uniq -d > $OutDir/reciprocal_best_hits.txt
  printf "Number of reciprocal blast pairs: "
  cat $OutDir/reciprocal_best_hits.txt | wc -l
  rm $OutDir/all_hits.txt
done
```

### Using Blast pipe:


```bash
for Assembly in $(ls assembly/genbank_03-19/P.*/*/*.fsa_nt); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  cat $Assembly | cut -f1 -d ' ' > ${Assembly%.fsa_nt}.fa
  echo "$Organism - $Strain"
  OutDir=analysis/genbank_03-19/$Organism/$Strain/blast_pipe
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Query=$(ls analysis/blast_homology/oomycete_avr_genes/appended_oomycete_effectors_cds_extended.fasta)
	qsub $ProgDir/blast_pipe.sh $Query dna ${Assembly%.fsa_nt}.fa $OutDir
done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
  for BlastHits in $(ls analysis/genbank_03-19/*/*/blast_pipe/*.csv); do
    Organism=$(echo $BlastHits | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
    Assembly=$(ls assembly/genbank_03-19/$Organism/$Strain/*.fa)
    HitsGff=${BlastHits%.csv}.gff
    HitsFa=${BlastHits%.csv}.fa
    Column2=Avr_homolog
    NumHits=1
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
    $ProgDir/sequence_extractor.py --coordinates_file $HitsGff --header_column 1 --start_column 4 --stop_column 5 --strand_column 7 --id_column 9 --fasta_file $Assembly \
    | tr -d '"' | tr -d ";" | sed "s/ID=//g" > $HitsFa
    ExpGff=${BlastHits%.csv}_expanded.gff
    ExpFa=${BlastHits%.csv}_expanded.fa
    Column2=Avr_homolog
    NumHits=1
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder
    $ProgDir/gffexpander.pl +- 500 $HitsGff > $ExpGff
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/sequence_extractor.py --coordinates_file $ExpGff --header_column 1 --start_column 4 --stop_column 5 --strand_column 7 --id_column 9 --fasta_file $Assembly \
    | tr -d '"' | tr -d ";" | sed "s/ID=//g" > $ExpFa
  done
```


Extract data for blast hits to Avh32

```bash
for BlastHits in $(ls analysis/genbank_03-19/*/*/blast_pipe/*.fa | grep -v 'expanded'); do
  Organism=$(echo $BlastHits | rev | cut -f4 -d '/' | rev)  
  Strain=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
  cat $BlastHits | grep -A1 'Avh32' | sed "s/>/>${Organism}_${Strain}_/g"
done > analysis/genbank_03-19/Avh32_best_blast_hits.fa
for BlastHits in $(ls analysis/genbank_03-19/*/*/blast_pipe/*.fa | grep -v 'expanded'); do
  Organism=$(echo $BlastHits | rev | cut -f4 -d '/' | rev)  
  Strain=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
  cat $BlastHits | grep -A1 'Avr3a' | sed "s/>/>${Organism}_${Strain}_/g"
done > analysis/genbank_03-19/Avr3a_best_blast_hits.fa
```

Extract expanded sequences

```bash
for BlastHits in $(ls analysis/genbank_03-19/*/*/blast_pipe/*.fa | grep 'expanded'); do
  Organism=$(echo $BlastHits | rev | cut -f4 -d '/' | rev)  
  Strain=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
  cat $BlastHits | grep -A1 'Avh32' | sed "s/>/>${Organism}_${Strain}_/g"
done > analysis/genbank_03-19/Avh32_best_blast_hits_expanded.fa
for BlastHits in $(ls analysis/genbank_03-19/*/*/blast_pipe/*.fa | grep 'expanded'); do
  Organism=$(echo $BlastHits | rev | cut -f4 -d '/' | rev)  
  Strain=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
  cat $BlastHits | grep -A1 'Avr3a' | sed "s/>/>${Organism}_${Strain}_/g"
done > analysis/genbank_03-19/Avr3a_best_blast_hits_expanded.fa
```


### The same analysis was repeated on assembled genomes:



```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -e 'P.cactorum' -e 'P.idaei' | grep -v -e '414' -e '10300'); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  OutDir=analysis/genbank_03-19/$Organism/$Strain/blast_pipe2
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Query=$(ls analysis/blast_homology/oomycete_avr_genes/appended_oomycete_effectors_cds_extended.fasta)
	qsub $ProgDir/blast_pipe.sh $Query dna $Assembly $OutDir
done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
  for BlastHits in $(ls analysis/genbank_03-19/*/*/blast_pipe2/*.csv); do
    Organism=$(echo $BlastHits | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
    Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/*_contigs_unmasked.fa)
    HitsGff=${BlastHits%.csv}.gff
    HitsFa=${BlastHits%.csv}.fa
    Column2=Avr_homolog
    NumHits=1
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
    $ProgDir/sequence_extractor.py --coordinates_file $HitsGff --header_column 1 --start_column 4 --stop_column 5 --strand_column 7 --id_column 9 --fasta_file $Assembly \
    | tr -d '"' | tr -d ";" | sed "s/ID=//g" > $HitsFa
    ExpGff=${BlastHits%.csv}_expanded.gff
    ExpFa=${BlastHits%.csv}_expanded.fa
    Column2=Avr_homolog
    NumHits=1
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder
    $ProgDir/gffexpander.pl +- 500 $HitsGff > $ExpGff
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/sequence_extractor.py --coordinates_file $ExpGff --header_column 1 --start_column 4 --stop_column 5 --strand_column 7 --id_column 9 --fasta_file $Assembly \
    | tr -d '"' | tr -d ";" | sed "s/ID=//g" > $ExpFa
  done
```

Extract expanded sequences

```bash
for BlastHits in $(ls analysis/genbank_03-19/*/*/blast_pipe2/*.fa | grep 'expanded'); do
  Organism=$(echo $BlastHits | rev | cut -f4 -d '/' | rev)  
  Strain=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
  cat $BlastHits | grep -A1 'Avh32' | sed "s/>/>${Organism}_${Strain}_/g"
done >> analysis/genbank_03-19/Avh32_best_blast_hits_expanded.fa
for BlastHits in $(ls analysis/genbank_03-19/*/*/blast_pipe2/*.fa | grep 'expanded'); do
  Organism=$(echo $BlastHits | rev | cut -f4 -d '/' | rev)  
  Strain=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
  cat $BlastHits | grep -A1 'Avr3a' | sed "s/>/>${Organism}_${Strain}_/g"
done >> analysis/genbank_03-19/Avr3a_best_blast_hits_expanded.fa
```

<!--
```bash
for Proteome in $(ls gene_pred/braker/P.cactorum/10300/P.cactorum/final_genes_Braker.cds.fasta ); do
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Query=$(ls analysis/blast_homology/oomycete_avr_genes/appended_oomycete_effectors_cds.fasta)
	qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
done
``` -->
