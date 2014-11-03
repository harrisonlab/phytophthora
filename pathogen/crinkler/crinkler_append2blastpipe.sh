#!/usr/bin/bash

#crinkler: append 2 blast

# grep -A1 '>' analysis/crinkler/P.cactorum/404/findmotif_LxLFLAK_HVLVVVP.fa | grep -v '\--' | sed 's/>/>P.cactorum_404_/' | less


for infile in $(ls analysis/crinkler/*/*/findmotif_LxLFLAK_HVLVVVP.fa); do
	ORGANISM=$(echo $infile | rev | cut -d "/" -f3 | rev)
	STRAIN=$(echo $infile | rev | cut -d "/" -f2 | rev)
	HEADER=">$ORGANISM""_$STRAIN""_"
	OUTNAME="$ORGANISM"_"$STRAIN""_LxLFLAK_HVLVVVP.fa"
	OUTFILE=$(echo $infile | sed "s/findmotif_LxLFLAK_HVLVVVP.fa/$OUTNAME/")
	grep -A1 '>' $infile | grep -v '\--' | sed "s/>/$HEADER/" | sed 's/Pcac10300_//' | sed 's/|/_/g' | sed 's/ Phytophthora fragariae strain 309.62 3runs_//g' | sed 's/, whole genome shotgun sequence//g' > $OUTFILE
done

cat analysis/crinkler/*/*/*_*_LxLFLAK_HVLVVVP.fa > analysis/blast_homology/crinkler/appended_crinklers.fa

for infile in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa); do
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/crinkler/appended_crinklers.fa dna $infile
done
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/crinkler/appended_crinklers.fa dna repeat_masked/P.ideai/371/371_assembly.41_repmask/371_contigs_softmasked.fa

/home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_differentials.pl analysis/blast_homology/P.cactorum/10300/10300_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.cactorum/404/404_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.cactorum/414/414_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.fragariae/JHVZ02/JHVZ02_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.ideai/371/371_appended_crinklers.fa_homologs.csv 

sort presence_tab.csv > presence_tab_sort.csv

mv *.csv analysis/blast_homology/crinkler/.


mkdir tblastx
cp analysis/blast_homology/P.cactorum/10300/10300_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.cactorum/404/404_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.cactorum/414/414_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.fragariae/JHVZ02/JHVZ02_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.ideai/371/371_appended_crinklers.fa_homologs.csv tblastx/.
mv tblastx analysis/blast_homology/crinkler/.


mkdir blastn
cp analysis/blast_homology/P.cactorum/10300/10300_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.cactorum/404/404_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.cactorum/414/414_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.fragariae/JHVZ02/JHVZ02_appended_crinklers.fa_homologs.csv analysis/blast_homology/P.ideai/371/371_appended_crinklers.fa_homologs.csv blastn/.
mv blastn analysis/blast_homology/crinkler/.
