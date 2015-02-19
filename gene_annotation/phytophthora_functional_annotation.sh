#!/bin/bash
#Functional analysis of phytophthora genomes

mkdir gene_pred/augustus/contaminants
mv gene_pred/augustus/P.cactorum/411 gene_pred/augustus/contaminants/.
mv gene_pred/augustus/P.fragariae/SCRP245 gene_pred/augustus/contaminants/.
SCRIPT_DIR=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan

for strain_path in $(ls -d gene_pred/augustus/P.*/*); do  
	echo $strain_path; 
	STRAIN=$(basename $strain_path); 
	echo $STRAIN
	$SCRIPT_DIR/sub_interproscan.sh $strain_path/"$STRAIN"_augustus_preds.aa
done
# append results once finished
SCRIPT_DIR=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for strain_path in $(ls -d gene_pred/interproscan/P.*/*); do
	STRAIN=$(basename $strain_path)
	ORGANISM=$(echo $strain_path | rev | cut -d "/" -f2 | rev)
	$SCRIPT_DIR/append_interpro.sh gene_pred/augustus/"$ORGANISM"/"$STRAIN"/"$STRAIN"_augustus_preds.aa gene_pred/interproscan/"$ORGANISM"/"$STRAIN"/raw
done

# Signal Peptides
SIGP_SCRIPTS=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
for FILEZ in $(ls -d gene_pred/interproscan/P.*/*/*_augustus_preds.aa_split_*); do   
echo $FILEZ;  
qsub $SIGP_SCRIPTS/pred_sigP.sh $FILEZ; 
done


for PATHZ in $(ls -d gene_pred/interproscan/P.*/*); do
STRAIN=$(echo $PATHZ | cut -d '/' -f4)
SPECIES=$(echo $PATHZ | cut -d '/' -f3)
echo $PATHZ
echo $STRAIN
IN_AA_STRING=''
IN_NEG_STRING=''
IN_TAB_STRING=''
IN_TXT_STRING=''
for GRP in $(ls -l gene_pred/interproscan/"$SPECIES"/"$STRAIN"/"$STRAIN"_augustus_preds* | cut -d '_' -f6 | sort -n); do  
echo $GRP
IN_AA_STRING="$IN_AA_STRING""gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_augustus_preds_split_$GRP""_sp.aa ";  
IN_NEG_STRING="$IN_NEG_STRING""gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_augustus_preds_split_$GRP""_sp_neg.aa ";  
IN_TAB_STRING="$IN_TAB_STRING""gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_augustus_preds_split_$GRP""_sp.tab "; 
IN_TXT_STRING="$IN_TXT_STRING""gene_pred/sigP/$SPECIES/$STRAIN/split/$STRAIN""_augustus_preds_split_$GRP""_sp.txt ";  
done
cat $IN_AA_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_sp.aa
cat $IN_NEG_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_neg_sp.aa
tail -n +2 -q $IN_TAB_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_sp.tabfile://localhost/.file/id=8321326.26/
cat $IN_TXT_STRING > gene_pred/sigP/$SPECIES/$STRAIN/"$STRAIN"_sp.txt
done







# Crinklers

# For predictions from nucleotide sequence data:

for strain_path in $(ls -d gene_pred/augustus/P.*/*); do  
	echo $strain_path; 
	STRAIN=$(basename $strain_path); 
	echo $STRAIN
	qsub /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/crinkler/sub_crinkler.sh $strain_path/"$STRAIN"_augustus_preds.codingseq LxLFLAK LYLAK HVLVVVP
done

MOTIF1=LxLFLAK
MOTIF2=LYLAK
MOTIF3=HVLVVVP

for infolder in $(ls -d analysis/crinkler/*/*); do
	STRAIN=$(echo $infolder | rev | cut -d "/" -f1 | rev)
	cat gene_pred/augustus/*/*/"$STRAIN"_augustus_preds.codingseq | grep -A1 -w "$(cat $infolder/findmotif_LxLFLAK.fa $infolder/findmotif_HVLVVVP.fa | grep '>' | cut -f1 | sed 's/ //g' | sort | uniq -d)" > $infolder/findmotif_LxLFLAK_HVLVVVP.fa
	cat gene_pred/augustus/*/*/"$STRAIN"_augustus_preds.codingseq | grep -A1 -w "$(cat $infolder/findmotif_LYLAK.fa $infolder/findmotif_HVLVVVP.fa | grep '>' | cut -f1 | sed 's/ //g' | sort | uniq -d)" > $infolder/findmotif_LYLAK_HVLVVVP.fa
done

# It was noticed that predictions were not particularly consistent between Phytophthora species. 
# Upon closer inspection grep results from amino acid sequence appeared to return more motifs.
# As such, sequences were linearised and motifs searched for.
# THis command linearises sequences:
# cat gene_pred/augustus/P.cactorum/404/404_augustus_preds.aa | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' 
for FILE in $(ls gene_pred/augustus/P.*/*/*_augustus_preds.aa); do
	echo "$FILE"
	echo "L.LFLAK"
	cat $FILE | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep 'L.LFLAK' | wc -l
	echo "HVLVVVP"
	cat $FILE | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep 'HVLVVVP' | wc -l
	echo "L.LFLAK & HVLVVVP"
	cat $FILE | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep 'L.LFLAK' | grep 'HVLVVVP' | wc -l
	echo "LYAK"
	cat $FILE | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep 'LYAK' | wc -l
	echo "LYAK & HVLVVVP"
	cat $FILE | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep 'LYAK' | grep 'HVLVVVP' | wc -l
done

# extract LxFLAK and HVLVVVP domain proteins from different isolates for alignment and analysis

for FILEZ in $(ls gene_pred/augustus/P.*/*/*_augustus_preds.aa); do
	SPECIES=$(printf $FILEZ | cut -f3 -d '/')
	STRAIN=$(printf $FILEZ | cut -f4 -d '/')
	mkdir -p analysis/motif_search/"$SPECIES"/"$STRAIN"
	OUTFILEZ=analysis/motif_search/"$SPECIES"/"$STRAIN"/"$STRAIN"_LxLFLAK_HVLVVVP.fa
	cat $FILEZ | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -B1 "L.LFLAK" | grep -B1 "HVLVVVP" | grep -v '\-\-' | sed "s/>/>$STRAIN\_/g" > $OUTFILEZ
done

# Use BLAST to identify if crinkler genes were present but were not predicted.
# First append all blast results
mkdir analysis/blast_homology/crinkler/blastx
cat analysis/motif_search/*/*/*_LxLFLAK_HVLVVVP.fa > analysis/blast_homology/crinkler/blastx/appended_LxLFLAK_HVLVVVP.fa
# Then perform BLASTx searches of proteins against each genome.
for GENOME in $(ls repeat_masked/P.*/*/*/*_contigs_unmasked.fa); do 
	echo $GENOME
	PATH_TO_BLAST=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	QUERY_FILE=analysis/blast_homology/crinkler/blastx
	qsub $PATH_TO_BLAST/blast_pipe.sh $QUERY_FILE/appended_LxLFLAK_HVLVVVP.fa protein $GENOME
done






#RxLR

# for strain_path in $(ls -d gene_pred/augustus/P.*/*); do  
# 	echo $strain_path; 
# 	STRAIN=$(basename $strain_path); 
# 	echo $STRAIN
# 	qsub /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/crinkler/sub_crinkler.sh $strain_path/"$STRAIN"_augustus_preds.codingseq RxLR
# done

for PATHZ in $(ls -d gene_pred/sigP/P.*/*); do 
echo $PATHZ; 
STRAIN=$(printf "$PATHZ" | cut -d '/' -f4); 
SPECIES=$(printf "$PATHZ" | cut -d '/' -f3); 
echo "strain: $STRAIN    species: $SPECIES";
mkdir -p analysis/sigP_rxlr/"$SPECIES"/"$STRAIN"
echo "the number of SigP gene is:"
cat gene_pred/sigP/"$SPECIES"/"$STRAIN"/"$STRAIN"_sp.aa | grep '>' | wc -l
echo "the number of SigP-RxLR genes are:"
cat gene_pred/sigP/"$SPECIES"/"$STRAIN"/"$STRAIN"_sp.aa | grep 'R.LR' | wc -l
cat gene_pred/sigP/"$SPECIES"/"$STRAIN"/"$STRAIN"_sp.aa | grep -B1 'R.LR' | grep -vw '\-\-' > analysis/sigP_rxlr/"$SPECIES"/"$STRAIN"/"$STRAIN"_sigP_RxLR.fa
done




# To ensure that all RxLR proteins are represented from each of the genomes the ATG.pl
# PexFinder script from the kamoun lab was run.
# This identified ORFs in the genome, used these to predict presence of SignalP
# and searched for the presence of an RxLR motif in these sequences.
for infile in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/path_pipe.sh $infile
done
# The number of SigP predicted for each genome was:
for FILEZ in $(ls analysis/rxlr_atg/P.*/*/*.sp.pve); do  
echo $FILEZ;  
cat $FILEZ | grep '>' | wc -l; 
done
# The number of SigP/RxLRs predicted for each genome was:
for FILEZ in $(ls analysis/rxlr_atg/P.*/*/*_sp_rxlr.fa); do 
echo $FILEZ; 
cat $FILEZ | grep '>' | wc -l; 
done
# The number of SigP/RxLRs predicted for each genome was:
for FILEZ in $(ls analysis/rxlr_atg/P.*/*/*_sp_rxlr_nuc.fa); do 
echo $FILEZ; 
cat $FILEZ | grep '>' | wc -l; 
done
# These produce two different results. 
# This is because the path_pipe script has been using grep without the -W function
# The _sp_rxlr.fa file is correct

for FILE in $(ls analysis/rxlr_atg/P.*/*/*_sp_rxlr.fa); do
	echo "$FILE"
	echo "RxLR"
	cat $FILE | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep 'R.LR' | wc -l
done
# get gff from sigP? > analysis/atg_sp_rxlr_unmasked/P.cactorum/10300/10300_sp_rxlr.gff
# some atg.pl outputs contain multiple amino acid sequences to a single fasta header (max 300aa long)
cat analysis/rxlr_atg/P.cactorum/10300/10300_sp_rxlr.fa | grep '>' | cut -f1 | sed 's/>//g' > names_tmp.txt
printf "" > 10300_sp_rxlr.gff
while read line; do
	grep -w "$line" analysis/rxlr_atg/P.cactorum/10300/10300_ORF.gff >> 10300_sp_rxlr.gff
done<names_tmp.txt
rm id_tmp.txt
# The number of genes predicted in 10300 with overlaps from the atg.pl script were identified.
bedtools intersect -c -a gene_pred/augustus/P.cactorum/10300/10300_augustus_preds.gtf -b 10300_sp_rxlr.gff > 10300_overlap.bed
cat 10300_overlap.bed | grep 'gene' | cut -f10 | grep -v -w '0' |  wc -l
# this identifeid that 348 predicted genes had overlaps with atg.pl ORF fragments
# with SignalP and RXLRs
bedtools intersect -c -a 10300_sp_rxlr.gff -b gene_pred/augustus/P.cactorum/10300/10300_augustus_preds.gtf > 10300_no_overlap.bed
# This accounted for 422 of the predicted ORF fragments.
cat 10300_no_overlap.bed | cut -f10 | grep -w '0' | wc -l
# This left 830 of the 1252 ORF fragments with SSigP RxLRs that did not relate to predicted genes.



# The path pipe .gff output is wrong. gff features must start at base 1.
bedtools intersect -c -a 10300_sp_rxlr.gff -b 10300_sp_rxlr.gff > 10300_self_overlap.bed



# use bowtie to align predicted RxLRs against the genome and use bedtools intersect to identify overlaps.
# for STRAINZ in $(ls -d analysis/rxlr_atg/P.*/*/*_sp_rxlr.fa); do 
# echo $STRAINZ; 
# STRAIN=$(echo $STRAINZ | rev | cut -f2 -d '/' | rev);
# THIS_DIR=$(dirname $STRAINZ) 
# grep '>' $THIS_DIR/"$STRAIN"_sp_rxlr.fa | cut -f1 > id_tmp.txt
# printf "" > $THIS_DIR/"$STRAIN"_sp_rxlr_nuc.fa
# while read line; do
# 	grep -w -A1 "$line" $THIS_DIR/"$STRAIN"_nuc.fa| sed 's/--//g' >> $THIS_DIR/"$STRAIN"_sp_rxlr_nuc.fa
# done<id_tmp.txt
# rm id_tmp.txt
# done
# 
# for STRAINZ in $(ls -d analysis/rxlr_atg/P.*/*/*_sp_rxlr.fa); do 
# echo $STRAINZ; 
# STRAIN=$(echo $STRAINZ | rev | cut -f2 -d '/' | rev);
# THIS_DIR=$(dirname $STRAINZ) 
# grep '>' $THIS_DIR/"$STRAIN"_sp_rxlr_nuc.fa | grep '>' | sed 's/>//g' > id_tmp.txt
# printf "" > $THIS_DIR/"$STRAIN"_sp_rxlr_nuc.fa
# while read line; do
# 	cat $THIS_DIR/"$STRAIN"_ORF.gff | grep -w "$line" >> $THIS_DIR/"$STRAIN"_sp_rxlr.gff
# done<id_tmp.txt
# rm id_tmp.txt
# done
# 
# PATHZ=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/annotation_by_alignment/
# qsub $PATHZ/bowtie_for_fasta_features.sh repeat_masked/P.cactorum/10300/version1_repmask/10300_contigs_softmasked.fa analysis/rxlr_atg/P.cactorum/10300/10300_sp_rxlr_nuc.fa sp_rxlr
# 

# To ensure all crinklers were predicted in the genome all ORFs were searched for crinkler motifs
for FILE in $(ls analysis/old/rxlr_atg_02-15/P.*/*/*.aa_cat.fa ); do
	echo "$FILE"
	echo "HVLVVVP LxFLAK"
	cat $FILE | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -A1 'HVLVVVP' | sed 's/--//g' | grep -A1 'L.LFLAK' | sed 's/--//g' | grep '>' | wc -l
done


cat analysis/rxlr_atg/P.cactorum/10300/10300_sp_rxlr_nuc.fa | grep '>' | sed 's/>//g' > tmp_gene_list.txt

