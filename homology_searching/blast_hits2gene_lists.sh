# To run on output files from blast_pipe.sh using appended_oomycete_avr_cds.fasta 
# sequences as input.

# These are commands used on output using blastn searches 
# (requires editing blast_pipe.sh script)

mkdir -p analysis/blast_homology/oomycete_avr_genes/extracted_hits

for effector_name in $(cat analysis/blast_homology/P.cactorum/404/404_appended_oomycete_avr_cds.fasta_homologs.csv | cut -f1 | tail -n+2); do
	QUERY_SEQ=$(cat analysis/blast_homology/P.cactorum/404/404_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f9)
	printf ">$effector_name\n$QUERY_SEQ\n" > analysis/blast_homology/oomycete_avr_genes/extracted_hits/$effector_name.fa
done

for strain_path in $(ls -d analysis/blast_homology/P.*/*); do
	STRAIN=$(echo $strain_path | rev | cut -d '/' -f1 | rev)
	for effector_name in $(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | cut -f1 | tail -n+2); do 
		touch analysis/blast_homology/oomycete_avr_genes/extracted_hits/$effector_name.fa 
		HEADER=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 12)
		SEQ=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 20)
		START=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 18)
		STOP=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 19)
		if [ "$HEADER" != '' ]; then
			printf ">$STRAIN""_$HEADER""_$START""_$STOP""\n$SEQ\n" >> analysis/blast_homology/oomycete_avr_genes/extracted_hits/$effector_name.fa
		fi
	done
done

# These are commands used on output using tblastx searches 
# (requires editing blast_pipe.sh script)

mkdir -p analysis/blast_homology/oomycete_avr_genes/extracted_tblastx_hits

for effector_name in $(cat analysis/blast_homology/P.cactorum/404/404_appended_oomycete_avr_cds.fasta_homologs.csv | cut -f1 | tail -n+2); do
	QUERY_SEQ=$(cat analysis/blast_homology/P.cactorum/404/404_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f9)
	printf '' > analysis/blast_homology/oomycete_avr_genes/extracted_tblastx_hits/$effector_name.fa
done

for strain_path in $(ls -d analysis/blast_homology/P.*/*); do
	STRAIN=$(echo $strain_path | rev | cut -d '/' -f1 | rev)
	for effector_name in $(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | cut -f1 | tail -n+2); do 
		touch analysis/blast_homology/oomycete_avr_genes/extracted_tblastx_hits/$effector_name.fa 
		HEADER=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 12)
		SEQ=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 20)
		START=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 18)
		STOP=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 19)
		if [ "$HEADER" != '' ]; then
			printf ">$STRAIN""_$HEADER""_$START""_$STOP""\n$SEQ\n" >> analysis/blast_homology/oomycete_avr_genes/extracted_tblastx_hits/$effector_name.fa
		fi
	done
done

for strain_path in $(ls -d analysis/blast_homology/P.*/*); do
STRAIN=$(echo $strain_path | rev | cut -d '/' -f1 | rev)
for effector_name in $(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | cut -f1 | tail -n+2); do 
touch analysis/blast_homology/oomycete_avr_genes/extracted_tblastx_hits/$effector_name.fa 
HEADER=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 12)
SEQ=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 20)
START=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 18)
STOP=$(cat $strain_path/"$STRAIN"_appended_oomycete_avr_cds.fasta_homologs.csv | grep "$effector_name" | cut -f 19)
if [ "$HEADER" != '' ]; then
printf ">$STRAIN""_$HEADER""_$START""_$STOP""\n$SEQ\n" >> analysis/blast_homology/oomycete_avr_genes/extracted_tblastx_hits/$effector_name.fa
fi
done
done
