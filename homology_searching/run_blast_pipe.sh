grep '>' analysis/rxlr/P.cactorum/404/404_sp_rxlr.fa  | cut -f1 > tmp.fa
while read line; do 
grep -A1 -w "$line" analysis/rxlr/P.cactorum/404/404_nuc.fa >> analysis/rxlr/P.cactorum/404/404_sp_rxlr_nuc.fa; 
done<tmp.fa
grep -A1 '>' analysis/rxlr/P.cactorum/404/404_sp_rxlr_nuc.fa | cut -d '|' -f5 > tmp.fa
sed 's/\s/_/g' tmp.fa > tmp2.fa
sed 's/>NODE/>P.cactorum_404_NODE/g' tmp2.fa > analysis/rxlr/P.cactorum/404/404_sp_rxlr_nuc_parsed.fa

grep '>' analysis/rxlr/P.cactorum/414/414_sp_rxlr.fa  | cut -f1 > tmp.fa
while read line; do 
grep -A1 -w "$line" analysis/rxlr/P.cactorum/414/414_nuc.fa >> analysis/rxlr/P.cactorum/414/414_sp_rxlr_nuc.fa; 
done<tmp.fa 
grep -A1 '>' analysis/rxlr/P.cactorum/414/414_sp_rxlr_nuc.fa | cut -d '|' -f5 > tmp.fa
sed 's/\s/_/g' tmp.fa > tmp2.fa
sed 's/>NODE/>P.cactorum_414_NODE/g' tmp2.fa > analysis/rxlr/P.cactorum/414/414_sp_rxlr_nuc_parsed.fa


grep '>' analysis/rxlr/P.cactorum/10300/10300_sp_rxlr.fa  | cut -f1 > tmp.fa
while read line; do 
grep -A1 -w "$line" analysis/rxlr/P.cactorum/10300/10300_nuc.fa >> analysis/rxlr/P.cactorum/10300/10300_sp_rxlr_nuc.fa; 
done<tmp.fa 
grep -A1 '>' analysis/rxlr/P.cactorum/10300/10300_sp_rxlr_nuc.fa | cut -d '|' -f5 > tmp.fa
sed 's/\s/_/g' tmp.fa > tmp2.fa
sed 's/>NODE/>P.cactorum_10300_NODE/g' tmp2.fa > analysis/rxlr/P.cactorum/10300/10300_sp_rxlr_nuc_parsed.fa



grep '>' analysis/rxlr/P.ideai/371/371_sp_rxlr.fa  | cut -f1 > tmp.fa
while read line; do 
grep -A1 -w "$line" analysis/rxlr/P.ideai/371/371_nuc.fa >> analysis/rxlr/P.ideai/371/371_sp_rxlr_nuc.fa; 
done<tmp.fa
grep -A1 '>' analysis/rxlr/P.ideai/371/371_sp_rxlr_nuc.fa | cut -d '|' -f5 > tmp.fa
sed 's/\s/_/g' tmp.fa > tmp2.fa
sed 's/>NODE/>P.ideai_371_NODE/g' tmp2.fa > analysis/rxlr/P.ideai/371/371_sp_rxlr_nuc_parsed.fa


grep '>' analysis/rxlr/P.fragariae/JHVZ02/JHVZ02_sp_rxlr.fa  | cut -f1 > tmp.fa
while read line; do 
grep -A1 -w "$line" analysis/rxlr/P.fragariae/JHVZ02/JHVZ02_nuc.fa >> analysis/rxlr/P.fragariae/JHVZ02/JHVZ02_sp_rxlr_nuc.fa; 
done<tmp.fa 
grep -A1 '>' analysis/rxlr/P.fragariae/JHVZ02/JHVZ02_sp_rxlr_nuc.fa | cut -d '|' -f5 > tmp.fa
sed 's/\s/_/g' tmp.fa > tmp2.fa
sed 's/_Phytophthora/>P.fragariae_JHVZ02/g' tmp2.fa > analysis/rxlr/P.fragariae/JHVZ02/JHVZ02_sp_rxlr_nuc_parsed.fa


cat analysis/rxlr/P.cactorum/404/404_sp_rxlr_nuc_parsed.fa analysis/rxlr/P.cactorum/414/414_sp_rxlr_nuc_parsed.fa analysis/rxlr/P.cactorum/10300/10300_sp_rxlr_nuc_parsed.fa analysis/rxlr/P.ideai/371/371_sp_rxlr_nuc_parsed.fa analysis/rxlr/P.fragariae/JHVZ02/JHVZ02_sp_rxlr_nuc_parsed.fa > analysis/blast_homology/rxlr_sp/appended_rxlr_sp_nuc.fa 
cat analysis/blast_homology/rxlr_sp/appended_rxlr_sp_nuc.fa | grep '>' | sort | uniq > tmp.txt
sed 's/,//g' tmp.fa > analysis/blast_homology/rxlr_sp/appended_rxlr_sp_nuc.fa

cat analysis/blast_homology/rxlr_sp/appended_rxlr_sp_nuc.fa | grep '>' | sort | uniq -d

qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/rxlr_sp/appended_rxlr_sp_nuc.fa dna assembly/velvet/P.cactorum/414/414_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/rxlr_sp/appended_rxlr_sp_nuc.fa dna assembly/velvet/P.cactorum/404/404_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/rxlr_sp/appended_rxlr_sp_nuc.fa dna assembly/velvet/P.cactorum/10300/10300_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/rxlr_sp/appended_rxlr_sp_nuc.fa dna assembly/velvet/P.ideai/371/371_assembly.41/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/rxlr_sp/appended_rxlr_sp_nuc.fa dna assembly/genbank/P.fragariae/JHVZ02/assembly_version2/JHVZ02.1.fsa_nt
