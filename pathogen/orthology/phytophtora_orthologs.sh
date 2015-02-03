# Phytophthora orthologous genes
ls -d gene_pred/augustus/P*/*
# gene_pred/augustus/P.cactorum/10300  gene_pred/augustus/P.cactorum/414       gene_pred/augustus/P.ideai/371
# gene_pred/augustus/P.cactorum/404    gene_pred/augustus/P.fragariae/JHVZ02
# gene_pred/augustus/P.cactorum/411    gene_pred/augustus/P.fragariae/SCRP245

set -- 10300 404 414 JHVZ02 371
	for a; do 
	shift
	for b; do
	printf "%s - %s\n" "$a" "$b"
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/sub_inparanoid.sh gene_pred/augustus/*/$a/"$a"_augustus_preds.aa gene_pred/augustus/*/$b/"$b"_augustus_preds.aa gene_pred/augustus/*/$a/"$a"_augustus_preds.gtf gene_pred/augustus/*/$b/"$b"_augustus_preds.gtf
	done 
done
