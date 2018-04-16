
Analysis of DeSeq2 output


```bash
for UpFile in $(ls alignment/salmon/DeSeq2/*_up.txt); do
  DownFile=$(echo $UpFile | sed 's/_up.txt/_down.txt/g')
  DegFile=$(echo $UpFile | sed 's/_up.txt/_DEGs.txt/g')
  # echo $UpFile
  # echo $DownFile
  cat $UpFile $DownFile | grep -v 'baseMean' | cut -f1 | sort -u > $DegFile
  echo $DegFile
  cat $DegFile | wc -l
done

```

```bash
cd /data/scratch/armita/idris
DegFiles=$(ls alignment/salmon/DeSeq2/*_DEGs.txt | sed 's/txt/txt /' | tr -d '\n')
Conditions=$(ls alignment/salmon/DeSeq2/*_DEGs.txt | sed 's&alignment/salmon/DeSeq2/&&g' | sed 's&_DEGs.txt& &g' | tr -d '\n')
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/RNAseq/P414
$ProgDir/deseq2venn.py --input $DegFiles --conditions $Conditions > alignment/salmon/DeSeq2/DEG_table.tsv
$ProgDir/P414_5way_venn.r --inp alignment/salmon/DeSeq2/DEG_table.tsv --out alignment/salmon/DeSeq2/DEG_venn.pdf
```
raw_counts.txt
```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/RNAseq/P414
$ProgDir/counts_venn.py --input alignment/salmon/DeSeq2/raw_counts.txt > alignment/salmon/DeSeq2/fpkm_table_gt5.tsv
$ProgDir/P414_evidence_of_expression_venn.r --inp alignment/salmon/DeSeq2/fpkm_table_gt5.tsv --out alignment/salmon/DeSeq2/raw_counts_venn.pdf
```

fpkm_norm_counts
```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/RNAseq/P414
$ProgDir/counts_venn.py --input alignment/salmon/DeSeq2/fpkm_counts.txt > alignment/salmon/DeSeq2/fpkm_table_gt5.tsv
$ProgDir/P414_evidence_of_expression_venn.r --inp alignment/salmon/DeSeq2/fpkm_table_gt5.tsv --out alignment/salmon/DeSeq2/fpkm_venn.pdf
```

fpkm_norm_counts
```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/RNAseq/P414
$ProgDir/counts_venn.py --input alignment/salmon/DeSeq2/fpkm_norm_counts.txt > alignment/salmon/DeSeq2/fpkm_table_gt5.tsv
$ProgDir/P414_evidence_of_expression_venn.r --inp alignment/salmon/DeSeq2/fpkm_table_gt5.tsv --out alignment/salmon/DeSeq2/fpkm_norm_venn.pdf
```
