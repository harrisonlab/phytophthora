# Commands used to identify motifs in the promotors of effectors based on RNAseq

## Extract promotor regions:

A new script was written to extract promotor regions upstream of a gene:

```bash
Assembly=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked_wrapped.fa)
GeneGff=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.gff3)
OutDir=analysis/meme/promotor_regions
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/RNAseq/P414/motif_discovery
$ProgDir/extract_promotor.pl --fasta $Assembly --gff $GeneGff --prefix $OutDir --ranges 1:100 101:200 201:300 301:400 401:500
```





## motif enrichment in RxLR sequences

Extraction of RxLR promotor sequences:

```bash
OutDir=analysis/meme/RxLR
mkdir -p $OutDir

# extract RxLRs
AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
RxlrHeaders=$OutDir/P414_RxLR_headers.txt
cat $AnnotTab | cut -f1,14 | grep 'RxLR' | cut -f1 | sed "s/\.t.//g"> $RxlrHeaders

# Create fasta files of RxLR upstream regions
for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
  Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
  mkdir $OutDir/$Region
  RxlrPromotors=$OutDir/$Region/P414_RxLR_${Region}_promotors.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $RxlrHeaders > $RxlrPromotors
done
```


```bash
# Run MEME on RxLR promotor regions
# initally run on the webserver using:
# meme P414_RxLR_promotors.fa -dna -oc . -nostatus -time 18000 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 0
# mast meme.xml P414_RxLR_promotors.fa -oc . -nostatus

screen -a

qlogin -pe smp 4
cd /data/scratch/armita/idris
for Query in $(ls analysis/meme/RxLR/*/P414_RxLR_*_promotors.fa); do
  Region=$(basename ${Query%.fa} | sed 's/promotor_regions.upstream//g')
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/meme
  meme $Query -p 4 -dna -oc $OutDir/meme -nostatus -mod zoops -nmotifs 5 -minw 6 -maxw 20 -objfun classic -revcomp -markov_order 0
  ls $OutDir/meme/meme.txt
  mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
  mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
  mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

# Also initial run of DREME on the webserver:
# dreme -verbosity 1 -oc . -dna -p P414_RxLR_promotors.fa -t 18000 -e 0.05 -dfile description
# qlogin -pe smp 4
# cd /data/scratch/armita/idris
# OutDir=analysis/meme/RxLR
# RxlrPromotors=$(ls $OutDir/P414_RxLR_promotors.fa)
# mkdir -p $OutDir/dreme
# dreme -verbosity 1 -oc $OutDir/dreme -dna -p $RxlrPromotors -t 18000 -e 0.05

for Query in $(ls analysis/meme/RxLR/*/P414_RxLR_*_promotors.fa); do
  Region=$(basename ${Query%.fa} | sed 's/promotor_regions.upstream//g')
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/dreme
  dreme -verbosity 1 -oc $OutDir/dreme -dna -p $Query -e 0.05
  ls $OutDir/dreme/dreme.txt
  mast $OutDir/dreme/dreme.xml $Query -oc $OutDir/dreme -nostatus
  mv $OutDir/dreme/mast.txt $OutDir/dreme/${Region}_mast.txt
  mv $OutDir/dreme/mast.html $OutDir/dreme/${Region}_mast.html
done
```

## Motif enrichment in early expressed RxLRs

Extraction of early expressed RxLR promotor sequences:

```bash
Prefix=RxLR_early
OutDir=analysis/meme/$Prefix
mkdir -p $OutDir

AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
QueryHeaders=$OutDir/${Prefix}_headers.txt
cat $AnnotTab | grep 'early expressed' | cut -f1,14 | grep 'RxLR' | cut -f1 | sed "s/\.t.//g" > $QueryHeaders

for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
  Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
  mkdir $OutDir/$Region
  Promotors=$OutDir/$Region/${Prefix}_${Region}_promotors.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $QueryHeaders > $Promotors
done
```

Extraction of non-early expressed promotor sequences:

```bash
  Prefix=RxLR_early
  OutDir=analysis/meme/$Prefix
  mkdir -p $OutDir

  AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
  ControlHeaders=$OutDir/${Prefix}_control_headers.txt
  cat $AnnotTab | grep -v 'early expressed' | cut -f1,14 | grep 'RxLR' | cut -f1 | sed "s/\.t.//g" > $ControlHeaders

  for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
    Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
    mkdir $OutDir/$Region
    Promotors=$OutDir/$Region/${Prefix}_control_${Region}_promotors.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $ControlHeaders > $Promotors
  done
```

Run MEME to compare early expressed to the control group

```bash

screen -a

qlogin -pe smp 4
cd /data/scratch/armita/idris
Prefix=RxLR_early
OutDir=analysis/meme/$Prefix
for Query in $(ls analysis/meme/$Prefix/*/${Prefix}*_promotors.fa | grep -v 'control'); do
  Control=$(echo $Query | sed "s/${Prefix}_/${Prefix}_control_/g")
  Region=$(basename ${Query%.fa} | sed "s/${Prefix}_//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/meme
  meme $Query -p 4 -dna -oc $OutDir/meme -neg $Control -nostatus -mod zoops -nmotifs 5 -minw 6 -maxw 20 -objfun de -revcomp -markov_order 0
  ls $OutDir/meme/meme.txt
  mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
  mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
  mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

for Query in $(ls analysis/meme/$Prefix/*/${Prefix}*_promotors.fa | grep -v 'control'); do
  Control=$(echo $Query | sed "s/${Prefix}_/${Prefix}_control_/g")
  Region=$(basename ${Query%.fa} | sed "s/${Prefix}_//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/dreme
  dreme -verbosity 1 -oc $OutDir/dreme -dna -p $Query -n $Control -e 0.05
  mast $OutDir/dreme/dreme.xml $Query -oc $OutDir/dreme -nostatus
  mv $OutDir/dreme/mast.txt $OutDir/dreme/${Region}_mast.txt
  mv $OutDir/dreme/mast.html $OutDir/dreme/${Region}_mast.html
done

ls analysis/meme/$Prefix/*/dreme/dreme.html
ls analysis/meme/$Prefix/*/dreme/*_mast.html
# ^ no file or directory means no dreme hits
```


## Motif enrichment in early expressed DEGs

Extraction of early expressed DEG promotor sequences:

```bash
Prefix=DEG_early
OutDir=analysis/meme/$Prefix
mkdir -p $OutDir

AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
QueryHeaders=$OutDir/${Prefix}_headers.txt
cat $AnnotTab | grep 'DEG' | grep 'early expressed' | cut -f1 | sed "s/\.t.//g" > $QueryHeaders

for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
  Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
  mkdir $OutDir/$Region
  Promotors=$OutDir/$Region/${Prefix}_${Region}_promotors.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $QueryHeaders > $Promotors
done
```

Extraction of non-early expressed promotor sequences:

```bash
  Prefix=DEG_early
  OutDir=analysis/meme/$Prefix
  mkdir -p $OutDir

  AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
  ControlHeaders=$OutDir/${Prefix}_control_headers.txt
  cat $AnnotTab | grep 'DEG' | grep -v 'early expressed' | cut -f1 | cut -f1 | sed "s/\.t.//g" > $ControlHeaders

  for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
    Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
    mkdir $OutDir/$Region
    Promotors=$OutDir/$Region/${Prefix}_control_${Region}_promotors.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $ControlHeaders > $Promotors
  done
```

Run MEME to compare early expressed to the control group

```bash

screen -a

qlogin -pe smp 4
cd /data/scratch/armita/idris
Prefix=DEG_early
OutDir=analysis/meme/$Prefix
for Query in $(ls analysis/meme/$Prefix/*/${Prefix}*_promotors.fa | grep -v 'control'); do
  Control=$(echo $Query | sed "s/${Prefix}_/${Prefix}_control_/g")
  Region=$(basename ${Query%.fa} | sed "s/${Prefix}_//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/meme
  meme $Query -p 4 -dna -oc $OutDir/meme -neg $Control -nostatus -mod zoops -nmotifs 5 -minw 6 -maxw 20 -objfun de -revcomp -markov_order 0
  ls $OutDir/meme/meme.txt
  mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
  mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
  mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

for Query in $(ls analysis/meme/$Prefix/*/${Prefix}*_promotors.fa | grep -v 'control'); do
  Control=$(echo $Query | sed "s/${Prefix}_/${Prefix}_control_/g")
  Region=$(basename ${Query%.fa} | sed "s/${Prefix}_//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/dreme
  dreme -verbosity 1 -oc $OutDir/dreme -dna -p $Query -n $Control -e 0.05
  mast $OutDir/dreme/dreme.xml $Query -oc $OutDir/dreme -nostatus
  mv $OutDir/dreme/mast.txt $OutDir/dreme/${Region}_mast.txt
  mv $OutDir/dreme/mast.html $OutDir/dreme/${Region}_mast.html
done

ls analysis/meme/$Prefix/*/dreme/dreme.html
ls analysis/meme/$Prefix/*/dreme/*_mast.html
# ^ no file or directory means no dreme hits
```


## Motif enrichment in late expressed DEGs

Extraction of late expressed DEG promotor sequences:


```bash
Prefix=DEG_late
OutDir=analysis/meme/$Prefix
mkdir -p $OutDir

AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
QueryHeaders=$OutDir/${Prefix}_headers.txt
cat $AnnotTab | grep 'DEG' | grep 'late expressed' | cut -f1 | sed "s/\.t.//g" > $QueryHeaders

for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
  Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
  mkdir $OutDir/$Region
  Promotors=$OutDir/$Region/${Prefix}_${Region}_promotors.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $QueryHeaders > $Promotors
done
```

Extraction of non-early expressed promotor sequences:

```bash
  Prefix=DEG_late
  OutDir=analysis/meme/$Prefix
  mkdir -p $OutDir

  AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
  ControlHeaders=$OutDir/${Prefix}_control_headers.txt
  cat $AnnotTab | grep 'DEG' | grep -v 'late expressed' | cut -f1 | cut -f1 | sed "s/\.t.//g" > $ControlHeaders

  for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
    Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
    mkdir $OutDir/$Region
    Promotors=$OutDir/$Region/${Prefix}_control_${Region}_promotors.fa
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $ControlHeaders > $Promotors
  done
```

Run MEME to compare early expressed to the control group

```bash

screen -a

qlogin -pe smp 4
cd /data/scratch/armita/idris
Prefix=DEG_late
OutDir=analysis/meme/$Prefix
for Query in $(ls analysis/meme/$Prefix/*/${Prefix}*_promotors.fa | grep -v 'control'); do
  Control=$(echo $Query | sed "s/${Prefix}_/${Prefix}_control_/g")
  Region=$(basename ${Query%.fa} | sed "s/${Prefix}_//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/meme
  meme $Query -p 4 -dna -oc $OutDir/meme -neg $Control -nostatus -mod zoops -nmotifs 5 -minw 6 -maxw 20 -objfun de -revcomp -markov_order 0
  ls $OutDir/meme/meme.txt
  mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
  mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
  mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

for Query in $(ls analysis/meme/$Prefix/*/${Prefix}*_promotors.fa | grep -v 'control'); do
  Control=$(echo $Query | sed "s/${Prefix}_/${Prefix}_control_/g")
  Region=$(basename ${Query%.fa} | sed "s/${Prefix}_//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/dreme
  dreme -verbosity 1 -oc $OutDir/dreme -dna -p $Query -n $Control -e 0.05
  mast $OutDir/dreme/dreme.xml $Query -oc $OutDir/dreme -nostatus
  mv $OutDir/dreme/mast.txt $OutDir/dreme/${Region}_mast.txt
  mv $OutDir/dreme/mast.html $OutDir/dreme/${Region}_mast.html
done

ls analysis/meme/$Prefix/*/dreme/dreme.html
ls analysis/meme/$Prefix/*/dreme/*_mast.html
# ^ no file or directory means no dreme hits
```

## Early vs Late

Extraction of early expressed DEG promotor sequences:

```bash
Prefix=DEG_early_vs_late
OutDir=analysis/meme/$Prefix
mkdir -p $OutDir

AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
QueryHeaders=$OutDir/${Prefix}_headers.txt
cat $AnnotTab | grep 'DEG' | grep 'early expressed' | cut -f1 | sed "s/\.t.//g" > $QueryHeaders

for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
  Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
  mkdir $OutDir/$Region
  Promotors=$OutDir/$Region/${Prefix}_${Region}_promotors.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $QueryHeaders > $Promotors
done
```

Extraction of late expressed DEG promotor sequences:

```bash
Prefix=DEG_early_vs_late
OutDir=analysis/meme/$Prefix
mkdir -p $OutDir

AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
QueryHeaders=$OutDir/${Prefix}_headers.txt
cat $AnnotTab | grep 'DEG' | grep 'late expressed' | cut -f1 | sed "s/\.t.//g" > $QueryHeaders

for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
  Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
  mkdir $OutDir/$Region
  Promotors=$OutDir/$Region/${Prefix}_control_${Region}_promotors.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $QueryHeaders > $Promotors
done
```

Run MEME to compare early expressed to the late expressed control group

```bash
screen -a

qlogin -pe smp 4
cd /data/scratch/armita/idris
Prefix=DEG_early_vs_late
OutDir=analysis/meme/$Prefix
for Query in $(ls analysis/meme/$Prefix/*/${Prefix}*_promotors.fa | grep -v 'control'); do
  Control=$(echo $Query | sed "s/${Prefix}_/${Prefix}_control_/g")
  Region=$(basename ${Query%.fa} | sed "s/${Prefix}_//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/meme
  meme $Query -p 4 -dna -oc $OutDir/meme -neg $Control -nostatus -mod zoops -nmotifs 5 -minw 6 -maxw 20 -objfun de -revcomp -markov_order 0
  ls $OutDir/meme/meme.txt
  mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
  mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
  mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

for Query in $(ls analysis/meme/$Prefix/*/${Prefix}*_promotors.fa | grep -v 'control'); do
  Control=$(echo $Query | sed "s/${Prefix}_/${Prefix}_control_/g")
  Region=$(basename ${Query%.fa} | sed "s/${Prefix}_//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/dreme
  dreme -verbosity 1 -oc $OutDir/dreme -dna -p $Query -n $Control -e 0.05
  mast $OutDir/dreme/dreme.xml $Query -oc $OutDir/dreme -nostatus
  mv $OutDir/dreme/mast.txt $OutDir/dreme/${Region}_mast.txt
  mv $OutDir/dreme/mast.html $OutDir/dreme/${Region}_mast.html
done

ls analysis/meme/$Prefix/*/dreme/dreme.html
ls analysis/meme/$Prefix/*/dreme/*_mast.html
# ^ no file or directory means no dreme hits
```


## Late vs Early DEGs

Extraction of late expressed DEG promotor sequences:

```bash
Prefix=DEG_late_vs_early
OutDir=analysis/meme/$Prefix
mkdir -p $OutDir

AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
QueryHeaders=$OutDir/${Prefix}_headers.txt
cat $AnnotTab | grep 'DEG' | grep 'late expressed' | cut -f1 | sed "s/\.t.//g" > $QueryHeaders

for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
  Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
  mkdir $OutDir/$Region
  Promotors=$OutDir/$Region/${Prefix}_${Region}_promotors.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $QueryHeaders > $Promotors
done
```

Extraction of early expressed DEG promotor sequences:

```bash
Prefix=DEG_late_vs_early
OutDir=analysis/meme/$Prefix
mkdir -p $OutDir

AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
QueryHeaders=$OutDir/${Prefix}_headers.txt
cat $AnnotTab | grep 'DEG' | grep 'early expressed' | cut -f1 | sed "s/\.t.//g" > $QueryHeaders

for Upstream in $(ls analysis/meme/promotor_regions.upstream*.fasta); do
  Region=$(basename ${Upstream%.fasta} | sed 's/promotor_regions.upstream//g')
  mkdir $OutDir/$Region
  Promotors=$OutDir/$Region/${Prefix}_control_${Region}_promotors.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Upstream --headers $QueryHeaders > $Promotors
done
```


Run MEME to compare early expressed to the late expressed control group

```bash
screen -a

qlogin -pe smp 4
cd /data/scratch/armita/idris
Prefix=DEG_late_vs_early
OutDir=analysis/meme/$Prefix
for Query in $(ls analysis/meme/$Prefix/*/${Prefix}*_promotors.fa | grep -v 'control'); do
  Control=$(echo $Query | sed "s/${Prefix}_/${Prefix}_control_/g")
  Region=$(basename ${Query%.fa} | sed "s/${Prefix}_//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/meme
  meme $Query -p 4 -dna -oc $OutDir/meme -neg $Control -nostatus -mod zoops -nmotifs 5 -minw 6 -maxw 20 -objfun de -revcomp -markov_order 0
  ls $OutDir/meme/meme.txt
  mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus
  mv $OutDir/meme/mast.txt $OutDir/meme/${Region}_mast.txt
  mv $OutDir/meme/mast.html $OutDir/meme/${Region}_mast.html
done

for Query in $(ls analysis/meme/$Prefix/*/${Prefix}*_promotors.fa | grep -v 'control'); do
  Control=$(echo $Query | sed "s/${Prefix}_/${Prefix}_control_/g")
  Region=$(basename ${Query%.fa} | sed "s/${Prefix}_//g")
  echo $Region
  OutDir=$(dirname $Query)
  mkdir -p $OutDir/dreme
  dreme -verbosity 1 -oc $OutDir/dreme -dna -p $Query -n $Control -e 0.05
  mast $OutDir/dreme/dreme.xml $Query -oc $OutDir/dreme -nostatus
  mv $OutDir/dreme/mast.txt $OutDir/dreme/${Region}_mast.txt
  mv $OutDir/dreme/mast.html $OutDir/dreme/${Region}_mast.html
done

ls analysis/meme/$Prefix/*/dreme/dreme.html
ls analysis/meme/$Prefix/*/dreme/*_mast.html
# ^ no file or directory means no dreme hits
```

Early expressed RxLRs included the motif 'CAMWWKYTGCRCMAN' in
8 of the 9 Early expressed RxLR effectors.


## Motif enrichment in CRNs

Extraction of late expressed CRN promotor sequences:

showing numbers of constant or DEGs:
```bash
AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
cat $AnnotTab | grep -w 'CRN' | grep 'constant' | wc -l
cat $AnnotTab | grep -w 'CRN' | grep 'DEG' | wc -l
cat $AnnotTab | grep -w 'CRN' | grep 'downregulated' | wc -l
cat $AnnotTab | grep -w 'CRN' | grep 'upregulated' | wc -l
```
<!--
```bash
OutDir=analysis/meme/DEG_early
mkdir -p $OutDir

# extract RxLRs
AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
QueryHeaders=$OutDir/P414_early_RxLR_headers.txt
cat $AnnotTab | grep 'DEG' | grep 'early expressed' | cut -f1 | sed "s/\.t./_upstream3000/g" | grep -v 'g9215' > $QueryHeaders

# Create fasta files of RxLR upstream regions
UpstreamGenes=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.upstream3000.fasta)
QueryPromotors=$OutDir/P414_early_DEG_promotors.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $UpstreamGenes --headers $QueryHeaders > $QueryPromotors
```

Extraction of non-early expressed promotor sequences:

```bash
OutDir=analysis/meme/RxLR_early
mkdir -p $OutDir

# extract RxLRs
AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
ControlHeaders=$OutDir/P414_control_DEG_headers.txt
cat $AnnotTab | grep 'DEG' | grep -v 'early expressed' | cut -f1 | grep -v 'g2856.t1' | sed "s/\.t./_upstream3000/g" > $ControlHeaders

# Create fasta files of RxLR upstream regions
UpstreamGenes=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.upstream3000.fasta)
ControlPromotors=$OutDir/P414_control_DEG_promotors.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $UpstreamGenes --headers $ControlHeaders > $ControlPromotors
```

Run MEME to compare early expressed to the control group

```bash
screen -a

qlogin -pe smp 4
cd /data/scratch/armita/idris
OutDir=analysis/meme/DEG_early
Query=$(ls $OutDir/P414_early_DEG_promotors.fa)
Control=$(ls $OutDir/P414_control_DEG_promotors.fa)
mkdir -p $OutDir/meme
meme $Query -p 4 -dna -oc $OutDir/meme -nostatus -time 18000 -mod zoops -nmotifs 5 -minw 6 -maxw 50 -objfun de -revcomp -markov_order 0 -neg $Control

mast $OutDir/meme/meme.xml $Query -oc $OutDir/meme -nostatus

ls $PWD/$OutDir/meme/mast.html
# this was downloaded for viewing in a web browser
``` -->

## Search for known motifs vs RxLRs:


Extraction of RxLR promotor sequences:

```bash
OutDir=analysis/meme/RxLR_jaspar
mkdir -p $OutDir

# extract RxLRs
AnnotTab=$(ls gene_pred/annotation/P.cactorum/414/414_annotation_ncbi.tsv)
RxlrHeaders=$OutDir/P414_RxLR_headers.txt
cat $AnnotTab | cut -f1,14 | grep 'RxLR' | cut -f1 | sed "s/\.t./_upstream3000/g"> $RxlrHeaders

# Create fasta files of RxLR upstream regions
UpstreamGenes=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.upstream3000.fasta)
RxlrPromotors=$OutDir/P414_RxLR_promotors.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $UpstreamGenes --headers $RxlrHeaders > $RxlrPromotors
```

```bash
screen -a

qlogin -pe smp 4
cd /data/scratch/armita/idris
OutDir=analysis/meme/RxLR_jaspar
Promotors=$(ls /home/armita/prog/meme/databases/motif_databases/JASPAR/JASPAR2018_CORE_non-redundant.meme)
Query=$(ls $OutDir/P414_RxLR_promotors.fa)
mkdir -p $OutDir/mast

mast $Promotors $Query -oc $OutDir/mast -nostatus

```
