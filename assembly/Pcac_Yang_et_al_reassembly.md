
## Data extraction


## Data download

```bash
mkdir -p /home/groups/harrisonlab/project_files/idris
cd /home/groups/harrisonlab/project_files/idris
OutDir=raw_dna/pacbio/P.cactorum/D-1/downloaded
mkdir -p $OutDir
fastq-dump --split-files -A SRR3386345 --gzip --outdir $OutDir
```

Pacbio coverage was determined using:

```bash
for RawData in $(ls raw_dna/pacbio/P.cactorum/*/downloaded/*q.gz | grep 'D-1'); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
GenomeSz=65
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

86.44 coverage was obtained from pacbio sequencing
plus 76.49 illumina sequencing


## Assembly


### Canu assembly

### Read correction using Canu

```bash
Reads=$(ls raw_dna/pacbio/P.cactorum/*/downloaded/*q.gz | grep 'D-1')
GenomeSz="66m"
Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
OutDir=assembly/canu-1.6/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction.sh $Reads $GenomeSz $Strain $OutDir
```


### Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.6/*/*/*.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```


Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/P.cactorum/*/*_smartdenovo.dmo.lay.utg | grep 'D-1'); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

checking using busco

```bash
for Assembly in $(ls assembly/SMARTdenovo/P.cactorum/*/*_smartdenovo.dmo.lay.utg | grep 'D-1'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
OutDir=$(dirname $Assembly)
# OutDir=/data/scratch/armita/idris/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

Results were summarised using the commands:

```bash
  for File in $(ls assembly/SMARTdenovo/P.cactorum/*/*/short_summary_*.txt ); do
  Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
  Prefix=$(basename $File)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Single=$(cat $File | grep "(S)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Prefix\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
  done
```

```
P.cactorum	D-1	short_summary_D-1_smartdenovo.dmo.lay.txt	223	212	32	48	303
```

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	Genome=$(ls assembly/SMARTdenovo/P.cactorum/*/*_smartdenovo.dmo.lay.utg | grep 'D-1')
  OutDir=gene_pred/cegma/P.cactorum/D-1
  Prefix=D-1_smartdenovo
	qsub $ProgDir/sub_cegma.sh $Genome $Prefix $OutDir
```


# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  for Assembly in $(ls assembly/SMARTdenovo/P.cactorum/*/*_smartdenovo.dmo.lay.utg); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs_repmask
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```



The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.


```bash
for File in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked.fa | grep -w '414'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked.fa | grep -w '414'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```
Number of masked bases:
16513494
```
