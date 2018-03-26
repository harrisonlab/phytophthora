
## Data extraction


for P.cactorum data:
```bash
  cd /home/groups/harrisonlab/project_files/idris
  RawDatDir=/home/harrir/projects/pacbio_test/p_cact
  mkdir -p raw_dna/pacbio/P.cactorum/414
  cp -r $RawDatDir/A07_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/B07_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/G06_1 raw_dna/pacbio/P.cactorum/414/.
  cp -r $RawDatDir/H06_1 raw_dna/pacbio/P.cactorum/414/.
  OutDir=raw_dna/pacbio/P.cactorum/414/extracted
  mkdir -p $OutDir
  cat raw_dna/pacbio/P.cactorum/414/*/Analysis_Results/*.subreads.fastq | gzip -cf > $OutDir/concatenated_pacbio.fastq.gz
  #
  RawDat=/home/groups/harrisonlab/raw_data/raw_seq/pacbio/Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage.tar.gz
  mkdir -p raw_dna/pacbio/P.cactorum/414
  cp -r $RawDat raw_dna/pacbio/P.cactorum/414/.
  cd raw_dna/pacbio/P.cactorum/414
  tar -zxvf Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage.tar.gz
  # Data for 414 was contained in E02_1, F02_1 and G02_1
  cat \
  Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage/E02_1/Analysis_Results/*.subreads.fastq \
  Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage/F02_1/Analysis_Results/*.subreads.fastq \
  Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage/G02_1/Analysis_Results/*.subreads.fastq \
  | gzip -cf > extracted/concatenated_pacbio_extra_coverage.fastq.gz
  rm Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage.tar.gz
  rm -r Richard_Harrison_NEMR.RH.ENQ-933.C.02_extra_coverage
```

Pacbio coverage was determined using:


Data quality was visualised once again following trimming:

```bash
for RawData in $(ls raw_dna/pacbio/P.cactorum/414/extracted/*q.gz); do
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

<!--
for P. fragariae data (commands for tom to run)

```bash
  cd /home/groups/harrisonlab/project_files/phytophthora_fragariae
  RawDatDir=/home/harrir/projects/pacbio_test/p_frag
  mkdir -p raw_dna/pacbio/P.fragariae/Bc16
  cp -r $RawDatDir/C07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  cp -r $RawDatDir/D07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  cp -r $RawDatDir/E07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  cp -r $RawDatDir/F07_1 raw_dna/pacbio/P.fragariae/Bc16/.
  OutDir=raw_dna/pacbio/P.fragariae/Bc16/extracted
  mkdir -p $OutDir
  cat raw_dna/pacbio/P.fragariae/Bc16/*/Analysis_Results/*.subreads.fastq > $OutDir/concatenated_pacbio.fastq
``` -->

<!-- Sequencing depth and possible contamination was identified using kmer counting
kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

```bash
  for Reads in $(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq); do
    echo $Reads
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/kmc_kmer_counting.sh $Reads
  done
``` -->

## Assembly


### Canu assembly


### Read correction using Canu

```bash
Run1=$(ls ../../../../home/groups/harrisonlab/project_files/idris/raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq.gz)
Run2=$(ls ../../../../home/groups/harrisonlab/project_files/idris/raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz)
Reads=../../../../home/groups/harrisonlab/project_files/idris/raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_both.fastq.gz
cat $Run1 $Run2 > $Reads
GenomeSz="72m"
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
for Assembly in $(ls assembly/SMARTdenovo/P.cactorum/414/414_smartdenovo.dmo.lay.utg); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

checking using busco

```bash
for Assembly in $(ls assembly/SMARTdenovo/P.cactorum/414/414_smartdenovo.dmo.lay.utg); do
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
  for File in $(ls data/scratch/armita/idris/busco/run_contigs_min_500bp_filtered_renamed/short_summary_contigs_min_500bp_filtered_renamed.txt ); do
  Strain='P414'
  Organism='P.cactorum'
  Prefix=$(basename $File)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Single=$(cat $File | grep "(S)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Prefix\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
  done
```


The falcon assembly was polished using Pilon

```bash
for Assembly in $(ls assembly/SMARTdenovo/P.cactorum/414/414_smartdenovo.dmo.lay.utg); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
IlluminaDir=$(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
TrimF2_Read=$(ls $IlluminaDir/F/414_run2_F_trim.fq.gz);
TrimR2_Read=$(ls $IlluminaDir/R/414_run2_R_trim.fq.gz);
TrimF3_Read=$(ls $IlluminaDir/F/414_170210_F_trim.fq.gz);
TrimR3_Read=$(ls $IlluminaDir/R/414_170210_R_trim.fq.gz);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
OutDir=$(dirname $Assembly)/polished
# OutDir=assembly/falcon/P.cactorum/414/v2
Iterations='5'
Ploidy='diploid'
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon_3_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations $Ploidy
done
```


Contigs were renamed in accordance with ncbi recomendations.

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/SMARTdenovo/P.cactorum/414/polished/pilon_5.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/SMARTdenovo/P.cactorum/414/polished/contigs_min_500bp_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

checking using busco

```bash
for Assembly in $(ls assembly/SMARTdenovo/P.cactorum/414/polished/contigs_min_500bp_renamed.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls assembly/SMARTdenovo/P.cactorum/414/polished/*/short_summary_*.txt ); do
  Strain=$(echo $File| rev | cut -d '/' -f5 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f6 | rev)
  Prefix=$(basename $File)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Prefix\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```


<!--
```bash
  Reads=$(ls raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq.gz)
  Run1=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq.gz)
  Run2=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz)
  Reads=raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_both.fastq.gz
  cat $Run1 $Run2 > $Reads
  GenomeSz="65m"
  Strain=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Reads | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir=assembly/canu/$Organism/"$Strain"_3
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu.sh $Reads $GenomeSz $Prefix $OutDir
```

```bash
Run1=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq.gz)
GenomeSz="75m"
Strain=$(echo $Run1 | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Run1 | rev | cut -f4 -d '/' | rev)
Prefix="$Strain"_canu
OutDir=assembly/canu/$Organism/"$Strain"_run1only
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/submit_canu.sh $Run1 $GenomeSz $Prefix $OutDir

  Run2=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz)
  GenomeSz="75m"
  Strain=$(echo $Run2 | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Run2 | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir=assembly/canu/$Organism/"$Strain"_run2only
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu.sh $Run2 $GenomeSz $Prefix $OutDir
```

-->

### Falcon Assembly

FALCON is an assembler designed by Pacific Biosciences to assemble long-read data. It is also 'diploid-aware'. This file contains an example set of commands for running FALCON on PacBio data for the P414 strain of Phytophthora cactorum. This is run on NIABs triticum computer.


```bash
mkdir -p assembly/falcon/P.cactorum/414

ssh adarmitage@10.1.10.170
/bin/bash
```

```bash
DataDir=/data/projects/armita
mkdir -p $DataDir
scp armita@149.155.34.72:/home/groups/harrisonlab/project_files/idris/raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq.gz $DataDir/.
scp armita@149.155.34.72:/home/groups/harrisonlab/project_files/idris/raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz $DataDir/.
gunzip $DataDir/*.fastq.gz
# scp -r armita@149.155.34.72:/home/groups/harrisonlab/project_files/idris/raw_dna/paired/P.cactorum/414 $DataDir/.

```

The following lines must be in your bash profile
```bash
  # source /home/sobczm/bin/FALCON-integrate/env.sh
  export PATH=/home/sobczm/bin/cmake-3.8.0/bin:${PATH}
  export PATH=/home/sobczm/bin/gawk-4.1.4:${PATH}
  export PYTHONPATH=/data/software/smrtanalysis/install/smrtanalysis_2.3.0.140936/analysis/bin
  export PYTHONPATH="$PYTHONPATH:/data/software/smrtanalysis/install/smrtanalysis_2.3.0.140936/common/lib"
  export PYTHONPATH="$PYTHONPATH:/data/software/smrtanalysis/install/smrtanalysis_2.3.0.140936/analysis/lib/python2.7"
  export PYTHONPATH="$PYTHONPATH:/home/sobczm/usr/local/lib/python2.7/site-packages"
  export PYTHONPATH="$PYTHONPATH:/home/sobczm/bin/FALCON-integrate/fc_env/lib/python2.7/site-packages"
  export PYTHONPATH="$PYTHONPATH:/data/software/smrtanalysis/install/smrtanalysis_2.3.0.140936/analysis/lib"
  export PYTHONUSERBASE=/home/sobczm/bin/FALCON-integrate/fc_env
  export PATH=$PYTHONUSERBASE/bin:${PATH}
  export PATH=/home/sobczm/usr/local/bin:${PATH}
  export PATH=/home/sobczm/bin/pbh5tools/bin:${PATH}
  export PATH=/data/software/smrtanalysis/install/smrtanalysis_2.3.0.140936/analysis/bin:${PATH}
```

In order to run, FALCON needs two files to be available. One that tells it where the fasta files of reads are and one that specifies parameters to run with.

```bash
# mkdir -p falcon
# cd falcon
# Run1=$(ls ../../../data/projects/armita/*.fastq | head -n1 | tail -n1)
# Run2=$(ls ../../../data/projects/armita/*.fastq | head -n2 | tail -n1)
# printf "$Run1\n$Run2\n" > input.fofn
cd /data/projects/armita
cat concatenated_pacbio.fastq | sed -n '1~4s/^@/>/p;2~4p' > concatenated_pacbio.fasta
cat concatenated_pacbio_extra_coverage.fastq | sed -n '1~4s/^@/>/p;2~4p' > concatenated_pacbio_extra_coverage.fasta
# Run1=$(ls ../../../data/projects/armita/*.fasta | head -n1 | tail -n1)
# Run2=$(ls ../../../data/projects/armita/*.fasta | head -n2 | tail -n1)
# printf "$Run1\n$Run2\n" > input.fofn
echo "concatenated_pacbio.fasta" > Pcac.fofn
echo "concatenated_pacbio_extra_coverage.fasta" >> Pcac.fofn

printf \
"[General]
use_tmpdir = True
job_type = local

# list of fasta files
input_fofn = Pcac.fofn

#input type, raw or pre-assembled reads (preads, error corrected reads)
input_type = raw

# The length cutoff used for seed reads used for initial mapping during error correction
#E. coli: automatic calculation
# length_cutoff = -1
#fungal
length_cutoff = 6000


# The length cutoff used for seed reads used for assembly overlapping of preads
# "-1" indicates FALCON should calculate the cutoff using
# the user-defined genome length and coverage cut off
# otherwise, user can specify length cut off in bp (e.g. 2000)
###In a general sense, longer pread length cut offs will increase the
###contiguity (contig N50) in your assembly, but may result in shorter over all assembly length.


#fungal
# length_cutoff_pr = 3500
length_cutoff_pr = 5000
genome_size = 66000000
#seed_coverage = 30


## resource usage ## EMPTY FOR LOCAL USAGE
# grid settings for...
jobqueue = production
# daligner step of raw reads
sge_option_da =
# las-merging of raw reads
sge_option_la =
# consensus calling for preads
sge_option_pda =
# daligner on preads
sge_option_pla =
# las-merging on preads
sge_option_fc =
# final overlap/assembly
sge_option_cns =


# job concurrency settings for...
# all jobs
default_concurrent_jobs = 32
# preassembly
da_concurrent_jobs = 32
la_concurrent_jobs = 32
# consensus calling of preads
cns_concurrent_jobs = 32
# overlap detection
pda_concurrent_jobs = 32
pla_concurrent_jobs = 32

# daligner parameter options for...
# https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide/

##initial overlap of raw reads
#Lambda
#S argument controls the size of jobs - bigger number will result in a smaller number
#of longer-running jobs.
#fungal
pa_HPCdaligner_option =  -v -B128 -t16 -e0.75 -M24 -l3200 -k18 -h480 -w8 -s100

## overlap of preads
#fungal
ovlp_HPCdaligner_option = -v -B128 -M24 -k24 -h1024 -e.96 -l2500 -s100

## parameters for creation of dazzler database of...
## https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide/
## raw reads
#fungal
pa_DBsplit_option = -a -x500 -s200

#fungal
ovlp_DBsplit_option = -s200

## settings for consensus calling for preads
#fungal
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --max_n_read 200 --n_core 8

#fungal
overlap_filtering_setting = --max_diff 120 --max_cov 120 --min_cov 2 --n_core 12
" \
> Pcac_fc_run.cfg
```

Run falcon job itself in a screen session
```bash
screen -a
/bin/bash
# source /home/sobczm/bin/FALCON-integrate/fc_env/bin/activate
# source /home/sobczm/bin/FALCON-integrate/env.sh
# export PYTHONUSERBASE=/data/software/FALCON-integrate/fc_env
# export PATH=$PYTHONUSERBASE/bin:$PATH

# FALCON_WORKSPACE=$(pwd)
# PYTHONUSERBASE=$(pwd)/fc_env
# FALCON_PREFIX=${PYTHONUSERBASE}
# mkdir -p ${FALCON_PREFIX}/include
# mkdir -p ${FALCON_PREFIX}/bin
# mkdir -p ${FALCON_PREFIX}/lib


# FALCON_WORKSPACE=$(pwd)
# PYTHONUSERBASE=$(pwd)/fc_env
# FALCON_PREFIX=${PYTHONUSERBASE}
# PATH=${PYTHONUSERBASE}/bin:${FALCON_PREFIX}/bin:${PATH}
# export PYTHONUSERBASE
# export FALCON_WORKSPACE
# export FALCON_PREFIX
# export PATH
# mkdir -p ${FALCON_PREFIX}/include
# mkdir -p ${FALCON_PREFIX}/bin
# mkdir -p ${FALCON_PREFIX}/lib

rm -rf 0-rawreads/ 1-preads_ovl/ 2-asm-falcon/ all.log mypwatcher/ scripts/ sge_log
rm fc_run.log pypeflow.log logging.ini
fc_run.py Pcac_fc_run.cfg
# fc_run.py fc_run.cfg
```

Assess assembly Graph and Pread Overlaps
https://pb-falcon.readthedocs.io/en/latest/tutorial.html#assembly-graph-and-pread-overlaps

Assembly contiguity can be enhanced by adjusting a few parameters in the last stage of the assembly process. You can try a grid of pread length cut offs for the filtering of the final overlaps in the assembly graph. In a general sense, longer pread length cut offs will increase the contiguity (contig N50) in your assembly, but may result in shorter over all assembly length. To try different length cut off, rename your 2-asm-falcon dir, modify the config file, rename the log and mypwatcher directory, and restart FALCON:

```bash
cd 2-asm-falcon
fc_ovlp_stats --fofn ../1-preads_ovl/merge-gather/las.fofn > ovlp.stats
scp -r ovlp.stats armita@149.155.34.72:/home/groups/harrisonlab/project_files/idris/assembly/falcon/P.cactorum/414/.
cd ../
```

```R
library(ggplot2)
ovlp<-read.table("ovlp.stats", header=F)
colnames(ovlp)<-c("pread","length","fivePrimeOvlps","threePrimeOvlps")

pdf(file="OvlpHist.pdf", width=11, height=5)
par(oma=c(3,3,2,0), cex=1.6, las=1, mar=c(4,4,2,2), mfrow=c(1,2))
hist(ovlp$fivePrimeOvlps, breaks=100,
        xlab="number of overlaps between preads",
        ylab="count", main="Five Prime")
hist(ovlp$threePrimeOvlps, breaks=100,
        xlab="number of overlaps between preads",
        ylab="count", main="Three Prime")
hist(ovlp$fivePrimeOvlps, breaks=1000,
        xlim=c(1,100),
        xlab="number of overlaps between preads",
        ylab="count", main="Five Prime")
hist(ovlp$threePrimeOvlps, breaks=1000,
        xlim=c(1,100),
        xlab="number of overlaps between preads",
        ylab="count", main="Three Prime")

dev.off()
```

```bash
cp Pcac_fc_run.cfg 2-asm-falcon/.
scp -r 2-asm-falcon armita@149.155.34.72:/home/groups/harrisonlab/project_files/idris/assembly/falcon/P.cactorum/414/v2
```


```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/falcon/P.cactorum/414/*/p_ctg.fa); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

checking using busco

```bash
for Assembly in $(ls assembly/falcon/P.cactorum/414/*/p_ctg.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

The falcon assembly was polished using Pilon

```bash
for Assembly in $(ls ../../../../home/groups/harrisonlab/project_files/idris/assembly/falcon/P.cactorum/414/v2/p_ctg.fa); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
IlluminaDir=$(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
TrimF2_Read=$(ls $IlluminaDir/F/414_run2_F_trim.fq.gz);
TrimR2_Read=$(ls $IlluminaDir/R/414_run2_R_trim.fq.gz);
TrimF3_Read=$(ls $IlluminaDir/F/414_170210_F_trim.fq.gz);
TrimR3_Read=$(ls $IlluminaDir/R/414_170210_R_trim.fq.gz);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
# OutDir=$(dirname $Assembly)/polished
OutDir=assembly/falcon/P.cactorum/414/v2
Iterations='5'
Ploidy='diploid'
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon_3_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations $Ploidy
done
```


Contigs were renamed in accordance with ncbi recomendations.

```bash
  touch tmp.csv
  for Assembly in $(ls assembly/falcon/P.cactorum/414/2-asm-falcon/polished/pilon_5.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/falcon/P.cactorum/414/2-asm-falcon/polished/contigs_min_500bp_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

checking using busco

```bash
for Assembly in $(ls assembly/falcon/P.cactorum/414/2-asm-falcon/polished/*.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls assembly/falcon/P.cactorum/414/2-asm-falcon/polished/*/short_summary_*.txt ); do
  Strain=$(echo $File| rev | cut -d '/' -f5 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f6 | rev)
  Prefix=$(basename $File)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Prefix\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

```
  P.cactorum	414	short_summary_pilon_1.txt	285	4	14	303
  P.cactorum	414	short_summary_pilon_2.txt	285	4	14	303
  P.cactorum	414	short_summary_pilon_3.txt	285	4	14	303
  P.cactorum	414	short_summary_pilon_4.txt	285	4	14	303
  P.cactorum	414	short_summary_pilon_5.txt	285	4	14	303
  P.cactorum	414	short_summary_contigs_min_500bp_renamed.txt	285	4	14	303
```


Contigs were identified that had blast hits to non-phytophthora genomes

```bash
for Assembly in $(ls assembly/falcon/P.cactorum/414/2-asm-falcon/polished/contigs_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
echo "$Organism - $Strain"
# Exclude_db="bact,virus,hsref"
Exclude_db="paenibacillus"
Good_db="phytoph"
AssemblyDir=$(dirname $Assembly)
OutDir=$AssemblyDir/../deconseq_Paen
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
done
```

Results were summarised using the commands:

```bash
  for File in $(ls assembly/falcon/P.*/*/*/deconseq_Paen/log.txt | grep '414'); do
    Name=$(echo $File | rev | cut -f4 -d '/' | rev);
    Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
    Both=$(cat $File |cut -f2 | head -n2 | tail -n1);
    Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
    printf "$Name\t$Good\t$Both\t$Bad\n";
  done
```

```
  414_v2	222	0	5
```


Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls ../../../../home/groups/harrisonlab/project_files/idris/assembly/falcon/P.*/*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta | grep -e '414'); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
# OutDir=$(dirname $Assembly)
OutDir=/data/scratch/armita/idris/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

checking using busco

```bash
for Assembly in $(ls ../../../../home/groups/harrisonlab/project_files/idris/assembly/falcon/P.*/*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta | grep -e '414'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
# OutDir=$(dirname $Assembly)
OutDir=/data/scratch/armita/idris/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

Results were summarised using the commands:

```bash
  for File in $(ls data/scratch/armita/idris/busco/run_contigs_min_500bp_filtered_renamed/short_summary_contigs_min_500bp_filtered_renamed.txt ); do
  Strain='P414'
  Organism='P.cactorum'
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
P.cactorum	P414	short_summary_contigs_min_500bp_filtered_renamed.txt	265 258	5	33	303
```

<!--
### Canu assembly

```bash
  Run1=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq.gz)
  Run2=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz)
  GenomeSz="75m"
  Strain=$(echo $Run1 | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Run1 | rev | cut -f4 -d '/' | rev)
  Prefix="$Strain"_canu
  OutDir=assembly/canu/$Organism/"$Strain"_modified_script
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu_2lib.sh $Run1 $Run2 $GenomeSz $Prefix $OutDir
```
-->
<!--
Data quality was visualised once again following trimming:

```bash
for RawData in $(ls assembly/canu/P.cactorum/414_modified_script/414_canu.trimmedReads.fasta.gz); do
echo $RawData;
GenomeSz=65
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
``` -->
<!--
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/canu/*/*/*.contigs.fasta | grep -v 'old' | grep -w '414'); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=$(dirname $Assembly)
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

Assemblies were polished using Pilon

```bash
  for Assembly in $(ls assembly/canu/*/*/*.contigs.fasta | grep -v 'old' | grep -w '414'); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    IlluminaDir=$(ls -d qc_dna/paired/$Organism/$Strain)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
    TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
    TrimF2_Read=$(ls $IlluminaDir/F/414_run2_F_trim.fq.gz);
    TrimR2_Read=$(ls $IlluminaDir/R/414_run2_R_trim.fq.gz);
    TrimF3_Read=$(ls $IlluminaDir/F/414_170210_F_trim.fq.gz);
    TrimR3_Read=$(ls $IlluminaDir/R/414_170210_R_trim.fq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    echo $TrimF3_Read
    echo $TrimR3_Read
    OutDir=assembly/canu/$Organism/$Strain/polished
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
    qsub $ProgDir/sub_pilon_3_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
  done
```

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/canu/*/*/polished/pilon.fasta); do
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    OutDir=assembly/canu/$Organism/$Strain/polished
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```
-->


### Spades Assembly

For P. cactorum

```bash
# for PacBioDat in $(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_extra_coverage.fastq.gz); do
for PacBioDat in $(ls ../../../../home/groups/harrisonlab/project_files/idris/raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio_both.fastq.gz); do
echo $StrainPath
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades/multiple_libraries
Organism=$(echo $PacBioDat | rev | cut -f4 -d '/' | rev)
Strain=$(echo $PacBioDat | rev | cut -f3 -d '/' | rev)
IlluminaDir=$(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
TrimF2_Read=$(ls $IlluminaDir/F/414_run2_F_trim.fq.gz);
TrimR2_Read=$(ls $IlluminaDir/R/414_run2_R_trim.fq.gz);
TrimF3_Read=$(ls $IlluminaDir/F/414_170210_F_trim.fq.gz);
TrimR3_Read=$(ls $IlluminaDir/R/414_170210_R_trim.fq.gz);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
# CuttOff="50"
CuttOff="20"
# OutDir=assembly/spades_pacbio/$Organism/$Strain
# OutDir=assembly/spades_pacbio/$Organism/"$Strain"_50x
OutDir=assembly/spades_pacbio/$Organism/"$Strain"_20x
qsub $ProgDir/subSpades_3lib_pacbio.sh $PacBioDat $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $CuttOff
done
```

Contigs shorter than 500bp were removed from the assembly

```bash
for Contigs in $(ls assembly/spades_pacbio/*/*/contigs.fasta); do
AssemblyDir=$(dirname $Contigs)
mkdir $AssemblyDir/filtered_contigs
FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
$FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs/contigs_min_500bp.fasta
done
```

Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades_pacbio/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=assembly/spades_pacbio/$Organism/$Strain/filtered_contigs
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


checking using busco

```bash
for Assembly in $(ls assembly/spades_pacbio/*/*/filtered_contigs/contigs_min_500bp.fasta); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
  echo "$Organism - $Strain"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # BuscoDB="Fungal"
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls assembly/falcon/P.cactorum/414/2-asm-falcon/polished/*/short_summary_*.txt ); do
  Strain=$(echo $File| rev | cut -d '/' -f5 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f6 | rev)
  Prefix=$(basename $File)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Prefix\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```


## Merging pacbio and hybrid assemblies

```bash
for PacBioAssembly in $(ls assembly/SMARTdenovo/P.cactorum/414/polished/contigs_min_500bp_renamed.fasta); do
# for PacBioAssembly in $(ls ../../../../home/groups/harrisonlab/project_files/idris/assembly/falcon/P.cactorum/414/2-asm-falcon/polished/contigs_min_500bp_renamed.fasta); do
Organism=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $PacBioAssembly | rev | cut -f3 -d '/' | rev)
# HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
HybridAssembly=$(ls assembly/spades_pacbio/${Organism}/${Strain}_20x/contigs.fasta)
OutDir=assembly/merged_SMARTdenovo_spades/$Organism/$Strain
AnchorLength=500000
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
done
```

Quast

```bash

for Assembly in $(ls assembly/spades_pacbio/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


This merged assembly was polished using Pilon

```bash
# for Assembly in $(ls -d assembly/merged_canu_spades/P.*/*/filtered_contigs/contigs_min_500bp_renamed.fasta | grep -w -e '414_v2'); do
for Assembly in $(ls assembly/merged_SMARTdenovo_spades/P.cactorum/414/merged.fasta); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
IlluminaDir=$(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_dna/paired/$Organism/$Strain)
echo $Strain
echo $Organism
TrimF1_Read=$(ls $IlluminaDir/F/414_run1_F_trim.fq.gz);
TrimR1_Read=$(ls $IlluminaDir/R/414_run1_R_trim.fq.gz);
TrimF2_Read=$(ls $IlluminaDir/F/414_run2_F_trim.fq.gz);
TrimR2_Read=$(ls $IlluminaDir/R/414_run2_R_trim.fq.gz);
TrimF3_Read=$(ls $IlluminaDir/F/414_170210_F_trim.fq.gz);
TrimR3_Read=$(ls $IlluminaDir/R/414_170210_R_trim.fq.gz);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
OutDir=$(dirname $Assembly)/polished
Iterations='5'
Ploidy='diploid'
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon_3_libs.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir $Iterations $Ploidy
done
```


# Preliminary analysis

## Checking MiSeq coverage against P414 contigs

```bash
# for Assembly in $(ls assembly/falcon/P.*/*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta | grep -e '414'); do
Jobs=$(qstat | grep 'sub_pilon_' | wc -l)
while [ $Jobs -gt 0 ]; do
sleep 2
printf "."
Jobs=$(qstat | grep 'sub_pilon_' | wc -l)
done
printf "\n"
touch tmp.csv
for Assembly in $(ls assembly/merged_SMARTdenovo_spades/P.cactorum/414/polished/pilon_5.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
done
rm tmp.csv
# for Assembly in $(ls ../../../../home/groups/harrisonlab/project_files/idris/assembly/falcon/P.cactorum/414/2-asm-falcon/polished/contigs_min_500bp_renamed.fasta); do
for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/polished/contigs_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
# IlluminaDir=$(ls -d qc_dna/paired/$Organism/414)
IlluminaDir=$(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_dna/paired/$Organism/414)
echo "$Organism - $Strain"
TrimF1_Read=$(ls $IlluminaDir/F/*trim.fq.gz | head -n1 | tail -n1);
TrimR1_Read=$(ls $IlluminaDir/R/*trim.fq.gz | head -n1 | tail -n1);
TrimF2_Read=$(ls $IlluminaDir/F/*trim.fq.gz | head -n2 | tail -n1);
TrimR2_Read=$(ls $IlluminaDir/R/*trim.fq.gz | head -n2 | tail -n1);
TrimF3_Read=$(ls $IlluminaDir/F/*trim.fq.gz | head -n3 | tail -n1);
TrimR3_Read=$(ls $IlluminaDir/R/*trim.fq.gz | head -n3 | tail -n1);
echo $TrimF1_Read
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
InDir=$(dirname $Assembly)
# OutDir=$InDir/aligned_MiSeq
OutDir=$InDir/aligned_MiSeq
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie_3lib.sh $Assembly $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir
done
```

Determine read depth over each bp.

```bash
qlogin -pe smp 8
# cd /home/groups/harrisonlab/project_files/idris
cd /data/scratch/armita/idris
# for Bam in $(ls assembly/falcon/P.cactorum/414/2-asm-falcon/deconseq_Paen/aligned_MiSeq/contigs_min_500bp_filtered_renamed.fasta_aligned.bam); do
for Bam in $(ls assembly/merged_SMARTdenovo_spades/P.cactorum/414/polished/aligned_MiSeq/contigs_min_500bp_renamed.fasta_aligned.bam); do
# Strain=$(echo $Bam | rev | cut -f5 -d '/' | rev)
# Organism=$(echo $Bam | rev | cut -f6 -d '/' | rev)
# echo "$Organism - $Strain"
OutDir=$(dirname $Bam)/read_depth
# Region="5:500000-1100000" # <chr:from-to>
# OutDir=/data/scratch/armita/idris/read_depth
mkdir -p $OutDir
samtools sort -@ 8 -o $OutDir/P414_illumina_vs_P414_assembly_sorted.bam $Bam
samtools depth -aa $OutDir/P414_illumina_vs_P414_assembly_sorted.bam > $OutDir/P414_illumina_vs_P414_assembly_depth.tsv


ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
$ProgDir/cov_by_window.py --cov $OutDir/P414_illumina_vs_P414_assembly_depth.tsv | sed "s/$/\tP414/g" | sed 's/contig_//g'> $OutDir/P414_illumina_vs_P414_assembly_depth_10kb.tsv
$ProgDir/coverage_by_contig.py --cov $OutDir/P414_illumina_vs_P414_assembly_depth.tsv | sort -k3 -n > $OutDir/coverage_by_conitg.tsv
cat $OutDir/coverage_by_conitg.tsv | head
done
```

```
  contig_1	3595319	0
  contig_117	141214	0
  contig_16	1005393	0
  contig_2	2482273	0
  contig_167	27024	33
  contig_45	590457	37
  contig_33	674714	50
  contig_129	100751	53
  contig_18	986313	53
  contig_139	57947	54
```

```R
library(readr)
appended_df <- read_delim("assembly/merged_SMARTdenovo_spades/P.cactorum/414/polished/aligned_MiSeq/read_depth/P414_illumina_vs_P414_assembly_depth_10kb.tsv",
     "\t", escape_double = FALSE, col_names = FALSE,
     trim_ws = TRUE)
colnames(appended_df) <- c("contig","position", "depth", "strain")
appended_df$depth_mod <- ifelse(appended_df$depth > 150, 150, appended_df$depth)



# install.packages("ggplot2")
library(ggplot2)
require(scales)



for (i in 1:198){
  paste(i)
  contig = paste(i, sep = "")
  p0 <- ggplot(data=appended_df[appended_df$contig == contig, ], aes(x=`position`, y=`depth_mod`, group=1)) +
    geom_line() +
    labs(x = "Position", y = "Coverage") +
    scale_y_continuous(breaks=seq(0,150,25), limits=c(0,150))
  outfile = paste("contig", i, "cov.tiff", sep = "_")
  ggsave(outfile , plot = p0, device = 'tiff', path = '/data/scratch/armita/idris/assembly/merged_SMARTdenovo_spades/P.cactorum/414/polished/aligned_MiSeq/read_depth',
    scale = 1, width = 250, height = 100, units = 'mm',
    dpi = 150, limitsize = TRUE)
  }
```
```
for (i in 1:227){
  paste(i)
  # contig = paste("Chr", i, sep = "")
  contig = paste(i, sep = "")
  p0 <- ggplot(data=appended_df[appended_df$contig == contig, ], aes(x=`position`, y=`depth_mod`, group=1)) +
    geom_line() +
    labs(x = "Position", y = "Coverage") +
    scale_y_continuous(breaks=seq(0,150,25), limits=c(0,150))
  outfile = paste("contig", i, "cov.tiff", sep = "_")
  ggsave(outfile , plot = p0, device = 'tiff', path = '/data/scratch/armita/idris/alignment/read_depth/',
    scale = 1, width = 250, height = 100, units = 'mm',
    dpi = 150, limitsize = TRUE)
  }
```

Contigs with 0 coverage were removed and contigs renamed

```bash
Exclude_contig_lines="contig_1\t3595319\t0\ncontig_117\t141214\t0\ncontig_16\t1005393\t0\ncontig_2\t2482273\t0"
printf \
"Exclude:\nSequence name,\tlength,\tapparent source
$Exclude_contig_lines
Trim:\nSequence name,\tlength,\tspan(s),\tapparent source\n"  \
> tmp.csv
for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/polished/contigs_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)/../filtered
mkdir -p $OutDir
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/filtered_contigs_renamed.fasta --coord_file tmp.csv
done
rm tmp.csv
```


```bash

for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/filtered/filtered_contigs_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


<!--

## Merging pacbio and hybrid assemblies

```bash
for PacBioAssembly in $(ls ../../../../home/groups/harrisonlab/project_files/idris/assembly/falcon/P.cactorum/414/2-asm-falcon/polished/contigs_min_500bp_renamed.fasta); do
Organism=$(echo $PacBioAssembly | rev | cut -f5 -d '/' | rev)
Strain=$(echo $PacBioAssembly | rev | cut -f4 -d '/' | rev)
# HybridAssembly=$(ls assembly/spades_pacbio/$Organism/$Strain/contigs.fasta)
HybridAssembly=$(ls ../../../../home/groups/harrisonlab/project_files/idris/assembly/spades_pacbio/P.cactorum/414/filtered_contigs/contigs_min_500bp.fasta)
OutDir=assembly/merged_canu_spades/$Organism/$Strain
AnchorLength=500000
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $PacBioAssembly $HybridAssembly $OutDir $AnchorLength
done
```

Quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/merged_canu_spades/P.cactorum/414/merged.fasta); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
OutDir=$(dirname $Assembly)
# OutDir=/data/scratch/armita/idris/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

checking using busco

```bash
for Assembly in $(ls assembly/merged_canu_spades/P.cactorum/414/merged.fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
OutDir=$(dirname $Assembly)
# OutDir=/data/scratch/armita/idris/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


Contigs were identified that had blast hits to non-phytophthora genomes

```bash
for Assembly in $(ls assembly/merged_canu_spades/P.cactorum/414/merged.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
# Exclude_db="bact,virus,hsref"
Exclude_db="paenibacillus"
Good_db="phytoph"
AssemblyDir=$(dirname $Assembly)
OutDir=$AssemblyDir/../deconseq_Paen
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
done
```

Results were summarised using the commands:

```bash
  for File in $(ls assembly/falcon/P.*/*/*/deconseq_Paen/log.txt | grep '414'); do
    Name=$(echo $File | rev | cut -f4 -d '/' | rev);
    Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
    Both=$(cat $File |cut -f2 | head -n2 | tail -n1);
    Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
    printf "$Name\t$Good\t$Both\t$Bad\n";
  done
```

```
  414_v2	222	0	5
```


## Checking PacBio coverage against P414 contigs

The accuracy of PacBio assembly pipelines is currently unknown. To help identify
regions that may have been missassembled the pacbio reads were aligned back to
the assembled genome. Coverage was determined using bedtools genomecov and
regions with low coverage flagged using a python script flag_low_coverage.py.
These low coverage regions were visually inspected using IGV.

```bash
  for Assembly in $(ls assembly/falcon/P.*/*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta | grep -e '414'); do
    Reads=$(ls raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq)
    OutDir=analysis/genome_alignment/bwa/P.cactorum/414/vs_414
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
    qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
  done
```
-->

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  # for BestAss in $(ls assembly/falcon/P.*/*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta | grep -e '414'); do
  for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/filtered/filtered_contigs_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
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

```bash
for RepDir in $(ls -d repeat_masked/P.*/*/filtered_contigs_repmask | grep -w '414'); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
# printf "The number of bases masked by RepeatMasker:\t"
RepMaskerBp=$(sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
# printf "The number of bases masked by TransposonPSI:\t"
TpsiBp=$(sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
# printf "The total number of masked bases are:\t"
Total=$(cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
printf "$Organism\t$Strain\t$RepMaskerBp\t$TpsiBp\t$Total\n"
done
```

```
  P.cactorum	414	18468452	5361346	19269912
```

# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

## Pre-gene prediction
<!--
Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  for Genome in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414_v2'); do
    echo $Genome;
    qsub $ProgDir/sub_cegma.sh $Genome dna;
  done
```

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/*/*/*_dna_cegma.completeness_report | grep -w -e '414_v2'); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
``` -->

#Gene prediction

Gene prediction was performed for the P. cactorum genome. Two gene prediction
approaches were used:

Gene prediction using Braker1
Prediction of all putative ORFs in the genome using the ORF finder (atg.pl)
approach.


## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to Fusarium genomes.



#### Aligning published RNAseq data

* qc of RNA seq data was performed as part of sequencing the 10300 genome:


```bash
for Assembly in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w '414'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNA in $(ls /home/groups/harrisonlab/project_files/idris/qc_rna/raw_rna/genbank/*/*/*_trim.fq.gz); do
Timepoint=$(echo $RNA | rev | cut -f1 -d '/' | rev | sed 's/_trim.*//g')
echo "$Timepoint"
Prefix="$Tiimepoint"
OutDir=alignment/star/$Organism/"$Strain"/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star_unpaired.sh $Assembly $RNA $OutDir
done
done
```


### Aligning in house RNAseq data

#### Timecourse of infection

make symbolic links to timecourse data

```bash
# Create symbolic links for all F read files
for File in $(ls /home/groups/harrisonlab/raw_data/raw_seq/fragaria/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07/*/*.gz | grep 'R1.fastq.gz'); do
  # echo $File;
  Sample=$(echo $File | rev | cut -d '/' -f2 | rev)
  echo "$Sample"
  OutDir=qc_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07/$Sample/F
  mkdir -p "$OutDir"
  cp -s $File $OutDir/.
done
# Create symbolic links for all R read files
for File in $(ls /home/groups/harrisonlab/raw_data/raw_seq/fragaria/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07/*/*.gz | grep 'R2.fastq.gz'); do
  # echo $File;
  Sample=$(echo $File | rev | cut -d '/' -f2 | rev)
  echo "$Sample"
  OutDir=qc_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07/$Sample/R
  mkdir -p "$OutDir"
  cp -s $File $OutDir/.
done
```


Perform qc of RNAseq timecourse data
```bash
  for FilePath in $(ls -d raw_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07/*); do
    echo $FilePath
    FileNum=$(ls $FilePath/F/*.gz | wc -l)
    for num in $(seq 1 $FileNum); do
      FileF=$(ls $FilePath/F/*.gz | head -n $num | tail -n1)
      FileR=$(ls $FilePath/R/*.gz | head -n $num | tail -n1)
      echo $FileF
      echo $FileR
      Jobs=$(qstat | grep 'rna_qc' | grep 'qw' | wc -l)
      while [ $Jobs -gt 16 ]; do
        sleep 5m
        printf "."
        Jobs=$(qstat | grep 'rna_qc' | grep 'qw' | wc -l)
      done		
      printf "\n"
      IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
      qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
    done
  done
```

Maria performed alignment of RNAseq data vs the F. vesca and those reads that did not align were
used for alignment vs the P.cactorum genome.


make symbolic links to timecourse data

```bash
  for File in $(ls /home/sobczm/popgen/rnaseq/vesca_*/*.mate1.fq.gz); do
    Sample=$(echo $File | rev | cut -d '/' -f2 | rev | sed 's/vesca_//g')
    echo "$Sample"
    OutDir=qc_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07_no_vesca/$Sample/F
    mkdir -p "$OutDir"
    cp -s $File $OutDir/.
  done
  for File in $(ls /home/sobczm/popgen/rnaseq/vesca_*/*.mate2.fq.gz); do
    Sample=$(echo $File | rev | cut -d '/' -f2 | rev | sed 's/vesca_//g')
    echo "$Sample"
    OutDir=qc_rna/paired/Transcriptome_Emily_Fenella_Pcactorum-2017-04-07_no_vesca/$Sample/R
    mkdir -p "$OutDir"
    cp -s $File $OutDir/.
  done
```

#### Mycelium

make symbolic links to mycelium data

```bash
# Create symbolic links for all F read files
for File in $(ls /data/seq_data/external/20171206_C101HW17030405_novogene/raw_data/CN_P414_*_1.fq.gz); do
  # echo $File;
  Sample="mycelium"
  echo "$Sample"
  OutDir=raw_rna/paired/mycelium-2017-12-06/$Sample/F
  mkdir -p "$OutDir"
  cp -s $File $OutDir/.
done
# Create symbolic links for all R read files
for File in $(ls /data/seq_data/external/20171206_C101HW17030405_novogene/raw_data/CN_P414_*_2.fq.gz); do
  # echo $File;
  Sample="mycelium"
  echo "$Sample"
  OutDir=raw_rna/paired/mycelium-2017-12-06/$Sample/R
  mkdir -p "$OutDir"
  cp -s $File $OutDir/.
done
```


Perform qc of RNAseq timecourse data
```bash
  for FilePath in $(ls -d raw_rna/paired/mycelium-2017-12-06/*); do
    echo $FilePath
    FileNum=$(ls $FilePath/F/*.gz | wc -l)
    for num in $(seq 1 $FileNum); do
      FileF=$(ls $FilePath/F/*.gz | head -n $num | tail -n1)
      FileR=$(ls $FilePath/R/*.gz | head -n $num | tail -n1)
      echo $FileF
      echo $FileR
      Jobs=$(qstat | grep 'rna_qc' | grep 'qw' | wc -l)
      while [ $Jobs -gt 16 ]; do
        sleep 5m
        printf "."
        Jobs=$(qstat | grep 'rna_qc' | grep 'qw' | wc -l)
      done		
      printf "\n"
      IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
      qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
    done
  done
```


Performed alignment of RNAseq data vs the F. vesca genome and those reads that did not align were
used for alignment vs the P.cactorum genome.


```bash
for Assembly in $(ls ../../../../home/sobczm/popgen/rnaseq/fvesca_v1.1_all.fa); do
echo "$Assembly"
for RNADir in $(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_rna/paired/*/* | grep -e 'mycelium'); do
FileNum=$(ls $RNADir/F/*.fq.gz | wc -l)
for num in $(seq 1 $FileNum); do
printf "\n"
FileF=$(ls $RNADir/F/*.fq.gz | head -n $num | tail -n1)
FileR=$(ls $RNADir/R/*.fq.gz | head -n $num | tail -n1)
echo $FileF
echo $FileR
Prefix=$(echo $FileF | rev | cut -f1 -d '/' | rev | sed "s/_1_trim.fq.gz//g")
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint"
OutDir=alignment/star/fvesca/v1.1/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star_unmapped.sh $Assembly $FileF $FileR $OutDir
done
done
done
```

```bash
for File in $(ls -d alignment/star/fvesca/v1.1/*/*/*.mate* | grep 'mycelium'); do
echo $File
cat $File | gzip -cf > $File.fq.gz
done
```


```bash
for Assembly in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w '414'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d alignment/star/fvesca/v1.1/*/* | grep 'mycelium'); do
FileF=$(ls $RNADir/*.mate1.fq.gz)
FileR=$(ls $RNADir/*.mate2.fq.gz)
echo $FileF
echo $FileR
Prefix=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
Timepoint=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
echo "$Timepoint"
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
done
```

#### Aligning


```bash
for Assembly in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w '414'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
# for RNADir in $(ls -d qc_rna/paired/Transcriptome*/*); do
# for RNADir in $(ls -d qc_rna/paired/mycelium-2017-12-06*/*); do
for RNADir in $(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_rna/paired/*/* | grep -e 'mycelium' -e '_no_vesca' | grep 'mycelium'); do
FileNum=$(ls $RNADir/F/*.mate1.fq.gz | wc -l)
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
for num in $(seq 1 $FileNum); do
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
done
printf "\n"
FileF=$(ls $RNADir/F/*.mate1.fq.gz | head -n $num | tail -n1)
FileR=$(ls $RNADir/R/*.mate2.fq.gz | head -n $num | tail -n1)
echo $FileF
echo $FileR
# Prefix=$(echo $FileF | rev | cut -f2 -d '/' | rev | sed "s/_R.*_trim.fq.gz//g")
Prefix=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint"
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
done
done
```

Alignment stats were collected:

```bash
for File in $(ls alignment/star/P.cactorum/414/*/*/star_aligmentLog.final.out); do
Sample=$(echo $File | rev | cut -f2 -d '/' | rev);
ReadNumU=$(cat $File | grep 'Uniquely' | grep 'number' | cut -f2);
ReadPercU=$(cat $File | grep 'Uniquely' | grep '%' | cut -f2);
ReadNumM=$(cat $File | grep 'multiple' | grep 'Number' | cut -f2);
ReadPercM=$(cat $File | grep 'multiple' | grep '%' | cut -f2);
echo -e "$Sample""\t""$ReadNumU""\t""$ReadPercU""\t""$ReadNumM""\t""$ReadPercM";  
done
```

```
CN_P414_1	27465450	89.31%	1595938	5.19%
CN_P414_2	26290608	89.01%	1484087	5.02%
CN_P414_3	26550405	88.81%	1853630	6.20%
PRO1467_S10_totRNA_S15_L001	5465	0.24%	805	0.04%
PRO1467_S10_totRNA_S15_L002	5469	0.19%	718	0.03%
PRO1467_S10_totRNA_S15_L003	6162	0.21%	772	0.03%
PRO1467_S10_totRNA_S15_L004	6003	0.21%	777	0.03%
PRO1467_S11_totRNA_S9_L001	3245	0.15%	539	0.02%
PRO1467_S11_totRNA_S9_L002	2557	0.11%	337	0.01%
PRO1467_S11_totRNA_S9_L003	2341	0.10%	335	0.01%
PRO1467_S11_totRNA_S9_L004	2618	0.11%	349	0.01%
PRO1467_S12_totRNA_S3_L001	5396	0.17%	763	0.02%
PRO1467_S12_totRNA_S3_L002	3793	0.12%	475	0.02%
PRO1467_S12_totRNA_S3_L003	3977	0.12%	476	0.01%
PRO1467_S12_totRNA_S3_L004	4109	0.13%	509	0.02%
PRO1467_S14_totRNA_S14_L001	80526	4.16%	7401	0.38%
PRO1467_S14_totRNA_S14_L002	117497	4.83%	9677	0.40%
PRO1467_S14_totRNA_S14_L003	125985	4.97%	10116	0.40%
PRO1467_S14_totRNA_S14_L004	124981	4.99%	10236	0.41%
PRO1467_S15_totRNA_S8_L001	121857	4.03%	10418	0.34%
PRO1467_S15_totRNA_S8_L002	195953	6.32%	15048	0.49%
PRO1467_S15_totRNA_S8_L003	150655	4.84%	11445	0.37%
PRO1467_S15_totRNA_S8_L004	202324	6.49%	14983	0.48%
PRO1467_S17_totRNA_S13_L001	2512253	40.07%	468769	7.48%
PRO1467_S17_totRNA_S13_L002	3555861	43.05%	573642	6.94%
PRO1467_S17_totRNA_S13_L003	2848182	32.67%	449902	5.16%
PRO1467_S17_totRNA_S13_L004	3778638	43.56%	595166	6.86%
PRO1467_S18_totRNA_S7_L001	2191645	40.28%	234652	4.31%
PRO1467_S18_totRNA_S7_L002	2262486	44.04%	226904	4.42%
PRO1467_S18_totRNA_S7_L003	1763092	33.38%	176386	3.34%
PRO1467_S18_totRNA_S7_L004	2342978	44.58%	233515	4.44%
PRO1467_S1_totRNA_S17_L001	4299	0.14%	622	0.02%
PRO1467_S1_totRNA_S17_L002	3977	0.13%	438	0.01%
PRO1467_S1_totRNA_S17_L003	5313	0.16%	550	0.02%
PRO1467_S1_totRNA_S17_L004	2815	0.09%	350	0.01%
PRO1467_S20_totRNA_S10_L001	3151625	46.77%	269889	4.00%
PRO1467_S20_totRNA_S10_L002	1837296	25.00%	150088	2.04%
PRO1467_S20_totRNA_S10_L003	3841586	49.86%	306696	3.98%
PRO1467_S20_totRNA_S10_L004	3828035	49.98%	305682	3.99%
PRO1467_S22_totRNA_S2_L001	492223	12.92%	51032	1.34%
PRO1467_S22_totRNA_S2_L002	490105	14.25%	48522	1.41%
PRO1467_S22_totRNA_S2_L003	521262	14.29%	50337	1.38%
PRO1467_S22_totRNA_S2_L004	516600	14.47%	49881	1.40%
PRO1467_S28_totRNA_S1_L001	3793687	50.50%	385999	5.14%
PRO1467_S28_totRNA_S1_L002	3907469	52.48%	378923	5.09%
PRO1467_S28_totRNA_S1_L003	4180466	52.64%	400397	5.04%
PRO1467_S28_totRNA_S1_L004	3116156	39.90%	301905	3.87%
PRO1467_S2_totRNA_S12_L001	4556	0.21%	768	0.04%
PRO1467_S2_totRNA_S12_L002	4237	0.17%	626	0.03%
PRO1467_S2_totRNA_S12_L003	3166	0.11%	480	0.02%
PRO1467_S2_totRNA_S12_L004	4478	0.17%	722	0.03%
PRO1467_S3_totRNA_S6_L001	22783	0.81%	2424	0.09%
PRO1467_S3_totRNA_S6_L002	25983	1.03%	2337	0.09%
PRO1467_S3_totRNA_S6_L003	27596	1.04%	2512	0.09%
PRO1467_S3_totRNA_S6_L004	27962	1.07%	2434	0.09%
PRO1467_S4_totRNA_S18_L001	145209	4.41%	11657	0.35%
PRO1467_S4_totRNA_S18_L002	221597	6.46%	16007	0.47%
PRO1467_S4_totRNA_S18_L003	235106	6.52%	17104	0.47%
PRO1467_S4_totRNA_S18_L004	234879	6.62%	16730	0.47%
PRO1467_S5_totRNA_S11_L001	22302	0.91%	2239	0.09%
PRO1467_S5_totRNA_S11_L002	97935	4.22%	8374	0.36%
PRO1467_S5_totRNA_S11_L003	105174	4.29%	8724	0.36%
PRO1467_S5_totRNA_S11_L004	78941	3.29%	6931	0.29%
PRO1467_S6_totRNA_S5_L001	165264	6.16%	12503	0.47%
PRO1467_S6_totRNA_S5_L002	192282	6.71%	13564	0.47%
PRO1467_S6_totRNA_S5_L003	204377	6.79%	13921	0.46%
PRO1467_S6_totRNA_S5_L004	202085	6.84%	13881	0.47%
PRO1467_S7_totRNA_S16_L001	3112318	40.42%	373807	4.85%
PRO1467_S7_totRNA_S16_L002	4520087	43.02%	483469	4.60%
PRO1467_S7_totRNA_S16_L003	5027169	43.46%	523155	4.52%
PRO1467_S7_totRNA_S16_L004	4937234	43.46%	516243	4.54%
PRO1467_S9_totRNA_S4_L001	4316377	58.38%	391297	5.29%
PRO1467_S9_totRNA_S4_L002	5139390	61.43%	438245	5.24%
PRO1467_S9_totRNA_S4_L003	5471977	61.63%	461445	5.20%
PRO1467_S9_totRNA_S4_L004	5438415	61.83%	458678	5.21%
```

Alignments were concatenated prior to gene prediction

```bash
BamFiles=$(ls alignment/star/P.cactorum/414/*/*/star_aligmentAligned.sortedByCoord.out.bam | grep -v -e 'PRO1467_S1_' -e 'PRO1467_S2_' -e 'PRO1467_S3_' -e 'PRO1467_S10_' -e 'PRO1467_S11_' -e 'PRO1467_S12_' | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=alignment/star/P.cactorum/414/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/concatenated.bam $BamFiles
```

#### Braker prediction

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_softmasked_repeatmasker_TPSI_appended.fa | grep -w -e '414'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
OutDir=gene_pred/braker/$Organism/"$Strain"_braker_pacbio
GeneModelName="$Organism"_"$Strain"_braker_pacbio
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_pacbio
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker.sh $Assembly $OutDir $AcceptedHits $GeneModelName
# OutDir=gene_pred/braker/$Organism/"$Strain"_braker_pacbio_fungi
# qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_softmasked_repeatmasker_TPSI_appended.fa | grep -w -e '414'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
    mkdir -p $OutDir
    AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
  done
```

Secondly, genes were predicted using CodingQuary:

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_softmasked_repeatmasker_TPSI_appended.fa | grep -w -e '414'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/codingquary/$Organism/$Strain
CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```


The next step had problems with the masked pacbio genome. Bioperl could not read in
the fasta sequences. This was overcome by wrapping the unmasked genome and using this
fasta file.

```bash
Assembly=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
NewName=$(echo $Assembly | sed 's/_unmasked.fa/_unmasked_wrapped.fa/g')
cat $Assembly | fold > $NewName
```


Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
for BrakerGff in $(ls gene_pred/braker/*/*_braker_pacbio/*/augustus.gff3 | grep '414'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_pacbio//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_unmasked_wrapped.fa)
CodingQuaryGff=$(ls gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3)
PGNGff=$(ls gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3)
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/final/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
# -
# This section is edited
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuary_unspliced.gff3
$ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuary_unspliced.gff3 > $FinalDir/final_genes_CodingQuary.gff3
# -
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_appended_renamed.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_appended_renamed.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_appended_renamed.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_appended_renamed.upstream3000.fasta


GffBraker=$FinalDir/final_genes_Braker.gff3
GffQuary=$FinalDir/final_genes_CodingQuary.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended

# cat $BrakerGff $AddDir/additional_gene_parsed.gff3 | bedtools sort > $FinalGff
done
```


The final number of genes per isolate was observed using:
```bash
  for DirPath in $(ls -d gene_pred/final/P.*/*/final | grep '414'); do
    echo $DirPath;
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_appended_renamed_renamed.pep.fasta | grep '>' | wc -l;
    echo "";
  done
```


In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified


```bash
for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep '414'); do
Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FinalDir=gene_pred/final/$Organism/$Strain/final
GffFiltered=$FinalDir/filtered_duplicates.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
LogFile=$FinalDir/final_genes_appended_renamed.log
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
rm $GffFiltered
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_unmasked_wrapped.fa)
$ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed
# The proteins fasta file contains * instead of Xs for stop codons, these should
# be changed
sed -i 's/\*/X/g' gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta
done
```
```
P.cactorum - 414
Identifiied the following duplicated transcripts:
CUFF_3035_1_82.t2
CUFF_1672_1_4.t2
CUFF_345_1_149.t2
```

The final number of genes per isolate was observed using:
```bash
  for DirPath in $(ls -d gene_pred/final/P.*/*/final | grep '414'); do
    echo $DirPath;
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_appended_renamed.pep.fasta | grep '>' | wc -l;
    echo "";
  done
```

```
gene_pred/final/P.cactorum/414/final
17472
12419
29891
```

<!--
# Assessing gene space in predicted transcriptomes

```bash
  for Transcriptome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gene.fasta | grep '414_v2'); do
    Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    # BuscoDB="Fungal"
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
    OutDir=gene_pred/busco/$Organism/$Strain/genes
    qsub $ProgDir/sub_busco2.sh $Transcriptome $BuscoDB $OutDir
  done
```
```bash
  for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt | grep '414_v2'); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

```
P.cactorum	414_v2	277	4	22	303
``` -->

## Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing
Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w -e '414'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
  	qsub $ProgDir/run_ORF_finder.sh $Assembly
  done
```


The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
	ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
  for OrfGff in $(ls gene_pred/ORF_finder/P.*/*/*_ORF.gff | grep -v 'atg' | grep -w -e '414'); do
    echo "$OrfGff"
  	OrfGffMod=$(echo $OrfGff | sed 's/.gff/.gff3/g')
  	$ProgDir/gff_corrector.pl $OrfGff > $OrfGffMod
  done
```

#Genomic analysis

## RxLR genes

Putative RxLR genes were identified within Augustus gene models using a number
of approaches:

 * A) From Augustus gene models - Signal peptide & RxLR motif  
 * B) From Augustus gene models - Hmm evidence of WY domains  
 * C) From Augustus gene models - Hmm evidence of RxLR effectors
 * D) From Augustus gene models - Hmm evidence of CRN effectors  
 * E) From ORF fragments - Signal peptide & RxLR motif  
 * F) From ORF fragments - Hmm evidence of WY domains  
 * G) From ORF fragments - Hmm evidence of RxLR effectors


 ### A) From Augustus gene models - Signal peptide & RxLR motif

 Required programs:
  * SigP
  * biopython

#### A.1) Signal peptide prediction using SignalP 2.0

 Proteins that were predicted to contain signal peptides were identified using
 the following commands:

 ```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta); do
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/final_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_final
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_final_*); do
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
while [ $Jobs -gt 20 ]; do
sleep 1
printf "."
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
done
printf "\n"
echo $File
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
qsub $ProgDir/pred_sigP.sh $File
qsub $ProgDir/pred_sigP.sh $File signalp-3.0
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```


The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:

```bash
 for SplitDir in $(ls -d gene_pred/final_split/P.*/* | grep '414'); do
   Strain=$(echo $SplitDir | cut -d '/' -f4)
   Organism=$(echo $SplitDir | cut -d '/' -f3)
   echo "$Organism - $Strain"
   for SigpDir in $(ls -d gene_pred/final_sig* | cut -f2 -d'/'); do
     InStringAA=''
     InStringNeg=''
     InStringTab=''
     InStringTxt=''
     for GRP in $(ls -l $SplitDir/*_final_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
       InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.aa";
       InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp_neg.aa";  
       InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.tab";
       InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_$GRP""_sp.txt";
     done
     cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
     cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
     tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
     cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
   done
 done
```

#### B.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '414'); do
   Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   echo "$Organism - $Strain"
   OutDir=analysis/phobius/$Organism/$Strain
   mkdir -p $OutDir
   phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
   ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
   $ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
   cat $OutDir/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d '>' > $OutDir/"$Strain"_phobius_headers.txt
 done
```


Secreted proteins from different sources were combined into a single file:

```bash
  for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '414'); do
   Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   echo "$Organism - $Strain"
   OutDir=gene_pred/combined_sigP/$Organism/$Strain
   mkdir -p $OutDir
   echo "The following number of sequences were predicted as secreted:"
   cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
   cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | wc -l
   echo "This represented the following number of unique genes:"
   cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
   cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
  done
```

```
P.cactorum - 414
The following number of sequences were predicted as secreted:
10583
This represented the following number of unique genes:
3465
```

Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
  for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '414'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
    qsub $ProgDir/submit_TMHMM.sh $Proteome
  done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt | grep '414'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
NonTmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $NonTmHeaders
SigP=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_secreted.fa)
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $NonTmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
echo "Number of SigP proteins:"
cat $SigP | grep '>' | wc -l
echo "Number without transmembrane domains:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
echo "Number of gene models:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l

# A text file was also made containing headers of proteins testing +ve
PosFile=$(ls gene_pred/trans_mem/$Organism/$Strain/"$Strain"_TM_genes_pos.txt)
TmHeaders=$(echo $PosFile | sed 's/.txt/_headers.txt/g')
cat $PosFile | cut -f1 > $TmHeaders

done
```

```
P.cactorum - 414
Number of SigP proteins:
3465
Number without transmembrane domains:
2399
Number of gene models:
2386
```

Proteins containing GPI anchors were also removed using GPIsom


These proteins were identified through submitting the combined protein file to
the webserver at: http://gpi.unibe.ch


An output directory was made to download the file to:
"GPI anchored (C&N-term signal) (SignalP):"

```bash
 for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '414'); do
   Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   echo "$Organism - $Strain"
   OutDir=gene_pred/trans_mem/$Organism/$Strain/GPIsom
   mkdir -p $OutDir
 done
```

The link to the proteins results was followed embedded in the text of
"GPI anchored (C&N-term signal) (SignalP):" thes +GPI anchor proteins
were copied and pasted into:

<!-- The link to the results was followed embedded in the text of
"Seqs with C-terminal signal (GPI-SOM):" thes +GPI anchor proteins
were copied and pasted into: -->

```bash
 nano gene_pred/trans_mem/P.cactorum/414/GPIsom/GPI_pos.fa
```

Those proteins with GPI anchors were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/GPIsom/GPI_pos.fa | grep '414'); do
Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
TmHeaders=$(echo "$File" | sed 's/.fa/.txt/g')
cat $File | grep '>' | cut -f1 -d ' ' | sed 's/>//g' > $TmHeaders
SigP=$(ls gene_pred/combined_sigP/$Organism/$Strain/*_sp_no_trans_mem.aa)
SigPHeaders=gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_sp_no_trans_mem_headers.txt
cat $SigP | grep '>' | cut -f1 | sed 's/>//g'> $SigPHeaders
GoodHeaders=$(echo "$File" | sed 's/_pos.fa/_neg.txt/g')
cat $SigPHeaders | grep -v -f $TmHeaders > $GoodHeaders
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# cat $SigP | grep -v -A1 -f $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $GoodHeaders  > $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa
echo "Number of SigP proteins:"
cat $SigP | grep '>' |wc -l
echo "Number with GPI anchors in entire proteome:"
cat $TmHeaders | wc -l
echo "Number without GPI anchors:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa | grep '>' | wc -l
echo "Number of gene models:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l
done
```

```
  P.cactorum - 414
  Number of SigP proteins:
  2399
  Number with GPI anchors in entire proteome:
  509
  Number without GPI anchors:
  2129
  Number of gene models:
  2118
```



### C) From Augustus gene models - Effector-like structure identification using EffectorP

Required programs:
* EffectorP.py

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '414'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
BaseName="$Organism"_"$Strain"_EffectorP
OutDir=analysis/effectorP/$Organism/$Strain
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
done
```

Those genes that were predicted as secreted and tested positive by effectorP
were identified:

Note - this doesnt exclude proteins with TM domains or GPI anchors

```bash
 for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt | grep '414'); do
   Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
   Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
   echo "$Organism - $Strain"
   Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
   cat $File | grep 'Effector' | cut -f1 > $Headers
   printf "EffectorP headers:\t"
   cat $Headers | wc -l
   Secretome=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_secreted.fa)
   OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
   OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
   cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
   printf "Secreted effectorP headers:\t"
   cat $OutFileHeaders | wc -l
   Gff=$(ls gene_pred/final/$Organism/$Strain/*/final_genes_appended_renamed.gff3)
   EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
 done
```
```
P.cactorum - 414
EffectorP headers:	19602
Secreted effectorP headers:	925
```

## D) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '414'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/CAZY/$Organism/$Strain
mkdir -p $OutDir
Prefix="$Strain"_CAZY
CazyHmm=../../dbCAN/dbCAN-fam-HMMs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep '414'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $File)
echo "$Organism - $Strain"
ProgDir=/home/groups/harrisonlab/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
printf "number of CAZY genes identified:\t"
cat $CazyHeaders | wc -l
Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
CazyGff=$OutDir/"$Strain"_CAZY.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

SecretedProts=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_secreted.fa)
# SecretedProts=$(ls gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem_no_GPI.aa)
SecretedHeaders=$(echo $SecretedProts | sed 's/.fa/_headers.txt/g' | sed 's/.aa/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
printf "number of Secreted CAZY genes identified:\t"
cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
done
```

```
P.cactorum - 10300
number of CAZY genes identified:	697
number of Secreted CAZY genes identified:	372
```

or, if those secreted proteins with TM domains or GPI anchors are removed:

```
P.cactorum - 10300
number of CAZY genes identified:	697
number of Secreted CAZY genes identified:	236
```


Note - the CAZY genes identified may need further filtering based on e value and
cuttoff length - see below:

Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage

* For fungi, use E-value < 1e-17 and coverage > 0.45

* The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)

### Summary of CAZY families by organism


```bash
for CAZY in $(ls gene_pred/CAZY/*/*/*_CAZY.out.dm.ps | grep '10300'); do
 Strain=$(echo $CAZY | rev | cut -f2 -d '/' | rev)
 Organism=$(echo $CAZY | rev | cut -f3 -d '/' | rev)
 OutDir=$(dirname $CAZY)
 echo "$Organism - $Strain"
 Secreted=$(ls gene_pred/combined_sigP/$Organism/$Strain/*_sp_no_trans_mem_no_GPI_headers.txt)
 Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
 ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/CAZY
 $ProgDir/summarise_CAZY.py --cazy $CAZY --inp_secreted $Secreted --inp_gff $Gff --summarise_family --trim_gene_id 2 --kubicek_2014
done
```

```
 P.cactorum - 10300
 Polygalacturonase - 17
 A-Arabinosidases - 5
 Xylanases - 12
 Polygalacturonate lyases - 34
 B-Glycosidases - 16
 Cellulases - 26
 other - 126
```

The following commands were used to identify b-glucan binding proteins as described in:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3627985/

```bash
 for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep '10300'); do
   cat $File | grep -w -e 'GH5' -e 'GH55' -e 'GH81' -e 'GH16' -e 'GH3' -e 'GH17' -e 'GH72' | sed -r "s/\s+/\t/g" | cut -f1,4 |  sort | uniq > tmp1.txt
   cat tmp1.txt | cut -f2 > tmp2.txt
   cat tmp1.txt | grep -f tmp3.txt | cut -f1 | sort | uniq -c
 done
```

```
 2 GH16.hmm
 18 GH17.hmm
 14 GH3.hmm
 14 GH5.hmm
 11 GH81.hmm
```

# D) Prediction of RxLRs from - Augustus/CodingQuary gene models

The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep '414'); do
Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
Proteome=$(ls gene_pred/final/$Organism/$Strain/*/final_genes_appended_renamed.pep.fasta)
Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain"
mkdir -p $OutDir;
printf "\nstrain: $Strain\tspecies: $Organism\n";
printf "the total number of SigP gene is:\t";
cat $Secretome | grep '>' | wc -l;
printf "the number of unique SigP gene is:\t";
cat $Secretome | grep '>' | cut -f1 | tr -d ' '| sort | uniq | wc -l;

printf "the number of SigP-RxLR genes are:\t";
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_all_secreted_RxLR_regex.fa;
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $OutDir/"$Strain"_RxLR_regex.txt
cat $OutDir/"$Strain"_RxLR_regex.txt | wc -l

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa


printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $OutDir/"$Strain"_RxLR_EER_regex.txt
cat $OutDir/"$Strain"_RxLR_EER_regex.txt | wc -l
$ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa

printf "\n"

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_regex.txt
sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_EER_regex.txt

cat $Gff | grep -w -f $OutDir/"$Strain"_RxLR_regex.txt> $OutDir/"$Strain"_RxLR_regex.gff3
cat $Gff | grep -w -f $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.gff3
done
```

```
  strain: 414	species: P.cactorum
  the total number of SigP gene is:	10583
  the number of unique SigP gene is:	3465
  the number of SigP-RxLR genes are:	313
  the number of SigP-RxLR-EER genes are:	133
```


### G) From Predicted gene models - Hmm evidence of RxLR effectors

```bash
 for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep '414'); do
   ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
   HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
   Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
   mkdir -p $OutDir
   HmmResults="$Strain"_RxLR_hmmer.txt
   hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
   echo "$Organism $Strain"
   cat $OutDir/$HmmResults | grep 'Initial search space'
   cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
   HmmFasta="$Strain"_RxLR_hmmer.fa
   $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
   Headers="$Strain"_RxLR_hmmer_headers.txt
   cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' | sort | uniq > $OutDir/$Headers
   Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
   cat $Gff | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
 done
```

```
P.cactorum 414
Initial search space (Z):              29888  [actual number of targets]
Domain search space  (domZ):             147  [number of targets reported over threshold]
```


### F) Combining RxLRs from Regex and hmm searches


The total RxLRs are

```bash
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_RxLR_EER_regex.txt | grep -v -e 'Aug' -e 'ORF' | grep '414'); do
Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
Proteome=$(ls gene_pred/final/$Organism/$Strain/*/final_genes_appended_renamed.pep.fasta)
HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer_headers.txt | grep -v -e 'Aug' -e 'ORF')
echo "$Organism - $Strain"
echo "Number of RxLRs identified by Regex:"
cat $RegexRxLR | sort | uniq | wc -l
echo "Number of RxLRs identified by Hmm:"
cat $HmmRxLR | sort | uniq | wc -l
echo "Number of RxLRs in combined dataset:"
cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
# echo "Number of RxLRs in both datasets:"
# cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
echo ""
# echo "Extracting RxLRs from datasets"
OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
mkdir -p $OutDir
cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_RxLR_headers.txt
Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
cat $Gff | grep -w -f $OutDir/"$Strain"_total_RxLR_headers.txt > $OutDir/"$Strain"_total_RxLR.gff
echo "Number of genes in the extracted gff file:"
cat $OutDir/"$Strain"_total_RxLR.gff | grep -w 'gene' | wc -l
done
```

```
P.cactorum - 10300
P.cactorum - 414
Number of RxLRs identified by Regex:
133
Number of RxLRs identified by Hmm:
147
Number of RxLRs in combined dataset:
170
Number of genes in the extracted gff file:
170
```



```bash
for Secretome in $(ls gene_pred/combined_sigP/*/*/*_all_secreted.fa | grep '414'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_Aug_WY_hmmer.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_Aug_WY_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_Aug_WY_hmmer_headers.txt
# cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' | sort | uniq > $OutDir/$Headers
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | cut -f1 | sort | uniq > $OutDir/$Headers
Gff=$(ls gene_pred/*/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
cat $Gff | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_WY_hmmer.gff
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/$Headers $Gff $HmmModel Name Augustus > $OutDir/"$Strain"_Aug_WY_hmmer.gff
done
```
```
P.cactorum 414
Initial search space (Z):              10583  [actual number of targets]
Domain search space  (domZ):             408  [number of targets reported over threshold]
```

<!-- 72 of 137 RxLRs were in the predicted 797 effectorP genes.
```bash
 cat analysis/RxLR_effectors/combined_evidence/P.idaei/SCRP370/SCRP370_total_RxLR_headers.txt $OutFileHeaders | sort | cut -f1 -d '.' | uniq -d | wc -l
``` -->

### D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in Augustus gene models. This was done with the following commands:


```bash
 HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
 LFLAK_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm)
 DWL_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm)
 for Proteome in $(ls gene_pred/final/*/*/*/final_genes_appended_renamed.pep.fasta | grep -w -e '414'); do
   Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
   mkdir -p $OutDir
   echo "$Organism - $Strain"
   # Run hmm searches LFLAK domains
   CrinklerProts_LFLAK=$OutDir/"$Strain"_pub_CRN_LFLAK_hmm.txt
   hmmsearch -T0 $LFLAK_hmm $Proteome > $CrinklerProts_LFLAK
   cat $CrinklerProts_LFLAK | grep 'Initial search space'
   cat $CrinklerProts_LFLAK | grep 'number of targets reported over threshold'
   ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
   $ProgDir/hmmer2fasta.pl $CrinklerProts_LFLAK $Proteome > $OutDir/"$Strain"_pub_CRN_LFLAK_hmm.fa
   # Run hmm searches DWL domains
   CrinklerProts_DWL=$OutDir/"$Strain"_pub_CRN_DWL_hmm.txt
   hmmsearch -T0 $DWL_hmm $Proteome > $CrinklerProts_DWL
   cat $CrinklerProts_DWL | grep 'Initial search space'
   cat $CrinklerProts_DWL | grep 'number of targets reported over threshold'
   ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
   $ProgDir/hmmer2fasta.pl $CrinklerProts_DWL $Proteome > $OutDir/"$Strain"_pub_CRN_DWL_hmm.fa
   # Identify the genes detected in both models
   cat $OutDir/"$Strain"_pub_CRN_LFLAK_hmm.fa $OutDir/"$Strain"_pub_CRN_DWL_hmm.fa | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $OutDir/"$Strain"_pub_CRN_LFLAK_DWL.txt
   cat $OutDir/"$Strain"_pub_CRN_LFLAK_DWL.txt | wc -l
 done
```


```
P.cactorum - 414
Initial search space (Z):              29888  [actual number of targets]
Domain search space  (domZ):             171  [number of targets reported over threshold]
Initial search space (Z):              29888  [actual number of targets]
Domain search space  (domZ):             159  [number of targets reported over threshold]
123
```
<!--
23 of 92 P.idaei CRNs were in the effectorP dataset.
```bash
cat analysis/CRN_effectors/hmmer_CRN/P.idaei/SCRP370/SCRP370_pub_CRN_LFLAK_DWL.txt $OutFileHeaders | sort | cut -f1 -d '.' | uniq -d | wc -l
``` -->

Extract gff annotations for Crinklers:

```bash
 for CRNlist in $(ls analysis/CRN_effectors/hmmer_CRN/*/*/*_pub_CRN_LFLAK_DWL.txt | grep -e '414'); do
   Strain=$(echo $CRNlist | rev | cut -f2 -d '/' | rev)
   Organism=$(echo $CRNlist | rev | cut -f3 -d '/' | rev)
   OutName=$(echo $CRNlist | sed 's/.txt/.gff/g')
   echo "$Organism - $Strain"
   Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
   cat $CRNlist | sed -r 's/\.t.$//g' > tmp.txt
   cat $Gff | grep -w -f tmp.txt > $OutName
   rm tmp.txt
 done
```



### E) From ORF gene models - Signal peptide & RxLR motif

Required programs:
* SigP
* Phobius
* biopython


#### E.1) Prediction using SignalP
Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep '414'); do
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
SplitDir=gene_pred/ORF_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_ORF_preds
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_ORF_preds_*); do
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
while [ $Jobs -gt 15 ]; do
sleep 1
printf "."
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
done		
printf "\n"
echo $File
qsub $ProgDir/pred_sigP.sh $File
qsub $ProgDir/pred_sigP.sh $File signalp-3.0
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```


The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:

```bash
for SplitDir in $(ls -d gene_pred/ORF_split/*/* | grep '414'); do
Strain=$(echo $SplitDir | cut -d '/' -f4)
Organism=$(echo $SplitDir | cut -d '/' -f3)
echo "$Organism - $Strain"
for SigpDir in $(ls -d gene_pred/ORF_sig* | cut -f2 -d'/'); do
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
for GRP in $(ls -l $SplitDir/*_ORF_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.aa";
InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp_neg.aa";
InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.tab";
InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.txt";
done
cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_ORF_sp.aa
cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_ORF_neg_sp.aa
tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_ORF_sp.tab
cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_ORF_sp.txt
done
done
```


#### E.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
	for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep '414'); do
   Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   echo "$Organism - $Strain"
   OutDir=analysis/phobius/$Organism/$Strain
   mkdir -p $OutDir
 	phobius.pl $Proteome > $OutDir/"$Strain"_phobius_ORF.txt
 	cat $OutDir/"$Strain"_phobius_ORF.txt | grep -B1 'SIGNAL' | grep 'ID' | sed s'/ID.*g/g/g' > $OutDir/"$Strain"_phobius_headers_ORF.txt
 done
```

Secreted proteins from different sources were combined into a single file:

```bash
 for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa | grep '414'); do
   Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   echo "$Organism - $Strain"
   OutDir=gene_pred/combined_sigP_ORF/$Organism/$Strain
   mkdir -p $OutDir
   echo "The following number of sequences were predicted as secreted:"
   cat gene_pred/ORF_sig*/$Organism/$Strain/*_ORF_sp.aa > $OutDir/"$Strain"_all_secreted.fa
   cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | tr -d '>' | tr -d ' ' | sed "s/HMM_score\t/HMM_score=\t/g" > $OutDir/"$Strain"_all_secreted_headers.txt
   cat $OutDir/"$Strain"_all_secreted_headers.txt | wc -l
   echo "This represented the following number of unique genes:"
   cat gene_pred/ORF_sig*/$Organism/$Strain/*_ORF_sp.aa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
   cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
 done
```

```
P.cactorum - 414
The following number of sequences were predicted as secreted:
66113
This represented the following number of unique genes:
32834
```



#### E.3) Prediction of RxLRs


Names of ORFs containing signal peptides were extracted from fasta files. This
included information on the position and hmm score of RxLRs.

```bash
	FastaFile=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp.aa
	SigP_headers=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_names.txt
	cat $FastaFile | grep '>' | sed -r 's/>//g' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | sed 's/--//g' > $SigP_headers
```

Due to the nature of predicting ORFs, some features overlapped with one another.
A single ORF was selected from each set of overlapped ORFs. This was was
selected on the basis of its SignalP Hmm score. Biopython was used to identify
overlaps and identify the ORF with the best signalP score.

Due to the nature of predicting ORFs, some features overlapped with one another.
A single ORF was selected from each set of overlapped ORFs. This was was
selected on the basis of its SignalP Hmm score. Biopython was used to identify
overlaps and identify the ORF with the best signalP score.

```bash
 for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff3 | grep -v -e 'atg' | grep '414'); do
   Organism=$(echo $ORF_Gff | rev |  cut -d '/' -f3 | rev) ;
   Strain=$(echo $ORF_Gff | rev | cut -d '/' -f2 | rev);
   OutDir=$(ls -d gene_pred/combined_sigP_ORF/$Organism/$Strain)
   echo "$Organism - $Strain"
   # SigP_fasta=$(ls $OutDir/"$Strain"_all_secreted.fa)
   SigP_headers=$(ls $OutDir/"$Strain"_all_secreted_headers.txt)
   ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)

   SigP_Gff=$OutDir/"$Strain"_all_secreted_unmerged.gff
   SigP_Merged_Gff=$OutDir/"$Strain"_all_secreted_merged.gff
   SigP_Merged_txt=$OutDir/"$Strain"_all_secreted_merged.txt
   SigP_Merged_AA=$OutDir/"$Strain"_all_secreted_merged.aa

   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/extract_gff_for_sigP_hits.pl $SigP_headers $ORF_Gff SigP Name > $SigP_Gff
   ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
   $ProgDir/make_gff_database.py --inp $SigP_Gff --db sigP_ORF.db
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
   $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF.db --id sigP_ORF --out sigP_ORF_merged.db --gff > $SigP_Merged_Gff
   cat $SigP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d'=' | rev > $SigP_Merged_txt
   # $ProgDir/extract_from_fasta.py --fasta $SigP_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
   $ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
 done
```

The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep '414'); do
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev);
Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev) ;
OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
mkdir -p $OutDir;
printf "\nstrain: $Strain\tspecies: $Organism\n";
printf "the number of SigP gene is:\t";
cat $Secretome | grep '>' | cut -f1 | tr -d ' ' | sort | uniq | wc -l
printf "the number of SigP-RxLR genes are:\t";
$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa;
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt
cat $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l
printf "the number of SigP-RxLR-EER genes are:\t";
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' '> $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt
cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt  $SigP_Merged_Gff 	RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
SigP_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt $SigP_Gff	RxLR_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt  $SigP_Gff	RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.gff
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt $SigP_Merged_Gff RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex.gff
RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_regex_merged.gff
RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_regex_merged.txt
RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_regex_merged.aa
RxLR_EER_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.gff
RxLR_EER_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.txt
RxLR_EER_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff --db sigP_ORF_RxLR.db
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.gff --db sigP_ORF_RxLR_EER.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR.db --id sigP_ORF_RxLR --out sigP_ORF_RxLR_merged.db --gff > $RxLR_Merged_Gff
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR_EER.db --id sigP_ORF_RxLR_EER --out sigP_ORF_RxLR_EER_merged.db --gff > $RxLR_EER_Merged_Gff
cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
cat $RxLR_EER_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_EER_Merged_txt
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_EER_Merged_txt > $RxLR_EER_Merged_AA
printf "Merged RxLR regex proteins:\t"
cat $RxLR_Merged_AA | grep '>' | wc -l
printf "Merged RxLR-EER regex proteins:\t"
cat $RxLR_EER_Merged_AA | grep '>' | wc -l
printf "\n"
done
```

```
strain: 414	species: P.cactorum
the number of SigP gene is:	32834
the number of SigP-RxLR genes are:	1924
the number of SigP-RxLR-EER genes are:	223
Merged RxLR regex proteins:	1613
Merged RxLR-EER regex proteins:	196
```

Quantification of ORF RxLRs was also performed.


#### Aligning


```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FileF=$(ls qc_rna/raw_rna/genbank/P.cactorum/F/SRR1206032_trim.fq.gz)
FileR=$(ls qc_rna/raw_rna/genbank/P.cactorum/R/SRR1206033_trim.fq.gz)
echo $FileF
echo $FileR
Prefix="genbank"
Timepoint="treatment"
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
```

<!--


```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
InGff=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain*/*_ORF_RxLR_regex_merged.gff)
InGff_features=$(echo $InGff | sed 's/.gff/_features.gff/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/add_ORF_features.pl $InGff $Assembly >> $InGff_features
for BamFile in $(ls alignment/star/P.cactorum/10300/treatment/genbank/star_aligmentAligned.sortedByCoord.out.bam); do
# OutDir=alignment/star/$Organism/$Strain/featureCounts_genbank
OutDir=$(dirname $BamFile)
mkdir -p $OutDir
# Prefix=$(echo $BamFile | rev | cut -f3 -d '/' | rev | sed 's/vesca_//g')
Prefix=$Strain
echo $Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
# echo "$BamFile $Gff $OutDir $Prefix"
qsub $ProgDir/sub_featureCounts.sh $BamFile $InGff_features $OutDir $Prefix
done
done
```

Those RxLRs with evidence of expression in planta (fpkm >5) were included in
ORF RxLR gene models:

```bash
 Organism="P.cactorum"
 Strain="10300"
 OutDir=$(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain)
 for FeatureCounts in $(ls alignment/star/$Organism/$Strain/treatment/genbank/*_featurecounts.txt); do
   OutFile=$(echo $FeatureCounts | sed 's/.txt/_fpkm_>5.txt/g')
 # cat $FeatureCounts | tail -n+3 | cut -f1,7 | grep -v -e '\s0$' -e '\s1$' -e '\s2$' -e '\s3$' -e '\s4$' > $OutFile
   cat $FeatureCounts | tail -n+3 | awk '$7>5 {print $0}' > $OutFile
 cat $OutFile
 done | cut -f1 | sort | uniq > $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers.txt
 cat $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers.txt | wc -l
 SigP_Gff=$(ls gene_pred/combined_sigP_ORF/P.cactorum/10300/10300_all_secreted_unmerged.gff)
 cat $SigP_Gff | grep -f $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers.txt | cut -f4 -d '=' | cut -f1 -d ';' > $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers_renamed.txt
 cat $OutDir/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers_renamed.txt | wc -l
```

The number of Aditional RxLRs from those  secreted ORF RxLRs with a fpkm >5 was:

```
 904
``` -->



### E5) From ORF gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search ORFs predicted with atg.pl. These were run with the following commands:


```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/P.*/*/*_all_secreted_merged.aa | grep '414'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_ORF_WY_hmmer.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_ORF_WY_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_ORF_WY_hmmer_headers.txt
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
SigP_Merged_Gff=$(echo $Secretome | sed 's/.aa/.gff/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_WY_hmmer.gff
done
```

P.cactorum 414
Initial search space (Z):              24058  [actual number of targets]
Domain search space  (domZ):             111  [number of targets reported over threshold]


### E6) From ORF gene models - Hmm evidence of RxLR effectors

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep '414'); do
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
mkdir -p $OutDir
HmmResults="$Strain"_ORF_RxLR_hmmer_unmerged.txt
hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
echo "$Organism $Strain"
cat $OutDir/$HmmResults | grep 'Initial search space'
cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
Headers="$Strain"_ORF_RxLR_hmmer_headers_unmerged.txt
cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
SigP_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_RxLR_hmmer_unmerged.gff3
RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.gff
RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.txt
RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_hmmer_unmerged.gff3 --db sigP_ORF_RxLR_hmm.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR_hmm.db --id sigP_ORF_RxLR_hmm --out sigP_ORF_RxLR_hmm_merged.db --gff > $RxLR_Merged_Gff
cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
printf "Merged RxLR-EER Hmm proteins:\t"
cat $RxLR_Merged_AA | grep '>' | wc -l
done
```

```
P.cactorum 414
Initial search space (Z):              66113  [actual number of targets]
Domain search space  (domZ):             481  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:	150
```

### E7) Combining RxLRs from Regex and hmm searches


The total ORF RxLRs are

```bash
for RegexRxLREER in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_ORF_RxLR_EER_regex_merged.txt | grep '414'); do
Organism=$(echo $RegexRxLREER | rev |  cut -d '/' -f3 | rev)
Strain=$(echo $RegexRxLREER | rev | cut -d '/' -f2 | rev)
Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF.gff3)
Proteome=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
# RegexRxLRfpkm=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Organism/$Strain/"$Strain"_ORF_RxLR_regex_merged_fpkm_5_headers_renamed.txt)
HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_ORF_RxLR_hmm_merged.txt)
echo "$Organism - $Strain"
echo "Number of RxLR EERs identified by Regex:"
cat $RegexRxLREER | sort | uniq | wc -l
# echo "Number of RxLRs identified by Regex with fpkm > 5"
# cat $RegexRxLRfpkm | sort | uniq | wc -l
echo "Number of RxLRs identified by Hmm:"
cat $HmmRxLR | sort | uniq | wc -l
echo "Number of RxLRs in combined dataset:"
# cat $RegexRxLREER $HmmRxLR $RegexRxLRfpkm | sort | uniq | wc -l
cat $RegexRxLREER $HmmRxLR | sort | uniq | wc -l
echo ""
# echo "Extracting RxLRs from datasets"
OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
mkdir -p $OutDir
# cat $RegexRxLREER $RegexRxLRfpkm $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
cat $RegexRxLREER $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_total_ORF_RxLR_headers.txt $Gff ORF_RxLR Name Augustus > $OutDir/"$Strain"_total_ORF_RxLR.gff
echo "Number of genes in the extracted gff file:"
cat $OutDir/"$Strain"_total_ORF_RxLR.gff | grep -w 'gene' | wc -l
done
```

```
P.cactorum - 414
Number of RxLR EERs identified by Regex:
196
Number of RxLRs identified by Hmm:
150
Number of RxLRs in combined dataset:
217

Number of genes in the extracted gff file:
217
```

## 4.2.c Analysis of RxLR effectors - merger of Augustus / published genes with ORFs

Intersection between the coodinates of putative RxLRs from gene models and ORFs
were identified to determine the total number of RxLRs predicted in these
genomes.

The RxLR effectors from both Gene models and ORF finding approaches were
combined into a single file.

This step was complicated by the inconsistency in downloaded gff files for gene
models.


```bash
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/*/* | grep '414'); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$MergeDir/"$Strain"_total_RxLR.gff
AugTxt=$MergeDir/"$Strain"_total_RxLR_headers.txt
AugFa=$(ls gene_pred/final/"$Species"/"$Strain"/final/final_genes_appended_renamed.pep.fasta)

ORFGff=$(ls $MergeDir/"$Strain"_total_ORF_RxLR.gff)
ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
ORFsTxt=$(ls $MergeDir/"$Strain"_total_ORF_RxLR_headers.txt)

ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_RxLR_EER_motif_hmm.gff
AugInORFs=$MergeDir/"$Strain"_AugInORFs_RxLR_EER_motif_hmm.gff
ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_RxLR_EER_motif_hmm.gff
AugUniq=$MergeDir/"$Strain"_Aug_Uniq_RxLR_EER_motif_hmm.gff
TotalRxLRsTxt=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.txt
TotalRxLRsGff=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.gff

bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq

echo "$Species - $Strain"
echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
cat $ORFsInAug | grep -w 'gene' | wc -l
echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
cat $AugInORFs | grep -w 'gene' | wc -l
echo "The number of RxLRs unique to ORF models:"
cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l
# cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
echo "The number of RxLRs unique to Augustus models:"
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
echo "The total number of putative RxLRs are:"
cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalRxLRsTxt
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
cat $ORFsUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f3 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
cat $TotalRxLRsTxt | wc -l
cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalRxLRsTxt > $TotalRxLRsGff

RxLRsFa=$MergeDir/"$Strain"_final_RxLR_EER.fa
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
# $ProgDir/unwrap_fasta.py --inp_fasta $AugFa | grep -A1 -w -f $AugTxt | grep -v -E '^--$' > $RxLRsFa
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalRxLRsTxt > $RxLRsFa
# echo "$Strain"
$ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalRxLRsTxt >> $RxLRsFa
# echo "$Strain"
echo "The number of sequences extracted is"
cat $RxLRsFa | grep '>' | wc -l
done
```

```
  P.cactorum - 414
  The number of ORF RxLRs overlapping Augustus RxLRs:
  150
  The number of Augustus RxLRs overlapping ORF RxLRs:
  150
  The number of RxLRs unique to ORF models:
  67
  The number of RxLRs unique to Augustus models:
  20
  The total number of putative RxLRs are:
  237
  The number of sequences extracted is
  237
```



### H) From ORF gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers
in ORF gene models. This was done with the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/*/*/*.aa_cat.fa | grep '414'); do
# Setting variables
Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
mkdir -p $OutDir
# Hmmer variables
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
# Searches for LFLAK domain
LFLAK_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm
HmmResultsLFLAK="$Strain"_ORF_CRN_LFLAK_unmerged_hmmer.txt
hmmsearch -T 0 $LFLAK_hmm $Proteome > $OutDir/$HmmResultsLFLAK
echo "Searching for LFLAK domains in: $Organism $Strain"
cat $OutDir/$HmmResultsLFLAK | grep 'Initial search space'
cat $OutDir/$HmmResultsLFLAK | grep 'number of targets reported over threshold'
HmmFastaLFLAK="$Strain"_ORF_CRN_LFLAK_unmerged_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResultsLFLAK $Proteome > $OutDir/$HmmFastaLFLAK
# Searches for DWL domain
DWL_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm
HmmResultsDWL="$Strain"_ORF_CRN_DWL_unmerged_hmmer.txt
hmmsearch -T 0 $DWL_hmm $Proteome > $OutDir/$HmmResultsDWL
echo "Searching for DWL domains in: $Organism $Strain"
cat $OutDir/$HmmResultsDWL | grep 'Initial search space'
cat $OutDir/$HmmResultsDWL | grep 'number of targets reported over threshold'
HmmFastaDWL="$Strain"_ORF_CRN_DWL_unmerged_hmmer.fa
$ProgDir/hmmer2fasta.pl $OutDir/$HmmResultsDWL $Proteome > $OutDir/$HmmFastaDWL
# Identify ORFs found by both models
CommonHeaders=$OutDir/"$Strain"_ORF_CRN_DWL_LFLAK_unmerged_headers.txt
cat $OutDir/$HmmFastaLFLAK $OutDir/$HmmFastaDWL | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $CommonHeaders
echo "The number of CRNs common to both models are:"
cat $CommonHeaders | wc -l
# The sequences will be merged based upon the strength of their DWL domain score
# For this reason headers as they appear in the DWL fasta file were extracted
Headers=$OutDir/"$Strain"_CRN_hmmer_unmerged_headers.txt
cat $OutDir/$HmmFastaDWL | grep '>' | grep -w -f $CommonHeaders | tr -d '>' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | tr -d '-' | sed 's/hmm_score/HMM_score/g' > $Headers
# As we are dealing with JGI and Broad sequences, some features need formatting:
ORF_Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/*_ORF.gff3 | grep -v '_atg_')
# Gff features were extracted for each header
CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
# cat $Headers | cut -f1 > tmp.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $Headers $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
# $ProgDir/extract_gff_for_sigP_hits.pl tmp.txt $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
# Gff features were merged based upon the DWL hmm score
DbDir=analysis/databases/$Organism/$Strain
mkdir -p $DbDir
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $CRN_unmerged_Gff --db $DbDir/CRN_ORF.db
CRN_Merged_Gff=$OutDir/"$Strain"_CRN_merged_hmmer.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp $DbDir/CRN_ORF.db --id LFLAK_DWL_CRN --out $DbDir/CRN_ORF_merged.db --gff > $CRN_Merged_Gff
# Final results are reported:
echo "Number of CRN ORFs after merging:"
cat $CRN_Merged_Gff | grep 'gene' | wc -l
done
```

```
  Searching for LFLAK domains in: P.cactorum 414
  Initial search space (Z):             552687  [actual number of targets]
  Domain search space  (domZ):             337  [number of targets reported over threshold]
  Searching for DWL domains in: P.cactorum 414
  Initial search space (Z):             552687  [actual number of targets]
  Domain search space  (domZ):             459  [number of targets reported over threshold]
  The number of CRNs common to both models are:
  223
  Number of CRN ORFs after merging:
  153
```

Extract crinklers from ORFs and Braker/codingquary gene models


```bash
for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/*/* | grep '414'); do
Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
AugGff=$(ls $MergeDir/"$Strain"_pub_CRN_LFLAK_DWL.gff)
AugFa=$(ls gene_pred/final/"$Species"/"$Strain"/final/final_genes_appended_renamed.pep.fasta)
ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer.gff3
ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_CRN_hmmer.bed
AugInORFs=$MergeDir/"$Strain"_AugInORFs_CRN_hmmer.bed
ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_CRN_hmmer.bed
AugUniq=$MergeDir/"$Strain"_Aug_Uniq_CRN_hmmer.bed
TotalCRNsTxt=$MergeDir/"$Strain"_final_CRN.txt
TotalCRNsGff=$MergeDir/"$Strain"_final_CRN.gff
TotalCRNsHeaders=$MergeDir/"$Strain"_Total_CRN_headers.txt
bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq
echo "$Species - $Strain"

echo "The number of ORF CRNs overlapping Augustus CRNs:"
cat $ORFsInAug | grep -w -e 'transcript' -e 'mRNA' | wc -l
echo "The number of Augustus CRNs overlapping ORF CRNs:"
cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA' | wc -l
cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalCRNsTxt
echo "The number of CRNs unique to ORF models:"
cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' | wc -l
cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt
echo "The number of CRNs unique to Augustus models:"
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt

cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalCRNsTxt > $TotalCRNsGff

CRNsFa=$MergeDir/"$Strain"_final_CRN.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalCRNsTxt > $CRNsFa
$ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalCRNsTxt >> $CRNsFa
echo "The number of sequences extracted is"
cat $CRNsFa | grep '>' | wc -l
done
```

```
  P.cactorum - 414
  The number of ORF CRNs overlapping Augustus CRNs:
  118
  The number of Augustus CRNs overlapping ORF CRNs:
  120
  The number of CRNs unique to ORF models:
  35
  The number of CRNs unique to Augustus models:
  3
  The number of sequences extracted is
  158
```



## C.ii) SSCP

Small secreted cysteine rich proteins were identified within secretomes. These
proteins may be identified by EffectorP, but this approach allows direct control
over what constitutes a SSCP.

```bash

for Secretome in $(ls gene_pred/combined_sigP/*/*/*_secreted.fa | grep -v 'all_secreted' |  grep '414'); do
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/sscp/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_all_results.fa
cat $OutDir/"$Strain"_sscp_all_results.fa | grep 'Yes' > $OutDir/"$Strain"_sscp.fa
printf "number of SSC-rich genes:\t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' | cut -f1 -d '.' | sort | uniq | wc -l
printf "Number of effectors predicted by EffectorP:\t"
EffectorP=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP_secreted_headers.txt)
cat $EffectorP | wc -l
printf "Number of SSCPs predicted by both effectorP and this approach: \t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' > $OutDir/"$Strain"_sscp_headers.txt
cat $OutDir/"$Strain"_sscp_headers.txt $EffectorP | cut -f1 | sort | uniq -d | wc -l
echo ""
done
```

```
P.cactorum - 414
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	402
number of SSC-rich genes:	401
Number of effectors predicted by EffectorP:	925
Number of SSCPs predicted by both effectorP and this approach: 	235
```
<!--
```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa | grep '10300'); do
ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"

OutDir=analysis/sscp_ORF/$Organism/$Strain
mkdir -p $OutDir
printf "the number of SigP gene is:\t";
cat $Secretome | grep '>' | cut -f1 | sort | uniq | wc -l
printf "the number of SSCP genes are:\t";
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_ORF.fa

cat $OutDir/"$Strain"_sscp_ORF.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/"$Strain"_sscp_ORF_headers_unmerged.txt
cat $OutDir/"$Strain"_sscp_ORF_headers_unmerged.txt | tr -d ' ' | sort | uniq | wc -l

SigP_Gff=$(ls gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_sscp_ORF_headers_unmerged.txt  $SigP_Gff	sscp_filter.py Name Augustus > $OutDir/"$Strain"_sscp_ORF_unmerged.gff

SSCP_Merged_Gff=$OutDir/"$Strain"_ORF_sscp_merged.gff
SSCP_Merged_txt=$OutDir/"$Strain"_ORF_sscp_merged.txt
SSCP_Merged_AA=$OutDir/"$Strain"_ORF_sscp_merged.aa
ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
$ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_sscp_ORF_unmerged.gff --db sigP_ORF_sscp.db
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_sscp.db --id sigP_ORF_sscp --out sigP_ORF_sscp_merged.db --gff > $SSCP_Merged_Gff
cat $SSCP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $SSCP_Merged_txt

ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
$ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $SSCP_Merged_txt > $SSCP_Merged_AA
printf "Merged SSCP proteins:\t"
cat $SSCP_Merged_AA | grep '>' | wc -l
printf "\n"
done
```
-->



# Making a combined file of Braker, Coding quary genes as well as additional ORF effector candidates


A gff file containing the combined Braker and CodingQuary genes as well as the
additional CRN and RxLR genes predicted by ORF analysis was made.

```bash
for GeneGff in $(ls gene_pred/final/P*/*/final/final_genes_appended_renamed.gff3 | grep '414'); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
GffOrfRxLR=$(ls analysis/RxLR_effectors/combined_evidence/$Organism/$Strain/*_ORFsUniq_RxLR_EER_motif_hmm.gff)
GffOrfCRN=$(ls analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/*_ORFsUniq_CRN_hmmer.bed)
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked_wrapped.fa)
OutDir=gene_pred/final_incl_ORF/$Organism/$Strain
mkdir -p $OutDir
cat $GeneGff > $OutDir/${Strain}_genes_incl_ORFeffectors.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/add_ORF_features.pl $GffOrfRxLR $Assembly >> $OutDir/${Strain}_genes_incl_ORFeffectors.gff3
$ProgDir/add_ORF_features.pl $GffOrfCRN $Assembly >> $OutDir/${Strain}_genes_incl_ORFeffectors.gff3
# Make gene models from gff files.
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked_wrapped.fa)
$ProgDir/gff2fasta.pl $Assembly $OutDir/${Strain}_genes_incl_ORFeffectors.gff3 $OutDir/${Strain}_genes_incl_ORFeffectors
done
```

<!-- High confidence genes overlapped by ORFs to be inserted were manually inspected
to see whether the original gene model or the ORF should be kept. Preference was
given to the ORF effector candidate. Before this could be done bedtools was
used to identify which genes were intersected: -->

High confidence genes overlapped by ORFs to be inserted were manually inspected
to see whether the original gene model or the ORF should be kept. Manual inspection.
indicated that predicted gene models were, in general, represented the true gene.
As this approach was to be used for 20 genomes, not all conflicts were assessed
manually. Bedtools was used to identify which genes were intersected:

```bash
# GeneGff=$(ls gene_pred/final/P.cactorum/10300/final/final_genes_appended.gff3)
for GeneGff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3); do
  Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=gene_pred/final_incl_ORF/$Organism/$Strain
  GffOrfRxLR=$(ls analysis/RxLR_effectors/combined_evidence/$Organism/$Strain/${Strain}_ORFsUniq_RxLR_EER_motif_hmm.gff)
  GffOrfCRN=$(ls analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain/${Strain}_ORFsUniq_CRN_hmmer.bed)
  Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked_wrapped.fa)
  # Identify RxLR ORFs intersecting non-RxLR gene models
  bedtools intersect -wo -a $GeneGff -b $GffOrfRxLR | grep -e "AUGUSTUS.gene" | grep "ORF_RxLR.gene"| cut -f1,9,18 | sed 's/ID=//g' | tr -d ';' > $OutDir/RxLR_ORFs_intersecting_non-RxLR_genes.txt
  # Identify CRN ORFs intersecting non-CRN gene models
  bedtools intersect -wo -a $GeneGff -b $GffOrfCRN | grep -e "AUGUSTUS.gene" | grep "CRN_HMM.gene"| cut -f1,9,18 | sed 's/ID=//g' | tr -d ';' > $OutDir/CRN_ORFs_intersecting_non-CRN_genes.txt
done

```

```
RxLR ORFs intersecting non-RxLR gene models
contig_1	g482	ORF_F_4497
contig_3	g1581	ORF_R_20817
contig_7	g4433	ORF_R_40635
contig_8	g4929	ORF_R_45860
contig_11	g6412	ORF_F_59902
contig_15	g8356	ORF_F_77083
contig_19	g9801	ORF_R_91763
contig_21	g10605	ORF_F_98304
contig_21	g10731	ORF_F_99312
contig_22	g10924	ORF_R_102635
contig_22	g10932	ORF_R_102487
contig_22	g10934	ORF_F_101323
contig_24	g11734	ORF_F_109274
contig_26	g12148	ORF_F_113182
contig_28	g12824	ORF_R_120577
contig_30	g13381	ORF_R_126202
contig_32	g13844	ORF_F_129292
contig_33	g14172	ORF_F_132408
contig_37	g15429	ORF_R_145005
contig_42	g16718	ORF_F_155678
contig_45	g17601	ORF_R_163887
contig_45	g17602	ORF_R_163887
contig_49	g18521	ORF_F_172351
contig_50	g18819	ORF_R_174418
contig_56	g20001	ORF_R_188375
contig_67	g22241	ORF_F_206962
contig_67	g22247	ORF_R_207819
contig_71	g22914	ORF_R_213951
contig_73	g23139	ORF_F_215655
contig_81	g24333	ORF_F_227132
contig_82	g24502	ORF_R_229202
contig_85	g24836	ORF_R_233563
contig_89	g25437	ORF_F_237194
contig_92	g25699	ORF_R_241479
contig_93	g25832	ORF_F_241018
contig_96	g26073	ORF_F_243324
contig_107	g27054	ORF_F_251751
contig_109	g27211	ORF_F_252862
contig_111	g27413	ORF_R_256682
contig_132	g28518	ORF_F_264890
contig_136	g28623	ORF_F_265849

CRN ORFs intersecting non-CRN gene models
contig_1	g180	ORF_R_5744
contig_6	g3550	ORF_R_37812
contig_9	g5161	ORF_R_53282
contig_9	g5438	ORF_R_50977
contig_12	g6859	ORF_R_63534
contig_12	g7052	ORF_F_65649
contig_15	g8156	ORF_F_75118
contig_16	g8748	ORF_R_78832
contig_23	g11215	ORF_F_104161
contig_25	g11846	ORF_R_112513
contig_36	g15092	ORF_F_141132
contig_46	g17850	ORF_R_166178
contig_49	g18550	ORF_F_172548
contig_57	g20253	ORF_F_188445
contig_60	g20951	ORF_R_195366
contig_64	g21705	ORF_F_201259
contig_66	g22126	ORF_F_205962
contig_74	g23413	ORF_F_218328
contig_78	g23977	ORF_F_223976
contig_88	g25210	ORF_R_237592
contig_88	g25213	ORF_R_237550
contig_94	g25873	ORF_R_243407
contig_99	g26391	ORF_R_247624
contig_112	g27429	ORF_R_257685
contig_114	g27567	ORF_R_258773
contig_122	g28079	ORF_F_260966
contig_122	g28112	ORF_R_262649
contig_144	g28804	ORF_R_269850
contig_145	g28843	ORF_R_270035
contig_146	g28854	ORF_R_270302
contig_158	g29049	ORF_R_272357
```

Some ORF effectors were noted to overlap AugustusCodingQuary gene models. These
were manually inspected in geneious and the following genes removed:

```bash
for CombinedGff in $(ls gene_pred/final_incl_ORF/P.cactorum/414/414_genes_incl_ORFeffectors.gff3); do
OutDir=$(dirname $CombinedGff)
cat $OutDir/RxLR_ORFs_intersecting_non-RxLR_genes.txt $OutDir/CRN_ORFs_intersecting_non-CRN_genes.txt | cut -f3 > $OutDir/exclude_list.txt
echo "genes before filtering"
cat $CombinedGff | grep -w 'gene' | wc -l
echo "genes in list"
cat $OutDir/exclude_list.txt | wc -l
echo "Unique genes in list"
cat $OutDir/exclude_list.txt | sort | uniq | wc -l
echo "Number of genes removed"
cat $CombinedGff | grep -w -f $OutDir/exclude_list.txt | grep -w 'gene' | wc -l
cat $CombinedGff | grep -w -v -f $OutDir/exclude_list.txt > $OutDir/${Strain}_genes_incl_ORFeffectors_filtered.gff3
echo "final number of genes"
cat $OutDir/${Strain}_genes_incl_ORFeffectors_filtered.gff3 | grep -w 'gene' | wc -l
done
```
```
genes before filtering
29629
genes in list
72
Unique genes in list
71
Number of genes removed
71
final number of genes
29558
```

```bash
 for GffAppended in $(ls gene_pred/final_incl_ORF/*/*/*_genes_incl_ORFeffectors_filtered.gff3 | grep '414'); do
   Strain=$(echo $GffAppended | rev | cut -d '/' -f2 | rev)
   Organism=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
   echo "$Organism - $Strain"
   FinalDir=$(dirname $GffAppended)
   GffFiltered=$FinalDir/filtered_duplicates.gff
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
   $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
   GffRenamed=$FinalDir/final_genes_genes_incl_ORFeffectors_renamed.gff3
   LogFile=$FinalDir/final_genes_appended_renamed.log
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
   $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
   rm $GffFiltered
  Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked_wrapped.fa)
 $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_genes_incl_ORFeffectors_renamed
   # The proteins fasta file contains * instead of Xs for stop codons, these should
   # be changed
   sed -i 's/\*/X/g' $FinalDir/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta
 done
```

No duplicate genes were found.



# Assessing gene space in predicted transcriptomes

```bash
for Transcriptome in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.gene.fasta); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Transcriptome | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/eukaryota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/genes
qsub $ProgDir/sub_busco3.sh $Transcriptome $BuscoDB $OutDir
done
```

```
INFO    281 Complete BUSCOs (C)
INFO    272 Complete and single-copy BUSCOs (S)
INFO    9 Complete and duplicated BUSCOs (D)
INFO    7 Fragmented BUSCOs (F)
INFO    15 Missing BUSCOs (M)
INFO    303 Total BUSCO groups searched
```

```bash
 for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt); do
 Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
 Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
 Complete=$(cat $File | grep "(C)" | cut -f2)
 Single=$(cat $File | grep "(S)" | cut -f2)
 Fragmented=$(cat $File | grep "(F)" | cut -f2)
 Missing=$(cat $File | grep "(M)" | cut -f2)
 Total=$(cat $File | grep "Total" | cut -f2)
 echo -e "$Organism\t$Strain\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
 done
```



Make a list of BUSCO genes:

```bash
BuscoHits=$(ls gene_pred/busco/P.cactorum/414_v2/genes/run_final_genes_combined.gene/single_copy_busco_sequences/EOG0937*.fna)
OutDir=$(ls -d gene_pred/busco/P.cactorum/414_v2/genes/run_final_genes_combined.gene)
cat $BuscoHits | grep '>' | cut -f3 -d ':' > $OutDir/busco_single_copy_gene_headers.txt
```


# Re-running Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	# cd /home/groups/harrisonlab/project_files/idris
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta | grep '414'); do
  FileName=$(basename $Genes)
  SymLinkDir=$(dirname $Genes)/symlink_directory
  mkdir $SymLinkDir
  cp -s $PWD/$Genes $SymLinkDir/$FileName
	echo $Genes
	$ProgDir/sub_interproscan.sh $SymLinkDir/$FileName
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteome in $(ls gene_pred/final_incl_ORF/P.*/*/symlink_directory/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteome $InterProRaw
done
```


## B) SwissProt

```bash
 for Proteome in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.pep.fasta | grep '414'); do
   Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   OutDir=gene_pred/swissprot/$Organism/$Strain
   SwissDbDir=../../../../home/groups/harrisonlab/uniprot/swissprot
   SwissDbName=uniprot_sprot
   ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
   qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
 done
```


## RNAseq
<!--
```bash
for Transcriptome in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta | grep '414'); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d alignment/star/fvesca/v1.1/*/* | grep 'mycelium'); do
FileF=$(ls $RNADir/*.mate1.fq.gz)
FileR=$(ls $RNADir/*.mate2.fq.gz)
echo $FileF
echo $FileR
Prefix=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
Timepoint=$(echo $RNADir | rev | cut -f2 -d '/' | rev)
echo "$Timepoint"
OutDir=alignment/salmon/$Organism/"$Strain"/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_salmon.sh $Transcriptome $FileF $FileR $OutDir
done
done
``` -->


```bash
for Transcriptome in $(ls gene_pred/final_incl_ORF/*/*/final_genes_genes_incl_ORFeffectors_renamed.cds.fasta | grep '414'); do
Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_rna/paired/*/* | grep -e 'mycelium' -e '_no_vesca'); do
FileNum=$(ls $RNADir/F/*.mate1.fq.gz | wc -l)
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
for num in $(seq 1 $FileNum); do
printf "\n"
FileF=$(ls $RNADir/F/*.mate1.fq.gz | head -n $num | tail -n1)
FileR=$(ls $RNADir/R/*.mate2.fq.gz | head -n $num | tail -n1)
echo $FileF
echo $FileR
Prefix=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint"
OutDir=alignment/salmon/$Organism/"$Strain"/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_salmon.sh $Transcriptome $FileF $FileR $OutDir
done
done
done
```


```bash
for RNADir in $(ls -d ../../../../home/groups/harrisonlab/project_files/idris/qc_rna/paired/*/* | grep -e 'mycelium' -e '_no_vesca' | grep 'mycelium'); do
FileNum=$(ls $RNADir/F/*.mate1.fq.gz | wc -l)
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
for num in $(seq 1 $FileNum); do
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
done
printf "\n"
FileF=$(ls $RNADir/F/*.mate1.fq.gz | head -n $num | tail -n1)
FileR=$(ls $RNADir/R/*.mate2.fq.gz | head -n $num | tail -n1)
echo $FileF
echo $FileR
# Prefix=$(echo $FileF | rev | cut -f2 -d '/' | rev | sed "s/_R.*_trim.fq.gz//g")
Prefix=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint"
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
done
done
```
