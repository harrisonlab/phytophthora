# Submission Commands

Submisison of annotations with an assembly appears to be a complex process.
If a genome is to be submitted without annotation then all that is needed is the
fasta file containing the assembled contigs. If an annotated genome is to be
submitted then a number of processing steps are required before submission. The
fasta file of contigs and the gff file of annotations must be combined to form a
.asn file. The program that does this conversion (tbl2asn) requires the fasta
files and gff files to be formatted correctly. In the case of the gff file, this
means parsing it to a .tbl file.

The commands used to parse these files and prepare the F. oxysporum f. sp.
cepae genome for submisson are shown below.


# Preliminary submission

A Bioproject and biosample number was prepared for the genome submission at:
https://submit.ncbi.nlm.nih.gov

A preliminary submission was made for the .fasta assembly to check if
any contigs needed to be split. This step was performed early in the annotation
process (prior to gene prediction) to ensure that annotation did not have to
be repeated at the end of the project.


The following note was provided in the WGS submission page on NCBI in the box
labeled "Private comments to NCBI staff":

```
I have been advised to submit my assemblies to NCBI early in my submission process to ensure that my contigs pass the contamination screen. This assembly will be revised as appropriate, including renaming of contigs where needed. Please allow me to modify this submission at a later date, including upload of the final gene models.

'For future submissions, you could send us the fasta files early
in the submission process so we can run them through our foreign
contamination screen. We will let you know if we find any
sequences to exclude or trim before you generate your final
WGS submission.'...'*IMPORTANT* Include a comment that you are submitting
the fasta files to be screened by the contamination screen
prior to creating your final annotated submission.'
```


# Submission of sequence data to SRA
<!--
Reads were submitted to the SRA at https://submit.ncbi.nlm.nih.gov/subs/sra/ .
To do this, a metadata file was provided detailing each of the files in the
bioproject. The file was downloaded in excel format and edited manually. A copy
of the edited file and the final .tsv file is present at:

```bash
  # For genmic reads:
  ls genome_submission/SRA_metadata_acc.txt genome_submission/SRA_metadata_acc.xlsx
  # For RNAseq reads:
  ls genome_submission/SRA_metadata_acc.txt genome_submission/RNAseq_SRA_metadata_acc.txt
```

As these files included a file > 500 Mb, a presubmission folder was requested.
This aids submission of large data files. This file was created on the ftp server
at ftp-private.ncbi.nlm.nih.gov, with a private folder named
uploads/andrew.armitage@emr.ac.uk_6L2oakBI. Ncbi provided a username a password.
Files were uploaded into a folder created within my preload folder using ftp.

For genomic reads:
```bash
  # Bioproject="PRJNA338236"
  SubFolder="FoC_PRJNA338256"
  mkdir $SubFolder
  for Read in $(ls raw_dna/paired/F.*/*/*/*.fastq.gz | grep -w -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e '414' | grep -v s_6_*_sequence.fastq.gz); do
    echo $Read;
    cp $Read $SubFolder/.
  done
  cp raw_dna/pacbio/P.cactorum/414/extracted/concatenated_pacbio.fastq $SubFolder/.
  cd $SubFolder
  gzip concatenated_pacbio.fastq
  ftp ftp-private.ncbi.nlm.nih.gov
  cd uploads/andrew.armitage@emr.ac.uk_6L2oakBI
  mkdir FoC_PRJNA338256_2
  cd FoC_PRJNA338256_2
  # put FoN_PRJNA338236
  prompt
  mput *
  bye
  cd ../
  rm -r $SubFolder
```

For RNAseq Reads:
```bash
  SubFolder="FoC_RNAseq_PRJNA338256"
  mkdir $SubFolder
  for Read in $(ls qc_rna/paired/F.*/*/*/*_trim.fq.gz); do
    echo $Read;
    cp $Read $SubFolder/.
  done
  cd $SubFolder
  ftp ftp-private.ncbi.nlm.nih.gov
  cd uploads/andrew.armitage@emr.ac.uk_6L2oakBI
  mkdir FoC_RNAseq_PRJNA338256
  cd FoC_RNAseq_PRJNA338256
  prompt
  mput *
  bye
  cd ../
  rm -r $SubFolder
``` -->

## Making a table for locus tags:

locus tags were provided by ncbi when the bioproject was registered.

Upon request the bioproject team chged the locus tag to "Pcac1"
<!--
A table detailing their relationship to the strain was made manually. This could
be searched later to extract locus tags for particular strains.

```bash
mkdir -p genome_submission/
printf \
"BFJ65 SAMN05529097 FoC_414
BFJ66 SAMN05529098 FoC_125
BFJ67 SAMN05529099 FoC_A23
BFJ68 SAMN05529100 Fo_A28
BFJ69 SAMN05529101 Fo_A13
BFJ70 SAMN05529102 Fo_PG
BFJ71 SAMN05529103 Fo_CB3
BFJ72 SAMN05529104 Fp_A8" \
> genome_submission/FoC_PRJNA338256_locus_tags.txt
``` -->

# Final Submission

These commands were used in the final submission of the Pcac P414 genome:


## Output directory
An output and working directory was made for genome submission:

```bash
  ProjDir=/data/scratch/armita/idris/
  cd $ProjDir
  OutDir="genome_submission/P.cactorum/414"
  mkdir -p $OutDir
```

## SbtFile
The genbank submission template tool was used at:
http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
This produce a template file '.sbt' detailing the submission.

## Setting varibales
Vairables containing locations of files and options for scripts were set:

```bash
# Program locations:
AnnieDir="/home/armita/prog/annie/genomeannotation-annie-c1e848b"
ProgDir="/home/armita/git_repos/emr_repos/tools/genbank_submission"
# File locations:
SbtFile="/home/groups/harrisonlab/project_files/idris/genome_submission/P.cactorum/414_v2/template.sbt"
Assembly=$(ls repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
InterProTab=$(ls gene_pred/interproscan/P.cactorum/414/414_interproscan.tsv)
SwissProtBlast=$(ls gene_pred/swissprot/P.cactorum/414/swissprot_vJul2016_tophit_parsed.tbl)
SwissProtFasta=$(ls /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot.fasta)
GffFile=$(ls gene_pred/final_incl_ORF/P.cactorum/414/final_genes_genes_incl_ORFeffectors_renamed.gff3)
# tbl2asn options:
Organism="Phytophthora cactorum"
Strain="P414"
# ncbi_tbl_corrector script options:
SubmissionID="Pcac1"
LabID="NIAB-EMR"
# GeneSource='ab initio prediction:Braker:1.9, CodingQuary:2.0'
# IDSource='similar to AA sequence:SwissProt:2016_07'
# IDSource='similar to AA sequence:UniProtKB/Swiss-Prot'
# Final submisison file name:
FinalName="Pcac_P414_2018"
```

<!-- ## Preparing Gff input file

Parse the Augustus Gff file.
Transcripts should be renamed as mRNA features. Exons should be added to the
Gff and unique IDs should be given to all features in the file.

```bash
cat $GffFile | sed 's/transcript/mRNA/g' > $OutDir/GffMRNA.gff
# $ProgDir/generate_tbl_file/exon_generator.pl $OutDir/GffMRNA.gff > $OutDir/corrected_exons.gff
# $ProgDir/generate_tbl_file/gff_add_id.py --inp_gff $OutDir/corrected_exons.gff --out_gff $OutDir/corrected_exons_id.gff
``` -->

## Generating .tbl file (GAG)

The Genome Annotation Generator (GAG.py) can be used to convert gff files into
.tbl format, for use by tbl2asn.

It can also add annotations to features as provided by Annie the Annotation
extractor.

### Extracting annotations (Annie)

Interproscan and Swissprot annotations were extracted using annie, the
ANNotation Information Extractor. The output of Annie was filtered to
keep only annotations with references to ncbi approved databases.
Note - It is important that transcripts have been re-labelled as mRNA by this
point.

```bash
  export PYTHONPATH="/home/armita/.local/lib/python3.5/site-packages"
  python3 $AnnieDir/annie.py -ipr $InterProTab -g $GffFile -b $SwissProtBlast -db $SwissProtFasta -o $OutDir/annie_output.csv --fix_bad_products
  $ProgDir/edit_tbl_file/annie_corrector.py --inp_csv $OutDir/annie_output.csv --out_csv $OutDir/annie_corrected_output.csv
  # re-load your python path variable
  . ~/.profile
```

### Running GAG

Gag was run using the modified gff file as well as the annie annotation file.
Gag was noted to output database references incorrectly, so these were modified.

```bash
mkdir -p $OutDir/gag/round1
gag.py -f $Assembly -g $GffFile -a $OutDir/annie_corrected_output.csv --fix_start_stop -o $OutDir/gag/round1 2>&1 | tee $OutDir/gag_log1.txt
sed -i 's/Dbxref/db_xref/g' $OutDir/gag/round1/genome.tbl
```

## manual edits

<!-- The gene NS_04463 was found to use the same start and stop codon as predicted
gene CUFF_4598_1_205. Both of these genes were predicted by codingquary. Neither
of these genes were predicted as having alternative splicing. As such the gene
NS_04463 was removed. The same was found for genes CUFF_11067_2_85 and
CUFF_11065_1_82 and as a result CUFF_11067_2_85 was removed.

```bash
  nano $OutDir/gag/round1/genome.tbl
``` -->

## tbl2asn round 1

tbl2asn was run an initial time to collect error reports on the current
formatting of the .tbl file.
Note - all input files for tbl2asn need to be in the same directory and have the
same basename.

```bash
cp $Assembly $OutDir/gag/round1/genome.fsa  
cp $SbtFile $OutDir/gag/round1/genome.sbt
mkdir -p $OutDir/tbl2asn/round1
tbl2asn -p $OutDir/gag/round1/. -t $OutDir/gag/round1/genome.sbt -r $OutDir/tbl2asn/round1 -M n -X E -Z $OutDir/gag/round1/discrep.txt -j "[organism=$Organism] [strain=$Strain]"
```

## Editing .tbl file

The tbl2asn .val output files were observed and errors corrected. This was done
with an in house script. The .val file indicated that some cds had premature
stops, so these were marked as pseudogenes ('pseudo' - SEQ_FEAT.InternalStop)
and that some genes had cds coordinates that did not match the end of the gene
if the protein was hanging off a contig ('stop' - SEQ_FEAT.NoStop).
Furthermore a number of other edits were made to bring the .tbl file in line
with ncbi guidelines. This included: Marking the source of gene
predictions and annotations ('add_inference'); Correcting locus_tags to use the
given ncbi_id ('locus_tag'); Correcting the protein and transcript_ids to
include the locus_tag and reference to submitter/lab id ('lab_id'), removal of
annotated names of genes if you don't have high confidence in their validity
(--gene_id 'remove'). If 5'-UTR and 3'-UTR were not predicted during gene
annotation then genes, mRNA and exon features need to reflect this by marking
them as incomplete ('unknown_UTR').
<!--
```bash
  mkdir -p $OutDir/gag/edited
  $ProgDir/edit_tbl_file/ncbi_tbl_corrector.py --inp_tbl $OutDir/gag/round1/genome.tbl --inp_val $OutDir/tbl2asn/round1/genome.val --locus_tag $SubmissionID --lab_id $LabID --gene_id "remove" --add_inference "$IDSource" --edits stop pseudo unknown_UTR --out_tbl $OutDir/gag/edited/genome.tbl
  # $ProgDir/edit_tbl_file/ncbi_tbl_corrector.py --inp_tbl $OutDir/gag/round1/genome.tbl --inp_val $OutDir/tbl2asn/round1/genome.val --locus_tag $SubmissionID --lab_id $LabID --gene_id "remove" --add_inference "$GeneSource" "$IDSource" --edits stop pseudo unknown_UTR --out_tbl $OutDir/gag/edited/genome.tbl
``` -->

```bash
  mkdir -p $OutDir/gag/edited
  $ProgDir/edit_tbl_file/ncbi_tbl_corrector.py --inp_tbl $OutDir/gag/round1/genome.tbl --inp_val $OutDir/tbl2asn/round1/genome.val --locus_tag $SubmissionID --lab_id $LabID --gene_id "remove" --edits stop pseudo unknown_UTR correct_partial --remove_product_locus_tags "True" --out_tbl $OutDir/gag/edited/genome.tbl
```


## Generating a structured comment detailing annotation methods

```bash
printf "StructuredCommentPrefix\t##Genome-Annotation-Data-START##
Annotation Provider\tHarrison Lab NIAB-EMR
Annotation Date\tMarch-2018
Annotation Version\tRelease 1.00
Annotation Method\tAb initio gene prediction: Braker 1.9 and CodingQuary 2.0; Functional annotation: Swissprot (July 2016 release) and Interproscan 5.18-57.0" \
> $OutDir/gag/edited/annotation_methods.strcmt.txt

```
## Final run of tbl2asn

Following correction of the GAG .tbl file, tbl2asn was re-run to provide the
final genbank submission file.

The options -l paired-ends -a r10k inform how to handle runs of Ns in the
sequence, these options show that paired-ends have been used to estimate gaps
and that runs of N's longer than 10 bp should be labelled as gaps.

```bash
  cp $Assembly $OutDir/gag/edited/genome.fsa
  cp $SbtFile $OutDir/gag/edited/genome.sbt
  mkdir $OutDir/tbl2asn/final
  tbl2asn -p $OutDir/gag/edited/. -t $OutDir/gag/edited/genome.sbt -r $OutDir/tbl2asn/final -M n -X E -Z $OutDir/tbl2asn/final/discrep.txt -j "[organism=$Organism] [strain=$Strain]" -l paired-ends -a r10k -w $OutDir/gag/edited/annotation_methods.strcmt.txt
  cat $OutDir/tbl2asn/final/genome.sqn | sed 's/_pilon//g' | sed 's/\. subunit/kDa subunit/g' | sed 's/, mitochondrial//g' > $OutDir/tbl2asn/final/$FinalName.sqn
```
