# Submission Commands

Submisison of annotations with an assembly appears to be a complex process.
If a genome is to be submitted without annotation then all that is needed is the
fasta file containing the assembled contigs. If an annotated genome is to be
submitted then a number of processing steps are required before submission. The
fasta file of contigs and the gff file of annotations must be combined to form a
.asn file. The program that does this conversion (tbl2asn) requires the fasta
files and gff files to be formatted correctly. In the case of the gff file, this
means parsing it to a .tbl file.

The commands used to parse these files and prepare the Alternaria spp. genomes
for submisson are shown below.


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

## Making a table for locus tags:

locus tags were provided by ncbi when the bioproject was registered.

A table detailing their relationship to the strain was made manually. This could
be searched later to extract locus tags for particular strains.

```bash
mkdir -p genome_submission/
printf \
"PC110 SUB1952873 10300
" \
> genome_submission/Pcac_PRJNA343457_locus_tags.txt
```

# Final Submission

These commands were used in the final submission of Alternaria spp. genomes:


## Output directory
An output and working directory was made for genome submission:

```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
  echo "$Organism - $Strain"
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  OutDir="genome_submission/$Organism/$Strain"
  mkdir -p $OutDir
done
```

## SbtFile
The genbank submission template tool was used at:
http://www.ncbi.nlm.nih.gov/WebSub/template.cgi
This produce a template file detailing the submission.

## Setting varibales
Vairables containing locations of files and options for scripts were set:

```bash
# Program locations:
AnnieDir="/home/armita/prog/annie/genomeannotation-annie-c1e848b"
ProgDir="/home/armita/git_repos/emr_repos/tools/genbank_submission"
# File locations:
SbtFile="genome_submission/P.cactorum/10300/template.sbt"
LabID="ArmitageEMR"
```

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
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  OutDir="genome_submission/$Organism/$Strain"
  GffFile=$(ls gene_pred/final_incl_ORF/$Organism/"$Strain"/final_genes_genes_incl_ORFeffectors_renamed.gff3)

  InterProTab=$(ls gene_pred/interproscan/$Organism/"$Strain"*/"$Strain"*_interproscan.tsv)
  SwissProtBlast=$(ls gene_pred/swissprot/$Organism/"$Strain"*/swissprot_vJul2016_tophit_parsed.tbl)
  SwissProtFasta=$(ls /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot.fasta)
  python3 $AnnieDir/annie.py -ipr $InterProTab -g $GffFile -b $SwissProtBlast -db $SwissProtFasta -o $OutDir/annie_output.csv --fix_bad_products
  ProgDir=/home/armita/git_repos/emr_repos/tools/genbank_submission
  $ProgDir/edit_tbl_file/annie_corrector.py --inp_csv $OutDir/annie_output.csv --out_csv $OutDir/annie_corrected_output.csv
done
```

### Running GAG

Gag was run using the modified gff file as well as the annie annotation file.
Gag was noted to output database references incorrectly, so these were modified.

```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir="genome_submission/$Organism/$Strain"
GffFile=$(ls gene_pred/final_incl_ORF/$Organism/"$Strain"/final_genes_genes_incl_ORFeffectors_renamed.gff3)
mkdir -p $OutDir/gag/round1
gag.py -f $Assembly -g $GffFile -a $OutDir/annie_corrected_output.csv --fix_start_stop -o $OutDir/gag/round1 2>&1 | tee $OutDir/gag_log1.txt
sed -i 's/Dbxref/db_xref/g' $OutDir/gag/round1/genome.tbl
done
```

<!-- ## manual edits

The gene NS_04463 was found to use the same start and stop codon as predicted
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
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir="genome_submission/$Organism/$Strain"

cp $Assembly $OutDir/gag/round1/genome.fsa
SbtFile="genome_submission/P.cactorum/10300/template.sbt"
cp $SbtFile $OutDir/gag/round1/genome.sbt
mkdir -p $OutDir/tbl2asn/round1
tbl2asn -p $OutDir/gag/round1/. -t $OutDir/gag/round1/genome.sbt -r $OutDir/tbl2asn/round1 -M n -X E -Z $OutDir/gag/round1/discrep.txt -j "[organism=$Organism] [strain=$Strain]"
done
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

```bash
  for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    OutDir="genome_submission/$Organism/$Strain"
    SubmissionID=$(cat genome_submission/Pcac_PRJNA343457_locus_tags.txt | grep "$Strain" | cut -f1 -d ' ')
    echo $SubmissionID
    mkdir -p $OutDir/gag/edited
    ProgDir=/home/armita/git_repos/emr_repos/tools/genbank_submission
    $ProgDir/edit_tbl_file/ncbi_tbl_corrector.py --inp_tbl $OutDir/gag/round1/genome.tbl --inp_val $OutDir/tbl2asn/round1/genome.val --locus_tag $SubmissionID --lab_id $LabID --gene_id "remove" --edits stop pseudo unknown_UTR correct_partial --remove_product_locus_tags "True" --del_name_from_prod "True" --out_tbl $OutDir/gag/edited/genome.tbl
  done
```


## Generating a structured comment detailing annotation methods

```bash
  for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    OutDir="genome_submission/$Organism/$Strain"
    printf "StructuredCommentPrefix\t##Genome-Annotation-Data-START##
    Annotation Provider\tHarrison Lab NIAB-EMR
    Annotation Date\tAUG-2017
    Annotation Version\tRelease 1.01
    Annotation Method\tAb initio gene prediction: Braker 1.9, CodingQuary 2.0 and ORF finding (see publication); Functional annotation: Swissprot (July 2016 release) and Interproscan 5.18-57.0" \
    > $OutDir/gag/edited/annotation_methods.strcmt.txt
  done
```

## Final run of tbl2asn

Following correction of the GAG .tbl file, tbl2asn was re-run to provide the
final genbank submission file.

The options -l paired-ends -a r10k inform how to handle runs of Ns in the
sequence, these options show that paired-ends have been used to estimate gaps
and that runs of N's longer than 10 bp should be labelled as gaps.

```bash
for Assembly in $(ls repeat_masked/P.cactorum/10300/10300_abyss_53_repmask/10300_contigs_unmasked.fa); do
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  echo "$Organism - $Strain"
  OutDir="genome_submission/$Organism/$Strain"
  FinalName="$Organism"_"$Strain"_Armitage_2017
  cp $Assembly $OutDir/gag/edited/genome.fsa
  cp $SbtFile $OutDir/gag/edited/genome.sbt
  mkdir $OutDir/tbl2asn/final
  tbl2asn -p $OutDir/gag/edited/. -t $OutDir/gag/edited/genome.sbt -r $OutDir/tbl2asn/final -M n -X E -Z $OutDir/tbl2asn/final/discrep.txt -j "[organism=$Organism] [strain=$Strain]" -l paired-ends -a r10k -w $OutDir/gag/edited/annotation_methods.strcmt.txt
  cat $OutDir/tbl2asn/final/genome.sqn | sed 's/_pilon//g' | sed 's/title "Saccharopine dehydrogenase.*/title "Saccharopine dehydrogenase/g' | sed 's/"Saccharopine dehydrogenase.*"/"Saccharopine dehydrogenase"/g' > $OutDir/tbl2asn/final/$FinalName.sqn
done
```
