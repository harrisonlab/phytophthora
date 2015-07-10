<!-- #!/bin/bash -->

# For a comparison between 414 and T30-4

```bash
  WorkDir=/home/groups/harrisonlab/project_files/idris/orthomcl_tmp
  mkdir -p $WorkDir
  cd $WorkDir
```


# Format fasta files

## for 414
```bash
  Taxon_code=P414
  Fasta_file=../gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
```

## for T30-4
```bash
  Taxon_code=T304
  Fasta_file=../assembly/external_group/P.infestans/T30-4/pep/Phytophthora_infestans.ASM14294v1.26.pep.all.fa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
```

# Filter proteins into good and poor sets.

```bash
  Input_dir=./
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=goodProteins.fasta
  Poor_proteins_file=poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

# Perform an all-vs-all blast of the proteins
<!--
blastall -d goodProteins.fasta -p blastp -v 100000 -b 100000 -e 1e-5 -m 8 # -F 'm S' -z protein_database_size
-->
```bash
  BlastDB=my_prot_blast_db

  makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB

  # blastp -db my_prot_blast_db -query goodProteins.fasta -outfmt 6 -out all-vs-all.tsv -num_threads 16
  mkdir blastall
  cd blastall
  BlastOut=all-vs-all_results.tsv
  blastall -d ../$BlastDB -p blastp -i ../$Good_proteins_file -v 100000 -b 100000 -e 1e-5 -m 8 -F 'm S' -1 16 -o $BlastOut
  Fasta_files_dir=compliant_files
  mkdir $Fasta_files_dir
  cp ../$Good_proteins_file $Fasta_files_dir/.
  SimilarSequences=similarSequences.txt
  orthomclBlastParser $BlastOut $Fasta_files_dir >> $SimilarSequences
  ls -lh $SimilarSequences # The database will be 5x the size of this file

  OrthoConfig=~/testing/armita_orthomcl/orthomcl.config
  orthomclLoadBlast $OrthoConfig $SimilarSequences
  Log_file=orthoMCL.log
  orthomclPairs $OrthoConfig $Log_file cleanup=yes #<startAfter=TAG>
  orthomclDumpPairsFiles $OrthoConfig
  MclInput=mclInput
  MclOutput=mclOutput
  mcl $MclInput --abc -I 1.5 -o $MclOutput
  Groups_output=groups.txt
  orthomclMclToGroups $MclOutput $Groups_output
  GitDir=~/git_repos/emr_repos/tools/pathogen/orthology/OrthoMCL
  GroupMatrix=groups_matrix.tab
  $GitDir/orthoMCLgroups2tab $Fasta_files_dir/$Good_proteins_file $Groups_output > $GroupMatrix
```

# Plot venn diagrams:

```R
  mydata <- read.table("groups_matrix.tab")
  transposedata <- t(mydata)
  summary(transposedata)
  # sum(transposedata[, "P414"])
  # sum(transposedata[, "T304"])
  area1=sum(transposedata[, 1])
  area2=sum(transposedata[, 2])
  colname1 <- paste(colnames(transposedata)[1])
  colname2 <- paste(colnames(transposedata)[2])
  label1 <- paste(colname1, ' (', area1, ')', sep="" )
  label2 <- paste(colname2, ' (', area2, ')', sep="" )
  n12=nrow(subset(transposedata, 414==1 & T304==1))
  n12=nrow(subset(transposedata, transposedata[,1]==1 & transposedata[,2]==1))
  pdf('414_vs_T30-4_venn.pdf')
  draw.pairwise.venn(area1, area2, n12,
      category = c(label1, label2), euler.d = TRUE, scaled = TRUE, inverted = FALSE, ext.text = TRUE, ext.percent = rep(0.05, 3),
      lwd = rep(2, 2), lty = rep("solid", 2), col = rep("black", 2), fill = NULL,
      alpha = rep(0.5, 2), label.col = rep("black", 3), cex = rep(1, 3),
      fontface = rep("plain", 3), fontfamily = rep("serif", 3), cat.pos = c(-50, 50),
      cat.dist = rep(0.025, 2), cat.cex = rep(1, 2), cat.col = rep("black", 2),
      cat.fontface = rep("plain", 2), cat.fontfamily = rep("serif", 2),
      cat.just = rep(list(c(0.5, 0.5)), 2), cat.default.pos = "outer",
      cat.prompts = FALSE,
      ext.pos = rep(0, 2), ext.dist = rep(0, 2), ext.line.lty = "solid",
      ext.length = rep(0.95, 2), ext.line.lwd = 1, rotation.degree = 0,
      rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist = 0.05, offset = 0
  )
  dev.off()
  q()
```
