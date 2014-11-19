#!/usr/bin/bash

# crinkler pipe fasta files sometimes contain additional lines of '--' between accessions.
# Need to identify why this happens. For now grep them out.


for directory in $(ls -d analysis/rxlr_hardmasked/P.*/*); do
	STRAIN=$(printf $directory | rev | cut -d '/' -f1 | rev)
	SPECIES=$(printf $directory | rev | cut -d '/' -f2 | rev)
	OUTDIR='analysis/orthology/crinkler-rxlr'
	mkdir -p "$OUTDIR"
	cat analysis/rxlr_hardmasked/$SPECIES/$STRAIN/"$STRAIN"_sp_rxlr.fa | grep '>' > $OUTDIR/"$STRAIN"_crinkler_RxLR_appended.fa
	cat analysis/crinkler/$SPECIES/$STRAIN/findmotif_LxLFLAK_HVLVVVP.fa | grep -v '-' >> $OUTDIR/"$STRAIN"_crinkler_RxLR_appended.fa
done

rm analysis/orthology/crinkler-rxlr/411_crinkler_RxLR_appended.fa 
rm analysis/orthology/crinkler-rxlr/SCRP245_crinkler_RxLR_appended.fa
