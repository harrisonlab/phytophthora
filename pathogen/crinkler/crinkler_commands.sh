#!/usr/bin/bash

#qsub /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/crinkler/sub_crinkler.sh analysis/rxlr/P.cactorum/10300/10300_nuc.fa LxLFLAK LYLAK HVLVVVP

for inpath in $(ls -d analysis/rxlr_unmasked/*/*); do
	STRAIN=$(echo $inpath | rev | cut -d "/" -f1 | rev)
	qsub /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/crinkler/sub_crinkler.sh $inpath/"$STRAIN"_nuc.fa LxLFLAK LYLAK HVLVVVP
done


MOTIF1=LxLFLAK
MOTIF2=LYLAK
MOTIF3=HVLVVVP

for infolder in $(ls -d analysis/crinkler/*/*); do
	STRAIN=$(echo $infolder | rev | cut -d "/" -f1 | rev)
	cat analysis/rxlr_unmasked/*/*/"$STRAIN"_nuc.fa | grep -A1 -w "$(cat $infolder/findmotif_LxLFLAK.fa $infolder/findmotif_HVLVVVP.fa | grep '>' | cut -f1 | sed 's/ //g' | sort | uniq -d)" > $infolder/findmotif_LxLFLAK_HVLVVVP.fa
	cat analysis/rxlr_unmasked/*/*/"$STRAIN"_nuc.fa | grep -A1 -w "$(cat $infolder/findmotif_LYLAK.fa $infolder/findmotif_HVLVVVP.fa | grep '>' | cut -f1 | sed 's/ //g' | sort | uniq -d)" > $infolder/findmotif_LYLAK_HVLVVVP.fa
done


grep '>' -c analysis/crinkler/*/*/findmotif_*.fa 

#--	unmasked contigs at kmer 41bp assembly
# analysis/crinkler/P.cactorum/10300/findmotif_HVLVVVP.fa:53
# analysis/crinkler/P.cactorum/10300/findmotif_LxLFLAK.fa:9
# analysis/crinkler/P.cactorum/10300/findmotif_LxLFLAK_HVLVVVP.fa:411514
# analysis/crinkler/P.cactorum/10300/findmotif_LYLAK.fa:119
# analysis/crinkler/P.cactorum/10300/findmotif_LYLAK_HVLVVVP.fa:411514
# 
# analysis/crinkler/P.cactorum/404/findmotif_HVLVVVP.fa:11
# analysis/crinkler/P.cactorum/404/findmotif_LxLFLAK.fa:6
# analysis/crinkler/P.cactorum/404/findmotif_LxLFLAK_HVLVVVP.fa:381647
# analysis/crinkler/P.cactorum/404/findmotif_LYLAK.fa:116
# analysis/crinkler/P.cactorum/404/findmotif_LYLAK_HVLVVVP.fa:381647
# 
# analysis/crinkler/P.cactorum/411/findmotif_HVLVVVP.fa:6
# analysis/crinkler/P.cactorum/411/findmotif_LxLFLAK.fa:6
# analysis/crinkler/P.cactorum/411/findmotif_LxLFLAK_HVLVVVP.fa:0
# analysis/crinkler/P.cactorum/411/findmotif_LYLAK.fa:16
# analysis/crinkler/P.cactorum/411/findmotif_LYLAK_HVLVVVP.fa:208235
# 
# analysis/crinkler/P.cactorum/414/findmotif_HVLVVVP.fa:29
# analysis/crinkler/P.cactorum/414/findmotif_LxLFLAK.fa:6
# analysis/crinkler/P.cactorum/414/findmotif_LxLFLAK_HVLVVVP.fa:382810
# analysis/crinkler/P.cactorum/414/findmotif_LYLAK.fa:84
# analysis/crinkler/P.cactorum/414/findmotif_LYLAK_HVLVVVP.fa:382810
# 
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_HVLVVVP.fa:116
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LxLFLAK.fa:122
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LxLFLAK_HVLVVVP.fa:41
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LYAK_HVLVVVP.fa:118139
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LYLAK.fa:120
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LYLAK_HVLVVVP.fa:4
# 
# analysis/crinkler/P.fragariae/SCRP245/findmotif_HVLVVVP.fa:0
# analysis/crinkler/P.fragariae/SCRP245/findmotif_LxLFLAK.fa:20
# analysis/crinkler/P.fragariae/SCRP245/findmotif_LxLFLAK_HVLVVVP.fa:260148
# analysis/crinkler/P.fragariae/SCRP245/findmotif_LYLAK.fa:1
# analysis/crinkler/P.fragariae/SCRP245/findmotif_LYLAK_HVLVVVP.fa:260148
# 
# analysis/crinkler/P.ideai/371/findmotif_HVLVVVP.fa:16
# analysis/crinkler/P.ideai/371/findmotif_LxLFLAK.fa:5
# analysis/crinkler/P.ideai/371/findmotif_LxLFLAK_HVLVVVP.fa:378319
# analysis/crinkler/P.ideai/371/findmotif_LYLAK.fa:70
# analysis/crinkler/P.ideai/371/findmotif_LYLAK_HVLVVVP.fa:378319
# 
# analysis/crinkler/rxlr/P.cactorum/findmotif_LxLFLAK_HVLVVVP.fa:0
# analysis/crinkler/rxlr/P.cactorum/findmotif_LYLAK_HVLVVVP.fa:0


# ---	unmasked contigs at best assembly
# analysis/crinkler/P.cactorum/10300/findmotif_HVLVVVP.fa:70
# analysis/crinkler/P.cactorum/10300/findmotif_LxLFLAK.fa:70
# analysis/crinkler/P.cactorum/10300/findmotif_LxLFLAK_HVLVVVP.fa:19
# analysis/crinkler/P.cactorum/10300/findmotif_LYLAK.fa:123
# analysis/crinkler/P.cactorum/10300/findmotif_LYLAK_HVLVVVP.fa:411514
# 
# analysis/crinkler/P.cactorum/404/findmotif_HVLVVVP.fa:64
# analysis/crinkler/P.cactorum/404/findmotif_LxLFLAK.fa:66
# analysis/crinkler/P.cactorum/404/findmotif_LxLFLAK_HVLVVVP.fa:10
# analysis/crinkler/P.cactorum/404/findmotif_LYLAK.fa:122
# analysis/crinkler/P.cactorum/404/findmotif_LYLAK_HVLVVVP.fa:381647
# 
# analysis/crinkler/P.cactorum/414/findmotif_HVLVVVP.fa:107
# analysis/crinkler/P.cactorum/414/findmotif_LxLFLAK.fa:95
# analysis/crinkler/P.cactorum/414/findmotif_LxLFLAK_HVLVVVP.fa:26
# analysis/crinkler/P.cactorum/414/findmotif_LYLAK.fa:89
# analysis/crinkler/P.cactorum/414/findmotif_LYLAK_HVLVVVP.fa:382810
# 
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_HVLVVVP.fa:0
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LxLFLAK.fa:0
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LxLFLAK_HVLVVVP.fa:527393
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LYAK_HVLVVVP.fa:118139
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LYLAK.fa:0
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LYLAK_HVLVVVP.fa:527393
# 
# analysis/crinkler/P.ideai/371/findmotif_HVLVVVP.fa:50
# analysis/crinkler/P.ideai/371/findmotif_LxLFLAK.fa:44
# analysis/crinkler/P.ideai/371/findmotif_LxLFLAK_HVLVVVP.fa:6
# analysis/crinkler/P.ideai/371/findmotif_LYLAK.fa:75
# analysis/crinkler/P.ideai/371/findmotif_LYLAK_HVLVVVP.fa:5
# 
# analysis/crinkler/P.fragariae/SCRP245/findmotif_HVLVVVP.fa:0
# analysis/crinkler/P.fragariae/SCRP245/findmotif_LxLFLAK.fa:20
# analysis/crinkler/P.fragariae/SCRP245/findmotif_LxLFLAK_HVLVVVP.fa:260148
# analysis/crinkler/P.fragariae/SCRP245/findmotif_LYLAK.fa:1
# analysis/crinkler/P.fragariae/SCRP245/findmotif_LYLAK_HVLVVVP.fa:260148
# 
# analysis/crinkler/P.cactorum/411/findmotif_HVLVVVP.fa:6
# analysis/crinkler/P.cactorum/411/findmotif_LxLFLAK.fa:6
# analysis/crinkler/P.cactorum/411/findmotif_LxLFLAK_HVLVVVP.fa:1
# analysis/crinkler/P.cactorum/411/findmotif_LYLAK.fa:16
# analysis/crinkler/P.cactorum/411/findmotif_LYLAK_HVLVVVP.fa:208235
