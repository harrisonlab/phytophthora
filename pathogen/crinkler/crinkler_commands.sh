#!/usr/bin/bash

#qsub /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/crinkler/sub_crinkler.sh analysis/rxlr/P.cactorum/10300/10300_nuc.fa LxLFLAK LYLAK HVLVVVP

for inpath in $(ls -d analysis/rxlr/*/*); do
	STRAIN=$(echo $inpath | rev | cut -d "/" -f1 | rev)
	qsub /home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/crinkler/sub_crinkler.sh $inpath/"$STRAIN"_nuc.fa LxLFLAK LYLAK HVLVVVP
done


MOTIF1=LxLFLAK
MOTIF2=LYLAK
MOTIF3=HVLVVVP

for infolder in $(ls -d analysis/crinkler/*/*); do
	STRAIN=$(echo $infolder | rev | cut -d "/" -f1 | rev)
	cat analysis/rxlr/*/*/"$STRAIN"_nuc.fa | grep -A1 -w "$(cat $infolder/findmotif_LxLFLAK.fa $infolder/findmotif_HVLVVVP.fa | grep '>' | cut -f1 | sort | uniq -d)" > $infolder/findmotif_LxLFLAK_HVLVVVP.fa
	cat analysis/rxlr/*/*/"$STRAIN"_nuc.fa | grep -A1 -w "$(cat $infolder/findmotif_LYLAK.fa $infolder/findmotif_HVLVVVP.fa | grep '>' | cut -f1 | sort | uniq -d)" > $infolder/findmotif_LYLAK_HVLVVVP.fa
done


grep '>' -c analysis/crinkler/*/*/findmotif_*.fa 

# analysis/crinkler/P.cactorum/10300/findmotif_HVLVVVP.fa:70
# analysis/crinkler/P.cactorum/10300/findmotif_LxLFLAK.fa:70
# analysis/crinkler/P.cactorum/10300/findmotif_LxLFLAK_HVLVVVP.fa:20
# analysis/crinkler/P.cactorum/10300/findmotif_LYLAK.fa:123
# analysis/crinkler/P.cactorum/10300/findmotif_LYLAK_HVLVVVP.fa:456656
# 
# analysis/crinkler/P.cactorum/404/findmotif_HVLVVVP.fa:54
# analysis/crinkler/P.cactorum/404/findmotif_LxLFLAK.fa:59
# analysis/crinkler/P.cactorum/404/findmotif_LxLFLAK_HVLVVVP.fa:14
# analysis/crinkler/P.cactorum/404/findmotif_LYLAK.fa:114
# analysis/crinkler/P.cactorum/404/findmotif_LYLAK_HVLVVVP.fa:421018
# 
# analysis/crinkler/P.cactorum/414/findmotif_HVLVVVP.fa:66
# analysis/crinkler/P.cactorum/414/findmotif_LxLFLAK.fa:64
# analysis/crinkler/P.cactorum/414/findmotif_LxLFLAK_HVLVVVP.fa:18
# analysis/crinkler/P.cactorum/414/findmotif_LYLAK.fa:94
# analysis/crinkler/P.cactorum/414/findmotif_LYLAK_HVLVVVP.fa:425387
# 
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_HVLVVVP.fa:117
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LxLFLAK.fa:122
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LxLFLAK_HVLVVVP.fa:41
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LYLAK.fa:134
# analysis/crinkler/P.fragariae/JHVZ02/findmotif_LYLAK_HVLVVVP.fa:12573
# 
# analysis/crinkler/P.ideai/371/findmotif_HVLVVVP.fa:50
# analysis/crinkler/P.ideai/371/findmotif_LxLFLAK.fa:44
# analysis/crinkler/P.ideai/371/findmotif_LxLFLAK_HVLVVVP.fa:10
# analysis/crinkler/P.ideai/371/findmotif_LYLAK.fa:75
# analysis/crinkler/P.ideai/371/findmotif_LYLAK_HVLVVVP.fa:5
