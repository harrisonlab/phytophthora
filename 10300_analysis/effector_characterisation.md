## 4. 2 Ananlysis of RxLR effectors

Due to RxLR effectors being predicted from a number of sources the number of
unique RxLRs were identified from motif and Hmm searches within gene models.
221 RxLR effectors were predicted in total from the P.cactorum genome. Of these,
90 were shared between both datasets.

```bash
  for InDir in $(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*); do
    Strain=$(echo "$InDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$InDir" | rev | cut -f2 -d '/' | rev)
    Source="pub"
    if [ $Strain == '10300' ]; then
      Source="Aug"
    fi
    RxLR_motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Species/$Strain/"$Strain"_"$Source"_RxLR_EER_regex.txt)
    RxLR_hmm=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Species/$Strain/"$Strain"_"$Source"_RxLR_hmmer_headers.txt)
    WY_hmm=$(ls analysis/RxLR_effectors/hmmer_WY/$Species/$Strain/"$Strain"_"$Source"_WY_hmmer_headers.txt)
    echo "$Species - $Strain"
    echo "Total number of RxLRs in predicted genes:"
    cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq | wc -l
    echo "Total number of RxLRs shared between prediction sources:"
    cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
    OutDir=analysis/RxLR_effectors/combined_evidence/$Species/$Strain
    mkdir -p $OutDir
    cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq > $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt
    echo "The number of combined RxLR containing proteins containing WY domains are:"
    cat $OutDir/"$Strain"_"$Source"_RxLR_EER_motif_hmm_headers.txt $WY_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
    echo ""
  done
```

```
P.cactorum - 10300
Total number of RxLRs in predicted genes:
145
Total number of RxLRs shared between prediction sources:
88
The number of combined RxLR containing proteins containing WY domains are:
69

P.capsici - LT1534
Total number of RxLRs in predicted genes:
125
Total number of RxLRs shared between prediction sources:
67
The number of combined RxLR containing proteins containing WY domains are:
42

P.infestans - T30-4
Total number of RxLRs in predicted genes:
398
Total number of RxLRs shared between prediction sources:
249
The number of combined RxLR containing proteins containing WY domains are:
157

P.parisitica - 310
Total number of RxLRs in predicted genes:
267
Total number of RxLRs shared between prediction sources:
156
The number of combined RxLR containing proteins containing WY domains are:
108

P.sojae - P6497
Total number of RxLRs in predicted genes:
344
Total number of RxLRs shared between prediction sources:
176
The number of combined RxLR containing proteins containing WY domains are:
132
```

A similar analysis was perfromed with ORF predicted RxLRs:

```bash
  for InDir in $(ls -d analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*); do
    Strain=$(echo "$InDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$InDir" | rev | cut -f2 -d '/' | rev)
    RxLR_motif=$(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Species/$Strain/"$Strain"_ORF_RxLR_EER_regex.txt)
    RxLR_hmm=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Species/$Strain/"$Strain"_ORF_RxLR_hmmer_headers.txt)
    WY_hmm=$(ls analysis/RxLR_effectors/hmmer_WY/$Species/$Strain/"$Strain"_ORF_WY_hmmer_headers.txt)
    echo "$Species - $Strain"
    echo "Total number of RxLRs in predicted ORFs:"
    cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq | wc -l
    echo "Total number of RxLRs shared between prediction sources:"
    cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
    OutDir=analysis/RxLR_effectors/combined_evidence/$Species/$Strain
    mkdir -p $OutDir
    cat $RxLR_motif $RxLR_hmm | cut -f1 -d' ' | sort | uniq > $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm_headers.txt
    echo "The number of combined RxLR containing ORFs containing WY domains are:"
    cat $OutDir/"$Strain"_ORF_RxLR_EER_motif_hmm_headers.txt $WY_hmm | cut -f1 -d' ' | sort | uniq -d | wc -l
    echo ""
  done
```

```
P.cactorum - 10300
Total number of RxLRs in predicted ORFs:
194
Total number of RxLRs shared between prediction sources:
121
The number of combined RxLR containing ORFs containing WY domains are:
74

P.capsici - LT1534
Total number of RxLRs in predicted ORFs:
260
Total number of RxLRs shared between prediction sources:
159
The number of combined RxLR containing ORFs containing WY domains are:
81

P.infestans - T30-4
Total number of RxLRs in predicted ORFs:
447
Total number of RxLRs shared between prediction sources:
254
The number of combined RxLR containing ORFs containing WY domains are:
151

P.parisitica - 310
Total number of RxLRs in predicted ORFs:
342
Total number of RxLRs shared between prediction sources:
211
The number of combined RxLR containing ORFs containing WY domains are:
123

P.sojae - 67593
Total number of RxLRs in predicted ORFs:
322
Total number of RxLRs shared between prediction sources:
188
The number of combined RxLR containing ORFs containing WY domains are:
118
```
