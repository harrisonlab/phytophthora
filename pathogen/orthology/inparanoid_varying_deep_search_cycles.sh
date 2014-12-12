mkdir analysis/inparanoid/summary_tables
cp tmp2/P.cact_orthology.csv analysis/inparanoid/summary_tables/.
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_orthology_tab.pl analysis/inparanoid/summary_tables/all_P.cac_RxLR.txt analysis/inparanoid/10300-414/sqltable.10300-414 analysis/inparanoid/10300-404/sqltable.10300-404 analysis/inparanoid/404-414/sqltable.404-414 > analysis/inparanoid/summary_tables/P.cact_ortho_groups.csv
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/analyse_orthology_tab.pl analysis/inparanoid/summary_tables/P.cact_ortho_groups.csv -file analysis/inparanoid/summary_tables/all_P.cac_RxLR.txt -deep -print_once
mv 10300_Pcac10300_* analysis/inparanoid/summary_tables/orthogroups/.
mv 4* analysis/inparanoid/summary_tables/orthogroups/.
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_final_ortho_tab.pl analysis/inparanoid/summary_tables/orthogroups/ 404 414 10300 > analysis/inparanoid/summary_tables/P.cact_orthology.csv 

# output of analyse orthology tab does not match the original number of RxLRs
# There are 1201 orthology groups in the final outfile including 177 unique genes
# and 1204 orthologs.



# run with 2 deep deep searching cycles
mkdir orthogroups3
cd orthogroups3
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/analyse_orthology_tab.pl ../P.cact_ortho_groups.csv -file ../all_P.cac_RxLR.txt -deep -print_once
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_final_ortho_tab.pl ./ 404 414 10300 > ../new_final_tab.csv
# output is 1199 orthogroups

# run with 3 deep deep searching cycles
mkdir orthogroups3
cd orthogroups3
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/analyse_orthology_tab.pl ../P.cact_ortho_groups.csv -file ../all_P.cac_RxLR.txt -deep -print_once
/home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/build_final_ortho_tab.pl ./ 404 414 10300 > ../new_final_tab.csv
# output stabilised at 1199 orthogroups