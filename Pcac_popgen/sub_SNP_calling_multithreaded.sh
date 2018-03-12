#$ -S /bin/bash
#$ -cwd
#$ -pe smp 24
#$ -l virtual_free=1.25G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

# Testing parallelisation of GATk HaplotypeCaller - may crash. (It did not! Resulted in 2x speedup)
# NOTE: this is a haploid organism. For diploid organism, change "ploidy" argument to 2.
# Changes required in the script:
# VARIABLES
# Reference - the genome reference used in read mapping.
# INSIDE THE GATK command:
# To specify which BAM mapping files (Out1 from pre_SNP_calling_cleanup.sh, RefName ending with "_rg" -> that is, with
# read group added) are to be used in SNP calling, use the -I argument with full path to each file following after that.
# Each new BAM file has to be specified after a separate -I

# Project=/home/groups/harrisonlab/project_files/idris
Project=/scratch/data/armita/idris
OutDir=analysis/popgen/SNP_calling
# Reference=$(ls /home/groups/harrisonlab/project_files/idris/repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)
Reference=$(ls /data/scratch/armita/idris/repeat_masked/P.cactorum/414/filtered_contigs_repmask/414_contigs_unmasked.fa)

RefName=$(basename "$Reference")
Out1="${RefName%.*}_temp.vcf"
Out2="${RefName%.*}.vcf"

ProgDir=/home/sobczm/bin/GenomeAnalysisTK-3.6

java -jar $ProgDir/GenomeAnalysisTK.jar \
     -R $Reference \
     -T HaplotypeCaller \
     -ploidy 2 \
     -nct 24 \
     --allow_potentially_misencoded_quality_scores \
     -I $Project/analysis/popgen/P.cactorum/11-40/11-40_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/12420/12420_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/15_13/15_13_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/15_7/15_7_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/17-21/17-21_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/2003_3/2003_3_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/4032/4032_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/4040/4040_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/404/404_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/414/414_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/415/415_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/416/416_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/62471/62471_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/P295/P295_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/P421/P421_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/PC13_15/PC13_15_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.cactorum/R36_14/R36_14_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.idaei/371/371_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.idaei/SCRP370/SCRP370_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.idaei/SCRP376/SCRP376_vs_414_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -o $Out1

#Break down complex SNPs into primitive ones with VariantsToAllelicPrimitives
#This tool will take an MNP (e.g. ACCCA -> TCCCG) and break it up into separate records for each component part (A-T and A->G).
#This tool modifies only bi-allelic variants.

java -jar $ProgDir/GenomeAnalysisTK.jar \
   -T VariantsToAllelicPrimitives \
   -R $Reference \
   -V $Out1 \
   -o $Out2 \


#####################################
# Notes on GATK parallelisation
#####################################
# http://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
