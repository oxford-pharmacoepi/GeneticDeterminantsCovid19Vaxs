#!/bin/sh

# It takes about an hour to run

# How to Run:
# ./partB_Merge.sh in the command line

# What this file does:
# Merge genotype calls

directory="MAH/one_dose_validation" #Path, example: MAH/cohorts
phenotype="one_dose_validation" #Example: one_dose_cohort

run_merge="cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* . ;\
        ls *.bed | sed -e 's/.bed//g'> files_to_merge.txt;\
        plink --merge-list files_to_merge.txt --make-bed\
        --autosome-xy --out ukb22418_c1_22_v2_merged;\
        rm files_to_merge.txt;"

dx run swiss-army-knife -iin="/${directory}/${phenotype}.phe" \
   -icmd="${run_merge}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/${directory}" --brief --yes --name="StepB_${phenotype}"
