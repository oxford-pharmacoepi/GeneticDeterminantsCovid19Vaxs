#!/bin/sh

# How to Run:
# ./partB.sh on the command line

# What this file does:
# Merge genotype calls

directory_input="Whole_genome/Breakthrough_gwas"
directory_output="Whole_genome/Breakthrough_gwas/Merged_files" 
phenotype="oneDose" 

run_merge="cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* . ;\
        ls *.bed | sed -e 's/.bed//g'> files_to_merge.txt;\
        plink --merge-list files_to_merge.txt --make-bed\
        --autosome-xy --out ukb22418_c1_22_v2_merged;\
        rm files_to_merge.txt;"

dx run swiss-army-knife -iin="/${directory_input}/Initial_input_${phenotype}.phe"\
   -icmd="${run_merge}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/${directory_output}" --brief --yes --name="StepB_${phenotype}"
