#!/bin/sh

directory_input="Whole_genome/Breakthrough_gwas"
directory_output="Whole_genome/Breakthrough_gwas/Intermediary_files"
phenotype="breakthroughSeverity"
outcome="severity_index"

run_regenie_step1="regenie --step 1\
 --lowmem --out StepD-${phenotype}_results --bed ukb22418_c1_22_v2_merged\
 --phenoFile Initial_input_${phenotype}.phe --covarFile Initial_input_${phenotype}.phe\
 --extract snps_qc_pass_${phenotype}.snplist --phenoCol ${outcome}\
 --covarCol Sex\
 --covarCol Age\
 --covarCol Genetic_batch\
 --covarCol PC{1:10}\
 --bsize 1000 --bt --loocv --gz --threads 16"

dx run swiss-army-knife -iin="/${directory_input}/Merged_files/ukb22418_c1_22_v2_merged.bed" \
   -iin="/${directory_input}/Merged_files/ukb22418_c1_22_v2_merged.bim" \
   -iin="/${directory_input}/Merged_files/ukb22418_c1_22_v2_merged.fam"\
   -iin="/${directory_input}/Intermediary_files/snps_qc_pass_${phenotype}.snplist"\
   -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe" \
   --name="StepD_"${phenotype}\
   -icmd="${run_regenie_step1}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/${directory_output}" --brief --yes
