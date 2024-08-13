#!/bin/sh

# The time it takes to run depends on the chromosome you are extracting the snps from

# How to Run:
# ./partD_Step1_Regenie_SelectSNPs.sh in the command line

# What this file does:
# Extracts the SNPs listed in validation.txt from the Imputation files


directory_input="Whole_genome/Breakthrough_gwas"
directory_output="Whole_genome/Breakthrough_gwas/Intermediary_files"
imputed_file_dir="/Bulk/Imputation/UKB imputation from genotype"
phenotype="breakthroughSeverity_validation"
outcome="severity_index"

for chr in {1,2,3,6,10,19}; do
run_regenie_step1="plink2 --bgen ukb22828_c${chr}_b0_v3.bgen ref-first\
        --sample ukb22828_c${chr}_b0_v3.sample\
	-extract validation_${chr}.txt\
	--make-bed\
	--out StepE-selected_snps_${phenotype}_${chr}\
	--no-pheno\
	--keep Initial_input_${phenotype}.phe"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.bgen" \
   -iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.sample"\
   -iin="/${directory_input}/validation_${chr}.txt"\
   -iin="/${directory_input}/${phenotype}.phe"\
   -icmd="${run_regenie_step1}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --name="StepE_${phenotype}_chr${chr}"\
   --destination="/${directory}/" --brief --yes
 
   
done
