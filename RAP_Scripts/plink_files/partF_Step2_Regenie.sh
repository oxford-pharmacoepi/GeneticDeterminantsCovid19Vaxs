#!/bin/bash


directory_input="Whole_genome/Breakthrough_gwas"
directory_output="Whole_genome/Breakthrough_gwas/Intermediary_files"
imputed_file_dir="Whole_genome/Breakthrough_gwas/Plink_files"
phenotype="oneDose"
outcome="immuneResponse"

for chr in {18..22}; do
  run_regenie_cmd="regenie --step 2 \
    --bed plink_files_c${chr}\
    --out ${phenotype}_assoc.c${chr}\
    --phenoFile Initial_input_${phenotype}.phe --covarFile Initial_input_${phenotype}.phe\
    --bt --approx --firth-se --firth --extract c${chr}_stepE_${phenotype}.snplist\
    --phenoCol ${outcome}\
    --covarCol Sex\
    --covarCol Age\
    --covarCol Genetic_batch\
    --covarCol PC{1:10}\
    --pred ${phenotype}_results_pred.list --bsize 200\
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"
    
  dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${chr}.bed"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.bim"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.fam"\
   -iin="/${directory_output}/c${chr}_stepE_${phenotype}.snplist"\
   -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe"\
   -iin="/${directory_output}/${phenotype}_results_pred.list"\
   -iin="/${directory_output}/${phenotype}_results_1.loco.gz"\
   -icmd="${run_regenie_cmd}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16"\
   --name="StepF_chr${chr}_${phenotype}"\
   --destination="${project}:/${directory_output}/" --brief --yes
done


for chr in {3..17}; do
  run_regenie_cmd="regenie --step 2 \
    --bed plink_files_c${chr}\
    --out ${phenotype}_assoc.c${chr}\
    --phenoFile Initial_input_${phenotype}.phe --covarFile Initial_input_${phenotype}.phe\
    --bt --approx --firth-se --firth --extract c${chr}_stepE_${phenotype}.snplist\
    --phenoCol ${outcome}\
    --covarCol Sex\
    --covarCol Age\
    --covarCol Genetic_batch\
    --covarCol PC{1:10}\
    --pred ${phenotype}_results_pred.list --bsize 200\
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"
    
  dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${chr}.bed"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.bim"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.fam"\
   -iin="/${directory_output}/c${chr}_stepE_${phenotype}.snplist"\
   -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe"\
   -iin="/${directory_output}/${phenotype}_results_pred.list"\
   -iin="/${directory_output}/${phenotype}_results_1.loco.gz"\
   -icmd="${run_regenie_cmd}" --tag="Step2" --instance-type "mem1_ssd1_v2_x36"\
   --name="StepF_chr${chr}_${phenotype}"\
   --destination="${project}:/${directory_output}/" --brief --yes
done

for chr in {1..2}; do
  run_regenie_cmd="regenie --step 2 \
    --bed plink_files_c${chr}\
    --out ${phenotype}_assoc.c${chr}\
    --phenoFile Initial_input_${phenotype}.phe --covarFile Initial_input_${phenotype}.phe\
    --bt --approx --firth-se --firth --extract c${chr}_stepE_${phenotype}.snplist\
    --phenoCol ${outcome}\
    --covarCol Sex\
    --covarCol Age\
    --covarCol Genetic_batch\
    --covarCol PC{1:10}\
    --pred ${phenotype}_results_pred.list --bsize 200\
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"
    
  dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${chr}.bed"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.bim"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.fam"\
   -iin="/${directory_output}/c${chr}_stepE_${phenotype}.snplist"\
   -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe"\
   -iin="/${directory_output}/${phenotype}_results_pred.list"\
   -iin="/${directory_output}/${phenotype}_results_1.loco.gz"\
   -icmd="${run_regenie_cmd}" --tag="Step2" --instance-type "mem1_ssd1_v2_x72"\
   --name="StepF_chr${chr}_${phenotype}"\
   --destination="${project}:/${directory_output}/" --brief --yes
done