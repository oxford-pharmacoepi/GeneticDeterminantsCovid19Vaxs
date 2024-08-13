#!/bin/bash

directory_input="Whole_genome/Breakthrough_gwas"
directory_output="Whole_genome/Breakthrough_gwas/Intermediary_files"
imputed_file_dir="Whole_genome/Breakthrough_gwas/Plink_files"
phenotype="oneDose"
outcome="immuneResponse"

for i in {18..22}; do
    run_plink_wes="plink2 --bfile plink_files_c${i}\
      --no-pheno --keep Initial_input_${phenotype}.phe --autosome\
      --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_stepE_${phenotype}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${i}.bim"\
     -iin="${imputed_file_dir}/plink_files_c${i}.bed" \
     -iin="${imputed_file_dir}/plink_files_c${i}.fam" \
     -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe" \
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done

for i in {3..17}; do
    run_plink_wes="plink2 --bfile plink_files_c${i}\
      --no-pheno --keep Initial_input_${phenotype}.phe --autosome\
      --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_stepE_${phenotype}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${i}.bim"\
     -iin="${imputed_file_dir}/plink_files_c${i}.bed" \
     -iin="${imputed_file_dir}/plink_files_c${i}.fam" \
     -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe" \
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x36"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done

for i in {1..2}; do
    run_plink_wes="plink2 --bfile plink_files_c${i}\
      --no-pheno --keep Initial_input_${phenotype}.phe --autosome\
      --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_stepE_${phenotype}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${i}.bim"\
     -iin="${imputed_file_dir}/plink_files_c${i}.bed" \
     -iin="${imputed_file_dir}/plink_files_c${i}.fam" \
     -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe" \
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x72"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done
