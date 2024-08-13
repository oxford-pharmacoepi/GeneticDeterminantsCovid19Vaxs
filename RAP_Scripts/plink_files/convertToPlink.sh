#!/bin/bash


imputed_file_dir="/Bulk/Imputation/UKB imputation from genotype"
directory_output="Whole_genome/Breakthrough_gwas/Plink_files/"

for i in {19..22}; do
run_plink="plink2 --bgen ukb22828_c${i}_b0_v3.bgen ref-first --sample ukb22828_c${i}_b0_v3.sample\
	--make-bed --out plink_files_c${i}"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
     -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
     -icmd="${run_plink}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16"\
     --name "Bgen_to_plink_c${i}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done

for i in {9..18}; do
run_plink="plink2 --bgen ukb22828_c${i}_b0_v3.bgen ref-first --sample ukb22828_c${i}_b0_v3.sample\
	--make-bed --out plink_files_c${i}"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
     -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
     -icmd="${run_plink}" --tag="Step2" --instance-type "mem1_ssd1_v2_x36"\
     --name "Bgen_to_plink_c${i}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done

for i in {3..8}; do
run_plink="plink2 --bgen ukb22828_c${i}_b0_v3.bgen ref-first --sample ukb22828_c${i}_b0_v3.sample\
	--make-bed --out plink_files_c${i}"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
     -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
     -icmd="${run_plink}" --tag="Step2" --instance-type "mem1_ssd1_v2_x72"\
     --name "Bgen_to_plink_c${i}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done

for i in {1..2}; do
run_plink="plink2 --bgen ukb22828_c${i}_b0_v3.bgen ref-first --sample ukb22828_c${i}_b0_v3.sample\
	--make-bed --out plink_files_c${i}"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
     -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
     -icmd="${run_plink}" --tag="Step2" --instance-type "mem1_ssd1_v2_x72"\
     --name "Bgen_to_plink_c${i}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done



