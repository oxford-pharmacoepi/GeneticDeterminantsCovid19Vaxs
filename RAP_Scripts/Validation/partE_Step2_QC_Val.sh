#!/bin/bash

directory_input="Whole_genome/Breakthrough_gwas"
directory_output="Whole_genome/Breakthrough_gwas/Intermediary_files"
imputed_file_dir="/Bulk/Imputation/UKB imputation from genotype"
phenotype="breakthroughSusceptibility"
outcome="bt_infection"


for chr in {1,2,3,6,10,19; do
        run_plink_wes="plink2 --bfile StepE-selected_snps_${phenotype}_${chr}\
            --no-pheno\
            --keep ${phenotype}.phe\
            --autosome\
            --maf 0.01\
            --mac 20\
            --geno 0.1\
            --hwe 1e-15\
            --mind 0.1\
            --write-snplist\
            --write-samples\
            --no-id-header\
            --rm-dup force-first\
            --out c${chr}_snps_qc_pass_${phenotype}"

        dx run swiss-army-knife -iin="/${directory}/selected_snps_${phenotype}_${chr}.bed" \
            -iin="/${directory}/selected_snps_${phenotype}_${chr}.bim" \
            -iin="/${directory}/selected_snps_${phenotype}_${chr}.fam"\
            -iin="/${directory}/${phenotype}.phe" \
            -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16"\
            --name "StepE_${phenotype}_chr${chr}"\
            --destination="${directory}/" --brief --yes
done




