#!/bin/bash

imputed_file_dir="/Auxiliary_genetic_dataset/Imputation"
directory="MAH/one_dose_validation" #Output directory, example: MAH/cohorts
phenotype="one_dose_validation" #Example: one_dose_cohort

    for chr in 6; do
        run_regenie_cmd="regenie --step 2\
            --bed selected_snps_${phenotype}_${chr}\
            --out ${phenotype}_assoc.c${chr}\
            --phenoFile ${phenotype}.phe\
            --covarFile ${phenotype}.phe\
            --bt\
            --approx\
            --firth-se\
            --firth\
            --phenoCol out\
            --covarCol Sex\
            --covarCol Age\
            --covarCol Genetic_batch\
            --covarCol PC1\
            --covarCol PC2\
            --covarCol PC3\
            --covarCol PC4\
            --covarCol PC5\
            --covarCol PC6\
            --covarCol PC7\
            --covarCol PC8\
            --covarCol PC9\
            --covarCol PC10\
            --pred ${phenotype}_results_pred.list\
            --bsize 200\
            --pThresh 0.05\
            --minMAC 3\
            --threads 16\
            --gz"
            
        dx run swiss-army-knife -iin="/${directory}/selected_snps_${phenotype}_${chr}.bed"\
            -iin="/${directory}/selected_snps_${phenotype}_${chr}.bim"\
            -iin="/${directory}/selected_snps_${phenotype}_${chr}.fam"\
            -iin="/${directory}/${phenotype}.phe"\
            -iin="/${directory}/${phenotype}_results_pred.list"\
            -iin="/${directory}/${phenotype}_results_1.loco.gz"\
            -icmd="${run_regenie_cmd}" --tag="Step2" --instance-type "mem1_hdd1_v2_x16"\
            --name "StepF_${phenotype}_chr${chr}"\
            --destination="${project}:/${directory}/" --brief --yes
    done






