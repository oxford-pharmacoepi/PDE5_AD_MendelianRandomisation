#!/bin/bash


imputed_file_dir="/Bulk/Imputation/UKB imputation from genotype"
directory_output="GWAS/Plink_files/"

run_plink="plink2 --bgen ukb22828_c4_b0_v3.bgen ref-first --sample ukb22828_c4_b0_v3.sample\
	--make-bed --out plink_files_c4"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c4_b0_v3.bgen" \
     -iin="${imputed_file_dir}/ukb22828_c4_b0_v3.sample" \
     -icmd="${run_plink}" --tag="Step2" --instance-type "mem1_ssd1_v2_x72"\
     --name "Bgen_to_plink_c4"\
     --destination="${project}:/${directory_output}/" --brief --yes
