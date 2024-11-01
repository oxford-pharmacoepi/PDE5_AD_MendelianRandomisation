#!/bin/bash

directory_output="GWAS/Plink_files/"

run_plink="plink --bfile plink_files_c4 --recode --out ped_files_c4"

dx run swiss-army-knife -iin="GWAS/Plink_files/plink_files_c4.bim" \
     -iin="GWAS/Plink_files/plink_files_c4.bed" \
     -iin="GWAS/Plink_files/plink_files_c4.fam" \
     -icmd="${run_plink}" --tag="Step2" --instance-type "mem1_ssd1_v2_x72"\
     --name "PlinkToPed_c4"\
     --destination="${project}:/${directory_output}/" --brief --yes