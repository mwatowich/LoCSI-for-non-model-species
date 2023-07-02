#!/bin/bash

module load bcftools/1.10.2
module load tabix/0.2.6

coverage=`sed -n ${SLURM_ARRAY_TASK_ID}p coverages` 

## Merge all samples per chromosome 
bcftools merge imputed_maf/*_${coverage}_${chr}.vcf.gz | bgzip -c > imputed_merge/all_${coverage}_${chr}.vcf.gz
tabix -p vcf imputed_merge/all_${coverage}_${chr}.vcf.gz

## Merge across all chromosomes
bcftools concat $(echo imputed_merge/all_${coverage}_{1..20}.vcf.gz) | bgzip > imputed_merge/all_imputed_${coverage}.vcf.gz
tabix imputed_merge/all_imputed_${coverage}.vcf.gz
