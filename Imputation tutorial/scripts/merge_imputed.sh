#!/bin/bash

module load bcftools/1.10.2
module load tabix/0.2.6

## Merge all samples per chromosome 
bcftools merge imputed/*_${chr}.vcf.gz | bgzip -c > imputed_merge/all_${chr}.vcf.gz
tabix -p vcf imputed_merge/all_${chr}.vcf.gz

## Merge across all chromosomes
bcftools concat $(echo imputed_merge/all_{1..20}.vcf.gz) | bgzip > imputed_merge/all_imputed.vcf.gz
tabix imputed_merge/all_imputed.vcf.gz
