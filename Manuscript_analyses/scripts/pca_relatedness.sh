#!/bin/bash

coverage=`sed -n ${SLURM_ARRAY_TASK_ID}p coverages` 

module load plink/1.9.0
module load vcftools/0.1.12b

## Make plink files 
gunzip imputed_merge/all_imputed_${coverage}.vcf
vcftools --vcf imputed_merge/all_imputed_${coverage}.vcf --plink --out ped/imputed_${coverage}

## PCA
plink --file ped/imputed_${coverage} --cluster --pca 50 --out ped/results/${coverage}_pca

## KING relatedness
vcftools --vcf all_imputed_${coverage}.vcf --relatedness2 --out all_imputed_${coverage}
