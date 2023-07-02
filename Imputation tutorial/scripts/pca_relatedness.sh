#!/bin/bash

module load plink/1.9.0
module load vcftools/0.1.12b

## Make plink files 
gunzip imputed_merge/all_imputed.vcf
vcftools --vcf imputed_merge/all_imputed.vcf --plink --out ped/imputed

## PCA
plink --file ped/imputed --cluster --pca 50 --out ped/results/gelada_${cov}_pca

## KING relatedness
vcftools --vcf all_imputed.vcf --relatedness2 --out all_imputed
