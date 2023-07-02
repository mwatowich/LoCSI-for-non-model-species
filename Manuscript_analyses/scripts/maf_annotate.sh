#!/bin/bash
## Annotate imputed files with reference VCF AF

module load bcftools/1.10.2
module load tabix/0.2.6

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p samples`
coverage=${1}

cat chroms | parallel --verbose -j 20 "bcftools annotate -a maf/chr{1}.AF.gz -h AFs.hdr \
    -c CHROM,POS,POS,REF_AF imputed/${sample}_${coverage}_{1}.vcf.gz | \
    bgzip -c > imputed/${sample}_${coverage}_{1}_MAF.vcf.gz"
cat chroms | parallel --verbose -j 20 "tabix -p vcf imputed_maf/${sample}_${coverage}_{1}.vcf.gz"
