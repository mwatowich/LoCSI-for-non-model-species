#!/bin/bash                                                                     

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p samples`     

odule load samtools/1.9                                                        
module load bcftools/1.10.2                                                     
module load tabix/0.2.6                                                         

## Make reference VCF that excludes the sample we are imputing 
cat chroms | awk '{print $1}' | parallel --verbose -j 21 "bcftools view -s ^${sample} --threads 24 \
    reference/chr{1}.phased.filt.vcf.gz -Oz -o LOO_impute_reference/chr{1}.phased.filt.LO${sample}.vcf.gz"

cat chroms | awk '{print $1}' | parallel --verbose -j 21 "tabix LOO_impute_reference/chr{1}.phased.filt.LO${sample}.vcf.gz"
