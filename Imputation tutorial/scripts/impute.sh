#!/bin/bash
## Genotype imputation using loimpute

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p samples`

module load speedseq/0.1.2
module load samtools/1.9
module load singularity/3.6.3
module load bcftools/1.10.2
module load tabix/0.2.6

cat chroms | parallel --verbose -j 20 "singularity exec -B /path loimpute_latest.sif loimpute \
    -id ${sample} \
    -i pileups/${sample}.{1}.gz \
    -h reference/chr{1}.phased.filt.vcf.gz \
    -o imputed/${sample}.{1}"

cat chroms | parallel --verbose -j 20 "tabix -p vcf imputed/${sample}.{1}.vcf.gz"
