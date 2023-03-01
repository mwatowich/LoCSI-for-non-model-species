#!/bin/bash
## low coverage imputation using loimpute

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p samples`
coverage=${1}
chrom=${2}

module load speedseq/0.1.2
module load samtools/1.9
module load singularity/3.6.3
module load bcftools/1.10.2
module load tabix/0.2.6

singularity exec -B /scratch/mwatowic/impute loimpute_latest.sif loimpute \
    -ne 1000 \
    -id ${sample} \
    -i pileups/${sample}.${coverage}.${chrom}.gz \
    -h mgap_new/mgap.chr${chrom}.phased.AFfilt.vcf.gz \
    -o imputed/${sample}.${coverage}.${chrom}

tabix -p vcf imputed/${sample}.${coverage}.${chrom}.vcf.gz
