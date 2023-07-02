#!/bin/bash
## Generate samtools pileup files - necessary input to loimpute 

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p samples`
coverage=${1}                                   #Coverage level (e.g., 10x, 5x, 1x)
subsample=subsamp/${sample}.${coverage}.bam     #path to subsample file

module load speedseq/0.1.2
module load samtools/1.9

cat chroms | parallel --verbose -j 20 "samtools mpileup -r {1} -l reference/chr{1}.phased.filt.vcf.gz \
    ${subsample} | gzip -c > pileups/${sample}.${coverage}.{1}.gz"
