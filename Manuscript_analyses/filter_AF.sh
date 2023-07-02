#!/bin/bash
## Filter reference panel to exclude invariant and multiallelic sites
module load bcftools/1.12.0
module load tabix/0.2.6

chrom=`sed -n ${SLURM_ARRAY_TASK_ID}p chroms | awk '{print $1}'`

bcftools filter -e 'INFO/AF = 0 | INFO/AF = 1' chr${chrom}.phased.vcf.gz | \
	bcftools view -m2 -M2 -Oz -o chr${chrom}.phased.filt.vcf.gz

tabix chr${chrom}.phased.filt.vcf.gz
