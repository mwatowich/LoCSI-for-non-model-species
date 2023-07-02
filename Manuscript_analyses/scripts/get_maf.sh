#!/bin/bash
## Extract AF from reference VCF, change any AF>0.5 to MAF (1-AF)

chrom=`sed -n ${SLURM_ARRAY_TASK_ID}p chroms`

module load bcftools/1.10.2
module load tabix/0.2.6

bcftools query -f '%CHROM\t%POS\t%POS\t%AF\n' reference/chr${chrom}.phased.filt.vcf.gz | \
    awk '{if ($4 < 0.5) print $0; else print $1,$2,$3,$5=1-$4;}' OFS='\t' | \
    bgzip -c > maf/chr${chrom}.AF.gz
tabix -s1 -b2 -e2 maf/chr${chrom}.AF.gz
