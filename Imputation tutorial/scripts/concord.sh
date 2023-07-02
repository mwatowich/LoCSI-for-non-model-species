#!/bin/bash

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p samples`

module load bcftools/1.10.2
module load tabix/0.2.6

bcftools stats -s ${sample} --af-bins 0.01,0.05,0.1,1.0 --af-tag AF \
    imputed_merge/all_imputed.vcf.gz reference/allChroms.phased.filt.vcf.gz | \
    grep 'GCsAF\|GCTs' > concord/unfiltered/${sample}.concord.txt

bcftools stats -f .,PASS -s ${sample} --af-bins 0.01,0.05,0.1,1.0 --af-tag AF \
    imputed_merge/all_imputed.vcf.gz reference/allChroms.phased.filt.vcf.gz | \
    grep 'GCsAF\|GCTs' > concord/filtered/${sample}.concord.txt
