#!/bin/bash

module load bcftools/1.12.0

bcftools stats -s - --threads 24 reference.chr${SLURM_ARRAY_TASK_ID}.vcf.gz > chr${SLURM_ARRAY_TASK_ID}.stats
