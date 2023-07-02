#!/bin/bash

module load bcftools/1.12.0
bcftools stats -S samples --threads 24 reference.chr1.vcf.gz > chr1.stats
