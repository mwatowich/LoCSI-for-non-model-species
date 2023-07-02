#!/bin/bash

## Filter for sites recapitulating RADseq approach 

module load bcftools/1.10.2
module load tabix/0.2.6

# Filter for sites in the RADseq approach
bcftools view -R SphI_digest_100bp_PE.bed ${coverage}.vcf.gz > ${coverage}_RAD.vcf
