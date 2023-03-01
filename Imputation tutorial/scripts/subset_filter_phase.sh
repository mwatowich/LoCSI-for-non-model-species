#!/bin/bash                                                                     
chrom=${SLURM_ARRAY_TASK_ID}                                                    
module load htslib/1.9.0                                                        
module load bcftools/1.10.2
module load tabix/0.2.6
module load beagle/4.0

# Filter out samples with high missingness and sites with high fraction missing 
bcftools view -r ${chrom} -S ref_samples_lowMiss --threads 24 \
    --force-samples /data/reference.vcf.gz | \
    bcftools view --threads 24 -i 'F_MISSING<0.1' | \
    bcftools +fill-tags -Oz -o chr${chrom}.filtered.vcf.gz

tabix chr${chrom}.filtered.vcf.gz

# Phase                                                                       
path_beagle=/packages/7x/beagle/4.0/beagle.r1399.jar
java -Xmx50g -jar $path_beagle nthreads=24 impute=false gt=chr${chrom}.filtered.vcf.gz out=mgap.chr${chrom}.phased chrom=${chrom}
tabix mgap.chr${chrom}.phased.vcf.gz
