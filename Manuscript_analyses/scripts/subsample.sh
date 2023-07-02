#!/bin/bash
##Downsample high coverage files

module load samtools/1.9

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p metadata/cov_x | awk '{print $1}'`
multiple10x=`sed -n ${SLURM_ARRAY_TASK_ID}p metadata/cov_x | awk '{print $3}'`
multiple3x=`sed -n ${SLURM_ARRAY_TASK_ID}p metadata/cov_x | awk '{print $4}'`
multiple1x=`sed -n ${SLURM_ARRAY_TASK_ID}p metadata/cov_x | awk '{print $5}'`
multiple0_5x=`sed -n ${SLURM_ARRAY_TASK_ID}p metadata/cov_x | awk '{print $6}'`
multiple0_1x=`sed -n ${SLURM_ARRAY_TASK_ID}p metadata/cov_x | awk '{print $7}'`

## Subsample to 10x 
samtools view -@ 20 -bh -s ${multiple10x} ${sample}.bam > subsamp/${sample}.10x.bam
samtools index subsamp/${sample}.10x.bam

## Subsample to 3x
samtools view -@ 20 -bh -s ${multiple3x} ${sample}.bam > subsamp/${sample}.3x.bam
samtools index subsamp/${sample}.3x.bam

## Subsample to 1x                                                             
samtools view -@ 20 -bh -s ${multiple1x} ${sample}.bam > subsamp/${sample}.1x.bam
samtools index subsamp/${sample}.1x.bam

## Subsample to 0.5x                                                             
samtools view -@ 20 -bh -s ${multiple0_5x} ${sample}.bam > subsamp/${sample}.0.5x.bam
samtools index subsamp/${sample}.0.5x.bam

## Subsample to 0.1x                                                             
samtools view -@ 20 -bh -s ${multiple0_1x} ${sample}.bam > subsamp/${sample}.0.1x.bam
samtools index subsamp/${sample}.0.1x.bam
