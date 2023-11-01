## Pipeline for implementing low-coverage sequencing and imputation in natural populations (Watowich et al. 2023). 


### Identify a reference panel for the species of interest
* Based on our analyses in Watowich et al. 2023, we suggest a reference panel of at least 50 unrelated individuals, sequenced to at least 10x coverage. However, special considerations of your species' genetic architecture should be accounted for when developing the reference panel.
* If the quality of genotype imputation from the reference panel is unknown, we suggest testing imputation quality using a leave-one-out approach. Please see the section "Validating reference panel using down-sampling and leave-one-out" below for more information. 
</br>

### Prepare high-confidence reference panel
#### For these steps, we assume a VCF file containing all reference samples is available. Our scripts assume separate VCFs per chromosome, though this format is not necessary

1. Remove samples with high-per-sample missingness (per_sample_missingness.sh)
   * Check the distribution of per-sample missingness from stats output files. Make file of reference individuals with low-missingness (example: data/ref_samples_lowMiss)
   * "High-missingness" will depend on the distribution of missingness in the reference, thus we do not provide a recommended threshold

2. Remove samples with high missingness, remove sites missing data in >10% of samples, phase (subset_filter_phase.sh)
   * We use 10%, but set the threshold accordingly to your reference panel
   * We phase with Beagle (Browning & Browning 2007), but any phasing program of similar performance can be used

3. Optional: Perform standard quality control steps for genetic analyses. For example, removing invariant and multiallelic sites, sites out of HWE, etc. (filter_AF.sh)
</br>

### Impute low-coverage data
1. Create a file of sample names to be used for calling individual file names (example: data/samples)

2. Impute: impute.sh
   * We use loimpue (Wasik 2021; Copyright © 2019-2020, Gencove, Inc.) but use any imputation program of your choice

4. Optional: annotate with MAF from the reference panel. We perform this step in our analyses for Watowich et al 2023 (in prep), as our reference panels were larger than the test imputation datasets. We suggest researchers use BCFtools fill-tags to calculate MAF of imputed data or annotate with reference panel MAF, depending on their population and dataset.
   * get_maf.sh #Note that if multiallelic sites are in the reference panel, this only keeps the first of multi-allelic alleles
   * Make header file (example: AFs.hdr)
   * maf_annotate.sh to annotate imputed files with reference panel MAF

5. Merge imputed VCFs across chromosomes and individuals: merge_imputed.sh

6. Optional: perform any QC (e.g., removing singletons, HWE or LD-filtering)
</br>

### Imputed data are generated. Perform downstream genetic analyses
1. pca_relatedness.sh
* Example 1: population structure analysis using PCA in plink
* Example 2: Calculate relatedness using KING from VCFtools
</br>

### OPTIONAL: Validating reference panel using down-sampling and leave-one-out approach
#### NOTE: we assume bam files for each individual per chromosome

1. Generate in silico low-coverage sequencing
   * Make file of the necessary multiplier to achieve each test coverage level (see example: cov_x)
   * Down-sample bams: subsample.sh
   * Pileups of down-sampled bams: pileup.sh

2. Impute using leave-one out approach
   * Make reference phased VCF for all LOO iterations: LOO_impute_ref.sh
   * Impute data for the left-out individual, using the reference file where they are removed

3. Merge imputed VCFs, perform QC steps as described above

4. Test concordance: concord.sh
   * Compare imputed data to 'truth' VCF, the high-coverage VCF of that animal
   * Note: the 'truth' VCF should only include sites in the reference dataset, i.e., sites that were imputed, otherwise estimates of imputation concordance will be biased
   * If imputed data are concatenated across chromosomes, ensure 'truth' VCF is as well
