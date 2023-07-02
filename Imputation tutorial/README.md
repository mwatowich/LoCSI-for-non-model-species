## Pipeline for implementing low-coverage sequencing and imputation in natural populations (Watowich 2023, under review). 


### Identify a reference panel for the species of interest
* Based on our analyses in Watowich 2023 (in prep), we suggest a reference panel of at least 50 unrelated individuals, sequenced to at least 10x coverage. However, special considerations of your species' genetic architecture should be accounted for when developing the reference panel.
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

3. Optional: Perform standard quality control steps for genetic analyses (remove invariant and multiallelic sites, sites out of HWE, etc.)
</br>

### Impute low-coverage data
1. Create a file of sample names to be used for calling individual file names (example: data/samples)

2. Impute (impute_single_chrom.sh). We use loimpue (Wasik 2021) but use any imputation program of your choice

3. Optional: annotate with MAF from the reference panel. We perform this step in our analyses for Watowich et al 2023 (in prep), as our reference panels were larger than the test imputation datasets. We suggest researchers use BCFtools fill-tags to calculate MAF of imputed data or annotate with reference panel MAF, depending on their specific population and dataset. 
   * get_maf.sh #Note that if multiallelic sites are in the reference panel, this only keeps the first of multi-allelic alleles
   * Make AFs.hdr, see file for example: ##INFO=<ID=REF_AF,Number=1,Type=Float,Description="Allele frequency in reference genotypes">
   * Run maf_annotate.sh: bcftools annotate imputed files with gelada ref MAF 

4. Remove singletons

5. Concat files / merge 
</br>

### Imputed data are generated! Perform downstream genetic analyses
* Example 1: population structure analyses using PCA (make_plink_files/sh and run_plink.sh)

* Example 2: Calculate relatedness using KING from VCFtools (relatedness.sh)
</br>

### OPTIONAL: Validating reference panel using down-sampling and leave-one-out approach
#### NOTE: requires bam files for each individual per chromosome

1. Generate in silico low-coverage sequencing
   * Make file of the necessary multiplier to achieve each test coverage level (see example: XXXX)
   * Down-sample bams: subsample.sh
   * Pileups of down-sampled bams: pileup.sh

2. Impute using leave-one out approach
   * Make reference phased VCF for all LOO iterations: LOO_impute_ref.sh.
   * Impute data for the left-out individual, using the reference file where they are removed: impute_single_chrom_LOO.sh

3. Remove singletons: 

4. Concatenate and merge VCFs: concat_vcfs.sh

5. Test concordance. Compare imputed to 'truth' VCF, the high-coverage VCF of that animal

8. Analyze concordance per coverage level 
