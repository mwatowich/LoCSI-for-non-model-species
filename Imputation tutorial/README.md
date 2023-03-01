## Pipeline for implementing low-coverage sequencing and imputation in natural populations (Watowich 2023, under review). 


### Identify a reference panel for the species of interest
* Based on our analyses in Watowich 2023 (in prep), we suggest a reference panel of at least 50 unrelated individuals, sequenced to at least 10x coverage. However, special considerations of your species' genetic architecture hould be accounted for when developing the reference panel.
* If the quality of genotype imputation from the reference panel is unknown, we suggest testing imputation quality using a leave-one-out approach. Please see "Validating reference panel using down-sampling and leave-one-out" below for more information. 


### Prepare high-confidence reference panel
* For these steps, we assume a VCF file containing all reference samples is available. Our scripts assume separate VCFs per chromosome, though this format is not necessary

1. Remove samples with high-per-sample missingness (per_sample_missingness.sh)
  * Check the distribution of per-sample missingness from stats output files. Make file of reference individuals with low-missingness (example: data/ref_samples_lowMiss)
  * "High-missingness" will depend on the distribution of missingness in the reference, thus we do not provide a recommended threshold

2. Remove samples with high missingness, remove sites missing data in >10% of samples, phase (subset_filter_phase.sh)
  * We use 10%, but set the threshold accordingly to your reference panel
  * We phase with Beagle (Browning & Browning 2007), but any phasing program of similar performance can be used

3. Optional: Perform standard quality control steps for genetic analyses (remove invariant and multiallelic sites, sites out of HWE, etc.)


### Impute low-coverage data
4. Create a file of sample names to be used for calling individual file names (example: data/samples)

5. Impute (impute_single_chrom.sh). We use loimpue (Wasik 2021) but use any imputation program of your choice

6. Optional: annotate with MAF from the reference panel. We perform this step in our analyses for Watowich et al 2023 (in prep), as our reference panels were larger than the test imputation datasets. We suggest researchers use BCFtools fill-tags to calculate MAF of imputed data or annotate with reference panel MAF, depending on their specific population and dataset. 
  * get_maf.sh #Note that if multiallelic sites are in the reference panel, this only keeps the first of multi-allelic alleles
  * Make AFs.hdr, see file for example: ##INFO=<ID=REF_AF,Number=1,Type=Float,Description="Allele frequency in reference genotypes">
  * Run maf_annotate.sh: bcftools annotate imputed files with gelada ref MAF 

7. Concat files / merge 

### Imputed data are generated! Time for downstream genetic analyses. 
8. Example analysis 1: population structure analyses using PCA (make plink files/run plink)

9. Calculate relatedness using KING from VCFtools (relatedness.sh)

### OPTIONAL: Validating reference panel using down-sampling and leave-one-out
* NOTE: requires bam files for each individual per chromosome

1. Make file of the necessary multiplier to achieve each test coverage level (see example: XXXX)

2. Down-sample bams: subsample.sh

3. Pileups of down-sampled bams: pileup.sh

4. Impute using leave-one out approach.
  * Make reference phased VCF for all LOO iterations: LOO_impute_ref.sh.
  * Impute the animal, using the reference file where they are removed: impute_single_chrom_LOO.sh

5. Concatenate VCFs: concat_vcfs.sh

6. Test concordance. Compare imputed to 'truth' VCF, the high-coverage VCF of that animal

7. Analyze concordance per coverage level 
