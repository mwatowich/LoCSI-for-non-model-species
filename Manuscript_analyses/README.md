## Imputation and analysis pipeline in Watowich et al. 2023 (under review)

### Filter reference panel
1. Remove samples with high missingness, remove sites missing data in >10% of samples, phase (subset_filter_phase.sh)

2. Remove invariant and multiallelic sites (filter_AF.sh)

3. For the gelada samples only, make leave-one-out files: LOO_impute_ref.sh


### Down-sample and impute low-coverage data
1. Make file of the necessary multiplier to achieve each test coverage level
	* Rhesus: cov_x_rhesus
	* Gelada: cov_x_gelada

2. Down-sample bams: subsample.sh

3. Pileups of down-sampled bams: pileup.sh

4. Impute: impute.sh
	* For gelada, use leave-one-out reference files

5. Annotate with MAF from the reference panel
	* get_maf.sh
	* maf_annotate.sh

6. Merge files: merge_imputed.sh

7. Test concordance of imputed data against 'truth' VCF (i.e., the high-coverage VCF for a given sample): concord.sh

8. Population-level analyses using PCA and relatedness: pca_relatedness.sh
	* Do for both the reference and imputed data

9. Compare in silico RADseq to imputed data 
	* inSilico_RADseq.sh
	* downsample_VCFs.sh

10. Analyze concordance per coverage level: imputation_analysis.R
