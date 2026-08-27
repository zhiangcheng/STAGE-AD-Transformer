# Overview
  Generally, this directory follows the steps outlined
  [here.](https://github.com/hakyimlab/summary-gwas-imputation/wiki/GWAS-Harmonization-And-Imputation#gwas-summary-statistics-harmonization)
  That is, GWAS files are: 
    1. Parsed (liftover happens here if necessary)
    2. Summary imputation takes place to impute variants for use in
    splicing/expression models.
    3. Post processing to combine the parsed and imputed results into GWAS data
    that is ready for use in the rest of the study.
*********************************************************************************************************
010_run_metal_bcac-white_ukb.pbs
  * Run metal to combine GWAS of UKB and BCAC.
  * Uses `scripts/metal_bcac-white_ukb.txt`
011_parse_metal_bcac-white_ukb_ovr_brca.pbs
  * Parse the METAL results output from above.

020_parse_bcac_white_gwas_ovr_brca.pbs
  * Parse BCAC Overall Breast Cancer GWAS.
030_parse_bcac_white_gwas_erpos_brca.pbs
  * Parse BCAC ER+ Breast Cancer GWAS.
040_run_metal_bcac-white_erneg_cimba.pbs
  * Run metal to ER- GWAS of BCAC/Cimba.
041_parse_bcac_white_gwas_erneg_cimba_brca.pbs
  * Parse the METAL result output from above.
  * Uses `scripts/metal_bcac-white_erneg_cimba.txt`

050_parse_intrinsic_subtype_1.pbs
050_parse_intrinsic_subtype_2.pbs
050_parse_intrinsic_subtype_3.pbs
050_parse_intrinsic_subtype_4.pbs
050_parse_intrinsic_subtype_5.pbs

090_run_summary_imputation.pbs
#090_run_summary_imputation.sh
#090_summary_imputation.jinja
#090_summary_imputation.yaml
#090_summary_imputation_memory.yaml
  * `source 090*.sh` Uses [badger](https://github.com/hakyimlab/badger) to create and submit multiple jobs to impute GWAS data for every study. 
  # badger 主要用于集群中进行密集数据处理！！

100_post_process_metal_bcac-white_ukb_ovr_brca.pbs
110_post_process_bcac-white_ovr_brca.pbs
120_post_process_bcac-white_erpos_brca.pbs
130_post_process_metal_bcac-white_erneg_cimba.pbs

140_post_process_intrinsic_subtype_1.pbs
140_post_process_intrinsic_subtype_2.pbs
140_post_process_intrinsic_subtype_3.pbs
140_post_process_intrinsic_subtype_4.pbs
140_post_process_intrinsic_subtype_5.pbs
  * These files create the final GWAS ready for use further on in the analysis.

999_add_rsid_col_to_ukb_bcac_metal_results.R
999_add_rsid_col_to_ukb_bcac_metal_results.pbs
  * Add an rsid column to the metal results of UKB/BCAC Overall Breast Cancer.
  * Used as a inspection check, not part of final analysis.


*********************************************************************************************************