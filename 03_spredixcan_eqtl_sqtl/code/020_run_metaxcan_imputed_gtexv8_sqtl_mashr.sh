#!/bin/bash

gwas_imputed_files=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/03_spredixcan_eqtl_sqtl/temp/gwas_imputed_files)
echo ${gwas_imputed_files}

mashr_sqtl_files=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/03_spredixcan_eqtl_sqtl/temp/mashr_sqtl_files)
echo ${mashr_sqtl_files}

# Run the SPrediXcan analysis
for STUDY in ${gwas_imputed_files}
do
	for MODEL in ${mashr_sqtl_files}
	do
		python3 $HOME/TWAS/Breast-cancer-Example/MetaXcan/software/SPrediXcan.py \
		--gwas_file /media/desk15/iy2120/Project2026/Myproject_AD/00_prep_inputs_for_metaxcan/output/metaxcan_gwas_imputed/${STUDY}.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--model_db_path /media/desk15/iy2120/TWAS/Breast-cancer-Example/data/prediction_model/gtex_v8_sqtl_dbs_mashr/mashr_${MODEL}.db \
		--covariance /media/desk15/iy2120/TWAS/Breast-cancer-Example/data/prediction_model/gtex_v8_sqtl_dbs_mashr/mashr_${MODEL}.txt.gz \
		--keep_non_rsid \
		--additional_output \
		--model_db_snp_key varID \
		--output_file /media/desk15/iy2120/Project2026/Myproject_AD/03_spredixcan_eqtl_sqtl/output/spredixcan_sqtl_mashr/spredixcan_igwas_gtexmashrv8_${STUDY}__PM__${MODEL}.csv
	done
done
