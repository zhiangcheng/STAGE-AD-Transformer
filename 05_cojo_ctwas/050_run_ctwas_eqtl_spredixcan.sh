#!/bin/bash

cond_gwas_files=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/spredixcan/cond_gwas_files)
echo ${cond_gwas_files}

tissues=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/tissue_lists/14tiss)
echo ${tissues}

# Run the SPrediXcan analysis
for STUDY in ${cond_gwas_files}
do
	for MODEL in ${tissues}
	do
		python3 $HOME/TWAS/Breast-cancer-Example/MetaXcan/software/SPrediXcan.py \
		--gwas_file /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/condTWAS_parse/final_metal_ad_kunkle_pgcalz_ukb/${STUDY}.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--model_db_path /media/desk15/iy2120/TWAS/Breast-cancer-Example/data/prediction_model/gtex_v8_eqtl_dbs_mashr/mashr_${MODEL}.db \
		--covariance /media/desk15/iy2120/TWAS/Breast-cancer-Example/data/prediction_model/gtex_v8_eqtl_dbs_mashr/mashr_${MODEL}.txt.gz \
		--keep_non_rsid \
		--additional_output \
		--model_db_snp_key varID \
		--output_file /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/spredixcan/spredixcan_eqtl_mashr/spredixcan_igwas_gtexmashrv8_${STUDY}__PM__${MODEL}.csv
	done
done
