#!/bin/bash

tiss_list=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/tiss_list)
echo ${tiss_list}

study_list=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/acat/study_list)
echo ${study_list}

# Run the acat_sqtl analysis
for TISSUE in ${tiss_list}
do
	for STUDY in ${study_list}
	do
		python3 /media/desk15/iy2120/Genomic_tools/MetaXcan/software/SMulTiXcanByFeature.py \
		--gwas_file /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/condTWAS_parse/final_metal_ad_kunkle_pgcalz_ukb/${STUDY}.txt.gz \
		--tiss_list /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/tissue_lists/${TISSUE} \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--model_db_path /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/combine_sqtl.db \
		--keep_non_rsid --model_db_snp_key varID \
		--covariance /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/combine_cross_intron_covar.txt.gz \
		--grouping /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/combine_phenotype_groups.txt.gz GTEx_sQTL \
		--associations /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/acat/acat_prep_sqtl/${STUDY}_sqtl_acat_input.csv SPrediXcan \
		--cutoff_condition_number 30 \
		--verbosity 10 \
		--acat \
		--output /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/acat/acat_sqtl/${TISSUE}__${STUDY}_sqtl_acat_results.tsv
	done
done
