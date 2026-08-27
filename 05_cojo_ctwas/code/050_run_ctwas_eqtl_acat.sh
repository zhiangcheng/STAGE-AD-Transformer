#!/bin/bash

tiss_list=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/tiss_list)
echo ${tiss_list}

study_list=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/acat/study_list)
echo ${study_list}

# Run the ACAT-O analysis
for TISSUE in ${tiss_list}
do
	for STUDY in ${study_list}
	do
		python3 /media/desk15/iy2120/Genomic_tools/MetaXcan/software/spredixcan_expression_acat.py \
		--results_dir /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/spredixcan/spredixcan_eqtl_mashr/ \
		--study_pattern "spredixcan_igwas_gtexmashrv8_${STUDY}*" \
		--tiss_list /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/input/tissue_lists/${TISSUE} \
		--output /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/acat/acat_eqtl/${TISSUE}__${STUDY}_eqtl_acat_results.tsv
	done
done
