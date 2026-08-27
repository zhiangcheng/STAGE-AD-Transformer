#!/bin/bash


study_list=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/acat/study_list)
echo ${study_list}


# Run the acat_prep_sqtl analysis
for STUDY in ${study_list}
do
	python3 /media/desk15/iy2120/Genomic_tools/MetaXcan/software/spredixcan_splicing_merge.py \
	--spred_sqtl_pattern "/media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/spredixcan/spredixcan_sqtl_mashr/spredixcan_igwas_gtexmashrv8_${STUDY}*.csv" \
	--output /media/desk15/iy2120/Project2026/Myproject_AD/05_cojo_ctwas/output/intermediate_data/acat/acat_prep_sqtl/${STUDY}_sqtl_acat_input.csv
done
