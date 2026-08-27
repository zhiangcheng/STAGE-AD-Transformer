#!/bin/bash

STUDY_NAME=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/06_focus/input/STUDY_NAME)
echo ${STUDY_NAME}

for STUDY in ${STUDY_NAME}
do
	for CHR in {1..22}
	do
		focus finemap "/media/desk15/iy2120/Project2026/Myproject_AD/06_focus/input/munged_${STUDY}.sumstats.gz" \
		"/media/desk15/iy2120/TWAS/Breast-cancer-Example/pipeline/05_focus/input/1000G_EUR_Phase3_plink_FOCUS/1000G.EUR.QC.${CHR}" \
		"/media/desk15/iy2120/Project2026/Myproject_AD/06_focus/input/combined_dbs/14tiss_eqtl.db" \
		--verbose \
		--tissue Brain_Hippocampus \
		--p-threshold 1e-3 \
		--chr ${CHR} \
		--out "/media/desk15/iy2120/Project2026/Myproject_AD/06_focus/output/eqtl_focus_14tiss_${STUDY}_Chr${CHR}"
	done
done
