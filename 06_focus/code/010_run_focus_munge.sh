#!/bin/bash

STUDY_NAME=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/06_focus/input/STUDY_NAME)
echo ${STUDY_NAME}

# Case/control numbers double checked by Charlie, used to calculate effective sample sizes for COJO/CTWAS
for STUDY_NAME in ${STUDY_NAME}
do
	focus munge "/media/desk15/iy2120/Project2026/Myproject_AD/06_focus/input/focus_ready_gwas/${STUDY_NAME}.txt.gz" \
	  --N 47525 \
	  --N-cas 1121785 \
	  --output "/media/desk15/iy2120/Project2026/Myproject_AD/06_focus/input/munged_${STUDY_NAME}"
done
