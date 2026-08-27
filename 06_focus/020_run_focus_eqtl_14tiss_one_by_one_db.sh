tissue_list=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/06_focus/input/14tiss.txt)
echo ${tissue_list}

for TISSUE in ${tissue_list}
do
	focus import /media/desk15/iy2120/TWAS/Breast-cancer-Example/pipeline/05_focus/input/patched/eqtl_dbs/mashr_${TISSUE}.db predixcan \
	--tissue ${TISSUE} \
	--name GTEx \
	--assay rnaseq \
	--predixcan-method mashr \
	--output /media/desk15/iy2120/Project2026/Myproject_AD/06_focus/input/combined_dbs/14tiss_eqtl
done
