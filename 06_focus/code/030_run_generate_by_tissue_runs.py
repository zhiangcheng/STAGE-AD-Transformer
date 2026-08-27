#!/bin/python

TISSUES = [
'Brain_Amygdala',
'Brain_Anterior_cingulate_cortex_BA24',
'Brain_Caudate_basal_ganglia',
'Brain_Cerebellar_Hemisphere',
'Brain_Cerebellum',
'Brain_Cortex',
'Brain_Frontal_Cortex_BA9',
'Brain_Hippocampus',
'Brain_Hypothalamus',
'Brain_Nucleus_accumbens_basal_ganglia',
'Brain_Putamen_basal_ganglia',
'Brain_Spinal_cord_cervical_c-1',
'Brain_Substantia_nigra',
'Nerve_Tibial'
]

PATTERN="""#!/bin/bash
#SBATCH --job-name=f_{eqtl_sqtl}
#SBATCH --array=1-22%4
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --tasks-per-node=1
#SBATCH --mem=16GB
#SBATCH --output=logs/focus_{eqtl_sqtl}_{cur_tissue}.out

#module load python
#module load pigz
#source activate /project2/guiming/Julian/software/envs/focus

#export PATH=`python -m site --user-base`/bin/:$PATH

cd /media/desk15/iy2120/Project2026/Myproject_AD/06_focus/code

STUDY_NAME=$(cat /media/desk15/iy2120/Project2026/Myproject_AD/06_focus/input/STUDY_NAME)
echo ${{STUDY_NAME}}

for STUDY in ${{STUDY_NAME}}
do
	for CHR in {{1..22}}
	do
		focus finemap "/media/desk15/iy2120/Project2026/Myproject_AD/06_focus/input/munged_${{STUDY}}.sumstats.gz" \\
		"/media/desk15/iy2120/TWAS/Breast-cancer-Example/pipeline/05_focus/input/1000G_EUR_Phase3_plink_FOCUS/1000G.EUR.QC.${{CHR}}" \\
		"/media/desk15/iy2120/TWAS/Breast-cancer-Example/pipeline/05_focus/input/eqtl_sqtl_db/{eqtl_sqtl}_dbs/mashr_{cur_tissue}.db" \\
		--verbose \\
		--tissue {cur_tissue} \\
		--p-threshold 1e-3 \\
		--chr ${{CHR}} \\
		--out "/media/desk15/iy2120/Project2026/Myproject_AD/06_focus/output/by_tissue/{eqtl_sqtl}_focus_14tiss_${{STUDY}}_Chr${{CHR}}__PM__{cur_tissue}"
	done
done
"""

import os
os.makedirs("/media/desk15/iy2120/Project2026/Myproject_AD/06_focus/code/by_tissue_runs", exist_ok=True)
os.makedirs("/media/desk15/iy2120/Project2026/Myproject_AD/06_focus/output/by_tissue", exist_ok=True)
if __name__ == '__main__':
    for cur_tissue in TISSUES:
        for cur_eqtl_sqtl in ['eqtl', 'sqtl']:
            cur_pattern = PATTERN.format(cur_tissue=cur_tissue,
              eqtl_sqtl=cur_eqtl_sqtl) 
            with open('../code/by_tissue_runs/080_focus_{}_{}.sh'.format(cur_eqtl_sqtl, cur_tissue), 'w') as f:
                f.write(cur_pattern)
