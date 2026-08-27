#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l mem=20GB
#PBS -l ncpus=5

# This script uses the GCTA package to create a sparse GRM from the full GRM, which is then used to run GWAS in fastGWA
# --make-bK-sparse sets the cutoff threshold - entries below this value are set to 0 (so only makes pairs of individuals whose entries in the GRM are greater than this value)
# Default threshold is 0.05. We used threshold of 0.03 as higher thresholds resulted in problems of inflation.

### Environment ###

module load GCTA/1.94.1

### Submit script ###

cd /working/directory/

gcta-1.94.1 \
 --grm GRM_autosomes_R11_Euro_GSAchip_QCed_LDpruned \
 --make-bK-sparse 0.03 \
 --out sparse_GRM_0.03_autosomes_R11_Euro_GSAchip_QCed_LDpruned



