#!/bin/bash
#PBS -N fastGWA_GWAS_binary_ChrX
#PBS -l walltime=01:00:00
#PBS -l mem=20GB
#PBS -l ncpus=10

# This script uses fastGWA in the GCTA package running a GWAS using the created sparse GRM
# Subset your pheno file to males only and females only. Then run analysis for females only and males only. (just subsetting pheno file is enough - can use genotype input and qcovar with all info)


### Environment ###

module load GCTA/1.94.1

### Preamble ###

# Edit below variables for your specific GWAS

directory=/path/to/working/directory/ # output files will go here

reference=/path/to/reference/ 

GRM= # name of sparse GRM using - e.g sparse_GRM_0.03_autosomes_R11_Euro

pheno= #name of yor phenotype file (FID IID phenotype(cases = 1, controls = 0). No headers.)

qcovar= # name of your quantitative covariate file (FID IID quantitative covars, e.g. 10 PCs. No headers)

output= # name of output

# Also update the --bfile for each of the X and XY genotype files below
# Use covar if you have any categorial covariates

### Submit script ###

cd ${directory}

#  Save the estimated fastGWA model parameters from an analysis with the autosomes
gcta-1.94.1 --mbfile ${reference}/filenames_autosomes_genotype.txt \
 --grm-sparse ${reference}${GRM} \
 --model-only \
 --fastGWA-mlm-binary \
 --pheno ${pheno} \
 --qcovar ${qcovar} \
 --thread-num 10 \
 --out ${output}_modelonly

# Load the saved model above to run association tests for autosomes
gcta-1.94.1 --mbfile ${reference}/filenames_autosomes_genotype.txt \
 --grm-sparse ${reference}sparse_GRM_0.03_autosomes_R11_Euro_GSAchip_QCed_LDpruned \
 --load-model ${output}_modelonly.fastGWA \
 --thread-num 10 \
 --out ${output}

# Load the saved model above to run association tests for ChrX nPAR
gcta-1.94.1 --bfile /path/chrX_nPAR_genotype_file \
 --grm-sparse ${reference}sparse_GRM_0.03_autosomes_R11_Euro_GSAchip_QCed_LDpruned \
 --load-model ${output}_modelonly.fastGWA \
 --thread-num 10 \
 --out ${output}_nPAR


#  Load the saved model above to run association tests for ChrX PAR
gcta-1.94.1 --bfile /path/chrX_PAR_genotype_file \
 --grm-sparse ${reference}sparse_GRM_0.03_autosomes_R11_Euro_GSAchip_QCed_LDpruned \
 --load-model ${output}_modelonly.fastGWA \
 --thread-num 10 \
 --out ${output}_PAR


