#!/bin/bash
#PBS -N fastGWA_GEI_binary_autosomes
#PBS -l walltime=01:00:00
#PBS -l mem=80GB
#PBS -l ncpus=10

# This script uses fastGWA in the GCTA package running a Genome wide gene by environment interaction analysis on a binary trait using the sparse GRM for the autosomes (chr 1 - 22)
# Uses fastGWA-mlm which is designed for continuous traits, but no specific method developed for binary traits and should be OK (though not ideal) for binary traits

### Environment ###

module load GCTA/1.94.1

### Preamble ###

# Edit below variables for your specific GWAS

directory=/path/to/working/directory/ # output files will go here

reference=/path/to/reference/ 

GRM= # name of sparse GRM using - e.g sparse_GRM_0.03_autosomes_R11_Euro

pheno= #name of yor phenotype file (FID IID phenotype(cases = 1, controls = 0). No headers.

qcovar= # name of your quantitative covariate file (FID IID quantitative covars, e.g. 10 PCs. No headers)

envir= # name of your environmental variable file (FID IID enviro covars as continuous or 0/1 if binary. No headers. male = 0, female = 1)

output= # name of output

# use covar if have any categorical covariates

### Submit script ###

cd ${directory}


#  Save the estimated fastGWA model parameters from an analysis with the autosomes
gcta-1.94.1 --mbfile ${reference}/fastGWA/filenames_autosomes_genotype.txt  \
 --grm-sparse ${reference}fastGWA/${GRM} \
 --fastGWA-mlm \
 --envir ${envir} \
 --pheno ${pheno} \
 --qcovar ${qcovar} \
 --model-only \
 --thread-num 10 \
 --out ${output}_modelonly
 
# Load the saved model above to run association tests for autosomes
gcta-1.94.1 --mbfile ${reference}fastGWA/filenames_autosomes_genotype.txt \
 --grm-sparse ${reference}fastGWA/${GRM} \
 --load-model ${output}_modelonly.fastGWA \
 --thread-num 10 \
 --out ${output}

# Load the saved model above to run association tests for ChrX nPAR: full dosage compensation
gcta-1.94.1 --bfile /path/chrX_nPAR_genotype_file \
 --grm-sparse ${reference}fastGWA/${GRM} \
 --load-model ${output}_modelonly.fastGWA \
 --dc 1 \
 --thread-num 10 \
 --out ${output}_nPAR_fulldc

# Load the saved model above to run association tests for ChrX nPAR: no dosage compensation
gcta-1.94.1 --bfile /path/chrX_nPAR_genotype_file \
 --grm-sparse ${reference}fastGWA/${GRM} \
 --load-model ${output}_modelonly.fastGWA \
 --dc 0 \
 --thread-num 10 \
 --out ${output}_nPAR_nodc

#  Load the saved model above to run association tests for ChrX PAR
gcta-1.94.1 --bfile /path/chrX_PAR_genotype_file \
 --grm-sparse ${reference}fastGWA/${GRM} \
 --load-model ${output}_modelonly.fastGWA \
 --thread-num 10 \
 --out ${output}_PAR
