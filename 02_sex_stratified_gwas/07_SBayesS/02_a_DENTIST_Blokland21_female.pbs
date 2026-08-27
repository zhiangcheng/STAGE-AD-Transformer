#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l mem=50GB
#PBS -l ncpus=10
#PBS -J 1-22

# This script runs DENTIST to make a list of problematic variants in sumstats
# Ran for Blokland female sumstats as wasn't converging


### Environment ###

module load DENTIST/1.3.0.0


### Preamble ###

directory=/path/12_SBayesS_h2/

chr=${PBS_ARRAY_INDEX}

input=Blokland21_MDD_female_sumstats_formatted_SBayes.ma

reference=/path/Meta-analysis/1000G_Phase3_v5/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_rsID_no_missing_dups



### Submit script ###
cd ${directory}


DENTIST --bfile ${reference} \
 --gwas-summary ${input} \
 --thread-num 10 \
 --out ${directory}/Blokland21_MDD_female_sumstats_chr${chr}


# No --keep option to only keep European IDs
# The .DENTIST.full.txt contains the statistics for all the tested variants
# The DENTIST.short.txt contains the rsIDs for variants that cannot pass DENTIST QC and are suggested for removal
