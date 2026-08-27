#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l mem=15GB
#PBS -l ncpus=1
#PBS -J 0-6

# This script creates Manhattan and QQ Plots from the Meta-analysis results (after QC and formatting)


### Environment ###

module load R/4.2.0


### Preamble ###

directory=/path/03_Metal/Males/ # This is where the plots will be output

p_threshold=5e-08 # p-value threshold for GWAS (where horizontal line will be placedon manhattan plot): 5e-08 for only common variants or 5e-09 for rare and common variants

reference=/path/Meta-analysis/

files=(/path/03_Metal/Males/Metaanalysis_MDD_male_*_QCed_rsID.txt)

name=$(echo ${files[${PBS_ARRAY_INDEX}]} | sed 's|/path/03_Metal/Males/||' | sed 's|_QCed_rsID.txt||')



### Run script ###

cd ${directory}

Rscript --vanilla ${reference}Plots_Man_QQ_all_SNPs.R ${files[${PBS_ARRAY_INDEX}]} ${name}_QCed ${p_threshold} 

