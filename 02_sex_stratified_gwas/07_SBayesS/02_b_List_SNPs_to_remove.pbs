#!/bin/bash
#PBS -l walltime=00:10:00
#PBS -l mem=10GB
#PBS -l ncpus=1

# This script combines all the DENTIST .short outputs from each chromosome to make a full list of SNPs to exclude


### Environment ###



### Preamble ###

directory=/path/12_SBayesS_h2/




### Submit script ###

cd ${directory}


# Take contents from all files with SNPs to remove and concatenate

cat Blokland21_MDD_female_sumstats_chr*.DENTIST.short.txt >> All_chromosomes_SNPs_to_remove_Blokland21_female.txt

