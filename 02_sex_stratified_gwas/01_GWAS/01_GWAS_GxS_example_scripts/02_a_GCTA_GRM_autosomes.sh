#!/bin/bash
#PBS -l walltime=05:00:00
#PBS -l mem=100GB
#PBS -l ncpus=5

# This script uses the GCTA package to create a GRM using a set of directly genotyped SNPs from the autosomes that have been filtered and LD pruned and using European individuals only
# Uses genotype data of only observed SNPs (directly genotyped) that are on the autosomes
# Filtered and LD pruned (Filtered on MAF of 0.01, geno 0.02, mind 0.02, hwe 0.0000000001 and LD pruning:  window size = 1500kb; step size (variant ct) = 150; r^2 threshold = 0.2)
# A list of individuals with European ancestry is input to only include these individuals when calculating the GRM
# (GCTA GRM calculation is not suitable for cross-ancestry data (inflated relationship coefficients will arise in this case))


### Environment ###

module load GCTA/1.94.1

### Submit script ###

cd /path/

# Then create GRM using this LD pruned set of SNPs

gcta-1.94.1 \
 --bfile /path/directly_genotyped_SNPs \
 --autosome  \
 --keep /path/European_FullIDs.txt \
 --extract /path/Observed_genotypes_GSAchip_autosomes_MAF0.01_Geno0.02_Mind0.02_HWE0.0000000001_LD0.2.prune.in \
 --make-grm  \
 --thread-num 5 \
 --out GRM_autosomes_R11_Euro_GSAchip_QCed_LDpruned

