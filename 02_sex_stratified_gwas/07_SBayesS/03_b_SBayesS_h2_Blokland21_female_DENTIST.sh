#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l mem=350GB
#PBS -l ncpus=15

# This script uses SBayesS in GCTB to calculate SNP-based heritability 
# Using sumstats with problem SNPs identified from DENTIST removed for Blokland21_females (which had problems converging)


### Environment ###

module load gctb/2.5.2


### Preamble ###

directory=/path/12_SBayesS_h2/

LD_matrix=/reference/genepi/public_reference_panels/ldm_ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.mldmlist #chr 1 - 22 doesnt include X chr

GWAS_sumstats=Blokland21_MDD_female_sumstats_formatted_SBayes.ma # need to be in the GCTA-COJO .ma format = header row of SNP A1 A2 freq b se p N (SNP identifier (rsID), the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size)

name=Blokland21_MDD_female_sumstats

remove=All_chromosomes_SNPs_to_remove_Blokland21_female.txt


### Submit script ###

cd ${directory}


gctb --sbayes S \
 --mldm ${LD_matrix} \
 --gwas-summary ${GWAS_sumstats} \
 --thread 15 \
 --chain-length 25000 \
 --burn-in 5000 \
 --num-chains 4 \
 --thin 10 \
 --write-mcmc-txt \
 --seed 37904 \
 --exclude ${remove} \
 --out ${name}_SBayesS_DENTISTremoved_redo
 
# *.parRes gives overall results
# Pi = polygenicity
# hsq = h2 = SNP-based heritability
# S = relationship between MAF and effect size (A negative relationship between effect size and MAF is a signature of negative (or purifying) selection)

# redo is adding flag --write-mcmc-txt to get the MCMC samples of h2
