#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l mem=5GB
#PBS -l ncpus=1

# Run filtering and LD pruning on directly genotyped SNPs on autosome
# Done autosomes only because when include X chr get Warning: --hwe observation counts vary by more than 10%, due to the X chromosome. You may want to use a less stringent --hwe p-value threshol>



### Environment ###

module load plink/1.90b7



### Submit script ###

plink --bfile /path/directly_genotyped_SNPs \
 --extract SNPs_on_autosomes.txt \
 --maf 0.01 \
 --geno 0.02 \
 --mind 0.02 \
 --hwe 0.0000000001 \
 --indep-pairwise 1500 150 0.2 \
 --out Observed_genotypes_GSAchip_autosomes_MAF0.01_Geno0.02_Mind0.02_HWE0.0000000001_LD0.2
